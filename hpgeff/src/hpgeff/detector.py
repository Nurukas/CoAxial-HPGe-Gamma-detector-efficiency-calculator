"""Parametric coaxial HPGe detector in exact CSG.

Coordinate convention: the detector axis is z; the detector looks up (+z);
``z = 0`` is the outer front face of the outermost layer (the end-cap window),
so samples placed at ``z >= 0`` sit on or above the end cap.

The active germanium is ONE OpenMC cell. OpenMC's pulse-height tally keeps a
separate energy sum per cell, so splitting the active volume into several
cells would split each photon's deposited energy and destroy the full-energy
peak. The region of that single cell can be as complex as needed (rounded front
edge, core hole, dead layers).

Layers (holder cup, window, end cap, ...) are concentric "cups" listed from the
crystal outward; each has a top thickness, a side thickness and the vacuum gaps
to the previous layer. A layer with ``side: 0`` is a flat disk (e.g. a window
foil); one with ``top: 0`` is a tube.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from .materials import MaterialLibrary
from .spec import ConfigError, Spec


@dataclass
class Detector:
    cells: list
    active_cell: object
    envelope: object  # region of the whole detector (outermost cup, solid)
    envelope_bbox: tuple
    active_bbox: tuple
    materials: list
    esa: list | None  # [(r_in, r_out, z0, z1, role, material)] or None if not axisymmetric-tubes
    info: dict = field(default_factory=dict)


def build_detector(s: Spec, k: float, lib: MaterialLibrary) -> Detector:
    import openmc

    cs = s.sub("crystal")
    R = cs.num("radius", positive=True) * k
    L = cs.num("length", positive=True) * k
    rho = cs.num("edge_rounding", 0.0, nonneg=True) * k
    ge_name = cs.get("material", "germanium")
    core = cs.sub("core", {})
    rc = core.num("radius", 0.0, nonneg=True) * k
    depth = core.num("depth", (L / k) if rc > 0 else 0.0, nonneg=True) * k
    bottom = core.choice("bottom", ("flat", "round"), "flat")
    core.done()
    dl = cs.sub("dead_layers", {})
    tf = dl.num("front", 0.0, nonneg=True) * k
    ts = dl.num("side", 0.0, nonneg=True) * k
    tb = dl.num("back", 0.0, nonneg=True) * k
    tc = dl.num("core", 0.0, nonneg=True) * k
    dl.done()
    cs.done()

    # ---------------- sanity checks -----------------------------------------
    true_coax = rc > 0 and depth >= L * (1 - 1e-12)
    closed = rc > 0 and not true_coax
    if rc > 0 and rc + tc >= R - ts:
        raise ConfigError("detector: core radius + core dead layer must be smaller than "
                          "crystal radius - side dead layer")
    if tf + tb >= L:
        raise ConfigError("detector: front + back dead layers exceed the crystal length")
    if closed and depth + tc >= L - tf:
        raise ConfigError("detector: closed-end core (depth + core dead layer) reaches the front dead layer")
    if closed and bottom == "round" and depth < rc:
        raise ConfigError("detector: a round core bottom needs depth >= core radius")
    if rho > min(R, L):
        raise ConfigError("detector: edge_rounding larger than crystal radius/length")
    if rho > 0 and rc > 0 and R - rho < rc + tc:
        raise ConfigError("detector: edge_rounding reaches the core")

    # ---------------- layers: positions --------------------------------------
    raw_layers = s.get("layers", [])
    if not isinstance(raw_layers, list):
        raise ConfigError("detector.layers must be a list")
    layers = []
    for i, raw in enumerate(raw_layers):
        ls = Spec(raw, f"detector.layers[{i}]")
        lay = dict(name=str(ls.get("name", f"layer{i}")), material=ls.get("material"),
                   top=ls.num("top", 0.0, nonneg=True) * k, side=ls.num("side", 0.0, nonneg=True) * k,
                   gap_top=ls.num("gap_top", 0.0, nonneg=True) * k,
                   gap_side=ls.num("gap_side", 0.0, nonneg=True) * k)
        ls.done()
        if lay["top"] == 0 and lay["side"] == 0:
            raise ConfigError(f"detector.layers[{i}] ('{lay['name']}') has zero top and side thickness")
        layers.append(lay)
    if not layers:
        raise ConfigError("detector.layers: at least the end cap must be given (z=0 is its outer face)")
    gap_mat_name = s.get("gap_material", "vacuum")
    back_ext = s.num("back_extension", 20.0, nonneg=True) * k

    pin_raw = s.get("contact_pin", None)  # null / absent -> no pin
    pin_spec = Spec(pin_raw, "detector.contact_pin") if pin_raw else None
    pin = None
    if pin_spec is not None:
        rp = pin_spec.num("radius", positive=True) * k
        lp = pin_spec.num("length", positive=True) * k
        pmat = pin_spec.get("material", "brass")
        pin_spec.done()
        if rc <= 0:
            raise ConfigError("detector.contact_pin given but the crystal has no core")
        if rp > rc:
            raise ConfigError(
                f"detector.contact_pin.radius ({rp / k:g}) is larger than the core radius ({rc / k:g}). "
                f"If the spec sheet value is a diameter, use radius = {rp / k / 2:g}.")
        hole_top = depth if not true_coax else L
        if bottom == "round" and closed:
            hole_top = depth - rc
        if lp > hole_top + 1e-12:
            raise ConfigError("detector.contact_pin.length exceeds the core depth")
        pin = dict(r=rp, length=lp, material=pmat)
    s.done()

    T = sum(l["gap_top"] + l["top"] for l in layers)
    z_ct = -T
    z_cb = z_ct - L
    z_f = z_cb - back_ext
    zprev, rprev = z_ct, R
    for lay in layers:
        lay["z_in"] = zprev + lay["gap_top"]
        lay["z_out"] = lay["z_in"] + lay["top"]
        lay["r_in"] = rprev + lay["gap_side"]
        lay["r_out"] = lay["r_in"] + lay["side"]
        zprev, rprev = lay["z_out"], lay["r_out"]
    layers[-1]["z_out"] = 0.0  # exact

    # ---------------- materials -------------------------------------------------
    mats = []
    ge = lib.make_openmc(ge_name, "detector_germanium")
    mats.append(ge)
    gapmat = lib.make_openmc(gap_mat_name, "detector_gap") if not lib.is_vacuum(gap_mat_name) else None
    if gapmat is not None:
        mats.append(gapmat)

    # ---------------- crystal surfaces & regions -------------------------------
    cyl_R = openmc.ZCylinder(r=R, name="crystal outer radius")
    p_ct = openmc.ZPlane(z0=z_ct, name="crystal front face")
    p_cb = openmc.ZPlane(z0=z_cb, name="crystal back face")
    if rho > 0:
        p_rz = openmc.ZPlane(z0=z_ct - rho)
        cyl_rr = openmc.ZCylinder(r=R - rho)
        tor = openmc.ZTorus(z0=z_ct - rho, a=R - rho, b=rho, c=rho)
        C_out = ((-cyl_R & +p_cb & -p_rz) | (-cyl_rr & +p_rz & -p_ct) | (-tor & +p_rz & +cyl_rr))
    else:
        C_out = -cyl_R & +p_cb & -p_ct
    cyl_Ra = cyl_R if ts == 0 else openmc.ZCylinder(r=R - ts, name="active outer radius")
    p_at = p_ct if tf == 0 else openmc.ZPlane(z0=z_ct - tf, name="active front")
    p_ab = p_cb if tb == 0 else openmc.ZPlane(z0=z_cb + tb, name="active back")
    if rho > ts and tf < rho:
        tor_a = tor if ts == 0 else openmc.ZTorus(z0=z_ct - rho, a=R - rho, b=rho - ts, c=rho - ts)
        C_in = ((-cyl_Ra & +p_ab & -p_rz) | (-cyl_rr & +p_rz & -p_at) | (-tor_a & +p_rz & +cyl_rr & -p_at))
    else:
        C_in = -cyl_Ra & +p_ab & -p_at

    H = H_dl = None
    if rc > 0:
        cyl_rc = openmc.ZCylinder(r=rc, name="core radius")
        cyl_rcd = cyl_rc if tc == 0 else openmc.ZCylinder(r=rc + tc, name="core dead-layer radius")
        if true_coax:
            H, H_dl = -cyl_rc, -cyl_rcd
            M = C_out & +cyl_rc
            A = C_in & +cyl_rcd
        elif bottom == "flat":
            p_h = openmc.ZPlane(z0=z_cb + depth, name="core bottom")
            p_hd = p_h if tc == 0 else openmc.ZPlane(z0=z_cb + depth + tc)
            H, H_dl = (-cyl_rc & -p_h), (-cyl_rcd & -p_hd)
            M = C_out & (+cyl_rc | +p_h)
            A = C_in & (+cyl_rcd | +p_hd)
        else:
            z_s = z_cb + depth - rc
            p_s = openmc.ZPlane(z0=z_s)
            sph = openmc.Sphere(z0=z_s, r=rc)
            sph_d = sph if tc == 0 else openmc.Sphere(z0=z_s, r=rc + tc)
            H = (-cyl_rc & -p_s) | (-sph & +p_s)
            H_dl = (-cyl_rcd & -p_s) | (-sph_d & +p_s)
            M = C_out & ~H
            A = C_in & ~H_dl
    else:
        M, A = C_out, C_in

    D = M & ~A
    cells = []
    active = openmc.Cell(name="active germanium", fill=ge, region=A)
    dead = openmc.Cell(name="dead germanium", fill=ge, region=D)
    cells += [active, dead]

    exclude = [M]
    if pin is not None:
        pmat = lib.make_openmc(pin["material"], "detector_contact_pin")
        if pmat is not None:
            mats.append(pmat)
        cyl_p = openmc.ZCylinder(r=pin["r"], name="contact pin")
        p_pt = openmc.ZPlane(z0=z_cb + pin["length"])
        P = -cyl_p & +p_cb & -p_pt
        cells.append(openmc.Cell(name="contact pin", fill=pmat, region=P))
        exclude.append(P)

    # ---------------- layers ----------------------------------------------------
    p_f = openmc.ZPlane(z0=z_f, name="detector model floor")
    esa_layers = []
    for lay in layers:
        cyl_in = openmc.ZCylinder(r=lay["r_in"])
        p_in = openmc.ZPlane(z0=lay["z_in"])
        p_out = openmc.ZPlane(z0=lay["z_out"])
        if lay["side"] > 0:
            cyl_out = openmc.ZCylinder(r=lay["r_out"])
        if lay["side"] > 0 and lay["top"] > 0:
            reg = (-cyl_out & +p_f & -p_out) & (+cyl_in | +p_in)
            esa_layers += [(0.0, lay["r_out"], lay["z_in"], lay["z_out"]),
                           (lay["r_in"], lay["r_out"], z_f, lay["z_in"])]
        elif lay["side"] > 0:
            reg = +cyl_in & -cyl_out & +p_f & -p_in
            esa_layers += [(lay["r_in"], lay["r_out"], z_f, lay["z_in"])]
        else:
            reg = -cyl_in & +p_in & -p_out
            esa_layers += [(0.0, lay["r_in"], lay["z_in"], lay["z_out"])]
        lay["region"] = reg
        lmat = lib.make_openmc(lay["material"], f"detector_{lay['name']}")
        if lmat is not None:
            mats.append(lmat)
        lay["openmc_material"] = lmat
        cells.append(openmc.Cell(name=f"detector {lay['name']}", fill=lmat, region=reg))
        exclude.append(reg)
        lay["esa"] = esa_layers[-2:] if (lay["side"] > 0 and lay["top"] > 0) else esa_layers[-1:]
    last = layers[-1]
    r_env = max(last["r_out"], last["r_in"])
    cyl_env = openmc.ZCylinder(r=r_env, name="detector envelope")
    p_0 = openmc.ZPlane(z0=0.0, name="end-cap outer face (z=0)")
    envelope = -cyl_env & +p_f & -p_0
    gap_region = envelope
    for r in exclude:
        gap_region = gap_region & ~r
    cells.append(openmc.Cell(name="detector gaps", fill=gapmat, region=gap_region))

    # ---------------- axisymmetric description for the cross-check -------------
    esa = None
    if rho == 0 and not (closed and bottom == "round"):
        gdef = ge._hpgeff_def
        esa = []
        z_h = z_cb + depth if closed else z_ct
        if rc > 0 and true_coax:
            esa.append((rc + tc, R - ts, z_cb + tb, z_ct - tf, "active", gdef))
        elif rc > 0:
            esa.append((rc + tc, R - ts, z_cb + tb, z_h + tc, "active", gdef))
            esa.append((0.0, R - ts, z_h + tc, z_ct - tf, "active", gdef))
        else:
            esa.append((0.0, R - ts, z_cb + tb, z_ct - tf, "active", gdef))
        r_front_in = rc if true_coax else 0.0
        for seg in [(r_front_in, R, z_ct - tf, z_ct), (R - ts, R, z_cb, z_ct - tf),
                    (rc + tc if rc > 0 else 0.0, R - ts, z_cb, z_cb + tb)]:
            if seg[1] > seg[0] and seg[3] > seg[2]:
                esa.append(seg + ("dead", gdef))
        if rc > 0 and tc > 0:
            top_core = (z_ct - tf) if true_coax else z_h
            esa.append((rc, rc + tc, z_cb, top_core, "dead", gdef))
            if closed:
                esa.append((0.0, rc + tc, z_h, z_h + tc, "dead", gdef))
        if pin is not None:
            pdef = lib.resolve(pin["material"])
            if pdef is not None:
                esa.append((0.0, pin["r"], z_cb, z_cb + pin["length"], "pin", pdef))
        for lay in layers:
            ldef = lib.resolve(lay["material"])
            if ldef is None:
                continue
            for seg in lay["esa"]:
                esa.append(tuple(seg) + (f"layer:{lay['name']}", ldef))
        if gapmat is not None:
            esa = None  # gas-filled gaps are not described by tube segments; skip cross-check

    info = dict(R=R, L=L, rc=rc, depth=depth, true_coax=true_coax, edge_rounding=rho,
                dead=dict(front=tf, side=ts, back=tb, core=tc), z_front=z_ct, z_back=z_cb,
                z_floor=z_f, r_envelope=r_env,
                layers=[{kk: v for kk, v in l.items() if kk in ("name", "material", "top", "side",
                                                                "gap_top", "gap_side", "z_in", "z_out",
                                                                "r_in", "r_out")} for l in layers],
                active_volume_cm3=_active_volume(R, L, rc, depth, true_coax, bottom, rho, tf, ts, tb, tc),
                pin=pin)
    return Detector(cells=cells, active_cell=active, envelope=envelope,
                    envelope_bbox=(np.array([-r_env, -r_env, z_f]), np.array([r_env, r_env, 0.0])),
                    active_bbox=(np.array([-(R - ts), -(R - ts), z_cb + tb]),
                                 np.array([R - ts, R - ts, z_ct - tf])),
                    materials=mats, esa=esa, info=info)


def _active_volume(R, L, rc, depth, true_coax, bottom, rho, tf, ts, tb, tc):
    """Analytic active volume (cm3), used to verify the CSG region."""
    Ra, La = R - ts, L - tf - tb
    v = math.pi * Ra * Ra * La
    # rounded front edge: remove the ring that is outside the torus
    if rho > ts and tf < rho:
        a = rho - ts  # minor radius of active corner
        Rm = R - rho  # major radius
        # corner block between r in [Rm, Rm+a], z in [zc, zc+a] minus quarter-torus section
        block = math.pi * ((Rm + a) ** 2 - Rm ** 2) * a
        # quarter torus volume (Pappus): area pi a^2/4, centroid radius Rm + 4a/(3pi)
        qt = (math.pi * a * a / 4.0) * 2 * math.pi * (Rm + 4 * a / (3 * math.pi))
        # the block extends up to z_ct - rho + a = z_ct - ts; the cylinder top is z_ct - tf
        # only valid when tf <= ts (block fully inside cylinder); otherwise approximate
        if tf <= ts:
            v -= (block - qt)
            # band between z_ct-ts and z_ct-tf at r in [Rm, Ra] is outside torus
            v -= math.pi * (Ra ** 2 - Rm ** 2) * (ts - tf)
        else:
            return None
    if rc > 0:
        r = rc + tc
        if true_coax:
            v -= math.pi * r * r * La
        elif bottom == "flat":
            h = depth + tc - tb
            v -= math.pi * r * r * h
        else:
            h_cyl = depth - rc - tb
            v -= math.pi * r * r * h_cyl + (2.0 / 3.0) * math.pi * r ** 3
    return v
