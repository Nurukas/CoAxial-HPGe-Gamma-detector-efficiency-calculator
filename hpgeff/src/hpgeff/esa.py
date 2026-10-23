"""Independent cross-check: effective solid angle (Moens et al., NIM 187 (1981) 451).

For straight (unscattered) photon paths the probability that a photon emitted
at point T in direction u has its first non-coherent interaction inside the
active germanium is

    P(T, u) = sum_k exp(-tau(t_k,in)) * [1 - exp(-mu_Ge * l_k)]

(t_k,in: entry distance of the k-th crossing of the active volume, l_k its
chord, tau: optical depth with coherent scattering excluded). The effective
solid angle is Omega_bar = integral of P over directions, averaged over the
source. Moens' assumption — the "virtual peak-to-total" eps_p / (Omega_bar/4pi)
is the same for every counting geometry — gives the transfer

    eps_x / eps_ref = Omega_bar_x / Omega_bar_ref     (Moens eq. 2)

This module evaluates Omega_bar by Monte Carlo integration over exact ray /
tube intersections — no particle transport, no OpenMC geometry — using the
same photon cross sections. Agreement of the Monte Carlo transfer factor with
this ratio (within the few-% accuracy of the Moens assumption) checks the whole
Monte Carlo set-up independently: positions, units, materials, densities,
source sampling.

Restricted to axisymmetric set-ups built from tube segments (no rounded crystal
edge, flat core bottom, exact cylinders/cups/tubes; point/disk/ring/volume
sources).
"""
from __future__ import annotations

import math
import os
import xml.etree.ElementTree as ET
from functools import lru_cache
from pathlib import Path

import numpy as np

from .spec import ConfigError

N_A = 6.02214076e23


@lru_cache(maxsize=None)
def _photon_lib():
    xs = os.environ.get("OPENMC_CROSS_SECTIONS")
    if not xs:
        raise ConfigError("OPENMC_CROSS_SECTIONS not set (needed for attenuation coefficients)")
    root = ET.parse(xs).getroot()
    base = Path(xs).parent
    out = {}
    for lib in root:
        if lib.get("type") == "photon":
            p = Path(lib.get("path"))
            out[lib.get("materials")] = p if p.is_absolute() else base / p
    return out


@lru_cache(maxsize=None)
def _element_xs(el: str):
    import openmc.data
    ph = openmc.data.IncidentPhoton.from_hdf5(str(_photon_lib()[el]))
    parts = []
    for mt in (504, 515, 517):
        if mt in ph.reactions:
            parts.append(ph.reactions[mt].xs)
    if 522 in ph.reactions:
        parts.append(ph.reactions[522].xs)
    else:
        parts += [ph.reactions[m].xs for m in ph.reactions if 534 <= m <= 599]
    return parts


def mu_no_coherent(mdef, E_eV) -> np.ndarray:
    """Linear attenuation coefficient without coherent scattering (1/cm)."""
    import openmc.data
    E = np.atleast_1d(np.asarray(E_eV, float))
    if mdef is None:
        return np.zeros_like(E)
    mu = np.zeros_like(E)
    for el, w in mdef.mass_fractions.items():
        n = mdef.density * w * N_A / openmc.data.atomic_weight(el) * 1e-24  # atoms / (b cm)
        sig = sum(np.asarray(f(E), float) for f in _element_xs(el))
        mu += n * sig
    return mu


# ----------------------------------------------------------------------------
# Ray / tube-segment intersection
# ----------------------------------------------------------------------------
def _tube_intervals(P, U, r0, r1, z0, z1):
    """Up to two [t_in, t_out] intervals (t >= 0) of rays inside a tube segment."""
    px, py, pz = P[:, 0], P[:, 1], P[:, 2]
    ux, uy, uz = U[:, 0], U[:, 1], U[:, 2]
    with np.errstate(divide="ignore", invalid="ignore"):
        ta = (z0 - pz) / uz
        tb = (z1 - pz) / uz
    zl = np.where(uz != 0, np.minimum(ta, tb), np.where((pz >= z0) & (pz <= z1), -np.inf, np.inf))
    zh = np.where(uz != 0, np.maximum(ta, tb), np.where((pz >= z0) & (pz <= z1), np.inf, -np.inf))
    a = ux * ux + uy * uy
    b = px * ux + py * uy
    c = px * px + py * py

    def roots(R):
        disc = b * b - a * (c - R * R)
        ok = (disc >= 0) & (a > 0)
        sq = np.sqrt(np.where(ok, disc, 0.0))
        with np.errstate(divide="ignore", invalid="ignore"):
            t1 = np.where(ok, (-b - sq) / a, np.inf)
            t2 = np.where(ok, (-b + sq) / a, -np.inf)
        vertical_in = (a == 0) & (c <= R * R)
        t1 = np.where(vertical_in, -np.inf, t1)
        t2 = np.where(vertical_in, np.inf, t2)
        return t1, t2

    o1, o2 = roots(r1)
    lo = np.maximum(np.maximum(zl, o1), 0.0)
    hi = np.minimum(zh, o2)
    if r0 > 0:
        i1, i2 = roots(r0)
        # a line that never enters the inner cylinder excludes nothing: make (i1, i2) empty
        miss = i1 > i2
        i1 = np.where(miss, np.inf, i1)
        i2 = np.where(miss, np.inf, i2)
        A0, A1 = lo, np.minimum(hi, i1)
        B0, B1 = np.maximum(lo, i2), hi
    else:
        A0, A1 = lo, hi
        B0 = B1 = np.zeros_like(lo)
    okA, okB = A1 > A0, B1 > B0  # empty intervals -> [0, 0] (avoids inf - inf)
    A0, A1 = np.where(okA, A0, 0.0), np.where(okA, A1, 0.0)
    B0, B1 = np.where(okB, B0, 0.0), np.where(okB, B1, 0.0)
    return (A0, A1), (B0, B1)


def _len_before(iv, t):
    (A0, A1), (B0, B1) = iv
    return np.clip(np.minimum(A1, t) - A0, 0, None) + np.clip(np.minimum(B1, t) - B0, 0, None)


# ----------------------------------------------------------------------------
# Setup description
# ----------------------------------------------------------------------------
def describe(built):
    """Tube segments (r0, r1, z0, z1, mdef, role) for detector + objects, or raise."""
    det = built.detector
    if det.esa is None:
        raise ConfigError("cross-check needs a detector without edge rounding / round core bottom / "
                          "gas-filled gaps")
    segs = [(a, b, c, d, mdef, role) for (a, b, c, d, role, mdef) in det.esa]
    for o in built.objects:
        if o.kind != "csg" or o.shape is None or o.shape.esa is None:
            raise ConfigError(f"cross-check: object '{o.name}' is not an axisymmetric tube composition")
        md = o.material._hpgeff_def.scaled(o.material.density) if o.material is not None else None
        for (a, b, c, d) in o.shape.esa:
            segs.append((a, b, c, d, md, f"object:{o.name}"))
    env = (0.0, det.info["r_envelope"], det.info["z_floor"], 0.0)
    return segs, env


def _sample_source(built, raw, n, rng):
    k = built.cfg.scale
    typ = raw.get("type")
    if typ == "point":
        p = np.asarray(raw["position"], float) * k
        return np.tile(p, (n, 1))
    if typ in ("disk", "ring"):
        R = float(raw["radius"]) * k
        r0 = float(raw.get("inner_radius", 0.0)) * k if typ == "disk" else R
        c = np.asarray(raw.get("center", [0, 0]), float) * k
        r = np.sqrt(rng.uniform(r0 * r0, R * R, n))
        ph = rng.uniform(0, 2 * math.pi, n)
        return np.column_stack([c[0] + r * np.cos(ph), c[1] + r * np.sin(ph), np.full(n, float(raw["z"]) * k)])
    if typ == "volume":
        o = built.obj(raw["object"])
        if o.shape is None or o.shape.esa is None:
            raise ConfigError("cross-check: volume source object must be an exact axisymmetric shape")
        segs = o.shape.esa
        vols = np.array([math.pi * (b * b - a * a) * (d - c) for a, b, c, d in segs])
        idx = rng.choice(len(segs), n, p=vols / vols.sum())
        out = np.empty((n, 3))
        for i, (a, b, c, d) in enumerate(segs):
            m = idx == i
            nm = int(m.sum())
            r = np.sqrt(rng.uniform(a * a, b * b, nm))
            ph = rng.uniform(0, 2 * math.pi, nm)
            out[m] = np.column_stack([r * np.cos(ph), r * np.sin(ph), rng.uniform(c, d, nm)])
        return out
    raise ConfigError(f"cross-check: source type '{typ}' not supported (point/disk/ring/volume)")


def effective_solid_angle(built, energies_keV, n_rays=2_000_000, chunk=200_000, seed=7):
    """Omega_bar / (4 pi) and its statistical uncertainty for each energy."""
    from .sources import angular_bias_cone
    segs, env = describe(built)
    E = np.asarray(energies_keV, float) * 1e3
    det = built.detector
    amb_def = built.lib.resolve(built.ambient_name)
    gap_def = None
    mu_amb = mu_no_coherent(amb_def, E)
    mu_gap = mu_no_coherent(gap_def, E)
    mus = [mu_no_coherent(s[4], E) for s in segs]
    inside_env = [s[5] in ("active", "dead", "pin") or s[5].startswith("layer:") for s in segs]
    bias = angular_bias_cone(built)
    theta = bias[0] if bias else math.pi
    fac = bias[1] if bias else 1.0
    rng = np.random.default_rng(seed)
    acc = np.zeros(len(E))
    acc2 = np.zeros(len(E))
    done = 0
    while done < n_rays:
        m = min(chunk, n_rays - done)
        P = _sample_source(built, built.source_spec, m, rng)
        mu_c = rng.uniform(math.cos(theta), 1.0, m)
        ph = rng.uniform(0, 2 * math.pi, m)
        st = np.sqrt(1 - mu_c * mu_c)
        U = np.column_stack([st * np.cos(ph), st * np.sin(ph), -mu_c])
        ivs = [_tube_intervals(P, U, *s[:4]) for s in segs]
        env_iv = _tube_intervals(P, U, *env)
        act = [i for i, s in enumerate(segs) if s[5] == "active"]

        def tau_at(t):
            # optical depth along [0, t] for every energy: shape (m, nE)
            tt = t[:, None] * mu_amb[None, :]
            tt = tt + _len_before(env_iv, t)[:, None] * (mu_gap - mu_amb)[None, :]
            for i, iv in enumerate(ivs):
                base = mu_gap if inside_env[i] else mu_amb
                tt = tt + _len_before(iv, t)[:, None] * (mus[i] - base)[None, :]
            return tt

        Pk = np.zeros((m, len(E)))
        for i in act:
            for (t0, t1) in ivs[i]:
                ell = t1 - t0
                has = ell > 0
                if not has.any():
                    continue
                tau = tau_at(np.where(has, t0, 0.0))
                Pk += np.where(has[:, None], np.exp(-tau) * (1 - np.exp(-mus[i][None, :] * ell[:, None])), 0.0)
        acc += Pk.sum(axis=0)
        acc2 += (Pk ** 2).sum(axis=0)
        done += m
    mean = acc / done
    var = np.maximum(acc2 / done - mean ** 2, 0.0) / done
    return fac * mean, fac * np.sqrt(var)
