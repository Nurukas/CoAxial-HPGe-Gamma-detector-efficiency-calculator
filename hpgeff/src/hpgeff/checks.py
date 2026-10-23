"""Pre-flight checks: nuclear data, volumes, overlaps, source sampling, plots."""
from __future__ import annotations

import os
import re
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np

from .regions import contains, mc_volume, sample_box


class CheckReport:
    def __init__(self):
        self.lines: list[str] = []
        self.errors: list[str] = []
        self.warnings: list[str] = []

    def info(self, msg):
        self.lines.append(msg)

    def warn(self, msg):
        self.warnings.append(msg)
        self.lines.append("WARNING: " + msg)

    def error(self, msg):
        self.errors.append(msg)
        self.lines.append("ERROR: " + msg)

    def text(self):
        return "\n".join(self.lines)


def check_nuclear_data(materials, rep: CheckReport):
    xs = os.environ.get("OPENMC_CROSS_SECTIONS")
    if not xs or not Path(xs).is_file():
        rep.error("OPENMC_CROSS_SECTIONS is not set to a cross_sections.xml file")
        return
    root = ET.parse(xs).getroot()
    base = Path(xs).parent
    photon, neutron = {}, {}
    for lib in root:
        typ, path = lib.get("type"), lib.get("path")
        for m in (lib.get("materials") or "").split():
            p = Path(path) if Path(path).is_absolute() else base / path
            (photon if typ == "photon" else neutron if typ == "neutron" else {})[m] = p
    missing_p, missing_n = set(), set()
    for mat in materials:
        if mat is None:
            continue
        for nuc in mat.get_nuclides():
            el = re.match(r"[A-Z][a-z]?", nuc).group(0)
            if el not in photon or not photon[el].is_file():
                missing_p.add(el)
            if nuc not in neutron or not neutron[nuc].is_file():
                missing_n.add(nuc)
    if missing_p:
        rep.error(f"photon data missing for elements {sorted(missing_p)} in {xs}")
    if missing_n:
        rep.error(f"OpenMC also requires neutron data for every nuclide, even for photon-only runs; "
                  f"missing: {sorted(missing_n)}")
    if not missing_p and not missing_n:
        rep.info(f"nuclear data: OK ({xs})")


def transport_smoke(built, rep, out_dir: Path, n=20_000, energy_keV=661.66):
    """Short run with OpenMC geometry debugging (-g): overlapping cells and lost particles."""
    import copy
    from .run import count_lost, read_result, run_openmc
    sim = copy.copy(built.cfg.sim)
    sim.particles, sim.batches, sim.rel_unc_target = n, 2, None
    model, meta = built.make_model(energy_keV, sim)
    wd = Path(out_dir) / f"_transport_{built.name}"
    try:
        sp, log = run_openmc(model, wd, sim, extra_args=["-g"])
    except RuntimeError as exc:
        rep.error(f"transport test (geometry debug) failed: {str(exc)[-600:]}")
        return
    lost = count_lost(log, wd)
    r = read_result(sp, meta)
    msg = (f"transport test: {n} histories at {energy_keV} keV with geometry debugging, "
           f"{lost} lost, FEP ~ {r['fep_eff']:.4g}, histories above E+window: {r['above_peak_counts']}")
    (rep.error if (lost or r["above_peak_counts"]) else rep.info)(msg)
    import shutil
    shutil.rmtree(wd, ignore_errors=True)


def check_setup(built, n_points=600_000, plots_dir: Path | None = None, rep: CheckReport | None = None,
                transport: bool = False):
    rep = rep or CheckReport()
    rng = np.random.default_rng(99)
    det = built.detector
    d = det.info
    k10 = 10.0  # cm -> mm for display
    rep.info(f"=== setup '{built.name}' ===")
    rep.info(f"crystal front face at z = {d['z_front'] * k10:.3f} mm (end-cap outer face is z = 0); "
             f"{'true (open-ended) coaxial' if d['true_coax'] else 'closed-end coaxial' if d['rc'] > 0 else 'no core'}")
    for lay in d["layers"]:
        rep.info(f"  layer {lay['name']:>12s} ({lay['material']}): front {lay['z_in'] * k10:+.3f}..{lay['z_out'] * k10:+.3f} mm, "
                 f"side r {lay['r_in'] * k10:.3f}..{lay['r_out'] * k10:.3f} mm")
    va = d.get("active_volume_cm3")
    ll, ur = det.active_bbox
    vmc, dv = mc_volume(det.active_cell.region, ll - 1e-3, ur + 1e-3, n=n_points, rng=rng)
    if va:
        dev = (vmc - va) / va
        msg = f"active Ge volume: analytic {va:.3f} cm3, sampled {vmc:.3f} ± {dv:.3f} cm3"
        (rep.info if abs(vmc - va) < 4 * dv + 1e-3 * va else rep.error)(msg)
    else:
        rep.info(f"active Ge volume (sampled): {vmc:.3f} ± {dv:.3f} cm3")

    csg_objs = [o for o in built.objects if o.region is not None]
    for o in built.objects:
        rep.info(f"object '{o.name}' [{o.kind}] material={o.material_name} density={o.density} g/cm3 "
                 f"volume={o.volume:.5g} cm3 ({o.volume_source})"
                 + (f" mass={o.density * o.volume:.4f} g" if o.density else ""))
        for n in o.notes:
            rep.info("    " + n)
        if o.region is not None and o.volume_source in ("analytic", "cad"):
            lo, hi = np.asarray(o.bbox[0]), np.asarray(o.bbox[1])
            pad = 1e-3 * (1 + np.max(hi - lo))
            v, dv = mc_volume(o.region, lo - pad, hi + pad, n=n_points, rng=rng)
            ok = abs(v - o.volume) < 4 * dv + 2e-3 * o.volume
            (rep.info if ok else rep.error)(
                f"    region volume check: sampled {v:.5g} ± {dv:.2g} vs {o.volume:.5g} cm3")

    # ---- overlaps ------------------------------------------------------------
    regions = [("detector", det.envelope, det.envelope_bbox)] + [(o.name, o.region, o.bbox) for o in csg_objs]
    for i, (na, ra, ba) in enumerate(regions):
        for nb, rb, bb in regions[i + 1:]:
            lo = np.maximum(ba[0], bb[0])
            hi = np.minimum(ba[1], bb[1])
            if np.any(hi <= lo):
                continue
            pts = sample_box(lo, hi, n_points // 4, rng)
            both = contains(ra, pts) & contains(rb, pts)
            vol = both.mean() * float(np.prod(hi - lo))
            if both.any():
                rep.error(f"overlap between '{na}' and '{nb}': ~{vol:.3g} cm3")
    stl = [o for o in built.objects if o.kind == "stl"]
    if stl:
        from .cad_stl import dagmc_available
        import trimesh
        if not dagmc_available():
            rep.error("STL objects need an OpenMC build with DAGMC (openmc.lib.dagmc_enabled is False). "
                      "Use STEP -> CSG or CSG shapes, or install a DAGMC-enabled OpenMC.")
        for o in stl:
            pts = trimesh.sample.volume_mesh(o.mesh, 50_000)
            for na, ra, _ in regions:
                if len(pts) and contains(ra, pts).any():
                    rep.error(f"STL object '{o.name}' overlaps '{na}'")
    # ---- source ----------------------------------------------------------------
    try:
        model, meta = built.make_model(1000.0)
        rep.info(f"source: {meta.get('source_type')} {meta.get('object', '')} "
                 f"sampling={meta.get('sampling', 'exact')}")
    except Exception as exc:  # pragma: no cover
        rep.error(f"source definition: {exc}")
        model = None
    for w in built.warnings:
        rep.warn(w)
    if transport and model is not None:
        transport_smoke(built, rep, out_dir=plots_dir or built.cfg.output_dir)
    if plots_dir is not None and model is not None:
        try:
            files = make_plots(built, model, plots_dir)
            rep.info(f"plots: {', '.join(str(f) for f in files)}")
        except Exception as exc:
            rep.warn(f"plotting failed: {exc}")
    return rep


def make_plots(built, model, out: Path, pixels=900):
    """Slice plots coloured by material, rendered by OpenMC itself."""
    import shutil
    import subprocess
    import openmc
    out.mkdir(parents=True, exist_ok=True)
    lo = np.min([built.detector.envelope_bbox[0]] + [np.asarray(o.bbox[0]) for o in built.objects], axis=0)
    hi = np.max([built.detector.envelope_bbox[1]] + [np.asarray(o.bbox[1]) for o in built.objects], axis=0)
    zc = float(np.clip(0.0, lo[2], hi[2]))
    ctr = (lo + hi) / 2
    span = float(np.max(hi - lo)) * 1.1
    colors = _material_colors(model.materials)
    plots = []
    p1 = openmc.Plot(name="xz")
    p1.basis, p1.origin, p1.width = "xz", (0.0, 0.0, float(ctr[2])), (span, span)
    p1.pixels = (pixels, pixels)
    p1.color_by, p1.colors, p1.filename = "material", colors, f"{built.name}_xz"
    plots.append(p1)
    zs = [float(o.bbox[0][2] + o.bbox[1][2]) / 2 for o in built.objects] or [zc]
    p2 = openmc.Plot(name="xy")
    p2.basis, p2.origin = "xy", (0.0, 0.0, zs[-1])
    w = float(max(hi[0] - lo[0], hi[1] - lo[1])) * 1.1
    p2.width, p2.pixels = (w, w), (pixels, pixels)
    p2.color_by, p2.colors, p2.filename = "material", colors, f"{built.name}_xy_z{zs[-1] * 10:.1f}mm"
    plots.append(p2)
    model.plots = openmc.Plots(plots)
    wd = out / f"_plot_{built.name}"
    wd.mkdir(parents=True, exist_ok=True)
    model.export_to_model_xml(str(wd / "model.xml"))
    proc = subprocess.run([built.cfg.sim.openmc_exec, "-p", "model.xml"], cwd=wd, stdout=subprocess.PIPE,
                          stderr=subprocess.STDOUT, text=True)
    files = []
    for p in plots:
        for ext in (".png", ".ppm"):
            f = wd / (p.filename + ext)
            if f.exists():
                dest = out / f.name
                shutil.move(str(f), dest)
                files.append(dest)
    if not files:
        raise RuntimeError(proc.stdout[-2000:])
    _legend(model.materials, colors, out / f"{built.name}_legend.png")
    shutil.rmtree(wd, ignore_errors=True)
    return files


_PALETTE = [(31, 119, 180), (255, 127, 14), (44, 160, 44), (214, 39, 40), (148, 103, 189),
            (140, 86, 75), (227, 119, 194), (127, 127, 127), (188, 189, 34), (23, 190, 207)]


def _material_colors(materials):
    cols = {}
    for i, m in enumerate(materials):
        name = (m.name or "").lower()
        if "germanium" in name:
            cols[m] = (70, 70, 200)
        elif name == "ambient":
            cols[m] = (240, 240, 240)
        else:
            cols[m] = _PALETTE[i % len(_PALETTE)]
    return cols


def _legend(materials, colors, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    fig, ax = plt.subplots(figsize=(3.2, 0.3 * len(materials) + 0.6))
    ax.axis("off")
    handles = [Patch(color=np.array(colors[m]) / 255, label=m.name) for m in materials]
    handles.append(Patch(facecolor="white", edgecolor="black", label="vacuum / void"))
    ax.legend(handles=handles, loc="center", frameon=False)
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)
