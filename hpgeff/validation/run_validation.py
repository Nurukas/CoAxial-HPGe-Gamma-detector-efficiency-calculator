"""Validation studies for hpgeff (writes validation/results/*.json and figures).

    python validation/run_validation.py [study ...]

Studies
  esa       Monte Carlo transfer factors vs the independent Moens effective-solid-angle integral
  step      STEP -> exact CSG (GEOUNED) vs native CSG, same container and sample
  facets    triangulation ("FED") study: cylinders replaced by inscribed n-gons, with and
            without mass-conserving density
  window    full-energy window: exact deposition vs +-0.1 / 0.5 / 1 keV
  bias      angular biasing on/off
  model     detector-model assumptions: absolute efficiency vs transfer factor
"""
from __future__ import annotations

import copy
import json
import math
import os
import sys
import tempfile
import time
from pathlib import Path

import numpy as np
import yaml

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT / "src"))

from hpgeff.config import parse_config  # noqa: E402
from hpgeff.model import build_setup  # noqa: E402
from hpgeff.run import read_result, run_openmc  # noqa: E402

OUT = HERE / "results"
WORK = Path(os.environ.get("HPGEFF_VALIDATION_WORK", Path(tempfile.gettempdir()) / "hpgeff_validation"))
BASE_YAML = ROOT / "examples" / "demo_detector2.yaml"
E5 = [46.54, 59.54, 122.06, 661.66, 1332.49]


def base_data():
    d = yaml.safe_load(BASE_YAML.read_text())
    d.pop("reference", None)
    return d


def run_point(data, setup, energies, particles=1_000_000, seed=1, tag="x", extra_windows=None):
    """Run one setup at several energies; returns list of result dicts."""
    data = copy.deepcopy(data)
    data["simulation"].update(particles=int(particles), batches=20, seed=seed)
    data["output_dir"] = str(WORK / tag)
    cfg = parse_config(data, BASE_YAML.parent)
    built = build_setup(cfg, setup)
    out = []
    for E in energies:
        model, meta = built.make_model(E)
        if extra_windows:
            import openmc
            for w in extra_windows:
                t = openmc.Tally(name=f"fep_w{w}")
                t.filters = [openmc.CellFilter([built.detector.active_cell]),
                             openmc.EnergyFilter([E * 1e3 - w * 1e3, E * 1e3 + w * 1e3])]
                t.scores = ["pulse-height"]
                model.tallies.append(t)
        sp, log = run_openmc(model, WORK / tag / f"{setup}_{E}", cfg.sim)
        r = read_result(sp, meta)
        if extra_windows:
            import openmc
            import warnings
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                with openmc.StatePoint(str(sp)) as f:
                    for w in extra_windows:
                        r[f"fep_w{w}"] = float(f.get_tally(name=f"fep_w{w}").mean.ravel()[0])
        r["setup"] = setup
        out.append(r)
        print(f"   {tag:>22s} {setup:>20s} {E:8.2f} keV  eff={r['fep_eff']:.5e} ± {r['fep_unc'] / r['fep_eff'] * 100:.2f}%",
              flush=True)
    return out


def run_paired(data, setup, energies, particles, batches, seed, tag):
    """Like run_point, but keeps every batch's FEP value (statepoint per batch) so that runs
    with the same seed can be compared batch by batch (correlated sampling)."""
    import glob
    import warnings
    import openmc
    data = copy.deepcopy(data)
    data["simulation"].update(particles=int(particles), batches=int(batches), seed=seed)
    data["output_dir"] = str(WORK / tag)
    cfg = parse_config(data, BASE_YAML.parent)
    built = build_setup(cfg, setup)
    out = []
    for E in energies:
        model, meta = built.make_model(E)
        model.settings.statepoint = {"batches": list(range(1, int(batches) + 1))}
        wd = WORK / tag / f"{setup}_{E}"
        for f in glob.glob(str(wd / "statepoint.*.h5")):
            os.remove(f)
        sp, log = run_openmc(model, wd, cfg.sim)
        r = read_result(sp, meta)
        sums = []
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            for b in range(1, int(batches) + 1):
                f = wd / f"statepoint.{b:0{len(str(batches))}d}.h5"
                if not f.exists():
                    f = Path(sorted(glob.glob(str(wd / f"statepoint.*{b}.h5")))[0])
                with openmc.StatePoint(str(f)) as spb:
                    sums.append(float(spb.get_tally(name="fep").sum.ravel()[0]))
        per_batch = np.diff([0.0] + sums)
        r["per_batch"] = per_batch.tolist()
        r["setup"] = setup
        out.append(r)
        print(f"   {tag:>22s} {E:8.2f} keV  eff={r['fep_eff']:.5e}", flush=True)
    return out


def save(name, obj):
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / f"{name}.json").write_text(json.dumps(obj, indent=1, default=float))


def load(name):
    p = OUT / f"{name}.json"
    return json.loads(p.read_text()) if p.exists() else None


# ----------------------------------------------------------------------------- studies
def study_esa():
    from hpgeff.esa import effective_solid_angle
    from hpgeff.run import read_csv
    res_dir = ROOT / "examples" / "results_demo2"
    data = base_data()
    cfg = parse_config(data, BASE_YAML.parent)
    rows = {}
    for setup in ("reference_point_5mm", "glass_in_pmma_cup"):
        mc = read_csv(res_dir / setup / "mc_efficiency.csv")
        E = [r["energy_keV"] for r in mc]
        om, du = effective_solid_angle(build_setup(cfg, setup), E, n_rays=4_000_000)
        rows[setup] = [dict(E=e, mc=r["fep_eff"], mc_u=r["fep_unc"], esa=o, esa_u=d)
                       for e, r, o, d in zip(E, mc, om, du)]
    out = []
    for a, b in zip(rows["glass_in_pmma_cup"], rows["reference_point_5mm"]):
        tm, te = a["mc"] / b["mc"], a["esa"] / b["esa"]
        utm = tm * math.hypot(a["mc_u"] / a["mc"], b["mc_u"] / b["mc"])
        out.append(dict(E=a["E"], transfer_mc=tm, transfer_mc_u=utm, transfer_esa=te,
                        ratio=tm / te, vpt_sample=a["mc"] / a["esa"], vpt_ref=b["mc"] / b["esa"]))
    save("esa", out)


def study_step():
    fc = os.environ.get("HPGEFF_FREECAD_PYTHON")
    if not fc:
        print("   skipped: set HPGEFF_FREECAD_PYTHON to a Python that can import FreeCAD and geouned")
        return
    data = base_data()
    data["setups"] = {
        "csg": copy.deepcopy(data["setups"]["glass_in_pmma_cup"]),
        "step": {"objects": [{"name": "cad", "shape": "step", "file": "cad/cup_sample.step",
                              "solids": [{"name": "cup", "material": "pmma"},
                                         {"name": "glass_sample", "material": "glass", "mass_g": 390.864}],
                              "placement": {"align": "bottom_on_endcap"},
                              "converter": {"python": fc}}],
                 "source": {"type": "volume", "object": "glass_sample"}},
    }
    a = run_point(data, "csg", E5, 2_000_000, seed=11, tag="step_csg")
    b = run_point(data, "step", E5, 2_000_000, seed=12, tag="step_step")
    save("step", dict(csg=a, step=b))


def study_facets():
    data = base_data()
    res = {}
    exact = copy.deepcopy(data)
    res["exact"] = run_point(exact, "glass_in_pmma_cup", E5, 1_000_000, seed=21, tag="facets_exact")
    for n in (8, 16, 32, 64):
        for corr in (False, True):
            d = copy.deepcopy(data)
            cup = d["setups"]["glass_in_pmma_cup"]["objects"][0]
            cup["facets"] = n
            if not corr:  # nominal density: the triangulated sample simply has less glass
                cup["fill"].pop("mass_g")
                cup["fill"]["density"] = 2.43
            key = f"n{n}_{'mass' if corr else 'nominal'}"
            res[key] = run_point(d, "glass_in_pmma_cup", E5, 1_000_000, seed=21, tag=f"facets_{key}")
    save("facets", res)


def study_facets_paired():
    """Correlated-sampling faceting study: every configuration samples the source by the same
    rejection scheme with the same seed, so per-batch differences cancel most statistical noise."""
    data = base_data()
    data["setups"]["glass_in_pmma_cup"]["source"]["sampling"] = "rejection"
    E = [59.54, 661.66, 1332.49]
    res = load("facets_paired") or {}
    configs = [("exact", None, True)] + [(f"n{n}_{c}", n, c == "mass") for n in (16, 32, 64, 128)
                                         for c in ("nominal", "mass")]
    for key, n, mass in configs:
        if key in res:
            continue
        d = copy.deepcopy(data)
        cup = d["setups"]["glass_in_pmma_cup"]["objects"][0]
        if n:
            cup["facets"] = n
        if not mass:
            cup["fill"].pop("mass_g")
            cup["fill"]["density"] = 2.43
        res[key] = run_paired(d, "glass_in_pmma_cup", E, 2_000_000, 20, 21, f"fp_{key}")
        save("facets_paired", res)


def study_window():
    data = base_data()
    ws = [0.1, 0.5, 1.0]
    E = [46.54, 59.54, 88.03, 122.06, 661.66]
    a = run_point(data, "reference_point_5mm", E, 2_000_000, seed=31, tag="win_ref", extra_windows=ws)
    b = run_point(data, "glass_in_pmma_cup", E, 2_000_000, seed=32, tag="win_glass", extra_windows=ws)
    save("window", dict(ref=a, glass=b, windows=ws))


def study_bias():
    data = base_data()
    res = {}
    for bias in ("off", "auto"):
        d = copy.deepcopy(data)
        d["simulation"]["angular_bias"] = bias
        for setup in ("reference_point_5mm", "glass_in_pmma_cup"):
            res[f"{setup}_{bias}"] = run_point(d, setup, [59.54, 661.66], 2_000_000,
                                               seed=41 if bias == "off" else 42, tag=f"bias_{bias}")
    save("bias", res)


MODEL_VARIANTS = {
    "nominal": [],
    "closed_end_55mm": [("detector.crystal.core.depth", 55.0)],
    "no_BeO": [("detector.layers", "DROP_BEO")],
    "no_pin": [("detector.contact_pin", None)],
    "side_dead_1.0mm": [("detector.crystal.dead_layers.side", 1.0)],
    "front_dead_0.3mm": [("detector.crystal.dead_layers.front", 0.3)],
    "endcap_gap_6mm": [("detector.layers.2.gap_top", 6.0)],
}


def _apply(data, changes):
    from hpgeff.spec import set_by_path
    d = copy.deepcopy(data)
    for key, val in changes:
        if val == "DROP_BEO":
            d["detector"]["layers"] = [l for l in d["detector"]["layers"] if l["name"] != "beo_window"]
        elif val is None:
            d["detector"].pop(key.split(".")[-1])
        else:
            set_by_path(d, key, val)
    return d


def study_model():
    """Detector-model variants vs nominal, with correlated sampling (same seeds, per-batch pairs)."""
    data = base_data()
    res = load("model_paired") or {}
    for name, ch in MODEL_VARIANTS.items():
        if name in res:
            continue
        d = _apply(data, ch)
        res[name] = dict(
            ref=run_paired(d, "reference_point_5mm", E5, 1_000_000, 20, 51, f"mp_{name}_ref"),
            glass=run_paired(d, "glass_in_pmma_cup", E5, 1_000_000, 20, 52, f"mp_{name}_glass"))
        save("model_paired", res)


STUDIES = dict(esa=study_esa, step=study_step, facets=study_facets, facets_paired=study_facets_paired,
               window=study_window, bias=study_bias, model=study_model)

if __name__ == "__main__":
    todo = sys.argv[1:] or list(STUDIES)
    for s in todo:
        t0 = time.time()
        print(f"=== {s} ===", flush=True)
        STUDIES[s]()
        print(f"=== {s} done in {time.time() - t0:.0f} s ===", flush=True)
