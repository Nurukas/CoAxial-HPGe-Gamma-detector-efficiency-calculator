"""Run OpenMC for each energy and extract full-energy-peak efficiencies."""
from __future__ import annotations

import csv
import glob
import json
import math
import os
import re
import shutil
import subprocess
import time
from pathlib import Path

import numpy as np

from .config import Config, SimSettings, setup_energies
from .model import BuiltSetup, build_setup

FIELDS = ["energy_keV", "fep_eff", "fep_unc", "fep_unc_batch", "total_eff", "total_unc",
          "peak_to_total", "above_peak_counts", "histories", "bias_factor", "lost_particles",
          "runtime_s"]


def run_openmc(model, workdir: Path, sim: SimSettings, extra_args=()) -> tuple[Path, str]:
    """Export XML, run the OpenMC executable, return (statepoint path, log text)."""
    workdir.mkdir(parents=True, exist_ok=True)
    for f in glob.glob(str(workdir / "statepoint.*.h5")) + glob.glob(str(workdir / "particle_*.h5")):
        os.remove(f)
    model.export_to_model_xml(str(workdir / "model.xml"))
    cmd = [sim.openmc_exec]
    if sim.threads:
        cmd += ["-s", str(sim.threads)]
    cmd += list(extra_args) + ["model.xml"]
    env = dict(os.environ)
    proc = subprocess.run(cmd, cwd=workdir, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                          text=True, env=env)
    (workdir / "openmc.log").write_text(proc.stdout)
    if proc.returncode != 0:
        tail = "\n".join(proc.stdout.strip().splitlines()[-25:])
        raise RuntimeError(f"OpenMC failed in {workdir} (exit {proc.returncode}):\n{tail}")
    sps = sorted(glob.glob(str(workdir / "statepoint.*.h5")), key=os.path.getmtime)
    if not sps:
        raise RuntimeError(f"OpenMC produced no statepoint in {workdir}")
    return Path(sps[-1]), proc.stdout


_LOST = re.compile(r"(could not be located|Could not find the cell|negative distance|"
                   r"Couldn't find particle|Cannot reflect|could not find)", re.I)


def count_lost(log: str, workdir: Path | None = None) -> int:
    """Lost particles: OpenMC writes one particle_*.h5 restart file per lost particle."""
    n_files = len(glob.glob(str(workdir / "particle_*.h5"))) if workdir is not None else 0
    n_msgs = len(_LOST.findall(log))
    return max(n_files, n_msgs)


def read_result(sp_path: Path, meta: dict) -> dict:
    import warnings
    import openmc
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # ID clashes with live Python objects are harmless here
        sp = openmc.StatePoint(str(sp_path))
    with sp:
        n = int(sp.n_particles) * int(sp.n_realizations)
        fep = sp.get_tally(name="fep")
        spec = sp.get_tally(name="spectrum")
        p_raw = float(fep.mean.ravel()[0])
        sd_batch = float(fep.std_dev.ravel()[0])
        bins = spec.mean.ravel()
        runtime = float(sp.runtime.get("total", float("nan"))) if hasattr(sp, "runtime") else float("nan")
    f = meta.get("bias_factor", 1.0)
    p = p_raw
    eff = f * p
    unc = f * math.sqrt(max(p * (1 - p), 0.0) / n) if n > 0 else float("nan")
    total = f * (bins[0] + bins[1])
    tot_unc = f * math.sqrt(max((bins[0] + bins[1]) * (1 - bins[0] - bins[1]), 0.0) / n)
    return dict(energy_keV=meta["energy_keV"], fep_eff=eff, fep_unc=unc, fep_unc_batch=f * sd_batch,
                total_eff=total, total_unc=tot_unc,
                peak_to_total=(bins[1] / (bins[0] + bins[1])) if (bins[0] + bins[1]) > 0 else float("nan"),
                above_peak_counts=int(round(bins[2] * n)), histories=n, bias_factor=f,
                runtime_s=runtime)


def run_setup(cfg: Config, name: str, energies=None, out_root: Path | None = None,
              built: BuiltSetup | None = None, log=print) -> list[dict]:
    built = built or build_setup(cfg, name)
    energies = energies or setup_energies(cfg, name)
    out = (out_root or cfg.output_dir) / name
    out.mkdir(parents=True, exist_ok=True)
    results = []
    for w in built.warnings:
        log(f"  WARNING: {w}")
    for E in energies:
        t0 = time.time()
        model, meta = built.make_model(E)
        wd = out / f"run_{E:.3f}keV"
        sp, text = run_openmc(model, wd, cfg.sim)
        res = read_result(sp, meta)
        res["lost_particles"] = count_lost(text, wd)
        target = cfg.sim.rel_unc_target
        if target and res["fep_eff"] > 0:
            # histories needed from the pilot estimate (binomial), then one top-up run with a
            # new seed. Deciding from p (not from the noisy observed error) avoids the bias of
            # stopping whenever a fluctuation happens to look good.
            p = res["fep_eff"] / res["bias_factor"]
            n_need = int(math.ceil(1.1 * (1 - p) / (p * target ** 2)))
            cap = int(cfg.sim.max_batches * cfg.sim.particles / cfg.sim.batches)
            n_add = min(n_need, cap) - res["histories"]
            if n_add > 0:
                import copy
                sim2 = copy.copy(cfg.sim)
                sim2.seed = cfg.sim.seed + 7919
                sim2.particles = n_add
                sim2.batches = max(cfg.sim.batches, 1)
                model2, meta2 = built.make_model(E, sim2)
                sp2, text2 = run_openmc(model2, wd / "topup", sim2)
                res = _combine(res, read_result(sp2, meta2))
                res["lost_particles"] += count_lost(text2, wd / "topup")
        if res["lost_particles"]:
            log(f"  WARNING: {res['lost_particles']} lost particle(s) at {E} keV — see {wd}/openmc.log")
        if not math.isfinite(res["runtime_s"]):
            res["runtime_s"] = time.time() - t0
        if res["above_peak_counts"] > 0:
            log(f"  WARNING: {res['above_peak_counts']} histories deposited more than E + window "
                f"at {E} keV — energy bookkeeping problem, check the model")
        results.append(res)
        rel = res["fep_unc"] / res["fep_eff"] if res["fep_eff"] > 0 else float("nan")
        if cfg.sim.rel_unc_target and rel > 1.25 * cfg.sim.rel_unc_target:
            log(f"  WARNING: binomial uncertainty {100 * rel:.2f}% above the target "
                f"{100 * cfg.sim.rel_unc_target:.2f}% (capped by simulation.max_batches)")
        log(f"  {name:>16s}  E = {E:9.3f} keV   FEP eff = {res['fep_eff']:.5e} ± {100 * rel:.2f}%"
            f"   (N = {res['histories']:.2e}, {time.time() - t0:.1f} s)")
        if not cfg.sim.keep_run_files:
            shutil.rmtree(wd, ignore_errors=True)
    write_csv(out / "mc_efficiency.csv", results, FIELDS)
    with open(out / "mc_efficiency.json", "w") as fh:
        json.dump(dict(setup=name, meta=setup_summary(built), results=results), fh, indent=2,
                  default=_jsonable)
    return results


def _combine(a: dict, b: dict) -> dict:
    """Pool two independent runs of the same model (counts add)."""
    na, nb = a["histories"], b["histories"]
    n = na + nb
    f = a["bias_factor"]
    out = dict(a)
    for key in ("fep_eff", "total_eff"):
        out[key] = (a[key] * na + b[key] * nb) / n
    p, t = out["fep_eff"] / f, out["total_eff"] / f
    out["fep_unc"] = f * math.sqrt(max(p * (1 - p), 0.0) / n)
    out["total_unc"] = f * math.sqrt(max(t * (1 - t), 0.0) / n)
    out["fep_unc_batch"] = math.sqrt((a["fep_unc_batch"] * na) ** 2 + (b["fep_unc_batch"] * nb) ** 2) / n
    out["peak_to_total"] = p / t if t > 0 else float("nan")
    out["above_peak_counts"] = a["above_peak_counts"] + b["above_peak_counts"]
    out["histories"] = n
    out["runtime_s"] = a["runtime_s"] + b["runtime_s"]
    return out


def write_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k) for k in fields})


def read_csv(path: Path) -> list[dict]:
    with open(path) as fh:
        rows = list(csv.DictReader(fh))
    out = []
    for r in rows:
        out.append({k: (float(v) if _isnum(v) else v) for k, v in r.items()})
    return out


def _isnum(v):
    try:
        float(v)
        return True
    except (TypeError, ValueError):
        return False


def _jsonable(o):
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, Path):
        return str(o)
    return str(o)


def setup_summary(built: BuiltSetup) -> dict:
    d = built.detector.info
    objs = []
    for o in built.objects:
        objs.append(dict(name=o.name, kind=o.kind, material=o.material_name,
                         density_g_cm3=o.density, volume_cm3=o.volume, volume_source=o.volume_source,
                         mass_g=(o.density * o.volume) if (o.density and o.volume) else None))
    return dict(detector=d, objects=objs, ambient=built.ambient_name,
                warnings=built.warnings)
