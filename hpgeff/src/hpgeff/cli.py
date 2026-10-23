"""Command line interface: hpgeff {check, run, transfer, all, crosscheck, sensitivity, stlcheck}."""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

from .spec import ConfigError


def _cfg(args):
    from .config import load_config
    cfg = load_config(args.config, args.set)
    if getattr(args, "particles", None):
        cfg.sim.particles = int(args.particles)
    if getattr(args, "threads", None):
        cfg.sim.threads = int(args.threads)
    if getattr(args, "output", None):
        cfg.output_dir = Path(args.output).resolve()
    return cfg


def _setups(cfg, args):
    names = args.setup or list(cfg.setups)
    for n in names:
        if n not in cfg.setups:
            raise ConfigError(f"unknown setup '{n}' (have {list(cfg.setups)})")
    return names


def cmd_check(args):
    from .checks import CheckReport, check_nuclear_data, check_setup
    from .model import build_setup
    cfg = _cfg(args)
    rep = CheckReport()
    out = cfg.output_dir / "checks"
    out.mkdir(parents=True, exist_ok=True)
    mats_all = []
    for name in _setups(cfg, args):
        built = build_setup(cfg, name)
        mats_all += built.materials
        check_setup(built, plots_dir=out if not args.no_plots else None, rep=rep,
                    transport=getattr(args, "transport", False))
    check_nuclear_data(mats_all, rep)
    text = rep.text()
    (out / "check_report.txt").write_text(text + "\n")
    print(text)
    print(f"\n{len(rep.errors)} error(s), {len(rep.warnings)} warning(s). Report: {out / 'check_report.txt'}")
    return 1 if rep.errors else 0


def cmd_run(args):
    from .config import setup_energies
    from .report import plot_mc
    from .run import run_setup
    cfg = _cfg(args)
    energies = [float(e) for e in args.energies.split(",")] if args.energies else None
    results = {}
    for name in _setups(cfg, args):
        print(f"[{name}]")
        results[name] = run_setup(cfg, name, energies=energies or setup_energies(cfg, name))
    plot_mc(results, cfg.output_dir / "mc_efficiency.png")
    print(f"results in {cfg.output_dir}")
    return 0


def cmd_transfer(args):
    from .report import plot_transfer
    from .run import read_csv
    from .transfer import reference_from_config, transfer, write_transfer_csv
    cfg = _cfg(args)
    if not cfg.reference:
        raise ConfigError("no 'reference' section in the configuration")
    ref = cfg.reference["setup"]
    ref_rows = read_csv(cfg.output_dir / ref / "mc_efficiency.csv")
    curve = reference_from_config(cfg)
    print(f"reference fit ({curve.method}): reduced chi2 = {curve.chi2_red:.2f}")
    out = {}
    for name in _setups(cfg, args):
        if name == ref:
            continue
        p = cfg.output_dir / name / "mc_efficiency.csv"
        if not p.exists():
            print(f"  skip {name}: no Monte Carlo results ({p})")
            continue
        rows = transfer(read_csv(p), ref_rows, curve)
        write_transfer_csv(cfg.output_dir / name / "efficiency_transferred.csv", rows)
        out[name] = rows
        print(f"[{name}]  E (keV)    final efficiency        transfer factor   meas/MC(ref)")
        for r in rows:
            flag = "  (extrapolated reference)" if r["extrapolated"] else ""
            print(f"  {r['energy_keV']:9.3f}   {r['eff_final']:.5e} ± {r['u_final_rel_pct']:.2f}%   "
                  f"{r['transfer_factor']:.5f}         {r['model_ratio_ref']:.3f}{flag}")
    if out:
        plot_transfer(curve, out, ref, cfg.output_dir / "efficiency_transferred.png")
    return 0


def cmd_all(args):
    rc = cmd_check(args)
    if rc and not args.force:
        print("checks failed — fix the errors above (or use --force)")
        return rc
    cmd_run(args)
    cfg = _cfg(args)
    if cfg.reference:
        cmd_transfer(args)
    return 0


def cmd_crosscheck(args):
    """Monte Carlo vs effective-solid-angle (Moens) transfer factors."""
    from .esa import effective_solid_angle
    from .model import build_setup
    from .run import read_csv
    cfg = _cfg(args)
    names = _setups(cfg, args)
    ref = args.reference or (cfg.reference or {}).get("setup") or names[0]
    data = {}
    for name in set(names) | {ref}:
        rows = read_csv(cfg.output_dir / name / "mc_efficiency.csv")
        E = [r["energy_keV"] for r in rows]
        om, du = effective_solid_angle(build_setup(cfg, name), E, n_rays=args.rays)
        data[name] = {r["energy_keV"]: (r["fep_eff"], r["fep_unc"], o, d) for r, o, d in zip(rows, om, du)}
    out_rows = []
    print("setup               E(keV)   MC transfer   ESA transfer   MC/ESA   virtual P/T")
    for name in names:
        if name == ref:
            continue
        for E, (m, um, o, do) in sorted(data[name].items()):
            if E not in data[ref]:
                continue
            mr, umr, orr, dor = data[ref][E]
            tm, te = m / mr, o / orr
            out_rows.append(dict(setup=name, energy_keV=E, transfer_mc=tm, transfer_esa=te,
                                 ratio=tm / te, virtual_pt_setup=m / o, virtual_pt_ref=mr / orr))
            print(f"{name:18s} {E:8.2f}   {tm:.5f}      {te:.5f}       {tm / te:.3f}    {m / o:.3f} / {mr / orr:.3f}")
    p = cfg.output_dir / "crosscheck_esa.csv"
    with open(p, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out_rows[0]) if out_rows else ["setup"])
        w.writeheader()
        w.writerows(out_rows)
    print(f"written {p}")
    return 0


def cmd_sensitivity(args):
    """One-at-a-time variations of model parameters: absolute efficiency vs transfer factor."""
    import math
    from .config import load_config, setup_energies
    from .run import run_setup
    base = _cfg(args)
    if not base.reference:
        raise ConfigError("sensitivity needs a 'reference' section (the transfer is what is compared)")
    ref = base.reference["setup"]
    targets = [n for n in _setups(base, args) if n != ref]
    energies = ([float(e) for e in args.energies.split(",")] if args.energies
                else setup_energies(base, targets[0]))
    variants = [("baseline", None, None)]
    for item in args.vary:
        key, vals = item.split("=", 1)
        for v in vals.split(","):
            variants.append((f"{key}={v}", key, v))
    root = base.output_dir / "sensitivity"
    table = {}
    for label, key, val in variants:
        ov = list(args.set) + ([f"{key}={val}"] if key else [])
        cfg = load_config(args.config, ov)
        if args.particles:
            cfg.sim.particles = int(args.particles)
        if args.threads:
            cfg.sim.threads = int(args.threads)
        cfg.output_dir = root / label.replace("/", "_")
        print(f"[{label}]")
        res = {n: run_setup(cfg, n, energies=energies) for n in [ref] + targets}
        table[label] = res
    rows = []
    b = table["baseline"]
    print("\nvariant                                setup              E(keV)   d(eff)/eff   d(T)/T     (T = transfer factor)")
    for label, res in table.items():
        if label == "baseline":
            continue
        for t in targets:
            for r0, r, rr0, rr in zip(b[t], res[t], b[ref], res[ref]):
                de = r["fep_eff"] / r0["fep_eff"] - 1
                T0, T = r0["fep_eff"] / rr0["fep_eff"], r["fep_eff"] / rr["fep_eff"]
                dT = T / T0 - 1
                u = math.sqrt(sum((x["fep_unc"] / x["fep_eff"]) ** 2 for x in (r0, r, rr0, rr)))
                rows.append(dict(variant=label, setup=t, energy_keV=r["energy_keV"], rel_change_eff=de,
                                 rel_change_transfer=dT, mc_unc=u))
                print(f"{label:38s} {t:18s} {r['energy_keV']:8.2f}   {100 * de:+7.2f}%    {100 * dT:+7.2f}%   (± {100 * u:.2f}%)")
    p = root / "sensitivity.csv"
    with open(p, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]) if rows else ["variant"])
        w.writeheader()
        w.writerows(rows)
    print(f"written {p}")
    return 0


def cmd_stlcheck(args):
    import trimesh
    from .cad_stl import UNIT_TO_CM, mesh_report
    m = trimesh.load(args.file, force="mesh", process=True)
    m.merge_vertices()
    rep = mesh_report(m, args.crease)
    u = args.units
    print(json.dumps({k: v for k, v in rep.items()}, indent=2))
    print(f"\nmax chord deviation ~ {rep['max_chord_deviation']:.4g} {u}; "
          f"volume deficit vs the smooth surface ~ {100 * rep['volume_deficit_rel_est']:.3f}%")
    ok = rep["watertight"] and rep["winding_consistent"]
    print("mesh integrity: " + ("OK" if ok else "NOT usable (must be watertight with consistent winding)"))
    return 0 if ok else 1


def main(argv=None):
    p = argparse.ArgumentParser(prog="hpgeff", description=__doc__)
    sub = p.add_subparsers(dest="cmd", required=True)

    def common(sp, run=False):
        sp.add_argument("config")
        sp.add_argument("--setup", action="append", help="limit to this setup (repeatable)")
        sp.add_argument("--set", action="append", default=[], metavar="KEY=VALUE",
                        help="override a config value, e.g. --set detector.crystal.core.depth=55")
        sp.add_argument("--output", help="output directory (overrides output_dir)")
        if run:
            sp.add_argument("--energies", help="comma-separated keV list (overrides config)")
            sp.add_argument("--particles", type=float, help="histories per energy")
            sp.add_argument("--threads", type=int)

    sp = sub.add_parser("check", help="validate geometry, data and sources; make plots")
    common(sp)
    sp.add_argument("--no-plots", action="store_true")
    sp.add_argument("--transport", action="store_true",
                    help="also run a short transport with OpenMC geometry debugging (overlaps, lost particles)")
    sp.set_defaults(func=cmd_check)
    sp = sub.add_parser("run", help="Monte Carlo FEP efficiencies")
    common(sp, run=True)
    sp.set_defaults(func=cmd_run)
    sp = sub.add_parser("transfer", help="efficiency transfer from the measured reference curve")
    common(sp)
    sp.set_defaults(func=cmd_transfer)
    sp = sub.add_parser("all", help="check + run + transfer")
    common(sp, run=True)
    sp.add_argument("--no-plots", action="store_true")
    sp.add_argument("--transport", action="store_true", default=True, help=argparse.SUPPRESS)
    sp.add_argument("--force", action="store_true", help="run even if checks report errors")
    sp.set_defaults(func=cmd_all)
    sp = sub.add_parser("crosscheck", help="compare MC transfer factors with the Moens effective solid angle")
    common(sp)
    sp.add_argument("--reference", help="reference setup (default: reference.setup)")
    sp.add_argument("--rays", type=int, default=2_000_000)
    sp.set_defaults(func=cmd_crosscheck)
    sp = sub.add_parser("sensitivity", help="vary model parameters one at a time; compare absolute "
                                            "efficiency and transfer-factor changes")
    common(sp, run=True)
    sp.add_argument("--vary", action="append", required=True, metavar="KEY=V1,V2",
                    help="e.g. --vary detector.crystal.core.depth=55 --vary detector.crystal.dead_layers.side=1.0")
    sp.set_defaults(func=cmd_sensitivity)
    sp = sub.add_parser("stlcheck", help="integrity and faceting report for an STL file")
    sp.add_argument("file")
    sp.add_argument("--units", default="mm")
    sp.add_argument("--crease", type=float, default=30.0, help="dihedral angle (deg) treated as a real edge")
    sp.set_defaults(func=cmd_stlcheck)
    args = p.parse_args(argv)
    try:
        return args.func(args)
    except ConfigError as exc:
        print(f"configuration error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
