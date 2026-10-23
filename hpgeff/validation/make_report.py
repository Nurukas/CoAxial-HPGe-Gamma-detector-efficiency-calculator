"""Turn validation/results/*.json into tables (printed as markdown) and figures."""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
RES = HERE / "results"
FIG = HERE.parent / "docs" / "figures"


def load(name):
    p = RES / f"{name}.json"
    return json.loads(p.read_text()) if p.exists() else None


def _plt():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def ratio(a, b):
    r = a["fep_eff"] / b["fep_eff"]
    u = r * math.hypot(a["fep_unc"] / a["fep_eff"], b["fep_unc"] / b["fep_eff"])
    return r, u


def esa_table():
    d = load("esa")
    if not d:
        return ""
    s = ["| E (keV) | transfer, Monte Carlo | transfer, Moens ESA | MC / ESA | virtual P/T sample / reference |",
         "|---:|---:|---:|---:|---:|"]
    for r in d:
        s.append(f"| {r['E']:.2f} | {r['transfer_mc']:.4f} ± {r['transfer_mc_u']:.4f} | {r['transfer_esa']:.4f} | "
                 f"{r['ratio']:.3f} | {r['vpt_sample']:.3f} / {r['vpt_ref']:.3f} |")
    return "\n".join(s)


def step_table():
    d = load("step")
    if not d:
        return ""
    s = ["| E (keV) | built-in CSG | STEP → GEOUNED CSG | STEP / CSG |", "|---:|---:|---:|---:|"]
    for a, b in zip(d["csg"], d["step"]):
        r, u = ratio(b, a)
        s.append(f"| {a['energy_keV']:.2f} | {a['fep_eff']:.5f} | {b['fep_eff']:.5f} | {r:.4f} ± {u:.4f} |")
    return "\n".join(s)


def facets_table_and_fig():
    """Correlated-sampling faceting study (facets_paired.json)."""
    d = load("facets_paired")
    if not d or "exact" not in d:
        return ""
    ex = d["exact"]
    E = [r["energy_keV"] for r in ex]
    R = 32.0
    rows = ["| facets n | sagitta s (mm) | volume deficit | " + " | ".join(f"{e:g} keV" for e in E) + " |",
            "|---|---:|---:|" + "---:|" * len(E)]
    plt = _plt()
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True, sharex=True)
    for k, (corr, title) in enumerate((("nominal", "density fixed (2.43 g/cm3)"),
                                        ("mass fixed", "mass fixed (mass_g -> density)"))):
        tag = "nominal" if corr == "nominal" else "mass"
        for i, e in enumerate(E):
            xs, ys, us = [], [], []
            for n in (16, 32, 64, 128):
                key = f"n{n}_{tag}"
                if key not in d:
                    continue
                m, u = _paired(d[key][i], ex[i])
                xs.append(R * (1 - math.cos(math.pi / n)))
                ys.append(100 * m)
                us.append(100 * u)
            if xs:
                axes[k].errorbar(xs, ys, yerr=us, marker="o", ms=4, capsize=2, label=f"{e:g} keV")
        axes[k].axhline(0, color="#888", lw=0.8)
        axes[k].set_xscale("log")
        axes[k].set_title(title, fontsize=10)
        axes[k].set_xlabel("sagitta s of the facets (mm), R = 32 mm")
        axes[k].grid(alpha=0.3)
        axes[k].axvline(0.05, color="#c33", lw=0.8, ls=":")
        axes[k].text(0.052, 0.95, "default max_deviation", color="#c33", fontsize=7,
                     transform=axes[k].get_xaxis_transform(), va="top")
        ticks = [R * (1 - math.cos(math.pi / n)) for n in (128, 64, 32, 16)]
        axes[k].set_xticks(ticks)
        axes[k].set_xticklabels([f"{t:.2g}\n(n={n})" for t, n in zip(ticks, (128, 64, 32, 16))], fontsize=8)
        axes[k].minorticks_off()
    axes[0].set_ylabel("FEP efficiency change vs exact cylinders (%)")
    axes[1].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / "faceting_effect.png", dpi=150)
    plt.close(fig)
    for tag, lab in (("nominal", "density fixed"), ("mass", "mass fixed")):
        for n in (16, 32, 64, 128):
            key = f"n{n}_{tag}"
            if key not in d:
                continue
            sag = R * (1 - math.cos(math.pi / n))
            deficit = 1 - (n / (2 * math.pi)) * math.sin(2 * math.pi / n)
            cells = []
            for i in range(len(E)):
                m, u = _paired(d[key][i], ex[i])
                cells.append(f"{100 * m:+.2f} ± {100 * u:.2f} %")
            rows.append(f"| {n} ({lab}) | {sag:.3f} | {100 * deficit:.2f} % | " + " | ".join(cells) + " |")
    return "\n".join(rows)


def _paired(r, r0):
    a = np.asarray(r["per_batch"], float)
    b = np.asarray(r0["per_batch"], float)
    rel = (a - b) / b.mean()
    return float(rel.mean()), float(rel.std(ddof=1) / math.sqrt(len(rel)))


def window_table():
    d = load("window")
    if not d:
        return ""
    ws = d["windows"]
    s = ["| E (keV) | setup | " + " | ".join(f"±{w:g} keV vs exact" for w in ws) + " | effect on transfer (±0.5 keV) |",
         "|---:|---|" + "---:|" * len(ws) + "---:|"]
    for a, b in zip(d["ref"], d["glass"]):
        for r, lab in ((a, "point 5 mm"), (b, "glass")):
            cells = [f"{100 * (r[f'fep_w{w}'] / r['fep_eff'] - 1):+.2f} %" for w in ws]
            extra = ""
            if lab == "glass":
                t0 = b["fep_eff"] / a["fep_eff"]
                t1 = b["fep_w0.5"] / a["fep_w0.5"]
                extra = f"{100 * (t1 / t0 - 1):+.2f} %"
            s.append(f"| {r['energy_keV']:.2f} | {lab} | " + " | ".join(cells) + f" | {extra} |")
    return "\n".join(s)


def bias_table():
    d = load("bias")
    if not d:
        return ""
    s = ["| setup | E (keV) | analog | angular bias (auto) | ratio |", "|---|---:|---:|---:|---:|"]
    for setup in ("reference_point_5mm", "glass_in_pmma_cup"):
        for a, b in zip(d[f"{setup}_off"], d[f"{setup}_auto"]):
            r, u = ratio(b, a)
            s.append(f"| {setup} | {a['energy_keV']:.2f} | {a['fep_eff']:.5f} | {b['fep_eff']:.5f} "
                     f"(f = {b['bias_factor']:.3f}) | {r:.4f} ± {u:.4f} |")
    return "\n".join(s)


def model_table_and_fig():
    d = load("model_paired")
    if not d or "nominal" not in d:
        return ""
    nom = d["nominal"]
    E = [r["energy_keV"] for r in nom["glass"]]
    s = ["| variant | E (keV) | point source, absolute | glass sample, absolute | transfer factor |",
         "|---|---:|---:|---:|---:|"]
    plt = _plt()
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
    for name, v in d.items():
        if name == "nominal":
            continue
        dabs, dT, uT_all = [], [], []
        for i, e in enumerate(E):
            mr_, ur_ = _paired(v["ref"][i], nom["ref"][i])
            mg, ug = _paired(v["glass"][i], nom["glass"][i])
            a = np.asarray(v["glass"][i]["per_batch"]); a0 = np.asarray(nom["glass"][i]["per_batch"])
            b = np.asarray(v["ref"][i]["per_batch"]); b0 = np.asarray(nom["ref"][i]["per_batch"])
            dt = (a - a0) / a0.mean() - (b - b0) / b0.mean()
            mT = (1 + mg) / (1 + mr_) - 1
            uT = float(dt.std(ddof=1) / math.sqrt(len(dt)))
            dabs.append(100 * mg)
            dT.append(100 * mT)
            uT_all.append(100 * uT)
            s.append(f"| {name} | {e:.2f} | {100 * mr_:+.2f} ± {100 * ur_:.2f} % | {100 * mg:+.2f} ± {100 * ug:.2f} % | "
                     f"{100 * mT:+.2f} ± {100 * uT:.2f} % |")
        axes[0].plot(E, dabs, "o-", ms=4, label=name)
        axes[1].errorbar(E, dT, yerr=uT_all, fmt="o-", ms=4, capsize=2, label=name)
    for ax, t in zip(axes, ("absolute MC efficiency, glass sample", "transfer factor glass / point (5 mm)")):
        ax.axhline(0, color="#888", lw=0.8)
        ax.set_xscale("log")
        ax.set_title(t, fontsize=10)
        ax.set_xlabel("E (keV)")
        ax.grid(alpha=0.3)
    axes[0].set_ylabel("change vs nominal model (%)")
    axes[1].legend(fontsize=7, frameon=False)
    fig.tight_layout()
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / "model_sensitivity.png", dpi=150)
    plt.close(fig)
    return "\n".join(s)


def stats_summary():
    vals = []
    for name in ("step", "facets", "facets_paired", "window", "bias", "model_paired"):
        d = load(name)
        if not d:
            continue
        stack = [d]
        while stack:
            x = stack.pop()
            if isinstance(x, dict):
                if "fep_unc_batch" in x and x.get("fep_unc"):
                    vals.append(x["fep_unc_batch"] / x["fep_unc"])
                else:
                    stack.extend(x.values())
            elif isinstance(x, list):
                stack.extend(x)
    if not vals:
        return ""
    v = np.array(vals)
    return (f"{len(v)} runs (20 batches each): batch/binomial uncertainty ratio mean {v.mean():.3f}, "
            f"standard deviation {v.std():.3f} (expected spread for 20 batches ≈ {1 / math.sqrt(2 * 19):.3f})")


if __name__ == "__main__":
    for title, fn in (("ESA", esa_table), ("STEP", step_table), ("FACETS", facets_table_and_fig),
                      ("WINDOW", window_table), ("BIAS", bias_table), ("MODEL", model_table_and_fig),
                      ("STATS", stats_summary)):
        print(f"\n### {title}\n")
        print(fn())
