"""Plots of efficiency curves and transfer diagnostics."""
from __future__ import annotations

from pathlib import Path

import numpy as np


def _plt():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


COLORS = ["#2a6fdb", "#d1495b", "#2e9e6b", "#e09f3e", "#7b5ea7", "#3a3a3a"]


def plot_mc(results: dict, path: Path):
    """results: {setup: [rows]} -> log-log FEP efficiency vs energy."""
    plt = _plt()
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for i, (name, rows) in enumerate(results.items()):
        E = np.array([r["energy_keV"] for r in rows], float)
        e = np.array([r["fep_eff"] for r in rows], float)
        u = np.array([r["fep_unc"] for r in rows], float)
        ax.errorbar(E, e, yerr=u, fmt="o-", ms=4, lw=1.2, capsize=2, color=COLORS[i % 6], label=name)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Photon energy (keV)")
    ax.set_ylabel("Full-energy peak efficiency (Monte Carlo)")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_transfer(curve, transfers: dict, ref_name: str, path: Path):
    plt = _plt()
    fig, (ax, ax2) = plt.subplots(2, 1, figsize=(7, 7), sharex=True,
                                  gridspec_kw={"height_ratios": [2.2, 1]})
    Ef = np.geomspace(curve.E.min() * 0.9, curve.E.max() * 1.1, 300)
    v, u = curve(Ef)
    ax.fill_between(Ef, v - u, v + u, color="#999999", alpha=0.25, lw=0)
    ax.plot(Ef, v, color="#555555", lw=1, label=f"{ref_name}: reference fit")
    ax.errorbar(curve.E, curve.eff, yerr=curve.u, fmt="s", ms=4, color="#222222",
                label=f"{ref_name}: measured")
    for i, (name, rows) in enumerate(transfers.items()):
        E = np.array([r["energy_keV"] for r in rows])
        f = np.array([r["eff_final"] for r in rows])
        uf = np.array([r["u_final"] for r in rows])
        ax.errorbar(E, f, yerr=uf, fmt="o-", ms=4, lw=1.2, capsize=2, color=COLORS[i % 6],
                    label=f"{name}: transferred")
        mr = np.array([r["model_ratio_ref"] for r in rows])
        umr = np.array([r["u_model_ratio_ref"] for r in rows])
        if i == 0:
            ax2.errorbar(E, mr, yerr=umr, fmt="o-", ms=4, color="#222222", lw=1)
    ax2.axhline(1.0, color="#888888", lw=0.8, ls="--")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylabel("FEP efficiency")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(frameon=False, fontsize=8)
    ax2.set_xlabel("Photon energy (keV)")
    ax2.set_ylabel("measured / MC\n(reference geometry)")
    ax2.grid(True, which="both", alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
