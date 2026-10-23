"""Efficiency transfer from a measured reference curve.

    eps_x(E) = eps_ref,meas(E) * eps_x,MC(E) / eps_ref,MC(E)

The ratio of two Monte Carlo efficiencies computed with the *same* detector
model is much less sensitive to errors in that model (dead layers, crystal
position, ...) than either efficiency alone — the idea of Moens et al. (1981),
eq. (2), with full Monte Carlo replacing the effective solid angle.
"""
from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np

from .spec import ConfigError, Spec


# ----------------------------------------------------------------------------
# Reference data
# ----------------------------------------------------------------------------
def read_reference(path: Path):
    """CSV/whitespace table: energy_keV, efficiency, uncertainty.

    The uncertainty column is absolute unless its header contains '%' (percent)
    or 'rel' (relative fraction).
    """
    rows, header = [], None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = [p for p in line.replace(";", ",").replace("\t", ",").split(",") if p.strip()]
            if len(parts) == 1:
                parts = line.split()
            try:
                vals = [float(p) for p in parts]
            except ValueError:
                header = [p.strip().lower() for p in parts]
                continue
            rows.append(vals)
    if not rows:
        raise ConfigError(f"reference file {path} has no numeric rows")
    a = np.asarray(rows, float)
    if a.shape[1] < 2:
        raise ConfigError("reference file needs at least energy and efficiency columns")
    E, eff = a[:, 0], a[:, 1]
    if a.shape[1] >= 3:
        u = a[:, 2]
        h = header[2] if header and len(header) > 2 else ""
        if "%" in h:
            u = u / 100.0 * eff
        elif "rel" in h:
            u = u * eff
    else:
        u = 0.02 * eff  # unknown: assume 2 %
    order = np.argsort(E)
    E, eff, u = E[order], eff[order], u[order]
    if np.any(eff <= 0) or np.any(u <= 0):
        raise ConfigError("reference efficiencies and uncertainties must be positive")
    return E, eff, u


class ReferenceCurve:
    """Smooth model of the measured reference efficiency with uncertainty."""

    def __init__(self, E, eff, u, method="loglog_poly", degree=4, use_measured_at_match=False,
                 match_tol_keV=0.05):
        self.E, self.eff, self.u = np.asarray(E), np.asarray(eff), np.asarray(u)
        self.method = method
        self.use_measured = use_measured_at_match
        self.tol = match_tol_keV
        self.Ec = float(np.exp(np.mean(np.log(self.E))))
        if method == "loglog_poly":
            n = len(self.E)
            if degree > n - 2:
                raise ConfigError(f"reference fit: degree {degree} needs at least {degree + 2} points "
                                  f"(have {n})")
            self.degree = degree
            x = np.log(self.E / self.Ec)
            y = np.log(self.eff)
            sy = self.u / self.eff
            A = np.vander(x, degree + 1, increasing=True)
            W = 1.0 / sy ** 2
            AtW = A.T * W
            cov = np.linalg.inv(AtW @ A)
            coef = cov @ (AtW @ y)
            res = y - A @ coef
            dof = max(n - degree - 1, 1)
            self.chi2_red = float(np.sum(res ** 2 * W) / dof)
            self.cov = cov * max(1.0, self.chi2_red)  # inflate if scatter exceeds stated errors
            self.coef = coef
            self.residuals_rel = np.exp(res) - 1.0
        elif method == "loglog_linear":
            self.chi2_red = float("nan")
        else:
            raise ConfigError("reference.fit.method must be loglog_poly or loglog_linear")

    def __call__(self, E):
        E = np.atleast_1d(np.asarray(E, float))
        val, unc = np.empty_like(E), np.empty_like(E)
        for i, e in enumerate(E):
            j = int(np.argmin(np.abs(self.E - e)))
            if self.use_measured and abs(self.E[j] - e) <= self.tol:
                val[i], unc[i] = self.eff[j], self.u[j]
                continue
            if self.method == "loglog_poly":
                xv = np.vander([math.log(e / self.Ec)], self.degree + 1, increasing=True)[0]
                y = float(xv @ self.coef)
                sy = math.sqrt(max(float(xv @ self.cov @ xv), 0.0))
                val[i], unc[i] = math.exp(y), math.exp(y) * sy
            else:
                lx, ly = np.log(self.E), np.log(self.eff)
                val[i] = math.exp(np.interp(math.log(e), lx, ly))
                unc[i] = val[i] * np.interp(math.log(e), lx, self.u / self.eff)
        return val, unc

    def extrapolated(self, E):
        E = np.atleast_1d(E)
        return (E < self.E.min() * 0.999) | (E > self.E.max() * 1.001)


def reference_from_config(cfg) -> ReferenceCurve:
    rs = Spec(cfg.reference, "reference")
    path = cfg.resolve_path(rs.get("data"))
    fit = Spec(rs.get("fit", {}) or {}, "reference.fit")
    method = fit.choice("method", ("loglog_poly", "loglog_linear"), "loglog_poly")
    deg = int(fit.num("degree", 4, positive=True))
    use_m = bool(fit.get("use_measured_at_match", False))
    fit.done()
    E, eff, u = read_reference(path)
    return ReferenceCurve(E, eff, u, method=method, degree=deg, use_measured_at_match=use_m)


# ----------------------------------------------------------------------------
# Transfer
# ----------------------------------------------------------------------------
TRANSFER_FIELDS = ["energy_keV", "eff_final", "u_final", "u_final_rel_pct", "transfer_factor",
                   "u_transfer", "eff_mc", "u_mc", "eff_ref_mc", "u_ref_mc", "eff_ref_meas",
                   "u_ref_meas", "model_ratio_ref", "u_model_ratio_ref", "extrapolated"]


def _match(rows, E, tol=1e-3):
    for r in rows:
        if abs(float(r["energy_keV"]) - E) <= tol:
            return r
    return None


def transfer(sample_rows, ref_rows, curve: ReferenceCurve) -> list[dict]:
    out = []
    for r in sample_rows:
        E = float(r["energy_keV"])
        rr = _match(ref_rows, E)
        if rr is None:
            raise ConfigError(f"no reference-setup Monte Carlo result at {E} keV — run the reference "
                              f"setup at the same energies")
        em, um = float(r["fep_eff"]), float(r["fep_unc"])
        er, ur = float(rr["fep_eff"]), float(rr["fep_unc"])
        if em <= 0 or er <= 0:
            raise ConfigError(f"zero Monte Carlo efficiency at {E} keV — increase histories")
        T = em / er
        uT = T * math.hypot(um / em, ur / er)
        ex, uex = curve(E)
        ex, uex = float(ex[0]), float(uex[0])
        final = ex * T
        ufin = final * math.hypot(uex / ex, uT / T)
        mr = ex / er
        umr = mr * math.hypot(uex / ex, ur / er)
        out.append(dict(energy_keV=E, eff_final=final, u_final=ufin, u_final_rel_pct=100 * ufin / final,
                        transfer_factor=T, u_transfer=uT, eff_mc=em, u_mc=um, eff_ref_mc=er,
                        u_ref_mc=ur, eff_ref_meas=ex, u_ref_meas=uex, model_ratio_ref=mr,
                        u_model_ratio_ref=umr, extrapolated=bool(curve.extrapolated(E)[0])))
    return out


def write_transfer_csv(path: Path, rows):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=TRANSFER_FIELDS)
        w.writeheader()
        for r in rows:
            w.writerow(r)
