"""Compose docs/VALIDATION.md from validation/results (run after run_validation.py)."""
from pathlib import Path

import make_report as mr

ROOT = Path(__file__).resolve().parent.parent

TEMPLATE = """# Validation

All runs: OpenMC 0.15.2 built from source, ENDF/B-VII.1 HDF5 photon data (the same version and
library as yours), detector and samples of `examples/demo_detector2.yaml` (your detector #2,
glass in the PMMA cup on the end cap, point-source reference 5 mm above it). 10⁶–2·10⁶
histories per point unless stated; uncertainties are 1σ statistical.

DAGMC could not be built in this environment (its MOAB dependency could not be downloaded), so
the STL route was exercised up to the DAGMC file and the OpenMC model export, not transported.
The *physics* of triangulation was measured instead inside exact CSG (§3), which isolates the
faceting effect from DAGMC's numerical robustness.

## 1. Independent cross-check: Moens effective solid angle vs full Monte Carlo

Transfer factor glass/point computed two ways: full transport (OpenMC) and the Moens
effective-solid-angle integral evaluated by exact ray–cylinder intersections with no transport
and no OpenMC geometry (THEORY.md §4). They share only the configuration file and the cross
sections.

{esa}

From 46 keV to 1.3 MeV the two agree within 1.3 %; at 1.8 MeV within 2.4 %. That residual is the
accuracy of Moens' assumption that the "virtual peak-to-total" is the same for both geometries
(he quotes a few per cent): the virtual P/T values in the last column differ by the same
amount. Positions, units, materials, densities and source sampling are therefore consistent
between two unrelated calculations — and the ESA transfer alone would already meet a 5 %
target for this sample.

## 2. STEP → exact CSG (GEOUNED) vs built-in CSG

The same cup and glass sample drawn in FreeCAD, exported as STEP and converted by GEOUNED
1.6.4, against the built-in `cup` + `fill` description (2·10⁶ histories each, different seeds):

{step}

Identical within statistics. Converted cell volumes matched the CAD volumes (13.631 and
160.850 cm³) within the sampling uncertainty. Found and fixed while testing: GEOUNED writes a
`vacuum` boundary on its last surface — here the inner wall of the cup — which would have
killed every photon crossing it; the importer now resets all boundaries and renumbers surfaces
(unit test `test_geouned_xml_is_sanitised`).

## 3. Triangles ("FED"): how much faceting changes the efficiency

Cup and glass cylinders replaced by inscribed regular n-gons — exactly what an STL export does —
everything else unchanged. "density fixed": the faceted sample simply contains less glass;
"mass fixed": `mass_g` sets density = 390.864 g / faceted volume (what the STL route does).

Correlated sampling: every configuration (including the exact reference) samples the source by
the same rejection scheme with the same seed, so each history follows the same random numbers
until the changed geometry makes it diverge; the uncertainty is that of the batch-by-batch
differences (2·10⁶ histories, 20 batches). A first, uncorrelated version of this study could not
resolve effects below ~1.5 % — the paired design resolves 0.02 %.

{facets}

![faceting effect](figures/faceting_effect.png)

The effect scales with the volume deficit (∝ 1/n²): about +0.4 × deficit with the density
fixed (less glass, less self-absorption, emitters closer to the axis) and about −0.3 × deficit
with the mass fixed (slightly denser glass) at 59.5 keV, less at high energy. Even a coarse
export with 0.6 mm chords changes the efficiency by ≤ 1 %; at the default `max_deviation` of
0.05 mm (n ≥ 57 for R = 32 mm) the change is below 0.1 % — negligible against a 5 % budget.
Triangles therefore do not break an FEP efficiency through faceting; what breaks DAGMC models is
mesh topology (leaks, non-matching shared surfaces), which the STL route rejects or flags.

Runtime is the other cost: 2·10⁶ histories took 17 s with exact cylinders and ~240 s with the
same geometry as 64-facet polygons in CSG (every plane is tested at every step; the 128-facet
mass-fixed case was skipped for this reason). DAGMC's bounding-volume hierarchy avoids that
scaling; exact surfaces avoid it altogether.

The mesh-only estimator (THEORY.md §5) against the exact values of FreeCAD exports:

| mesh | estimated chord deviation | estimated volume deficit | true deficit |
|---|---:|---:|---:|
| glass sample, 0.5 mm deflection (100 triangles) | 0.233 mm | 0.979 % | 0.970 % |
| glass sample, 0.05 mm deflection (316 triangles) | 0.0247 mm | 0.103 % | 0.103 % |
| PMMA cup, 0.05 mm (640 triangles; convex + concave faces) | 0.0248 mm | 0.071 % | 0.070 % |

The 0.5 mm export is rejected by the default check (0.233 mm > 0.05 mm).

## 4. Definition of "full energy": exact deposition vs a finite window

FEP efficiency counted with a ±0.1 / 0.5 / 1 keV window, relative to exact full deposition
(±1 eV, the default):

{window}

A wider window adds photons that Compton-scattered by small angles before reaching the
crystal (THEORY.md §2). The effect is negligible at 662 keV but reaches several per cent below
~100 keV in the 50 mm glass column, and only partly cancels in the transfer (+2.9 % at
46.5 keV for ±0.5 keV). This is a real systematic for the low-energy lines, set by how the
peak areas are measured: the exact definition matches peak fits that put the small-angle
scatter tail into the background (step) function; if the fitted peak absorbs part of that
tail, set `peak_window_keV` to roughly FWHM/2 at those energies (≈ 0.3–0.4 keV for a typical
coaxial HPGe at 50–60 keV) and quote the difference between the two as an uncertainty
component.

## 5. Angular biasing (optional)

{bias}

## 6. Detector-model assumptions: absolute efficiency vs transfer factor

One parameter changed at a time from the nominal model in `demo_detector2.yaml`, with the same
random-number streams for nominal and variant (correlated sampling); uncertainties come from
batch-by-batch differences, which is why they are much smaller than the statistical
uncertainty of either efficiency.

{model}

![model sensitivity](figures/model_sensitivity.png)

Compare the glass column (what an absolute Monte Carlo result would be off by) with the
transfer column (what the transferred result is off by). For the parameters with large effects
the transfer helps by factors of 1.5 to over 100: a 0.3 mm front dead layer (the spec says
none) removes 53 % of the absolute efficiency at 46.5 keV but changes the transferred value by
2.5 %; a 1 mm error in the crystal-to-end-cap distance moves the absolute efficiency by 3–4 % and
the transferred one by 0.7–1.4 %; leaving out the BeO window changes them by 4–9 % and 1.1–2.4 %.
The transfer does not help for details that matter only to the on-axis reference source: the
core bore and the pin change the point-source efficiency at low energy by 0.3–0.7 % while
hardly touching the glass sample, so that difference enters the transferred result instead of
cancelling — small, but it is why the bore geometry still matters.

For your open spec-sheet questions the transferred efficiency of the glass sample moves by:
closed-end instead of true coaxial 0.4–2 % (largest at high energy), the 0.8 mm BeO window
1.1–2.4 % (mostly because it sets the crystal 0.8 mm further from the end cap), the contact pin
≤ 0.4 %. Adding a ±1 mm crystal-position and a +0.3 mm side-dead-layer allowance, all in
quadrature with a 2 % reference uncertainty and the Monte Carlo statistics, gives about 2.8 % at
46.5 keV, 2.7 % at 59.5 keV, 3.8 % at 662 keV and 4.2 % at 1332 keV (the peak-definition term of
§4 comes on top at the lowest energies). That fits a 5 % budget, but the BeO / crystal-position
question is the one most worth settling with the manufacturer.

## 7. Statistics and transport integrity

{stats}. Across all validation runs: 0 lost particles, no overlapping cells, and no history
deposited more than the source energy.

## 8. Failure modes deliberately provoked

| mistake | caught by |
|---|---|
| cup pushed 1 mm into the end cap | sampled overlap 3.42 cm³ (exact π·3.3²·0.1 = 3.42 cm³) and OpenMC geometry debugging: "Overlapping cells detected" |
| STL exported with 0.5 mm deflection | rejected: 0.233 mm chords, 0.98 % volume deficit |
| GEOUNED vacuum boundary on an internal surface | removed on import (unit test) |
| contact pin radius 7.775 mm > 5.35 mm bore | error suggesting the value is a diameter |
| misspelt key (`heigth`) | rejected — unknown keys are errors |
| deposits above E + window (energy bookkeeping) | counted in every run; zero in all runs here |
| two-batch trigger stop | replaced OpenMC triggers by a pilot + top-up scheme decided from the binomial estimate |
"""


def main():
    md = TEMPLATE.format(esa=mr.esa_table(), step=mr.step_table(), facets=mr.facets_table_and_fig(),
                         window=mr.window_table(), bias=mr.bias_table(), model=mr.model_table_and_fig(),
                         stats=mr.stats_summary())
    (ROOT / "docs" / "VALIDATION.md").write_text(md)
    print(md)


if __name__ == "__main__":
    main()
