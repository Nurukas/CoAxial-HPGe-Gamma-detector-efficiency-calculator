# hpgeff — full-energy-peak efficiency of coaxial HPGe detectors with OpenMC

`hpgeff` computes full-energy-peak (FEP) efficiencies of a coaxial HPGe detector for arbitrary
sample geometries and materials, with OpenMC photon transport (pulse-height tally), and
optionally transfers them to your measured reference curve.

- **Geometry without triangles by default.** Detector, containers and samples are exact CSG
  (quadric surfaces). CAD comes in as STEP and is converted to *exact* CSG by GEOUNED. STL
  meshes are supported through DAGMC only as a guarded fallback. See
  [docs/CAD_AND_TRIANGLES.md](docs/CAD_AND_TRIANGLES.md) — the answer to "will triangles (FED)
  break the calculation?" with measured numbers.
- **Any material**: formulas, element or compound mass fractions, mixtures of materials;
  density from a weighed mass.
- **Source models**: point, disk, ring, uniform volume, surface contamination (chosen faces),
  depth profiles (sedimentation), user point lists, mixtures.
- **Efficiency transfer** (Moens 1981 idea with full Monte Carlo): final efficiency =
  measured reference × MC(sample)/MC(reference), with uncertainty propagation.
- **Checks before you trust a number**: volumes vs analytic/CAD values, overlaps (sampling and
  OpenMC geometry debugging), lost particles, energy bookkeeping, nuclear data, slice plots,
  and an independent effective-solid-angle cross-check.

Derivations: [docs/THEORY.md](docs/THEORY.md). Validation: [docs/VALIDATION.md](docs/VALIDATION.md).

## Install (WSL / Ubuntu, in your existing OpenMC 0.15.x environment)

```bash
cd hpgeff
pip install -e .                 # numpy, scipy, h5py, pyyaml, matplotlib
pip install -e ".[stl]"          # optional: STL route (trimesh, manifold3d, cad_to_dagmc)
export OPENMC_CROSS_SECTIONS=/path/to/endfb-vii.1-hdf5/cross_sections.xml
```

OpenMC needs **neutron data for every nuclide even in photon-only runs** (your full
ENDF/B-VII.1 library has them) plus the photon files. `hpgeff check` verifies both.

## Quick start (your detector #2, glass sample in the PMMA cup)

```bash
cd examples
hpgeff check demo_detector2.yaml --transport   # geometry, volumes, overlaps, plots, short debug run
hpgeff run demo_detector2.yaml                 # MC efficiencies, both setups, all energies
hpgeff transfer demo_detector2.yaml            # transfer to the measured point-source curve
hpgeff crosscheck demo_detector2.yaml          # independent Moens effective-solid-angle check
```

Outputs (in `output_dir`):

| file | content |
|---|---|
| `checks/check_report.txt`, `checks/*_xz.png` | pre-flight report and OpenMC slice plots |
| `<setup>/mc_efficiency.csv` | MC FEP efficiency, binomial and batch uncertainties, total efficiency, P/T, lost particles |
| `<setup>/efficiency_transferred.csv` | **final efficiency** with total uncertainty, transfer factor, reference value used, measured/MC ratio of the reference |
| `mc_efficiency.png`, `efficiency_transferred.png` | plots |

`reference_curve_SYNTHETIC.csv` in `examples/` is made-up test data so the pipeline runs end
to end — replace it with the experimentalist's curve (`energy_keV, efficiency, uncertainty`;
an uncertainty column header containing `%` or `rel` is read as relative).

Change any value from the command line, e.g. to test the closed-end reading of the drawing:

```bash
hpgeff run demo_detector2.yaml --set detector.crystal.core.depth=55 --output results_closed_end
```

## Configuration in one page

```yaml
units: mm                     # all lengths; z = 0 is the outer face of the end cap, detector looks up (+z)
materials:                    # densities g/cm3; built-ins: air, germanium, aluminium, beryllium(_oxide),
  glass: {density: 2.43, compounds: {SiO2: 75, Na2O: 11.7, CaO: 5}}   # copper, brass, lead, pmma, ...
  steel: {density: 7.9, elements: {Fe: 0.7, Cr: 0.19, Ni: 0.11}}     # (mass fractions; fractions: atom)
  pmma2: {base: pmma, density: 1.18}
  slurry: {density: 1.3, mixture: {soil: 0.7, water: 0.3}}
detector:
  crystal: {radius, length, edge_rounding, material,
            core: {radius, depth, bottom: flat|round},        # depth >= length: true coaxial
            dead_layers: {front, side, back, core}}
  contact_pin: {radius, length, material}                      # optional, inside the core
  layers:                                                      # concentric cups, crystal -> out
    - {name, material, top, side, gap_top, gap_side}           # side: 0 -> disk; top: 0 -> tube
  gap_material: vacuum
setups:
  <name>:
    objects:
      - {name, shape: cylinder|tube|cup|marinelli|box|sphere|cone|polycone|polygon_prism|csg|step|stl,
         material, mass_g | density, <shape parameters>, rotate: [ax, ay, az]}
    source: {type: point|disk|ring|volume|surface|depth_profile|points|mixture, ...}
energies_keV: [...]
simulation: {particles, batches, rel_unc_target, peak_window_keV, threads, seed, angular_bias: off|auto}
reference: {setup, data, fit: {method: loglog_poly|loglog_linear, degree}}
```

Shape parameters (mm): `cylinder` radius, height, z, center · `tube` inner_radius,
outer_radius, height · `cup` inner_radius, wall, bottom, inner_height, lid {thickness},
**fill {name, material, height, mass_g}** (creates the sample sharing the cup's exact inner
surfaces) · `marinelli` inner_radius, wall, bottom, inner_height, well_radius, well_depth,
seat_z, fill · `polycone` planes [[z, r_in, r_out], ...] · `csg` op union|intersection|difference,
shapes [...] · `step` file, solids [{name, material, mass_g}], placement, converter {python} ·
`stl` file, units, material, mass_g, placement, checks {max_deviation}. `facets: n` on
cylinder/tube/cup emulates an STL export inside exact CSG. Typos are rejected: unknown keys
are errors, so a misspelt dimension cannot silently change a result.

[examples/shapes_gallery.yaml](examples/shapes_gallery.yaml) shows every shape and source
model; `cad_step_example.yaml` and `cad_stl_example.yaml` show the two CAD routes.

## Your detector: items still to confirm

`demo_detector2.yaml` marks four assumptions (`ASSUMPTION` comments): true vs closed-end core,
the position of the 0.8 mm BeO window, the contact pin (spec "radius 7.775 mm" > 5.35 mm bore
— read as a diameter; length assumed), and the batch-vs-glass composition. VALIDATION.md §6
shows how much each moves the absolute efficiency and how much less it moves the transferred
one; `hpgeff sensitivity` repeats that for any parameter:

```bash
hpgeff sensitivity demo_detector2.yaml --vary detector.crystal.core.depth=55 \
       --vary detector.layers.2.gap_top=6.0 --energies 59.54,661.66
```

## Accuracy in one paragraph

Against an independent effective-solid-angle calculation the Monte Carlo transfer factors
agree within 1.3 % from 46 keV to 1.3 MeV. Transferring from the measured reference cuts the
large detector-model errors (dead layers, crystal position, window) by factors of 1.5 to over
100, so the open spec-sheet questions cost about 0.3–2.4 % each on the transferred glass
efficiency and the estimated total stays under 5 %.
The one systematic that does not shrink much is the definition of "full energy" below
~100 keV for thick samples (small-angle Compton scatter, up to ~3 % at 46.5 keV) — match
`peak_window_keV` to how the peak areas are fitted. Details: docs/VALIDATION.md.

## Python API

```python
from hpgeff.config import load_config
from hpgeff.model import build_setup
from hpgeff.run import run_setup
cfg = load_config("demo_detector2.yaml", ["simulation.particles=2e6"])
rows = run_setup(cfg, "glass_in_pmma_cup", energies=[661.66])
built = build_setup(cfg, "glass_in_pmma_cup")
model, meta = built.make_model(661.66)      # a plain openmc.Model, to inspect or extend
```

## Tests

`pytest -q` (unit tests: shapes, detector partition, materials, ray/tube intersections,
reference fit, STL faceting estimator, attenuation coefficients). `validation/run_validation.py`
reproduces VALIDATION.md.
