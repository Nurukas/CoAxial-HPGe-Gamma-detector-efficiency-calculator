# CAD, STL and triangles ("FED") — how they are kept from breaking the calculation

Short answer: **nothing that matters is ever triangulated unless you explicitly choose an STL.**

| Part of the model | Representation | Triangles? |
|---|---|---|
| Detector (crystal, dead layers, core, pin, holder, window, end cap) | exact CSG (planes, cylinders, torus, sphere) | never |
| Containers and samples from built-in shapes (`cup`, `marinelli`, `cylinder`, `tube`, `box`, `sphere`, `cone`, `polycone`, `polygon_prism`, boolean `csg`) | exact CSG | never |
| `shape: step` (any CAD solid made of planes / cylinders / cones / spheres / tori) | GEOUNED → exact CSG | never |
| `shape: stl` | DAGMC triangle mesh | yes — guarded |

## Why triangles hurt in neutronics, and which of those reasons apply here

1. **Lost particles.** A mesh with holes, or two volumes whose shared surface was triangulated
   twice in slightly different ways, leaves gaps/overlaps that the ray tracer cannot resolve.
   This is the real problem in neutronics CAD models. It is a *topology* problem, not a physics
   one, and it applies equally to photons.
2. **Faceting (chord) error.** Flat facets cut inside curved surfaces: volume and mass go down,
   surfaces move inwards by the sagitta. In neutronics this biases reaction rates and masses.
   For an FEP efficiency it biases the source distribution, the self-attenuation and the
   source–crystal distance.
3. **Speed and memory** (big meshes, BVH). Not an issue for a container.

## What the code does about each

**Route 1 — built-in exact shapes.** Most laboratory containers are cylinders, cups, Marinelli
beakers, bottles (`polycone`), boxes, or booleans of these. Measured dimensions go straight in.
Shapes that share a boundary (a cup and its contents via `fill:`) are built from the *same*
surface objects, so they meet exactly — the CSG equivalent of an imprinted, merged CAD model.

**Route 2 — STEP → exact CSG (GEOUNED).** When the shape is easier to draw than to describe,
draw it in any CAD program and export STEP. GEOUNED (developed for ITER/F4E neutronics
models) decomposes each solid into half-spaces of its *exact* analytic surfaces.
Guards in `hpgeff.cad_step`:
- solids with spline (NURBS) faces are refused rather than approximated;
- the vacuum boundary GEOUNED writes on its last surface is removed (left in place it would
  silently kill photons crossing an ordinary container surface — found while testing);
- surface IDs are renumbered so they cannot collide with the rest of the model;
- every converted cell's volume is sampled and must equal the CAD kernel's exact B-rep volume
  (default tolerance 0.5 %), otherwise the conversion is rejected;
- `mass_g` sets the density from the exact CAD volume.
Validated: the same cup + glass sample as built-in CSG and as STEP → GEOUNED gives the same
efficiency within statistics (VALIDATION.md, §2).

**Route 3 — STL through DAGMC, guarded** (only when a mesh is all you have — a scan, or a
3D-printed part, for which the STL *is* the printed geometry):
- non-watertight or inconsistently wound meshes are rejected;
- the chord deviation and volume deficit are estimated from the dihedral angles between
  neighbouring triangles (derivation: THEORY.md §5; reproduced the exact values of the
  FreeCAD exports to 1 %); above `checks.max_deviation` (default 0.05 mm) the mesh is rejected;
- `mass_g` → density = mass / mesh volume (mass-conserving);
- overlaps of meshes with each other, with CSG objects and with the detector are measured;
  touching meshes (separately triangulated shared boundary) trigger a warning;
- the DAGMC universe sits in a bounding cell from which the detector and CSG objects are carved
  out, so a bounding-box margin can never cover them;
- OpenMC's lost-particle limit stays strict and lost particles are counted in every run.

`hpgeff stlcheck part.stl` prints the integrity and faceting report for any mesh.

## How large is the faceting effect on an efficiency?

Measured with the real geometry of this project (glass sample in the PMMA cup, on the end cap)
by replacing the cylinders with inscribed n-gons — exactly what an STL export does, but inside
CSG so that nothing else changes — and comparing with correlated sampling (VALIDATION.md §3).

| facets per circle (R = 32 mm) | chord sagitta | volume deficit | FEP change at 59.5 keV, density fixed / mass fixed |
|---:|---:|---:|---:|
| 16 | 0.62 mm | 2.6 % | +1.0 % / −0.7 % |
| 32 | 0.15 mm | 0.64 % | +0.2 % / −0.2 % |
| 64 | 0.04 mm | 0.16 % | +0.08 % / −0.03 % |

The effect scales with the volume deficit (∝ 1/n²) and is smaller at high energy. The default
`max_deviation` of 0.05 mm (n ≥ 57 for R = 32 mm) keeps it below 0.1 %. So faceting itself does
not break an FEP efficiency at a 5 % target; topology problems (leaks, mismatched shared
surfaces) would, and those are what the STL guards reject.

## Installing the optional pieces (WSL / Ubuntu)

```bash
# GEOUNED + FreeCAD in their own environment (FreeCAD pins its own Python):
conda create -n freecad -c conda-forge freecad python=3.11
conda run -n freecad pip install geouned
# then in the config:  converter: {python: ~/miniforge3/envs/freecad/bin/python}

# STL route: a DAGMC-enabled OpenMC and the mesh tools
python -c "import openmc.lib; print(openmc.lib.dagmc_enabled)"   # must print True
pip install trimesh rtree manifold3d cad_to_dagmc
```
