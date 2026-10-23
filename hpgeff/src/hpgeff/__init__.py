"""hpgeff — full-energy-peak (FEP) efficiency of coaxial HPGe detectors with OpenMC.

Geometry routes, in order of preference:

1. Exact CSG (detector, containers, samples built from quadric surfaces) — no triangles.
2. STEP CAD -> exact CSG through GEOUNED (planes, cylinders, cones, spheres, tori) — no triangles.
3. STL triangle meshes through DAGMC — guarded: watertightness, faceting deviation,
   volume/mass, overlap and lost-particle checks; needs an OpenMC build with DAGMC.

The efficiency is obtained from an analog OpenMC pulse-height tally on a single
active-germanium cell; optionally transferred to a measured reference curve
(efficiency transfer), which cancels most detector-model errors.
"""

__version__ = "0.1.0"


def load_config(path, overrides=None):
    """Load and validate a YAML configuration (see :func:`hpgeff.config.load_config`)."""
    from .config import load_config as _lc
    return _lc(path, overrides)
