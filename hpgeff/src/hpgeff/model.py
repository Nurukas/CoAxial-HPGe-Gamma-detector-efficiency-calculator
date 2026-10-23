"""Assemble OpenMC models: detector + objects + ambient + source + tallies."""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from .config import Config, SimSettings
from .detector import Detector, build_detector
from .materials import MaterialLibrary
from .regions import mc_volume
from .shapes import Shape, build_csg_object
from .spec import ConfigError, Spec


@dataclass
class BuiltObject:
    name: str
    kind: str                      # 'csg' | 'step' | 'stl'
    material_name: str | None
    material: object               # openmc.Material or None (vacuum)
    cell: object | None            # openmc.Cell for CSG/STEP objects
    region: object | None          # world-frame CSG region (CSG/STEP), None for STL
    bbox: tuple
    volume: float | None           # cm3
    volume_unc: float = 0.0
    volume_source: str = ""        # 'analytic' | 'cad' | 'mesh' | 'monte-carlo'
    mass_g: float | None = None
    density: float | None = None
    shape: Shape | None = None     # CSG shape (for exact source sampling, ESA)
    mesh: object | None = None     # trimesh mesh in cm (STL)
    notes: list = field(default_factory=list)


@dataclass
class BuiltSetup:
    name: str
    cfg: Config
    lib: MaterialLibrary
    detector: Detector
    objects: list
    ambient_material: object
    ambient_name: str
    world_radius: float
    cells: list
    materials: list
    dagmc: object | None = None            # (DAGMCUniverse, bounding cell) if any STL
    source_spec: dict = field(default_factory=dict)
    warnings: list = field(default_factory=list)

    def obj(self, name: str) -> BuiltObject:
        for o in self.objects:
            if o.name == name:
                return o
        raise ConfigError(f"setup '{self.name}': no object named '{name}' "
                          f"(objects: {[o.name for o in self.objects]})")

    # ------------------------------------------------------------------
    def make_model(self, energy_keV: float, sim: SimSettings | None = None):
        """Full OpenMC model for one photon energy. Returns (model, meta)."""
        import openmc
        from .sources import build_sources, angular_bias_cone
        sim = sim or self.cfg.sim
        E = float(energy_keV) * 1e3  # eV
        w = sim.peak_window_keV * 1e3
        if w >= E:
            raise ConfigError("peak window wider than the photon energy")

        bias = None
        if sim.angular_bias == "auto":
            bias = angular_bias_cone(self)
        sources, src_meta = build_sources(self.source_spec, self, E, bias)

        settings = openmc.Settings()
        settings.run_mode = "fixed source"
        settings.photon_transport = True
        settings.electron_treatment = sim.electron_treatment
        settings.cutoff = {"energy_photon": sim.photon_cutoff_keV * 1e3}
        settings.survival_biasing = False  # pulse-height scores ignore particle weight
        settings.seed = sim.seed
        settings.source = sources
        settings.max_lost_particles = sim.max_lost_particles
        settings.rel_max_lost_particles = sim.rel_max_lost_particles
        settings.output = {"tallies": False, "summary": False}
        per_batch = max(1, int(round(sim.particles / sim.batches)))
        settings.particles = per_batch
        settings.batches = sim.batches

        active = self.detector.active_cell
        fep = openmc.Tally(name="fep")
        fep.filters = [openmc.CellFilter([active]), openmc.EnergyFilter([E - w, E + w])]
        fep.scores = ["pulse-height"]
        spec = openmc.Tally(name="spectrum")
        spec.filters = [openmc.CellFilter([active]),
                        openmc.EnergyFilter([1.0, E - w, E + w, E + w + 10e3])]
        spec.scores = ["pulse-height"]

        geom = openmc.Geometry(self.cells)
        model = openmc.Model(geometry=geom, materials=openmc.Materials(self.materials),
                             settings=settings, tallies=openmc.Tallies([fep, spec]))
        meta = dict(energy_keV=float(energy_keV), window_keV=sim.peak_window_keV,
                    bias_factor=bias[1] if bias else 1.0,
                    bias_half_angle_deg=math.degrees(bias[0]) if bias else 180.0, **src_meta)
        return model, meta


# ----------------------------------------------------------------------------
# Building a setup
# ----------------------------------------------------------------------------
def build_setup(cfg: Config, name: str) -> BuiltSetup:
    import openmc
    if name not in cfg.setups:
        raise ConfigError(f"unknown setup '{name}' (have: {list(cfg.setups)})")
    lib = MaterialLibrary(cfg.materials)
    det = build_detector(Spec(cfg.detector, "detector"), cfg.scale, lib)
    st = cfg.setups[name]
    objects: list[BuiltObject] = []
    warnings: list[str] = []
    dagmc_items = []
    for i, raw in enumerate(st.get("objects") or []):
        s = Spec(raw, f"setups.{name}.objects[{i}]")
        oname = str(s.get("name", f"object{i}"))
        kind = s.get("shape", None)
        if kind == "step":
            from .cad_step import build_step_objects
            objs = build_step_objects(oname, s, cfg, lib)
            objects += objs
        elif kind == "stl":
            from .cad_stl import load_stl_object
            obj = load_stl_object(oname, s, cfg, lib)
            objects.append(obj)
            dagmc_items.append(obj)
        else:
            objects += _build_csg_objects(oname, s, cfg, lib)
    names = [o.name for o in objects]
    dup = {n for n in names if names.count(n) > 1}
    if dup:
        raise ConfigError(f"setup '{name}': duplicate object names {sorted(dup)}")

    # ---------------- world & ambient ------------------------------------------
    corners = [det.envelope_bbox[0], det.envelope_bbox[1]]
    for o in objects:
        corners += [np.asarray(o.bbox[0]), np.asarray(o.bbox[1])]
    rmax = max(float(np.linalg.norm(c)) for c in corners)
    world_r = 1.5 * rmax + 10.0
    world = openmc.Sphere(r=world_r, boundary_type="vacuum", name="world")
    amb_name = cfg.ambient
    amb_mat = lib.make_openmc(amb_name, "ambient")
    cells = list(det.cells)
    materials = [m for m in det.materials if m is not None]
    amb_region = -world & ~det.envelope
    for o in objects:
        if o.cell is not None:
            cells.append(o.cell)
            amb_region = amb_region & ~o.region
        if o.material is not None:
            materials.append(o.material)
    dagmc = None
    if dagmc_items:
        from .cad_stl import build_dagmc_universe
        dag_univ, bound_region = build_dagmc_universe(dagmc_items, cfg, amb_name, setup_name=name,
                                                      warnings=warnings)
        # the DAGMC universe is visible only inside this cell: carve out the detector and the
        # CSG objects so the box margin can never cover them (overlaps are checked separately)
        bound_region = bound_region & ~det.envelope
        for o in objects:
            if o.cell is not None:
                bound_region = bound_region & ~o.region
        bcell = openmc.Cell(name="dagmc bounding cell", fill=dag_univ, region=bound_region)
        cells.append(bcell)
        amb_region = amb_region & ~bound_region
        dagmc = (dag_univ, bcell)
    if amb_mat is not None:
        materials.append(amb_mat)
    cells.append(openmc.Cell(name="ambient", fill=amb_mat, region=amb_region))

    for o in objects:
        warnings += [n for n in o.notes if "differs" in n or "WARNING" in n]
    built = BuiltSetup(name=name, cfg=cfg, lib=lib, detector=det, objects=objects,
                       ambient_material=amb_mat, ambient_name=amb_name, world_radius=world_r,
                       cells=cells, materials=materials, dagmc=dagmc,
                       source_spec=st.get("source"), warnings=warnings)
    return built


def _density_from_mass(name, mass, volume, vol_unc, nominal, warnings_list):
    rho = mass / volume
    if nominal and abs(rho / nominal - 1) > 0.25:
        warnings_list.append(
            f"object '{name}': density from mass/volume = {rho:.4g} g/cm3 differs from the "
            f"material's nominal {nominal:.4g} g/cm3 by {100 * (rho / nominal - 1):+.0f}% "
            f"(check units/mass) — using the mass-derived value")
    return rho


def _build_csg_objects(oname: str, s: Spec, cfg: Config, lib: MaterialLibrary) -> list[BuiltObject]:
    import openmc
    mat_name = s.get("material", None)
    mass = s.get("mass_g", None)
    dens = s.get("density", None)
    shapes = build_csg_object(oname, s, cfg.scale)
    s.done()
    out = []
    for sh in shapes:
        if sh is shapes[0]:
            m_name, m_mass, m_dens = mat_name, mass, dens
        else:  # cup/marinelli contents or separate lid
            fd = sh.params.get("fill") or sh.params
            m_name, m_mass, m_dens = fd.get("material"), fd.get("mass_g"), fd.get("density")
        if m_name is None:
            raise ConfigError(f"object '{sh.name}': 'material' is required")
        vol, vunc, vsrc = sh.volume, 0.0, "analytic"
        if vol is None:
            ll, ur = sh.bbox
            vol, vunc = mc_volume(sh.region, ll, ur, n=2_000_000)
            vsrc = "monte-carlo"
        md = lib.resolve(m_name)
        notes = []
        rho = None
        if m_mass is not None:
            if md is None:
                raise ConfigError(f"object '{sh.name}': mass given for a vacuum object")
            rho = _density_from_mass(sh.name, float(m_mass), vol, vunc, md.density, notes)
        elif m_dens is not None:
            rho = float(m_dens)
        mat = lib.make_openmc(m_name, sh.name, density=rho)
        cell = openmc.Cell(name=f"object {sh.name}", fill=mat, region=sh.region)
        out.append(BuiltObject(
            name=sh.name, kind="csg", material_name=m_name, material=mat, cell=cell,
            region=sh.region, bbox=sh.bbox, volume=vol, volume_unc=vunc, volume_source=vsrc,
            mass_g=None if m_mass is None else float(m_mass),
            density=(mat.density if mat is not None else None), shape=sh, notes=notes))
    return out
