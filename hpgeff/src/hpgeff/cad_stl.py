"""STL triangle meshes through DAGMC — the guarded fallback route.

Prefer exact CSG (built-in shapes) or STEP -> CSG (GEOUNED). Use STL only when a
mesh is the only description of an object (3D scans, 3D-printed parts — for a
printed part the STL *is* the geometry that was printed).

What can go wrong with triangles, and what is done about it:

1. Leaky / non-manifold meshes -> lost particles. Every mesh must be watertight
   with consistent winding and positive volume, or it is rejected. OpenMC's
   lost-particle limit stays at its strict default and lost particles are
   counted in every run.
2. Faceting (chord) error: a curved surface replaced by flat triangles moves
   the surface inwards by the sagitta s and removes volume. The sagitta and the
   volume deficit are estimated from the dihedral angles between neighbouring
   triangles; ``max_deviation`` (default 0.05 mm) is enforced.
3. Mass error from (2): ``mass_g`` sets density = mass / mesh volume, so the
   attenuating mass is exactly right even if the volume is slightly off.
4. Touching meshes have non-matching triangulations (gaps/overlaps at the
   interface). Pairwise overlap volumes are measured and must be ~0; objects
   that must touch exactly should come from STEP (shared exact surfaces) instead.
5. Overlap with the detector or CSG objects is checked by sampling.

Running needs an OpenMC build with DAGMC (``openmc.lib.dagmc_enabled``), e.g.
conda-forge's ``openmc`` DAGMC variant.
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np

from .spec import ConfigError, Spec

UNIT_TO_CM = {"mm": 0.1, "cm": 1.0, "m": 100.0, "in": 2.54}


def mesh_report(mesh, crease_deg: float = 30.0) -> dict:
    """Integrity and faceting diagnostics of a closed triangle mesh (lengths in mesh units)."""
    rep = dict(n_faces=int(len(mesh.faces)), n_vertices=int(len(mesh.vertices)),
               watertight=bool(mesh.is_watertight), winding_consistent=bool(mesh.is_winding_consistent),
               volume=float(mesh.volume) if mesh.is_watertight else float("nan"),
               area=float(mesh.area), bodies=int(mesh.body_count),
               euler_number=int(mesh.euler_number))
    areas = mesh.area_faces
    rep["degenerate_faces"] = int(np.sum(areas < 1e-14 * max(mesh.area, 1e-30)))
    # faceting deviation from dihedral angles
    ang = np.asarray(mesh.face_adjacency_angles)
    pairs = np.asarray(mesh.face_adjacency)
    edges = np.asarray(mesh.face_adjacency_edges)
    if len(ang):
        elen = np.linalg.norm(mesh.vertices[edges[:, 0]] - mesh.vertices[edges[:, 1]], axis=1)
        elen = np.maximum(elen, 1e-300)
        h1 = 2 * areas[pairs[:, 0]] / elen
        h2 = 2 * areas[pairs[:, 1]] / elen
        c = 0.5 * (h1 + h2)
        smooth = (ang > 1e-9) & (ang < math.radians(crease_deg))
        s = np.where(smooth, 0.5 * c * np.tan(ang / 4.0), 0.0)
        # convex edges (outer surfaces) lose material, concave ones (cavities) gain it
        sign = np.where(np.asarray(mesh.face_adjacency_convex), 1.0, -1.0)
        dV = np.where(smooth, sign * (2.0 / 3.0) * c * s * elen, 0.0)
        rep["max_chord_deviation"] = float(s.max()) if s.size else 0.0
        rep["mean_chord_deviation"] = float(s[smooth].mean()) if np.any(smooth) else 0.0
        rep["volume_deficit_est"] = float(dV.sum())
        rep["volume_deficit_rel_est"] = float(dV.sum() / rep["volume"]) if rep["watertight"] and rep["volume"] else float("nan")
        rep["smooth_edges"] = int(smooth.sum())
    else:
        rep.update(max_chord_deviation=0.0, mean_chord_deviation=0.0, volume_deficit_est=0.0,
                   volume_deficit_rel_est=0.0, smooth_edges=0)
    return rep


def _rotation(ang):
    from .shapes import _rotation_matrix
    return _rotation_matrix(ang)


def load_stl_object(oname: str, s: Spec, cfg, lib):
    try:
        import trimesh
    except ImportError:
        raise ConfigError("STL support needs 'trimesh' (pip install trimesh rtree)") from None
    from .model import BuiltObject, _density_from_mass
    path = cfg.resolve_path(s.get("file"))
    units = s.choice("units", tuple(UNIT_TO_CM), "mm")
    mname = s.get("material")
    mass = s.get("mass_g", None)
    dens = s.get("density", None)
    placement = Spec(s.get("placement", {}) or {}, f"objects.{oname}.placement")
    checks = Spec(s.get("checks", {}) or {}, f"objects.{oname}.checks")
    s.done()
    max_dev = checks.num("max_deviation", 0.05, positive=True) * cfg.scale  # config units -> cm
    crease = checks.num("crease_angle_deg", 30.0, positive=True)
    checks.done()
    if not Path(path).is_file():
        raise ConfigError(f"STL file not found: {path}")
    mesh = trimesh.load(str(path), force="mesh", process=True)
    mesh.merge_vertices()
    mesh.apply_scale(UNIT_TO_CM[units])
    notes = []
    if not mesh.is_watertight:
        raise ConfigError(
            f"STL '{path.name}' is not watertight (holes or non-manifold edges). A leaky mesh "
            f"loses particles in DAGMC. Repair it (e.g. MeshLab/Blender 3D-print toolbox) and re-export.")
    if not mesh.is_winding_consistent:
        raise ConfigError(f"STL '{path.name}' has inconsistent triangle winding — repair before use.")
    if mesh.volume < 0:
        mesh.invert()
        notes.append(f"STL '{path.name}': normals pointed inwards; inverted")
    if mesh.body_count != 1:
        notes.append(f"STL '{path.name}' contains {mesh.body_count} separate bodies; all get material '{mname}'")
    rot = placement.vec("rotate", 3, [0.0, 0.0, 0.0])
    pivot = np.asarray(placement.vec("pivot", 3, [0.0, 0.0, 0.0])) * cfg.scale
    trans = np.asarray(placement.vec("translate", 3, [0.0, 0.0, 0.0])) * cfg.scale
    align = placement.choice("align", ("none", "bottom_on_endcap"), "none")
    placement.done()
    if any(rot):
        T = np.eye(4)
        R = _rotation(rot)
        T[:3, :3] = R
        T[:3, 3] = pivot - R @ pivot
        mesh.apply_transform(T)
    shift = trans.copy()
    if align == "bottom_on_endcap":
        shift[2] -= mesh.bounds[0][2]
    if np.any(shift):
        mesh.apply_translation(shift)
    rep = mesh_report(mesh, crease)
    if rep["max_chord_deviation"] > max_dev:
        raise ConfigError(
            f"STL '{path.name}': estimated faceting deviation {10 * rep['max_chord_deviation']:.3f} mm "
            f"exceeds checks.max_deviation ({10 * max_dev:.3f} mm); estimated volume deficit "
            f"{100 * rep['volume_deficit_rel_est']:.2f}%. Re-export with a finer chord tolerance, or "
            f"use STEP -> CSG.")
    notes.append(f"STL '{path.name}': {rep['n_faces']} triangles, max chord deviation "
                 f"~{10 * rep['max_chord_deviation']:.4f} mm, volume deficit ~"
                 f"{100 * rep['volume_deficit_rel_est']:.3f}%")
    md = lib.resolve(mname)
    vol = float(mesh.volume)
    rho = None
    if mass is not None:
        if md is None:
            raise ConfigError(f"object '{oname}': mass given for vacuum")
        rho = _density_from_mass(oname, float(mass), vol, 0.0, md.density, notes)
        notes.append(f"density from mass / mesh volume = {rho:.5g} g/cm3 (mass-conserving)")
    elif dens is not None:
        rho = float(dens)
    mat = lib.make_openmc(mname, oname, density=rho)
    if mat is None:
        raise ConfigError(f"STL object '{oname}' cannot be vacuum")
    obj = BuiltObject(name=oname, kind="stl", material_name=mname, material=mat, cell=None, region=None,
                      bbox=(np.asarray(mesh.bounds[0]), np.asarray(mesh.bounds[1])), volume=vol,
                      volume_source="mesh", mass_g=None if mass is None else float(mass),
                      density=mat.density, mesh=mesh, notes=notes)
    obj.mesh_report = rep
    return obj


def stl_pair_overlap(a, b) -> float:
    """Overlap volume (cm3) of two closed meshes (manifold3d boolean, else sampling)."""
    lo = np.maximum(a.mesh.bounds[0], b.mesh.bounds[0])
    hi = np.minimum(a.mesh.bounds[1], b.mesh.bounds[1])
    if np.any(hi <= lo):
        return 0.0
    try:
        inter = a.mesh.intersection(b.mesh, engine="manifold")
        return float(abs(inter.volume)) if inter is not None and len(inter.faces) else 0.0
    except Exception:
        import trimesh
        pts = trimesh.sample.volume_mesh(a.mesh, 20000)
        frac = float(np.mean(b.mesh.contains(pts))) if len(pts) else 0.0
        return frac * a.volume


def stl_min_distance(a, b, n=20000) -> float:
    """Smallest distance (cm) between the surfaces of two meshes (sampled)."""
    import trimesh
    pts, _ = trimesh.sample.sample_surface(a.mesh, n, seed=1)
    _, dist, _ = trimesh.proximity.closest_point(b.mesh, pts)
    return float(np.min(dist))


def build_dagmc_universe(stl_objects, cfg, ambient_name: str, setup_name: str, warnings=None):
    """Write all STL objects into one DAGMC .h5m and return (universe, bounding region)."""
    import openmc
    for i, a in enumerate(stl_objects):
        for b in stl_objects[i + 1:]:
            ov = stl_pair_overlap(a, b)
            if ov > 1e-4 * min(a.volume, b.volume):
                raise ConfigError(f"STL objects '{a.name}' and '{b.name}' overlap by {ov:.4g} cm3")
            dmin = stl_min_distance(a, b)
            if dmin < 1e-3 and warnings is not None:
                warnings.append(
                    f"STL objects '{a.name}' and '{b.name}' touch (gap {dmin * 1e4:.1f} um) with separately "
                    f"triangulated surfaces — the classic source of DAGMC lost particles. Watch the "
                    f"lost-particle count, or use STEP -> CSG for parts that share a boundary.")
    h5m = cfg.output_dir / "_dagmc" / f"{setup_name}.h5m"
    h5m.parent.mkdir(parents=True, exist_ok=True)
    verts, tri_sets, tags = [], {}, []
    off = 0
    for sid, o in enumerate(stl_objects, start=1):
        v = np.asarray(o.mesh.vertices, float)
        f = np.asarray(o.mesh.faces, int) + off
        verts.append(v)
        tri_sets[sid] = {sid: f.tolist()}
        tags.append(o.material.name)
        off += len(v)
    vertices = [tuple(p) for p in np.concatenate(verts)]
    try:
        from cad_to_dagmc.core import vertices_to_h5m
    except Exception as exc:
        raise ConfigError(f"writing DAGMC files needs 'cad_to_dagmc' (pip install cad_to_dagmc): {exc}") from None
    try:
        import pymoab  # noqa: F401
        method = "pymoab"
    except Exception:
        method = "h5py"
    amb = "void" if ambient_name in ("vacuum", "void") else "ambient"
    vertices_to_h5m(vertices=vertices, triangles_by_solid_by_face=tri_sets, material_tags=tags,
                    h5m_filename=str(h5m), implicit_complement_material_tag=amb, method=method)
    univ = openmc.DAGMCUniverse(filename=str(h5m), auto_geom_ids=True, name=f"stl_{setup_name}")
    lo = np.min([o.bbox[0] for o in stl_objects], axis=0) - 0.05
    hi = np.max([o.bbox[1] for o in stl_objects], axis=0) + 0.05
    box = openmc.model.RectangularParallelepiped(lo[0], hi[0], lo[1], hi[1], lo[2], hi[2])
    return univ, -box


def dagmc_available() -> bool:
    try:
        import openmc.lib
        return bool(openmc.lib.dagmc_enabled)
    except Exception:
        return False
