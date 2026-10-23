"""Photon source models.

Types (lengths in the configuration unit):

* ``point``          — ``position: [x, y, z]``
* ``disk``           — ``radius``, ``z``, optional ``inner_radius`` (annulus), ``center``
* ``ring``           — thin circular line source: ``radius``, ``z``, ``center``
* ``volume``         — uniform in an object: ``object: <name>``
* ``surface``        — uniform on object faces: ``object``, ``faces: [inner_side, ...]`` or ``all``
                       (CSG faces; STL objects: the whole mesh surface)
* ``depth_profile``  — volume source whose activity varies with height (sedimentation, ...):
                       ``object``, ``profile: {type: exponential, length: 5, from: bottom}``,
                       ``{type: linear, bottom: 1, top: 0}`` or ``{z: [...], weight: [...]}``
                       (z measured from the object's lowest point)
* ``points``         — user-defined emission points: ``file: points.csv`` (x, y, z[, weight])
                       or ``positions: [[x,y,z], ...]``
* ``mixture``        — ``components: [{weight: 0.7, type: volume, ...}, {weight: 0.3, ...}]``

Uniform sampling uses exact distributions where the shape allows it, otherwise
OpenMC rejection sampling inside the object's bounding box. If rejection would
accept fewer than ~10 % of samples (OpenMC aborts below 5 %), a point cloud
drawn from the exact region is used instead.
"""
from __future__ import annotations

import math

import numpy as np

from .regions import rejection_points
from .spec import ConfigError, Spec

TWO_PI = 2.0 * math.pi
_CLOUD_POINTS = 200_000


def _isotropic(bias):
    import openmc.stats as st
    if bias is None:
        return st.Isotropic()
    theta, _ = bias
    return st.PolarAzimuthal(mu=st.Uniform(math.cos(theta), 1.0), phi=st.Uniform(0.0, TWO_PI),
                             reference_uvw=(0.0, 0.0, -1.0))


def _src(space, E, strength, bias, constraints=None):
    import openmc
    import openmc.stats as st
    kw = {}
    if constraints:
        kw["constraints"] = constraints
    return openmc.IndependentSource(space=space, angle=_isotropic(bias),
                                    energy=st.Discrete([E], [1.0]), particle="photon",
                                    strength=float(strength), **kw)


def build_sources(raw: dict, built, E: float, bias=None):
    """Return (list of openmc.IndependentSource, meta) for one energy E (eV)."""
    if raw is None:
        raise ConfigError(f"setup '{built.name}': 'source' is required")
    pieces, meta = _build(Spec(raw, f"setups.{built.name}.source"), built, E, bias, 1.0)
    tot = sum(p.strength for p in pieces)
    for p in pieces:
        p.strength = p.strength / tot
    return pieces, meta


def _build(s: Spec, built, E, bias, weight):
    import openmc.stats as st
    k = built.cfg.scale
    typ = s.choice("type", ("point", "disk", "ring", "volume", "surface", "depth_profile",
                            "points", "mixture"))
    meta = {"source_type": typ}
    if typ == "point":
        pos = [v * k for v in s.vec("position", 3)]
        s.done()
        return [_src(st.Point(pos), E, weight, bias)], meta
    if typ in ("disk", "ring"):
        R = s.num("radius", positive=True) * k
        z = s.num("z") * k
        c = [v * k for v in s.vec("center", 2, [0.0, 0.0])]
        r0 = s.num("inner_radius", 0.0, nonneg=True) * k if typ == "disk" else R
        s.done()
        rdist = st.PowerLaw(r0, R, 1) if typ == "disk" else st.Discrete([R], [1.0])
        sp = st.CylindricalIndependent(r=rdist, phi=st.Uniform(0.0, TWO_PI),
                                       z=st.Discrete([z], [1.0]), origin=(c[0], c[1], 0.0))
        return [_src(sp, E, weight, bias)], meta
    if typ == "volume":
        obj = built.obj(s.get("object"))
        npts = int(s.num("cloud_points", _CLOUD_POINTS, positive=True))
        how = s.choice("sampling", ("auto", "rejection"), "auto")
        s.done()
        meta["object"] = obj.name
        return _volume(obj, built, E, bias, weight, npts, meta, force_rejection=(how == "rejection")), meta
    if typ == "surface":
        obj = built.obj(s.get("object"))
        faces = s.get("faces", "all")
        npts = int(s.num("cloud_points", _CLOUD_POINTS, positive=True))
        s.done()
        meta["object"] = obj.name
        return _surface(obj, faces, E, bias, weight, npts), meta
    if typ == "depth_profile":
        obj = built.obj(s.get("object"))
        prof = s.sub("profile")
        npts = int(s.num("cloud_points", _CLOUD_POINTS, positive=True))
        s.done()
        meta["object"] = obj.name
        meta["sampling"] = "z-profile + rejection by cell/material"
        return _depth(obj, prof, built, E, bias, weight, npts), meta
    if typ == "points":
        if s.has("file"):
            path = built.cfg.resolve_path(s.get("file"))
            arr = _read_points(path)
        else:
            arr = np.asarray(s.get("positions"), float)
        s.done()
        if arr.ndim != 2 or arr.shape[1] not in (3, 4):
            raise ConfigError("points source needs rows of x, y, z[, weight]")
        pos = arr[:, :3] * k
        wts = arr[:, 3] if arr.shape[1] == 4 else np.ones(len(arr))
        if np.any(wts < 0) or wts.sum() <= 0:
            raise ConfigError("points source weights must be non-negative with a positive sum")
        sp = st.PointCloud(pos, strengths=wts)
        return [_src(sp, E, weight, bias)], meta
    # mixture
    comps = s.get("components")
    s.done()
    if not isinstance(comps, list) or not comps:
        raise ConfigError("mixture source needs a non-empty 'components' list")
    out = []
    ws = []
    for i, c in enumerate(comps):
        cs = Spec(c, f"{s.path}.components[{i}]")
        w = cs.num("weight", positive=True)
        ws.append(w)
    wsum = sum(ws)
    for i, c in enumerate(comps):
        c2 = dict(c)
        c2.pop("weight")
        sub, _ = _build(Spec(c2, f"{s.path}.components[{i}]"), built, E, bias, 1.0)
        tot = sum(p.strength for p in sub)
        for p in sub:
            p.strength = weight * (ws[i] / wsum) * p.strength / tot
        out += sub
    return out, meta


def _read_points(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.replace(",", " ").replace(";", " ").split()
            try:
                rows.append([float(p) for p in parts])
            except ValueError:
                continue  # header
    return np.asarray(rows, float)


def _acceptance(obj):
    ll, ur = obj.bbox
    vb = float(np.prod(np.asarray(ur) - np.asarray(ll)))
    return (obj.volume / vb) if (obj.volume and vb > 0) else None


def _volume(obj, built, E, bias, weight, npts, meta, force_rejection=False):
    import openmc.stats as st
    sh = obj.shape
    if sh is not None and sh.volume_parts and not force_rejection:
        tot = sum(w for _, w in sh.volume_parts)
        meta["sampling"] = "exact"
        return [_src(b(), E, weight * w / tot, bias) for b, w in sh.volume_parts]
    acc = _acceptance(obj)
    ll, ur = obj.bbox
    if obj.kind == "stl":
        if acc is not None and acc < 0.10:
            pts = _stl_volume_points(obj, npts)
            meta["sampling"] = f"point cloud ({len(pts)} pts in mesh)"
            return [_src(st.PointCloud(pts), E, weight, bias)]
        meta["sampling"] = f"rejection by material (acceptance {acc:.2f})"
        return [_src(st.Box(list(ll), list(ur)), E, weight, bias,
                     constraints={"domains": [obj.material]})]
    if acc is not None and acc < 0.10:
        pts = rejection_points(obj.region, ll, ur, npts, np.random.default_rng(2024))
        meta["sampling"] = f"point cloud ({len(pts)} pts, acceptance {acc:.3f})"
        return [_src(st.PointCloud(pts), E, weight, bias)]
    meta["sampling"] = f"rejection by cell (acceptance {acc:.2f})" if acc else "rejection by cell"
    return [_src(st.Box(list(ll), list(ur)), E, weight, bias,
                 constraints={"domains": [obj.cell]})]


def _surface(obj, faces, E, bias, weight, npts):
    import openmc.stats as st
    if obj.kind == "stl":
        import trimesh
        pts, _ = trimesh.sample.sample_surface(obj.mesh, npts, seed=2024)
        return [_src(st.PointCloud(np.asarray(pts)), E, weight, bias)]
    sh = obj.shape
    if sh is None or not sh.faces:
        raise ConfigError(
            f"surface source on '{obj.name}': this object has no analytic faces "
            f"(available for exact cylinder/tube/cup/box/sphere shapes and STL meshes). "
            f"Model a thin contamination layer as its own object and use a volume source.")
    names = list(sh.faces) if faces in ("all", None) else list(faces)
    bad = [f for f in names if f not in sh.faces]
    if bad:
        raise ConfigError(f"surface source on '{obj.name}': unknown face(s) {bad}; "
                          f"available: {list(sh.faces)}")
    areas = [sh.faces[f][1] for f in names]
    tot = sum(areas)
    return [_src(sh.faces[f][0](), E, weight * a / tot, bias) for f, a in zip(names, areas)]


def _profile_weights(prof: Spec, z0, z1, k):
    n = 801
    z = np.linspace(z0, z1, n)
    if prof.has("z"):
        zz = np.asarray(prof.get("z"), float) * k + z0
        ww = np.asarray(prof.get("weight"), float)
        prof.done()
        if len(zz) != len(ww) or len(zz) < 2:
            raise ConfigError("depth profile table needs matching z/weight lists (>= 2 points)")
        w = np.interp(z, zz, ww, left=0.0, right=0.0)
    else:
        typ = prof.choice("type", ("exponential", "linear", "uniform"))
        if typ == "exponential":
            lam = prof.num("length", positive=True) * k
            frm = prof.choice("from", ("bottom", "top"), "bottom")
            w = np.exp(-(z - z0) / lam) if frm == "bottom" else np.exp(-(z1 - z) / lam)
        elif typ == "linear":
            wb, wt = prof.num("bottom", nonneg=True), prof.num("top", nonneg=True)
            w = wb + (wt - wb) * (z - z0) / (z1 - z0)
        else:
            w = np.ones_like(z)
        prof.done()
    if np.any(w < 0) or w.sum() <= 0:
        raise ConfigError("depth profile must be non-negative and not all zero")
    return z, w


def _depth(obj, prof, built, E, bias, weight, npts):
    import openmc.stats as st
    ll, ur = obj.bbox
    z0, z1 = float(ll[2]), float(ur[2])
    z, w = _profile_weights(prof, z0, z1, built.cfg.scale)
    acc = _acceptance(obj)
    if acc is not None and acc < 0.10:
        rng = np.random.default_rng(2025)
        if obj.kind == "stl":
            pts = _stl_volume_points(obj, npts * 3)
        else:
            pts = rejection_points(obj.region, ll, ur, npts * 3, rng)
        wz = np.interp(pts[:, 2], z, w)
        keep = rng.random(len(pts)) < wz / wz.max()
        return [_src(st.PointCloud(pts[keep]), E, weight, bias)]
    space = st.CartesianIndependent(st.Uniform(float(ll[0]), float(ur[0])),
                                    st.Uniform(float(ll[1]), float(ur[1])),
                                    st.Tabular(z, w, interpolation="linear-linear"))
    dom = obj.material if obj.kind == "stl" else obj.cell
    return [_src(space, E, weight, bias, constraints={"domains": [dom]})]


def _stl_volume_points(obj, n):
    import trimesh
    out, got = [], 0
    for _ in range(100):
        pts = trimesh.sample.volume_mesh(obj.mesh, max(4 * (n - got), 10000))
        out.append(pts)
        got += len(pts)
        if got >= n:
            break
    return np.concatenate(out)[:n]


# ----------------------------------------------------------------------------
# Angular biasing (optional, off by default)
# ----------------------------------------------------------------------------
def source_bbox(raw, built):
    """Conservative bounding box of all emission points of a source spec (cm)."""
    k = built.cfg.scale
    s = dict(raw)
    typ = s.get("type")
    if typ == "point":
        p = np.asarray(s["position"], float) * k
        return p, p
    if typ in ("disk", "ring"):
        R = float(s["radius"]) * k
        c = np.asarray(s.get("center", [0, 0]), float) * k
        z = float(s["z"]) * k
        return np.array([c[0] - R, c[1] - R, z]), np.array([c[0] + R, c[1] + R, z])
    if typ in ("volume", "surface", "depth_profile"):
        o = built.obj(s["object"])
        return np.asarray(o.bbox[0], float), np.asarray(o.bbox[1], float)
    if typ == "points":
        arr = (_read_points(built.cfg.resolve_path(s["file"])) if "file" in s
               else np.asarray(s["positions"], float))
        p = arr[:, :3] * k
        return p.min(axis=0), p.max(axis=0)
    if typ == "mixture":
        boxes = [source_bbox({kk: v for kk, v in c.items() if kk != "weight"}, built)
                 for c in s["components"]]
        return np.min([b[0] for b in boxes], axis=0), np.max([b[1] for b in boxes], axis=0)
    raise ConfigError(f"cannot bound source type {typ}")


def angular_bias_cone(built, margin_deg: float = 2.0):
    """Half-angle (rad) of a downward cone containing every straight line from any
    emission point to any point of the crystal, and the solid-angle factor.

    Returns None when no useful cone exists (source not entirely above the crystal).
    Exact for unscattered photons; it drops only photons that would have to scatter
    elastically (coherent) back towards the detector — documented as negligible
    and checked in the validation.
    """
    sll, sur = source_bbox(built.source_spec, built)
    det = built.detector
    R = det.info["R"]
    cll = np.array([-R, -R, det.info["z_back"]])
    cur = np.array([R, R, det.info["z_front"]])
    dz_min = sll[2] - cur[2]  # smallest height of a source point above the crystal
    if dz_min <= 0:
        return None
    dx = max(abs(cur[0] - sll[0]), abs(cll[0] - sur[0]))
    dy = max(abs(cur[1] - sll[1]), abs(cll[1] - sur[1]))
    theta = math.atan2(math.hypot(dx, dy), dz_min) + math.radians(margin_deg)
    if theta >= math.radians(89.0):
        theta = math.pi / 2  # hemisphere is still exact when the crystal lies wholly below
    return theta, (1.0 - math.cos(theta)) / 2.0
