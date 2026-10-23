"""Exact CSG shapes for containers, samples and other objects.

All inputs are in the configuration length unit (mm by default) and are
converted to cm (OpenMC's unit) here. Every shape knows

* its OpenMC region (quadric surfaces — no triangles),
* its bounding box,
* its exact volume when that has a closed form,
* exact uniform source distributions for its volume and named faces when they
  exist (otherwise the source module falls back to rejection sampling),
* an axisymmetric "tube segment" decomposition used by the independent
  attenuation cross-check (:mod:`hpgeff.esa`).

Shapes that share a boundary (a cup and its contents) are built from the *same*
surface objects, so the two cells meet exactly: there is no gap or overlap,
the CSG equivalent of a conformal mesh.

``facets: n`` on cylinders/tubes/cups replaces each circle by the inscribed
regular n-gon — exactly what an STL export does — so the effect of
triangulation can be measured inside pure CSG.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Callable

import numpy as np

from .spec import ConfigError, Spec

TWO_PI = 2.0 * math.pi


@dataclass
class Shape:
    name: str
    kind: str
    region: object  # openmc.Region
    bbox: tuple  # (ll[3], ur[3]) in cm
    volume: float | None = None  # cm3 (exact if known)
    esa: list | None = None  # [(r_in, r_out, z0, z1)] (cm) about the z axis through (0,0)
    volume_parts: list | None = None  # [(spatial_builder, weight)] exact uniform sampling
    faces: dict = field(default_factory=dict)  # name -> (spatial_builder, area cm2)
    z_range: tuple | None = None  # (zmin, zmax) for depth profiles
    params: dict = field(default_factory=dict)

    def bbox_volume(self) -> float:
        ll, ur = self.bbox
        return float(np.prod(np.asarray(ur) - np.asarray(ll)))


# ----------------------------------------------------------------------------
# Radial boundaries: exact circle or inscribed regular polygon ("faceted")
# ----------------------------------------------------------------------------
class Radial:
    """Boundary r = R around (cx, cy): a ZCylinder, or an inscribed n-gon of planes."""

    def __init__(self, R: float, cx: float = 0.0, cy: float = 0.0, facets: int = 0,
                 phase: float = 0.0):
        import openmc
        self.R, self.cx, self.cy, self.n = float(R), float(cx), float(cy), int(facets or 0)
        if self.n:
            if self.n < 3:
                raise ConfigError("facets must be 0 (exact) or >= 3")
            a = self.R * math.cos(math.pi / self.n)  # apothem
            self.planes = []
            for k in range(self.n):
                ang = phase + TWO_PI * (k + 0.5) / self.n  # outward normal direction
                nx, ny = math.cos(ang), math.sin(ang)
                self.planes.append(openmc.Plane(a=nx, b=ny, c=0.0, d=a + nx * self.cx + ny * self.cy))
            self.surface = None
        else:
            self.surface = openmc.ZCylinder(x0=self.cx, y0=self.cy, r=self.R)

    @property
    def area(self) -> float:
        if self.n:
            return 0.5 * self.n * self.R ** 2 * math.sin(TWO_PI / self.n)
        return math.pi * self.R ** 2

    def inside(self):
        import openmc
        if self.n:
            return openmc.Intersection([-p for p in self.planes])
        return -self.surface

    def outside(self):
        import openmc
        if self.n:
            return openmc.Union([+p for p in self.planes])
        return +self.surface


# ----------------------------------------------------------------------------
# Source-distribution builders (return callables so OpenMC objects are fresh)
# ----------------------------------------------------------------------------
def _cyl_volume_space(r0, r1, z0, z1, cx, cy):
    def build():
        import openmc.stats as st
        return st.CylindricalIndependent(
            r=st.PowerLaw(r0, r1, 1), phi=st.Uniform(0.0, TWO_PI), z=st.Uniform(z0, z1),
            origin=(cx, cy, 0.0))
    return build


def _disk_space(r0, r1, z, cx, cy):
    def build():
        import openmc.stats as st
        return st.CylindricalIndependent(
            r=st.PowerLaw(r0, r1, 1), phi=st.Uniform(0.0, TWO_PI), z=st.Discrete([z], [1.0]),
            origin=(cx, cy, 0.0))
    return build


def _lateral_space(R, z0, z1, cx, cy):
    def build():
        import openmc.stats as st
        return st.CylindricalIndependent(
            r=st.Discrete([R], [1.0]), phi=st.Uniform(0.0, TWO_PI), z=st.Uniform(z0, z1),
            origin=(cx, cy, 0.0))
    return build


def _box_space(ll, ur):
    def build():
        import openmc.stats as st
        return st.Box(list(ll), list(ur))
    return build


def _sphere_volume_space(R, c):
    def build():
        import openmc.stats as st
        return st.SphericalIndependent(r=st.PowerLaw(0.0, R, 2), cos_theta=st.Uniform(-1.0, 1.0),
                                       phi=st.Uniform(0.0, TWO_PI), origin=tuple(c))
    return build


def _sphere_surface_space(R, c):
    def build():
        import openmc.stats as st
        return st.SphericalIndependent(r=st.Discrete([R], [1.0]), cos_theta=st.Uniform(-1.0, 1.0),
                                       phi=st.Uniform(0.0, TWO_PI), origin=tuple(c))
    return build


def _bbox_cyl(R, z0, z1, cx, cy):
    return (np.array([cx - R, cy - R, z0]), np.array([cx + R, cy + R, z1]))


# ----------------------------------------------------------------------------
# Shape builders
# ----------------------------------------------------------------------------
class _Ctx:
    def __init__(self, s: Spec, scale: float):
        self.s, self.k = s, scale

    def L(self, key, default=..., positive=False, nonneg=False):
        if default is ...:
            v = self.s.num(key, positive=positive, nonneg=nonneg)
        else:
            v = self.s.num(key, default, positive=positive, nonneg=nonneg)
        return None if v is None else v * self.k

    def center2(self):
        c = self.s.vec("center", 2, [0.0, 0.0])
        return c[0] * self.k, c[1] * self.k


def build_cylinder(name, s: Spec, k: float) -> list[Shape]:
    import openmc
    c = _Ctx(s, k)
    R, H = c.L("radius", positive=True), c.L("height", positive=True)
    z0 = c.L("z", 0.0)
    cx, cy = c.center2()
    nf = int(s.num("facets", 0))
    rad = Radial(R, cx, cy, nf)
    p0, p1 = openmc.ZPlane(z0=z0), openmc.ZPlane(z0=z0 + H)
    sh = Shape(name, "cylinder", rad.inside() & +p0 & -p1, _bbox_cyl(R, z0, z0 + H, cx, cy),
               volume=rad.area * H, z_range=(z0, z0 + H),
               params=dict(R=R, H=H, z0=z0, cx=cx, cy=cy, facets=nf))
    if not nf:
        sh.volume_parts = [(_cyl_volume_space(0.0, R, z0, z0 + H, cx, cy), 1.0)]
        sh.faces = {"top": (_disk_space(0.0, R, z0 + H, cx, cy), math.pi * R * R),
                    "bottom": (_disk_space(0.0, R, z0, cx, cy), math.pi * R * R),
                    "side": (_lateral_space(R, z0, z0 + H, cx, cy), TWO_PI * R * H)}
        if cx == 0 and cy == 0:
            sh.esa = [(0.0, R, z0, z0 + H)]
    return [sh]


def build_tube(name, s: Spec, k: float) -> list[Shape]:
    import openmc
    c = _Ctx(s, k)
    Ri, Ro = c.L("inner_radius", nonneg=True), c.L("outer_radius", positive=True)
    H, z0 = c.L("height", positive=True), c.L("z", 0.0)
    cx, cy = c.center2()
    nf = int(s.num("facets", 0))
    if not Ro > Ri:
        raise ConfigError(f"object '{name}': outer_radius must exceed inner_radius")
    ro, ri = Radial(Ro, cx, cy, nf), Radial(Ri, cx, cy, nf) if Ri > 0 else None
    p0, p1 = openmc.ZPlane(z0=z0), openmc.ZPlane(z0=z0 + H)
    reg = ro.inside() & +p0 & -p1
    if ri is not None:
        reg = reg & ri.outside()
    vol = (ro.area - (ri.area if ri else 0.0)) * H
    sh = Shape(name, "tube", reg, _bbox_cyl(Ro, z0, z0 + H, cx, cy), volume=vol, z_range=(z0, z0 + H),
               params=dict(Ri=Ri, Ro=Ro, H=H, z0=z0, cx=cx, cy=cy, facets=nf))
    if not nf:
        sh.volume_parts = [(_cyl_volume_space(Ri, Ro, z0, z0 + H, cx, cy), 1.0)]
        sh.faces = {"top": (_disk_space(Ri, Ro, z0 + H, cx, cy), math.pi * (Ro ** 2 - Ri ** 2)),
                    "bottom": (_disk_space(Ri, Ro, z0, cx, cy), math.pi * (Ro ** 2 - Ri ** 2)),
                    "outer_side": (_lateral_space(Ro, z0, z0 + H, cx, cy), TWO_PI * Ro * H)}
        if Ri > 0:
            sh.faces["inner_side"] = (_lateral_space(Ri, z0, z0 + H, cx, cy), TWO_PI * Ri * H)
        if cx == 0 and cy == 0:
            sh.esa = [(Ri, Ro, z0, z0 + H)]
    return [sh]


def build_cup(name, s: Spec, k: float, fill_name_default: str) -> list[Shape]:
    """Container: flat bottom + cylindrical wall, optional lid, optional contents.

    ``z`` is the outer bottom face. ``fill`` creates a second object (the sample)
    that shares the cup's inner surfaces exactly.
    """
    import openmc
    c = _Ctx(s, k)
    Ri = c.L("inner_radius", positive=True)
    t, b = c.L("wall", positive=True), c.L("bottom", nonneg=True)
    Hin = c.L("inner_height", positive=True)
    z0 = c.L("z", 0.0)
    cx, cy = c.center2()
    nf = int(s.num("facets", 0))
    Ro = Ri + t
    zb, zt = z0 + b, z0 + b + Hin
    ro, ri = Radial(Ro, cx, cy, nf), Radial(Ri, cx, cy, nf)
    p0, pb, pt = openmc.ZPlane(z0=z0), openmc.ZPlane(z0=zb), openmc.ZPlane(z0=zt)
    wall = ri.outside() & ro.inside() & +pb & -pt
    pieces = [wall]
    esa = [(Ri, Ro, zb, zt)]
    vol = (ro.area - ri.area) * Hin
    if b > 0:
        pieces.append(ro.inside() & +p0 & -pb)
        esa.append((0.0, Ro, z0, zb))
        vol += ro.area * b
    ztop = zt
    lid = s.sub("lid", None) if s.has("lid") else None
    out: list[Shape] = []
    if lid is not None:
        lt = lid.num("thickness", positive=True) * k
        lid_mat = lid.get("material", None)
        lid_mass = lid.get("mass_g", None)
        lid.done()
        plid = openmc.ZPlane(z0=zt + lt)
        lid_reg = ro.inside() & +pt & -plid
        ztop = zt + lt
        if lid_mat is None and lid_mass is None:
            pieces.append(lid_reg)
            esa.append((0.0, Ro, zt, zt + lt))
            vol += ro.area * lt
        else:
            ztop = zt  # separate lid object; the cup itself ends at the rim
            lsh = Shape(f"{name}_lid", "lid", lid_reg, _bbox_cyl(Ro, zt, zt + lt, cx, cy),
                        volume=ro.area * lt, z_range=(zt, zt + lt),
                        esa=[(0.0, Ro, zt, zt + lt)] if (cx == 0 and cy == 0 and not nf) else None)
            lsh.params = {"material": lid_mat, "mass_g": lid_mass}
            out.append(lsh)
    region = openmc.Union(pieces) if len(pieces) > 1 else pieces[0]
    cup = Shape(name, "cup", region, _bbox_cyl(Ro, z0, ztop, cx, cy), volume=vol, z_range=(z0, ztop),
                params=dict(Ri=Ri, Ro=Ro, wall=t, bottom=b, Hin=Hin, z0=z0, cx=cx, cy=cy, facets=nf))
    if not nf:
        cup.faces = {
            "inner_side": (_lateral_space(Ri, zb, zt, cx, cy), TWO_PI * Ri * Hin),
            "inner_bottom": (_disk_space(0.0, Ri, zb, cx, cy), math.pi * Ri * Ri),
            "outer_side": (_lateral_space(Ro, z0, zt, cx, cy), TWO_PI * Ro * (zt - z0)),
            "outer_bottom": (_disk_space(0.0, Ro, z0, cx, cy), math.pi * Ro * Ro),
            "rim": (_disk_space(Ri, Ro, zt, cx, cy), math.pi * (Ro * Ro - Ri * Ri)),
        }
        if cx == 0 and cy == 0:
            cup.esa = esa
    out.insert(0, cup)

    if s.has("fill"):
        f = s.sub("fill")
        fname = str(f.get("name", fill_name_default))
        hf = f.num("height", Hin / k, positive=True) * k
        if hf > Hin * (1 + 1e-12):
            raise ConfigError(f"object '{name}': fill height {hf / k} exceeds inner_height {Hin / k}")
        pf = pt if abs(hf - Hin) < 1e-12 else openmc.ZPlane(z0=zb + hf)
        freg = ri.inside() & +pb & -pf
        fsh = Shape(fname, "cup_fill", freg, _bbox_cyl(Ri, zb, zb + hf, cx, cy), volume=ri.area * hf,
                    z_range=(zb, zb + hf), params=dict(R=Ri, H=hf, z0=zb, cx=cx, cy=cy, facets=nf,
                                                        fill=f.raw()))
        for key in ("material", "mass_g", "density"):
            f.get(key, None)
        f.done()
        if not nf:
            fsh.volume_parts = [(_cyl_volume_space(0.0, Ri, zb, zb + hf, cx, cy), 1.0)]
            fsh.faces = {"top": (_disk_space(0.0, Ri, zb + hf, cx, cy), math.pi * Ri * Ri),
                         "bottom": (_disk_space(0.0, Ri, zb, cx, cy), math.pi * Ri * Ri),
                         "side": (_lateral_space(Ri, zb, zb + hf, cx, cy), TWO_PI * Ri * hf)}
            if cx == 0 and cy == 0:
                fsh.esa = [(0.0, Ri, zb, zb + hf)]
        out.append(fsh)
    return out


def build_marinelli(name, s: Spec, k: float, fill_name_default: str) -> list[Shape]:
    """Marinelli (re-entrant) beaker. ``seat_z`` is the inner surface of the well
    ceiling, i.e. where the beaker rests on the detector end cap (default 0)."""
    import openmc
    c = _Ctx(s, k)
    Ri = c.L("inner_radius", positive=True)
    t = c.L("wall", positive=True)
    b = c.L("bottom", t / k)
    Hin = c.L("inner_height", positive=True)
    rw = c.L("well_radius", positive=True)
    dw = c.L("well_depth", positive=True)
    tw = c.L("well_wall", t / k)
    seat = c.L("seat_z", 0.0)
    cx, cy = c.center2()
    Ro, rwo = Ri + t, rw + tw
    zbo = seat - dw
    zb, zt, zwc = zbo + b, zbo + b + Hin, seat + tw
    if not rwo < Ri:
        raise ConfigError(f"object '{name}': well_radius + well_wall must be < inner_radius")
    if not seat > zb:
        raise ConfigError(f"object '{name}': well_depth must exceed the bottom thickness")
    if not zwc < zt:
        raise ConfigError(f"object '{name}': well ceiling is above the rim (check inner_height)")
    cyl = {r: openmc.ZCylinder(x0=cx, y0=cy, r=r) for r in (Ro, Ri, rwo, rw)}
    pl = {z: openmc.ZPlane(z0=z) for z in (zbo, zb, zt, seat, zwc)}
    pieces = [
        +cyl[Ri] & -cyl[Ro] & +pl[zbo] & -pl[zt],             # outer wall
        +cyl[rwo] & -cyl[Ri] & +pl[zbo] & -pl[zb],            # bottom ring
        +cyl[rw] & -cyl[rwo] & +pl[zbo] & -pl[zwc],           # well side wall
        -cyl[rw] & +pl[seat] & -pl[zwc],                      # well ceiling
    ]
    vol = (math.pi * (Ro ** 2 - Ri ** 2) * (zt - zbo) + math.pi * (Ri ** 2 - rwo ** 2) * b
           + math.pi * (rwo ** 2 - rw ** 2) * (zwc - zbo) + math.pi * rw ** 2 * tw)
    beaker = Shape(name, "marinelli", openmc.Union(pieces), _bbox_cyl(Ro, zbo, zt, cx, cy),
                   volume=vol, z_range=(zbo, zt),
                   params=dict(Ri=Ri, Ro=Ro, rw=rw, rwo=rwo, zbo=zbo, zb=zb, zt=zt, seat=seat, zwc=zwc))
    if cx == 0 and cy == 0:
        beaker.esa = [(Ri, Ro, zbo, zt), (rwo, Ri, zbo, zb), (rw, rwo, zbo, zwc), (0.0, rw, seat, zwc)]
    out = [beaker]
    if s.has("fill"):
        f = s.sub("fill")
        fname = str(f.get("name", fill_name_default))
        hf = f.num("height", Hin / k, positive=True) * k
        if hf > Hin * (1 + 1e-12):
            raise ConfigError(f"object '{name}': fill height exceeds inner_height")
        zf = zb + hf
        pf = pl[zt] if abs(zf - zt) < 1e-12 else (pl.get(zf) or openmc.ZPlane(z0=zf))
        parts, vparts, esa = [], [], []
        if zf <= zwc:
            parts.append(+cyl[rwo] & -cyl[Ri] & +pl[zb] & -pf)
            vparts.append((_cyl_volume_space(rwo, Ri, zb, zf, cx, cy), math.pi * (Ri ** 2 - rwo ** 2) * (zf - zb)))
            esa.append((rwo, Ri, zb, zf))
        else:
            parts.append(+cyl[rwo] & -cyl[Ri] & +pl[zb] & -pl[zwc])
            parts.append(-cyl[Ri] & +pl[zwc] & -pf)
            v1 = math.pi * (Ri ** 2 - rwo ** 2) * (zwc - zb)
            v2 = math.pi * Ri ** 2 * (zf - zwc)
            vparts += [(_cyl_volume_space(rwo, Ri, zb, zwc, cx, cy), v1),
                       (_cyl_volume_space(0.0, Ri, zwc, zf, cx, cy), v2)]
            esa += [(rwo, Ri, zb, zwc), (0.0, Ri, zwc, zf)]
        vtot = sum(w for _, w in vparts)
        fsh = Shape(fname, "marinelli_fill", openmc.Union(parts) if len(parts) > 1 else parts[0],
                    _bbox_cyl(Ri, zb, zf, cx, cy), volume=vtot, z_range=(zb, zf),
                    volume_parts=vparts, esa=esa if (cx == 0 and cy == 0) else None,
                    params=dict(fill=f.raw()))
        for key in ("material", "mass_g", "density"):
            f.get(key, None)
        f.done()
        out.append(fsh)
    return out


def build_box(name, s: Spec, k: float) -> list[Shape]:
    import openmc
    c = _Ctx(s, k)
    size = [v * k for v in s.vec("size", 3)]
    if min(size) <= 0:
        raise ConfigError(f"object '{name}': box size must be positive")
    cx, cy = c.center2()
    z0 = c.L("z", 0.0)
    ll = np.array([cx - size[0] / 2, cy - size[1] / 2, z0])
    ur = np.array([cx + size[0] / 2, cy + size[1] / 2, z0 + size[2]])
    xs = [openmc.XPlane(ll[0]), openmc.XPlane(ur[0])]
    ys = [openmc.YPlane(ll[1]), openmc.YPlane(ur[1])]
    zs = [openmc.ZPlane(ll[2]), openmc.ZPlane(ur[2])]
    reg = +xs[0] & -xs[1] & +ys[0] & -ys[1] & +zs[0] & -zs[1]
    sh = Shape(name, "box", reg, (ll, ur), volume=float(np.prod(size)), z_range=(ll[2], ur[2]),
               volume_parts=[(_box_space(ll, ur), 1.0)])
    sx, sy, sz = size
    sh.faces = {
        "top": (_box_space([ll[0], ll[1], ur[2]], [ur[0], ur[1], ur[2]]), sx * sy),
        "bottom": (_box_space([ll[0], ll[1], ll[2]], [ur[0], ur[1], ll[2]]), sx * sy),
        "xmin": (_box_space([ll[0], ll[1], ll[2]], [ll[0], ur[1], ur[2]]), sy * sz),
        "xmax": (_box_space([ur[0], ll[1], ll[2]], [ur[0], ur[1], ur[2]]), sy * sz),
        "ymin": (_box_space([ll[0], ll[1], ll[2]], [ur[0], ll[1], ur[2]]), sx * sz),
        "ymax": (_box_space([ll[0], ur[1], ll[2]], [ur[0], ur[1], ur[2]]), sx * sz),
    }
    return [sh]


def build_sphere(name, s: Spec, k: float) -> list[Shape]:
    import openmc
    R = s.num("radius", positive=True) * k
    ctr = [v * k for v in s.vec("center", 3, [0.0, 0.0, 0.0])]
    sph = openmc.Sphere(x0=ctr[0], y0=ctr[1], z0=ctr[2], r=R)
    c = np.asarray(ctr)
    sh = Shape(name, "sphere", -sph, (c - R, c + R), volume=4.0 / 3.0 * math.pi * R ** 3,
               z_range=(c[2] - R, c[2] + R), volume_parts=[(_sphere_volume_space(R, ctr), 1.0)])
    sh.faces = {"surface": (_sphere_surface_space(R, ctr), 4.0 * math.pi * R * R)}
    return [sh]


def _frustum_region(r_bot, r_top, z0, z1, cx, cy, p0, p1):
    """Region inside a frustum between planes p0 (z0) and p1 (z1)."""
    import openmc
    if abs(r_bot - r_top) < 1e-12:
        if r_bot <= 0:
            return None, None
        cylr = openmc.ZCylinder(x0=cx, y0=cy, r=r_bot)
        return -cylr & +p0 & -p1, cylr
    slope = (r_top - r_bot) / (z1 - z0)
    zapex = z0 - r_bot / slope
    cone = openmc.ZCone(x0=cx, y0=cy, z0=zapex, r2=slope * slope)
    return -cone & +p0 & -p1, cone


def build_cone(name, s: Spec, k: float) -> list[Shape]:
    import openmc
    c = _Ctx(s, k)
    rb, rt = c.L("radius_bottom", nonneg=True), c.L("radius_top", nonneg=True)
    H, z0 = c.L("height", positive=True), c.L("z", 0.0)
    cx, cy = c.center2()
    if rb <= 0 and rt <= 0:
        raise ConfigError(f"object '{name}': cone needs a positive radius")
    p0, p1 = openmc.ZPlane(z0=z0), openmc.ZPlane(z0=z0 + H)
    reg, _ = _frustum_region(rb, rt, z0, z0 + H, cx, cy, p0, p1)
    R = max(rb, rt)
    vol = math.pi * H * (rb * rb + rb * rt + rt * rt) / 3.0
    return [Shape(name, "cone", reg, _bbox_cyl(R, z0, z0 + H, cx, cy), volume=vol, z_range=(z0, z0 + H))]


def build_polycone(name, s: Spec, k: float) -> list[Shape]:
    """Axisymmetric solid defined like Geant4's G4Polycone:
    ``planes: [[z, r_inner, r_outer], ...]`` with linear (conical) variation
    between consecutive planes. ``z`` offsets all plane z values."""
    import openmc
    c = _Ctx(s, k)
    z0 = c.L("z", 0.0)
    cx, cy = c.center2()
    raw = s.get("planes")
    try:
        pts = [(float(p[0]) * k + z0, float(p[1]) * k, float(p[2]) * k) for p in raw]
    except Exception:
        raise ConfigError(f"object '{name}': planes must be a list of [z, r_inner, r_outer]") from None
    if len(pts) < 2:
        raise ConfigError(f"object '{name}': polycone needs at least two planes")
    for (za, _, _), (zb_, _, _) in zip(pts, pts[1:]):
        if not zb_ >= za:
            raise ConfigError(f"object '{name}': polycone z values must be non-decreasing")
    planes = {}
    for z, _, _ in pts:
        planes.setdefault(round(z, 12), openmc.ZPlane(z0=z))
    pieces, vol, esa, all_cyl = [], 0.0, [], True
    for (za, ria, roa), (zb_, rib, rob) in zip(pts, pts[1:]):
        if zb_ - za <= 0:
            continue
        if roa < ria or rob < rib:
            raise ConfigError(f"object '{name}': r_outer < r_inner in polycone")
        p0, p1 = planes[round(za, 12)], planes[round(zb_, 12)]
        outer, _ = _frustum_region(roa, rob, za, zb_, cx, cy, p0, p1)
        if outer is None:
            continue
        reg = outer
        if ria > 0 or rib > 0:
            inner, _ = _frustum_region(ria, rib, za, zb_, cx, cy, p0, p1)
            if inner is not None:
                reg = reg & ~inner
        pieces.append(reg)
        h = zb_ - za
        vol += math.pi * h * (roa ** 2 + roa * rob + rob ** 2 - ria ** 2 - ria * rib - rib ** 2) / 3.0
        if abs(roa - rob) > 1e-12 or abs(ria - rib) > 1e-12:
            all_cyl = False
        else:
            esa.append((ria, roa, za, zb_))
    R = max(max(p[2] for p in pts), 1e-9)
    sh = Shape(name, "polycone", openmc.Union(pieces) if len(pieces) > 1 else pieces[0],
               _bbox_cyl(R, pts[0][0], pts[-1][0], cx, cy), volume=vol, z_range=(pts[0][0], pts[-1][0]))
    if all_cyl and cx == 0 and cy == 0:
        sh.esa = esa
    return [sh]


def build_polygon_prism(name, s: Spec, k: float) -> list[Shape]:
    import openmc
    c = _Ctx(s, k)
    n = int(s.num("sides"))
    R = c.L("radius", positive=True)  # circumradius
    H, z0 = c.L("height", positive=True), c.L("z", 0.0)
    cx, cy = c.center2()
    phase = math.radians(s.num("rotation_deg", 0.0))
    rad = Radial(R, cx, cy, n, phase)
    p0, p1 = openmc.ZPlane(z0=z0), openmc.ZPlane(z0=z0 + H)
    return [Shape(name, "polygon_prism", rad.inside() & +p0 & -p1, _bbox_cyl(R, z0, z0 + H, cx, cy),
                  volume=rad.area * H, z_range=(z0, z0 + H))]


def build_csg(name, s: Spec, k: float) -> list[Shape]:
    import openmc
    op = s.choice("op", ("union", "intersection", "difference"))
    subs = s.get("shapes")
    if not isinstance(subs, list) or len(subs) < 2:
        raise ConfigError(f"object '{name}': csg needs a list of at least two 'shapes'")
    built = []
    for i, sub in enumerate(subs):
        ss = Spec(sub, f"objects.{name}.shapes[{i}]")
        kind = ss.get("shape")
        shp = _PRIMITIVES.get(kind)
        if shp is None:
            raise ConfigError(f"object '{name}': csg sub-shape '{kind}' not supported "
                              f"(use {sorted(_PRIMITIVES)})")
        res = shp(f"{name}#{i}", ss, k)
        ss.done()
        built.append(res[0])
    regs = [b.region for b in built]
    if op == "union":
        reg = openmc.Union(regs)
        ll = np.min([b.bbox[0] for b in built], axis=0)
        ur = np.max([b.bbox[1] for b in built], axis=0)
    elif op == "intersection":
        reg = openmc.Intersection(regs)
        ll = np.max([b.bbox[0] for b in built], axis=0)
        ur = np.min([b.bbox[1] for b in built], axis=0)
    else:
        reg = regs[0]
        for r in regs[1:]:
            reg = reg & ~r
        ll, ur = built[0].bbox
    return [Shape(name, "csg", reg, (np.asarray(ll), np.asarray(ur)), volume=None,
                  z_range=(float(ll[2]), float(ur[2])))]


_PRIMITIVES: dict[str, Callable] = {
    "cylinder": build_cylinder,
    "tube": build_tube,
    "box": build_box,
    "sphere": build_sphere,
    "cone": build_cone,
    "polycone": build_polycone,
    "polygon_prism": build_polygon_prism,
}

CSG_KINDS = sorted(list(_PRIMITIVES) + ["cup", "marinelli", "csg"])


def build_csg_object(name: str, s: Spec, k: float) -> list[Shape]:
    """Build one configured object (may return several shapes: a cup and its fill)."""
    kind = s.get("shape")
    if kind == "cup":
        shapes = build_cup(name, s, k, f"{name}_fill")
    elif kind == "marinelli":
        shapes = build_marinelli(name, s, k, f"{name}_fill")
    elif kind == "csg":
        shapes = build_csg(name, s, k)
    elif kind in _PRIMITIVES:
        shapes = _PRIMITIVES[kind](name, s, k)
    else:
        raise ConfigError(f"object '{name}': unknown shape '{kind}'. CSG shapes: {CSG_KINDS}; "
                          f"CAD: 'step' (exact CSG via GEOUNED) or 'stl' (DAGMC)")
    if s.has("rotate"):
        ang = s.vec("rotate", 3)
        pivot = [v * k for v in s.vec("pivot", 3, [0.0, 0.0, 0.0])]
        memo: dict = {}
        for sh in shapes:
            sh.region = sh.region.rotate(ang, pivot=pivot, memo=memo)
            sh.bbox = _rotated_bbox(sh.bbox, ang, pivot)
            sh.esa = None
            sh.volume_parts = None
            sh.faces = {}
            sh.z_range = (float(sh.bbox[0][2]), float(sh.bbox[1][2]))
    return shapes


def _rotation_matrix(ang_deg):
    phi, theta, psi = (math.radians(a) for a in ang_deg)
    # same convention as openmc (rotation about x, then y, then z)
    cx, sx, cy, sy, cz, sz = (math.cos(phi), math.sin(phi), math.cos(theta), math.sin(theta),
                              math.cos(psi), math.sin(psi))
    Rx = np.array([[1, 0, 0], [0, cx, -sx], [0, sx, cx]])
    Ry = np.array([[cy, 0, sy], [0, 1, 0], [-sy, 0, cy]])
    Rz = np.array([[cz, -sz, 0], [sz, cz, 0], [0, 0, 1]])
    return Rz @ Ry @ Rx


def _rotated_bbox(bbox, ang, pivot):
    ll, ur = (np.asarray(b, float) for b in bbox)
    corners = np.array([[x, y, z] for x in (ll[0], ur[0]) for y in (ll[1], ur[1]) for z in (ll[2], ur[2])])
    R = _rotation_matrix(ang)
    p = np.asarray(pivot, float)
    rc = (corners - p) @ R.T + p
    return rc.min(axis=0), rc.max(axis=0)
