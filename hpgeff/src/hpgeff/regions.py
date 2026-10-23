"""Vectorised point-in-region evaluation for OpenMC CSG regions.

Used for overlap checks, volume checks and for building point-cloud sources
when rejection sampling inside OpenMC would accept too few points. It walks the
same region tree OpenMC transports with, so it checks the model actually run.
"""
from __future__ import annotations

import numpy as np


def eval_surface(surf, X, Y, Z):
    import openmc
    if isinstance(surf, openmc.XPlane):
        return X - surf.x0
    if isinstance(surf, openmc.YPlane):
        return Y - surf.y0
    if isinstance(surf, openmc.ZPlane):
        return Z - surf.z0
    if isinstance(surf, openmc.Plane):
        return surf.a * X + surf.b * Y + surf.c * Z - surf.d
    if isinstance(surf, openmc.ZCylinder):
        return (X - surf.x0) ** 2 + (Y - surf.y0) ** 2 - surf.r ** 2
    if isinstance(surf, openmc.XCylinder):
        return (Y - surf.y0) ** 2 + (Z - surf.z0) ** 2 - surf.r ** 2
    if isinstance(surf, openmc.YCylinder):
        return (X - surf.x0) ** 2 + (Z - surf.z0) ** 2 - surf.r ** 2
    if isinstance(surf, openmc.Sphere):
        return (X - surf.x0) ** 2 + (Y - surf.y0) ** 2 + (Z - surf.z0) ** 2 - surf.r ** 2
    if isinstance(surf, openmc.ZCone):
        return (X - surf.x0) ** 2 + (Y - surf.y0) ** 2 - surf.r2 * (Z - surf.z0) ** 2
    if isinstance(surf, openmc.XCone):
        return (Y - surf.y0) ** 2 + (Z - surf.z0) ** 2 - surf.r2 * (X - surf.x0) ** 2
    if isinstance(surf, openmc.YCone):
        return (X - surf.x0) ** 2 + (Z - surf.z0) ** 2 - surf.r2 * (Y - surf.y0) ** 2
    if isinstance(surf, (openmc.ZTorus, openmc.XTorus, openmc.YTorus)):
        dx, dy, dz = X - surf.x0, Y - surf.y0, Z - surf.z0
        if isinstance(surf, openmc.ZTorus):
            axial, rad = dz, np.sqrt(dx * dx + dy * dy)
        elif isinstance(surf, openmc.XTorus):
            axial, rad = dx, np.sqrt(dy * dy + dz * dz)
        else:
            axial, rad = dy, np.sqrt(dx * dx + dz * dz)
        return axial ** 2 / surf.b ** 2 + (rad - surf.a) ** 2 / surf.c ** 2 - 1.0
    # generic quadric (Quadric, Cylinder, Cone, rotated surfaces)
    if hasattr(surf, "_get_base_coeffs"):
        co = surf._get_base_coeffs()
        if len(co) == 4:
            a, b, c, d = co
            return a * X + b * Y + c * Z - d
        A, B, C, D, E, F, G, H, J, K = co
        return (A * X * X + B * Y * Y + C * Z * Z + D * X * Y + E * Y * Z + F * X * Z
                + G * X + H * Y + J * Z + K)
    raise TypeError(f"cannot evaluate surface type {type(surf).__name__}")


def contains(region, pts: np.ndarray) -> np.ndarray:
    """Boolean mask: which points (N,3) lie in ``region``."""
    import openmc
    X, Y, Z = pts[:, 0], pts[:, 1], pts[:, 2]
    cache: dict[int, np.ndarray] = {}

    def f(surf):
        key = id(surf)
        if key not in cache:
            cache[key] = eval_surface(surf, X, Y, Z)
        return cache[key]

    def rec(node):
        if isinstance(node, openmc.Halfspace):
            v = f(node.surface)
            return v > 0 if node.side == "+" else v < 0
        if isinstance(node, openmc.Intersection):
            out = np.ones(len(X), bool)
            for n in node:
                out &= rec(n)
            return out
        if isinstance(node, openmc.Union):
            out = np.zeros(len(X), bool)
            for n in node:
                out |= rec(n)
            return out
        if isinstance(node, openmc.Complement):
            return ~rec(node.node)
        raise TypeError(f"unsupported region node {type(node).__name__}")

    if region is None:
        return np.ones(len(X), bool)
    return rec(region)


def sample_box(ll, ur, n, rng) -> np.ndarray:
    ll, ur = np.asarray(ll, float), np.asarray(ur, float)
    return ll + (ur - ll) * rng.random((int(n), 3))


def mc_volume(region, ll, ur, n=400_000, rng=None):
    """Monte Carlo volume (cm3) and 1-sigma uncertainty of a region inside a box."""
    rng = rng or np.random.default_rng(12345)
    pts = sample_box(ll, ur, n, rng)
    m = contains(region, pts)
    vbox = float(np.prod(np.asarray(ur) - np.asarray(ll)))
    p = m.mean()
    return vbox * p, vbox * np.sqrt(max(p * (1 - p), 1e-300) / n)


def rejection_points(region, ll, ur, n, rng, max_tries=200):
    """Uniform points inside ``region`` (rejection in a box), shape (n, 3)."""
    out, got = [], 0
    for _ in range(max_tries):
        pts = sample_box(ll, ur, max(4 * (n - got), 10_000), rng)
        pts = pts[contains(region, pts)]
        out.append(pts)
        got += len(pts)
        if got >= n:
            break
    if got < n:
        raise RuntimeError("region too small for rejection sampling in its bounding box")
    return np.concatenate(out)[:n]
