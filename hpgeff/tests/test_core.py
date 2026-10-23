"""Fast unit tests (no transport). Run: pytest -q"""
import math
import os

import numpy as np
import pytest

from hpgeff.config import parse_config
from hpgeff.detector import build_detector
from hpgeff.materials import MaterialLibrary, parse_formula
from hpgeff.regions import contains, mc_volume
from hpgeff.shapes import build_csg_object
from hpgeff.spec import ConfigError, Spec

K = 0.1  # mm -> cm


# ----------------------------------------------------------------------------- materials
def test_formula_parser():
    assert parse_formula("Ca(OH)2") == {"Ca": 1, "O": 2, "H": 2}
    assert parse_formula("CaSO4*2H2O") == {"Ca": 1, "S": 1, "O": 6, "H": 4}
    assert parse_formula("Na0.5Cl0.5") == {"Na": 0.5, "Cl": 0.5}
    assert parse_formula("Al2(SO4)3") == {"Al": 2, "S": 3, "O": 12}
    with pytest.raises(ConfigError):
        parse_formula("Ca(OH")


def test_material_definitions():
    lib = MaterialLibrary({
        "glass": {"density": 2.43, "compounds": {"SiO2": 75, "Na2O": 11.7, "CaO": 5}},
        "mix": {"density": 1.0, "mixture": {"water": 0.5, "pmma": 0.5}},
        "pmma2": {"base": "pmma", "density": 1.18},
    })
    g = lib.resolve("glass")
    assert abs(sum(g.mass_fractions.values()) - 1) < 1e-12
    # SiO2 mass fraction of Si = 28.0855/60.0843
    assert abs(g.mass_fractions["Si"] - 75 / 91.7 * 28.0855 / 60.0843) < 1e-3
    assert lib.resolve("PMMA").density == 1.19  # alias / case-insensitive
    assert lib.resolve("pmma2").density == 1.18
    assert lib.resolve("vacuum") is None
    with pytest.raises(ConfigError):
        MaterialLibrary({"bad": {"density": 1, "formula": "H2O", "typo": 1}}).resolve("bad")
    with pytest.raises(ConfigError):
        lib.resolve("unobtainium")


# ----------------------------------------------------------------------------- shapes
SHAPES = {
    "cup": {"shape": "cup", "inner_radius": 32, "wall": 1, "bottom": 1, "inner_height": 50,
            "fill": {"name": "glass", "height": 40}},
    "marinelli": {"shape": "marinelli", "inner_radius": 60, "wall": 2, "inner_height": 120,
                  "well_radius": 40, "well_depth": 80, "fill": {"name": "water", "height": 110}},
    "polycone": {"shape": "polycone", "planes": [[0, 0, 20], [2, 0, 20], [2, 18, 20], [40, 18, 25], [42, 0, 25]]},
    "cone": {"shape": "cone", "radius_bottom": 10, "radius_top": 5, "height": 30},
    "sphere": {"shape": "sphere", "radius": 12, "center": [1, 2, 20]},
    "tube": {"shape": "tube", "inner_radius": 5, "outer_radius": 9, "height": 10, "center": [3, 4]},
    "cupfacets": {"shape": "cup", "inner_radius": 32, "wall": 1, "bottom": 1, "inner_height": 50,
                  "facets": 12, "fill": {"name": "gf"}},
}


@pytest.mark.parametrize("name", list(SHAPES))
def test_shape_volumes_and_no_overlap(name):
    d = SHAPES[name]
    s = Spec(d, name)
    shapes = build_csg_object(name, s, K)
    rng = np.random.default_rng(1)
    for sh in shapes:
        ll, ur = sh.bbox
        v, dv = mc_volume(sh.region, np.asarray(ll) - 1e-3, np.asarray(ur) + 1e-3, n=1_500_000, rng=rng)
        assert abs(v - sh.volume) < 5 * dv + 1e-6, (sh.name, v, dv, sh.volume)
    if len(shapes) > 1:
        ll = np.min([x.bbox[0] for x in shapes], axis=0)
        ur = np.max([x.bbox[1] for x in shapes], axis=0)
        pts = ll + (ur - ll) * rng.random((500_000, 3))
        n_in = np.sum([contains(x.region, pts) for x in shapes], axis=0)
        assert np.all(n_in <= 1)


def test_faceted_cylinder_area():
    s = Spec({"shape": "cylinder", "radius": 10, "height": 10, "facets": 8}, "c")
    sh = build_csg_object("c", s, K)[0]
    R = 1.0
    assert abs(sh.volume - 0.5 * 8 * R * R * math.sin(2 * math.pi / 8) * 1.0) < 1e-12


def test_unknown_shape_key_rejected():
    s = Spec({"shape": "cylinder", "radius": 10, "height": 10, "heigth": 3}, "c")
    build_csg_object("c", s, K)
    with pytest.raises(ConfigError):
        s.done()


# ----------------------------------------------------------------------------- detector
DET = {
    "crystal": {"radius": 29.5, "length": 72.5, "core": {"radius": 5.35, "depth": 55},
                "dead_layers": {"front": 0.3, "side": 0.7, "back": 0.5, "core": 0.3}},
    "layers": [{"name": "holder", "material": "aluminium", "top": 0.05, "side": 1.0},
               {"name": "endcap", "material": "aluminium", "top": 1.5, "side": 1.0,
                "gap_top": 5, "gap_side": 5}],
    "contact_pin": {"radius": 3, "length": 40},
}


@pytest.mark.parametrize("variant", [{}, {"edge_rounding": 8.0}, {"core": {"radius": 5.35, "depth": 72.5}},
                                     {"core": {"radius": 5.35, "depth": 55, "bottom": "round"}}])
def test_detector_active_volume_and_partition(variant):
    d = {**DET, "crystal": {**DET["crystal"], **variant}}
    det = build_detector(Spec(d, "detector"), K, MaterialLibrary({}))
    ll, ur = det.active_bbox
    vals = [mc_volume(det.active_cell.region, ll, ur, n=1_000_000, rng=np.random.default_rng(s))[0]
            for s in range(4)]
    va = det.info["active_volume_cm3"]
    assert abs(np.mean(vals) - va) / va < 2e-3
    # cells partition the envelope exactly
    ll, ur = det.envelope_bbox
    pts = ll + (ur - ll) * np.random.default_rng(5).random((600_000, 3))
    n_in = np.sum([contains(c.region, pts) for c in det.cells], axis=0)
    inside = contains(det.envelope, pts)
    assert np.all(n_in[inside] == 1) and np.all(n_in[~inside] == 0)
    assert abs(det.info["z_front"] - (-(0.05 + 5 + 1.5) * K)) < 1e-12


def test_pin_larger_than_core_is_explained():
    d = {**DET, "contact_pin": {"radius": 7.775, "length": 40}}
    with pytest.raises(ConfigError, match="diameter"):
        build_detector(Spec(d, "detector"), K, MaterialLibrary({}))


# ----------------------------------------------------------------------------- ESA geometry
def test_ray_tube_intervals_against_brute_force():
    from hpgeff.esa import _len_before, _tube_intervals
    rng = np.random.default_rng(0)
    for trial in range(400):
        P = rng.uniform([-4, -4, -2], [4, 4, 6])
        if trial % 5 == 0:
            P[:2] = 0
        u = rng.normal(size=3)
        if trial % 7 == 0:
            u = np.array([0.0, 0.0, 1.0])
        u /= np.linalg.norm(u)
        seg = (rng.choice([0.0, rng.uniform(0, 1.5)]), rng.uniform(1.6, 3.5), rng.uniform(-3, 0), rng.uniform(0.1, 3))
        iv = _tube_intervals(P[None, :], u[None, :], *seg)
        L = _len_before(iv, np.array([1e9]))[0]
        t = np.linspace(0, 30, 300_001)
        x, y, z = P[0] + t * u[0], P[1] + t * u[1], P[2] + t * u[2]
        r = np.hypot(x, y)
        Lb = np.sum((r >= seg[0]) & (r <= seg[1]) & (z >= seg[2]) & (z <= seg[3])) * (t[1] - t[0])
        assert abs(L - Lb) < 2e-3


# ----------------------------------------------------------------------------- transfer
def test_reference_fit_and_transfer(tmp_path):
    from hpgeff.transfer import ReferenceCurve, read_reference, transfer
    E = np.array([60, 88, 122, 166, 279, 392, 662, 898, 1173, 1332, 1836.0])
    eff = 0.05 * (E / 100.0) ** -0.8 * (1 - np.exp(-(E / 60.0) ** 2))
    p = tmp_path / "ref.csv"
    p.write_text("energy_keV,efficiency,unc %\n" + "\n".join(f"{a},{b},2.0" for a, b in zip(E, eff)))
    E2, e2, u2 = read_reference(p)
    assert np.allclose(u2, 0.02 * e2)
    c = ReferenceCurve(E2, e2, u2, degree=4)
    v, _ = c([150.0, 700.0])
    true = 0.05 * (np.array([150, 700]) / 100.0) ** -0.8 * (1 - np.exp(-(np.array([150, 700]) / 60.0) ** 2))
    assert np.all(np.abs(v / true - 1) < 0.02)
    sample = [dict(energy_keV=662.0, fep_eff=0.01, fep_unc=0.0001)]
    ref = [dict(energy_keV=662.0, fep_eff=0.05, fep_unc=0.0005)]
    out = transfer(sample, ref, c)[0]
    assert abs(out["transfer_factor"] - 0.2) < 1e-12
    assert abs(out["eff_final"] - 0.2 * c([662.0])[0][0]) < 1e-15


# ----------------------------------------------------------------------------- config
def test_config_rejects_unknown_keys(tmp_path):
    base = {"detector": DET, "setups": {"a": {"objects": [], "source": {"type": "point", "position": [0, 0, 5]}}},
            "energies_keV": [661.66]}
    parse_config(base, tmp_path)
    with pytest.raises(ConfigError):
        parse_config({**base, "simulaton": {}}, tmp_path)


# ----------------------------------------------------------------------------- STL checks
def test_mesh_faceting_estimate():
    trimesh = pytest.importorskip("trimesh")
    from hpgeff.cad_stl import mesh_report
    R, H = 32.0, 50.0
    for n in (16, 32, 64, 128):
        m = trimesh.creation.cylinder(radius=R, height=H, sections=n)
        rep = mesh_report(m)
        true_def = 1 - (n / (2 * math.pi)) * math.sin(2 * math.pi / n)
        sag = R * (1 - math.cos(math.pi / n))
        assert rep["watertight"] and rep["winding_consistent"]
        assert abs(rep["volume_deficit_rel_est"] - true_def) / true_def < 0.05
        assert abs(rep["max_chord_deviation"] - sag) / sag < 0.05


@pytest.mark.skipif(not os.environ.get("OPENMC_CROSS_SECTIONS"), reason="needs photon data")
def test_attenuation_coefficients_against_xcom():
    from hpgeff.esa import mu_no_coherent
    lib = MaterialLibrary({})
    # XCOM (total minus coherent), cm2/g: Ge 662 keV 0.0694, Al 662 keV 0.0745, water 662 keV 0.0856
    for name, ref in (("germanium", 0.0694), ("aluminium", 0.0745), ("water", 0.0856)):
        md = lib.resolve(name)
        mu = mu_no_coherent(md, [661657.0])[0] / md.density
        assert abs(mu / ref - 1) < 0.01, (name, mu)


# ----------------------------------------------------------------------------- STEP import
def test_geouned_xml_is_sanitised():
    """GEOUNED marks its last surface 'vacuum' and uses IDs 1..n: both must be neutralised."""
    from pathlib import Path
    from hpgeff.cad_step import load_geouned_xml
    xml = Path(__file__).parent / "data" / "geouned_cup_sample.xml"
    assert 'boundary="vacuum"' in xml.read_text()
    cells = load_geouned_xml(xml)
    assert set(cells) == {1, 2}
    surfs = {}
    for reg in cells.values():
        surfs.update(reg.get_surfaces())
    assert all(s.boundary_type == "transmission" for s in surfs.values())
    # converted volumes equal the exact cup (13.631 cm3) and glass (160.85 cm3) volumes
    rng = np.random.default_rng(3)
    v1, d1 = mc_volume(cells[1], [-3.31, -3.31, -0.01], [3.31, 3.31, 5.11], n=1_000_000, rng=rng)
    v2, d2 = mc_volume(cells[2], [-3.21, -3.21, 0.09], [3.21, 3.21, 5.11], n=1_000_000, rng=rng)
    assert abs(v1 - math.pi * (3.3 ** 2 * 0.1 + (3.3 ** 2 - 3.2 ** 2) * 5.0)) < 5 * d1
    assert abs(v2 - math.pi * 3.2 ** 2 * 5.0) < 5 * d2


def test_pooling_of_two_runs():
    from hpgeff.run import _combine
    a = dict(fep_eff=0.10, total_eff=0.3, fep_unc=0, total_unc=0, fep_unc_batch=0.001, histories=1000,
             bias_factor=1.0, above_peak_counts=0, runtime_s=1.0, peak_to_total=0)
    b = dict(a, fep_eff=0.12, histories=3000, fep_unc_batch=0.0005)
    c = _combine(a, b)
    assert c["histories"] == 4000
    assert abs(c["fep_eff"] - (0.1 * 1000 + 0.12 * 3000) / 4000) < 1e-15
    assert abs(c["fep_unc"] - math.sqrt(c["fep_eff"] * (1 - c["fep_eff"]) / 4000)) < 1e-15
