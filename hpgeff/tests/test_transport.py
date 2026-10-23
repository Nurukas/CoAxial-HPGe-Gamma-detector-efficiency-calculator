"""Short OpenMC transport tests (skipped without the openmc executable and photon data)."""
import os
import shutil

import pytest

pytestmark = pytest.mark.skipif(
    shutil.which("openmc") is None or not os.environ.get("OPENMC_CROSS_SECTIONS"),
    reason="needs the openmc executable and OPENMC_CROSS_SECTIONS")

CFG = {
    "units": "mm",
    "detector": {
        "crystal": {"radius": 29.5, "length": 72.5, "core": {"radius": 5.35, "depth": 72.5},
                    "dead_layers": {"side": 0.7}},
        "layers": [{"name": "endcap", "material": "aluminium", "top": 1.5, "side": 1.0,
                    "gap_top": 5.0, "gap_side": 5.0}],
    },
    "energies_keV": [661.66],
    "setups": {
        "point": {"objects": [], "source": {"type": "point", "position": [0, 0, 50]}},
        "cup": {"objects": [{"name": "cup", "shape": "cup", "material": "pmma", "inner_radius": 32,
                             "wall": 1, "bottom": 1, "inner_height": 50,
                             "fill": {"name": "water", "material": "water", "height": 40}}],
                "source": {"type": "volume", "object": "water"}},
    },
    "simulation": {"particles": 40000, "batches": 4},
}


def test_energy_bookkeeping_and_bounds(tmp_path):
    from hpgeff.config import parse_config
    from hpgeff.run import run_setup
    cfg = parse_config(CFG, tmp_path)
    cfg.output_dir = tmp_path / "out"
    for name in ("point", "cup"):
        rows = run_setup(cfg, name, log=lambda *a: None)
        r = rows[0]
        assert r["above_peak_counts"] == 0          # nothing deposits more than E
        assert r["lost_particles"] == 0
        assert 0 < r["fep_eff"] < r["total_eff"] < 1
        assert abs(r["fep_unc_batch"] / r["fep_unc"] - 1) < 0.9


def test_angular_bias_is_unbiased(tmp_path):
    """Cone biasing must reproduce the analog result (point source 50 mm above the end cap)."""
    import copy
    from hpgeff.config import parse_config
    from hpgeff.run import run_setup
    d = copy.deepcopy(CFG)
    d["simulation"]["particles"] = 200000
    a = parse_config(d, tmp_path)
    a.output_dir = tmp_path / "a"
    d["simulation"]["angular_bias"] = "auto"
    d["simulation"]["seed"] = 2
    b = parse_config(d, tmp_path)
    b.output_dir = tmp_path / "b"
    ra = run_setup(a, "point", log=lambda *x: None)[0]
    rb = run_setup(b, "point", log=lambda *x: None)[0]
    assert rb["bias_factor"] < 0.3
    diff = rb["fep_eff"] - ra["fep_eff"]
    assert abs(diff) < 4 * (ra["fep_unc"] ** 2 + rb["fep_unc"] ** 2) ** 0.5
