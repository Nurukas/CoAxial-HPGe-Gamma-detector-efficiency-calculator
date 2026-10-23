"""Configuration loading (YAML) with strict validation and command-line overrides."""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import yaml

from .spec import ConfigError, Spec, set_by_path

LENGTH_UNITS = {"mm": 0.1, "cm": 1.0, "m": 100.0}


@dataclass
class SimSettings:
    particles: int = 1_000_000          # histories per energy (fixed-length runs)
    batches: int = 20
    threads: int | None = None
    seed: int = 1
    peak_window_keV: float = 0.001      # half-width of the full-energy window
    photon_cutoff_keV: float = 1.0
    electron_treatment: str = "ttb"
    rel_unc_target: float | None = None  # e.g. 0.005 -> run until 0.5 % (OpenMC trigger)
    max_batches: int = 2000
    angular_bias: str = "off"           # off | auto  (see docs: exact for unscattered photons)
    max_lost_particles: int = 10
    rel_max_lost_particles: float = 1e-6
    keep_run_files: bool = False
    openmc_exec: str = "openmc"


@dataclass
class Config:
    base_dir: Path
    scale: float                         # config length unit -> cm
    unit_name: str
    materials: dict
    detector: dict
    ambient: str
    setups: dict
    energies_keV: list
    sim: SimSettings
    reference: dict | None
    output_dir: Path
    raw: dict = field(default_factory=dict)

    def resolve_path(self, p) -> Path:
        p = Path(p).expanduser()
        return p if p.is_absolute() else (self.base_dir / p)


def _parse_sim(s: Spec) -> SimSettings:
    d = SimSettings()
    d.particles = int(s.num("particles", d.particles, positive=True))
    d.batches = int(s.num("batches", d.batches, positive=True))
    thr = s.get("threads", None)
    d.threads = None if thr in (None, "auto") else int(thr)
    d.seed = int(s.num("seed", d.seed, positive=True))
    d.peak_window_keV = s.num("peak_window_keV", d.peak_window_keV, positive=True)
    d.photon_cutoff_keV = s.num("photon_cutoff_keV", d.photon_cutoff_keV, positive=True)
    d.electron_treatment = s.choice("electron_treatment", ("ttb", "led"), d.electron_treatment)
    rut = s.get("rel_unc_target", None)
    d.rel_unc_target = None if rut is None else float(rut)
    d.max_batches = int(s.num("max_batches", d.max_batches, positive=True))
    d.angular_bias = s.choice("angular_bias", ("off", "auto"), d.angular_bias)
    d.max_lost_particles = int(s.num("max_lost_particles", d.max_lost_particles, positive=True))
    d.rel_max_lost_particles = s.num("rel_max_lost_particles", d.rel_max_lost_particles, positive=True)
    d.keep_run_files = bool(s.get("keep_run_files", d.keep_run_files))
    d.openmc_exec = str(s.get("openmc_exec", d.openmc_exec))
    s.done()
    if d.particles < d.batches:
        raise ConfigError("simulation.particles must be >= simulation.batches")
    if d.peak_window_keV > 5:
        raise ConfigError("simulation.peak_window_keV is a half-width in keV; values > 5 keV "
                          "would merge Compton/escape features into the peak")
    return d


def parse_config(data: dict, base_dir: Path) -> Config:
    s = Spec(data)
    unit = s.choice("units", tuple(LENGTH_UNITS), "mm")
    materials = s.get("materials", {}) or {}
    detector = s.get("detector")
    ambient = s.get("ambient", "air")
    setups = s.get("setups")
    if not isinstance(setups, dict) or not setups:
        raise ConfigError("'setups' must map setup names to {objects, source}")
    energies = s.get("energies_keV", None)
    sim = _parse_sim(s.sub("simulation", {}))
    reference = s.get("reference", None)
    out = s.get("output_dir", "results")
    s.done()
    if energies is not None:
        energies = [float(e) for e in energies]
        if any(e <= 0 for e in energies):
            raise ConfigError("energies_keV must be positive")
    for name, st in setups.items():
        ss = Spec(st, f"setups.{name}")
        ss.get("objects", [])
        ss.get("source")
        ss.get("energies_keV", None)
        ss.get("description", None)
        ss.done()
        if st.get("energies_keV") is None and energies is None:
            raise ConfigError(f"setup '{name}': no energies (set top-level energies_keV)")
    if reference is not None:
        rs = Spec(reference, "reference")
        if rs.get("setup") not in setups:
            raise ConfigError(f"reference.setup '{reference.get('setup')}' is not one of the setups")
        rs.get("data")
        rs.get("fit", None)
        rs.done()
    cfg = Config(base_dir=base_dir, scale=LENGTH_UNITS[unit], unit_name=unit, materials=materials,
                 detector=detector, ambient=str(ambient), setups=setups,
                 energies_keV=energies or [], sim=sim, reference=reference,
                 output_dir=Path(out), raw=data)
    cfg.output_dir = cfg.resolve_path(out)
    return cfg


def load_config(path, overrides: list[str] | None = None) -> Config:
    path = Path(path)
    with open(path) as fh:
        data = yaml.safe_load(fh) or {}
    for ov in overrides or []:
        if "=" not in ov:
            raise ConfigError(f"override must look like key.path=value, got '{ov}'")
        key, val = ov.split("=", 1)
        set_by_path(data, key.strip(), yaml.safe_load(val))
    return parse_config(data, path.resolve().parent)


def setup_energies(cfg: Config, name: str) -> list[float]:
    e = cfg.setups[name].get("energies_keV")
    return [float(x) for x in (e if e is not None else cfg.energies_keV)]
