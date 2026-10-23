"""Materials: a small built-in library plus user definitions.

Every material is reduced to (density, elemental mass fractions). That single
representation feeds both OpenMC (``add_element(..., 'wo')``) and the
independent attenuation calculation in :mod:`hpgeff.esa`, so the two can never
disagree about what a material is.

User definitions (YAML, densities in g/cm3)::

    materials:
      glass:   {density: 2.43, compounds: {SiO2: 81.8, Na2O: 12.8, CaO: 5.45}}
      pmma:    {density: 1.19, formula: C5H8O2}
      brass:   {density: 8.5,  elements: {Cu: 0.63, Zn: 0.37}}          # mass fractions
      gas:     {density: 0.0018, elements: {Ar: 1}, fractions: atom}
      wetsoil: {density: 1.6,  mixture: {soil: 0.8, water: 0.2}}       # mass fractions of materials
      pmma2:   {base: pmma, density: 1.20}                               # override a density
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field

from .spec import ConfigError, Spec

# ----------------------------------------------------------------------------
# Built-in library. Compositions: NIST/PNNL compendium style (mass fractions or
# formulas); densities in g/cm3. Users should override anything they know better.
# ----------------------------------------------------------------------------
BUILTIN: dict[str, dict | None] = {
    "vacuum": None,
    "air": {"density": 0.001205,  # dry air near sea level (NIST)
            "elements": {"C": 0.000124, "N": 0.755268, "O": 0.231781, "Ar": 0.012827}},
    "germanium": {"density": 5.323, "formula": "Ge"},
    "aluminium": {"density": 2.699, "formula": "Al"},
    "beryllium": {"density": 1.848, "formula": "Be"},
    "beryllium_oxide": {"density": 3.01, "formula": "BeO"},
    "copper": {"density": 8.96, "formula": "Cu"},
    "brass": {"density": 8.47, "elements": {"Cu": 0.65, "Zn": 0.35}},
    "lead": {"density": 11.35, "formula": "Pb"},
    "tin": {"density": 7.31, "formula": "Sn"},
    "iron": {"density": 7.874, "formula": "Fe"},
    "stainless_steel_304": {"density": 8.00,
                            "elements": {"Fe": 0.695, "Cr": 0.19, "Ni": 0.095, "Mn": 0.02}},
    "pmma": {"density": 1.19, "formula": "C5H8O2"},
    "polyethylene": {"density": 0.94, "formula": "C2H4"},
    "polypropylene": {"density": 0.90, "formula": "C3H6"},
    "polystyrene": {"density": 1.06, "formula": "C8H8"},
    "ptfe": {"density": 2.20, "formula": "C2F4"},
    "mylar": {"density": 1.40, "formula": "C10H8O4"},
    "kapton": {"density": 1.42, "formula": "C22H10N2O5"},
    "water": {"density": 1.00, "formula": "H2O"},
    "quartz": {"density": 2.65, "formula": "SiO2"},
    "fused_silica": {"density": 2.20, "formula": "SiO2"},
    # A generic commercial soda-lime glass (illustrative; define your own glass).
    "soda_lime_glass": {"density": 2.50,
                        "compounds": {"SiO2": 72.6, "Na2O": 13.9, "CaO": 8.4, "MgO": 3.9,
                                      "Al2O3": 1.1, "K2O": 0.1}},
}

ALIASES = {
    "aluminum": "aluminium", "al": "aluminium", "ge": "germanium", "be": "beryllium",
    "beo": "beryllium_oxide", "cu": "copper", "pb": "lead", "acrylic": "pmma",
    "plexiglass": "pmma", "plexiglas": "pmma", "perspex": "pmma", "pe": "polyethylene",
    "hdpe": "polyethylene", "pp": "polypropylene", "ps": "polystyrene", "teflon": "ptfe",
    "pet": "mylar", "void": "vacuum", "none": "vacuum",
}


# ----------------------------------------------------------------------------
# Chemical formula parsing
# ----------------------------------------------------------------------------
_TOKEN = re.compile(r"([A-Z][a-z]?)|(\d+(?:\.\d+)?)|(\()|(\))")


def _atomic_weight(sym: str) -> float:
    import openmc.data
    try:
        return openmc.data.atomic_weight(sym)
    except Exception as exc:  # pragma: no cover - depends on openmc tables
        raise ConfigError(f"unknown element or element without natural abundance: '{sym}'") from exc


def parse_formula(formula: str) -> dict[str, float]:
    """Atom counts of a chemical formula. Supports parentheses and hydrates.

    >>> parse_formula("Ca(OH)2")  # {'Ca': 1, 'O': 2, 'H': 2}
    """
    total: dict[str, float] = {}
    for part in _split_hydrate(formula.replace(" ", "")):
        if not part:
            continue
        m = re.match(r"^(\d+(?:\.\d+)?)(.*)$", part)
        mult = 1.0
        if m and m.group(2) and m.group(2)[0].isupper():
            mult, part = float(m.group(1)), m.group(2)
        for el, n in _parse_group(part, formula).items():
            total[el] = total.get(el, 0.0) + mult * n
    if not total:
        raise ConfigError(f"could not parse chemical formula '{formula}'")
    return total


def _split_hydrate(s: str) -> list[str]:
    """Split 'CaSO4*2H2O' (or 'CaSO4·2H2O') into ['CaSO4', '2H2O'].

    A '.' is always a decimal point ('Na0.5Cl0.5'); write hydrates with '*' or '·'.
    """
    return re.split(r"[·*]", s)


def _parse_group(s: str, original: str) -> dict[str, float]:
    pos = 0
    stack: list[dict[str, float]] = [{}]
    last: dict[str, float] | None = None
    while pos < len(s):
        m = _TOKEN.match(s, pos)
        if not m:
            raise ConfigError(f"invalid character in formula '{original}' at '{s[pos:]}'")
        el, num, lp, rp = m.groups()
        pos = m.end()
        if el:
            stack[-1][el] = stack[-1].get(el, 0.0) + 1.0
            last = {el: 1.0}
        elif num:
            if last is None:
                raise ConfigError(f"number without element in formula '{original}'")
            k = float(num)
            for e, n in last.items():
                stack[-1][e] += n * (k - 1.0)
            last = None
        elif lp:
            stack.append({})
            last = None
        elif rp:
            if len(stack) == 1:
                raise ConfigError(f"unbalanced ')' in formula '{original}'")
            grp = stack.pop()
            for e, n in grp.items():
                stack[-1][e] = stack[-1].get(e, 0.0) + n
            last = grp
    if len(stack) != 1:
        raise ConfigError(f"unbalanced '(' in formula '{original}'")
    return stack[0]


def formula_mass_fractions(formula: str) -> dict[str, float]:
    atoms = parse_formula(formula)
    masses = {e: n * _atomic_weight(e) for e, n in atoms.items()}
    tot = sum(masses.values())
    return {e: m / tot for e, m in masses.items()}


def _normalise(d: dict[str, float], what: str) -> dict[str, float]:
    tot = float(sum(d.values()))
    if tot <= 0:
        raise ConfigError(f"{what}: fractions must sum to a positive number")
    if any(v < 0 for v in d.values()):
        raise ConfigError(f"{what}: negative fraction")
    return {k: float(v) / tot for k, v in d.items()}


# ----------------------------------------------------------------------------
# Material definitions
# ----------------------------------------------------------------------------
@dataclass
class MaterialDef:
    name: str
    density: float  # g/cm3
    mass_fractions: dict[str, float] = field(default_factory=dict)  # element -> w

    def scaled(self, density: float) -> "MaterialDef":
        return MaterialDef(self.name, float(density), dict(self.mass_fractions))


class MaterialLibrary:
    """Resolves material names (built-in or user-defined) to :class:`MaterialDef`."""

    def __init__(self, user: dict | None = None):
        self._raw: dict[str, dict | None | str] = {}
        for k, v in BUILTIN.items():
            self._raw[k] = v
        for k, v in (user or {}).items():
            self._raw[str(k)] = v
        self._cache: dict[str, MaterialDef | None] = {}
        self._counter = 0

    def canonical(self, name: str) -> str:
        key = str(name)
        if key in self._raw:
            return key
        low = key.lower()
        if low in self._raw:
            return low
        if low in ALIASES:
            return ALIASES[low]
        raise ConfigError(
            f"unknown material '{name}'. Define it under 'materials:' or use one of: "
            f"{sorted(self._raw)}")

    def is_vacuum(self, name: str) -> bool:
        return self.resolve(name) is None

    def resolve(self, name: str, _stack: tuple = ()) -> MaterialDef | None:
        key = self.canonical(name)
        if key in self._cache:
            return self._cache[key]
        if key in _stack:
            raise ConfigError(f"circular material definition: {' -> '.join(_stack + (key,))}")
        raw = self._raw[key]
        if raw is None or (isinstance(raw, str) and raw.lower() in ("vacuum", "void")):
            self._cache[key] = None
            return None
        if isinstance(raw, str):  # alias to another material
            md = self.resolve(raw, _stack + (key,))
            self._cache[key] = md
            return md
        md = self._build(key, Spec(raw, f"materials.{key}"), _stack + (key,))
        self._cache[key] = md
        return md

    def _build(self, key: str, s: Spec, stack: tuple) -> MaterialDef:
        base = s.get("base", None)
        density = s.num("density", None, positive=True)
        kinds = [k for k in ("formula", "elements", "compounds", "mixture") if s.has(k)]
        if base is not None:
            bdef = self.resolve(base, stack)
            if bdef is None:
                raise ConfigError(f"materials.{key}: base material '{base}' is vacuum")
            if kinds:
                raise ConfigError(f"materials.{key}: use either 'base' or a composition, not both")
            s.done()
            return MaterialDef(key, density if density is not None else bdef.density,
                               dict(bdef.mass_fractions))
        if len(kinds) != 1:
            raise ConfigError(f"materials.{key}: give exactly one of formula/elements/compounds/"
                              f"mixture (got {kinds or 'none'})")
        if density is None:
            raise ConfigError(f"materials.{key}: 'density' (g/cm3) is required")
        kind = kinds[0]
        if kind == "formula":
            w = formula_mass_fractions(str(s.get("formula")))
        elif kind == "elements":
            fr = s.choice("fractions", ("mass", "atom"), "mass")
            els = {str(k): float(v) for k, v in dict(s.get("elements")).items()}
            for e in els:
                _atomic_weight(e)
            if fr == "atom":
                els = {e: n * _atomic_weight(e) for e, n in els.items()}
            w = _normalise(els, f"materials.{key}.elements")
        elif kind == "compounds":
            comps = _normalise({str(k): float(v) for k, v in dict(s.get("compounds")).items()},
                               f"materials.{key}.compounds")
            w = {}
            for comp, frac in comps.items():
                for e, we in formula_mass_fractions(comp).items():
                    w[e] = w.get(e, 0.0) + frac * we
        else:  # mixture of other materials by mass
            mix = _normalise({str(k): float(v) for k, v in dict(s.get("mixture")).items()},
                             f"materials.{key}.mixture")
            w = {}
            for mname, frac in mix.items():
                md = self.resolve(mname, stack)
                if md is None:
                    continue  # vacuum contributes no mass
                for e, we in md.mass_fractions.items():
                    w[e] = w.get(e, 0.0) + frac * we
            w = _normalise(w, f"materials.{key}.mixture")
        s.done()
        return MaterialDef(key, density, w)

    # ------------------------------------------------------------------
    def make_openmc(self, name: str, instance_name: str | None = None,
                    density: float | None = None):
        """Return a new ``openmc.Material`` (or None for vacuum).

        Each call creates a distinct material so that source constraints by
        material can single out one object.
        """
        import openmc
        md = self.resolve(name)
        if md is None:
            return None
        if density is not None:
            md = md.scaled(density)
        mat = openmc.Material(name=instance_name or md.name)
        for el, w in sorted(md.mass_fractions.items()):
            if w > 0:
                mat.add_element(el, w, percent_type="wo")
        mat.set_density("g/cm3", md.density)
        mat._hpgeff_def = md  # used by the independent attenuation calculation
        return mat


def element_set(defs) -> set[str]:
    out: set[str] = set()
    for d in defs:
        if d is not None:
            out |= set(d.mass_fractions)
    return out
