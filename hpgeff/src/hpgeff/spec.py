"""Small strict-dictionary helper used to parse the YAML configuration.

Every section is read through :class:`Spec`, which remembers which keys were
consumed. Calling :meth:`Spec.done` raises if keys are left over, so a typo such
as ``radiu: 32`` is an error instead of being silently ignored (a silently
ignored dimension would quietly change an efficiency).
"""
from __future__ import annotations

from typing import Any


class ConfigError(ValueError):
    """Raised for invalid or inconsistent configuration input."""


_MISSING = object()


class Spec:
    def __init__(self, data: Any, path: str = ""):
        if data is None:
            data = {}
        if not isinstance(data, dict):
            raise ConfigError(f"'{path or '<root>'}' must be a mapping, got {type(data).__name__}")
        self._data = dict(data)
        self._used: set[str] = set()
        self.path = path

    # -- access -------------------------------------------------------------
    def _key(self, key: str) -> str:
        return f"{self.path}.{key}" if self.path else key

    def has(self, key: str) -> bool:
        return key in self._data

    def get(self, key: str, default: Any = _MISSING) -> Any:
        if key in self._data:
            self._used.add(key)
            return self._data[key]
        if default is _MISSING:
            raise ConfigError(f"missing required key '{self._key(key)}'")
        return default

    def num(self, key: str, default: Any = _MISSING, *, positive: bool = False,
            nonneg: bool = False) -> float:
        val = self.get(key, default)
        if val is None:
            return None  # type: ignore[return-value]
        try:
            x = float(val)
        except (TypeError, ValueError):
            raise ConfigError(f"'{self._key(key)}' must be a number, got {val!r}") from None
        if positive and not x > 0:
            raise ConfigError(f"'{self._key(key)}' must be > 0, got {x}")
        if nonneg and x < 0:
            raise ConfigError(f"'{self._key(key)}' must be >= 0, got {x}")
        return x

    def vec(self, key: str, n: int, default: Any = _MISSING) -> list[float]:
        val = self.get(key, default)
        if val is None:
            return None  # type: ignore[return-value]
        try:
            out = [float(v) for v in val]
        except (TypeError, ValueError):
            raise ConfigError(f"'{self._key(key)}' must be a list of {n} numbers, got {val!r}") from None
        if len(out) != n:
            raise ConfigError(f"'{self._key(key)}' must have {n} values, got {len(out)}")
        return out

    def sub(self, key: str, default: Any = _MISSING) -> "Spec":
        val = self.get(key, default)
        return Spec(val, self._key(key))

    def choice(self, key: str, choices, default: Any = _MISSING):
        val = self.get(key, default)
        if val not in choices:
            raise ConfigError(f"'{self._key(key)}' must be one of {list(choices)}, got {val!r}")
        return val

    def raw(self) -> dict:
        return dict(self._data)

    def keys(self):
        return self._data.keys()

    def done(self) -> None:
        extra = set(self._data) - self._used
        if extra:
            where = self.path or "<root>"
            raise ConfigError(
                f"unknown key(s) in '{where}': {sorted(extra)} — check spelling "
                f"(unknown keys are rejected so typos cannot silently change a result)"
            )


def set_by_path(data: dict, dotted: str, value: Any) -> None:
    """Set ``data['a']['b'][0]['c'] = value`` from ``'a.b.0.c'`` (creates dicts)."""
    parts = dotted.split(".")
    cur: Any = data
    for i, part in enumerate(parts[:-1]):
        nxt = parts[i + 1]
        if isinstance(cur, list):
            cur = cur[int(part)]
            continue
        if part not in cur or cur[part] is None:
            cur[part] = [] if nxt.isdigit() else {}
        cur = cur[part]
    last = parts[-1]
    if isinstance(cur, list):
        cur[int(last)] = value
    else:
        cur[last] = value
