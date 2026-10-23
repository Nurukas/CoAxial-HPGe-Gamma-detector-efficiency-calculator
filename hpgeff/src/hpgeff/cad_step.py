"""STEP CAD -> exact CSG through GEOUNED (no triangles).

GEOUNED (UNED / Fusion for Energy, used for ITER models) decomposes each CAD
solid bounded by planes, cylinders, cones, spheres and tori into half-spaces of
those *exact* surfaces and writes OpenMC CSG. Solids with B-spline/NURBS faces
cannot be converted exactly and are rejected (spline_surfaces='stop') rather
than approximated.

GEOUNED needs FreeCAD, which usually lives in its own environment. Tell hpgeff
which Python can ``import FreeCAD, geouned`` with ``converter.python`` (or the
``HPGEFF_FREECAD_PYTHON`` environment variable); the conversion then runs in a
subprocess and its result is cached.

Guards applied to every converted solid:

* all boundary conditions in the GEOUNED output are reset to transmission
  (GEOUNED marks its last surface 'vacuum', which would silently kill photons
  crossing an ordinary surface of your container);
* the volume of each CSG cell is measured by Monte Carlo and must match the
  exact B-rep volume from the CAD kernel, otherwise the conversion is rejected;
* mass_g sets the density from the exact CAD volume.
"""
from __future__ import annotations

import hashlib
import json
import os
import subprocess
import sys
import textwrap
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np

from .regions import mc_volume
from .spec import ConfigError, Spec

_RUNNER = textwrap.dedent("""
    import json, sys, os
    import geouned
    step, out_dir, tol = sys.argv[1], sys.argv[2], sys.argv[3]
    os.chdir(out_dir)
    s = geouned.Settings(voidGen=False, debug=False, compSolids=False)
    geo = geouned.CadToCsg(settings=s)
    geo.load_step_file(filename=step, spline_surfaces='stop')
    geo.start()
    geo.export_csg(title='hpgeff', geometryName='geouned_csg', outFormat=('openmc_xml',),
                   cellSummaryFile=False, volCARD=False)
    info = []
    for m in geo.meta_list:
        if m.__id__ is None or m.NullCell or m.Volume is None:
            continue
        bb = m.BoundBox
        info.append(dict(id=int(m.__id__), label=str(m.Comments).strip(), volume_cm3=m.Volume / 1000.0,
                         bbox_cm=[bb.XMin / 10, bb.YMin / 10, bb.ZMin / 10,
                                  bb.XMax / 10, bb.YMax / 10, bb.ZMax / 10]))
    json.dump(info, open('solids.json', 'w'), indent=1)
    print('HPGEFF_OK', len(info))
""")


def _freecad_python(python: str | None) -> list[str]:
    py = python or os.environ.get("HPGEFF_FREECAD_PYTHON")
    if py:
        return [str(Path(py).expanduser())]
    try:  # FreeCAD importable in this interpreter?
        import FreeCAD  # noqa: F401
        import geouned  # noqa: F401
        return [sys.executable]
    except Exception:
        raise ConfigError(
            "STEP import needs GEOUNED + FreeCAD. Install them (e.g. `conda create -n freecad -c "
            "conda-forge freecad` then `pip install geouned` in that env) and set "
            "converter.python: /path/to/that/python (or HPGEFF_FREECAD_PYTHON).") from None


def convert_step(step: Path, cache_root: Path, python: str | None = None) -> tuple[Path, list[dict]]:
    """Run GEOUNED (cached). Returns (geometry xml, solids info)."""
    step = Path(step).resolve()
    if not step.is_file():
        raise ConfigError(f"STEP file not found: {step}")
    h = hashlib.sha256(step.read_bytes()).hexdigest()[:16]
    out = cache_root / f"geouned_{step.stem}_{h}"
    xml, info_file = out / "geouned_csg.xml", out / "solids.json"
    if xml.is_file() and info_file.is_file():
        return xml, json.loads(info_file.read_text())
    out.mkdir(parents=True, exist_ok=True)
    runner = out / "run_geouned.py"
    runner.write_text(_RUNNER)
    cmd = _freecad_python(python) + [str(runner), str(step), str(out), "0"]
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    (out / "geouned.log").write_text(proc.stdout)
    if proc.returncode != 0 or "HPGEFF_OK" not in proc.stdout:
        tail = "\n".join(proc.stdout.strip().splitlines()[-20:])
        hint = ""
        if "spline" in proc.stdout.lower():
            hint = ("\nThe CAD contains spline (NURBS) faces, which have no exact CSG form. Rebuild "
                    "those features with planes/cylinders/cones/spheres/tori, or use the STL route.")
        raise ConfigError(f"GEOUNED conversion of {step.name} failed:\n{tail}{hint}")
    return xml, json.loads(info_file.read_text())


def load_geouned_xml(xml: Path):
    """Parse GEOUNED's OpenMC XML into {cell_id: region}; boundaries reset."""
    import openmc
    root = ET.parse(xml).getroot()
    surfs = {}
    import warnings
    for el in root.findall("surface"):
        el.attrib.pop("boundary", None)  # GEOUNED tags its last surface 'vacuum'
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            srf = openmc.Surface.from_xml_element(el)
        srf.id = None  # fresh ID: GEOUNED's IDs would collide with the rest of the model
        srf.boundary_type = "transmission"
        surfs[int(el.get("id"))] = srf
    cells = {}
    for el in root.findall("cell"):
        expr = el.get("region")
        if not expr:
            continue
        cells[int(el.get("id"))] = openmc.Region.from_expression(expr, surfs)
    return cells


def build_step_objects(oname: str, s: Spec, cfg, lib):
    import openmc
    from .model import BuiltObject, _density_from_mass
    path = cfg.resolve_path(s.get("file"))
    solids_spec = s.get("solids")
    placement = Spec(s.get("placement", {}) or {}, f"objects.{oname}.placement")
    conv = Spec(s.get("converter", {}) or {}, f"objects.{oname}.converter")
    tol = float(s.num("volume_tolerance", 0.005, positive=True))
    s.done()
    python = conv.get("python", None)
    conv.done()
    xml, info = convert_step(path, cfg.output_dir / "_cad_cache", python)
    regions = load_geouned_xml(xml)
    info = [i for i in info if i["id"] in regions]
    if not isinstance(solids_spec, list) or len(solids_spec) != len(info):
        listing = "\n".join(f"  [{n}] '{i['label']}'  volume {i['volume_cm3']:.4g} cm3  "
                            f"bbox(cm) {np.round(i['bbox_cm'], 3).tolist()}" for n, i in enumerate(info))
        raise ConfigError(f"object '{oname}': the STEP file has {len(info)} solid(s); give one "
                          f"'solids' entry per solid, in this order:\n{listing}")

    k = cfg.scale
    rot = placement.vec("rotate", 3, [0.0, 0.0, 0.0])
    pivot = [v * k for v in placement.vec("pivot", 3, [0.0, 0.0, 0.0])]
    trans = np.asarray(placement.vec("translate", 3, [0.0, 0.0, 0.0])) * k
    align = placement.choice("align", ("none", "bottom_on_endcap"), "none")
    placement.done()

    from .shapes import _rotated_bbox
    bbs = []
    for i in info:
        b = np.asarray(i["bbox_cm"], float)
        bb = (b[:3], b[3:])
        if any(rot):
            bb = _rotated_bbox(bb, rot, pivot)
        bbs.append(bb)
    zmin_all = min(float(b[0][2]) for b in bbs)
    shift = trans.copy()
    if align == "bottom_on_endcap":
        shift[2] += -zmin_all
    memo: dict = {}
    out = []
    for n, (i, raw) in enumerate(zip(info, solids_spec)):
        ss = Spec(raw, f"objects.{oname}.solids[{n}]")
        name = str(ss.get("name", f"{oname}_{n}"))
        mname = ss.get("material")
        mass = ss.get("mass_g", None)
        dens = ss.get("density", None)
        ss.done()
        reg = regions[i["id"]]
        if any(rot):
            reg = reg.rotate(rot, pivot=pivot, memo=memo)
        if np.any(shift != 0):
            reg = reg.translate(shift, memo=memo)
        ll, ur = bbs[n][0] + shift, bbs[n][1] + shift
        pad = 1e-3 + 1e-3 * np.max(ur - ll)
        vmc, dv = mc_volume(reg, ll - pad, ur + pad, n=2_000_000)
        vcad = float(i["volume_cm3"])
        rel = abs(vmc - vcad) / vcad
        if rel > tol + 4 * dv / vcad:
            raise ConfigError(
                f"STEP solid '{i['label']}' ({name}): CSG volume {vmc:.5g} ± {dv:.2g} cm3 differs from "
                f"the CAD volume {vcad:.5g} cm3 by {100 * rel:.2f}% — the conversion is not faithful. "
                f"Simplify the solid or check it in FreeCAD.")
        md = lib.resolve(mname)
        notes = [f"STEP solid '{i['label']}' -> exact CSG; volume check {vmc:.5g} vs CAD {vcad:.5g} cm3"]
        rho = None
        if mass is not None:
            if md is None:
                raise ConfigError(f"object '{name}': mass given for vacuum")
            rho = _density_from_mass(name, float(mass), vcad, 0.0, md.density, notes)
        elif dens is not None:
            rho = float(dens)
        mat = lib.make_openmc(mname, name, density=rho)
        cell = openmc.Cell(name=f"object {name} (STEP)", fill=mat, region=reg)
        out.append(BuiltObject(name=name, kind="step", material_name=mname, material=mat, cell=cell,
                               region=reg, bbox=(ll, ur), volume=vcad, volume_source="cad",
                               mass_g=None if mass is None else float(mass),
                               density=mat.density if mat is not None else None, notes=notes))
    return out
