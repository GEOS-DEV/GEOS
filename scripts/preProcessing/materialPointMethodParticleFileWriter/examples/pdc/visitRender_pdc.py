#!/usr/bin/env python3
"""VisIt renderer for the PDC cutter example.

The generic examples renderer is deliberately simple.  This PDC example needs a
three-dimensional view and two particle-region pseudocolor plots because the
cutter and substrate fields are stored under different ParticleRegion paths.

Saved view copied from the interactive VisIt GUI:
  viewNormal    = (-0.340733810949161, 0.4135483524293768, -0.8443211693392261)
  focus         = (20.0, 14.82962941761171, 5.0)
  viewUp        = (0.1648063063473661, 0.910428673249851, 0.3794186504544214)
  parallelScale = 26.769
  near/far      = -53.538, 53.538
  imagePan      = (0.05201177520211786, -0.02212855592370033)
  perspective   = on
"""
from __future__ import annotations

import argparse
import os
import re
import sys
from pathlib import Path

PDC_RENDERER_VERSION = 1

PDC_VIEW = {
    "viewNormal": (-0.340733810949161, 0.4135483524293768, -0.8443211693392261),
    "focus": (20.0, 14.82962941761171, 5.0),
    "viewUp": (0.1648063063473661, 0.910428673249851, 0.3794186504544214),
    "viewAngle": 30.0,
    "parallelScale": 26.769,
    "nearPlane": -53.538,
    "farPlane": 53.538,
    "imagePan": (0.05201177520211786, -0.02212855592370033),
    "imageZoom": 1.0,
    "perspective": 1,
}


def parse_args(argv=None):
    p = argparse.ArgumentParser(description="Render PDC cutter states with the stored 3D VisIt view")
    p.add_argument("--run-dir", default=".")
    p.add_argument("--output-dir", default=None)
    p.add_argument("--case-name", default="pdc")
    p.add_argument("--states", default="initial,final")
    p.add_argument("--width", type=int, default=1600)
    p.add_argument("--height", type=int, default=1000)
    p.add_argument("--database", default=None)
    p.add_argument("--mesh", default="CellRegion1")
    p.add_argument("--colortable", default="hot_desaturated")
    p.add_argument("--list-databases", action="store_true")
    p.add_argument("--no-annotations", action="store_true")
    p.add_argument("--dry-run", action="store_true")
    args, unknown = p.parse_known_args(sys.argv[1:] if argv is None else argv)
    if unknown:
        print("Ignoring unknown renderer arguments: " + " ".join(unknown))
    return args


def expand_path(text, base=None):
    value = os.path.expandvars(os.path.expanduser(str(text)))
    suffix = ""
    if value.endswith(" database"):
        value = value[:-len(" database")]
        suffix = " database"
    if base is not None and not os.path.isabs(value) and ":" not in value:
        value = str(Path(base) / value)
    return value + suffix


def natural_key(path):
    name = Path(path).name
    return [int(p) if p.isdigit() else p for p in re.split(r"(\d+)", name)]


def call_visit(name, *args):
    fn = globals().get(name)
    if fn is None:
        try:
            import __main__
            fn = getattr(__main__, name)
        except Exception:
            fn = None
    if fn is None:
        raise RuntimeError("VisIt function %s is not available; run with visit -cli -s" % name)
    return fn(*args)


def discover_databases(run_dir):
    candidates = []
    silo = run_dir / "siloFiles"
    patterns = []
    if silo.is_dir():
        patterns.extend([silo / "mpm_cpdi_*", silo / "mpm_*", silo / "*.visit", silo / "*.silo"])
    patterns.extend([run_dir / "mpm_cpdi_*", run_dir / "*.visit", run_dir / "*.silo"])
    for pat in patterns:
        matches = sorted((p for p in pat.parent.glob(pat.name) if p.is_file() and "restart" not in p.name.lower()), key=natural_key)
        if not matches:
            continue
        if len(matches) > 1 and "*" in pat.name:
            db = str(pat) + " database"
        else:
            db = str(matches[0])
        candidates.append((str(pat.relative_to(run_dir)) if pat.is_absolute() else str(pat), db, matches))
    return candidates


def write_visit_file(run_dir, label, matches):
    path = run_dir / ("pdc_" + re.sub(r"[^A-Za-z0-9_.-]+", "_", label).strip("_") + ".visit")
    with path.open("w") as f:
        for m in matches:
            try:
                f.write(str(m.relative_to(run_dir)) + "\n")
            except ValueError:
                f.write(str(m) + "\n")
    return str(path)


def open_database(run_dir, explicit, list_databases=False):
    if explicit:
        database = expand_path(explicit, run_dir)
        print("Opening database: " + database)
        if not call_visit("OpenDatabase", database):
            raise RuntimeError("OpenDatabase failed for " + database)
        return database
    candidates = discover_databases(run_dir)
    if list_databases:
        print("Discovered VisIt database candidates:")
        for label, db, matches in candidates:
            print("  %s: %s" % (label, db))
            for m in matches[:8]:
                print("    - " + str(m))
    errors = []
    for label, db, matches in candidates:
        attempts = [db]
        if db.endswith(" database") and len(matches) > 1:
            attempts.append(write_visit_file(run_dir, label, matches))
        for candidate in attempts:
            print("Opening database: " + candidate)
            try:
                result = call_visit("OpenDatabase", candidate)
            except Exception as exc:
                print("OpenDatabase failed: %s" % exc)
                result = 0
            if result:
                print("Matched %d files in selected database family." % len(matches))
                return candidate
            errors.append(candidate)
    raise RuntimeError("Could not open any PDC plot database. Tried: " + ", ".join(errors))


def safe_metadata_name(value):
    try:
        return str(value.name)
    except Exception:
        try:
            return str(value)
        except Exception:
            return ""


def metadata_names(md, count_method, get_method):
    out = []
    n_fn = getattr(md, count_method, None)
    get_fn = getattr(md, get_method, None)
    if not n_fn or not get_fn:
        return out
    try:
        n = int(n_fn())
    except Exception:
        return out
    for i in range(n):
        try:
            name = safe_metadata_name(get_fn(i))
            if name:
                out.append(name)
        except Exception:
            pass
    return out


def read_metadata(database):
    try:
        md = call_visit("GetMetaData", database)
    except Exception:
        md = call_visit("GetMetaData")
    scalars = metadata_names(md, "GetNumScalars", "GetScalars")
    meshes = metadata_names(md, "GetNumMeshes", "GetMeshes")
    return scalars, meshes


def base(region):
    return "ParticleRegion%d_ParticleDomains_ParticleFields" % region


def first_existing(scalars, names):
    scal = set(scalars)
    for name in names:
        if name in scal:
            return name
    low = {s.lower(): s for s in scalars}
    for name in names:
        if name.lower() in low:
            return low[name.lower()]
    return None


def define_scalar_expression(name, expression):
    try:
        call_visit("DefineScalarExpression", name, expression)
        print("Defined scalar expression %s = %s" % (name, expression))
        return name
    except Exception as exc:
        print("Could not define scalar expression %s: %s" % (name, exc))
        return None


def plastic_strain_var(scalars, region):
    b = base(region)
    direct = first_existing(scalars, [
        b + "/particlePlasticStrain_magnitude",
        b + "/particlePlasticStrainMagnitude",
        b + "/plasticStrainMagnitude",
    ])
    if direct:
        return direct
    comps = []
    for suffix in ["0", "1", "2", "3", "4", "5"]:
        name = b + "/particlePlasticStrain_" + suffix
        if name in scalars:
            comps.append(name)
    if not comps:
        # Some model-specific fields include material-name prefixes.
        for s in scalars:
            if s.startswith(b + "/") and "plasticstrain" in s.lower() and re.search(r"_[0-5]$", s):
                comps.append(s)
        comps = sorted(comps)[:6]
    if comps:
        expr = "sqrt(" + " + ".join("<%s>*<%s>" % (c, c) for c in comps) + ")"
        return define_scalar_expression("pdc_region%d_particlePlasticStrainMagnitude" % region, expr)
    return None


def material_type_var(scalars, region):
    b = base(region)
    return first_existing(scalars, [b + "/particleMaterialType", b + "/particleGroup", b + "/particleColor", b + "/particleRank"])


def add_mesh(meshes, requested):
    mesh = requested
    if mesh == "auto":
        mesh = "CellRegion1" if "CellRegion1" in meshes else (meshes[0] if meshes else "")
    if not mesh or mesh.lower() == "none":
        return False
    try:
        call_visit("AddPlot", "Mesh", mesh)
        print("Added mesh plot: " + mesh)
        return True
    except Exception as exc:
        print("Could not add mesh plot %s: %s" % (mesh, exc))
        return False


def set_pseudocolor_attributes(colortable, range_mode):
    PseudocolorAttributes = globals().get("PseudocolorAttributes")
    SetPlotOptions = globals().get("SetPlotOptions")
    if PseudocolorAttributes is None or SetPlotOptions is None:
        import __main__
        PseudocolorAttributes = getattr(__main__, "PseudocolorAttributes")
        SetPlotOptions = getattr(__main__, "SetPlotOptions")
    attrs = PseudocolorAttributes()
    if hasattr(attrs, "colorTableName"):
        attrs.colorTableName = colortable
    if range_mode == "unit":
        if hasattr(attrs, "minFlag"):
            attrs.minFlag = 1
        if hasattr(attrs, "maxFlag"):
            attrs.maxFlag = 1
        if hasattr(attrs, "min"):
            attrs.min = 0.0
        if hasattr(attrs, "max"):
            attrs.max = 1.0
    if hasattr(attrs, "pointSizePixels"):
        attrs.pointSizePixels = 4
    SetPlotOptions(attrs)


def add_pseudocolor(varname, colortable, range_mode):
    if not varname:
        return False
    try:
        call_visit("AddPlot", "Pseudocolor", varname)
        set_pseudocolor_attributes(colortable, range_mode)
        print("Added pseudocolor plot: " + varname)
        return True
    except Exception as exc:
        print("Could not add pseudocolor plot %s: %s" % (varname, exc))
        return False


def configure_view():
    try:
        view = call_visit("GetView3D")
    except Exception as exc:
        print("Could not get 3D view: %s" % exc)
        return
    for key, value in PDC_VIEW.items():
        if hasattr(view, key):
            setattr(view, key, value)
    try:
        call_visit("SetView3D", view)
        print("Applied PDC saved 3D view: viewNormal=%s focus=%s viewUp=%s parallelScale=%s perspective=%s" % (
            PDC_VIEW["viewNormal"], PDC_VIEW["focus"], PDC_VIEW["viewUp"], PDC_VIEW["parallelScale"], PDC_VIEW["perspective"]))
    except Exception as exc:
        print("Could not set 3D view: %s" % exc)


def configure_annotations(enabled=True):
    AnnotationAttributes = globals().get("AnnotationAttributes")
    SetAnnotationAttributes = globals().get("SetAnnotationAttributes")
    if AnnotationAttributes is None or SetAnnotationAttributes is None:
        import __main__
        AnnotationAttributes = getattr(__main__, "AnnotationAttributes")
        SetAnnotationAttributes = getattr(__main__, "SetAnnotationAttributes")
    annot = AnnotationAttributes()
    if not enabled:
        for attr in ["userInfoFlag", "databaseInfoFlag", "legendInfoFlag", "timeInfoFlag"]:
            if hasattr(annot, attr):
                setattr(annot, attr, 0)
    else:
        if hasattr(annot, "userInfoFlag"):
            annot.userInfoFlag = 1
        if hasattr(annot, "databaseInfoFlag"):
            annot.databaseInfoFlag = 1
        if hasattr(annot, "legendInfoFlag"):
            annot.legendInfoFlag = 1
        if hasattr(annot, "timeInfoFlag"):
            annot.timeInfoFlag = 1
    SetAnnotationAttributes(annot)


def add_time_slider(enabled=True):
    if not enabled:
        return
    try:
        slider = call_visit("CreateAnnotationObject", "TimeSlider")
        if hasattr(slider, "height"):
            slider.height = 0.07
        if hasattr(slider, "position"):
            slider.position = (0.02, 0.02)
        print("Added TimeSlider annotation")
    except Exception as exc:
        print("Could not add TimeSlider annotation: %s" % exc)


def save_window(output_dir, filename, width, height):
    SaveWindowAttributes = globals().get("SaveWindowAttributes")
    SetSaveWindowAttributes = globals().get("SetSaveWindowAttributes")
    SaveWindow = globals().get("SaveWindow")
    if SaveWindowAttributes is None or SetSaveWindowAttributes is None or SaveWindow is None:
        import __main__
        SaveWindowAttributes = getattr(__main__, "SaveWindowAttributes")
        SetSaveWindowAttributes = getattr(__main__, "SetSaveWindowAttributes")
        SaveWindow = getattr(__main__, "SaveWindow")
    attrs = SaveWindowAttributes()
    attrs.outputDirectory = str(output_dir)
    attrs.fileName = filename
    attrs.family = 0
    attrs.format = attrs.PNG
    attrs.width = width
    attrs.height = height
    if hasattr(attrs, "screenCapture"):
        attrs.screenCapture = 0
    SetSaveWindowAttributes(attrs)
    SaveWindow()
    print("Saved " + str(output_dir / (filename + ".png")))


def resolve_states(spec, n_states):
    out = []
    for raw in [s.strip() for s in spec.split(",") if s.strip()]:
        key = raw.lower()
        if key == "initial":
            out.append(("initial", 0))
        elif key == "final":
            out.append(("final", max(n_states - 1, 0)))
        else:
            idx = int(raw)
            if idx < 0:
                idx = max(n_states + idx, 0)
            out.append(("state%d" % idx, min(idx, max(n_states - 1, 0))))
    return out or [("initial", 0), ("final", max(n_states - 1, 0))]


def make_scene(label, scalars, meshes, args):
    try:
        call_visit("DeleteAllPlots")
    except Exception:
        pass
    add_mesh(meshes, args.mesh)
    if label == "initial":
        variables = [material_type_var(scalars, 1), material_type_var(scalars, 2)]
        filename = args.case_name + "_initial_particleMaterialType_regions"
        range_mode = "auto"
    else:
        variables = [plastic_strain_var(scalars, 1), plastic_strain_var(scalars, 2)]
        filename = args.case_name + "_final_particlePlasticStrainMagnitude_regions"
        range_mode = "auto"
    added = 0
    for v in variables:
        if add_pseudocolor(v, args.colortable, range_mode):
            added += 1
    if added == 0:
        raise RuntimeError("Could not add any PDC pseudocolor plots for state " + label)
    call_visit("DrawPlots")
    configure_view()
    configure_annotations(not args.no_annotations)
    add_time_slider(not args.no_annotations)
    return filename


def main(argv=None):
    args = parse_args(argv)
    run_dir = Path(expand_path(args.run_dir)).resolve()
    output_dir = Path(expand_path(args.output_dir, run_dir)).resolve() if args.output_dir else run_dir / "visit_frames"
    output_dir.mkdir(parents=True, exist_ok=True)
    print("Run directory: " + str(run_dir))
    print("Output directory: " + str(output_dir))
    database = open_database(run_dir, args.database, args.list_databases)
    scalars, meshes = read_metadata(database)
    print("Scalar variables: " + (str(scalars) if scalars else "<none reported>"))
    print("Meshes: " + (str(meshes) if meshes else "<none reported>"))
    if args.dry_run:
        return 0
    try:
        n_states = int(call_visit("TimeSliderGetNStates"))
    except Exception:
        n_states = 1
    states = resolve_states(args.states, n_states)
    print("Rendering PDC states: " + str(states))
    for label, state in states:
        try:
            call_visit("SetTimeSliderState", state)
        except Exception:
            pass
        filename = make_scene(label, scalars, meshes, args)
        save_window(output_dir, filename, args.width, args.height)
    try:
        call_visit("DeleteAllPlots")
    except Exception:
        pass
    try:
        call_visit("CloseDatabase", database)
    except Exception as exc:
        print("CloseDatabase warning: " + str(exc))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
