#!/usr/bin/env python3
"""Case-aware VisIt renderer for GEOS-MPM example cases.

This script is copied into the staged Lustre run directory and executed with
``visit -nowin -cli -s`` by examples/mpm_example_postprocess.py.  It prefers the
GEOS-MPM plot family ``siloFiles/mpm_cpdi_*`` and ignores restart/root files.

The scene is selected from --case-name / script filename:
  borehole_Ghareb: pressure -> plastic-strain magnitude, xy view
  collidingDisks: velocity magnitude -> velocity magnitude, xy view
  elasticDisk: velocity magnitude -> density, xy view
  pdc: velocity magnitude -> plastic-strain magnitude, saved 3-D view
"""
from __future__ import annotations

import argparse
import os
import re
import sys
from pathlib import Path

EXAMPLE_VISIT_RENDERER_VERSION = 3

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
    p = argparse.ArgumentParser(description="Render selected GEOS-MPM example states")
    p.add_argument("--run-dir", default=".")
    p.add_argument("--output-dir", default=None)
    p.add_argument("--case-name", default=None)
    p.add_argument("--states", default="initial,final")
    p.add_argument("--mesh", default="CellRegion1")
    p.add_argument("--database", default=None)
    p.add_argument("--width", type=int, default=1600)
    p.add_argument("--height", type=int, default=1100)
    p.add_argument("--colortable", default="hot_desaturated")
    p.add_argument("--list-databases", action="store_true")
    p.add_argument("--no-annotations", action="store_true")
    p.add_argument("--dry-run", action="store_true")
    args, unknown = p.parse_known_args(sys.argv[1:] if argv is None else argv)
    if unknown:
        print("Ignoring unknown renderer arguments: " + " ".join(unknown))
    return args


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


def natural_key(path):
    name = Path(path).name
    return [int(p) if p.isdigit() else p for p in re.split(r"(\d+)", name)]


def expand_path(text, base=None):
    value = os.path.expandvars(os.path.expanduser(str(text)))
    suffix = ""
    if value.endswith(" database"):
        value = value[:-len(" database")]
        suffix = " database"
    if base is not None and not os.path.isabs(value) and ":" not in value:
        value = str(Path(base) / value)
    return value + suffix


def discover_databases(run_dir):
    candidates = []
    silo = run_dir / "siloFiles"
    patterns = []
    if silo.is_dir():
        patterns.extend([silo / "mpm_cpdi_*", silo / "mpm_*", silo / "*.visit", silo / "*.silo"])
    patterns.extend([run_dir / "mpm_cpdi_*", run_dir / "*.visit", run_dir / "*.silo", run_dir / "*.pvd", run_dir / "*.vtk"])
    for pat in patterns:
        try:
            matches = sorted(
                (p for p in pat.parent.glob(pat.name)
                 if p.is_file() and "restart" not in p.name.lower() and not p.name.lower().endswith(".root")),
                key=natural_key,
            )
        except Exception:
            matches = []
        if not matches:
            continue
        if "*" in pat.name and len(matches) > 1:
            db = str(pat) + " database"
        else:
            db = str(matches[0])
        try:
            label = str(pat.relative_to(run_dir))
        except Exception:
            label = str(pat)
        candidates.append((label, db, matches))
    return candidates


def write_visit_file(run_dir, label, matches):
    safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", label).strip("_") or "database"
    path = run_dir / (run_dir.name + "_" + safe + ".visit")
    with path.open("w") as f:
        for m in matches:
            try:
                f.write(str(m.relative_to(run_dir)) + "\n")
            except ValueError:
                f.write(str(m) + "\n")
    return str(path)


def open_database(run_dir, explicit=None, list_databases=False):
    if explicit:
        database = expand_path(explicit, run_dir)
        print("Opening explicit database: " + database)
        if not call_visit("OpenDatabase", database):
            raise RuntimeError("OpenDatabase failed for " + database)
        return database
    candidates = discover_databases(run_dir)
    if list_databases:
        print("Discovered VisIt database candidates:")
        for label, db, matches in candidates:
            print("  %s: %s" % (label, db))
            for m in matches[:10]:
                print("    - " + str(m))
            if len(matches) > 10:
                print("    - ... %d more" % (len(matches) - 10))
    errors = []
    for label, db, matches in candidates:
        attempts = [db]
        if db.endswith(" database") and len(matches) > 1:
            attempts.append(write_visit_file(run_dir, label, matches))
        for candidate in attempts:
            print("Opening database: " + candidate)
            try:
                ok = call_visit("OpenDatabase", candidate)
            except Exception as exc:
                print("OpenDatabase raised: %s" % exc)
                ok = 0
            if ok:
                print("Matched %d files in selected database family." % len(matches))
                return candidate
            errors.append(candidate)
    raise RuntimeError("No GEOS-MPM VisIt database could be opened. Tried: " + ", ".join(errors))


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
    vectors = metadata_names(md, "GetNumVectors", "GetVectors")
    meshes = metadata_names(md, "GetNumMeshes", "GetMeshes")
    print("Scalar variables: %d" % len(scalars))
    print("Vector variables: %d" % len(vectors))
    print("Meshes: " + (", ".join(meshes) if meshes else "<none>"))
    return scalars, vectors, meshes


def region_numbers(names):
    nums = set()
    for name in names:
        m = re.search(r"ParticleRegion(\d+)_", name)
        if m:
            nums.add(int(m.group(1)))
    return sorted(nums) or [1]


def base(region):
    return "ParticleRegion%d_ParticleDomains_ParticleFields" % region


def first_existing(names, candidates):
    exact = set(names)
    for c in candidates:
        if c in exact:
            return c
    lowered = {n.lower(): n for n in names}
    for c in candidates:
        if c.lower() in lowered:
            return lowered[c.lower()]
    return None


def define_scalar_expression(name, expression):
    try:
        call_visit("DefineScalarExpression", name, expression)
        print("Defined scalar expression %s = %s" % (name, expression))
        return name
    except Exception as exc:
        print("WARNING: could not define scalar expression %s: %s" % (name, exc))
        return None


def component_label_sets(stem, nmax):
    """Preferred grouped VisIt scalar component labels for native Silo fields."""
    low = str(stem).lower()
    numeric = [str(i) for i in range(nmax)]
    labels = []
    if any(token in low for token in ("stress", "strain", "tensor")):
        if nmax == 9:
            labels.append(["xx", "xy", "xz", "yx", "yy", "yz", "zx", "zy", "zz"])
        else:
            labels.append(["xx", "yy", "zz", "yz", "xz", "xy"][:nmax])
    elif nmax == 3:
        labels.append(["x", "y", "z"])
    if numeric not in labels:
        labels.append(numeric)
    return labels


def find_component_name(names, candidate):
    exact = set(names)
    if candidate in exact:
        return candidate
    lowered = {n.lower(): n for n in names}
    return lowered.get(candidate.lower())


def natural_string_key(text):
    return [int(part) if part.isdigit() else part for part in re.split(r"(\d+)", str(text))]


def component_sort_key(name):
    text = str(name)
    leaf = text.rsplit("/", 1)[-1].lower()
    component_order = {
        "x": 0, "y": 1, "z": 2,
        "xx": 0, "yy": 1, "zz": 2, "yz": 3, "xz": 4, "xy": 5,
        "yx": 3, "zx": 6, "zy": 7,
    }
    if leaf in component_order:
        return (0, component_order[leaf], natural_string_key(text))
    match = re.search(r"_0_([0-9]+)$", text) or re.search(r"_([0-9]+)$", text)
    if match:
        return (1, int(match.group(1)), natural_string_key(text))
    if leaf.isdigit():
        return (2, int(leaf), natural_string_key(text))
    return (9, 9999, natural_string_key(text))


def component_series(names, prefix, stem, nmax):
    """Return components in physical order for both grouped and legacy layouts."""
    for labels in component_label_sets(stem, nmax):
        candidates = [prefix + "/" + stem + "/" + label for label in labels]
        found = [find_component_name(names, candidate) for candidate in candidates]
        if all(found):
            return found
    flat = [prefix + "/" + stem + "_" + str(i) for i in range(nmax)]
    found = [find_component_name(names, candidate) for candidate in flat]
    if all(found):
        return found

    out = []
    for labels in component_label_sets(stem, nmax):
        for label in labels:
            candidate = find_component_name(names, prefix + "/" + stem + "/" + label)
            if candidate and candidate not in out:
                out.append(candidate)
    for candidate in flat:
        found = find_component_name(names, candidate)
        if found and found not in out:
            out.append(found)
    return out


def material_component_series(names, prefix, token, nmax):
    # Match one material-specific field family such as quartz_stress/xx or quartz_stress_0_0.
    groups = {}
    allowed_labels = set()
    label_order = {}
    for labels in component_label_sets(token, nmax):
        for i, label in enumerate(labels):
            allowed_labels.add(label.lower())
            label_order.setdefault(label.lower(), i)

    def add_component(base_name, order, full_name):
        groups.setdefault(base_name, []).append((order, full_name))

    for name in names:
        if not name.startswith(prefix + "/"):
            continue
        leaf = name[len(prefix) + 1:]
        if token.lower() not in leaf.lower():
            continue
        if "/" in leaf:
            base_name, grouped_label = leaf.rsplit("/", 1)
            grouped_label = grouped_label.lower()
            if grouped_label in allowed_labels:
                add_component(base_name, label_order.get(grouped_label, nmax), name)
                continue
        match = re.search(r"_0_([0-9]+)$", leaf) or re.search(r"_([0-9]+)$", leaf)
        if match:
            order = int(match.group(1))
            base_name = leaf[:match.start()]
            add_component(base_name, order, name)

    if not groups:
        return []
    best_base, best_components = sorted(
        groups.items(),
        key=lambda item: (-len(item[1]), natural_string_key(item[0])),
    )[0]
    return [name for _, name in sorted(best_components, key=lambda item: (item[0], natural_string_key(item[1])))[:nmax]]


def velocity_variable(scalars, vectors, region):
    p = base(region)
    direct = first_existing(scalars + vectors, [p + "/particleVelocity_magnitude", p + "/particleVelocityMagnitude"])
    if direct:
        return direct
    comps = component_series(scalars, p, "particleVelocity", 3)
    if len(comps) == 3:
        expr = "sqrt(" + " + ".join("<%s>*<%s>" % (c, c) for c in comps) + ")"
        return define_scalar_expression("region%d_velocityMagnitude" % region, expr)
    # Some databases expose velocity as a vector.
    vec = first_existing(vectors, [p + "/particleVelocity", p + "/velocity"])
    if vec:
        return define_scalar_expression("region%d_velocityMagnitude" % region, "magnitude(<%s>)" % vec)
    return None


def plastic_variable(scalars, region):
    p = base(region)
    direct = first_existing(scalars, [
        p + "/particlePlasticStrain_magnitude",
        p + "/particlePlasticStrainMagnitude",
        p + "/plasticStrainMagnitude",
    ])
    if direct:
        return direct
    comps = component_series(scalars, p, "particlePlasticStrain", 6)
    if not comps:
        comps = material_component_series(scalars, p, "plasticStrain", 6)
    if comps:
        expr = "sqrt(" + " + ".join("<%s>*<%s>" % (c, c) for c in comps) + ")"
        return define_scalar_expression("region%d_plasticStrainMagnitude" % region, expr)
    return None


def pressure_variable(scalars, region):
    p = base(region)
    direct = first_existing(scalars, [p + "/particlePressure", p + "/pressure", p + "/Pressure"])
    if direct:
        return direct
    comps = component_series(scalars, p, "particleStress", 3)
    if len(comps) >= 3:
        expr = "-(<%s> + <%s> + <%s>)/3.0" % (comps[0], comps[1], comps[2])
        return define_scalar_expression("region%d_pressure" % region, expr)
    comps = material_component_series(scalars, p, "stress", 6)
    if len(comps) >= 3:
        expr = "-(<%s> + <%s> + <%s>)/3.0" % (comps[0], comps[1], comps[2])
        return define_scalar_expression("region%d_pressure" % region, expr)
    # Fall back to cell/material stress if particles do not carry stress.
    cell = first_existing(scalars, [
        "CellRegion1_Solid_MaterialFields/stress_11",
        "CellRegion1_Solid_MaterialFields/stress_22",
        "CellRegion1_Solid_MaterialFields/stress_33",
    ])
    if cell:
        sx = first_existing(scalars, ["CellRegion1_Solid_MaterialFields/stress_11"])
        sy = first_existing(scalars, ["CellRegion1_Solid_MaterialFields/stress_22"])
        sz = first_existing(scalars, ["CellRegion1_Solid_MaterialFields/stress_33"])
        if sx and sy and sz:
            return define_scalar_expression("cell_pressure", "-(<%s> + <%s> + <%s>)/3.0" % (sx, sy, sz))
    return None


def density_variable(scalars, region):
    p = base(region)
    direct = first_existing(scalars, [
        p + "/particleDensity",
        p + "/density",
        p + "/particleMassDensity",
    ])
    if direct:
        return direct
    # Prefer particle material density over cell density when present.
    for name in scalars:
        if name.startswith(p + "/") and "density" in name.lower() and "null" not in name.lower():
            return name
    return first_existing(scalars, [
        "CellRegion1_Solid_MaterialFields/density",
        "CellRegion1_ElementFields/density",
    ])


def variable_for_key(key, scalars, vectors, region):
    k = key.lower()
    if k in ("velocity", "velocitymagnitude", "speed"):
        return velocity_variable(scalars, vectors, region)
    if k in ("plastic", "plasticstrain", "plasticstrainmagnitude"):
        return plastic_variable(scalars, region)
    if k in ("pressure", "meanstress"):
        return pressure_variable(scalars, region)
    if k in ("density", "rho"):
        return density_variable(scalars, region)
    return first_existing(scalars + vectors, [key])


def add_mesh(meshes, requested):
    if not requested or requested.lower() == "none":
        return
    mesh = requested
    if mesh == "auto":
        mesh = "CellRegion1" if "CellRegion1" in meshes else (meshes[0] if meshes else "")
    if not mesh:
        return
    try:
        call_visit("AddPlot", "Mesh", mesh)
        print("Added mesh overlay: " + mesh)
    except Exception as exc:
        print("WARNING: could not add mesh %s: %s" % (mesh, exc))


def set_pc_atts(colortable, unit_range=False):
    PseudocolorAttributes = globals().get("PseudocolorAttributes")
    if PseudocolorAttributes is None:
        return
    try:
        a = PseudocolorAttributes()
        if hasattr(a, "colorTableName"):
            a.colorTableName = colortable
        if hasattr(a, "minFlag"):
            a.minFlag = 1 if unit_range else 0
        if hasattr(a, "maxFlag"):
            a.maxFlag = 1 if unit_range else 0
        if unit_range:
            if hasattr(a, "min"):
                a.min = 0.0
            if hasattr(a, "max"):
                a.max = 1.0
        if hasattr(a, "pointType"):
            try:
                a.pointType = a.Point
            except Exception:
                pass
        if hasattr(a, "pointSizePixels"):
            a.pointSizePixels = 4
        call_visit("SetPlotOptions", a)
    except Exception as exc:
        print("WARNING: could not set PseudocolorAttributes: %s" % exc)


def add_pseudocolors(var_key, scalars, vectors, regions, colortable, prefer_all_regions=False):
    added = []
    # For density, pressure in single-region examples, use first available region/cell field.
    selected_regions = regions if prefer_all_regions else regions[:]
    for region in selected_regions:
        var = variable_for_key(var_key, scalars, vectors, region)
        if not var:
            if not prefer_all_regions:
                continue
            print("No %s variable found for ParticleRegion%d" % (var_key, region))
            continue
        try:
            call_visit("AddPlot", "Pseudocolor", var)
            set_pc_atts(colortable, unit_range=False)
            added.append(var)
            print("Added pseudocolor: " + var)
        except Exception as exc:
            print("WARNING: could not add pseudocolor %s: %s" % (var, exc))
        if not prefer_all_regions and added:
            break
    if not added and var_key.lower() == "density":
        var = density_variable(scalars, regions[0] if regions else 1)
        if var:
            call_visit("AddPlot", "Pseudocolor", var)
            set_pc_atts(colortable, unit_range=False)
            added.append(var)
    if not added:
        raise RuntimeError("No variable found for requested scene key %s" % var_key)
    return added


def configure_annotations(enabled=True):
    try:
        AnnotationAttributes = globals().get("AnnotationAttributes")
        if AnnotationAttributes is not None:
            a = AnnotationAttributes()
            a.userInfoFlag = 1 if enabled else 0
            a.databaseInfoFlag = 1 if enabled else 0
            if hasattr(a, "timeInfoFlag"):
                a.timeInfoFlag = 1 if enabled else 0
            call_visit("SetAnnotationAttributes", a)
    except Exception:
        pass
    if enabled:
        try:
            call_visit("CreateAnnotationObject", "TimeSlider")
        except Exception:
            pass


def configure_view(case_key):
    try:
        call_visit("ResetView")
    except Exception:
        pass
    View3DAttributes = globals().get("View3DAttributes")
    if View3DAttributes is None:
        return
    try:
        v = call_visit("GetView3D")
        if "pdc" in case_key:
            for key, value in PDC_VIEW.items():
                if hasattr(v, key):
                    setattr(v, key, value)
        else:
            v.viewNormal = (0.0, 0.0, 1.0)
            v.viewUp = (0.0, 1.0, 0.0)
            if hasattr(v, "perspective"):
                v.perspective = 0
        call_visit("SetView3D", v)
        print("Configured view for " + case_key)
    except Exception as exc:
        print("WARNING: could not configure 3D view: %s" % exc)


def scene_for_case(case_name):
    key = re.sub(r"[^A-Za-z0-9]+", "_", (case_name or "").strip()).strip("_").lower()
    if "borehole" in key:
        return {
            "case_key": "borehole",
            "initial": ("pressure", "pressure", False),
            "final": ("plasticStrainMagnitude", "plasticStrainMagnitude", False),
        }
    if "colliding" in key:
        return {
            "case_key": "collidingDisks",
            "initial": ("velocityMagnitude", "velocityMagnitude", True),
            "final": ("velocityMagnitude", "velocityMagnitude", True),
        }
    if "elasticdisk" in key or "elastic_disk" in key:
        return {
            "case_key": "elasticDisk",
            "initial": ("velocityMagnitude", "velocityMagnitude", False),
            "final": ("density", "density", False),
        }
    if key == "pdc" or "pdc" in key:
        return {
            "case_key": "pdc",
            "initial": ("velocity", "velocity", True),
            "final": ("plasticStrainMagnitude", "plasticStrainMagnitude", True),
        }
    return {
        "case_key": key or "mpm",
        "initial": ("velocityMagnitude", "velocityMagnitude", False),
        "final": ("damage", "damage", False),
    }


def num_states():
    try:
        return int(call_visit("TimeSliderGetNStates"))
    except Exception:
        try:
            return int(call_visit("GetDatabaseNStates"))
        except Exception:
            return 1


def set_state(index):
    try:
        call_visit("SetTimeSliderState", int(index))
    except Exception as exc:
        print("WARNING: could not set time slider state %s: %s" % (index, exc))


def safe_name(text):
    text = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(text)).strip("_")
    return text or "field"


def save_window(output_dir, file_prefix, width, height):
    SaveWindowAttributes = globals().get("SaveWindowAttributes")
    if SaveWindowAttributes is None:
        raise RuntimeError("SaveWindowAttributes is unavailable")
    output_dir.mkdir(parents=True, exist_ok=True)
    a = SaveWindowAttributes()
    a.outputDirectory = str(output_dir)
    a.fileName = file_prefix
    a.family = 0
    a.width = int(width)
    a.height = int(height)
    try:
        a.format = a.PNG
    except Exception:
        pass
    call_visit("SetSaveWindowAttributes", a)
    path = call_visit("SaveWindow")
    print("Saved " + str(Path(output_dir) / (file_prefix + ".png")))
    return path


def render_state(args, scene, state_label, state_index, scalars, vectors, meshes, regions):
    call_visit("DeleteAllPlots")
    key, token, all_regions = scene[state_label]
    added = add_pseudocolors(key, scalars, vectors, regions, args.colortable, prefer_all_regions=all_regions)
    add_mesh(meshes, args.mesh)
    call_visit("DrawPlots")
    configure_view(scene["case_key"])
    configure_annotations(not args.no_annotations)
    set_state(state_index)
    # Draw again after state changes so expressions update.
    try:
        call_visit("DrawPlots")
    except Exception:
        pass
    if not args.dry_run:
        prefix = safe_name(args.case_name or scene["case_key"]) + "_" + state_label + "_" + safe_name(token)
        save_window(Path(args.output_dir), prefix, args.width, args.height)
    return added


def main():
    args = parse_args()
    run_dir = Path(os.path.expanduser(os.path.expandvars(args.run_dir))).resolve()
    output_dir = Path(os.path.expanduser(os.path.expandvars(args.output_dir or str(run_dir / "visit_frames")))).resolve()
    args.output_dir = str(output_dir)
    case_name = args.case_name or Path(sys.argv[0]).stem.replace("visitRender_", "")
    args.case_name = case_name
    scene = scene_for_case(case_name)
    print("Run directory: " + str(run_dir))
    print("Output directory: " + str(output_dir))
    print("Case scene: " + scene["case_key"])
    database = open_database(run_dir, args.database, args.list_databases)
    scalars, vectors, meshes = read_metadata(database)
    regions = region_numbers(scalars + vectors)
    print("Particle regions: " + ", ".join(str(r) for r in regions))
    n = max(num_states(), 1)
    states = []
    for s in [p.strip() for p in args.states.split(",") if p.strip()]:
        if s == "initial":
            states.append(("initial", 0))
        elif s == "final":
            states.append(("final", max(n - 1, 0)))
        else:
            states.append(("state" + s, int(s)))
    # Only the configured initial/final scenes are special.  Integer state aliases use final variable.
    for label, idx in states:
        scene_label = label if label in ("initial", "final") else "final"
        print("Rendering %s at state %s" % (label, idx))
        render_state(args, scene, scene_label, idx, scalars, vectors, meshes, regions)
    try:
        call_visit("DeleteAllPlots")
        call_visit("CloseDatabase", database)
    except Exception:
        pass


if __name__ == "__main__":
    main()
