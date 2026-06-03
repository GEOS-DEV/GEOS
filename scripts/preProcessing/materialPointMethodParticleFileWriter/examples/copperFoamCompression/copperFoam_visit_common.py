#!/usr/bin/env python3
"""Shared VisIt renderer for the copper foam compression examples."""
from __future__ import annotations

import argparse
import math
import os
import re
import sys
from pathlib import Path

COPPER_FOAM_RENDERER_VERSION = 1


def parse_args(argv, default_case_name, default_view, default_width, default_height):
    p = argparse.ArgumentParser(description="Render copper foam initial/final MPM states with VisIt")
    p.add_argument("--run-dir", default=".")
    p.add_argument("--database", default=None)
    p.add_argument("--output-dir", default=None)
    p.add_argument("--case-name", default=default_case_name)
    p.add_argument("--states", default="initial,final")
    p.add_argument("--width", type=int, default=default_width)
    p.add_argument("--height", type=int, default=default_height)
    p.add_argument("--view", choices=("auto", "xy", "xz", "yz"), default=default_view)
    p.add_argument("--view-padding", type=float, default=1.06)
    p.add_argument("--colortable", default="hot_desaturated")
    p.add_argument("--point-size-pixels", type=int, default=4)
    p.add_argument("--mesh", default="CellRegion1")
    p.add_argument("--region", default="ParticleRegion1")
    p.add_argument("--no-annotations", action="store_true")
    p.add_argument("--time-slider", dest="time_slider", action="store_true", default=True)
    p.add_argument("--no-time-slider", dest="time_slider", action="store_false")
    p.add_argument("--list-databases", action="store_true")
    p.add_argument("--list-variables", action="store_true")
    p.add_argument("--dry-run", action="store_true")
    args, unknown = p.parse_known_args(argv)
    if unknown:
        print("Ignoring unknown VisIt/script arguments: {0}".format(" ".join(unknown)))
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
        raise RuntimeError("VisIt function {0} is not available; run with visit -cli -s".format(name))
    return fn(*args)


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
    return [int(part) if part.isdigit() else part for part in re.split(r"(\d+)", name)]


def natural_string_key(text):
    return [int(part) if part.isdigit() else part for part in re.split(r"(\d+)", str(text))]


def existing_files(directory, pattern):
    try:
        return sorted([p for p in directory.glob(pattern) if p.is_file()], key=natural_key)
    except Exception:
        return []


def not_restart(path):
    lowered = str(path).lower()
    return "restart" not in lowered and not lowered.endswith(".root")


def make_candidate(label, pattern_path, matches):
    matches = [p for p in matches if not_restart(p)]
    if not matches:
        return None
    if len(matches) > 1:
        database = str(pattern_path) + " database"
        is_family = True
    else:
        database = str(matches[0])
        is_family = False
    return {"label": label, "database": database, "matches": matches, "is_family": is_family}


def discover_databases(run_dir):
    out = []
    silo = run_dir / "siloFiles"
    if silo.is_dir():
        for pattern in ("mpm_cpdi_*", "mpm_*", "*.visit", "*.silo", "*.h5", "*.hdf5"):
            matches = [p for p in existing_files(silo, pattern) if p.name != "data"]
            cand = make_candidate("siloFiles/" + pattern, silo / pattern, matches)
            if cand:
                out.append(cand)
    for pattern in ("mpm_cpdi_*", "mpm_*", "*.visit", "*.silo", "*.pvd", "*.vtk"):
        cand = make_candidate(pattern, run_dir / pattern, existing_files(run_dir, pattern))
        if cand:
            out.append(cand)
    return out


def write_visit_file(run_dir, candidate):
    safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", candidate.get("label", "database")).strip("_") or "database"
    visit_path = run_dir / (run_dir.name + "_" + safe + ".visit")
    with visit_path.open("w") as handle:
        for path in candidate.get("matches", []):
            try:
                handle.write(str(path.relative_to(run_dir)) + "\n")
            except ValueError:
                handle.write(str(path) + "\n")
    return str(visit_path)


def print_candidates(candidates):
    print("Discovered VisIt database candidates:")
    if not candidates:
        print("  <none>")
    for cand in candidates:
        print("  {0}: {1}".format(cand["label"], cand["database"]))
        for match in cand["matches"][:8]:
            print("    - {0}".format(match))
        if len(cand["matches"]) > 8:
            print("    - ... {0} more".format(len(cand["matches"]) - 8))


def try_open_database(database):
    print("Opening database: {0}".format(database))
    try:
        result = call_visit("OpenDatabase", database)
    except Exception as exc:
        print("OpenDatabase failed for {0}: {1}".format(database, exc))
        return False
    if result == 0:
        print("OpenDatabase returned failure for {0}".format(database))
        return False
    return True


def open_database(run_dir, explicit, list_databases=False):
    if explicit:
        database = expand_path(explicit, run_dir)
        if not try_open_database(database):
            raise RuntimeError("Could not open explicit database: " + database)
        return database
    candidates = discover_databases(run_dir)
    if list_databases:
        print_candidates(candidates)
    for cand in candidates:
        attempts = [cand["database"]]
        if cand.get("is_family"):
            attempts.append(write_visit_file(run_dir, cand))
        for database in attempts:
            if try_open_database(database):
                print("Selected VisIt database candidate: {0}".format(cand["label"]))
                print("Matched {0} files in this database family.".format(len(cand.get("matches", []))))
                return database
    raise FileNotFoundError("No GEOS-MPM plot database found. Expected <run-dir>/siloFiles/mpm_cpdi_*")


def metadata_name(value):
    try:
        name = getattr(value, "name")
        return str(name) if name else None
    except Exception:
        return None


def metadata_entries(md, count_method, get_method):
    out = []
    n_fn = getattr(md, count_method, None)
    get_fn = getattr(md, get_method, None)
    if n_fn is None or get_fn is None:
        return out
    try:
        n = int(n_fn())
    except Exception:
        return out
    for i in range(n):
        try:
            name = metadata_name(get_fn(i))
        except Exception:
            name = None
        if name:
            out.append(name)
    return out


def metadata_variables(database):
    try:
        md = call_visit("GetMetaData", database)
    except Exception as exc:
        print("Could not query metadata: {0}".format(exc))
        return [], [], []
    return (
        metadata_entries(md, "GetNumScalars", "GetScalars"),
        metadata_entries(md, "GetNumVectors", "GetVectors"),
        metadata_entries(md, "GetNumMeshes", "GetMeshes"),
    )


def compact(text):
    return re.sub(r"[^a-z0-9]+", "", str(text).lower())


def region_prefix(region):
    return str(region) + "_ParticleDomains_ParticleFields"


def field(region, name):
    return region_prefix(region) + "/" + str(name)


def quote_var(name):
    return "<" + str(name) + ">"


def square_sum_expression(components):
    terms = ["({0})*({0})".format(quote_var(name)) for name in components]
    return "sqrt(" + " + ".join(terms) + ")"


def define_scalar_expression(name, expression):
    print("Defining expression {0} = {1}".format(name, expression))
    call_visit("DefineScalarExpression", name, expression)
    return name


def find_direct_scalar(scalars, region, names, required_tokens=()):
    scalar_set = set(scalars)
    for name in names:
        exact = field(region, name)
        if exact in scalar_set:
            return exact
    r = compact(region)
    token_compacts = [compact(t) for t in required_tokens]
    name_compacts = [compact(n) for n in names]
    matches = []
    for scalar in scalars:
        c = compact(scalar)
        if r in c and all(t in c for t in token_compacts) and any(n in c for n in name_compacts):
            matches.append(scalar)
    return sorted(matches, key=natural_string_key)[0] if matches else None


def component_label_sets(basename, n_expected):
    if n_expected is None:
        return []
    low = str(basename).lower()
    numeric = [str(i) for i in range(n_expected)]
    labels = []
    if any(token in low for token in ("stress", "strain", "tensor")):
        if n_expected == 9:
            labels.append(["xx", "xy", "xz", "yx", "yy", "yz", "zx", "zy", "zz"])
        else:
            labels.append(["xx", "yy", "zz", "yz", "xz", "xy"][:n_expected])
    elif n_expected == 3:
        labels.append(["x", "y", "z"])
    if numeric not in labels:
        labels.append(numeric)
    return labels


def grouped_component_candidates(region, basename, n_expected):
    candidates = []
    for labels in component_label_sets(basename, n_expected):
        candidates.append([field(region, basename + "/" + label) for label in labels])
    return candidates


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


def choose_components(scalars, region, basename, n_expected=None, required_tokens=()):
    scalar_set = set(scalars)
    if n_expected is not None:
        for grouped in grouped_component_candidates(region, basename, n_expected):
            if all(name in scalar_set for name in grouped):
                return grouped
        exact = [field(region, basename + "_" + str(i)) for i in range(n_expected)]
        if all(name in scalar_set for name in exact):
            return exact
    comps = []
    r = compact(region)
    b = compact(basename)
    token_compacts = [compact(t) for t in required_tokens]
    for name in scalars:
        c = compact(name)
        if r in c and b in c and all(t in c for t in token_compacts):
            comps.append(name)
    return sorted(comps, key=component_sort_key)


def prepare_velocity_magnitude(scalars, vectors, region, expression_prefix):
    direct = find_direct_scalar(
        scalars,
        region,
        ["particleVelocity_magnitude", "particleVelocityMagnitude", "velocityMagnitude", "VelocityMagnitude"],
        required_tokens=("velocity",),
    )
    if direct:
        print("Using velocity magnitude scalar: {0}".format(direct))
        return direct
    comps = choose_components(scalars, region, "particleVelocity", n_expected=3)
    if len(comps) < 3:
        comps = choose_components(scalars, region, "velocity")[:3]
    if len(comps) >= 3:
        print("Velocity components: {0}".format(comps[:3]))
        return define_scalar_expression(expression_prefix + "_velocityMagnitude", square_sum_expression(comps[:3]))
    r = compact(region)
    for vector in vectors:
        c = compact(vector)
        if r in c and "velocity" in c:
            expr = "magnitude({0})".format(quote_var(vector))
            return define_scalar_expression(expression_prefix + "_velocityMagnitude", expr)
    raise RuntimeError("Could not find velocity components in {0}".format(region))


def prepare_plastic_strain_magnitude(scalars, region, expression_prefix):
    direct = find_direct_scalar(
        scalars,
        region,
        [
            "particlePlasticStrain_magnitude",
            "particlePlasticStrainMagnitude",
            "plasticStrainMagnitude",
            "equivalentPlasticStrain",
        ],
        required_tokens=("plastic",),
    )
    if direct:
        print("Using plastic strain magnitude scalar: {0}".format(direct))
        return direct
    comps = choose_components(scalars, region, "particlePlasticStrain", n_expected=6)
    if len(comps) < 1:
        comps = choose_components(scalars, region, "plasticStrain")
    if len(comps) < 1:
        raise RuntimeError("Could not find plastic-strain components in {0}".format(region))
    print("Plastic-strain components: {0}".format(comps))
    return define_scalar_expression(expression_prefix + "_plasticStrainMagnitude", square_sum_expression(comps))


def prepare_expressions(scalars, vectors, region, case_name):
    prefix = re.sub(r"[^A-Za-z0-9_]+", "_", case_name).strip("_") or "copperFoam"
    return {
        "velocityMagnitude": prepare_velocity_magnitude(scalars, vectors, region, prefix),
    }


def configure_pseudocolor(colortable, point_size_pixels):
    try:
        PseudocolorAttributes = globals().get("PseudocolorAttributes")
        SetPlotOptions = globals().get("SetPlotOptions")
        if PseudocolorAttributes is None or SetPlotOptions is None:
            import __main__
            PseudocolorAttributes = getattr(__main__, "PseudocolorAttributes")
            SetPlotOptions = getattr(__main__, "SetPlotOptions")
        attrs = PseudocolorAttributes()
        if colortable and hasattr(attrs, "colorTableName"):
            attrs.colorTableName = str(colortable)
        if hasattr(attrs, "pointSizePixels"):
            attrs.pointSizePixels = int(point_size_pixels)
        if hasattr(attrs, "pointSizeVarEnabled"):
            attrs.pointSizeVarEnabled = 0
        if hasattr(attrs, "legendFlag"):
            attrs.legendFlag = 1
        if hasattr(attrs, "minFlag"):
            attrs.minFlag = 0
        if hasattr(attrs, "maxFlag"):
            attrs.maxFlag = 0
        SetPlotOptions(attrs)
    except Exception as exc:
        print("Could not configure Pseudocolor attributes: {0}".format(exc))


def add_pseudocolor(variable, label, colortable, point_size_pixels):
    print("Adding {0} Pseudocolor plot: {1}".format(label, variable))
    call_visit("AddPlot", "Pseudocolor", variable)
    configure_pseudocolor(colortable, point_size_pixels)


def add_mesh(mesh_name):
    if not mesh_name or str(mesh_name).lower() in ("none", "off", "false", "0"):
        return
    try:
        print("Adding background mesh plot: {0}".format(mesh_name))
        call_visit("AddPlot", "Mesh", mesh_name)
        try:
            MeshAttributes = globals().get("MeshAttributes")
            SetPlotOptions = globals().get("SetPlotOptions")
            if MeshAttributes is None or SetPlotOptions is None:
                import __main__
                MeshAttributes = getattr(__main__, "MeshAttributes")
                SetPlotOptions = getattr(__main__, "SetPlotOptions")
            attrs = MeshAttributes()
            if hasattr(attrs, "legendFlag"):
                attrs.legendFlag = 0
            SetPlotOptions(attrs)
        except Exception:
            pass
    except Exception as exc:
        print("Could not add mesh plot {0}: {1}".format(mesh_name, exc))


def flatten_numbers(value):
    out = []
    if value is None:
        return out
    if isinstance(value, (int, float)):
        return [float(value)]
    if isinstance(value, str):
        for part in re.split(r"[,\s]+", value.strip()):
            if part:
                try:
                    out.append(float(part))
                except ValueError:
                    pass
        return out
    try:
        iterator = iter(value)
    except TypeError:
        return out
    for item in iterator:
        out.extend(flatten_numbers(item))
    return out


def query_extents():
    for query in ("SpatialExtents", "OriginalSpatialExtents"):
        try:
            call_visit("Query", query)
            value = call_visit("GetQueryOutputValue")
            nums = flatten_numbers(value)
            if len(nums) >= 6:
                ext = tuple(float(v) for v in nums[:6])
                print("{0}: {1}".format(query, ext))
                return ext
        except Exception as exc:
            print("Could not query {0}: {1}".format(query, exc))
    return None


def get_view3d():
    try:
        return call_visit("GetView3D")
    except Exception:
        View3DAttributes = globals().get("View3DAttributes")
        if View3DAttributes is None:
            import __main__
            View3DAttributes = getattr(__main__, "View3DAttributes")
        return View3DAttributes()


def view_axes(view_name):
    if view_name == "xy":
        return (0.0, 0.0, 1.0), (0.0, 1.0, 0.0), (0, 1), (2, 3), (4, 5)
    if view_name == "xz":
        return (0.0, -1.0, 0.0), (0.0, 0.0, 1.0), (0, 1), (4, 5), (2, 3)
    if view_name == "yz":
        return (1.0, 0.0, 0.0), (0.0, 0.0, 1.0), (2, 3), (4, 5), (0, 1)
    return None


def fit_view(view_name, extents, width, height, padding):
    if view_name == "auto" or extents is None:
        return
    axes = view_axes(view_name)
    if axes is None:
        return
    vals = list(float(v) for v in extents)
    for i in (0, 2, 4):
        if vals[i + 1] < vals[i]:
            vals[i], vals[i + 1] = vals[i + 1], vals[i]
    x0, x1, y0, y1, z0, z1 = vals
    center = ((x0 + x1) * 0.5, (y0 + y1) * 0.5, (z0 + z1) * 0.5)
    spans = [max(x1 - x0, 1.0e-9), max(y1 - y0, 1.0e-9), max(z1 - z0, 1.0e-9)]
    diag = math.sqrt(sum(s * s for s in spans))
    normal, up, h_pair, v_pair, _ = axes
    aspect = float(width) / float(height) if height else 1.0
    hspan = max(vals[h_pair[1]] - vals[h_pair[0]], 1.0e-9)
    vspan = max(vals[v_pair[1]] - vals[v_pair[0]], 1.0e-9)
    parallel_scale = 0.5 * max(vspan, hspan / max(aspect, 1.0e-9)) * padding

    v = get_view3d()
    v.viewNormal = normal
    v.viewUp = up
    v.focus = center
    v.parallelScale = parallel_scale
    v.nearPlane = -10.0 * max(diag, parallel_scale)
    v.farPlane = 10.0 * max(diag, parallel_scale)
    if hasattr(v, "viewAngle"):
        v.viewAngle = 30.0
    if hasattr(v, "perspective"):
        v.perspective = 0
    if hasattr(v, "imagePan"):
        v.imagePan = (0.0, 0.0)
    if hasattr(v, "imageZoom"):
        v.imageZoom = 1.0
    if hasattr(v, "centerOfRotationSet"):
        v.centerOfRotationSet = 1
    if hasattr(v, "centerOfRotation"):
        v.centerOfRotation = center
    call_visit("SetView3D", v)
    print("Configured {0} view: focus={1}, viewNormal={2}, viewUp={3}, parallelScale={4}".format(
        view_name, center, normal, up, parallel_scale
    ))


def configure_annotations(enabled):
    try:
        AnnotationAttributes = globals().get("AnnotationAttributes")
        SetAnnotationAttributes = globals().get("SetAnnotationAttributes")
        if AnnotationAttributes is None or SetAnnotationAttributes is None:
            import __main__
            AnnotationAttributes = getattr(__main__, "AnnotationAttributes")
            SetAnnotationAttributes = getattr(__main__, "SetAnnotationAttributes")
        annot = AnnotationAttributes()
        for attr, value in (("userInfoFlag", 0), ("databaseInfoFlag", 0), ("timeInfoFlag", 1 if enabled else 0), ("legendInfoFlag", 1 if enabled else 0)):
            if hasattr(annot, attr):
                setattr(annot, attr, value)
        try:
            annot.axes3D.visible = 1 if enabled else 0
        except Exception:
            pass
        try:
            annot.axes2D.visible = 1 if enabled else 0
        except Exception:
            pass
        SetAnnotationAttributes(annot)
    except Exception as exc:
        print("Could not configure annotations: {0}".format(exc))


def add_time_slider(enabled):
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
        print("Could not add TimeSlider annotation: {0}".format(exc))


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
    attrs.width = int(width)
    attrs.height = int(height)
    if hasattr(attrs, "screenCapture"):
        attrs.screenCapture = 0
    SetSaveWindowAttributes(attrs)
    SaveWindow()


def resolve_states(spec, n_states):
    out = []
    for raw in [s.strip() for s in str(spec).split(",") if s.strip()]:
        lower = raw.lower()
        if lower == "initial":
            out.append(("initial", 0))
        elif lower == "final":
            out.append(("final", max(n_states - 1, 0)))
        else:
            idx = int(raw)
            if idx < 0:
                idx = max(n_states + idx, 0)
            label = "state{0}".format(idx)
            out.append((label, min(idx, max(n_states - 1, 0))))
    return out or [("initial", 0), ("final", max(n_states - 1, 0))]


def render_scene(label, state, exprs, args, output_dir, case_name):
    try:
        call_visit("DeleteAllPlots")
    except Exception:
        pass
    try:
        call_visit("SetTimeSliderState", int(state))
    except Exception:
        pass

    variable = exprs["velocityMagnitude"]
    if label == "initial":
        suffix = "initial_velocityMagnitude"
    else:
        suffix = "final_velocityMagnitude"
    plot_label = "velocity magnitude"

    add_pseudocolor(variable, plot_label, args.colortable, args.point_size_pixels)
    add_mesh(args.mesh)
    call_visit("DrawPlots")
    extents = query_extents()
    fit_view(args.view, extents, args.width, args.height, args.view_padding)
    configure_annotations(not args.no_annotations)
    add_time_slider((not args.no_annotations) and args.time_slider)

    filename = case_name + "_" + suffix
    save_window(output_dir, filename, args.width, args.height)
    print("Saved {0}".format(output_dir / (filename + ".png")))


def main(argv=None, default_case_name="copperFoamCompression", default_view="xz", default_width=1200, default_height=1600):
    args = parse_args(sys.argv[1:] if argv is None else argv, default_case_name, default_view, default_width, default_height)
    run_dir = Path(expand_path(args.run_dir)).resolve()
    output_dir = Path(expand_path(args.output_dir)).resolve() if args.output_dir else run_dir / "visit_frames"
    output_dir.mkdir(parents=True, exist_ok=True)
    case_name = args.case_name or default_case_name

    print("Run directory: {0}".format(run_dir))
    print("Output directory: {0}".format(output_dir))
    database = open_database(run_dir, args.database, args.list_databases)

    scalars, vectors, meshes = metadata_variables(database)
    print("Scalar variables: {0}".format(scalars if scalars else "<none reported>"))
    print("Vector variables: {0}".format(vectors if vectors else "<none reported>"))
    print("Meshes: {0}".format(meshes if meshes else "<none reported>"))
    if args.list_variables:
        print("Variables for {0}:".format(args.region))
        for name in scalars + vectors:
            if args.region in name:
                print("  " + name)

    exprs = prepare_expressions(scalars, vectors, args.region, case_name)
    if args.dry_run:
        return 0

    try:
        n_states = int(call_visit("TimeSliderGetNStates"))
    except Exception:
        n_states = 1
    states = resolve_states(args.states, n_states)
    print("Rendering states: {0}".format(states))
    for label, state in states:
        render_scene("initial" if label == "initial" else "final", state, exprs, args, output_dir, case_name)

    try:
        call_visit("DeleteAllPlots")
    except Exception:
        pass
    try:
        call_visit("CloseDatabase", database)
    except Exception as exc:
        print("CloseDatabase warning: {0}".format(exc))
    return 0
