#!/usr/bin/env python3
"""Render the GEOS-MPM Benchy contact example with VisIt.

The Benchy example intentionally uses separate particle regions/fields for the
STL boat and the impactor ball.  This script therefore creates multiple
Pseudocolor plots in one window rather than assuming a single particle region.

Default scene contents
----------------------
initial state:
  * Mesh:        CellRegion1 background grid
  * Pseudocolor: ParticleRegion1 ghostRank, intended to show the boat particles
  * Pseudocolor: ParticleRegion2 velocity magnitude, intended to show the ball

final state:
  * Mesh:        CellRegion1 background grid
  * Pseudocolor: ParticleRegion1 plastic-strain magnitude, intended to show boat deformation
  * Pseudocolor: ParticleRegion2 velocity magnitude, intended to show the ball

The script prefers GEOS-MPM Silo plot files under:

    <run-dir>/siloFiles/mpm_cpdi_*

and ignores restart ROOT files unless --database is supplied explicitly.
"""

import argparse
import math
import os
import re
import sys
from pathlib import Path

BENCHY_VISIT_RENDER_VERSION = 1


def parse_args(argv):
    parser = argparse.ArgumentParser(description="Render Benchy multi-region MPM frames with VisIt")
    parser.add_argument("--run-dir", default=".")
    parser.add_argument("--database", default=None)
    parser.add_argument("--output-dir", default=None)
    parser.add_argument("--case-name", default=None)
    parser.add_argument("--states", default="initial,final", help="Comma-list of initial, final, or state ids")
    parser.add_argument("--width", type=int, default=1600)
    parser.add_argument("--height", type=int, default=1200)
    parser.add_argument("--view", choices=("auto", "iso", "xy", "xz", "yz"), default="iso")
    parser.add_argument("--view-padding", type=float, default=1.10)
    parser.add_argument("--colortable", default="hot_desaturated")
    parser.add_argument("--ghost-colortable", default=None, help="Color table for the boat ghostRank plot; default uses --colortable")
    parser.add_argument("--point-size-pixels", type=int, default=4)
    parser.add_argument("--mesh", default="CellRegion1")
    parser.add_argument("--boat-region", default="ParticleRegion1")
    parser.add_argument("--ball-region", default="ParticleRegion2")
    parser.add_argument("--no-annotations", action="store_true")
    parser.add_argument("--time-slider", dest="time_slider", action="store_true", default=True)
    parser.add_argument("--no-time-slider", dest="time_slider", action="store_false")
    parser.add_argument("--list-databases", action="store_true")
    parser.add_argument("--list-variables", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    args, unknown = parser.parse_known_args(argv)
    if unknown:
        print("Ignoring unknown VisIt/script arguments: {0}".format(unknown))
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
        raise RuntimeError("VisIt function {0} is not available. Run with visit -cli -s.".format(name))
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
    for pattern in ("mpm_cpdi_*", "*.visit", "*.silo", "*.pvd", "*.vtk"):
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


def find_exact_or_token(scalars, region, exact_name, tokens):
    exact = field(region, exact_name)
    if exact in scalars:
        return exact
    r = compact(region)
    token_compacts = [compact(t) for t in tokens]
    matches = []
    for name in scalars:
        c = compact(name)
        if r in c and all(t in c for t in token_compacts):
            matches.append(name)
    return sorted(matches)[0] if matches else None


def choose_components(scalars, region, basename, n_expected=None, required_tokens=None):
    required_tokens = required_tokens or []
    comps = []
    if n_expected is not None:
        exact = [field(region, basename + "_" + str(i)) for i in range(n_expected)]
        if all(name in scalars for name in exact):
            return exact
    r = compact(region)
    b = compact(basename)
    token_compacts = [compact(t) for t in required_tokens]
    for name in scalars:
        c = compact(name)
        if r in c and b in c and all(t in c for t in token_compacts):
            comps.append(name)
    return sorted(comps, key=natural_string_key)


def natural_string_key(text):
    return [int(part) if part.isdigit() else part for part in re.split(r"(\d+)", str(text))]


def quote_var(name):
    return "<" + str(name) + ">"


def square_sum_expression(components):
    terms = ["({0})*({0})".format(quote_var(name)) for name in components]
    return "sqrt(" + " + ".join(terms) + ")"


def define_scalar_expression(name, expression):
    print("Defining expression {0} = {1}".format(name, expression))
    call_visit("DefineScalarExpression", name, expression)
    return name


def prepare_expressions(scalars, boat_region, ball_region):
    boat_ghost = find_exact_or_token(scalars, boat_region, "ghostRank", ["ghost", "rank"])
    if boat_ghost is None:
        raise RuntimeError("Could not find boat ghostRank in {0}. Available region variables:\n{1}".format(
            boat_region, "\n".join([s for s in scalars if boat_region in s])
        ))

    ball_vel = choose_components(scalars, ball_region, "particleVelocity", n_expected=3)
    if len(ball_vel) < 3:
        # Fall back to any three velocity-like scalars in the ball region.
        ball_vel = choose_components(scalars, ball_region, "velocity")[:3]
    if len(ball_vel) < 3:
        raise RuntimeError("Could not find three ball velocity components in {0}. Available region variables:\n{1}".format(
            ball_region, "\n".join([s for s in scalars if ball_region in s])
        ))
    ball_speed = define_scalar_expression("benchy_ball_velocityMagnitude", square_sum_expression(ball_vel[:3]))

    boat_plastic = choose_components(scalars, boat_region, "particlePlasticStrain", n_expected=6)
    if len(boat_plastic) < 1:
        boat_plastic = choose_components(scalars, boat_region, "plasticStrain")
    if len(boat_plastic) < 1:
        raise RuntimeError("Could not find boat plastic-strain components in {0}. Available region variables:\n{1}".format(
            boat_region, "\n".join([s for s in scalars if boat_region in s])
        ))
    boat_plastic_mag = define_scalar_expression("benchy_boat_plasticStrainMagnitude", square_sum_expression(boat_plastic))

    print("Boat ghostRank variable: {0}".format(boat_ghost))
    print("Ball velocity components: {0}".format(ball_vel[:3]))
    print("Boat plastic-strain components: {0}".format(boat_plastic))
    return {
        "boat_ghost": boat_ghost,
        "ball_speed": ball_speed,
        "boat_plastic_mag": boat_plastic_mag,
    }


def configure_pseudocolor(colortable, point_size_pixels, autoscale=True):
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
        if autoscale:
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
    configure_pseudocolor(colortable, point_size_pixels, autoscale=True)


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
            if hasattr(attrs, "opaqueColorSource") and hasattr(attrs, "OpaqueCustom"):
                attrs.opaqueColorSource = attrs.OpaqueCustom
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


def normalized(vec):
    mag = math.sqrt(sum(float(v) * float(v) for v in vec))
    if mag <= 0.0:
        return tuple(vec)
    return tuple(float(v) / mag for v in vec)


def fit_view(view_name, extents, width, height, padding):
    if view_name == "auto" or extents is None:
        return
    vals = list(float(v) for v in extents)
    for i in (0, 2, 4):
        if vals[i + 1] < vals[i]:
            vals[i], vals[i + 1] = vals[i + 1], vals[i]
    x0, x1, y0, y1, z0, z1 = vals
    dx = max(x1 - x0, 1.0e-9)
    dy = max(y1 - y0, 1.0e-9)
    dz = max(z1 - z0, 1.0e-9)
    center = ((x0 + x1) * 0.5, (y0 + y1) * 0.5, (z0 + z1) * 0.5)
    diag = math.sqrt(dx * dx + dy * dy + dz * dz)

    if view_name == "xy":
        view_normal = (0.0, 0.0, 1.0)
        view_up = (0.0, 1.0, 0.0)
        aspect = float(width) / float(height) if height else 1.0
        parallel_scale = 0.5 * max(dy, dx / max(aspect, 1.0e-9)) * padding
    elif view_name == "xz":
        view_normal = (0.0, -1.0, 0.0)
        view_up = (0.0, 0.0, 1.0)
        aspect = float(width) / float(height) if height else 1.0
        parallel_scale = 0.5 * max(dz, dx / max(aspect, 1.0e-9)) * padding
    elif view_name == "yz":
        view_normal = (1.0, 0.0, 0.0)
        view_up = (0.0, 0.0, 1.0)
        aspect = float(width) / float(height) if height else 1.0
        parallel_scale = 0.5 * max(dz, dy / max(aspect, 1.0e-9)) * padding
    else:
        # Orthographic isometric-ish view for the 3D Benchy demo.
        view_normal = normalized((-0.65, -0.45, 0.62))
        view_up = (0.0, 0.0, 1.0)
        parallel_scale = 0.5 * diag * padding

    v = get_view3d()
    v.viewNormal = view_normal
    v.viewUp = view_up
    v.focus = center
    v.parallelScale = parallel_scale
    v.nearPlane = -10.0 * max(diag, parallel_scale)
    v.farPlane = 10.0 * max(diag, parallel_scale)
    if hasattr(v, "viewAngle"):
        v.viewAngle = 30.0
    if hasattr(v, "perspective"):
        # Keep the batch images orthographic/reproducible.  This is especially
        # important when the same renderer is used for 2D views.
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
    print("Configured {0} view: focus={1}, viewNormal={2}, viewUp={3}, parallelScale={4}, perspective=0".format(
        view_name, center, view_normal, view_up, parallel_scale
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
            out.append(("state{0}".format(idx), min(idx, max(n_states - 1, 0))))
    return out


def render_scene(label, state, exprs, args, output_dir, case_name):
    try:
        call_visit("DeleteAllPlots")
    except Exception:
        pass
    try:
        call_visit("SetTimeSliderState", int(state))
    except Exception:
        pass

    ghost_colortable = args.ghost_colortable or args.colortable
    if label == "initial":
        add_pseudocolor(exprs["boat_ghost"], "boat ghostRank", ghost_colortable, args.point_size_pixels)
        add_pseudocolor(exprs["ball_speed"], "ball velocity magnitude", args.colortable, args.point_size_pixels)
        add_mesh(args.mesh)
        suffix = "initial_boatGhostRank_ballVelocityMagnitude"
    else:
        add_pseudocolor(exprs["boat_plastic_mag"], "boat plastic-strain magnitude", args.colortable, args.point_size_pixels)
        add_pseudocolor(exprs["ball_speed"], "ball velocity magnitude", args.colortable, args.point_size_pixels)
        add_mesh(args.mesh)
        suffix = "final_boatPlasticStrainMagnitude_ballVelocityMagnitude"

    call_visit("DrawPlots")
    extents = query_extents()
    fit_view(args.view, extents, args.width, args.height, args.view_padding)
    configure_annotations(not args.no_annotations)
    add_time_slider((not args.no_annotations) and args.time_slider)

    filename = case_name + "_" + suffix
    save_window(output_dir, filename, args.width, args.height)
    print("Saved {0}".format(output_dir / (filename + ".png")))


def main(argv=None):
    args = parse_args(sys.argv[1:] if argv is None else argv)
    run_dir = Path(expand_path(args.run_dir)).resolve()
    output_dir = Path(expand_path(args.output_dir)).resolve() if args.output_dir else run_dir / "visit_frames"
    output_dir.mkdir(parents=True, exist_ok=True)
    case_name = args.case_name or run_dir.name

    print("Run directory: {0}".format(run_dir))
    print("Output directory: {0}".format(output_dir))
    database = open_database(run_dir, args.database, args.list_databases)

    scalars, vectors, meshes = metadata_variables(database)
    print("Scalar variables: {0}".format(scalars if scalars else "<none reported>"))
    print("Vector variables: {0}".format(vectors if vectors else "<none reported>"))
    print("Meshes: {0}".format(meshes if meshes else "<none reported>"))
    if args.list_variables:
        print("Variables for {0}:".format(args.boat_region))
        for name in scalars:
            if args.boat_region in name:
                print("  " + name)
        print("Variables for {0}:".format(args.ball_region))
        for name in scalars:
            if args.ball_region in name:
                print("  " + name)

    exprs = prepare_expressions(scalars, args.boat_region, args.ball_region)
    if args.dry_run:
        return 0

    try:
        n_states = int(call_visit("TimeSliderGetNStates"))
    except Exception:
        n_states = 1
    states = resolve_states(args.states, n_states)
    print("Rendering states: {0}".format(states))
    for label, state in states:
        if label not in ("initial", "final"):
            print("State label {0} uses final-style variables.".format(label))
            label_for_scene = "final"
        else:
            label_for_scene = label
        render_scene(label_for_scene, state, exprs, args, output_dir, case_name)

    try:
        call_visit("DeleteAllPlots")
    except Exception:
        pass
    try:
        call_visit("CloseDatabase", database)
    except Exception as exc:
        print("CloseDatabase warning: {0}".format(exc))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
