#!/usr/bin/env python3
"""
Render GEOS-MPM frames with VisIt.

For GEOS-MPM Silo output this script prefers the plot file family:

    <run-dir>/siloFiles/mpm_cpdi_*

Those files are the master Silo/HDF5 files that link to the payload below
siloFiles/data.  Restart ROOT files are intentionally ignored unless an explicit
--database path is supplied.
"""

import argparse
import math
import os
import re
import sys
from pathlib import Path

MPM_CPDI_RENDERER_VERSION = 9


def parse_args(argv):
    parser = argparse.ArgumentParser(description="Render GEOS-MPM frames using VisIt")
    parser.add_argument("--run-dir", default=".", help="GEOS-MPM run directory")
    parser.add_argument(
        "--database",
        default=None,
        help=("Explicit VisIt database path. Relative paths are resolved below "
              "--run-dir. VisIt wildcard database strings such as "
              "'siloFiles/mpm_cpdi_* database' are accepted."),
    )
    parser.add_argument("--output-dir", default=None, help="Output directory for PNG frames")
    parser.add_argument("--variable", default=None, help="Preferred pseudocolor variable to render")
    parser.add_argument(
        "--mesh",
        default="auto",
        help="Mesh overlay variable. Default 'auto' uses CellRegion1 when available; use 'none' to disable.",
    )
    parser.add_argument("--case-name", default=None, help="Prefix for PNG names")
    parser.add_argument("--states", default="initial,final", help="Comma-list of initial, middle, final, or integer state ids")
    parser.add_argument("--width", type=int, default=1600)
    parser.add_argument("--height", type=int, default=1200)
    parser.add_argument("--view", choices=("auto", "xy", "xz", "yz"), default="auto")
    parser.add_argument("--view-padding", type=float, default=1.04, help="Padding multiplier around fitted extents")
    parser.add_argument("--bounds", default=None, help="Override fitted extents as xmin,xmax,ymin,ymax,zmin,zmax")
    parser.add_argument("--focus", nargs=3, type=float, metavar=("X", "Y", "Z"), help="Override 3D view focus")
    parser.add_argument("--parallel-scale", type=float, help="Override 3D view parallelScale")
    parser.add_argument("--near", type=float, help="Override 3D view nearPlane")
    parser.add_argument("--far", type=float, help="Override 3D view farPlane")
    parser.add_argument("--view-angle", type=float, default=30.0)
    parser.add_argument("--no-fit-extents", action="store_true", help="Do not query/fix the orthographic view to extents")
    parser.add_argument("--perspective", action="store_true", help="Enable perspective. Default is off for explicit 2D views")
    parser.add_argument("--colortable", default="hot_desaturated", help="VisIt color table for pseudocolor plots")
    parser.add_argument("--point-size-pixels", type=int, default=5, help="Point size for particle pseudocolor plots")
    parser.add_argument("--range-mode", choices=("unit", "auto", "explicit"), default="unit", help="Pseudocolor scaling: unit fixes [0,1], auto uses VisIt/data range, explicit uses --color-min/--color-max")
    parser.add_argument("--color-min", type=float, default=None, help="Explicit pseudocolor minimum; implies --range-mode explicit when paired with --color-max")
    parser.add_argument("--color-max", type=float, default=None, help="Explicit pseudocolor maximum; implies --range-mode explicit when paired with --color-min")
    parser.add_argument("--strict-variable", action="store_true", help="Fail instead of falling back to unrelated variables when the requested variable cannot be plotted")
    parser.add_argument("--vector-variable", action="append", default=[], help="Vector variable to overlay. May be repeated.")
    parser.add_argument("--vector-scale", type=float, default=0.25, help="VisIt vector glyph scale for overlay vectors")
    parser.add_argument("--vector-count", type=int, default=500, help="Target vector glyph count when supported by VisIt")
    parser.add_argument("--vector-stride", type=int, default=1, help="Stride for vector glyphs when stride mode is enabled")
    parser.add_argument("--no-annotations", action="store_true")
    parser.add_argument("--time-slider", dest="time_slider", action="store_true", default=True, help="Add VisIt TimeSlider annotation")
    parser.add_argument("--no-time-slider", dest="time_slider", action="store_false", help="Disable TimeSlider annotation")
    parser.add_argument("--list-databases", action="store_true", help="Print database candidates before opening one")
    parser.add_argument("--dry-run", action="store_true", help="Discover databases/variables without saving windows")
    args, unknown = parser.parse_known_args(argv)
    if unknown:
        print("Ignoring unknown VisIt/script arguments: {0}".format(unknown))
    return args


def expand_local_path(path_text, base_dir=None):
    """Expand ~ and environment variables, preserving VisIt's database suffix."""
    text = os.path.expandvars(os.path.expanduser(str(path_text)))
    suffix = ""
    if text.endswith(" database"):
        text = text[:-len(" database")]
        suffix = " database"
    if base_dir is not None and not os.path.isabs(text) and ":" not in text:
        text = str(Path(base_dir) / text)
    return text + suffix


def natural_key(path):
    name = Path(path).name
    return [int(part) if part.isdigit() else part for part in re.split(r"(\d+)", name)]


def existing_files(directory, pattern):
    try:
        return sorted((p for p in directory.glob(pattern) if p.is_file()), key=natural_key)
    except Exception:
        return []


def not_restart_file(path):
    lowered = str(path).lower()
    return "restart" not in lowered and "restartoutput" not in lowered and not lowered.endswith(".root")


def make_candidate(label, pattern_path, matches, prefer_family=True):
    matches = [p for p in matches if p.is_file() and not_restart_file(p)]
    if not matches:
        return None
    if prefer_family and len(matches) > 1:
        database = str(pattern_path) + " database"
        is_family = True
    else:
        database = str(matches[0])
        is_family = False
    return {
        "label": label,
        "database": database,
        "matches": matches,
        "is_family": is_family,
    }


def discover_databases(run_dir):
    """Return database candidates, preferring GEOS-MPM plot output."""
    candidates = []
    silo_dir = run_dir / "siloFiles"

    if silo_dir.is_dir():
        for pattern in ("*.visit", "mpm_cpdi_*", "mpm_*", "*.silo", "*.h5", "*.hdf5"):
            matches = [p for p in existing_files(silo_dir, pattern) if p.name != "data"]
            candidate = make_candidate("siloFiles/" + pattern, silo_dir / pattern, matches)
            if candidate:
                candidates.append(candidate)

    for pattern in ("*.visit", "mpm_cpdi_*", "mpm_*.silo", "output*.silo", "*.silo", "*.pvd", "*.vtk"):
        candidate = make_candidate(pattern, run_dir / pattern, existing_files(run_dir, pattern))
        if candidate:
            candidates.append(candidate)

    return candidates


def write_visit_file(run_dir, candidate):
    """Create a .visit file for a matched file family and return its path string."""
    matches = candidate.get("matches") or []
    if len(matches) < 1:
        return None
    safe_label = re.sub(r"[^A-Za-z0-9_.-]+", "_", candidate.get("label", "database")).strip("_")
    visit_path = run_dir / (run_dir.name + "_" + safe_label + ".visit")
    with visit_path.open("w") as handle:
        for path in matches:
            try:
                handle.write(str(path.relative_to(run_dir)) + "\n")
            except ValueError:
                handle.write(str(path) + "\n")
    return str(visit_path)


def print_candidates(candidates):
    print("Discovered VisIt database candidates:")
    if not candidates:
        print("  <none>")
    for candidate in candidates:
        print("  {0}: {1}".format(candidate["label"], candidate["database"]))
        for match in candidate["matches"][:8]:
            print("    - {0}".format(match))
        if len(candidate["matches"]) > 8:
            print("    - ... {0} more".format(len(candidate["matches"]) - 8))


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


def open_database_from_candidates(run_dir, explicit, list_candidates=False):
    if explicit:
        database = expand_local_path(explicit, run_dir)
        if not try_open_database(database):
            raise RuntimeError("Could not open explicit VisIt database: {0}".format(database))
        return database

    candidates = discover_databases(run_dir)
    if list_candidates:
        print_candidates(candidates)

    errors = []
    for candidate in candidates:
        print("Selected VisIt database candidate: {0}".format(candidate["label"]))
        database_attempts = [candidate["database"]]
        if candidate.get("is_family"):
            visit_database = write_visit_file(run_dir, candidate)
            if visit_database:
                database_attempts.append(visit_database)
        for database in database_attempts:
            if try_open_database(database):
                if len(candidate["matches"]) > 1:
                    print("Matched {0} files in this database family.".format(len(candidate["matches"])))
                return database
            errors.append(database)

    if errors:
        raise RuntimeError("Could not open any discovered VisIt database. Tried: " + ", ".join(errors))
    raise FileNotFoundError(
        "No VisIt database found in {0}. Expected GEOS-MPM plot files such as {0}/siloFiles/mpm_cpdi_*".format(run_dir)
    )


def safe_metadata_name(value):
    try:
        name = getattr(value, "name")
    except Exception:
        return None
    try:
        text = str(name)
    except Exception:
        return None
    return text if text else None


def metadata_entries(md, count_method, get_method):
    out = []
    n_fn = getattr(md, count_method, None)
    get_fn = getattr(md, get_method, None)
    if n_fn is None or get_fn is None:
        return out
    try:
        count = int(n_fn())
    except Exception:
        return out
    for idx in range(count):
        try:
            name = safe_metadata_name(get_fn(idx))
        except Exception:
            name = None
        if name:
            out.append(name)
    return out


def metadata_variable_names(database):
    try:
        md = call_visit("GetMetaData", database)
    except Exception as exc:
        print("Could not query VisIt metadata: {0}".format(exc))
        return [], [], [], [], []
    scalars = metadata_entries(md, "GetNumScalars", "GetScalars")
    vectors = metadata_entries(md, "GetNumVectors", "GetVectors")
    meshes = metadata_entries(md, "GetNumMeshes", "GetMeshes")
    arrays = metadata_entries(md, "GetNumArrays", "GetArrays")
    tensors = metadata_entries(md, "GetNumTensors", "GetTensors")
    # Some VisIt versions split symmetric tensors into a separate metadata
    # category.  Query it opportunistically so expression construction can find
    # native tensor variables when they are present.
    tensors += [name for name in metadata_entries(md, "GetNumSymmetricTensors", "GetSymmetricTensors") if name not in tensors]
    return scalars, vectors, meshes, arrays, tensors


def compact_name(name):
    return re.sub(r"[^a-z0-9]+", "", str(name).lower())



# -----------------------------------------------------------------------------
# Derived scalar fields for example-suite visualization.
# -----------------------------------------------------------------------------
def _compact_name(text):
    return re.sub(r"[^a-z0-9]+", "", str(text).lower())


def _quote_visit_var(name):
    return "<" + str(name) + ">"


def _visit_square_sum_expression(components):
    terms = ["({0})*({0})".format(_quote_visit_var(name)) for name in components]
    return "sqrt(" + " + ".join(terms) + ")"


def _derived_priority(name):
    low = str(name).lower()
    score = 0
    if "particleregion" in low:
        score -= 200
    if "particledomains" in low:
        score -= 100
    if "particlefields" in low:
        score -= 100
    if "particle" in low:
        score -= 25
    if "ghareb" in low:
        score -= 20
    if "null_" in low or "/null" in low:
        score += 10000
    if "cellregion" in low:
        score += 100
    return (score, natural_key(name))


def _find_best_scalar(scalars, suffixes=(), contains=()):
    matches = []
    compact_suffixes = [_compact_name(s) for s in suffixes]
    compact_contains = [_compact_name(c) for c in contains]
    for name in scalars:
        cname = _compact_name(name)
        if compact_suffixes and not any(cname.endswith(s) for s in compact_suffixes):
            continue
        if compact_contains and not all(c in cname for c in compact_contains):
            continue
        matches.append(name)
    if not matches:
        return None
    return sorted(matches, key=_derived_priority)[0]


def _component_suffix_candidates(stem, label, legacy_index):
    """Suffixes for grouped Silo component expressions and legacy flat components."""
    return (
        "particle{0}/{1}".format(stem, label),
        "{0}/{1}".format(stem, label),
        "particle{0}_{1}".format(stem, legacy_index),
        "{0}_0_{1}".format(stem, legacy_index),
        "{0}_{1}".format(stem, legacy_index),
    )


def _define_scalar_expression(name, expression):
    print("Defining expression {0} = {1}".format(name, expression))
    call_visit("DefineScalarExpression", name, expression)
    return name


def _define_tensor_expression(name, expression):
    print("Defining tensor expression {0} = {1}".format(name, expression))
    call_visit("DefineTensorExpression", name, expression)
    return name


def _unique_keep_order(values):
    out = []
    seen = set()
    for value in values:
        if value and value not in seen:
            out.append(value)
            seen.add(value)
    return out


def _quoted_candidate(path):
    return _quote_visit_var(path)


def _indexed_candidate(path, index=0):
    # VisIt's [] operator extracts a scalar component from vector-valued
    # expressions.  This is needed for GEOS/Silo particleStress component leaves
    # on VisIt 3.5, where paths such as particleStress/XX are shown in scalar
    # menus but are evaluated by the engine as vector variables.
    return "{0}[{1}]".format(_quote_visit_var(path), int(index))


def _array_decompose_candidate(path, index):
    return "array_decompose({0},{1})".format(_quote_visit_var(path), int(index))


def _stress_base_candidates(all_names):
    candidates = [
        "ParticleRegion1_ParticleDomains_ParticleFields/particleStress",
        "particleStress",
    ]
    for name in all_names:
        cname = _compact_name(name)
        if cname.endswith("particlestress") or cname == "particlestress":
            candidates.append(name)
    return _unique_keep_order(candidates)


def _stress_leaf_candidates(all_names, label, legacy_index):
    labels = [label, label.lower(), label.upper(), str(legacy_index)]
    candidates = []
    for base in _stress_base_candidates(all_names):
        for leaf in labels:
            candidates.append("{0}/{1}".format(base, leaf))
        candidates.append("{0}_{1}".format(base, legacy_index))
    suffixes = _component_suffix_candidates("Stress", label.lower(), legacy_index)
    for name in all_names:
        cname = _compact_name(name)
        if "stress" not in cname:
            continue
        if any(cname.endswith(_compact_name(suffix)) for suffix in suffixes):
            candidates.append(name)
    return _unique_keep_order(candidates)


def _stress_component_expression_candidates(all_names, label, legacy_index):
    expressions = []

    # Preferred path: decompose the native GEOS Voigt stress array.  This creates
    # true scalar component expressions without adding any solver-side scalar
    # fields.
    for base in _stress_base_candidates(all_names):
        expressions.append(_array_decompose_candidate(base, legacy_index))

    # Fallback path for VisIt/Silo outputs that only expose component leaves.
    # The vector-index form handles component leaves that appear in scalar menus
    # but evaluate as vector variables.  The direct form remains useful for older
    # outputs where the leaves are genuine scalars.
    for leaf in _stress_leaf_candidates(all_names, label, legacy_index):
        expressions.append(_indexed_candidate(leaf, 0))
        expressions.append(_array_decompose_candidate(leaf, 0))
        expressions.append(_quoted_candidate(leaf))

    return _unique_keep_order(expressions)


def _define_stress_component_expression(name, all_names, label, legacy_index):
    names = []
    for i, expression in enumerate(_stress_component_expression_candidates(all_names, label, legacy_index)):
        expr_name = name if i == 0 else "{0}_candidate{1}".format(name, i)
        _define_scalar_expression(expr_name, expression)
        names.append(expr_name)
    return names


def _define_particle_stress_tensor_candidates(all_names):
    component_info = (
        ("XX", 0, "particleStressXX"),
        ("YY", 1, "particleStressYY"),
        ("ZZ", 2, "particleStressZZ"),
        ("YZ", 3, "particleStressYZ"),
        ("XZ", 4, "particleStressXZ"),
        ("XY", 5, "particleStressXY"),
    )
    component_candidates = {}
    for label, index, base_name in component_info:
        component_candidates[label] = _define_stress_component_expression(base_name, all_names, label, index)

    tensor_names = []
    max_candidates = max(len(v) for v in component_candidates.values())
    for i in range(max_candidates):
        def pick(label):
            candidates = component_candidates[label]
            return candidates[min(i, len(candidates) - 1)]

        sxx = pick("XX")
        syy = pick("YY")
        szz = pick("ZZ")
        syz = pick("YZ")
        sxz = pick("XZ")
        sxy = pick("XY")
        tensor_name = "particleStressTensor" if i == 0 else "particleStressTensor_candidate{0}".format(i)
        # Full symmetric 3x3 tensor in VisIt's row-major tensor expression
        # syntax.  VisIt's effective_tensor() then computes sqrt(3*J2), i.e.
        # the von-Mises equivalent stress, from this tensor-valued expression.
        tensor_expr = "{{{{{sxx},{sxy},{sxz}}},{{{sxy},{syy},{syz}}},{{{sxz},{syz},{szz}}}}}".format(
            sxx=sxx, syy=syy, szz=szz, syz=syz, sxz=sxz, sxy=sxy,
        )
        _define_tensor_expression(tensor_name, tensor_expr)
        tensor_names.append(tensor_name)
    return tensor_names


def _derive_particle_stress_component(all_names, requested):
    key = _compact_name(requested or "")
    component_map = {
        "particlestressxx": ("XX", 0),
        "stressxx": ("XX", 0),
        "sigmaxx": ("XX", 0),
        "particlestressyy": ("YY", 1),
        "stressyy": ("YY", 1),
        "sigmayy": ("YY", 1),
        "particlestresszz": ("ZZ", 2),
        "stresszz": ("ZZ", 2),
        "sigmazz": ("ZZ", 2),
        "particlestressyz": ("YZ", 3),
        "stressyz": ("YZ", 3),
        "sigmayz": ("YZ", 3),
        "particlestressxz": ("XZ", 4),
        "stressxz": ("XZ", 4),
        "sigmaxz": ("XZ", 4),
        "particlestressxy": ("XY", 5),
        "stressxy": ("XY", 5),
        "sigmaxy": ("XY", 5),
    }
    for compact, value in component_map.items():
        if key == compact or compact in key:
            label, index = value
            return _define_stress_component_expression("particleStress{0}".format(label), all_names, label, index)
    return None


def _derive_pressure(scalars):
    # Prefer particle stress components.  Support both grouped Silo components
    # such as particleStress/xx and legacy flat components such as particleStress_0.
    sxx = _find_best_scalar(
        scalars,
        suffixes=_component_suffix_candidates("Stress", "xx", 0) + ("stress_11",),
        contains=("stress",),
    )
    syy = _find_best_scalar(
        scalars,
        suffixes=_component_suffix_candidates("Stress", "yy", 1) + ("stress_22", "stress_1_1"),
        contains=("stress",),
    )
    szz = _find_best_scalar(
        scalars,
        suffixes=_component_suffix_candidates("Stress", "zz", 2) + ("stress_33", "stress_2_2"),
        contains=("stress",),
    )
    if not (sxx and syy and szz):
        print("Could not derive pressure; missing stress components. sxx={0}, syy={1}, szz={2}".format(sxx, syy, szz))
        return None
    expression = "-(({0}) + ({1}) + ({2}))/3.0".format(_quote_visit_var(sxx), _quote_visit_var(syy), _quote_visit_var(szz))
    return _define_scalar_expression("pressure", expression)


def _derive_particle_von_mises_stress(all_names):
    tensor_names = _define_particle_stress_tensor_candidates(all_names)
    if not tensor_names:
        print("Could not derive particleVonMisesStress; no particleStress array/component candidates found")
        return None
    vm_names = []
    for i, tensor_name in enumerate(tensor_names):
        name = "particleVonMisesStress" if i == 0 else "particleVonMisesStress_candidate{0}".format(i)
        vm_names.append(_define_scalar_expression(name, "effective_tensor({0})".format(tensor_name)))
    return vm_names


def _derive_plastic_strain_magnitude(scalars):
    components = []
    # Support grouped particlePlasticStrain/xx-style component expressions and
    # legacy particlePlasticStrain_0...5 or material-specific plasticStrain_0_0...0_5 fields.
    for i, label in enumerate(("xx", "yy", "zz", "yz", "xz", "xy")):
        comp = _find_best_scalar(
            scalars,
            suffixes=_component_suffix_candidates("PlasticStrain", label, i),
            contains=("plasticStrain",),
        )
        if comp:
            components.append(comp)
    if not components:
        print("Could not derive plasticStrainMagnitude; no plastic-strain components found")
        return None
    expression = _visit_square_sum_expression(components)
    return _define_scalar_expression("plasticStrainMagnitude", expression)


def _find_best_mesh(meshes, contains=()):
    compact_contains = [_compact_name(c) for c in contains]
    matches = []
    for name in meshes:
        cname = _compact_name(name)
        if compact_contains and not all(c in cname for c in compact_contains):
            continue
        matches.append(name)
    if not matches:
        return None
    return sorted(matches, key=_derived_priority)[0]


def _component_scalar(scalars, stem, component, field_index=None):
    label = ("x", "y", "z")[component]
    suffixes = [
        "{0}_{1}".format(stem, component),
        "{0}/{1}".format(stem, label),
        "{0}/{1}".format(stem, component),
    ]
    if field_index is not None:
        suffixes = [
            "{0}_{1}_{2}".format(stem, field_index, component),
            "{0}_{1}_{2}".format(stem, component, field_index),
        ] + suffixes
    return _find_best_scalar(scalars, suffixes=tuple(suffixes), contains=(stem,))


def _known_component_paths(stem, component, field_index=None):
    """Likely GEOS/VisIt component variable names for grouped Silo vectors.

    Some VisIt/Silo combinations expose vector components in the plot menu but
    do not report those leaves through metadata.GetScalars().  In that case a
    direct expression such as
    <ParticleRegion1_ParticleDomains_ParticleFields/particleVelocity/x> still
    plots correctly.  Prefer these deterministic component paths over falling
    back to damage or material fields when the user requested a velocity.
    """
    label = ("x", "y", "z")[component]
    paths = []
    if stem == "particleVelocity":
        paths.extend([
            "ParticleRegion1_ParticleDomains_ParticleFields/particleVelocity/{0}".format(label),
            "ParticleRegion1_ParticleDomains_ParticleFields/particleVelocity/{0}".format(component),
            "particleVelocity/{0}".format(label),
            "particleVelocity/{0}".format(component),
        ])
    elif stem == "gridVelocity":
        field = 0 if field_index is None else int(field_index)
        paths.extend([
            "CellRegion1/gridVelocity/{0}".format(label),
            "CellRegion1/gridVelocity/{0}".format(component),
            "gridVelocity/{0}".format(label),
            "gridVelocity/{0}".format(component),
            "gridVelocity_{0}_{1}".format(field, component),
            "gridVelocity_{0}_{1}".format(field, label),
        ])
    return paths


def _derive_velocity_component(scalars, stem, output_name, component, field_index=None):
    comp = _component_scalar(scalars, stem, component, field_index=field_index)
    candidates = []
    if comp:
        candidates.append(comp)
    for candidate in _known_component_paths(stem, component, field_index=field_index):
        if candidate not in candidates:
            candidates.append(candidate)
    if not candidates:
        print("Could not derive {0}; missing {1} component {2}".format(output_name, stem, component))
        return None

    # Try the best metadata-reported component first.  If metadata omitted the
    # component, fall back to common GEOS grouped-component paths.  VisIt does not
    # always validate expression definitions until the plot is added, so the main
    # candidate loop still controls success/failure.
    expression = _quote_visit_var(candidates[0])
    if len(candidates) > 1:
        print("Velocity component candidates for {0}: {1}".format(output_name, candidates))
    return _define_scalar_expression(output_name, expression)


def _derive_grid_displacement_magnitude(scalars):
    components = []
    for component in range(3):
        comp = _component_scalar(scalars, "gridDisplacement", component, field_index=0)
        if comp:
            components.append(comp)
    if len(components) < 2:
        print("Could not derive gridDisplacementMagnitude; missing gridDisplacement components: {0}".format(components))
        return None
    expression = _visit_square_sum_expression(components)
    return _define_scalar_expression("gridDisplacementMagnitude", expression)


def _derive_particle_displacement_magnitude(scalars, meshes):
    refs = []
    for component in range(3):
        comp = _component_scalar(scalars, "particleReferencePosition", component)
        if comp:
            refs.append(comp)
    if len(refs) < 2:
        print("Could not derive particleDisplacementMagnitude; missing particleReferencePosition components: {0}".format(refs))
        return None
    mesh = _find_best_mesh(meshes, contains=("ParticleRegion", "ParticleDomains")) or _find_best_mesh(meshes, contains=("Particle",))
    if not mesh:
        print("Could not derive particleDisplacementMagnitude; particle mesh not found")
        return None
    terms = []
    for component, ref in enumerate(refs):
        coord = "coord({0})[{1}]".format(_quote_visit_var(mesh), component)
        terms.append("(({0}) - ({1}))*(({0}) - ({1}))".format(coord, _quote_visit_var(ref)))
    expression = "sqrt(" + " + ".join(terms) + ")"
    return _define_scalar_expression("particleDisplacementMagnitude", expression)


def prepare_derived_requested_variable(scalars, requested, meshes=None, arrays=None, tensors=None, vectors=None):
    meshes = meshes or []
    arrays = arrays or []
    tensors = tensors or []
    vectors = vectors or []
    all_names = _unique_keep_order(list(scalars or []) + list(vectors or []) + list(arrays or []) + list(tensors or []))
    key = _compact_name(requested or "")
    stress_component = _derive_particle_stress_component(all_names, requested)
    if stress_component:
        return stress_component
    if key in ("pressure", "meanstress", "hydrostaticpressure"):
        return _derive_pressure(scalars)
    if key in ("particlevonmisesstress", "vonmisesstress", "vonmises", "equivalentstress", "misesstress"):
        return _derive_particle_von_mises_stress(all_names)
    if key in ("plasticstrainmagnitude", "plasticstrain", "plasticmagnitude", "equivalentplasticstrain"):
        return _derive_plastic_strain_magnitude(scalars)
    if key in ("griddisplacementmagnitude", "griddisplacement", "backgrounddisplacement"):
        return _derive_grid_displacement_magnitude(scalars)
    if key in ("particledisplacementmagnitude", "particledisplacement", "displacementmagnitude"):
        return _derive_particle_displacement_magnitude(scalars, meshes)
    if key in ("particlevelocityx", "particlexvelocity", "xparticlevelocity", "particlevx"):
        return _derive_velocity_component(scalars, "particleVelocity", "particleVelocityX", 0)
    if key in ("particlevelocityy", "particleyvelocity", "yparticlevelocity", "particlevy"):
        return _derive_velocity_component(scalars, "particleVelocity", "particleVelocityY", 1)
    if key in ("particlevelocityz", "particlezvelocity", "zparticlevelocity", "particlevz"):
        return _derive_velocity_component(scalars, "particleVelocity", "particleVelocityZ", 2)
    if key in ("gridvelocityx", "gridxvelocity", "xgridvelocity", "gridvx", "backgroundvelocityx"):
        return _derive_velocity_component(scalars, "gridVelocity", "gridVelocityX", 0, field_index=0)
    if key in ("gridvelocityy", "gridyvelocity", "ygridvelocity", "gridvy", "backgroundvelocityy"):
        return _derive_velocity_component(scalars, "gridVelocity", "gridVelocityY", 1, field_index=0)
    if key in ("gridvelocityz", "gridzvelocity", "zgridvelocity", "gridvz", "backgroundvelocityz"):
        return _derive_velocity_component(scalars, "gridVelocity", "gridVelocityZ", 2, field_index=0)
    return None


def score_variable(name, requested):
    compact = compact_name(name)
    req = compact_name(requested) if requested else ""

    if requested:
        if compact == req:
            return (0, 0, name)
        if req in compact:
            if "particle" in compact and "damage" in compact:
                return (0, 1, name)
            if "quartz" in compact and "damage" in compact:
                return (0, 2, name)
            return (0, 10, name)

    priority_tests = [
        ("particledamage", 10),
        ("quartzdamage", 11),
        ("damage", 12),
        ("particlematerialtype", 20),
        ("materialtype", 21),
        ("particlecolor", 22),
        ("particlecztag", 23),
        ("particlestrengthscale", 30),
        ("strengthscale", 31),
        ("particleplasticstrain", 40),
        ("plasticstrain", 41),
        ("particletemperature", 50),
        ("temperature", 51),
        ("particledensity", 80),
        ("density", 81),
        ("elementvolume", 90),
    ]
    for token, score in priority_tests:
        if token in compact:
            return (1, score, name)
    return (9, 999, name)


def ordered_scalar_candidates(scalars, requested):
    seen = set()
    candidates = []
    for name in sorted(scalars, key=lambda n: score_variable(n, requested)):
        if name not in seen:
            candidates.append(name)
            seen.add(name)

    extra = []
    if requested:
        extra.append(requested)
    extra.extend(["particleDamage", "Damage", "damage", "quartz_damage_0", "particleMaterialType", "particleColor", "particleDensity", "density"])
    for name in extra:
        if name and name not in seen:
            candidates.append(name)
            seen.add(name)
    return candidates


def configure_pseudocolor(point_size_pixels, colortable, range_mode="unit", color_min=None, color_max=None):
    try:
        PseudocolorAttributes = globals()["PseudocolorAttributes"]
        SetPlotOptions = globals()["SetPlotOptions"]
    except KeyError:
        try:
            import __main__
            PseudocolorAttributes = getattr(__main__, "PseudocolorAttributes")
            SetPlotOptions = getattr(__main__, "SetPlotOptions")
        except Exception:
            return
    try:
        attrs = PseudocolorAttributes()
        if hasattr(attrs, "colorTableName") and colortable:
            attrs.colorTableName = str(colortable)
        if hasattr(attrs, "pointSizePixels"):
            attrs.pointSizePixels = int(point_size_pixels)
        if hasattr(attrs, "pointSizeVarEnabled"):
            attrs.pointSizeVarEnabled = 0
        if hasattr(attrs, "legendFlag"):
            attrs.legendFlag = 1
        mode = str(range_mode).lower()
        if color_min is not None and color_max is not None:
            mode = "explicit"
        if mode == "auto":
            if hasattr(attrs, "minFlag"):
                attrs.minFlag = 0
            if hasattr(attrs, "maxFlag"):
                attrs.maxFlag = 0
        else:
            if hasattr(attrs, "minFlag"):
                attrs.minFlag = 1
            if hasattr(attrs, "maxFlag"):
                attrs.maxFlag = 1
            if hasattr(attrs, "min"):
                attrs.min = float(color_min) if mode == "explicit" and color_min is not None else 0.0
            if hasattr(attrs, "max"):
                attrs.max = float(color_max) if mode == "explicit" and color_max is not None else 1.0
        SetPlotOptions(attrs)
    except Exception as exc:
        print("Could not configure Pseudocolor attributes: {0}".format(exc))


def try_add_pseudocolor(variable, point_size_pixels, colortable, range_mode="unit", color_min=None, color_max=None):
    try:
        call_visit("AddPlot", "Pseudocolor", variable)
        configure_pseudocolor(point_size_pixels, colortable, range_mode, color_min=color_min, color_max=color_max)
        # VisIt often accepts an expression at AddPlot time and only discovers
        # rank/type errors when the engine draws the plot.  Validate each
        # candidate immediately so strict render requests can fall through to
        # the next stress-component expression instead of saving a blank/error
        # frame.
        result = call_visit("DrawPlots")
        if result == 0:
            raise RuntimeError("DrawPlots returned failure")
        return True
    except Exception as exc:
        print("Could not add Pseudocolor plot for {0}: {1}".format(variable, exc))
        try:
            call_visit("DeleteAllPlots")
        except Exception:
            pass
        return False


def choose_mesh(meshes, requested):
    if requested is None or str(requested).strip() == "" or str(requested).lower() == "auto":
        lower_to_name = {str(mesh).lower(): mesh for mesh in meshes}
        for name in ("CellRegion1", "CellRegion", "cellregion1", "cellregion"):
            if name.lower() in lower_to_name:
                return lower_to_name[name.lower()]
        for mesh in meshes:
            if str(mesh).lower().startswith("cellregion"):
                return mesh
        return None
    if str(requested).lower() in ("none", "off", "false", "0"):
        return None
    return requested


def try_add_mesh(mesh_variable):
    if not mesh_variable:
        return False
    try:
        call_visit("AddPlot", "Mesh", mesh_variable)
        try:
            MeshAttributes = globals().get("MeshAttributes")
            SetPlotOptions = globals().get("SetPlotOptions")
            if MeshAttributes is not None and SetPlotOptions is not None:
                attrs = MeshAttributes()
                if hasattr(attrs, "legendFlag"):
                    attrs.legendFlag = 0
                SetPlotOptions(attrs)
        except Exception:
            pass
        print("Added mesh overlay: {0}".format(mesh_variable))
        return True
    except Exception as exc:
        print("Could not add mesh overlay {0}: {1}".format(mesh_variable, exc))
        return False


def _find_best_vector(vectors, suffixes=(), contains=()):
    matches = []
    compact_suffixes = [_compact_name(s) for s in suffixes]
    compact_contains = [_compact_name(c) for c in contains]
    for name in vectors or []:
        cname = _compact_name(name)
        if compact_suffixes and not any(cname.endswith(s) for s in compact_suffixes):
            continue
        if compact_contains and not all(c in cname for c in compact_contains):
            continue
        matches.append(name)
    if not matches:
        return None
    return sorted(matches, key=_derived_priority)[0]


def _known_vector_paths(requested):
    key = _compact_name(requested or "")
    if key in ("particlesurfacenormal", "surfacenormal", "surfaceNormal"):
        return [
            "ParticleRegion1_ParticleDomains_ParticleFields/particleSurfaceNormal",
            "particleSurfaceNormal",
        ]
    if key in ("particlesurfaceposition", "surfaceposition", "surfacePosition"):
        return [
            "ParticleRegion1_ParticleDomains_ParticleFields/particleSurfacePosition",
            "particleSurfacePosition",
        ]
    if key in ("particlevelocity", "velocity"):
        return [
            "ParticleRegion1_ParticleDomains_ParticleFields/particleVelocity",
            "particleVelocity",
        ]
    return [str(requested)] if requested else []


def choose_vector(vectors, requested):
    requested = str(requested or "").strip()
    if not requested:
        return None
    candidates = []
    if requested in vectors:
        candidates.append(requested)
    found = _find_best_vector(vectors, suffixes=(requested,), contains=())
    if found and found not in candidates:
        candidates.append(found)
    found = _find_best_vector(vectors, contains=(requested,))
    if found and found not in candidates:
        candidates.append(found)
    for candidate in _known_vector_paths(requested):
        if candidate not in candidates:
            candidates.append(candidate)
    return candidates[0] if candidates else None


def configure_vector(scale, vector_count, stride):
    try:
        VectorAttributes = globals()["VectorAttributes"]
        SetPlotOptions = globals()["SetPlotOptions"]
    except KeyError:
        try:
            import __main__
            VectorAttributes = getattr(__main__, "VectorAttributes")
            SetPlotOptions = getattr(__main__, "SetPlotOptions")
        except Exception:
            return
    try:
        attrs = VectorAttributes()
        for attr, value in (
            ("scale", float(scale)),
            ("nVectors", int(vector_count)),
            ("stride", max(1, int(stride))),
            ("useStride", 0),
            ("autoScale", 1),
            ("scaleByMagnitude", 1),
            ("origOnly", 1),
            ("legendFlag", 1),
            ("useLegend", 1),
        ):
            if hasattr(attrs, attr):
                try:
                    setattr(attrs, attr, value)
                except Exception:
                    pass
        SetPlotOptions(attrs)
    except Exception as exc:
        print("Could not configure Vector attributes: {0}".format(exc))


def try_add_vector(vector_variable, vectors, scale, vector_count, stride):
    resolved = choose_vector(vectors, vector_variable)
    if not resolved:
        print("Could not resolve requested vector overlay: {0}".format(vector_variable))
        return False
    try:
        call_visit("AddPlot", "Vector", resolved)
        configure_vector(scale, vector_count, stride)
        print("Added vector overlay: {0}".format(resolved))
        return True
    except Exception as exc:
        print("Could not add Vector plot for {0}: {1}".format(resolved, exc))
        return False


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


def parse_bounds(bounds_text):
    if not bounds_text:
        return None
    values = flatten_numbers(bounds_text)
    if len(values) == 4:
        values = [values[0], values[1], values[2], values[3], 0.0, 0.0]
    if len(values) < 6:
        raise ValueError("--bounds must be xmin,xmax,ymin,ymax[,zmin,zmax]")
    return tuple(float(v) for v in values[:6])


def query_spatial_extents():
    for query_name in ("SpatialExtents", "OriginalSpatialExtents"):
        try:
            call_visit("Query", query_name)
            value = call_visit("GetQueryOutputValue")
            numbers = flatten_numbers(value)
            if len(numbers) >= 6:
                extents = tuple(numbers[:6])
                print("{0}: {1}".format(query_name, extents))
                return extents
            if len(numbers) == 4:
                extents = (numbers[0], numbers[1], numbers[2], numbers[3], 0.0, 0.0)
                print("{0}: {1}".format(query_name, extents))
                return extents
            print("{0} returned unusable value: {1}".format(query_name, value))
        except Exception as exc:
            print("Could not query {0}: {1}".format(query_name, exc))
    return None


def sanitized_extents(extents):
    if extents is None:
        return None
    vals = [float(v) for v in extents]
    for i in (0, 2, 4):
        if vals[i + 1] < vals[i]:
            vals[i], vals[i + 1] = vals[i + 1], vals[i]
        if not math.isfinite(vals[i]) or not math.isfinite(vals[i + 1]):
            return None
    return tuple(vals)


def view_axes(view):
    if view == "xy":
        return {
            "normal": (0.0, 0.0, 1.0),
            "up": (0.0, 1.0, 0.0),
            "h": (0, 1),
            "v": (2, 3),
            "d": (4, 5),
        }
    if view == "xz":
        return {
            "normal": (0.0, -1.0, 0.0),
            "up": (0.0, 0.0, 1.0),
            "h": (0, 1),
            "v": (4, 5),
            "d": (2, 3),
        }
    if view == "yz":
        return {
            "normal": (1.0, 0.0, 0.0),
            "up": (0.0, 0.0, 1.0),
            "h": (2, 3),
            "v": (4, 5),
            "d": (0, 1),
        }
    return None


def get_view3d_object():
    try:
        return call_visit("GetView3D")
    except Exception:
        try:
            View3DAttributes = globals().get("View3DAttributes")
            if View3DAttributes is None:
                import __main__
                View3DAttributes = getattr(__main__, "View3DAttributes")
            return View3DAttributes()
        except Exception:
            raise RuntimeError("View3DAttributes/GetView3D unavailable; cannot configure view")


def configure_fitted_view(view, extents, width, height, padding, perspective, focus_override=None,
                          parallel_scale_override=None, near_override=None, far_override=None, view_angle=30.0):
    if view == "auto":
        return
    axes = view_axes(view)
    if axes is None:
        return

    extents = sanitized_extents(extents)
    if extents is None:
        print("No valid spatial extents available; setting only view orientation/fixed overrides")
    x0, x1, y0, y1, z0, z1 = extents if extents else (-1.0, 1.0, -1.0, 1.0, -1.0, 1.0)
    center = tuple(float(v) for v in focus_override) if focus_override else ((x0 + x1) * 0.5, (y0 + y1) * 0.5, (z0 + z1) * 0.5)

    h0, h1 = axes["h"]
    v0, v1 = axes["v"]
    d0, d1 = axes["d"]
    hspan = max(abs(extents[h1] - extents[h0]), 1.0e-12) if extents else 2.0
    vspan = max(abs(extents[v1] - extents[v0]), 1.0e-12) if extents else 2.0
    dspan = max(abs(extents[d1] - extents[d0]), 1.0e-12) if extents else 2.0
    aspect = float(width) / float(height) if height else 1.0

    if parallel_scale_override is not None:
        parallel_scale = float(parallel_scale_override)
    else:
        parallel_scale = 0.5 * max(vspan, hspan / max(aspect, 1.0e-12)) * float(padding)

    clip = max(dspan, parallel_scale, hspan, vspan) * 10.0
    near_plane = float(near_override) if near_override is not None else -clip
    far_plane = float(far_override) if far_override is not None else clip

    v = get_view3d_object()
    v.viewNormal = axes["normal"]
    v.viewUp = axes["up"]
    v.focus = center
    if hasattr(v, "viewAngle"):
        v.viewAngle = float(view_angle)
    v.parallelScale = parallel_scale
    v.nearPlane = near_plane
    v.farPlane = far_plane
    if hasattr(v, "imagePan"):
        v.imagePan = (0.0, 0.0)
    if hasattr(v, "imageZoom"):
        v.imageZoom = 1.0
    if hasattr(v, "shear"):
        v.shear = (0.0, 0.0, 1.0)
    if hasattr(v, "perspective"):
        v.perspective = 1 if perspective else 0
    if hasattr(v, "eyeAngle"):
        v.eyeAngle = 2.0
    if hasattr(v, "centerOfRotationSet"):
        v.centerOfRotationSet = 1
    if hasattr(v, "centerOfRotation"):
        v.centerOfRotation = center
    if hasattr(v, "axis3DScaleFlag"):
        v.axis3DScaleFlag = 0
    if hasattr(v, "axis3DScales"):
        v.axis3DScales = (1.0, 1.0, 1.0)

    call_visit("SetView3D", v)
    print(
        "Configured {0} view: viewNormal={1}, focus={2}, viewUp={3}, horizontal_span={4}, "
        "vertical_span={5}, aspect={6}, parallelScale={7}, nearPlane={8}, farPlane={9}, perspective={10}".format(
            view, axes["normal"], center, axes["up"], hspan, vspan, aspect,
            parallel_scale, near_plane, far_plane, 1 if perspective else 0
        )
    )


def configure_annotations(enabled):
    try:
        AnnotationAttributes = globals().get("AnnotationAttributes")
        SetAnnotationAttributes = globals().get("SetAnnotationAttributes")
        if AnnotationAttributes is None or SetAnnotationAttributes is None:
            import __main__
            AnnotationAttributes = getattr(__main__, "AnnotationAttributes")
            SetAnnotationAttributes = getattr(__main__, "SetAnnotationAttributes")
    except Exception:
        return
    annot = AnnotationAttributes()
    try:
        annot.axes3D.visible = 1 if enabled else 0
    except Exception:
        pass
    try:
        annot.axes2D.visible = 1 if enabled else 0
    except Exception:
        pass
    for attr, value in (
        ("userInfoFlag", 0),
        ("databaseInfoFlag", 0),
        ("timeInfoFlag", 1 if enabled else 0),
        ("legendInfoFlag", 1 if enabled else 0),
    ):
        if hasattr(annot, attr):
            try:
                setattr(annot, attr, value)
            except Exception:
                pass
    try:
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
    if hasattr(attrs, "outputToCurrentDirectory"):
        attrs.outputToCurrentDirectory = 0
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


def resolve_states(state_spec, n_states):
    out = []
    for raw in [s.strip() for s in state_spec.split(",") if s.strip()]:
        key = raw.lower()
        max_state = max(n_states - 1, 0)
        if key == "initial":
            out.append(("initial", 0))
        elif key == "final":
            out.append(("final", max_state))
        elif key in ("middle", "mid", "half"):
            out.append(("middle", max_state // 2))
        elif key in ("quarter", "firstquarter", "q1"):
            out.append(("quarter", max_state // 4))
        elif key in ("threequarter", "threequarters", "thirdquarter", "q3"):
            out.append(("threequarter", (3 * max_state) // 4))
        elif key.endswith("%"):
            fraction = max(0.0, min(1.0, float(key[:-1]) / 100.0))
            out.append(("pct{0}".format(int(round(100.0 * fraction))), int(round(fraction * max_state))))
        else:
            idx = int(raw)
            if idx < 0:
                idx = max(n_states + idx, 0)
            out.append(("state{0}".format(idx), min(idx, max_state)))
    seen = set()
    unique = []
    for label, idx in out:
        key = (label, idx)
        if key not in seen:
            unique.append((label, idx))
            seen.add(key)
    return unique


def basename_of_variable(variable):
    text = str(variable)
    if "/" in text:
        text = text.rsplit("/", 1)[-1]
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", text).strip("_") or "variable"


def main(argv=None):
    args = parse_args(sys.argv[1:] if argv is None else argv)
    run_dir = Path(expand_local_path(args.run_dir)).resolve()
    output_dir = Path(expand_local_path(args.output_dir)).resolve() if args.output_dir else run_dir / "visit_frames"
    output_dir.mkdir(parents=True, exist_ok=True)
    case_name = args.case_name or run_dir.name

    print("Run directory: {0}".format(run_dir))
    print("Output directory: {0}".format(output_dir))
    database = open_database_from_candidates(run_dir, args.database, args.list_databases)

    scalars, vectors, meshes, arrays, tensors = metadata_variable_names(database)
    print("Scalar variables: {0}".format(scalars if scalars else "<none reported>"))
    print("Vector variables: {0}".format(vectors if vectors else "<none reported>"))
    print("Array variables: {0}".format(arrays if arrays else "<none reported>"))
    print("Tensor variables: {0}".format(tensors if tensors else "<none reported>"))
    print("Meshes: {0}".format(meshes if meshes else "<none reported>"))
    if args.dry_run:
        return 0

    try:
        call_visit("DeleteAllPlots")
    except Exception:
        pass

    derived_variable = prepare_derived_requested_variable(
        scalars,
        args.variable,
        meshes,
        arrays=arrays,
        tensors=tensors,
        vectors=vectors,
    )
    candidates = []
    if derived_variable:
        if isinstance(derived_variable, (list, tuple)):
            for variable in derived_variable:
                if variable not in candidates:
                    candidates.append(variable)
        else:
            candidates.append(derived_variable)
    if not args.strict_variable:
        for candidate in ordered_scalar_candidates(scalars, args.variable):
            if candidate not in candidates:
                candidates.append(candidate)
    elif args.variable:
        req = compact_name(args.variable)
        for candidate in scalars:
            cname = compact_name(candidate)
            if (req == cname or req in cname or cname in req) and candidate not in candidates:
                candidates.append(candidate)
        # A strict requested variable should fail loudly instead of silently
        # plotting damage or material type.  Add the literal request as a last
        # chance for exact VisIt variable names supplied by the caller.
        if args.variable not in candidates:
            candidates.append(args.variable)
    print("Scalar candidate order: {0}".format(candidates[:12]))

    plotted_variable = None
    for variable in candidates:
        if try_add_pseudocolor(
            variable,
            args.point_size_pixels,
            args.colortable,
            args.range_mode,
            color_min=args.color_min,
            color_max=args.color_max,
        ):
            plotted_variable = variable
            break
    if plotted_variable is None:
        if args.strict_variable:
            raise RuntimeError("Could not create a pseudocolor plot for requested variable {0}".format(args.variable))
        raise RuntimeError("Could not create a pseudocolor plot for any scalar variable candidate")

    mesh_variable = choose_mesh(meshes, args.mesh)
    if mesh_variable:
        try_add_mesh(mesh_variable)

    for vector_variable in args.vector_variable:
        try_add_vector(vector_variable, vectors, args.vector_scale, args.vector_count, args.vector_stride)

    call_visit("DrawPlots")

    extents = parse_bounds(args.bounds) if args.bounds else None
    if extents:
        print("Using explicit --bounds: {0}".format(extents))
    elif not args.no_fit_extents and args.view != "auto":
        extents = query_spatial_extents()

    if args.view != "auto":
        # Explicit planar views are orthographic/parallel by default; --perspective
        # must be requested to turn perspective back on.
        configure_fitted_view(
            args.view,
            extents,
            args.width,
            args.height,
            args.view_padding,
            perspective=args.perspective,
            focus_override=args.focus,
            parallel_scale_override=args.parallel_scale,
            near_override=args.near,
            far_override=args.far,
            view_angle=args.view_angle,
        )

    configure_annotations(not args.no_annotations)
    add_time_slider((not args.no_annotations) and args.time_slider)

    try:
        n_states = int(call_visit("TimeSliderGetNStates"))
    except Exception:
        n_states = 1
    states = resolve_states(args.states, n_states)
    print("Rendering states: {0}; variable={1}".format(states, plotted_variable))

    safe_variable = basename_of_variable(plotted_variable)
    for label, state in states:
        try:
            call_visit("SetTimeSliderState", state)
        except Exception:
            pass
        if args.view != "auto":
            configure_fitted_view(
                args.view,
                extents,
                args.width,
                args.height,
                args.view_padding,
                perspective=args.perspective,
                focus_override=args.focus,
                parallel_scale_override=args.parallel_scale,
                near_override=args.near,
                far_override=args.far,
                view_angle=args.view_angle,
            )
        filename = "{0}_{1}_{2}".format(case_name, label, safe_variable)
        save_window(output_dir, filename, args.width, args.height)
        print("Saved {0}".format(output_dir / (filename + ".png")))

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
