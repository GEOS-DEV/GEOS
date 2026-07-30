#!/usr/bin/env python3
"""VisIt renderer for weakTraceRotatedBar nodal interface diagnostics."""

import argparse
from pathlib import Path


SCALAR_QUANTITY_FIELDS = {
    "mechanism": ("gridInterfaceMechanism", False, 0.0, 5.0),
    "normalForce": ("gridInterfaceNormalForce", True, None, None),
    "tangentialForce": ("gridInterfaceTangentialForce", True, None, None),
    "normalVelocityJump": ("gridInterfaceNormalVelocityJump", True, None, None),
    "tangentialVelocityJump": ("gridInterfaceTangentialVelocityJump", True, None, None),
    "normalDisplacementJump": ("gridInterfaceNormalDisplacementJump", True, None, None),
    "tangentialDisplacementJump": ("gridInterfaceTangentialDisplacementJump", True, None, None),
    "traceActive": ("gridWeakInterfaceTraceActive", False, 0.0, 1.0),
    "traceSkipReason": ("gridWeakInterfaceTraceSkipReason", False, 0.0, 10.0),
    "traceContactSuppressed": ("gridWeakInterfaceTraceContactSuppressed", False, 0.0, 1.0),
}


VECTOR_QUANTITY_FIELDS = {
    "surfaceNormal": "gridSurfaceNormal",
    "explicitSurfaceNormal": "gridExplicitSurfaceNormal",
    "principalExplicitSurfaceNormal": "gridPrincipalExplicitSurfaceNormal",
    "surfacePosition": "gridSurfacePosition",
    "traceForce": "gridWeakInterfaceTraceForce",
    "tracePoint": "gridWeakInterfaceTracePoint",
    "traceSurfaceJump": "gridWeakInterfaceTraceSurfaceJump",
    "traceVelocityJump": "gridWeakInterfaceTraceVelocityJump",
    "traceVelocityJumpPost": "gridWeakInterfaceTraceVelocityJumpPost",
}


QUANTITY_NAMES = sorted(set(SCALAR_QUANTITY_FIELDS) | set(VECTOR_QUANTITY_FIELDS))


def _call(name, *args):
    fn = globals().get(name)
    if fn is None:
        import __main__
        fn = getattr(__main__, name)
    if not args:
        return fn
    return fn(*args)


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--run-dir", default=".")
    p.add_argument("--output", default="weak_trace_grid_mechanism.png")
    p.add_argument("--state", default="final")
    p.add_argument("--quantity", choices=QUANTITY_NAMES, default="mechanism")
    p.add_argument("--field-index", default="all", help="'all' or zero-based velocity-field index")
    p.add_argument("--field-count", type=int, default=2)
    p.add_argument("--color-min", type=float, default=None)
    p.add_argument("--color-max", type=float, default=None)
    p.add_argument("--vector-stride", type=int, default=1)
    p.add_argument("--vector-scale", type=float, default=0.04)
    return p.parse_known_args()[0]


def open_database(run_dir):
    run_dir = Path(run_dir)
    silo = run_dir / "siloFiles"
    candidates = []
    if silo.is_dir():
        matches = sorted(silo.glob("mpm_cpdi_*")) or sorted(silo.glob("mpm_*"))
        if len(matches) > 1:
            visit_file = run_dir / "weakTraceRotatedBar.visit"
            visit_file.write_text("\n".join(str(p.relative_to(run_dir)) for p in matches) + "\n")
            candidates.append(str(visit_file))
        candidates.extend(str(p) for p in matches[:1])
    for db in candidates:
        try:
            if _call("OpenDatabase", db):
                return True
        except Exception:
            pass
    return False


def set_time_state(state):
    try:
        nstates = _call("TimeSliderGetNStates")()
        if state == "final":
            target = max(0, nstates - 1)
        elif state == "middle":
            target = max(0, nstates // 2)
        elif state == "initial":
            target = 0
        elif state.startswith("fraction:"):
            fraction = min(1.0, max(0.0, float(state.split(":", 1)[1])))
            target = int(round(fraction * max(0, nstates - 1)))
        else:
            target = int(state)
        _call("SetTimeSliderState", target)
        print(f"Set VisIt time-slider state {target} of {nstates} for requested state {state}")
    except Exception as exc:
        print(f"Could not set VisIt time-slider state {state}: {exc}")


def metadata_variable_names():
    names = set()
    try:
        md = _call("GetMetaData")()
        for attr in ("scalars", "vectors", "tensors", "symmTensors"):
            values = getattr(md, attr, None)
            if values is None:
                continue
            for i in range(len(values)):
                try:
                    names.add(values[i].name)
                except Exception:
                    pass
    except Exception:
        return names
    return names


def grid_field_candidates(field_name, field_index):
    suffixes = [f"velocityField{field_index + 1}", str(field_index)]
    prefixes = ["CellRegion1_GridFields", "CellRegion1_NodalFields", "/CellRegion1_GridFields", "/CellRegion1_NodalFields", ""]
    out = []
    for prefix in prefixes:
        for suffix in suffixes:
            out.append(f"{prefix}/{field_name}_{suffix}" if prefix else f"{field_name}_{suffix}")
        out.append(f"{prefix}/{field_name}" if prefix else field_name)
    return out


def quote_var(name):
    return f"<{name}>"


def nested_max(expressions):
    if not expressions:
        return None
    out = expressions[0]
    for expr in expressions[1:]:
        out = f"max({out},{expr})"
    return out


def selected_grid_variables(field_name, field_index, field_count, metadata_names):
    indices = range(field_count) if field_index == "all" else [int(field_index)]
    variables = []
    seen = set()
    for idx in indices:
        candidates = grid_field_candidates(field_name, idx)
        if metadata_names:
            candidates = [candidate for candidate in candidates if candidate in metadata_names]
        if candidates:
            variable = candidates[0]
            if variable not in seen:
                variables.append(variable)
                seen.add(variable)
    return variables


def expression_safe_name(name):
    return "".join(ch if ch.isalnum() else "_" for ch in name).strip("_")


def define_vector_from_components(quantity, candidate, metadata_names):
    component_sets = [
        [f"{candidate}_{i}" for i in range(3)],
        [f"{candidate}/{axis}" for axis in ("X", "Y", "Z")],
    ]
    for components in component_sets:
        if metadata_names and not all(component in metadata_names for component in components):
            continue
        expr_name = f"weakTraceVector_{quantity}_{expression_safe_name(candidate)}"
        expr = "{" + ",".join(quote_var(component) for component in components) + "}"
        _call("DefineVectorExpression", expr_name, expr)
        print(f"Defined VisIt vector expression {expr_name} = {expr}")
        return expr_name
    return None


def diagnostic_scalar_variable(quantity, field_index, field_count):
    field_name, use_abs, default_min, default_max = SCALAR_QUANTITY_FIELDS[quantity]
    metadata_names = metadata_variable_names()
    variables = selected_grid_variables(field_name, field_index, field_count, metadata_names)
    if not variables:
        raise RuntimeError(f"could not find VisIt variables for {field_name}")

    terms = [f"abs({quote_var(v)})" if use_abs else quote_var(v) for v in variables]
    if len(terms) == 1:
        return variables[0], terms[0], variables, default_min, default_max

    expr_name = f"weakTraceGrid_{quantity}_{field_index}"
    expr = nested_max(terms)
    _call("DefineScalarExpression", expr_name, expr)
    print(f"Defined VisIt expression {expr_name} = {expr}")
    return expr_name, terms[0], variables, default_min, default_max


def diagnostic_vector_variables(quantity, field_index, field_count):
    field_name = VECTOR_QUANTITY_FIELDS[quantity]
    metadata_names = metadata_variable_names()
    variables = selected_grid_variables(field_name, field_index, field_count, metadata_names)
    if not variables:
        indices = range(field_count) if field_index == "all" else [int(field_index)]
        seen = set()
        for idx in indices:
            for candidate in grid_field_candidates(field_name, idx):
                expr_name = define_vector_from_components(quantity, candidate, metadata_names)
                if expr_name and expr_name not in seen:
                    variables.append(expr_name)
                    seen.add(expr_name)
                    break
    if not variables:
        raise RuntimeError(f"could not find VisIt vector variables for {field_name}")
    return variables


def _set_attr(obj, name, value):
    try:
        setattr(obj, name, value)
    except Exception:
        pass


def add_mesh_overlay():
    try:
        _call("AddPlot", "Mesh", "CellRegion1")
        ma = _call("MeshAttributes")()
        ma.legendFlag = 0
        ma.lineWidth = 1
        _call("SetPlotOptions", ma)
    except Exception:
        pass


def set_default_view():
    try:
        v = _call("View2DAttributes")()
        v.windowCoords = (-0.58, 0.58, -0.58, 0.58)
        _call("SetView2D", v)
    except Exception:
        pass


def render_scalar(args):
    plot_var, fallback_expr, variables, default_min, default_max = diagnostic_scalar_variable(
        args.quantity, args.field_index, args.field_count
    )
    _call("DeleteAllPlots")()
    try:
        _call("AddPlot", "Pseudocolor", plot_var)
    except Exception as exc:
        print(f"Could not add combined plot {plot_var}: {exc}")
        if not variables:
            raise
        plot_var = variables[0]
        _call("AddPlot", "Pseudocolor", plot_var)
    pc = _call("PseudocolorAttributes")()
    color_min = args.color_min if args.color_min is not None else default_min
    color_max = args.color_max if args.color_max is not None else default_max
    if color_min is not None:
        pc.minFlag = 1
        pc.min = color_min
    if color_max is not None:
        pc.maxFlag = 1
        pc.max = color_max
    pc.colorTableName = "viridis"
    try:
        pc.centering = pc.Nodal
    except Exception:
        pass
    _call("SetPlotOptions", pc)
    add_mesh_overlay()
    _call("DrawPlots")()
    set_default_view()


def render_vector(args):
    variables = diagnostic_vector_variables(args.quantity, args.field_index, args.field_count)
    _call("DeleteAllPlots")()
    for variable in variables:
        _call("AddPlot", "Vector", variable)
        try:
            va = _call("VectorAttributes")()
            _set_attr(va, "useStride", 1)
            _set_attr(va, "stride", max(1, args.vector_stride))
            _set_attr(va, "scale", args.vector_scale)
            _set_attr(va, "scaleByMagnitude", 0)
            _set_attr(va, "autoScale", 0)
            _set_attr(va, "colorByMagnitude", 1)
            _call("SetPlotOptions", va)
        except Exception:
            pass
    add_mesh_overlay()
    _call("DrawPlots")()
    set_default_view()


def save_window(output):
    out = Path(output)
    out.parent.mkdir(parents=True, exist_ok=True)
    swa = _call("SaveWindowAttributes")()
    try:
        swa.outputToCurrentDirectory = 0
    except Exception:
        pass
    swa.outputDirectory = str(out.parent)
    swa.fileName = out.stem
    swa.family = 0
    swa.format = swa.PNG
    swa.width = 1400
    swa.height = 1100
    _call("SetSaveWindowAttributes", swa)
    _call("SaveWindow")()


def main():
    args = parse_args()
    if not open_database(args.run_dir):
        raise SystemExit("could not open database")
    set_time_state(args.state)

    if args.quantity in VECTOR_QUANTITY_FIELDS:
        render_vector(args)
    else:
        render_scalar(args)
    save_window(args.output)


if __name__ == "__main__":
    main()
