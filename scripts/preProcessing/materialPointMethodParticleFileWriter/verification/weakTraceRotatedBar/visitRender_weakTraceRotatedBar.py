#!/usr/bin/env python3
"""Case-specific VisIt renderer for weakTraceRotatedBar.

This renderer does not require, request, or create any GEOS output field named
vonMisesStress.  The database is expected to contain the ordinary particle
stress components written by GEOS-MPM.  The von-Mises pseudocolor variable used
for the figure is a transient VisIt CLI scalar expression reconstructed from
those components separately for each particle region.
"""

import argparse
from pathlib import Path


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
    p.add_argument("--output", default="weak_trace_rotated_bar_vm.png")
    p.add_argument("--state", default="final")
    p.add_argument("--quantity", choices=["vm", "stressXX", "sigmaNN"], default="vm")
    p.add_argument("--color-min", type=float, default=None)
    p.add_argument("--color-max", type=float, default=None)
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


def quote_visit_var(name):
    return f"<{name}>"


def region_stress_array(region):
    return f"{region}_ParticleDomains_ParticleFields/particleStress"


def stress_component_expression(region, component):
    return "array_decompose({array},{index})".format(
        array=quote_visit_var(region_stress_array(region)),
        index=component_index(component),
    )


def component_index(component):
    return {"XX": 0, "YY": 1, "ZZ": 2, "YZ": 3, "XZ": 4, "XY": 5}[component.upper()]


def add_vm_expression(region):
    sxx = stress_component_expression(region, "XX")
    syy = stress_component_expression(region, "YY")
    szz = stress_component_expression(region, "ZZ")
    syz = stress_component_expression(region, "YZ")
    sxz = stress_component_expression(region, "XZ")
    sxy = stress_component_expression(region, "XY")

    expr_name = f"weakTraceVM_{region}"
    expr = (
        "sqrt(0.5*((({sxx})-({syy}))^2+(({syy})-({szz}))^2+(({szz})-({sxx}))^2)"
        "+3.0*((({sxy})^2)+(({sxz})^2)+(({syz})^2)))"
    ).format(sxx=sxx, syy=syy, szz=szz, syz=syz, sxz=sxz, sxy=sxy)
    _call("DefineScalarExpression", expr_name, expr)
    print(f"Defined VisIt expression {expr_name} from {region_stress_array(region)}")
    return expr_name


def add_sigma_nn_expression(region):
    expr_name = f"weakTraceSigmaNN_{region}"
    expr = "0.5*({sxx})+0.5*({syy})+({sxy})".format(
        sxx=stress_component_expression(region, "XX"),
        syy=stress_component_expression(region, "YY"),
        sxy=stress_component_expression(region, "XY"),
    )
    _call("DefineScalarExpression", expr_name, expr)
    print(f"Defined VisIt expression {expr_name} from {region_stress_array(region)}")
    return expr_name


def stress_xx_variable(region):
    expr_name = f"weakTraceStressXX_{region}"
    _call("DefineScalarExpression", expr_name, stress_component_expression(region, "XX"))
    print(f"Defined VisIt expression {expr_name} from {region_stress_array(region)}")
    return expr_name


def plot_variable(region, quantity):
    if quantity == "vm":
        return add_vm_expression(region)
    if quantity == "sigmaNN":
        return add_sigma_nn_expression(region)
    if quantity == "stressXX":
        return stress_xx_variable(region)
    return None


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


def main():
    args = parse_args()
    if not open_database(args.run_dir):
        raise SystemExit("could not open database")
    set_time_state(args.state)

    _call("DeleteAllPlots")()
    added = 0
    for region in ["ParticleRegion1", "ParticleRegion2"]:
        var = plot_variable(region, args.quantity)
        if not var:
            var = f"{region}_ParticleDomains_ParticleFields/particleDensity"
        try:
            _call("AddPlot", "Pseudocolor", var)
            try:
                _call("SetActivePlots", added)
            except Exception:
                pass
            pc = _call("PseudocolorAttributes")()
            if args.color_min is not None:
                pc.minFlag = 1
                pc.min = args.color_min
            else:
                pc.minFlag = 0
            if args.color_max is not None:
                pc.maxFlag = 1
                pc.max = args.color_max
            else:
                pc.maxFlag = 0
            pc.colorTableName = "hot_desaturated"
            try:
                pc.pointType = pc.Point
                pc.pointSizePixels = 5
            except Exception:
                pass
            _call("SetPlotOptions", pc)
            added += 1
        except Exception as exc:
            print(f"Skipping {region}: {exc}")
    try:
        _call("AddPlot", "Mesh", "CellRegion1")
        ma = _call("MeshAttributes")()
        ma.legendFlag = 0
        ma.lineWidth = 1
        _call("SetPlotOptions", ma)
    except Exception:
        pass
    if added == 0:
        raise SystemExit("no particle-region plots were added")
    _call("DrawPlots")()
    try:
        v = _call("View2DAttributes")()
        v.windowCoords = (-0.58, 0.58, -0.58, 0.58)
        _call("SetView2D", v)
    except Exception:
        pass
    out = Path(args.output)
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


if __name__ == "__main__":
    main()
