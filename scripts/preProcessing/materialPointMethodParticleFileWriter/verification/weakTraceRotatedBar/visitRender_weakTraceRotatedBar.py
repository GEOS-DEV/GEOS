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
    return fn(*args)


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--run-dir", default=".")
    p.add_argument("--output", default="weak_trace_rotated_bar_vm.png")
    p.add_argument("--state", default="final")
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


def region_stress_component_candidates(region, component):
    """Return likely stress-component variable paths for a particle region.

    GEOS/PFW branches have used both flat stressXX-style fields and grouped
    particleStress/xx-style fields.  The renderer tries both without changing
    the database output structure.
    """
    base = f"{region}_ParticleDomains_ParticleFields"
    lower = component.lower()
    upper = component.upper()
    return [
        f"{base}/stress{upper}",
        f"{base}/particleStress/{lower}",
        f"{base}/particleStress_{component_index(component)}",
    ]


def component_index(component):
    return {"XX": 0, "YY": 1, "ZZ": 2, "YZ": 3, "XZ": 4, "XY": 5}[component.upper()]


def variable_exists(var):
    """Best-effort VisIt variable existence check."""
    try:
        md = _call("GetMetaData")()
        names = set()
        for attr in ("scalars", "vectors", "tensors", "symmTensors"):
            values = getattr(md, attr, None)
            if values is None:
                continue
            for i in range(len(values)):
                try:
                    names.add(values[i].name)
                except Exception:
                    pass
        if names:
            return var in names
    except Exception:
        # Some VisIt builds do not expose metadata cleanly in -cli scripts.
        # In that case, optimistically let DefineScalarExpression validate.
        return True
    return False


def first_existing(candidates):
    for c in candidates:
        if variable_exists(c):
            return c
    return None


def add_vm_expression(region):
    comps = {}
    for component in ("XX", "YY", "ZZ", "YZ", "XZ", "XY"):
        comps[component] = first_existing(region_stress_component_candidates(region, component))
    missing = [k for k, v in comps.items() if not v]
    if missing:
        print(f"Could not create VisIt VM expression for {region}; missing stress components {missing}")
        return None

    sxx = quote_visit_var(comps["XX"])
    syy = quote_visit_var(comps["YY"])
    szz = quote_visit_var(comps["ZZ"])
    syz = quote_visit_var(comps["YZ"])
    sxz = quote_visit_var(comps["XZ"])
    sxy = quote_visit_var(comps["XY"])

    expr_name = f"weakTraceVM_{region}"
    expr = (
        "sqrt(0.5*((({sxx})-({syy}))^2+(({syy})-({szz}))^2+(({szz})-({sxx}))^2)"
        "+3.0*((({sxy})^2)+(({sxz})^2)+(({syz})^2)))"
    ).format(sxx=sxx, syy=syy, szz=szz, syz=syz, sxz=sxz, sxy=sxy)
    _call("DefineScalarExpression", expr_name, expr)
    print(f"Defined VisIt expression {expr_name} from stress components: {comps}")
    return expr_name


def main():
    args = parse_args()
    if not open_database(args.run_dir):
        raise SystemExit("could not open database")
    try:
        nstates = _call("TimeSliderGetNStates")()
        if args.state == "final":
            _call("SetTimeSliderState", max(0, nstates - 1))
        elif args.state == "middle":
            _call("SetTimeSliderState", max(0, nstates // 2))
        else:
            _call("SetTimeSliderState", int(args.state))
    except Exception:
        pass

    _call("DeleteAllPlots")()
    added = 0
    for region in ["ParticleRegion1", "ParticleRegion2"]:
        var = add_vm_expression(region)
        if not var:
            var = f"{region}_ParticleDomains_ParticleFields/particleDensity"
        try:
            _call("AddPlot", "Pseudocolor", var)
            pc = _call("PseudocolorAttributes")()
            pc.minFlag = 0
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
