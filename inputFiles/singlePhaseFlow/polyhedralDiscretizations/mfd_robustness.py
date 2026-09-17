#!/usr/bin/env python3
"""h-robustness, contrast-robustness and convergence of SinglePhaseMixedMFD on regular hexahedral and tetrahedral meshes.

Standard library only. Each case is a deck generated from mfd_robustness_template.xml.in, run serially
with geosx, and parsed from the log and the ASCII VTK output.

  python3 mfd_robustness.py --geosx /path/to/geosx                 # everything: h-sweep + contrast sweep, hex and tet
  python3 mfd_robustness.py --geosx geosx --mesh tet --levels 4,8,16
  python3 mfd_robustness.py --geosx geosx --contrasts 1,1e2,1e4,1e6 --pattern layersX --level 16
  python3 mfd_robustness.py --geosx geosx --mesh tet --mfd-percents 0,25,50,75,100 --level 8   # prescribed MFD fraction
  python3 mfd_robustness.py --geosx geosx --mfd-percent 50                                     # sweeps at a fixed 50 % MFD
  python3 mfd_robustness.py --geosx geosx --table                                              # mesh x mfd% x n x contrast table, n up to 64
  python3 mfd_robustness.py --geosx geosx --table --np 8 --table-levels 32,64 --no-output      # the large levels on 8 MPI ranks, iterations only
  python3 mfd_robustness.py --geosx geosx --table --mpirun "srun --overlap"                    # inside a Slurm step: every run through srun, serial included

Exact solution used for the convergence check (homogeneous permeability k, injection rate density q):
  p(x) = p0 + q mu / (2 rho k) x (1 - x)
"""
import argparse, glob, math, os, re, shlex, shutil, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
STATES = {"tpfa": "1.0e+20", "mixed": "0.1", "mfd": "1.0e-20"}
RHO, MU, K_REF, L = 1000.0, 1.0e-3, 1.0e-14, 1.0
VOLUME = 1.0 * 1.0 * 2.0

def linspace(a, b, m):
    return [a + (b - a) * i / m for i in range(m + 1)]

def mesh_params(mesh, n, pattern):
    """four cell blocks: layers of thickness 0.5 along z (parallel to the flow) or two bands along x (in series);
    returns the deck tokens and the cell edges along x, y, z"""
    elem = "C3D8" if mesh == "hex" else "C3D4"
    ye = linspace(0.0, 1.0, n)
    if pattern == "layersX":
        m = max(n // 2, 1)
        xc, nx = "0, 0.5, 1", f"{m}, {m}"
        zc, nz = "0, 2", str(2 * n)
        xe = [0.0] + [0.5 * b + 0.5 * i / m for b in range(2) for i in range(1, m + 1)]
        ze = linspace(0.0, 2.0, 2 * n)
    else:
        m = max(n // 2, 1)
        xc, nx = "0, 1", str(n)
        zc, nz = "0, 0.5, 1, 1.5, 2", ", ".join(str(m) for _ in range(4))
        xe = linspace(0.0, 1.0, n)
        ze = [0.0] + [0.5 * b + 0.5 * i / m for b in range(4) for i in range(1, m + 1)]
    nblocks = 2 if pattern == "layersX" else 4
    names = ", ".join(f"cb{i}" for i in range(nblocks))
    tokens = dict(elementType=elem, xCoords=xc, nx=nx, ny=str(n), zCoords=zc, nz=nz, cellBlockNames=names,
                  blocksA="cb0, cb2" if nblocks == 4 else "cb0", blocksB="cb1, cb3" if nblocks == 4 else "cb1")
    return tokens, (xe, ye, ze)

def prescription(edges, percent):
    """prescribedMfdFlag = 0 on every cell, then 1 on the first m cells in z-major order (m = percent of the
    cells, rounded), selected with up to three boxes; the resulting percentage is exact"""
    if percent is None:
        return "", ""
    xe, ye, ze = edges; nx, ny, nz = len(xe) - 1, len(ye) - 1, len(ze) - 1
    m = int(round(percent / 100.0 * nx * ny * nz))
    a, rem = divmod(m, nx * ny); b, c = divmod(rem, nx)
    eps = 1.0e-6
    boxes = []
    if a > 0: boxes.append((xe[0], xe[-1], ye[0], ye[-1], ze[0], ze[a]))
    if a < nz and b > 0: boxes.append((xe[0], xe[-1], ye[0], ye[b], ze[a], ze[a + 1]))
    if a < nz and b < ny and c > 0: boxes.append((xe[0], xe[c], ye[b], ye[b + 1], ze[a], ze[a + 1]))
    box_xml = "".join(f'    <Box name="mfdBox{i}" xMin="{{ {x0 - eps:.6f}, {y0 - eps:.6f}, {z0 - eps:.6f} }}" '
                      f'xMax="{{ {x1 + eps:.6f}, {y1 + eps:.6f}, {z1 + eps:.6f} }}"/>\n' for i, (x0, x1, y0, y1, z0, z1) in enumerate(boxes))
    spec_xml = ('    <FieldSpecification name="prescribeTpfa" initialCondition="1" setNames="{ all }" objectPath="ElementRegions" fieldName="prescribedMfdFlag" scale="0.0"/>\n'
                + "".join(f'    <FieldSpecification name="prescribeMfd{i}" initialCondition="1" setNames="{{ mfdBox{i} }}" objectPath="ElementRegions" fieldName="prescribedMfdFlag" scale="1.0"/>\n'
                          for i in range(len(boxes))))
    return box_xml, spec_xml

def write_case(template, workdir, mesh, n, pattern, contrast, state, rate, p0, percent=None, output=True):
    tokens, edges = mesh_params(mesh, n, pattern)
    tokens["prescriptionBoxes"], tokens["prescription"] = prescription(edges, percent)
    write_case.cells = tuple(len(e) - 1 for e in edges)      # cells per axis, for the MPI partitioning
    tokens.update(consistencyTolerance=STATES[state], permA=f"{K_REF:.6e}", permB=f"{K_REF / contrast:.6e}", sourceRate=f"{rate:.6e}", p0=f"{p0:.6e}")
    text = open(template).read()
    if not output:
        text = re.sub(r'\s*<PeriodicEvent name="output"[^>]*/>', "", text)
        text = re.sub(r'\s*<Outputs>.*?</Outputs>', "", text, flags=re.S)
    for k, v in tokens.items():
        text = text.replace(f"@{k}@", v)
    assert "@" not in text, "unreplaced token in the template"
    os.makedirs(workdir, exist_ok=True)
    with open(os.path.join(workdir, "deck.xml"), "w") as f:
        f.write(text)

SOLVE = re.compile(r"Linear Solver \| (\w+) \| Unknowns: ([\d,]+) \| Nonzeros: [\d,]+ \| Iterations: (\d+)"
                   r" \| Final Rel Res: \S+ \| Setup Time: ([\d.eE+-]+) s \| Solve Time: ([\d.eE+-]+) s")

def partition(np, cells):
    """split np ranks over the three axes of the internal mesh, each prime factor going to the axis
    with the most cells per rank, so that no axis gets more partitions than cells"""
    p = [1, 1, 1]; f = 2; m = np
    factors = []
    while m > 1:
        while m % f == 0: factors.append(f); m //= f
        f += 1
    for f in sorted(factors, reverse=True):
        axis = max(range(3), key=lambda i: cells[i] / p[i])
        p[axis] *= f
    return p

def run_case(geosx, workdir, timeout, np=1, launcher=None, cells=(1, 1, 1)):
    t0 = time.time()
    cmd = [geosx, "-i", "deck.xml"]
    if np > 1:
        px, py, pz = partition(np, cells)
        cmd += ["-x", str(px), "-y", str(py), "-z", str(pz)]
    if launcher:
        # "-n" is understood by mpirun, mpiexec and srun
        cmd = launcher + ["-n", str(np)] + cmd
    with open(os.path.join(workdir, "run.log"), "w") as log:
        try:
            subprocess.run(cmd, cwd=workdir, stdout=log, stderr=subprocess.STDOUT, timeout=timeout)
        except subprocess.TimeoutExpired:
            pass
    text = open(os.path.join(workdir, "run.log"), errors="replace").read()
    r = dict(wall=time.time() - t0, error="Error" in text or "Signal encountered" in text)
    m = re.search(r"eta = 1 on (\d+) / (\d+) cells", text)
    if m: r["mfdCells"], r["cells"] = int(m.group(1)), int(m.group(2))
    m = re.search(r"flux dofs (\d+) non-condensed \(saddle point\), (\d+) condensed", text)
    if m: r["live"], r["condensed"] = int(m.group(1)), int(m.group(2))
    s = SOLVE.findall(text)
    if s:
        r["status"], r["unknowns"], r["iterations"] = s[0][0], int(s[0][1].replace(",", "")), int(s[0][2])
        r["setupTime"], r["solveTime"] = float(s[0][3]), float(s[0][4])
        r["numSolves"] = len(s)
    else:
        r["error"] = True
        lines = [l for l in text.splitlines() if l.strip() and not l.startswith("-")]
        print(f"  no linear solve in {os.path.join(workdir, 'run.log')}: {lines[-1][:120] if lines else 'empty log'}", file=sys.stderr, flush=True)
    return r

def read_vtu_arrays(path, names):
    text = open(path).read()
    out = {}
    for name in names:
        m = re.search(r'<DataArray[^>]*Name="%s"[^>]*>(.*?)</DataArray>' % name, text, re.S)
        if m:
            body = re.sub(r"<InformationKey.*?</InformationKey>", "", m.group(1), flags=re.S)
            out[name] = [float(v) for v in body.split()]
    return out

def exact_series(x, q, contrast):
    """-(k p')' = q on [0,1], p(0) = p(1) = 0, k = K_REF on [0, 0.5), K_REF/contrast on [0.5, 1]"""
    k = [K_REF, K_REF / contrast]
    # u = -k p' = q x + c ; p(x) = -int_0^x (q s + c)/k(s) ds ; p(1) = 0 fixes c
    def integrals(x):
        i0 = i1 = 0.0; t0 = 0.0
        for j in range(2):
            t1 = min(0.5 * (j + 1), x)
            if t1 <= t0: break
            i0 += (t1 - t0) / k[j]; i1 += 0.5 * (t1 * t1 - t0 * t0) / k[j]; t0 = t1
        return i0, i1
    I0, I1 = integrals(1.0); c = -q * I1 / I0
    i0, i1 = integrals(x)
    return -(q * i1 + c * i0) * MU / RHO

def pressure_error(workdir, rate, contrast, p0, pattern):
    """relative L2 error of the cell pressures against the exact solution (homogeneous, or layers in series)"""
    files = sorted(glob.glob(os.path.join(workdir, "robustness", "*", "mesh", "Level0", "*", "rank_*.vtu")))
    files = [f for f in files if "/000000/" not in f]          # skip the initial-time output
    if not files or (contrast != 1.0 and pattern != "layersX"):
        return None, None
    q = -rate / VOLUME                                          # injected mass rate per unit volume
    num = den = 0.0
    for f in files:
        a = read_vtu_arrays(f, ["pressure", "elementCenter"])
        p, c = a.get("pressure", []), a.get("elementCenter", [])
        for i, pi in enumerate(p):
            x = c[3 * i]
            pe = exact_series(x, q, contrast) if contrast != 1.0 else q * MU / (2.0 * RHO * K_REF) * x * (L - x)
            num += (pi - p0 - pe) ** 2; den += pe ** 2
    return math.sqrt(num / den) if den > 0 else None, q * MU / (8.0 * RHO * K_REF)

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--geosx", default=os.environ.get("GEOSX", shutil.which("geosx") or ""), help="geosx executable")
    ap.add_argument("--template", default=os.path.join(HERE, "mfd_robustness_template.xml.in"))
    ap.add_argument("--workdir", default="mfd_robustness_runs")
    ap.add_argument("--mesh", default="hex,tet", help="hex, tet or both (comma separated)")
    ap.add_argument("--levels", default="4,8,16,32", help="cells per unit length for the h-sweep (even numbers)")
    ap.add_argument("--level", type=int, default=16, help="mesh size of the contrast sweep")
    ap.add_argument("--contrasts", default="1,1e2,1e4,1e6")
    ap.add_argument("--pattern", default="layersZ", choices=["layersZ", "layersX"], help="contrast layers parallel (Z) or in series (X) with the flow")
    ap.add_argument("--states", default="tpfa,mixed,mfd")
    ap.add_argument("--rate", type=float, default=-0.16, help="total mass rate [kg/s], negative = injection; -0.16 gives p_max - p0 = 1e6 Pa")
    ap.add_argument("--p0", type=float, default=1.0e7, help="reference pressure of the initial state and of the two Dirichlet faces")
    ap.add_argument("--mfd-percent", type=float, default=None, help="prescribe the MFD product on this percentage of the cells (a y-slab) instead of classifying; 100 = full MFD")
    ap.add_argument("--mfd-percents", default=None, help="extra sweep over prescribed MFD percentages at --level, e.g. 0,25,50,75,100")
    ap.add_argument("--table", action="store_true", help="4-D table mesh x mfd%% x n x contrast (layers in series), entries 'rel L2 error / iterations'")
    ap.add_argument("--table-percents", default="0,25,75,100"); ap.add_argument("--table-levels", default="2,4,8,16,32,64"); ap.add_argument("--table-contrasts", default="1,1e4,1e7")
    ap.add_argument("--no-output", action="store_true", help="no VTK output: iterations only, no pressure error (use for the large levels)")
    ap.add_argument("--np", type=int, default=1, help="number of MPI ranks (1 = run geosx directly unless --mpirun is given)")
    ap.add_argument("--mpirun", default=None, help="MPI launcher, options allowed, e.g. 'srun --overlap'; default: mpirun when --np > 1, "
                                                   "none when --np 1. Give it explicitly to launch the serial runs too (needed inside a Slurm step)")
    ap.add_argument("--timeout", type=int, default=3600)
    ap.add_argument("--write-only", action="store_true", help="generate the decks without running")
    a = ap.parse_args()
    # geosx runs with the case directory as working directory: resolve the paths given on the command line now
    if not a.write_only:
        geosx = os.path.abspath(a.geosx) if os.path.exists(a.geosx) else shutil.which(a.geosx or "geosx")
        if not geosx or not os.access(geosx, os.X_OK):
            sys.exit(f"geosx executable not found or not executable: '{a.geosx}' (pass --geosx with a valid path or set GEOSX)")
        a.geosx = geosx
        if a.np > 1 or a.mpirun:
            words = shlex.split(a.mpirun or "mpirun")
            exe = os.path.abspath(words[0]) if os.path.exists(words[0]) else shutil.which(words[0])
            if not exe:
                sys.exit(f"MPI launcher not found: '{words[0]}' (pass --mpirun)")
            a.mpirun = [exe] + words[1:]
    a.template = os.path.abspath(a.template)
    if not os.path.exists(a.template):
        sys.exit(f"template not found: {a.template}")

    def do(mesh, n, pattern, contrast, state, percent=None):
        tag = state if percent is None else f"mfd{percent:g}pct"
        d = os.path.join(a.workdir, f"{mesh}_n{n}_{pattern}_c{contrast:g}_{tag}")
        write_case(a.template, d, mesh, n, pattern, contrast, state, a.rate, a.p0, percent, not a.no_output)
        if a.write_only:
            return None
        r = run_case(a.geosx, d, a.timeout, a.np, a.mpirun, write_case.cells)
        r["l2"], r["pmax"] = pressure_error(d, a.rate, contrast, a.p0, pattern)
        return r

    if a.table:
        contrasts = [float(x) for x in a.table_contrasts.split(",")]
        tol = re.search(r'krylovTol="([^"]+)"', open(a.template).read())
        print(f"\nrelative L2 pressure error / GMRES iterations / linear solve time [s] (setup excluded) of the first solve;"
              f" krylovTol = {tol.group(1) if tol else '?'}")
        print(f"isotropic permeability in series along the flow (x): k = {K_REF:g} m^2 on x < 0.5, k = {K_REF:g} / contrast on x > 0.5,"
              f" interface on mesh faces; p = p0 on x = 0 and x = 1, uniform source")
        print("unknowns = cells + faces of the assembled matrix; reduced = cells + non-condensed flux dofs, the coupled system after the"
              " condensation (a condensed flux is a one-way closure row)")
        print(f"{'mesh':>4} {'mfd%':>5} {'n':>3} {'cells':>8} {'unknowns':>9} {'reduced':>9} | " + " | ".join(f"{'contrast ' + format(c, 'g'):>32}" for c in contrasts) + " | wall")
        for mesh in a.mesh.split(","):
            for percent in [float(x) for x in a.table_percents.split(",")]:
                for n in [int(x) for x in a.table_levels.split(",")]:
                    cells, walls = [], []
                    for contrast in contrasts:
                        r = do(mesh, n, "layersX", contrast, "mfd", percent)
                        if r is None:
                            continue
                        got = 100.0 * r.get("mfdCells", 0) / max(r.get("cells", 1), 1)
                        err = f"{r['l2']:.2e}" if r["l2"] is not None else "   n/a  "
                        solve = f"{r['solveTime']:9.3f}" if "solveTime" in r else "      n/a"
                        cells.append(f"{err} / {r.get('iterations', '-'):>3} / {solve}" + ("!" if r["error"] else " ") + (f"({got:.0f}%)" if abs(got - percent) > 0.5 else "    "))
                        walls.append(r["wall"])
                    if not cells:
                        continue
                    reduced = r["cells"] + r["live"] if "cells" in r and "live" in r else "-"
                    print(f"{mesh:>4} {percent:>5g} {n:>3} {r.get('cells', '-'):>8} {r.get('unknowns', '-'):>9} {reduced:>9} | "
                          + " | ".join(f"{c:>32}" for c in cells) + f" | {sum(walls):6.0f} s", flush=True)
        return

    for mesh in a.mesh.split(","):
        print(f"\n=== {mesh}: h-sweep, contrast 1 (iterations of the first solve; L2 error vs the quadratic solution, observed order) ===")
        print(f"{'state':>6} {'n':>4} {'cells':>8} {'unknowns':>9} {'mfd%':>5} {'its':>4} {'rel L2 err':>11} {'order':>6} {'s':>6}")
        states = a.states.split(",") if a.mfd_percent is None else [f"mfd{a.mfd_percent:g}pct"]
        for state in states:
            prev = None
            for n in [int(x) for x in a.levels.split(",")]:
                r = do(mesh, n, "layersZ", 1.0, "mfd", a.mfd_percent) if a.mfd_percent is not None else do(mesh, n, "layersZ", 1.0, state)
                if r is None: continue
                order = math.log(prev / r["l2"]) / math.log(2.0) if (prev and r["l2"]) else float("nan")
                prev = r["l2"]
                frac = 100.0 * r.get("mfdCells", 0) / max(r.get("cells", 1), 1)
                print(f"{state:>6} {n:>4} {r.get('cells','-'):>8} {r.get('unknowns','-'):>9} {frac:>5.0f} {r.get('iterations','-'):>4} "
                      f"{(r['l2'] if r['l2'] is not None else float('nan')):>11.3e} {order:>6.2f} {r['wall']:>6.1f}{'  ERROR' if r['error'] else ''}", flush=True)
        print(f"\n=== {mesh}: contrast sweep at n = {a.level}, pattern {a.pattern} ===")
        print(f"{'state':>6} {'contrast':>9} {'mfd%':>5} {'live':>8} {'its':>4} {'s':>6}")
        for state in states:
            for contrast in [float(x) for x in a.contrasts.split(",")]:
                r = do(mesh, a.level, a.pattern, contrast, "mfd", a.mfd_percent) if a.mfd_percent is not None else do(mesh, a.level, a.pattern, contrast, state)
                if r is None: continue
                frac = 100.0 * r.get("mfdCells", 0) / max(r.get("cells", 1), 1)
                print(f"{state:>6} {contrast:>9g} {frac:>5.0f} {r.get('live','-'):>8} {r.get('iterations','-'):>4} {r['wall']:>6.1f}{'  ERROR' if r['error'] else ''}", flush=True)
        if a.mfd_percents:
            print(f"\n=== {mesh}: prescribed MFD percentage sweep at n = {a.level}, contrast 1 (y-slab of MFD cells) ===")
            print(f"{'asked%':>7} {'mfd%':>5} {'live':>8} {'condensed':>9} {'its':>4} {'rel L2 err':>11} {'s':>6}")
            for percent in [float(x) for x in a.mfd_percents.split(",")]:
                r = do(mesh, a.level, "layersZ", 1.0, "mfd", percent)
                if r is None: continue
                frac = 100.0 * r.get("mfdCells", 0) / max(r.get("cells", 1), 1)
                print(f"{percent:>7g} {frac:>5.0f} {r.get('live','-'):>8} {r.get('condensed','-'):>9} {r.get('iterations','-'):>4} "
                      f"{(r['l2'] if r['l2'] is not None else float('nan')):>11.3e} {r['wall']:>6.1f}{'  ERROR' if r['error'] else ''}", flush=True)

if __name__ == "__main__":
    main()
