"""
Kirsch wellbore: finite elements and mixed virtual elements against the analytical solution while the
Poisson ratio approaches one half.

The shear modulus of kirschWellbore_base.xml is kept and the bulk modulus follows from each Poisson
ratio, K = 2 G (1 + nu) / (3 (1 - 2 nu)). The plane strain Kirsch stresses do not depend on the elastic
constants, so a stress error that grows with nu belongs to the discretization.

For each Poisson ratio the committed decks are copied into the work directory with the new bulk modulus,
geosx is run for both methods when their output is missing, and a profile figure is written. The summary
figure shows against nu the stress and displacement errors, max_E |sigma_h,ij(E)| / ||sigma_inf||_max for
ij in {yz, xz}, which is zero analytically (see kirschWellboreComparison.py), and the total number of
linear solver iterations. The committed decks are
never modified.

Usage, from the bin directory of a build:
    python3 ../../inputFiles/solidMechanics/polyhedral/kirschWellborePoissonSweep.py --geosx ./geosx --np 8
"""

import argparse
import os
import re
import subprocess
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kirschWellboreComparison as comparison  # noqa: E402

DECKS = {"FEM": "kirschWellbore_fem.xml", "mixed VEM": "kirschWellbore_mixedVEM.xml"}
FOLDERS = {"FEM": "fem", "mixed VEM": "mixedVEM"}


def bulkModulus(shearModulus, nu):
    return 2.0 * shearModulus * (1.0 + nu) / (3.0 * (1.0 - 2.0 * nu))


def writeDecks(caseDir, nu):
    # geosx resolves an included file against the including deck, so the path has to be absolute
    caseDir = os.path.abspath(caseDir)
    base = open(os.path.join(comparison.scriptDir, "kirschWellbore_base.xml")).read()
    shear = float(re.search(r'defaultShearModulus="([^"]+)"', base).group(1))
    K = bulkModulus(shear, nu)
    base = re.sub(r'defaultBulkModulus="[^"]+"', f'defaultBulkModulus="{K:.10e}"', base)

    os.makedirs(caseDir, exist_ok=True)
    basePath = os.path.join(caseDir, "kirschWellbore_base.xml")
    open(basePath, "w").write(base)

    decks = {}
    for label, name in DECKS.items():
        deck = open(os.path.join(comparison.scriptDir, name)).read()
        deck = deck.replace('name="./kirschWellbore_base.xml"', f'name="{basePath}"')
        decks[label] = os.path.join(caseDir, name)
        open(decks[label], "w").write(deck)
    return basePath, decks, K


def solverSummary(logPath):
    # the last linear solve of the run, as reported by the solver
    lines = [l for l in open(logPath, errors="replace") if "Linear Solver |" in l]
    if not lines:
        return "no linear solve reported", None
    status = re.search(r"Linear Solver \| (\w+)", lines[-1]).group(1)
    iterations = int(re.search(r"Iterations: (\d+)", lines[-1]).group(1))
    residual = float(re.search(r"Final Rel Res: ([0-9.eE+-]+)", lines[-1]).group(1))
    return f"{status}, {iterations} iterations, residual {residual:.2e}", status == "Success"


def run(geosx, np_, deck, outputDir, logPath):
    command = [geosx, "-i", deck, "-o", outputDir]
    if np_ > 1:
        # the wellbore mesh is generated internally, so the partition has to be given; split theta only
        command = ["mpirun", "-np", str(np_)] + command + ["-x", "1", "-y", str(np_), "-z", "1"]
    print("  running: " + " ".join(command), flush=True)
    with open(logPath, "w") as log:
        subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--nu", type=float, nargs="+", default=[0.4, 0.49, 0.499, 0.4999], help="Poisson ratios")
    parser.add_argument("--workDir", default="kirschPoissonSweep", help="decks, outputs and figures")
    parser.add_argument("--geosx", help="geosx executable, needed only for cases without output")
    parser.add_argument("--np", type=int, default=1, help="MPI ranks per run")
    parser.add_argument("--rerun", action="store_true", help="run again even when output exists")
    parser.add_argument("--theta", type=float, default=45.0, help="angle of the profile, degrees")
    parser.add_argument("--band", type=float, default=0.45, help="half width of the cell selection, degrees")
    parser.add_argument("--bins", type=int, default=30, help="log spaced radial bins for the markers")
    args = parser.parse_args()

    workDir = os.path.abspath(args.workDir)
    errors, shears, iterations, solves = {}, {}, {}, {}

    for nu in args.nu:
        caseDir = os.path.join(workDir, f"nu_{nu:g}")
        basePath, decks, K = writeDecks(caseDir, nu)
        print(f"\nnu = {nu:g}, K = {K:.4e}", flush=True)

        runs = []
        for label in comparison.METHODS:
            outputDir = os.path.join(caseDir, FOLDERS[label])
            logPath = outputDir + ".log"
            # a run that was interrupted leaves its folder behind, so only a logged solve counts as done
            finished = os.path.exists(logPath) and solverSummary(logPath)[1] is not None
            if args.rerun or not finished:
                if not args.geosx:
                    sys.exit(f"{outputDir} is missing, give --geosx to run it")
                run(os.path.abspath(args.geosx), args.np, decks[label], outputDir, logPath)
            summary, converged = solverSummary(logPath) if os.path.exists(logPath) else ("no log", None)
            solves[(nu, label)] = summary
            print(f"  {label:10s} {summary}")
            if not converged:
                print(f"  {label:10s} did not converge, its errors mix the solver and the discretization")
            runs.append((label, outputDir, logPath))

        params = comparison.readParameters(basePath, decks["FEM"])
        errors[nu], shears[nu], iterations[nu] = comparison.compareRuns(runs, params, args.theta, args.band, args.bins,
                                            os.path.join(caseDir, "kirschWellboreComparison.png"),
                                            title=rf"Kirsch wellbore, $\nu$ = {nu:g}, $\theta$ = {args.theta:g}$^\circ$")

    writeSummary(args.nu, errors, shears, iterations, solves, os.path.join(workDir, "kirschWellborePoissonSweep.png"))


def writeSummary(nus, errors, shears, iterations, solves, save):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    names = ["s_rr", "s_tt", "s_rt", "u_r", "u_t", "s_yz", "s_xz"]
    print("\nrelative L2 error and max out of plane shear / |sigma_inf| against nu")
    print(f"{'nu':>8s} {'method':10s} " + " ".join(f"{n:>9s}" for n in names) + f" {'lin its':>8s}   last linear solve")
    for nu in nus:
        for label in comparison.METHODS:
            its = iterations[nu][label]
            print(f"{nu:8g} {label:10s} " + " ".join(f"{e:9.2e}" for e in list(errors[nu][label]) + list(shears[nu][label])) +
                  f" {'-' if its is None else its:>8}   {solves[(nu, label)]}")

    fsize = 18
    cmap = plt.get_cmap("tab10")
    styles = {"FEM": ("-", "o"), "mixed VEM": ("--", "s")}
    # distance to the incompressible limit, so the four ratios spread evenly
    x = 0.5 - np.array(nus)
    fig, axes = plt.subplots(2, 2, figsize=(22, 18))
    axes = axes.ravel()

    for label in comparison.METHODS:
        line, marker = styles[label]
        table = np.array([errors[nu][label] for nu in nus])
        for c in range(3):
            axes[0].loglog(x, table[:, c], line + marker, lw=2.5, ms=10, mfc="none", mew=2, color=cmap(c),
                           label=comparison.STRESS_NAMES[c] + " " + label)
        for c in range(2):
            axes[1].loglog(x, table[:, 3 + c], line + marker, lw=2.5, ms=10, mfc="none", mew=2, color=cmap(c),
                           label=comparison.DISP_NAMES[c] + " " + label)
        shear = np.array([shears[nu][label] for nu in nus])
        for c in range(2):
            axes[2].loglog(x, shear[:, c], line + marker, lw=2.5, ms=10, mfc="none", mew=2, color=cmap(3 + c),
                           label=comparison.SHEAR_NAMES[c] + " " + label)
        its = [iterations[nu][label] for nu in nus]
        if all(i is not None for i in its):
            axes[3].semilogx(x, its, line + marker, lw=2.5, ms=10, mfc="none", mew=2, color="k", label=label)

    ylabels = ["relative L2 error", "relative L2 error", r"$\max_E |\sigma_{h,ij}(E)| \, / \, \|\sigma_\infty\|_{\max}$",
               "total linear solver iterations"]
    for a, title, ylabel in zip(axes, ["stress", "induced displacement", r"$\sigma_{h,ij}$, $ij \in \{yz, xz\}$ (exactly zero)",
                                       "linear solver"], ylabels):
        a.set_xticks(x)
        a.set_xticklabels([f"{nu:g}" for nu in nus])
        a.minorticks_off()
        a.invert_xaxis()
        a.set_xlabel(r"Poisson ratio $\nu$", size=fsize)
        a.set_ylabel(ylabel, size=fsize)
        a.set_title(title, size=fsize)
        a.grid(True, which="both", alpha=0.3)
        a.tick_params(labelsize=fsize * 0.8)
        a.legend(fontsize=fsize * 0.65, ncol=2)

    fig.suptitle("Kirsch wellbore: FEM (solid) and mixed VEM (dashed) as the Poisson ratio approaches 1/2",
                 size=fsize * 1.05)
    fig.tight_layout()
    fig.savefig(save, dpi=110)
    plt.close(fig)
    print(f"summary figure written to {save}")


if __name__ == "__main__":
    main()
