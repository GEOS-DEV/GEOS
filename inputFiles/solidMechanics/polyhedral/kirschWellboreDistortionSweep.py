"""Kirsch wellbore on tetrahedra with a thin inner shell of thickness delta = ratio * L, L the
tangential element size at the wall, so the split tets there have radius ratio ~ ratio. Sweeps the
ratio and the stabilization option and tabulates the mixed-VEM MGR iterations."""
import argparse, os, re, subprocess, sys

DECK = """<?xml version="1.0" ?>
<Problem>
  <Mesh>
    <InternalWellbore name="mesh1" elementTypes="{{ C3D4 }}" radius="{{ {rw}, {R} }}" theta="{{ 0, 90 }}"
      zCoords="{{ {z0}, {z1} }}" nr="{{ {nr} }}" nt="{{ {nt} }}" nz="{{ {nz} }}"
      hardRadialCoords="{{ {coords} }}"
      trajectory="{{ {{ 0.0, 0.0, {z0} }}, {{ 0.0, 0.0, {z1} }} }}"
      useCartesianOuterBoundary="0" cellBlockNames="{{ cb1 }}"/>
  </Mesh>
  <ElementRegions>
    <CellElementRegion name="Omega" cellBlocks="{{ * }}" materialList="{{ rock }}"/>
  </ElementRegions>
  <Constitutive>
    <ElasticIsotropic name="rock" defaultDensity="2700" defaultBulkModulus="5.0e8" defaultShearModulus="3.0e8"/>
  </Constitutive>
  <Solvers gravityVector="{{ 0.0, 0.0, 0.0 }}">
    <SolidMechanicsMixedVEM name="wellboreSolver" discretization="mixedVEM1" targetRegions="{{ Omega }}" logLevel="1">
      <NonlinearSolverParameters newtonTol="1.0e-3" newtonMaxIter="5"/>
      <LinearSolverParameters solverType="gmres" preconditionerType="mgr" krylovTol="1.0e-5" krylovMaxIter="2000"/>
    </SolidMechanicsMixedVEM>
  </Solvers>
  <NumericalMethods>
    <MixedVEM>
      <MixedVEMDiscretization name="mixedVEM1" hybridization="0" stabilizationLength="{stab}"/>
    </MixedVEM>
  </NumericalMethods>
  <Events maxTime="1.0">
    <PeriodicEvent name="solverApplications" forceDt="1" target="/Solvers/wellboreSolver"/>{outputEvent}
  </Events>{outputs}
  <FieldSpecifications>
    <Traction name="farField" objectPath="faceManager" setNames="{{ rpos }}" tractionType="stress" scale="1.0"
      inputStress="{{ -11.25e6, -9.0e6, -15.0e6, 0.0, 0.0, 0.0 }}"/>
    <Traction name="wellLoad" objectPath="faceManager" setNames="{{ rneg }}" tractionType="normal" scale="-2.0e6"/>
    <FieldSpecification name="symmetryY" objectPath="faceManager" fieldName="displacementTrace" component="1" scale="0.0" setNames="{{ tneg }}"/>
    <FieldSpecification name="symmetryX" objectPath="faceManager" fieldName="displacementTrace" component="0" scale="0.0" setNames="{{ tpos }}"/>
    <FieldSpecification name="planeStrain" objectPath="faceManager" fieldName="displacementTrace" component="2" scale="0.0" setNames="{{ zneg, zpos }}"/>
  </FieldSpecifications>
</Problem>
"""

OUTPUT_EVENT = """
    <PeriodicEvent name="outputs" timeFrequency="1" target="/Outputs/vtkOutput"/>"""
OUTPUTS = """
  <Outputs>
    <VTK name="vtkOutput" plotFileRoot="plot"/>
  </Outputs>"""

def radialCoords(rw, R, nr, ratio, nt):
    # geometric spacing from rw to R, plus one shell of thickness ratio * L at the wall
    L = rw * (0.5 * 3.141592653589793) / nt
    q = (R / rw) ** (1.0 / nr)
    coords = [rw * q ** i for i in range(nr + 1)]
    if ratio < 1.0:
        coords.insert(1, rw + ratio * L)
    return coords, L

def solverSummary(logPath):
    lines = [l for l in open(logPath, errors="replace") if "Linear Solver |" in l]
    if not lines:
        return "n/a", "failed"
    status = re.search(r"Linear Solver \| (\w+)", lines[-1]).group(1)
    its = int(re.search(r"Iterations: (\d+)", lines[-1]).group(1))
    res = float(re.search(r"Final Rel Res: ([0-9.eE+-]+)", lines[-1]).group(1))
    return f"{its} ({res:.1e})", status

def dofCount(logPath):
    m = re.findall(r"global dofs?\D*(\d+)|Number of global dofs\D*(\d+)", open(logPath, errors="replace").read())
    return next((a or b for a, b in m), "?")

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--geosx", required=True)
    p.add_argument("--np", type=int, default=1)
    p.add_argument("--nr", type=int, default=8)
    p.add_argument("--nt", type=int, default=8)
    p.add_argument("--nz", type=int, default=2)
    p.add_argument("--height", type=float, default=0.2)
    p.add_argument("--ratios", default="1,1e-1,1e-2,1e-3,1e-4")
    p.add_argument("--stab", default="1")
    p.add_argument("--out", default="runs")
    p.add_argument("--no-vtk", action="store_true", help="skip the VTK output of the solution")
    a = p.parse_args()
    geosx = os.path.abspath(a.geosx)
    rw, R = 0.1, 5.0
    os.makedirs(a.out, exist_ok=True)
    ratios = [float(r) for r in a.ratios.split(",")]
    stabs = [int(s) for s in a.stab.split(",")]
    print(f"{'ratio':>8} | " + " | ".join(f"stab {s:>16}" for s in stabs))
    for ratio in ratios:
        row = []
        for stab in stabs:
            coords, L = radialCoords(rw, R, a.nr, ratio, a.nt)
            tag = f"r{ratio:g}_s{stab}"
            case = os.path.join(a.out, tag)
            os.makedirs(case, exist_ok=True)
            deck = os.path.join(case, "kirsch.xml")
            open(deck, "w").write(DECK.format(rw=rw, R=R, z0=-a.height / 2, z1=a.height / 2, nr=len(coords) - 1,
                                              nt=a.nt, nz=a.nz, coords=", ".join(f"{c:.12g}" for c in coords), stab=stab,
                                              outputEvent="" if a.no_vtk else OUTPUT_EVENT,
                                              outputs="" if a.no_vtk else OUTPUTS))
            log = os.path.join(case, "log.txt")
            if not os.path.exists(log):
                cmd = ([geosx] if a.np == 1 else ["mpirun", "-np", str(a.np), geosx]) + ["-i", deck, "-o", case]
                with open(log, "w") as f:
                    subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, check=False)
            s, status = solverSummary(log)
            row.append(f"{s:>16}" + ("" if status == "Success" else " !"))
        print(f"{ratio:>8g} | " + " | ".join(row), flush=True)
    print(f"tets = {6 * (a.nr + 1) * a.nt * a.nz}, wall tangential size L = {L:.3g}, dz = {a.height / a.nz:.3g}")

if __name__ == "__main__":
    main()
