"""Full PFW initialization checks, not a GEOS or real-MPI mechanics test.

Use --communicator mpi on a machine with an MPI launcher. Emulated mode uses
independent worker processes and explicitly test-only file-backed collectives.
"""
from __future__ import annotations
import argparse
import getpass
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import numpy as np
from spherical_reference import Sphere


def check_static_reference():
    for ratio in (.25,1.,4.,100.):
        model=Sphere(inclusion_E=ratio)
        inner=model.evaluate(model.a,"inclusion")
        outer=model.evaluate(model.a,"matrix")
        assert abs(inner["displacement_r"]-outer["displacement_r"]) < 1e-14
        assert abs(inner["stress_rr"]-outer["stress_rr"]) < 1e-14
        assert abs(model.evaluate(model.b)["stress_rr"]+model.pressure) < 1e-14
        if ratio==1:
            assert abs(model.coefficients()[2]) < 1e-14
            assert abs(inner["stress_rr"]+model.pressure) < 1e-14


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--communicator",choices=["mpi","emulated"],default="mpi")
    parser.add_argument("--launcher",default="mpiexec")
    parser.add_argument("--ranks",nargs="+",type=int,default=[1,3,7])
    parser.add_argument("--modes",nargs="+",default=["cartesian","simplex"])
    parser.add_argument("--loadings",nargs="+",default=["pressure","ftableCube"])
    args=parser.parse_args()
    here=Path(__file__).resolve().parent
    pfw=here.parent.parent
    worker=pfw/"tests/_collective_test_worker.py"
    if args.communicator=="mpi" and shutil.which(args.launcher) is None:
        parser.error("MPI launcher unavailable; select emulated explicitly for the test-only API")
    check_static_reference()
    output=args.output.resolve(); output.mkdir(parents=True,exist_ok=True)
    results=[]
    for loading in args.loadings:
        for mode in args.modes:
            method_arrays=[]
            for variant in ("singleField","traceContactGroups"):
                reference=None
                for ranks in args.ranks:
                    label=f"{loading}_{mode}_{variant}_r{ranks}"
                    directory=output/label
                    if directory.exists(): shutil.rmtree(directory)
                    directory.mkdir()
                    shutil.copy2(here/"pfw_input_weakTraceSphericalInclusionOctant.py",directory/"pfw_input_case.py")
                    (directory/("userDefs_"+getpass.getuser()+".py")).write_text(
                        'geosPath="/bin/true"\npfwPath=%r\ndefaultBank="test"\n' % str(pfw))
                    env=os.environ.copy()
                    env.update(WEAK_TRACE_LOADING=loading,WEAK_TRACE_VARIANT=variant,
                               WEAK_TRACE_CPP="10",WEAK_TRACE_PPC="2",WEAK_TRACE_GAP_CELLS="2",
                               WEAK_TRACE_SURFACE_REFINEMENT_LEVEL="2",WEAK_TRACE_SURFACE_REFINEMENT_METHOD=mode,
                               VV_CASE_NAME="case",OPENBLAS_NUM_THREADS="1",OMP_NUM_THREADS="1")
                    # PFW process count is independent of the requested GEOS partitioning.
                    env.update(WEAK_TRACE_XPAR="1",WEAK_TRACE_YPAR="1",WEAK_TRACE_ZPAR="1")
                    logs=[]; processes=[]
                    try:
                        commands=([[args.launcher,"-n",str(ranks),sys.executable,str(pfw/"particleFileWriter.py"),"pfw_input_case.py"]]
                                  if args.communicator=="mpi" else
                                  [[sys.executable,str(worker),str(pfw/"particleFileWriter.py"),"pfw_input_case.py",
                                    str(rank),str(ranks),str(directory/"collectives")] for rank in range(ranks)])
                        for i,command in enumerate(commands):
                            logs.append((directory/f"rank{i}.log").open("w"))
                            processes.append(subprocess.Popen(command,cwd=directory,env=env,stdout=logs[-1],stderr=subprocess.STDOUT))
                        for process in processes:
                            if process.wait(timeout=180)!=0:
                                raise RuntimeError("PFW failed: "+str(directory))
                    finally:
                        for process in processes:
                            if process.poll() is None: process.kill()
                        for stream in logs: stream.close()
                    particle_file=directory/"mpmParticleFile_case"
                    contents=particle_file.read_bytes()
                    if reference is None: reference=contents
                    assert contents==reference,"Rank-dependent initialization: "+label
                    names=particle_file.open().readline().split()
                    data=np.loadtxt(particle_file,skiprows=1,ndmin=2)
                    assert np.isfinite(data).all()
                    assert np.array_equal(data[:,0],np.arange(1,len(data)+1))
                    ij=("XX","XY","XZ","YX","YY","YZ","ZX","ZY","ZZ")
                    rv=data[:,[names.index("RVector"+s) for s in ij]].reshape(-1,3,3)
                    volume=8*np.abs(np.linalg.det(rv))
                    assert (volume>0).all()
                    shape=data[:,names.index("ParticleReferenceGeometry")].astype(int)
                    assert set(shape)<=set((0,2)),"Expected boxes and tetrahedra in 3D"
                    gv=data[:,[names.index("GeometryRVector"+s) for s in ij]].reshape(-1,3,3)
                    geometry_volume=np.linalg.det(gv)*np.choose(shape,[8.,.5,1/6])
                    assert np.allclose(geometry_volume,volume,rtol=1e-9,atol=0)
                    xyz=data[:,[names.index(s) for s in ("PositionX","PositionY","PositionZ")]]
                    assert (xyz>=-1e-12).all()
                    tetra=shape==2
                    v0=xyz[tetra]-gv[tetra].sum(axis=1)/4
                    if tetra.any():
                        verts=np.stack((v0,v0+gv[tetra,0],v0+gv[tetra,1],v0+gv[tetra,2]),axis=1)
                        assert (verts>=-2e-12).all(),"A simplex crossed a symmetry plane"
                    assert set(data[:,names.index("MaterialType")].astype(int))=={0,1}
                    flags=data[:,names.index("SurfaceFlag")].astype(int)
                    assert (flags==5).any(),"No weak-interface samples"
                    xml=(directory/"mpm_case.xml").read_text()
                    assert 'planeStrain="0"' in xml
                    assert 'weakInterfaceTraceBoundaryCoupling="1"' in xml if variant=="traceContactGroups" else 'weakInterfaceTraceBoundaryCoupling="1"' not in xml
                    if ranks==args.ranks[0]:
                        # Both methods must use precisely the same integration points,
                        # masses, flags, normal/offset data, and native shapes.
                        method_arrays.append((names,np.delete(data,names.index("ContactGroup"),axis=1)))
                    results.append(dict(case=label,particles=len(data),volume=float(volume.sum()),
                        native_tetrahedra=int(tetra.sum()),reference_rank_count=args.ranks[0],byte_identical_to_reference=True,communicator=args.communicator))
                    print(label, len(data),"particles; reference-rank-identical",flush=True)
            assert method_arrays[0][0]==method_arrays[1][0]
            assert np.array_equal(method_arrays[0][1],method_arrays[1][1]),"Method comparison changed particles"
    (output/"initialization_results.json").write_text(json.dumps(dict(
        static_reference_checks="passed, four modulus ratios",cases=results,
        scope="Full PFW only; no GEOS physics or real MPI unless communicator=mpi"),indent=2)+"\n")
    print("PASS",len(results),"writer configurations and matched single-field/mortar particles")


if __name__=="__main__": main()
