"""Full surface-example PFW/XML emission; emulated 1/3-rank collectives, no GEOS run."""
from pathlib import Path
import argparse
import getpass
import os
import subprocess
import sys
import xml.etree.ElementTree as ET
import numpy as np
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
p.add_argument('--ranks',type=int,nargs='+',default=[1,3]);args=p.parse_args()
base=Path(__file__).resolve().parents[1];sys.path.insert(0,str(base))
import pfw_axisymmetric as axis
common=dict(xmin=0,xmax=2,computeParticleSurfaceNormalsAndPositions=1,computeNodalArea=1,enableSurfaceTension=1)
for pairs in ([[0,-1,1]],[[0,1,.1],[1,-1,0]],'{{0,-1,1}}'):
    axis.validate_options(dict(common,surfaceTensionPairs=pairs))
for pairs in ([],[[0,0,1]],[[.5,-1,1]],[[0,-2,1]],[[0,-1,-1]],[[0,-1,float('nan')]],[[0,-1]],[[0,float('inf'),1]]):
    try:axis.validate_options(dict(common,surfaceTensionPairs=pairs))
    except ValueError:pass
    else:raise AssertionError(pairs)
source=(base/'examples/axisymmetricSurfaces/pfw_input_axisymmetricSurfaces.py').read_text()
for case,interpolation in [('cylinder','mixed'),('cylinder','cpdi'),('cylinder','cpti'),('cylinder','bspline'),('sphere','mixed'),('sphere','bspline')]:
    reference=None
    for ranks in args.ranks:
        folder=args.output.resolve()/f'{case}_{interpolation}_r{ranks}';folder.mkdir(parents=True,exist_ok=False)
        env=dict(os.environ,AXIS_SURFACE_CASE=case,AXIS_SURFACE_INTERPOLATION=interpolation,
                 AXIS_SURFACE_RESOLUTION='12',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
        (folder/'pfw_input_case.py').write_text(source)
        (folder/f'userDefs_{getpass.getuser()}.py').write_text('geosPath="/bin/true"\npfwPath='+repr(str(base))+'\ndefaultBank="test"\n')
        logs=[];jobs=[]
        try:
            for rank in range(ranks):
                log=(folder/f'rank{rank}.log').open('w');logs.append(log)
                jobs.append(subprocess.Popen([sys.executable,str(base/'tests/_collective_test_worker.py'),
                    str(base/'particleFileWriter.py'),'pfw_input_case.py',str(rank),str(ranks),str(folder/'collectives')],
                    cwd=folder,env=env,stdout=log,stderr=subprocess.STDOUT))
            for job in jobs:assert job.wait(timeout=180)==0,'Inspect logs: '+str(folder)
        finally:
            for job in jobs:
                if job.poll() is None:job.kill();job.wait()
            for log in logs:log.close()
        content=(folder/'mpmParticleFile_case').read_bytes()
        if reference is None:reference=content
        assert content==reference,'Particle emission depends on writer rank count'
        names=content.splitlines()[0].decode().split();data=np.loadtxt(folder/'mpmParticleFile_case',skiprows=1,ndmin=2)
        assert np.all(np.isfinite(data));assert np.all(data[:,names.index('ParticleVolume')]>0)
        xml=ET.parse(folder/'mpm_case.xml');solver=xml.find('.//SolidMechanics_MPM')
        for key in ('computeParticleSurfaceNormalsAndPositions','computeNodalArea','enableSurfaceTension'):
            assert solver.get(key)=='1',(key,solver.attrib)
        assert solver.get('normalAndPositionMethod')=='DFGAndVolumeIntegration'
        assert solver.get('overlapCorrection')=='Off'
        assert solver.find('.//InitializeStress').get('pressure')==str(.01 if case=='cylinder' else .02)
        subprocess.run([sys.executable,str(base/'examples/axisymmetricSurfaces/postProcess.py')],cwd=folder,env=env,check=True)
        print(case,interpolation,ranks,'emulated ranks:',len(data),'particles',flush=True)
