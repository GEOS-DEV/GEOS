"""Full PFW emission for axisymmetric grid controls, with emulated MPI collectives.

Use --output NEW_DIR. This validates generated inputs and files, not GEOS runtime.
"""
from pathlib import Path
import argparse
import copy
import getpass
import os
import subprocess
import sys
import xml.etree.ElementTree as ET
import numpy as np

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output',required=True,type=Path)
parser.add_argument('--ranks',type=int,nargs='+',default=[1,3])
args=parser.parse_args()
base=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(base))
import pfw_axisymmetric as axis

valid=dict(xmin=0,xmax=1,ymin=0,ymax=1,boundaryConditionTypes=[1,2,1,2],
           prescribedBoundaryFTable=1,fTable=[[0,1,1,1],[1,.9,.8,1]])
for change in ({},dict(updateMethod='FMPM'),dict(stressControl=[1,0,0],stressTable=[[0,-100,0,0]]),
               dict(stressControl=[1,1,0],stressTable=[[0,-100,-50,0]],stressControlAdaptiveSeek=1),
               dict(xmin=1,xmax=2,boundaryConditionTypes=[3,3,1,2]),
               dict(fTable='{{0,1,1,1},{1,.9,.8,1}}'),
               dict(boundaryConditionTypes=None,prescribedBcTable=1,bcTable=[[0,1,2,1,2],[1,1,0,1,0]]),
               dict(prescribedBcTable=1,bcTable=[[0,1,2,1,2],[1,1,3,1,2]])):
    axis.validate_options(copy.deepcopy(dict(valid,**change)))
for change in (dict(prescribedFTable=1),dict(stressControl=[0,0,1]),
               dict(stressControl=[1,0,0],stressTable=[[0,0,0]]),
               dict(stressControl=[1,0,0],stressTable=[[1,0,0,0]]),
               dict(stressControl=[1,0,0],stressTable=[[0,0,0,0],[0,0,0,0]]),
               dict(stressControl=[2,0,0],stressTable=[[0,0,0,0]],stressControlAdaptiveSeek=1),
               dict(fTable=[[0,1,1,1],[1,.9,1,.9]]),dict(fTable=[[0,1,1,1],[1,float('nan'),1,1]]),
               dict(boundaryConditionTypes=[1,1,1,2]),dict(boundaryConditionTypes=[1,0,1,2]),
               dict(updateMethod='FMPM',boundaryConditionTypes=[1,3,1,2]),dict(periodic=[False,True]),
               dict(prescribedBcTable=1,bcTable=[[0,1,2,1,2],[1,0,2,1,2]]),
               dict(prescribedBcTable=1,bcTable=[[0,1,2,1,2],[1,1,3,1,2]],updateMethod='FMPM'),
               dict(enablePrescribedBoundaryTransverseVelocities=[0,1,0,0,0,0],
                    prescribedBoundaryTransverseVelocities=[[0,0],[0,1],[0,0],[0,0],[0,0],[0,0]])):
    try:axis.validate_options(copy.deepcopy(dict(valid,**change)))
    except ValueError:pass
    else:raise AssertionError('Unsupported control accepted: '+repr(change))
print('Driving input validation passed',flush=True)
source=(base/'examples/axisymmetricGridDriving/pfw_input_axisymmetricGridDriving.py').read_text()
for mode in ('stretch','stress','mixed','fmpm','contact','bc_switch','adaptive','limit'):
    reference=None
    for ranks in args.ranks:
        folder=args.output.resolve()/f'{mode}_r{ranks}';folder.mkdir(parents=True,exist_ok=False)
        extra=''
        env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',
                 AXIS_GRID_MODE=mode if mode in ('stretch','stress','mixed') else 'mixed',
                 AXIS_GRID_WALL='contact' if mode=='contact' else 'moving',
                 AXIS_GRID_UPDATE='FMPM' if mode=='fmpm' else 'FLIP')
        if mode=='bc_switch':extra="\npfw.update(prescribedBcTable=1,bcTable=[[0,1,2,1,2],[.1,1,3,1,2]])\n"
        if mode=='adaptive':extra="\npfw.update(stressControlAdaptiveSeek=1,stressControlMinStrainRate=.05)\n"
        if mode=='limit':extra="\npfw.update(stressControl=[2,0,0],fTable=[[0,1,1,1],[.2,.99,.995,1]])\n"
        (folder/'pfw_input_case.py').write_text(source+extra)
        (folder/f'userDefs_{getpass.getuser()}.py').write_text('geosPath="/bin/true"\npfwPath='+repr(str(base))+'\ndefaultBank="test"\n')
        logs=[];jobs=[]
        try:
            for rank in range(ranks):
                log=(folder/f'rank{rank}.log').open('w');logs.append(log)
                jobs.append(subprocess.Popen([sys.executable,str(base/'tests/_collective_test_worker.py'),
                    str(base/'particleFileWriter.py'),'pfw_input_case.py',str(rank),str(ranks),str(folder/'collectives')],
                    cwd=folder,env=env,stdout=log,stderr=subprocess.STDOUT))
            for job in jobs:
                assert job.wait(timeout=120)==0,'Inspect writer logs: '+str(folder)
        finally:
            for job in jobs:
                if job.poll() is None:job.kill();job.wait()
            for log in logs:log.close()
        data=(folder/'mpmParticleFile_case').read_bytes()
        if reference is None:reference=data
        assert data==reference
        names=data.splitlines()[0].decode().split();particles=np.loadtxt(folder/'mpmParticleFile_case',skiprows=1,ndmin=2)
        assert np.all(np.isfinite(particles))
        assert np.isclose(particles[:,names.index('ParticleVolume')].sum(),np.pi,rtol=1e-12)
        xml=ET.parse(folder/'mpm_case.xml');solver=xml.find('.//SolidMechanics_MPM')
        assert solver.get('formulation')=='Axisymmetric' and solver.get('prescribedBoundaryFTable')=='1'
        assert solver.get('fTableInterpType')=='Smoothstep'
        assert len(solver.get('boundaryConditionTypes').strip('{} ').split(',')) in (4,6)
        assert solver.get('updateMethod')==('FMPM' if mode=='fmpm' else 'FLIP')
        if mode!='stretch':assert solver.get('stressControl') is not None
        if mode=='bc_switch':assert solver.get('prescribedBcTable')=='1'
        if mode=='adaptive':assert solver.get('stressControlAdaptiveSeek')=='1'
        print(mode,ranks,'emulated ranks:',len(particles),'particles',flush=True)
