"""Generate complete axisymmetric CZ inputs using 1/3 emulated MPI ranks."""
import argparse
import getpass
import os
from pathlib import Path
import subprocess
import sys
import numpy as np

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output', type=Path, required=True)
parser.add_argument('--ranks', type=int, nargs='+', default=[1,3])
args=parser.parse_args()
base=Path(__file__).resolve().parents[1]
case=base/'examples/axisymmetricCZ'
sys.path.insert(0,str(base))
import pfw_axisymmetric as axis

data=dict(particleType=2,group=0,surfaceFlag=3,czTag=1,velocity=[0,0,0],
          surfaceNormal=[0,1,0],surfacePosition=[0,.5,0])
for group in (0,1):
    for flag in (3,4):
        axis.validate_particle([2,0,0],dict(data,group=group,surfaceFlag=flag),(3,np.eye(2)))
for change in (dict(group=2),dict(czTag=-1),dict(czTag=.3),dict(surfaceNormal=[0,0,1]),
               dict(surfaceNormal=[0,0,0]),dict(surfacePosition=[0,0,1])):
    try: axis.validate_particle([2,0,0],dict(data,**change),(3,np.eye(2)))
    except ValueError: pass
    else: raise AssertionError(change)
options=dict(xmin=0,xmax=3,useEvents=1,
             mpmEventsString='<ReferenceCohesiveZones computeNormalsAndPositions="0"/>',
             cohesiveZoneRegions='<CohesiveZoneRegion name="b" constitutiveModel="law" tag="1"/>',
             particleFileFields=['ContactGroup','SurfaceFlag','SurfaceNormal','SurfacePosition','CZTag'])
for model in ('UncoupledCohesiveZone','CoupledCohesiveZone','BicrystalCohesiveZone'):
    axis.validate_options(dict(options,materialPropertyString=f'<{model} name="law"/>'))
axis.validate_options(dict(options,useEvents='1'))
# Match GEOS's default and disabled event dispatcher: valid surface flags and
# a declared event cannot create a bond unless useEvents is explicitly one.
for event_switch in (None,0,'0',2):
    invalid=dict(options)
    if event_switch is None: invalid.pop('useEvents')
    else: invalid['useEvents']=event_switch
    try: axis.validate_options(invalid)
    except ValueError as error: assert 'useEvents=1' in str(error)
    else: raise AssertionError(('Disabled CZ activation accepted',event_switch))
for change in (dict(materialPropertyString='<SurfaceInformedPolymerCohesiveZone name="law"/>'),
               dict(materialPropertyString='<PolymerCohesiveZone name="law"/>'),
               dict(particleFileFields=[]),
               dict(mpmEventsString='<ReferenceCohesiveZones computeNormalsAndPositions="1"/>')):
    try: axis.validate_options(dict(options,**change))
    except ValueError: pass
    else: raise AssertionError(change)
print('Axisymmetric CZ input validation passed.',flush=True)

for mode in ('elastic','debond'):
    for update in ('FLIP','FMPM'):
        reference=None
        for size in args.ranks:
            folder=args.output.resolve()/f'{mode}_{update}_r{size}'
            folder.mkdir(parents=True,exist_ok=False)
            env=dict(os.environ,AXIS_CZ_CASE=mode,AXIS_CZ_UPDATE=update,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
            (folder/f'userDefs_{getpass.getuser()}.py').write_text('geosPath="/bin/true"\npfwPath='+repr(str(base))+'\ndefaultBank="test"\n')
            subprocess.run([sys.executable,str(case/'caseSettings.py'),'--write','case.json'],cwd=folder,env=env,check=True)
            logs=[];jobs=[]
            try:
                for rank in range(size):
                    log=(folder/f'rank{rank}.log').open('w');logs.append(log)
                    jobs.append(subprocess.Popen([sys.executable,str(base/'tests/_collective_test_worker.py'),
                        str(base/'particleFileWriter.py'),str(case/'pfw_input_axisymmetricCZ.py'),str(rank),str(size),str(folder/'collectives')],
                        cwd=folder,env=env,stdout=log,stderr=subprocess.STDOUT))
                for job in jobs:
                    if job.wait(timeout=120): raise RuntimeError('Writer failed: '+str(folder))
            finally:
                for job in jobs:
                    if job.poll() is None: job.kill();job.wait()
                for log in logs: log.close()
            subprocess.run([sys.executable,str(case/'checkParticles.py')],cwd=folder,env=env,check=True)
            content=(folder/'mpmParticleFile_axisymmetricCZ').read_bytes()
            if reference is None: reference=content
            assert content==reference,'Writer rank count changes the particle file'
            print(mode,update,size,'emulated ranks passed',flush=True)
