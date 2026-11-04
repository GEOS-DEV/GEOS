"""Full writer with explicitly emulated collectives, not an MPI runtime test.

Usage: python tests/run_axisymmetric_writer_tests.py --output /tmp/axisymmetric-pfw
Checks ordinary CPDI/CPTI emission and two-pass CPDI bulk/CPTI curved surfaces.
"""
from pathlib import Path
import argparse
import getpass
import json
import os
import shutil
import subprocess
import sys
import numpy as np

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output',type=Path,required=True)
parser.add_argument('--ranks',type=int,nargs='+',default=[1,3,7])
args=parser.parse_args()
b=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(b))
import pfw_axisymmetric as axis
axis.validate_particle([2,0,0],dict(particleType=1),(3,np.eye(2)))
try:
    axis.validate_particle([2,0,0],dict(particleType=1,czTag=1),(3,np.eye(2)))
except ValueError:
    pass
else:
    raise AssertionError('Q1 cohesive area operator must not accept tagged B-spline particles')

for case in ('cylinder','refined','bspline','bspline_fmpm'):
    reference=None
    for ranks in args.ranks:
        d=args.output.resolve()/f'{case}_r{ranks}'
        d.mkdir(parents=True,exist_ok=False)
        source=(b/'examples/axisymmetricCylinder/pfw_input_axisymmetricCylinder.py').read_text()
        if case=='refined':
            source+='''\npfw.update(surfaceRefinementMethod='simplex', surfaceRefinementLevel=2, simplexParticleType=3)
pfw['objects']=[geom.disk('torus',[0.6,0.6],0.31,mat=0,group=0,particleType=2)]
'''
        if case.startswith('bspline'):
            source += '''\nfor obj in pfw['objects']: obj.particleType=1
pfw.update(surfaceDetection=1, useEvents=1,
    mpmEventsString='<FrictionCoefficientSwap name="swap" startTime="0.001" frictionCoefficient="0.25"/>')
'''
            if case=='bspline_fmpm': source += "pfw.update(updateMethod='FMPM', updateOrder=2)\n"
        (d/'pfw_input_case.py').write_text(source)
        (d/f'userDefs_{getpass.getuser()}.py').write_text('geosPath="/bin/true"\npfwPath='+repr(str(b))+'\ndefaultBank="test"\n')
        files=[];jobs=[]
        try:
            for rank in range(ranks):
                log=(d/f'rank{rank}.log').open('w');files.append(log)
                jobs.append(subprocess.Popen([sys.executable,str(b/'tests/_collective_test_worker.py'),
                    str(b/'particleFileWriter.py'),'pfw_input_case.py',str(rank),str(ranks),str(d/'collectives')],
                    cwd=d,env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1'),
                    stdout=log,stderr=subprocess.STDOUT))
            for job in jobs:
                assert job.wait(timeout=120)==0, 'Inspect writer logs in '+str(d)
        finally:
            for job in jobs:
                if job.poll() is None:
                    job.kill();job.wait()
            for log in files:log.close()
        path=d/'mpmParticleFile_case'
        with path.open() as stream:names=stream.readline().split()
        a=np.loadtxt(path,skiprows=1,ndmin=2)
        radius=a[:,names.index('PositionX')];volume=a[:,names.index('ParticleVolume')]
        shape=a[:,names.index('GeometryShape')].astype(int)
        g=a[:,[names.index('Geometry'+ij) for ij in ('XX','XY','YX','YY')]].reshape(-1,2,2)
        expected=2*np.pi*radius*np.linalg.det(g)*np.where(shape==4,0.5,1.0)
        assert np.all(np.isfinite(a)) and np.all(radius>0) and np.all(volume>0)
        assert np.allclose(volume,expected,rtol=1e-12,atol=0)
        assert set(shape)==({3} if case.startswith('bspline') else {3,4})
        if case.startswith('bspline'): assert set(a[:,names.index('ParticleType')].astype(int))=={1}
        content=path.read_bytes()
        if reference is None:reference=content
        assert reference==content
        xml=(d/'mpm_case.xml').read_text()
        assert 'axisymmetric="1"' in xml and 'formulation="Axisymmetric"' in xml
        assert 'planeStrainThickness=' not in xml
        assert 'timeIntegrationOption="ExplicitDynamic"' in xml
        assert ('updateMethod="FMPM"' if case=='bspline_fmpm' else 'updateMethod="FLIP"') in xml
        if case.startswith('bspline'):
            assert 'FrictionCoefficientSwap' in xml and 'surfaceDetection="1"' in xml
        if case!='refined':assert np.isclose(volume.sum(),np.pi,rtol=1e-12)
        else:
            report=json.loads((d/'mpmParticleFile_case_refinement.json').read_text())
            assert np.isclose(report['particle_volume'],volume.sum(),rtol=1e-12)
            assert np.isclose(report['particle_layout']['total_volume'],volume.sum(),rtol=1e-12)
            assert np.isclose(volume.sum(),2*np.pi**2*0.6*0.31**2,rtol=.02)
        print(case,ranks,'emulated ranks:',len(a),'particles, full-ring volume',volume.sum(),flush=True)
