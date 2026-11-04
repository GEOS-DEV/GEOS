"""Full pressure-region PFW emission with 1/3 emulated MPI ranks.

Use --output NEW_DIR. Includes validation failures; does not run GEOS.
"""
from pathlib import Path
import argparse
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

common='startTime="0" endTime="1" startPressure="0" endPressure="10"'
def event(tag,attrs):
    return '<'+tag+' name="pressure" '+common+' '+attrs+'/>'
def validate(xml,**extra):
    options=dict(xmin=0,xmax=1,ymin=-1,ymax=1,boundaryConditionTypes=[1,0,0,0],
                 useEvents=1,mpmEventsString=xml)
    options.update(extra)
    axis.validate_options(options)

valid=[event('BoreholePressure','boreholeRadius="0.2"'),
       event('ConfiningPressure','confiningPressureBoxMin="{-1,-1,0}" confiningPressureBoxMax="{.8,1,0}"'),
       event('ConfiningPressure','confiningPressureSphereRadius=".8"'),
       event('ConfiningPressure','uniformPressure="1"')]
for xml in valid+[valid[0]+valid[1]]:validate(xml)
bad=[event('BoreholePressure',attr) for attr in ('','boreholeRadius="0"','boreholeRadius="-1"','boreholeRadius="nan"')]
bad += [event('ConfiningPressure',attr) for attr in ('','confiningPressureSphereRadius="-1"',
    'confiningPressureSphereRadius="nan"','uniformPressure="2"',
    'uniformPressure="1" confiningPressureSphereRadius="1"',
    'uniformPressure="1" confiningPressureBoxMin="{0,0,0}"',
    'confiningPressureBoxMin="{0,0,0}"',
    'confiningPressureBoxMin="{0,0}" confiningPressureBoxMax="{1,1,0}"',
    'confiningPressureBoxMin="{0,0,0}" confiningPressureBoxMax="{-1,1,0}"',
    'confiningPressureBoxMin="{0,0,0}" confiningPressureBoxMax="{1,nan,0}"')]
for xml in valid:
    bad.extend([xml.replace('endPressure="10"','endPressure="inf"'),
                xml.replace('endPressure="10"',''),xml.replace('/>',' interpType="3"/>')])
for xml in bad:
    try:validate(xml)
    except ValueError:pass
    else:raise AssertionError('Accepted invalid pressure event: '+xml)
for xml in valid:
    try:validate(xml,useEvents=0)
    except ValueError:pass
    else:raise AssertionError('Accepted disabled event processing')
print('Pressure-region input validation passed',flush=True)

source=(base/'examples/axisymmetricPressureRegions/pfw_input_axisymmetricPressureRegions.py').read_text()
for case in ('borehole','sphere','uniform'):
    reference=None
    for ranks in args.ranks:
        folder=args.output.resolve()/f'{case}_r{ranks}';folder.mkdir(parents=True,exist_ok=False)
        env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',
                 AXIS_PRESSURE_CASE='sphere' if case=='uniform' else case)
        extra=''
        if case=='uniform':
            extra="\npfw['mpmEventsString']=pfw['mpmEventsString'].replace('confiningPressureSphereRadius=\"0.0195\"','uniformPressure=\"1\"')\n"
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
        content=(folder/'mpmParticleFile_case').read_bytes()
        if reference is None:reference=content
        assert content==reference,'Writer rank count changed particle records'
        names=content.splitlines()[0].decode().split()
        particles=np.loadtxt(folder/'mpmParticleFile_case',skiprows=1,ndmin=2)
        assert np.all(np.isfinite(particles))
        volume=particles[:,names.index('ParticleVolume')].sum()
        expected=np.pi*(.02**2-.01**2)*.04 if case=='borehole' else 4*np.pi*.02**3/3
        assert np.isclose(volume,expected,rtol=1e-12 if case=='borehole' else .02)
        xml=ET.parse(folder/'mpm_case.xml');solver=xml.find('.//SolidMechanics_MPM')
        assert solver.get('formulation')=='Axisymmetric' and solver.get('useEvents')=='1'
        events=solver.find('MPMEvents');assert events is not None
        outer=events.find('ConfiningPressure');assert outer is not None
        assert outer.get('endPressure')=='20000' and outer.get('interpType')=='2'
        if case=='borehole':
            inner=events.find('BoreholePressure');assert inner is not None
            assert inner.get('boreholeRadius')=='0.0105' and inner.get('endPressure')=='100000'
            assert outer.get('confiningPressureBoxMax')=='{0.0195,1,0}'
            assert outer.get('uniformPressure') is None
        elif case=='sphere':assert outer.get('confiningPressureSphereRadius')=='0.0195'
        else:assert outer.get('uniformPressure')=='1'
        print(case,ranks,'emulated ranks:',len(particles),'particles; ring volume',volume,flush=True)
