"""Full native-2D PFW tests with emulated process collectives, never GEOS jobs."""
import argparse
import getpass
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import xml.etree.ElementTree as ET
import numpy as np


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    tests=Path(__file__).resolve().parent;writer=tests.parent/'particleFileWriter.py'
    output=args.output.resolve();output.mkdir(parents=True,exist_ok=True)
    fixture=(tests/'pfw_input_refinement_smoke.py').read_text()
    reports=[]
    for case in ('cpdi','bspline','cpti','mixed','legacy'):
      serial=None
      for ranks in (1,3,7):
        directory=output/(case+'_r'+str(ranks))
        if directory.exists():shutil.rmtree(directory)
        directory.mkdir()
        bulk=1 if case=='bspline' else (3 if case=='cpti' else 2)
        extra="""
pfw['particleGeometry']='fromRefinement'
pfw['surfaceRefinementMethod']=%r
pfw['surfaceRefinementLevel']=1
pfw['planeStrainThickness']=0.07
pfw['nK']=999 # Compatibility input is ignored; it creates no z cells.
%s
%s
original_factory=make_objects
def make_objects(rankxmin,rankxmax):
    objects=original_factory(rankxmin,rankxmax)
    for obj in objects:obj.particleType=%d
    return objects
""" % ('simplex' if case in ('mixed','cpti') else 'cartesian',
       "pfw['simplexParticleType']=3" if case=='mixed' else '',
       "pfw.pop('planeStrain');pfw['formulation']='PlaneStrain'" if case!='legacy' else "pfw.pop('planeStrainThickness');pfw['zmax']=0.07",
       bulk)
        (directory/'pfw_input_case.py').write_text(fixture+extra)
        (directory/('userDefs_'+getpass.getuser()+'.py')).write_text('geosPath="/bin/true"\npfwPath=%r\ndefaultBank="test"\n'%str(tests.parent))
        env=dict(os.environ,PFW_TEST_DIM='2',PFW_TEST_GEOMETRY='plane',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
        processes=[];handles=[]
        try:
          for rank in range(ranks):
            log=(directory/('rank%d.log'%rank)).open('w');handles.append(log)
            processes.append(subprocess.Popen([sys.executable,str(tests/'_collective_test_worker.py'),str(writer),'pfw_input_case.py',str(rank),str(ranks),str(directory/'collectives')],cwd=directory,env=env,stdout=log,stderr=subprocess.STDOUT))
          for process in processes:
            if process.wait(timeout=120):raise RuntimeError('Writer failed: '+str(directory))
        finally:
          for process in processes:
            if process.poll() is None:process.kill()
          for log in handles:log.close()
        particle=directory/'mpmParticleFile_case'; content=particle.read_bytes()
        names=particle.read_text().splitlines()[0].split();data=np.loadtxt(particle,skiprows=1,ndmin=2)
        index={key:i for i,key in enumerate(names)}
        shape=data[:,index['GeometryShape']].astype(int);ptype=data[:,index['ParticleType']].astype(int)
        assert np.all(data[:,index['PositionZ']]==0)
        assert not any(name in names for name in ('GeometryXZ','GeometryYZ','GeometryZX','GeometryZY','GeometryZZ'))
        g=data[:,[index['Geometry'+ij] for ij in ('XX','XY','YX','YY')]].reshape(-1,2,2)
        area=np.linalg.det(g)*np.where(shape==4,0.5,1)
        volume=data[:,index['ParticleVolume']]
        np.testing.assert_allclose(volume,area*0.07,rtol=2e-12,atol=0)
        assert np.all(area>0) and abs(volume.sum()-0.07)<2e-12
        assert np.all(shape[ptype==2]==3) and np.all(shape[ptype==3]==4)
        assert set(shape)<= {3,4}
        if case=='mixed':assert set(shape)=={3,4} and set(ptype)=={2,3}
        xml=ET.parse(directory/'mpm_case.xml');mesh=xml.find('.//InternalMesh');solver=xml.find('.//SolidMechanics_MPM')
        assert mesh.attrib['spatialDimension']=='2' and float(mesh.attrib['planeStrainThickness'])==0.07
        assert 'C2D4' in mesh.attrib['elementTypes'] and 'zCoords' not in mesh.attrib and 'nz' not in mesh.attrib
        assert solver.attrib['formulation']=='PlaneStrain' and 'planeStrain' not in solver.attrib
        reportfile=directory/'mpmParticleFile_case_refinement.json'
        if reportfile.exists():
          report=json.loads(reportfile.read_text())
          assert report['emitted_particles']==len(data)
          assert {r['geometry_shape'] for r in report['particle_layout']['by_type_and_shape']}<= {3,4}
        if serial is None:serial=content
        assert content==serial,'Native 2D output depends on rank count'
        record=dict(case=case,ranks=ranks,particles=len(data),volume=float(volume.sum()),communicator='emulated',serial_identical=True)
        reports.append(record);print(record,flush=True)
    (output/'results.json').write_text(json.dumps(reports,indent=2)+'\n')

if __name__=='__main__':main()
