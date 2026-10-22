"""Full PFW CPTI tests using explicitly emulated collectives (not real MPI).

Tests all-CPTI Cartesian/simplex generation and mixed CPDI/CPTI output with
1, 3 and 7 ranks, including zero-row ranks. Never submits a GEOS job.
"""
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
from run_writer_regression import inspect_particle_file


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    tests = Path(__file__).resolve().parent
    writer = tests.parent/'particleFileWriter.py'
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    fixture = (tests/'pfw_input_refinement_smoke.py').read_text()
    records = []
    for case in ('mixed', 'all_cartesian', 'all_simplex'):
        serial = None
        for ranks in (1, 3, 7):
            directory = output/(case+'_r'+str(ranks))
            if directory.exists():
                shutil.rmtree(directory)
            directory.mkdir()
            extra = '''
pfw['planeStrain'] = 0
pfw['particleGeometry'] = 'fromRefinement'
pfw['surfaceRefinementMethod'] = %r
pfw['surfaceRefinementLevel'] = 1
pfw['cpdiDomainScaling'] = 1
pfw['geometryReset'] = 'FOnly'
%s
original_factory = make_objects
def make_objects(rankxmin, rankxmax):
    objects = original_factory(rankxmin, rankxmax)
    for obj in objects:
        obj.particleType = %d
    return objects
''' % ('cartesian' if case == 'all_cartesian' else 'simplex',
       "pfw['simplexParticleType'] = 3" if case == 'mixed' else '', 2 if case == 'mixed' else 3)
            (directory/'pfw_input_case.py').write_text(fixture+extra)
            (directory/('userDefs_'+getpass.getuser()+'.py')).write_text(
                'geosPath="/bin/true"\npfwPath=%r\ndefaultBank="test"\n' % str(tests.parent))
            env = dict(os.environ, PFW_TEST_DIM='3', PFW_TEST_GEOMETRY='plane',
                       OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
            processes = []
            handles = []
            try:
                for rank in range(ranks):
                    handle = (directory/('rank%d.log' % rank)).open('w')
                    handles.append(handle)
                    processes.append(subprocess.Popen([sys.executable, str(tests/'_collective_test_worker.py'),
                        str(writer), 'pfw_input_case.py', str(rank), str(ranks), str(directory/'collectives')],
                        cwd=directory, env=env, stdout=handle, stderr=subprocess.STDOUT))
                for process in processes:
                    if process.wait(timeout=120):
                        raise RuntimeError('Writer failed: '+str(directory))
            finally:
                for process in processes:
                    if process.poll() is None:
                        process.kill()
                for handle in handles:
                    handle.close()
            particle_file = directory/'mpmParticleFile_case'
            count, volume = inspect_particle_file(particle_file)
            names = particle_file.read_text().splitlines()[0].split()
            data = np.loadtxt(particle_file, skiprows=1, ndmin=2)
            types = data[:, names.index('ParticleType')].astype(int)
            shapes = data[:, names.index('GeometryShape')].astype(int)
            assert np.all(shapes[types == 3] == 2)
            assert np.all(shapes[types == 2] == 0)
            assert set(types) == ({2, 3} if case == 'mixed' else {3})
            assert abs(volume-1.) < 2e-12
            xml = ET.parse(directory/'mpm_case.xml')
            mesh = xml.find('.//ParticleMesh')
            assert mesh is not None
            assert 'CPTI' in mesh.attrib['particleTypes']
            assert ('CPDI' in mesh.attrib['particleTypes']) == (case == 'mixed')
            assert 'simplexParticleType=' not in (directory/'mpm_case.xml').read_text()
            report_file = directory/'mpmParticleFile_case_refinement.json'
            if report_file.exists():
                report = json.loads(report_file.read_text())
                assert report['emitted_particles'] == count
                assert report['emitted_refinement_pieces'] <= count
            content = particle_file.read_bytes()
            if serial is None:
                serial = content
            assert content == serial, 'CPTI output depends on rank count'
            record = dict(case=case, ranks=ranks, particles=count, volume=volume,
                          communicator='emulated', serial_identical=True)
            records.append(record)
            print(record, flush=True)
    # Native 2D CPTI is covered by run_true_2d_writer_tests.py. Reject actual conflicts.
    for name, settings, message in (
        ('box_geometry', "pfw['particleGeometry']='box'\n", 'CPTI requires particleGeometry'),
        ('invalid_override', "pfw['simplexParticleType']=2\n", 'simplexParticleType must be'),
        ('override_without_simplex', "pfw['simplexParticleType']=3\npfw['surfaceRefinementMethod']='cartesian'\n", 'simplexParticleType requires')):
        directory = output/name
        directory.mkdir(exist_ok=True)
        (directory/'pfw_input_case.py').write_text(fixture+extra+settings)
        (directory/('userDefs_'+getpass.getuser()+'.py')).write_text(
            'geosPath="/bin/true"\npfwPath=%r\ndefaultBank="test"\n' % str(tests.parent))
        result = subprocess.run([sys.executable, str(tests/'_collective_test_worker.py'),
            str(writer), 'pfw_input_case.py', '0', '1', str(directory/'collectives')],
            cwd=directory, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=120)
        (directory/'writer.log').write_text(result.stdout)
        assert result.returncode != 0 and message in result.stdout, name
        print(name+': rejected as expected', flush=True)
    (output/'results.json').write_text(json.dumps(records, indent=2)+'\n')


if __name__ == '__main__':
    main()
