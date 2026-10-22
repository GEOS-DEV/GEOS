"""Complete hexTet writer tests; emulated collectives, never a GEOS/MPI job.

Exercise planar free surfaces, two-material interfaces, curved surfaces and
missing-plane fallback at 1/3/7 ranks. Compare the preset against explicit
simplex/native/CPTI settings and the summary against actual serialized records.
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


def run_case(tests, directory, text, ranks, expected_error=None):
    if directory.exists():
        shutil.rmtree(directory)
    directory.mkdir(parents=True)
    (directory/'pfw_input_case.py').write_text(text)
    (directory/('userDefs_'+getpass.getuser()+'.py')).write_text(
        'geosPath="/bin/true"\npfwPath=%r\ndefaultBank="test"\n' % str(tests.parent))
    env = dict(os.environ, PFW_TEST_DIM='3', PFW_TEST_GEOMETRY='plane',
               OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
    processes, handles = [], []
    try:
        for rank in range(ranks):
            handle = (directory/('rank%d.log' % rank)).open('w')
            handles.append(handle)
            processes.append(subprocess.Popen([sys.executable, str(tests/'_collective_test_worker.py'),
                str(tests.parent/'particleFileWriter.py'), 'pfw_input_case.py', str(rank),
                str(ranks), str(directory/'collectives')], cwd=directory, env=env,
                stdout=handle, stderr=subprocess.STDOUT))
        codes = [process.wait(timeout=120) for process in processes]
    finally:
        for process in processes:
            if process.poll() is None:
                process.kill()
        for handle in handles:
            handle.close()
    if expected_error:
        assert all(code != 0 for code in codes)
        for rank in range(ranks):
            assert expected_error in (directory/('rank%d.log' % rank)).read_text(), str(directory)
    elif any(codes):
        raise RuntimeError('Writer failed; inspect '+str(directory))


def inspect(directory, scene):
    path = directory/'mpmParticleFile_case'
    count, volume = inspect_particle_file(path)
    names = path.read_text().splitlines()[0].split()
    data = np.loadtxt(path, skiprows=1, ndmin=2)
    types = data[:, names.index('ParticleType')].astype(int)
    shapes = data[:, names.index('GeometryShape')].astype(int)
    volumes = data[:, names.index('ParticleVolume')]
    assert np.all(shapes[types == 2] == 0) and np.all(shapes[types == 3] == 2)
    assert set(types) == ({2} if scene in ('bulk', 'fallback') else {2, 3})
    summary = json.loads((directory/'mpmParticleFile_case_refinement.json').read_text())
    layout = summary['particle_layout']
    assert layout['total_count'] == count == summary['emitted_particles']
    assert summary['method'] == 'simplex' and summary['resolved_geometry'] == 'fromRefinement'
    assert summary['resolved_simplex_particle_type'] == 3
    assert abs(layout['total_volume']-volume) < 2e-12
    for group in layout['by_type_and_shape']:
        selected = volumes[(types == group['particle_type']) & (shapes == group['geometry_shape'])]
        assert len(selected) == group['count']
        for key, actual in (('volume', selected.sum()), ('min_particle_volume', selected.min()),
                            ('max_particle_volume', selected.max())):
            assert np.isclose(group[key], actual, rtol=1e-9, atol=1e-15), (key, group, actual)
    columns = [names.index('Geometry'+ij) for ij in ('XX','XY','XZ','YX','YY','YZ','ZX','ZY','ZZ')]
    generators = data[:, columns].reshape(-1, 3, 3)
    tet_g = generators[shapes == 2]
    # Recover all four vertices from centroid/full edges and check the physical
    # plane for accepted cuts; this also catches accidental equivalent boxes.
    points = data[:, [names.index(k) for k in ('PositionX', 'PositionY', 'PositionZ')]]
    if len(tet_g):
        v0 = points[shapes == 2] - tet_g.sum(axis=1)/4
        vertices = np.concatenate((v0[:, None, :], v0[:, None, :]+tet_g), axis=1)
        assert vertices.min() >= -1e-10 and vertices.max() <= 1+1e-10
        if scene in ('free', 'interface'):
            normal = np.array([1., .31, .19]); normal /= np.linalg.norm(normal)
            distance = vertices @ normal - .49
            mat = data[shapes == 2, names.index('MaterialType')].astype(int)
            assert np.all(distance[mat == 0] <= 1e-10)
            assert np.all(distance[mat == 1] >= -1e-10)
    if scene == 'interface':
        assert abs(volume-1) < 2e-12
    elif scene == 'free':
        # x <= .49*norm([1,.31,.19]) - .31*y - .19*z is inside [0,1]
        # throughout this unit yz square; integrate the affine height exactly.
        expected = .49*np.linalg.norm([1, .31, .19])-.25
        assert abs(volume-expected) < 2e-12
    elif scene == 'curved':
        expected = 4*np.pi*.3**3/3
        assert abs(volume-expected)/expected < .02
    elif scene == 'bulk':
        assert count == 60 and abs(volume-1) < 2e-12
        assert summary.get('candidate_domains', 0) == 0
    else:
        assert summary.get('fallback_no_usable_plane', 0) > 0
        assert not len(tet_g)
    mesh = ET.parse(directory/'mpm_case.xml').find('.//ParticleMesh')
    assert mesh is not None and 'CPDI' in mesh.attrib['particleTypes']
    assert 'CPTI' in mesh.attrib['particleTypes']  # declared even if fallback emits no tets
    xml = (directory/'mpm_case.xml').read_text()
    for option in ('surfaceRefinementMethod=', 'simplexParticleType=', 'particleGeometry='):
        assert option not in xml
    return dict(particles=count, volume=volume, layout=layout,
                fallbacks={k: v for k, v in summary.items() if k.startswith('fallback_')})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    # Collective paths must stay absolute when subprocesses change directory.
    args.output = args.output.resolve()
    tests = Path(__file__).resolve().parent
    fixture = (tests/'pfw_input_refinement_smoke.py').read_text()
    records = []
    for scene in ('bulk', 'free', 'interface', 'curved', 'fallback'):
        extra = '''
# The new preset supplies native geometry and CPTI cuts. Remove the fixture's
# deliberately explicit legacy equivalent-box setting to avoid a contradiction.
pfw.pop('particleReferenceGeometry')
pfw['surfaceRefinementMethod'] = 'hexTet'
pfw['surfaceRefinementLevel'] = 2
class SolidBox(PlanarSide):
    def isInterior(self, point, depth):
        return -1 if np.any(np.asarray(point) < 0) or np.any(np.asarray(point) >= 1) else 0
class UnflaggedSide(PlanarSide):
    def isInterior(self, point, depth):
        return min(0, super().isInterior(point, depth))
def make_objects(rankxmin, rankxmax):
    if %r == 'curved':
        return [geom.sphere('ball', [.5, .5, .5], .3, mat=0, group=0, particleType=2)]
    objects = [UnflaggedSide(1, 0)] if %r == 'fallback' else [PlanarSide(1, 0)]
    if %r == 'bulk':
        objects = [SolidBox(1, 0)]
    if %r == 'interface':
        objects.append(PlanarSide(-1, 1))
    for obj in objects:
        obj.particleType = 2
    return objects
''' % (scene, scene, scene, scene)
        # Resolve a curved object on a less anisotropic lattice; exact planar
        # tests above use only five x rows and exercise zero-row writer ranks.
        if scene == 'curved':
            extra += '\npfw.update(nI=10, nJ=10, nK=10)\n'
        text = fixture+extra
        serial = None
        for ranks in (1, 3, 7):
            directory = args.output/(scene+'_r'+str(ranks))
            run_case(tests, directory, text, ranks)
            result = inspect(directory, scene)
            content = (directory/'mpmParticleFile_case').read_bytes()
            if serial is None:
                serial = content
            assert content == serial, scene+' output changes with rank count'
            result.update(scene=scene, ranks=ranks, communicator='emulated', serial_identical=True)
            records.append(result)
            print(scene, ranks, result['particles'], result['volume'], result['fallbacks'], flush=True)
        explicit = args.output/(scene+'_explicit')
        run_case(tests, explicit, text+"\npfw.update(surfaceRefinementMethod='simplex', "
                 "particleGeometry='fromRefinement', simplexParticleType=3)\n", 1)
        inspect(explicit, scene)
        assert (explicit/'mpmParticleFile_case').read_bytes() == serial
        # Paths in comments can differ by run directory; compare parsed XML.
        def mesh_xml(directory):
            return ET.tostring(ET.parse(directory/'mpm_case.xml').find('.//ParticleMesh'))
        assert mesh_xml(explicit) == mesh_xml(args.output/(scene+'_r1'))
        if scene == 'fallback':
            run_case(tests, args.output/'strict_fallback', text+"\npfw['surfaceRefinementFallback']='error'\n",
                     1, 'no_usable_plane')
        if scene == 'interface':
            for name, setting, error in (
                ('box', "pfw['particleGeometry']='box'", 'hexTet requires'),
                ('strain', "pfw['planeStrain']=1\npfw['nK']=3", 'requires planeStrain=0'),
                ('bulk_tets', "_factory=make_objects\ndef make_objects(a,b):\n objects=_factory(a,b)\n for obj in objects: obj.particleType=3\n return objects", 'hexTet bulk objects')):
                run_case(tests, args.output/name, text+'\n'+setting+'\n', 1, error)
            # Only a non-root rank sees this invalid object locally. Global
            # metadata validation must reject it on every rank before writing.
            hidden_bad_type = """
_local_factory = make_objects
def make_objects(a, b):
    if a < .3:
        return []
    objects = _local_factory(a, b)
    for obj in objects:
        obj.particleType = 3
    return objects
"""
            run_case(tests, args.output/'rank_local_bad_type', text+hidden_bad_type, 3,
                     'hexTet bulk objects')
    (args.output/'results.json').write_text(json.dumps(records, indent=2)+'\n')
    print('15 preset runs, 5 explicit-equivalence runs and 5 rejection cases passed.', flush=True)


if __name__ == '__main__':
    main()
