"""Exercise the complete writer without running GEOS.

Real MPI:
  python tests/run_writer_regression.py --output /tmp/pfw-check --communicator mpi
  # Choose --launcher srun on an allocated Slurm node.
Without MPI (explicitly emulated collectives, NOT a real MPI test):
  python tests/run_writer_regression.py --output /tmp/pfw-check --communicator emulated

The merged particle file must be byte-identical for 1, 3, and 7 writer ranks;
7 ranks also tests ranks with zero owned rows on this five-row coarse lattice.
"""
import argparse
import getpass
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np


def inspect_particle_file(path):
    with path.open() as stream:
        names = stream.readline().split()
    data = np.loadtxt(path, skiprows=1, ndmin=2)
    assert np.all(np.isfinite(data)), 'Nonfinite particle data'
    assert np.all(data[:,names.index('GeometryVersion')]==2)
    shape = data[:, names.index('GeometryShape')].astype(int)
    assert set(shape).issubset({0, 1, 2})
    volume = data[:, names.index('ParticleVolume')]
    columns = [names.index('Geometry'+ij) for ij in ('XX','XY','XZ','YX','YY','YZ','ZX','ZY','ZZ')]
    measure = np.linalg.det(data[:, columns].reshape(-1, 3, 3)) * np.choose(shape, [1, .5, 1/6])
    assert np.allclose(measure, volume, rtol=1e-9, atol=0)
    assert not any(name.startswith(('RVector','GeometryRVector','IntegrationWidth')) for name in names)
    assert np.all(volume > 0), 'Nonpositive particle volume'
    assert np.array_equal(data[:, 0], np.arange(1, len(data) + 1)), 'Duplicate/missing IDs'
    return len(data), float(volume.sum())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--communicator', choices=['mpi', 'emulated'], default='mpi')
    parser.add_argument('--launcher', default='mpiexec')
    parser.add_argument('--dimensions', type=int, nargs='+', default=[2, 3], choices=[2, 3])
    parser.add_argument('--scenes', nargs='+', default=['plane', 'separated'], choices=['plane', 'separated'])
    parser.add_argument('--modes', nargs='+', default=['cartesian', 'simplex'], choices=['cartesian', 'simplex'])
    parser.add_argument('--reference-geometry', choices=['equivalentBox', 'fromRefinement'], default='equivalentBox')
    parser.add_argument('--geometry-format', choices=['legacy', 'unified'], default='legacy')
    parser.add_argument('--writer', type=Path, help='Optional baseline writer for byte comparisons')
    parser.add_argument('--ranks', type=int, nargs='+', default=[1, 3, 7])
    parser.add_argument('--level', type=int, default=2)
    args = parser.parse_args()
    if min(args.ranks) < 1 or args.level < 1:
        parser.error('Ranks and refinement level must be positive')
    tests = Path(__file__).resolve().parent
    writer = args.writer.resolve() if args.writer else tests.parent / 'particleFileWriter.py'
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=True)
    if args.communicator == 'mpi' and shutil.which(args.launcher) is None:
        parser.error('MPI launcher not found; emulated mode must be requested explicitly')
    results = []
    for dimension in args.dimensions:
        for scene in args.scenes:
            for mode in args.modes:
                reference = None
                for ranks in args.ranks:
                    name = '%dd_%s_%s_r%d' % (dimension, scene, mode, ranks)
                    directory = root / name
                    # Do not reuse old collective files or stale output.
                    if directory.exists():
                        shutil.rmtree(directory)
                    directory.mkdir()
                    shutil.copyfile(tests / 'pfw_input_refinement_smoke.py', directory / 'pfw_input_case.py')
                    (directory / ('userDefs_' + getpass.getuser() + '.py')).write_text(
                        'geosPath="/bin/true"\npfwPath=%r\ndefaultBank="test"\n' % str(tests.parent))
                    env = os.environ.copy()
                    env.update(PFW_TEST_DIM=str(dimension), PFW_TEST_GEOMETRY=scene,
                               PFW_TEST_MODE=mode, PFW_TEST_LEVEL=str(args.level),
                               PFW_TEST_REFERENCE_GEOMETRY=args.reference_geometry,
                               PFW_TEST_GEOMETRY_FORMAT=args.geometry_format,
                               OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
                    handles, processes = [], []
                    try:
                        if args.communicator == 'mpi':
                            handles.append((directory / 'mpi.log').open('w'))
                            processes.append(subprocess.Popen(
                                [args.launcher, '-n', str(ranks), sys.executable, str(writer), 'pfw_input_case.py'],
                                cwd=directory, env=env, stdout=handles[-1], stderr=subprocess.STDOUT))
                        else:
                            for rank in range(ranks):
                                handles.append((directory / ('rank%d.log' % rank)).open('w'))
                                processes.append(subprocess.Popen(
                                    [sys.executable, str(tests / '_collective_test_worker.py'), str(writer),
                                     'pfw_input_case.py', str(rank), str(ranks), str(directory / 'collectives')],
                                    cwd=directory, env=env, stdout=handles[-1], stderr=subprocess.STDOUT))
                        for process in processes:
                            if process.wait(timeout=120) != 0:
                                raise RuntimeError('Writer failed; inspect logs in ' + str(directory))
                    finally:
                        for process in processes:
                            if process.poll() is None:
                                process.kill()
                        for handle in handles:
                            handle.close()
                    particle_file = directory / 'mpmParticleFile_case'
                    count, volume = inspect_particle_file(particle_file)
                    content = particle_file.read_bytes()
                    if reference is None:
                        reference = content
                    assert content == reference, 'Rank-dependent particle file: ' + name
                    xml = (directory / 'mpm_case.xml').read_text()
                    assert 'materialList="{ mat0 }"' in xml and 'materialList="{ mat1 }"' in xml
                    assert 'surfaceRefinementMethod=' not in xml, 'Writer-only option leaked into solver XML'
                    assert 'particleGeometryFormat=' not in xml, 'File format option leaked into solver XML'
                    assert 'geometryReset="FOnly"' in xml
                    assert 'particleDomainOutput=' not in xml
                    if scene == 'plane':
                        assert abs(volume - 1.0) < 2e-12
                    report_path = directory / 'mpmParticleFile_case_refinement.json'
                    if report_path.exists():
                        report = json.loads(report_path.read_text())
                        assert abs(report['particle_volume'] - volume) < 2e-12
                    else:
                        assert args.level == 1 and mode == 'cartesian', 'Missing refinement report'
                    results.append(dict(case=name, communicator=args.communicator, reference_geometry=args.reference_geometry,
                                        geometry_format=args.geometry_format, level=args.level, particles=count,
                                        volume=volume, byte_identical_to_serial=True))
                    print(name, count, 'particles; serial-identical', flush=True)
    (root / 'writer_regression_results.json').write_text(json.dumps(results, indent=2) + '\n')
    print('All %d full-writer cases passed (%s collectives).' % (len(results), args.communicator))


if __name__ == '__main__':
    main()
