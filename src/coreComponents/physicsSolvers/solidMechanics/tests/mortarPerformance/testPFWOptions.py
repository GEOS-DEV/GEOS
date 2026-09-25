"""Full PFW metadata passthrough; explicitly emulated writer collectives, not MPI.

Optionally compare particle bytes with an unpatched PFW tree. No GEOS is run.
"""
import argparse
import getpass
import json
import os
from pathlib import Path
import subprocess
import sys
import xml.etree.ElementTree as ET


def root_path():
    for parent in Path(__file__).resolve().parents:
        if (parent / 'scripts/preProcessing/materialPointMethodParticleFileWriter').is_dir():
            return parent
    raise RuntimeError('Run from a GEOS source tree')


def generate(pfw, directory, dimension, mode, ranks, options):
    directory.mkdir(parents=True, exist_ok=False)
    text = (pfw / 'tests/pfw_input_refinement_smoke.py').read_text()
    text += '\n# Solver-only performance controls must not change the particles.\n'
    text += 'pfw.update(' + repr(options) + ')\n'
    (directory / 'pfw_input_case.py').write_text(text)
    (directory / ('userDefs_' + getpass.getuser() + '.py')).write_text(
        'geosPath="/bin/true"\npfwPath=%r\ndefaultBank="test"\n' % str(pfw))
    env = dict(os.environ, PFW_TEST_DIM=str(dimension), PFW_TEST_MODE=mode,
               PFW_TEST_GEOMETRY='plane', PFW_TEST_LEVEL='2',
               PFW_TEST_REFERENCE_GEOMETRY='fromRefinement', OPENBLAS_NUM_THREADS='1',
               OMP_NUM_THREADS='1')
    processes, handles = [], []
    try:
        for rank in range(ranks):
            handle = (directory / ('rank%d.log' % rank)).open('w')
            handles.append(handle)
            command = [sys.executable, str(pfw / 'tests/_collective_test_worker.py'),
                       str(pfw / 'particleFileWriter.py'), 'pfw_input_case.py',
                       str(rank), str(ranks), str(directory / 'collectives')]
            processes.append(subprocess.Popen(command, cwd=directory, env=env,
                                               stdout=handle, stderr=subprocess.STDOUT))
        for process in processes:
            if process.wait(timeout=120) != 0:
                raise RuntimeError('PFW failed; inspect ' + str(directory))
    finally:
        for process in processes:
            if process.poll() is None:
                process.kill()
                process.wait()
        for handle in handles:
            handle.close()
    xml = ET.parse(directory / 'mpm_case.xml')
    for key, value in options.items():
        assert any(node.get(key) == str(value) for node in xml.iter()), (key, value)
    return (directory / 'mpmParticleFile_case').read_bytes()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--baseline-pfw', type=Path)
    args = parser.parse_args()
    pfw = root_path() / 'scripts/preProcessing/materialPointMethodParticleFileWriter'
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=True)
    options = dict(weakInterfaceTraceQRBackend='native', weakInterfaceTraceComponentRanks=3,
                   weakInterfaceTraceFastMarkerSearch=1, weakInterfaceTracePerformanceLog=1)
    rows = []
    for dimension in (2, 3):
        for mode in ('cartesian', 'simplex'):
            stem = '%dd_%s' % (dimension, mode)
            base = None
            if args.baseline_pfw:
                base = generate(args.baseline_pfw.resolve(), root / (stem + '_baseline'),
                                dimension, mode, 1, {})
            reference = None
            for ranks in (1, 3, 7):
                data = generate(pfw, root / (stem + '_r%d' % ranks), dimension, mode, ranks, options)
                if reference is None:
                    reference = data
                assert reference == data, 'rank-dependent PFW output'
                if base is not None:
                    assert base == data, 'performance controls changed particle data'
                rows.append(dict(dimension=dimension, refinement=mode, ranks=ranks,
                                 option_attributes_verified=True, serial_identical=True,
                                 baseline_identical=None if base is None else True,
                                 communicator='explicitly emulated; not MPI'))
                print(stem, ranks, 'PASS', flush=True)
    (root / 'results.json').write_text(json.dumps(rows, indent=2) + '\n')


if __name__ == '__main__':
    main()
