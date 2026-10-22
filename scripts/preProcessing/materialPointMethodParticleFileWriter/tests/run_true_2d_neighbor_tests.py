"""Run the complete production neighbor routine with serial framework adapters.

Both CPU and GEOS_USE_DEVICE preprocessor paths are exercised serially. This
checks list construction and bounds, not real GPU memory or MPI communication.
An optional --source path permits reproduction against an earlier source file.
"""
from pathlib import Path
import argparse
import subprocess
import tempfile

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--source', type=Path)
args = parser.parse_args()
tests = Path(__file__).resolve().parent
root = next(p for p in tests.parents if (p / 'src/coreComponents').exists())
source = args.source or root / 'src/coreComponents/physicsSolvers/solidMechanics/SolidMechanicsMPM.cpp'
text = source.read_text()
start = text.index('real64 SolidMechanicsMPM::computeNeighborList( ParticleManager & particleManager )')
end, depth = text.index('{', start) + 1, 1
while depth:
    depth += (text[end] == '{') - (text[end] == '}')
    end += 1
program = (tests / 'neighbor_cpp/neighbor_test_support.hpp').read_text() + '\n' + text[start:end] + '\n' + (tests / 'neighbor_cpp/test_neighbor_list.cpp').read_text()
with tempfile.TemporaryDirectory(prefix='geos-neighbor-') as folder:
    folder = Path(folder)
    cpp = folder / 'test.cpp'
    cpp.write_text(program)
    for mode in ['host', 'device-adapter']:
        exe = folder / mode
        flags = ['-DMPM_DEBUG_GPU']
        if mode == 'device-adapter':
            flags.append('-DGEOS_USE_DEVICE')
        subprocess.run(['c++', '-std=c++17', '-O1', '-g', '-Wall', '-Wextra', '-Werror',
                        '-Wshadow', '-Wfloat-equal', '-pedantic', '-fsanitize=address,undefined',
                        '-fno-omit-frame-pointer', *flags, str(cpp), '-o', str(exe)], check=True)
        print('Testing ' + mode, flush=True)
        subprocess.run([str(exe)], check=True)
