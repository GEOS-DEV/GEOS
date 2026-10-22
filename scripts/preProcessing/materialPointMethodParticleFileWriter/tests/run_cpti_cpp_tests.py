"""Compile the production CPTI kernels with host adapters, not a GEOS build.

Both array-view-style outputs and raw on-demand arrays instantiate the same
production template. Run --sanitize for address/undefined-behavior checks.
"""
import argparse
from pathlib import Path
import subprocess
import tempfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sanitize', action='store_true')
    args = parser.parse_args()
    tests = Path(__file__).resolve().parent
    from production_mapping_test import program_prefix
    root,prefix=program_prefix()
    test=(tests/'native_geometry_cpp/test_cpti.cpp').read_text()
    test=test.replace('kernel.computeCPTIShapeFunctions(', 'kernel.computeParticleMapping<3>(ParticleType::CPTI,')
    program=prefix+test
    with tempfile.TemporaryDirectory(prefix='geos-cpti-') as folder:
        tmp = Path(folder)
        (tmp/'common').mkdir()
        (tmp/'common/DataTypes.hpp').write_text('#pragma once\n#define GEOS_HOST_DEVICE\n#define GEOS_FORCE_INLINE inline\n'
            'namespace geos { using real64=double; using integer=int; using localIndex=int; }\n')
        cpp = tmp/'test.cpp'
        cpp.write_text(program)
        exe = tmp/'test'
        flags = ['-fsanitize=address,undefined', '-fno-omit-frame-pointer'] if args.sanitize else []
        subprocess.run(['c++', '-std=c++17', '-O1', '-g', '-Wall', '-Wextra', '-Werror', '-Wshadow',
                        '-Wfloat-equal', '-pedantic', *flags, '-I'+str(tmp),
                        '-I'+str(root/'src/coreComponents'), str(cpp), '-o', str(exe)], check=True)
        subprocess.run([str(exe)], check=True)


if __name__ == '__main__':
    main()
