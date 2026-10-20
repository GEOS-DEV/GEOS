"""Compile the production geometry helper in a standalone CPU test.

Uses a temporary minimal DataTypes header containing only typedefs/macros. This
checks the actual geometry math, not GEOS registration, MPI, CUDA/HIP or linking.
Run from anywhere inside a full repository or the supplied files overlay.
"""
import argparse
from pathlib import Path
import subprocess
import tempfile


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cxx',default='c++')
    args=parser.parse_args()
    tests=Path(__file__).resolve().parent
    root=next(p for p in tests.parents if (p/'src/coreComponents/mesh/ParticleReferenceGeometry.hpp').is_file())
    with tempfile.TemporaryDirectory(prefix='geos-native-geometry-') as directory:
        tmp=Path(directory)
        (tmp/'common').mkdir()
        (tmp/'common/DataTypes.hpp').write_text(
            '#pragma once\n#define GEOS_HOST_DEVICE\n'
            'namespace geos { using real64=double; using integer=int; using localIndex=int; }\n')
        executable=tmp/'test_geometry'
        subprocess.run([args.cxx,'-std=c++17','-Wall','-Wextra','-Werror','-pedantic',
                        '-I',str(tmp),'-I',str(root/'src/coreComponents'),
                        str(tests/'native_geometry_cpp/test_geometry.cpp'),'-o',str(executable)],check=True)
        subprocess.run([str(executable)],check=True)


if __name__ == '__main__':
    main()
