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
    parser.add_argument('--writer-results', type=Path,
                        help='Optional run_writer_regression directory to validate with the production C++ codec')
    args=parser.parse_args()
    tests=Path(__file__).resolve().parent
    root=next(p for p in tests.parents if (p/'src/coreComponents/mesh/ParticleGeometry.hpp').is_file())
    with tempfile.TemporaryDirectory(prefix='geos-native-geometry-') as directory:
        tmp=Path(directory)
        (tmp/'common').mkdir()
        (tmp/'common/DataTypes.hpp').write_text(
            '#pragma once\n#define GEOS_HOST_DEVICE\n'
            'namespace geos { using real64=double; using integer=int; using localIndex=int; }\n')
        for name in ('test_unified_geometry',):
            executable=tmp/name
            subprocess.run([args.cxx,'-std=c++17','-Wall','-Wextra','-Werror','-pedantic',
                            '-I',str(tmp),'-I',str(root/'src/coreComponents'),
                            str(tests/'native_geometry_cpp'/(name+'.cpp')),'-o',str(executable)],check=True)
            subprocess.run([str(executable)],check=True)
        if args.writer_results:
            import numpy as np
            count=0
            for path in sorted(args.writer_results.glob('*_r1/mpmParticleFile_case')):
                names=path.read_text().splitlines()[0].split()
                data=np.loadtxt(path,skiprows=1,ndmin=2)
                columns=['GeometryShape','ParticleVolume']+['Geometry'+ij for ij in
                         ('XX','XY','XZ','YX','YY','YZ','ZX','ZY','ZZ')]
                rows=data[:,[names.index(c) for c in columns]]
                import io
                text=io.StringIO(); np.savetxt(text,rows,fmt='%.17g')
                run=subprocess.run([str(executable),'--decode'],input=text.getvalue(),
                                   text=True,capture_output=True,check=True)
                assert len(run.stdout.splitlines())==len(rows)
                count+=len(rows)
            assert count>0
            print('Production C++ input validator accepted %d generated PFW records.'%count)

if __name__ == '__main__':
    main()
