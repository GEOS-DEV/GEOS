"""Compile the production DIM kernels with scalar framework type aliases.

This is an arithmetic/geometry test, not a full GEOS, MPI or GPU build.
"""
from pathlib import Path
import subprocess
import tempfile

tests=Path(__file__).resolve().parent
root=next(p for p in tests.parents if (p/'src/coreComponents/mesh/ParticleGeometry.hpp').exists())
with tempfile.TemporaryDirectory(prefix='geos-native-2d-') as temp:
    temp=Path(temp)
    (temp/'common').mkdir()
    (temp/'common/DataTypes.hpp').write_text('#pragma once\n#define GEOS_HOST_DEVICE\nnamespace geos { using real64=double; using integer=int; using localIndex=int; }\n')
    executable=temp/'test'
    subprocess.run(['c++','-std=c++17','-O1','-g','-Wall','-Wextra','-Werror','-Wshadow','-Wfloat-equal','-pedantic',
                    '-fsanitize=address,undefined','-fno-omit-frame-pointer','-I'+str(temp),'-I'+str(root/'src/coreComponents'),
                    str(tests/'native_geometry_cpp/test_true_2d.cpp'),'-o',str(executable)],check=True)
    subprocess.run([str(executable)],check=True)
