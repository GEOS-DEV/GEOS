"""Test actual dimension-specialized CPDI projection and shared mapping bodies."""
from pathlib import Path
import subprocess
import tempfile
from production_mapping_test import program_prefix
root, prefix=program_prefix()
tests=Path(__file__).resolve().parent
program=prefix+(tests/'native_geometry_cpp/test_cpdi_geometry.cpp').read_text()
with tempfile.TemporaryDirectory(prefix='geos-cpdi-geometry-') as folder:
    tmp=Path(folder); (tmp/'common').mkdir()
    (tmp/'common/DataTypes.hpp').write_text('#pragma once\n#define GEOS_HOST_DEVICE\n#define GEOS_FORCE_INLINE inline\nnamespace geos { using real64=double; using integer=int; using localIndex=int; }\n')
    cpp=tmp/'test.cpp'; cpp.write_text(program); exe=tmp/'test'
    subprocess.run(['c++','-std=c++17','-O1','-g','-Wall','-Wextra','-Werror','-Wshadow','-Wfloat-equal','-pedantic',
                    '-fsanitize=address,undefined','-fno-omit-frame-pointer',
                    '-I'+str(tmp),'-I'+str(root/'src/coreComponents'),str(cpp),'-o',str(exe)],check=True)
    subprocess.run([str(exe)],check=True)
