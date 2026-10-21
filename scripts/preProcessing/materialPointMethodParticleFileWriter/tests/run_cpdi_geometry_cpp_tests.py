"""Compile production CPDI projection/mapping bodies with small CPU adapters.

This checks the actual changed kernels, without pretending to compile GEOS.
The test supplies flat array views and basic tensor operations (including a
Jacobi eigensolve); GEOS/LvArray integration, GPU and MPI still require a build.
"""
from pathlib import Path
import re
import subprocess
import tempfile

tests = Path(__file__).resolve().parent
root = next(p for p in tests.parents if (p/'src/coreComponents/mesh/ParticleGeometry.hpp').exists())
source = (root/'src/coreComponents/physicsSolvers/solidMechanics/SolidMechanicsMPM.cpp').read_text()


def function(signature):
    start = source.index(signature)
    a = source.index('{', start); b = a+1; depth = 1
    while depth:
        depth += (source[b] == '{') - (source[b] == '}'); b += 1
    return source[start:b]


start = source.index('struct VPHenckyDomainProjectionResult')
end = source.index('GEOS_HOST_DEVICE\nGEOS_FORCE_INLINE\nint countSubdivisionDirections', start)
free = source[start:end]
methods = [function(signature).replace('SolidMechanicsMPM::', '') for signature in (
    'int SolidMechanicsMPM::projectIntegrationGeometry(',
    'void SolidMechanicsMPM::invalidateShapeFunctionRow(',
    'int SolidMechanicsMPM::computeLocalCellIndexWithBoundaryTolerance(',
    'bool SolidMechanicsMPM::computeCPDIShapeFunctions(',
    'void SolidMechanicsMPM::computeCPDIParticleShapeFunctions(')]
program = (tests/'native_geometry_cpp/cpdi_test_support.hpp').read_text() + free + '''
struct Kernel {
  int m_cpdiDomainScaling=0, m_numDims=3, m_planeStrain=0;
  double m_neighborRadius=1.1, m_hEl[3]={1,1,1};
  mpm::CPDIDomainScalingTypeOption m_cpdiDomainScalingType=mpm::CPDIDomainScalingTypeOption::VPHencky;
''' + '\n'.join(methods) + '\n};\n' + (tests/'native_geometry_cpp/test_cpdi_geometry.cpp').read_text()
with tempfile.TemporaryDirectory(prefix='geos-cpdi-geometry-') as tmp:
    tmp = Path(tmp)
    (tmp/'common').mkdir()
    (tmp/'common/DataTypes.hpp').write_text('#pragma once\n#define GEOS_HOST_DEVICE\n#define GEOS_FORCE_INLINE inline\n'
        'namespace geos { using real64=double; using integer=int; using localIndex=int; }\n')
    cpp = tmp/'test.cpp'; cpp.write_text(program)
    exe = tmp/'test'
    subprocess.run(['c++','-std=c++17','-O2','-Wall','-Wextra','-Werror','-pedantic',
                    '-I'+str(tmp),'-I'+str(root/'src/coreComponents'),str(cpp),'-o',str(exe)],check=True)
    subprocess.run([str(exe)],check=True)
