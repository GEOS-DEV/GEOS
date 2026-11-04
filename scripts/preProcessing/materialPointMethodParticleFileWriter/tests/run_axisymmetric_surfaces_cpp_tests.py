"""Compile production axisymmetric surface kernels and solver excerpts on CPU.

Strict warnings and ASan/UBSan; no full GEOS, MPI, GPU or restart execution.
"""
from pathlib import Path
import re
import subprocess
import tempfile
from axisymmetric_test_support import program_prefix

root,prefix=program_prefix()
base=Path(__file__).resolve().parent
source=(root/'src/coreComponents/physicsSolvers/solidMechanics/SolidMechanicsMPM.cpp').read_text()
def block(text,marker):
    start=text.index(marker)
    masked=re.sub(r'//[^\n]*|/\*[\s\S]*?\*/|"(?:\\.|[^"\\])*"',lambda m:' '*len(m[0]),text)
    brace=masked.index('{',start);end=brace+1;depth=1
    while depth:depth+=(masked[end]=='{')-(masked[end]=='}');end+=1
    return text[start:end]

# Keep the production mapper with a larger 2D test grid; no 3D mapping is used.
prefix=prefix.replace('using arrayView3d=T(*)[9][9];','using arrayView3d=T(*)[65][1];')
prefix=prefix.replace('template<class T> using arraySlice1d=T*;','template<class T,int=0> using arraySlice1d=T*;')
program=prefix+'\n#include "physicsSolvers/solidMechanics/MPMAxisymmetricSurfaces.hpp"\n#include <vector>\n#include <array>\n#include <memory>\n'
program+='namespace mpm {enum class BoundaryConditionOption {Outflow,Symmetry};}\n'
program+='struct SurfaceKernel:Kernel {\n'
program+=(base/'native_geometry_cpp/axisymmetric_surface_test_adapter.hpp').read_text()
for marker in ('axisymmetricSurface::Support SolidMechanicsMPM::axisymmetricSurfaceSupport(',
               'localIndex SolidMechanicsMPM::axisymmetricSurfaceColorNode('):
    program+=block(source,marker).replace('SolidMechanicsMPM::','')+'\n'
program+='void assemble(Particle const & particle,int inputA,int inputB,double inputGamma) {\n'
program+=r"""
  int pp=0,pGroup=particle.group,numContactFlags=1,numContactGroups=2,numVelocityFields=2,numDims=2;
  double pVolume=particle.volume,gradientTolerance=1e-11/m_hEl[0],forceSign=1;
  double supportVolume=m_hEl[0]*m_hEl[1],invSupportVolume=1/supportVolume,volumeTolerance=1e-12*supportVolume;
  int effectiveMappedNodes[1][16]={},effectiveMappedFields[1][16]={};
  double weights[16]={},effectiveShapeFunctionGradientValues[1][16][3]={};
  int numberOfEffectiveMappedNodesPerParticle=particle.type==ParticleType::CPTI?12:
    (particle.type==ParticleType::SinglePoint?4:16);
  bool valid=computeParticleMapping<2>(particle.type,particle.x,particle.geometry,m_ijkMap,m_xLocalMin,m_hEl,m_nEl,
    effectiveMappedNodes[0],weights,effectiveShapeFunctionGradientValues[0]);
  assert(valid);
  for(int a=0;a<numberOfEffectiveMappedNodesPerParticle;++a)effectiveMappedFields[0][a]=particle.group;
  auto clamp01=[](double v){return std::clamp(v,0.,1.);};
"""
program+=block(source,'auto assembleSurfaceStressPair =')+';\nassembleSurfaceStressPair(inputA,inputB,inputGamma); }\n};\n'
program+=(base/'native_geometry_cpp/test_axisymmetric_surfaces.cpp').read_text()
program+='\nint main(){std::cout<<std::unitbuf;testMeasures();testRegression();testSurfaceStress();testGatherAndSphere();}\n'
with tempfile.TemporaryDirectory(prefix='geos-axis-surfaces-') as directory:
    tmp=Path(directory);(tmp/'common').mkdir()
    (tmp/'common/DataTypes.hpp').write_text('#pragma once\n#define GEOS_HOST_DEVICE\n#define GEOS_HOST\n#define GEOS_FORCE_INLINE inline\nnamespace geos {using real64=double;using integer=int;using localIndex=int;}\n')
    cpp=tmp/'test.cpp';cpp.write_text(program);exe=tmp/'test'
    subprocess.run(['c++','-std=c++17','-O1','-g','-Wall','-Wextra','-Werror','-Wshadow','-Wfloat-equal','-pedantic',
                    '-fsanitize=address,undefined','-fno-omit-frame-pointer',
                    '-I'+str(tmp),'-I'+str(root/'src/coreComponents'),str(cpp),'-o',str(exe)],check=True)
    subprocess.run([str(exe)],check=True)
