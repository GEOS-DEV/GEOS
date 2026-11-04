"""Compile axisymmetric boundary regressions against production source excerpts.

Exercises the shared mapper, the real cleanup corner loop, surface-detection
neighbor loop, and the event allow-list/swap handler. Repository arrays and
scheduler inputs have CPU adapters; this is not a full GEOS/MPI/GPU build.
"""
from pathlib import Path
import os
import re
import subprocess
import tempfile
from axisymmetric_test_support import program_prefix

root, prefix = program_prefix()
tests = Path(__file__).resolve().parent
source = (root/'src/coreComponents/physicsSolvers/solidMechanics/SolidMechanicsMPM.cpp').read_text()

def block(text, signature):
    text = text[text.index(signature):]
    masked = re.sub(r'//[^\n]*|/\*[\s\S]*?\*/|"(?:\\.|[^"\\])*"',
                    lambda m: ' '*len(m[0]), text)
    end, depth = masked.index('{')+1, 1
    while depth:
        depth += (masked[end] == '{') - (masked[end] == '}')
        end += 1
    return text[:end]

# Compile the actual cleanup body, not a restatement of its projection policy.
cleanup = source[source.index('void SolidMechanicsMPM::flagOutOfRangeParticles('):]
start = cleanup.index('          real64 integration[3][3];')
end = cleanup.index('\n        } );', start)
cleanup = cleanup[start:end]
method = r'''
  template<int DIM>
  bool cleanup(ParticleType type, double const (&center)[3], double const (&geometry)[3][3]) const {
    constexpr int p=0;
    int particleDeleteFlag[1]={},periodic[3]={};
    double particlePosition[1][3];double particleGeometryCurrent[1][3][3];
    for(int i=0;i<3;++i) {particlePosition[0][i]=center[i];
      for(int j=0;j<3;++j) particleGeometryCurrent[0][i][j]=geometry[i][j];}
    ParticleType const particleType=type;
    bool const isCPTI=type==ParticleType::CPTI;
    int numReportedOutOfRangeCPDIParticles=0;
    double globalMin[3],globalMax[3],tolerance[3];
    for(int i=0;i<3;++i) {globalMin[i]=m_xGlobalMin[i]-m_hEl[i];globalMax[i]=m_xGlobalMax[i]+m_hEl[i];
      tolerance[i]=64*DBL_EPSILON*std::max({m_hEl[i],std::abs(globalMin[i]),std::abs(globalMax[i])});}
''' + cleanup + r'''
    (void)numReportedOutOfRangeCPDIParticles;
    return particleDeleteFlag[0]!=0;
  }
'''
prefix = prefix[:-3] + method + '};\n'  # append to the extracted Kernel class

# Surface loop with repository-view adapters. Threshold/flag assignment also
# come from production, so a helper-only test cannot mask a wiring regression.
surface = source[source.index('void SolidMechanicsMPM::computeSurfaceFlags('):]
surface = block(surface, '      if( axisymmetric )')
program = prefix + r'''
#include "physicsSolvers/solidMechanics/MPMSurfaceDetection.hpp"
#include <array>
#include <vector>
#include <map>
#include <limits>
#include <stdexcept>
#define GEOS_ERROR_IF(condition,message) do {if(condition) throw std::runtime_error("input rejected");} while(false)
#define GEOS_LOG_RANK_0(message) do {} while(false)
#define GEOS_LOG_RANK_0_IF(condition,message) do {} while(false)
namespace mpm { enum struct SurfaceFlag { FullyDamaged=1 }; inline int toInteger(SurfaceFlag v){return int(v);} }
struct SurfaceSample {double r,y,volume;};
int detect(std::vector<SurfaceSample> const & samples,double queryRadius,double axial,double support,int flag=0) {
  int const p=0,numNeighbors=int(samples.size());
  std::vector<int> regionIndices(samples.size()),subRegionIndices(samples.size()),particleIndices(samples.size());
  std::vector<std::array<double,3>> positions(samples.size());std::vector<double> volumes(samples.size());
  for(int q=0;q<numNeighbors;++q) {particleIndices[q]=q;positions[q]={samples[q].r,samples[q].y,0};volumes[q]=samples[q].volume;}
  auto const *particlePositionView0=positions.data();auto const **particlePositionView1=&particlePositionView0;
  auto const ***particlePositionView=&particlePositionView1;
  auto const *particleVolumeView0=volumes.data();auto const **particleVolumeView1=&particleVolumeView0;
  auto const ***particleVolumeView=&particleVolumeView1;
  double const particlePosition[1][3]={{queryRadius,axial,0}},neighborRadius=support;
  int particleSurfaceFlag[1]={flag};bool const axisymmetric=true;
  auto evaluate=[&]() {
''' + surface + r'''
  };
  evaluate();return particleSurfaceFlag[0];
}
'''

# Actual guard and swap handler, with a minimal event object. Scheduler gating
# itself is unchanged; include its actual completed-event gate in this adapter.
event = source[source.index('void SolidMechanicsMPM::triggerEvents('):]
start = event.index('    GEOS_ERROR_IF(m_axisymmetric &&')
end = event.index('    // 1. Do not run completed events.', start)
guard = event[start:end]
completed = block(event, '    if( event.isComplete() )')
swap = block(event, '    if( event.getCatalogName() == "FrictionCoefficientSwap" )')
program += r'''
struct FrictionCoefficientSwapMPMEvent {
  std::string catalog="FrictionCoefficientSwap";double coefficient=.25;
  std::vector<double> table={.2,.3,.4};bool complete=false;
  std::string getCatalogName() const{return catalog;}
  double getFrictionCoefficient() const{return coefficient;}
  std::vector<double> getFrictionCoefficientTable() const{return table;}
  bool isComplete() const{return complete;}void setIsComplete(int v){complete=v!=0;}
};
template<class T> T & dynamicCast(FrictionCoefficientSwapMPMEvent &event){return static_cast<T&>(event);}
struct FrictionHarness {
  bool m_axisymmetric=true;double m_frictionCoefficient=0;std::vector<double> m_frictionCoefficientTable;
  int initialized=0;void initializeFrictionCoefficients(){++initialized;}
  void trigger(FrictionCoefficientSwapMPMEvent &event) {
''' + guard + completed + swap + r'''
  }
};
'''
program += (tests/'native_geometry_cpp/test_axisymmetric_boundary.cpp').read_text()
with tempfile.TemporaryDirectory(prefix='geos-axis-boundary-') as directory:
    tmp=Path(directory);(tmp/'common').mkdir()
    (tmp/'common/DataTypes.hpp').write_text('#pragma once\n#define GEOS_HOST_DEVICE\n#define GEOS_FORCE_INLINE inline\n'
       'namespace geos {using real64=double;using integer=int;using localIndex=int;}\n')
    cpp=tmp/'test.cpp';cpp.write_text(program);exe=tmp/'test'
    subprocess.run(['c++','-std=c++17','-O1','-g','-Wall','-Wextra','-Werror','-Wshadow','-Wfloat-equal','-pedantic',
       '-fsanitize=address,undefined','-fno-omit-frame-pointer','-I'+str(tmp),'-I'+str(root/'src/coreComponents'),str(cpp),'-o',str(exe)],check=True)
    env=dict(os.environ);env.setdefault('ASAN_OPTIONS','detect_leaks=0')
    subprocess.run([str(exe)],check=True,env=env)
