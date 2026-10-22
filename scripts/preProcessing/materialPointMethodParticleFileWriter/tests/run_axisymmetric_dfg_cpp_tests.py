"""Compile production DFG reconstruction/partition/separability and contact code.

Reuses CPU adapters from the ordinary contact test. Checks actual source
excerpts, not a full GEOS build, MPI run or GPU execution.
"""
import argparse
import os
from pathlib import Path
import re
import subprocess
import tempfile

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--particle-files', type=Path, nargs='*', default=[],
                    help='optional generated axisymmetricDFG files; check central-band mapping')
args = parser.parse_args()
tests = Path(__file__).resolve().parent
support = (tests/'run_axisymmetric_contact_cpp_tests.py').read_text()
scope = {'__file__': str(tests/'run_axisymmetric_contact_cpp_tests.py')}
# The setup before this marker only builds the production-code test program.
exec(compile(support[:support.index('\nwith tempfile.TemporaryDirectory')],
             str(tests/'run_axisymmetric_contact_cpp_tests.py'), 'exec'), scope)
root, source = scope['root'], scope['source']


def extract(text, signature):
    text = text[text.index(signature):]
    masked = re.sub(r'''//[^\n]*|/\*[\s\S]*?\*/|"(?:\\.|[^"\\])*"|'(?:\\.|[^'\\])*' ''',
                    lambda match: ' '*len(match[0]), text, flags=re.X)
    end, depth = masked.index('{')+1, 1
    while depth:
        depth += (masked[end] == '{') - (masked[end] == '}')
        end += 1
    return text[:end].replace('SolidMechanicsMPM::', '')

adapters = scope['adapters'].replace('int main() {mappedInterfaces();contactLaws();}', '')
program = scope['prefix'] + scope['definitions'] + adapters
program += '\n#include "physicsSolvers/solidMechanics/MPMDamageGradient.hpp"\n#include <vector>\n#include <fstream>\n#include <sstream>\n'
enums = scope['enums']
program += 'namespace mpm {\n' + extract(enums, 'enum struct SurfaceFlag') + ';\n'
program += extract(enums, 'inline integer toInteger( SurfaceFlag') + '\n}\n'
program += '#define GEOS_ERROR_IF(condition, message) do { if(condition) throw std::runtime_error(message); } while(false)\n'
start = source.index('constexpr real64 explicitSurfaceNormalTolerance')
end = source.index('\n}', source.index('bool hasExplicitPartitionSurfaceNormal', start))+2
program += source[start:end] + '\n'
for name, result in [('markSurfaceAsDamage','bool'), ('partitionField','localIndex'),
                     ('computeSurfaceQualityFromMappingNormalTensor','real64'),
                     ('evaluateSeparabilityCriterion','bool'), ('kernelDevice','real64')]:
    program += extract(source, result+' SolidMechanicsMPM::'+name+'(')+'\n'

# Compile the actual axisymmetric neighbor loop, including
# the surface flag convention. Adapters provide the nested repository views.
start = source.index('void SolidMechanicsMPM::computeDamageFieldGradient(')
start = source.index('      if( axisymmetric )', start)
block = extract(source[start:], '      if( axisymmetric )')[len('      if( axisymmetric )'):]
program += r'''
struct Sample {double r,y,area,d; int group=0,flag=0;};
Vec reconstruct(std::vector<Sample> const & samples,double queryRadius,double axial,double neighborRadius) {
  constexpr int count=1024; assert(samples.size()<=count);
  int const p=0,pp=0;
  int numNeighborsAll[1]={int(samples.size())};
  int neighborRegions[1][count]={},neighborSubRegions[1][count]={},neighborIndices[1][count]={};
  int particleSurfaceFlagView[1][1][count]={},particleCohesiveZoneFlagView[1][1][count]={};
  double particlePosition[1][3]={{queryRadius,axial,0}},particlePositionView[1][1][count][3]={};
  double particleDamageView[1][1][count]={},particleVolumeView[1][1][count]={},grad[3]={};
  for(int n=0;n<int(samples.size());++n) {
    auto const &q=samples[n];neighborIndices[0][n]=n;
    particleSurfaceFlagView[0][0][n]=q.flag;
    particlePositionView[0][0][n][0]=q.r;particlePositionView[0][0][n][1]=q.y;
    particleDamageView[0][0][n]=q.d;particleVolumeView[0][0][n]=axisymmetric::twoPi*q.r*q.area;
  }
''' + block + r'''
  return {grad[0],grad[1],grad[2]};
}
'''
program += (tests/'native_geometry_cpp/test_axisymmetric_dfg.cpp').read_text()

with tempfile.TemporaryDirectory(prefix='geos-axis-dfg-') as directory:
    tmp = Path(directory)
    (tmp/'common').mkdir()
    (tmp/'common/DataTypes.hpp').write_text(
        '#pragma once\n#define GEOS_HOST_DEVICE\n#define GEOS_FORCE_INLINE inline\n'
        'namespace geos { using real64=double; using integer=int; using localIndex=int; }\n')
    cpp = tmp/'test.cpp'
    cpp.write_text(program)
    exe = tmp/'test'
    subprocess.run(['c++', '-std=c++17', '-O1', '-g', '-Wall', '-Wextra', '-Werror',
                    '-Wshadow', '-Wfloat-equal', '-pedantic', '-fsanitize=address,undefined',
                    '-fno-omit-frame-pointer', '-I'+str(tmp), '-I'+str(root/'src/coreComponents'),
                    str(cpp), '-o', str(exe)], check=True)
    env = dict(os.environ)
    env.setdefault('ASAN_OPTIONS', 'detect_leaks=0')
    subprocess.run([str(exe), *map(str, args.particle_files)], check=True, env=env)
