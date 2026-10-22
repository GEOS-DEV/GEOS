"""Compile production CZ kernel/area/relative-F/reset/gather excerpts on CPU.

Adapters replace repository arrays and execution policies, not the numerical
algorithms under test. This does not compile or run the complete GEOS solver.
"""
from pathlib import Path
import re
import subprocess
import tempfile

base = Path(__file__).resolve().parent
root = next(p for p in base.parents if (p/'src/coreComponents').is_dir())
source = (root/'src/coreComponents/physicsSolvers/solidMechanics/SolidMechanicsMPM.cpp').read_text()
kernels = (root/'src/coreComponents/physicsSolvers/solidMechanics/kernels/ExplicitMPM.hpp').read_text()
uncoupled = (root/'src/coreComponents/constitutive/cohesiveZone/UncoupledCohesiveZone.hpp').read_text()

def block(text, start):
    # Braces in comments and strings are not C++ scope delimiters.
    masked = re.sub(r'//[^\n]*|/\*[\s\S]*?\*/|"(?:\\.|[^"\\])*"', lambda m: ' '*len(m[0]), text)
    begin = masked.index('{', start)
    end, depth = begin+1, 1
    while depth:
        depth += (masked[end] == '{') - (masked[end] == '}')
        end += 1
    return text[begin:end]

prefix = (base/'native_geometry_cpp/cz_test_support.hpp').read_text()
program = prefix + '\nstruct CohesiveZoneStateUpdateKernel ' + block(kernels, kernels.index('struct CohesiveZoneStateUpdateKernel'))+';\n'
program += '''struct Law {
  double m_normalForceConstant=200, m_shearForceConstant=100;
  double m_maxNormalDisplacement=.1, m_maxTangentialDisplacement=.1;
  double *m_damage;
  void saveConvergedState(int,int) const {}
  void jumpDisplacementUpdate(int k, double normalDisplacement, double tangentialDisplacement,
                              double &normalStress, double &shearStress) const
'''+block(uncoupled,uncoupled.index('void jumpDisplacementUpdate'))+'\n};\n'
a=source.index('void SolidMechanicsMPM::enforceCohesiveLaw(')
b=source.index('real64 deformationGradientCofactor[3][3] = {};',a)
e=source.index('real64 centerDisplacement[3]',b)
program += '''void relativeArea(double const (&particleDeformationGradient)[1][3][3],
    double const (&particleCohesiveReferenceDeformationGradient)[1][3][3], double (&out)[3][3]) {
    int const p=0;
'''+source[b:e]+'''LvArray::tensorOps::copy<3,3>(out,deformationGradientCofactor); }\n'''
start=source.index('        // Preserve F F_cz^{-1} across the reset:')
end=source.index('\n      }',start)
program+='''void rebase(double (&particleCohesiveReferenceDeformationGradient)[1][3][3],
 double const (&oldF)[3][3], double const (&newDeformationGradient)[3][3]) {
 int const p=0; double inverseOldDeformationGradient[3][3];
 LvArray::tensorOps::copy<3,3>(inverseOldDeformationGradient,oldF);
 LvArray::tensorOps::invert<3>(inverseOldDeformationGradient);
'''+source[start:end]+'\n}\n'
start=source.index('    // Map cohesive law back to particle',a)
start=source.index('          if( particleCZTag[p] == czTag',start)
program+='''void partialGather(double const (&damage)[2], int &flag, double &force) {
 int const p=0, czTag=1, numDims=2, mappedNodeCapacity=2;
 double const smallMass=1e-20;
 globalIndex const invalidCohesiveNodeIndex=99;
 int particleCZTag[1]={1}, particleSurfaceFlag[1]={3}, particleCohesiveZoneFlag[1]={1};
 globalIndex particleReferenceGlobalNodeIndex[1][2]={{0,1}};
 int particleCohesiveFieldMapping[1][2]={{0,0}};
 double particleReferenceShapeFunctionValues[1][2]={{.4,.6}}, particleMass[1]={1}, particleDamage[1]={0};
 double tempGridMassGlobal[2][2]={{.4,.4},{.6,.6}};
 double tempGridCohesiveTraction[2][2][3]={{{2,0,0},{-2,0,0}},{{3,0,0},{-3,0,0}}};
 double particleCohesiveForce[1][3]={};
 auto const &czDamage=damage;
'''+source[start:source.index('{',start)]+block(source,start)+'''
 flag=particleCohesiveZoneFlag[0];force=particleCohesiveForce[0][0];
}\n'''
start=source.index('              auto const particleDomain = particleGeometry::affine::Domain::fromGenerators(',source.index('void SolidMechanicsMPM::initializeCohesiveReferenceConfiguration('))
end=source.index('\n            }',start)
program += """void anchor(int shape, double const (&g)[2][2],double const (&n)[3],double (&out)[3]) {
  int const p=0; int geometryShape[1]={shape};double materialGeometry[1][2][2],particleSurfaceNormal[1][3];
  double particleCohesiveReferenceSurfacePosition[1][3]={};
  LvArray::tensorOps::copy<2,2>(materialGeometry[0],g);LvArray::tensorOps::copy<3>(particleSurfaceNormal[0],n);
"""+source[start:end]+"""
  LvArray::tensorOps::copy<3>(out,particleCohesiveReferenceSurfacePosition[0]);
}\n"""
program += (base/'native_geometry_cpp/test_axisymmetric_cz.cpp').read_text()
with tempfile.TemporaryDirectory(prefix='geos-axisym-cz-') as directory:
    tmp=Path(directory); (tmp/'common').mkdir()
    (tmp/'common/DataTypes.hpp').write_text('#pragma once\n#define GEOS_HOST_DEVICE\n#define GEOS_FORCE_INLINE inline\nnamespace geos { using real64=double; using integer=int; using localIndex=int; }\n')
    cpp=tmp/'test.cpp'; cpp.write_text(program); exe=tmp/'test'
    subprocess.run(['c++','-std=c++17','-O1','-g','-Wall','-Wextra','-Werror','-Wshadow','-Wfloat-equal','-pedantic',
                    '-fsanitize=address,undefined','-fno-omit-frame-pointer',
                    '-I'+str(tmp),'-I'+str(root/'src/coreComponents'),str(cpp),'-o',str(exe)],check=True)
    subprocess.run([str(exe)],check=True)
