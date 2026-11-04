"""Compile production wall, table, grid-map and controller code with CPU adapters.

No GEOS runtime is available here. Array/reduction/MPI and constitutive lookup
adapters replace dependencies; table interpolation, controller algebra (including
adaptive seek), wall target assembly, validation and resize bodies are production.
"""
from pathlib import Path
import os
import re
import subprocess
import tempfile

base=Path(__file__).resolve().parent
root=next(p for p in base.parents if (p/'src/coreComponents').exists())
source=(root/'src/coreComponents/physicsSolvers/solidMechanics/SolidMechanicsMPM.cpp').read_text()
header=(root/'src/coreComponents/physicsSolvers/solidMechanics/SolidMechanicsMPM.hpp').read_text()

def block(signature, text=source):
    start=text.index(signature)
    masked=re.sub(r'//[^\n]*|/\*[\s\S]*?\*/|"(?:\\.|[^"\\])*"',lambda m:' '*len(m[0]),text[start:])
    end=masked.index('{')+1;depth=1
    while depth:
        depth+=(masked[end]=='{')-(masked[end]=='}');end+=1
    return text[start:start+end]

def method(name):
    return block('void SolidMechanicsMPM::'+name+'(').replace('SolidMechanicsMPM::','')

program=(base/'native_geometry_cpp/driving_test_support.hpp').read_text()
for sig in ('real64 vectorInfinityNorm3(', 'void capVectorInfinityNorm3(', 'bool solveSmallLinearSystem3(', 'void dampedLeastSquaresSolve3('):
    program+=block(sig)+'\n'
a=source.index('enum : int\n{\n  mpmBCTargetNone')
b=source.index('real64 vectorInfinityNorm3(',a)
program+=source[a:b]
program+='struct Driving {\n'+r'''
 bool m_axisymmetric=true;int m_numDims=2,m_prescribedBoundaryFTable=1,m_prescribedFTable=0,m_prescribedBcTable=0;
 int m_boxAverageResizeWithDomain=0;
 double m_hEl[3]={.1,.1,1},m_xGlobalMin[3]={0,0,0},m_xGlobalMax[3]={1,1,0},m_domainExtent[3]={1,1,0};
 double m_xLocalMin[3]={-.1,-.1,0},m_xLocalMax[3]={1.1,1.1,0},m_xLocalMinNoGhost[3]={0,0,0},m_xLocalMaxNoGhost[3]={1,1,0};
 double m_partitionExtent[3]={1.2,1.2,0},m_boxAverageMin[3]={0,0,0},m_boxAverageMax[3]={1,1,0};
 double m_domainF[3]={1,1,1},m_domainL[3]={},m_domainIncrementalStretch[3]={1,1,1};
 Array<double> m_domainStress=Array<double>({0,0,0});
 int m_boundaryConditionTypes[6]={1,2,1,2,0,0},m_enablePrescribedBoundaryTransverseVelocities[6]={};
 double m_prescribedBoundaryTransverseVelocities[6][2]={};
 mpm::UpdateMethodOption m_updateMethod=mpm::UpdateMethodOption::FLIP;
 mpm::InterpolationOption m_fTableInterpType=mpm::InterpolationOption::Linear,m_stressTableInterpType=mpm::InterpolationOption::Linear;
 Table m_fTable{{0,1,1,1},{1,.9,.8,1}},m_stressTable{{0,0,0,0},{1,-100,0,0}},m_bcTable;
 double m_cflFactor=.2,m_planeStrainThickness=1;
'''
# Preserve production scalar controller defaults and member layout.
for line in header.splitlines():
    match=re.match(r'  (int|real64|array1d< int >|array1d< real64 >|array2d< real64 >) (m_stressControl\w*);',line)
    if not match:continue
    kind,name=match.groups()
    default=re.search(r'  '+name+r'\( ([^()]+) \)',source)
    if kind in ('int','real64'):
        program+=kind+' '+name+'='+default[1]+';\n'
    elif kind=='array2d< real64 >':program+='double '+name+'[3][3]={};\n'
    else:program+=('int ' if 'int' in kind else 'double ')+name+'[3]={};\n'
program+=block('  real64 physicalDomainVolume() const',header)+'\n'
program+='template<class T> double physicalVolume(T const &e) const{return e[0]*e[1]*(m_numDims==2?m_planeStrainThickness:e[2]);}\n'
for name in ('interpolateTable','interpolateFTable','interpolateStressTable','validateAxisymmetricBoundaryControls','boundaryConditionUpdate','computeBoxMetrics','adaptiveSeekStressControl','updatePrescribedKinematicsForExplicitStep','resizeGrid'):
    program+=method(name)+'\n'
# Only the constitutive lookup is replaced: adapters supply its bulk modulus.
s=method('stressControl')
a=s.index('  // arrayView1d< int const >');b=s.index('  // Could some numerical artifact')
s=s[:a]+'''  double targetStress[3];LvArray::tensorOps::copy<3>(targetStress,m_domainStress);
  Array<double> currentStress(3);double boxMaterialVolume=0;
  computeBoxMetrics(particleManager,currentStress,boxMaterialVolume);
  double const maximumBulkModulus=particleManager.bulk;
'''+s[b:]
program+=s+'\n};\n'

# Actual per-node wall target assembly and staging helpers, including legacy
# Cartesian branch, compiled under strict warnings. Outer device/reduce loops
# are replaced by one node and one material field.
s=block('void SolidMechanicsMPM::applyEssentialBCs(')
a=s.index('          // Contact faces prescribe');b=s.index('\n        } );',a)
program+='void Wall::apply() {\n'+s[a:b]+'\n}\n'
program+=(base/'native_geometry_cpp/test_axisymmetric_driving.cpp').read_text()
with tempfile.TemporaryDirectory(prefix='geos-axis-driving-') as directory:
    tmp=Path(directory);cpp=tmp/'test.cpp';cpp.write_text(program);exe=tmp/'test'
    subprocess.run(['c++','-std=c++17','-O1','-g','-Wall','-Wextra','-Werror','-Wshadow','-Wfloat-equal','-pedantic',
        '-fsanitize=address,undefined','-fno-omit-frame-pointer',str(cpp),'-o',str(exe)],check=True)
    env=dict(os.environ);env.setdefault('ASAN_OPTIONS','detect_leaks=0')
    subprocess.run([str(exe)],check=True,env=env)
