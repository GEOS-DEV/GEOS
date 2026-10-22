"""Compile production subdivision code with standalone host adapters.

No GEOS link, MPI, GPU runtime, or constitutive integration is exercised.
"""
import argparse
from pathlib import Path
import subprocess
import tempfile
import re



def solver_program(root, tests):
    source=(root/'src/coreComponents/physicsSolvers/solidMechanics/SolidMechanicsMPM.cpp').read_text()
    def function(signature):
        start=source.index(signature);a=source.index('{',start);b=a+1;depth=1
        while depth:
            depth+=(source[b]=='{')-(source[b]=='}');b+=1
        return source[start:b]
    splitter=function('bool SolidMechanicsMPM::subdivideParticles(')
    scalar_int={'particleSubdivideFlag','particleCopyFlag','particleDeleteFlag','particleCohesiveZoneFlag','particleDomainScaledFlag'}
    rank2={'particleReferencePosition','particleReferenceSurfaceTraction','particleCohesiveForce','particleReferenceSurfacePosition','particleEstimatedSurfacePosition'}
    rank3={'particleGeometryReference','particleDeformationGradient','particleVelocityGradient'}
    names=sorted(set(re.findall(r'getField< fields::mpm::(\w+) >',splitter)))
    tags=[];create=[]
    for name in names:
        rank=3 if name in rank3 else 2 if name in rank2 else 1
        dtype='int' if name in scalar_int else 'double'
        key='meltFlag' if name=='particleMeltFlag' else name
        atype=f'Array<{dtype},{rank}>'
        tags.append('struct '+name+' { using Type='+atype+'; static string key(){return "'+key+'";} };')
        create.append('r.add<'+atype+'>("'+key+'"'+',3'*(rank-1)+');')
    support=(tests/'subdivision_cpp/solver_test_support.hpp').read_text().replace('// FIELD_TAGS','namespace fields { namespace mpm {\n'+'\n'.join(tags)+'\n}}')
    helpers=[]
    for signature,prefix in [
        ('mpm::OversizedParticleTreatmentOption effectiveOversizedParticleTreatment(',''),
        ('bool exceedsOversizedParticleResetRankThreshold(',''),
        ('real64 particleDomainMeasure(','template<class R_VECTORS>\n'),
        ('real64 particleDomainMaximumFeasibleMeasure(',''),
        ('int countSubdivisionDirections(',''),
        ('int longestParticleDomainDirectionsMask(','template<class R_VECTORS>\n'),
        ('real64 particleDomainMaximumDiagonalAfterSubdivisionMask(','template<class R_VECTORS>\n')]:
        helpers.append(prefix+function(signature))
    solver=r"""
struct Solver {
 int m_numDims=3,m_planeStrain=0,m_numberOfSubRegions=1;
 int m_applyVelocityGradientWhenSplittingParticles=1,m_subdivideParticlesAutomatic=0,m_subdivideParticles=1;
 int m_cpdiDomainScaling=0,m_resetDefGradForOversizedParticles=0,m_oversizedParticleResetRankParticleCountThreshold=-1;
 int m_resetDefGradForFullyDamagedParticles=0,m_resetDefGradForMeltedParticles=0,m_resetDefGradForScaledSurfaceParticles=0;
 double m_neighborRadius=1.5,m_hEl[3]={1,1,1},m_defGradResetMaxParticleDomainToGridCellRatio=-1;
 mpm::CPDIDomainScalingTypeOption m_cpdiDomainScalingType=mpm::CPDIDomainScalingTypeOption::Homel;
 mpm::DomainResetTypeOption m_domainResetType=mpm::DomainResetTypeOption::Isotropic;
 mpm::OversizedParticleTreatmentOption m_oversizedParticleTreatment=mpm::OversizedParticleTreatmentOption::None;
 std::vector<Array<int,2>> m_mappedNodes,m_mappedFields,m_effectiveMappedNodes,m_effectiveMappedFields;
 std::vector<Array<double,2>> m_shapeFunctionValues,m_effectiveShapeFunctionValues;
 std::vector<Array<double,3>> m_shapeFunctionGradientValues,m_effectiveShapeFunctionGradientValues;
 std::vector<Array<int,1>> m_numEffectiveMappedNodes;
 struct viewKeyStruct{static string solidMaterialNamesString(){return "material";}};
 void moveParticleWrappersToHost(ParticleManager &){}
 void initialize(int n){m_numberOfSubRegions=n;m_mappedNodes.resize(n);m_mappedFields.resize(n);m_effectiveMappedNodes.resize(n);m_effectiveMappedFields.resize(n);m_shapeFunctionValues.resize(n);m_effectiveShapeFunctionValues.resize(n);m_shapeFunctionGradientValues.resize(n);m_effectiveShapeFunctionGradientValues.resize(n);m_numEffectiveMappedNodes.resize(n);}
"""+splitter.replace('SolidMechanicsMPM::','')+'\n};\n'
    checks=(tests/'subdivision_cpp/test_solver_subdivision.cpp').read_text().replace('// CREATE_FIELDS','\n'.join(create))
    return support+'\n#include "physicsSolvers/solidMechanics/MPMMappingKernels.hpp"\n'+'\n'.join(helpers)+'\n'+solver+checks


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cxx',default='c++')
    parser.add_argument('--sanitize',action='store_true')
    args=parser.parse_args()
    tests=Path(__file__).resolve().parent
    root=next(p for p in tests.parents if (p/'src/coreComponents/mesh/ParticleSubdivision.hpp').is_file())
    with tempfile.TemporaryDirectory(prefix='geos-subdivision-') as folder:
        tmp=Path(folder);(tmp/'common').mkdir()
        (tmp/'common/DataTypes.hpp').write_text('#pragma once\n#define GEOS_HOST_DEVICE\nnamespace geos { using real64=double;using integer=int;using localIndex=int; }\n')
        flags=['-std=c++17','-O1','-g','-Wall','-Wextra','-Werror','-Wshadow','-Wfloat-equal','-pedantic']
        if args.sanitize: flags+=['-fsanitize=address,undefined','-fno-omit-frame-pointer']
        binary=tmp/'geometry'
        subprocess.run([args.cxx,*flags,'-I',str(tmp),'-I',str(root/'src/coreComponents'),str(tests/'subdivision_cpp/test_subdivision.cpp'),'-o',str(binary)],check=True)
        subprocess.run([str(binary)],check=True)
        program=tmp/'solver.cpp';program.write_text(solver_program(root,tests))
        binary=tmp/'solver'
        subprocess.run([args.cxx,*flags,'-I',str(tmp),'-I',str(root/'src/coreComponents'),str(program),'-o',str(binary)],check=True)
        subprocess.run([str(binary)],check=True)

if __name__=='__main__':main()
