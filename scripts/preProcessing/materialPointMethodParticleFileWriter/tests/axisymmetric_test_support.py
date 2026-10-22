"""Extract the shared production mapper and projection bodies for CPU tests.

Array/tensor adapters replace GEOS dependencies; this is not a full GEOS build.
"""
from pathlib import Path

def program_prefix():
    tests = Path(__file__).resolve().parent
    root = next(p for p in tests.parents if (p/'src/coreComponents/mesh/ParticleGeometry.hpp').exists())
    source = (root/'src/coreComponents/physicsSolvers/solidMechanics/SolidMechanicsMPM.cpp').read_text()
    def function(signature):
        start=source.index(signature)
        end=source.index('{',start)+1; depth=1
        while depth:
            depth+=(source[end]=='{')-(source[end]=='}'); end+=1
        return source[start:end].replace('SolidMechanicsMPM::','')
    start=source.index('struct VPHenckyDomainProjectionResult')
    end=source.index('GEOS_HOST_DEVICE\nGEOS_FORCE_INLINE\nint countSubdivisionDirections',start)
    prefix=(tests/'native_geometry_cpp/cpdi_test_support.hpp').read_text()
    prefix+='\n#include "physicsSolvers/solidMechanics/MPMMappingKernels.hpp"\n#include "mesh/ParticleSubdivision.hpp"\n'
    prefix+='\n#include "mesh/AxisymmetricGeometry.hpp"\n'
    prefix+='enum class ParticleType {SinglePoint,CPDI,SinglePointBSpline,CPTI};\n'
    prefix+=source[start:end]
    prefix+='struct Kernel {\n bool m_axisymmetric=false;\n int m_cpdiDomainScaling=0;\n double m_neighborRadius=1.1,m_hEl[3]={1,1,1};\n'
    prefix+=' mpm::CPDIDomainScalingTypeOption m_cpdiDomainScalingType=mpm::CPDIDomainScalingTypeOption::VPHencky;\n'
    for name in ('projectIntegrationGeometry','projectCPTIIntegrationGeometry'):
        prefix+='template<int DIM=3>\n'+function('int SolidMechanicsMPM::'+name+'(')+'\n'
    prefix+='template<int DIM,class NODES,class VALUES,class GRADIENTS>\n'
    prefix+=function('bool SolidMechanicsMPM::computeParticleMapping(')+'\n};\n'
    return root,prefix
