# SPDX-License-Identifier: LGPL-2.1-only
#
# Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
# Copyright (c) 2018-2024 TotalEnergies
# Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
# Copyright (c) 2023-2024 Chevron
# Copyright (c) 2019-     GEOS/GEOSX Contributors
# All rights reserved
#
# See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
#
#--------------------------------------------------------------------------------------------------

#[[
Package: solidMechanicsSolvers

Contains:
#]]

# Specify solver headers
set( solidMechanicsSolvers_headers
     SolidMechanicsFields.hpp
     SolidMechanicsLagrangianFEM.hpp
     SolidMechanicsLagrangianFEM.hpp
     kernels/SolidMechanicsLagrangianFEMKernels.hpp
     SolidMechanicsMPM.hpp
     MPMSolverFields.hpp
     kernels/ExplicitFiniteStrain.hpp
     kernels/ExplicitFiniteStrain_impl.hpp
     kernels/ExplicitMPM.hpp
     kernels/ExplicitSmallStrain.hpp
     kernels/ExplicitSmallStrain_impl.hpp
     kernels/FixedStressThermoPoromechanics.hpp
     kernels/FixedStressThermoPoromechanics_impl.hpp
     kernels/ImplicitSmallStrainNewmark.hpp
     kernels/ImplicitSmallStrainNewmark_impl.hpp
     kernels/ImplicitSmallStrainQuasiStatic.hpp
     kernels/ImplicitSmallStrainQuasiStatic_impl.hpp
     SolidMechanicsStateReset.hpp
     SolidMechanicsStatistics.hpp
     contact/ContactSolverBase.hpp
     contact/FractureState.hpp
     contact/ContactFields.hpp
     contact/SolidMechanicsEmbeddedFractures.hpp
     contact/SolidMechanicsLagrangeContact.hpp
     contact/SolidMechanicsLagrangeContactBubbleStab.hpp
     contact/SolidMechanicsAugmentedLagrangianContact.hpp
     contact/kernels/SolidMechanicsConformingContactKernelsBase.hpp
     contact/kernels/SolidMechanicsDisplacementJumpUpdateKernels.hpp
     contact/kernels/SolidMechanicsEFEMKernelsBase.hpp
     contact/kernels/SolidMechanicsEFEMKernels.hpp
     contact/kernels/SolidMechanicsEFEMStaticCondensationKernels.hpp
     contact/kernels/SolidMechanicsEFEMKernelsHelper.hpp
     contact/kernels/SolidMechanicsALMKernelsBase.hpp
     contact/kernels/SolidMechanicsALMKernels.hpp
     contact/kernels/SolidMechanicsConformingContactKernelsHelper.hpp
     contact/kernels/SolidMechanicsContactFaceBubbleKernels.hpp
     contact/kernels/SolidMechanicsLagrangeContactKernels.hpp 
     contact/LogLevelsInfo.hpp )

# Specify solver sources
set( solidMechanicsSolvers_sources
     SolidMechanicsLagrangianFEM.cpp
     SolidMechanicsMPM.cpp
     SolidMechanicsStateReset.cpp
     SolidMechanicsStatistics.cpp
     contact/ContactSolverBase.cpp
     contact/SolidMechanicsEmbeddedFractures.cpp
     contact/SolidMechanicsLagrangeContact.cpp
     contact/SolidMechanicsLagrangeContactBubbleStab.cpp
     contact/SolidMechanicsAugmentedLagrangianContact.cpp )

#include( kernels/SolidMechanicsKernels.cmake)

set( dependencyList ${parallelDeps} physicsSolversBase )

geos_decorate_link_dependencies( LIST decoratedDependencies DEPENDENCIES ${dependencyList} )


file( READ "${CMAKE_CURRENT_LIST_DIR}/kernelSpecs.json" kernelSpecs )
set( kernelTemplateFileList "" )

  list( APPEND kernelTemplateFileList 
        kernels/SolidMechanicsKernels.cpp.template
        kernels/SolidMechanicsFixedStressThermoPoromechanicsKernels.cpp.template )


foreach( kernelTemplateFile ${kernelTemplateFileList} )
  get_filename_component( jsonKey ${kernelTemplateFile} NAME_WE )
  generateKernels( TEMPLATE ${kernelTemplateFile}
                   JSON kernelSpecs
                   KEY ${jsonKey}
                   HEADERS headerFiles
                   SOURCES sourceFiles )

list(APPEND solidMechanicsSolvers_headers ${headerFiles})
list(APPEND solidMechanicsSolvers_sources ${sourceFiles})
endforeach()


blt_add_library( NAME       solidMechanicsSolvers
                 SOURCES    ${solidMechanicsSolvers_sources}
                 HEADERS    ${solidMechanicsSolvers_headers}
                 DEPENDS_ON ${decoratedDependencies} ${externalComponentDeps}
                 OBJECT     ${GEOS_BUILD_OBJ_LIBS}
                 SHARED     ${GEOS_BUILD_SHARED_LIBS}
               )

target_include_directories( solidMechanicsSolvers PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )

install( TARGETS solidMechanicsSolvers LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib )

if( externalComponentDeps )
  target_include_directories( solidMechanicsSolvers PUBLIC ${CMAKE_SOURCE_DIR}/externalComponents )
endif()