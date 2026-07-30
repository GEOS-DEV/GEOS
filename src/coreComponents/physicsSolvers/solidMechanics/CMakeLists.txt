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

if( NOT GEOS_ENABLE_SOLIDMECHANICS_MPM AND NOT GEOS_ENABLE_SOLIDMECHANICS_FEM )
  message( FATAL_ERROR "GEOS_ENABLE_SOLIDMECHANICS requires at least one of GEOS_ENABLE_SOLIDMECHANICS_MPM or GEOS_ENABLE_SOLIDMECHANICS_FEM" )
endif()

# SolidMechanicsFields is shared by MPM and FEM. State reset and statistics
# currently depend on SolidMechanicsLagrangianFEM and are added with FEM below.
set( solidMechanicsSolvers_headers
     SolidMechanicsFields.hpp )

set( solidMechanicsSolvers_sources )

if( GEOS_ENABLE_SOLIDMECHANICS_MPM )
  list( APPEND solidMechanicsSolvers_headers
        SolidMechanicsMPM.hpp
        MPMSolverFields.hpp
        kernels/ExplicitMPM.hpp )

  list( APPEND solidMechanicsSolvers_sources
        SolidMechanicsMPM.cpp )
endif()

if( GEOS_ENABLE_SOLIDMECHANICS_FEM )
  list( APPEND solidMechanicsSolvers_headers
        SolidMechanicsLagrangianFEM.hpp
        SolidMechanicsStateReset.hpp
        SolidMechanicsStatistics.hpp
        kernels/SolidMechanicsLagrangianFEMKernels.hpp
        kernels/ExplicitFiniteStrain.hpp
        kernels/ExplicitFiniteStrain_impl.hpp
        kernels/ExplicitSmallStrain.hpp
        kernels/ExplicitSmallStrain_impl.hpp
        kernels/FixedStressThermoPoromechanics.hpp
        kernels/FixedStressThermoPoromechanics_impl.hpp
        kernels/ImplicitSmallStrainNewmark.hpp
        kernels/ImplicitSmallStrainNewmark_impl.hpp
        kernels/ImplicitSmallStrainQuasiStatic.hpp
        kernels/ImplicitSmallStrainQuasiStatic_impl.hpp )

  list( APPEND solidMechanicsSolvers_sources
        SolidMechanicsLagrangianFEM.cpp
        SolidMechanicsStateReset.cpp
        SolidMechanicsStatistics.cpp )
endif()

if( GEOS_ENABLE_CONTACT )
  if( NOT GEOS_ENABLE_SOLIDMECHANICS_FEM )
    message( FATAL_ERROR "GEOS_ENABLE_CONTACT requires GEOS_ENABLE_SOLIDMECHANICS_FEM=ON" )
  endif()

  list( APPEND solidMechanicsSolvers_headers
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

  list( APPEND solidMechanicsSolvers_sources
        contact/ContactSolverBase.cpp
        contact/SolidMechanicsEmbeddedFractures.cpp
        contact/SolidMechanicsLagrangeContact.cpp
        contact/SolidMechanicsLagrangeContactBubbleStab.cpp
        contact/SolidMechanicsAugmentedLagrangianContact.cpp )
endif()

set( dependencyList ${parallelDeps} physicsSolversBase )

geos_decorate_link_dependencies( LIST decoratedDependencies DEPENDENCIES ${dependencyList} )

if( GEOS_ENABLE_SOLIDMECHANICS_FEM )
  file( READ "${CMAKE_CURRENT_LIST_DIR}/kernelSpecs.json" kernelSpecs )
  set( kernelTemplateFileList
       kernels/SolidMechanicsKernels.cpp.template
       kernels/SolidMechanicsFixedStressThermoPoromechanicsKernels.cpp.template )

  foreach( kernelTemplateFile ${kernelTemplateFileList} )
    get_filename_component( jsonKey ${kernelTemplateFile} NAME_WE )
    generateKernels( TEMPLATE ${kernelTemplateFile}
                     JSON kernelSpecs
                     KEY ${jsonKey}
                     HEADERS headerFiles
                     SOURCES sourceFiles )

    list( APPEND solidMechanicsSolvers_headers ${headerFiles} )
    list( APPEND solidMechanicsSolvers_sources ${sourceFiles} )
  endforeach()
endif()

option(
  GEOS_MPM_FORCE_O0
  "Compile SolidMechanicsMPM.cpp with -O0 for GPU diagnosis"
  OFF
)

if( GEOS_MPM_FORCE_O0 AND GEOS_ENABLE_SOLIDMECHANICS_MPM )
  set_property(
    SOURCE SolidMechanicsMPM.cpp
    APPEND
    PROPERTY COMPILE_OPTIONS "-O0"
  )
endif()

blt_add_library( NAME       solidMechanicsSolvers
                 SOURCES    ${solidMechanicsSolvers_sources}
                 HEADERS    ${solidMechanicsSolvers_headers}
                 DEPENDS_ON ${decoratedDependencies} ${externalComponentDeps}
                 OBJECT     ${GEOS_BUILD_OBJ_LIBS}
                 SHARED     ${GEOS_BUILD_SHARED_LIBS}
               )

target_include_directories( solidMechanicsSolvers PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )

if( GEOS_ENABLE_MPM_GPU_DEBUG )
  # SolidMechanicsMPM.cpp instantiates the RAJA/HIP kernels defined in its
  # included headers.  Preserve source-level device information and avoid
  # optimization reordering while diagnosing a fault.
  target_compile_options( solidMechanicsSolvers PRIVATE
                          "$<$<COMPILE_LANGUAGE:CXX>:-O0;-ggdb;-fno-omit-frame-pointer>"
                          "$<$<COMPILE_LANGUAGE:HIP>:-O0;-ggdb;-fno-omit-frame-pointer>" )
endif()

install( TARGETS solidMechanicsSolvers LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib )

if( externalComponentDeps )
  target_include_directories( solidMechanicsSolvers PUBLIC ${CMAKE_SOURCE_DIR}/externalComponents )
endif()
