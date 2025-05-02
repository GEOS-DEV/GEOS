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
Package: physicsSolvers

Contains:
 - physics solvers base and manager classes.
 - implementations of different physics models.
 - physics solver wrapper designed for PyGEOSX interface.
#]]

#
# Specify solver headers
#
set( physicsSolversBase_headers
     LinearSolverParameters.hpp
     NonlinearSolverParameters.hpp
     PhysicsSolverBase.hpp
     PhysicsSolverBaseKernels.hpp
     SolverStatistics.hpp
     FieldStatisticsBase.hpp )
#
# Specify solver sources
#
set( physicsSolversBase_sources
     LinearSolverParameters.cpp
     NonlinearSolverParameters.cpp
     PhysicsSolverBase.cpp
     SolverStatistics.cpp )

set( dependencyList ${parallelDeps} fileIO discretizationMethods events linearAlgebra )

if( ENABLE_PYGEOSX )
  list( APPEND physicsSolversBase_headers
  python/PySolver.hpp
  python/PySolverType.hpp )
  list( APPEND physicsSolversBase_sources
        python/PySolver.cpp )
  list( APPEND dependencyList Python3::Python pylvarray )
endif()

geos_decorate_link_dependencies( LIST decoratedDependencies
                                 DEPENDENCIES ${dependencyList} )

blt_add_library( NAME       physicsSolversBase
                 SOURCES    ${physicsSolversBase_sources}
                 HEADERS    ${physicsSolversBase_headers}
                 DEPENDS_ON ${decoratedDependencies} ${externalComponentDeps}
                 OBJECT     ${GEOS_BUILD_OBJ_LIBS}
                 SHARED     ${GEOS_BUILD_SHARED_LIBS}
               )

target_include_directories( physicsSolversBase PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )

install( TARGETS physicsSolversBase LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib )

if( externalComponentDeps )
  target_include_directories( physicsSolversBase PUBLIC ${CMAKE_SOURCE_DIR}/externalComponents )
endif()

# PhysicsSolvers interface
#
# Specify solver headers
#
set( physicsSolvers_headers
     PhysicsSolverManager.hpp)
#
# Specify solver sources
#
set( physicsSolvers_sources
     PhysicsSolverManager.cpp )

set( kernelTemplateFileList )


set( dependencyList ${parallelDeps})
if( GEOS_ENABLE_FLUIDFLOW )
  add_subdirectory( fluidFlow )
  list( APPEND dependencyList fluidFlowSolvers )
  list( APPEND kernelTemplateFileList
        fluidFlow/SinglePhaseHybridFVMKernels.cpp.template
        fluidFlow/CompositionalMultiphaseHybridFVMKernels_upwinding.cpp.template
        fluidFlow/CompositionalMultiphaseHybridFVMKernels_assembly.cpp.template
        fluidFlow/CompositionalMultiphaseHybridFVMKernels_flux.cpp.template
        fluidFlow/ReactiveCompositionalMultiphaseOBLKernels.cpp.template )
endif()

if( GEOS_ENABLE_INDUCEDSEISMICITY )
  add_subdirectory( inducedSeismicity )
  list( APPEND dependencyList inducedSeismicitySolvers )
endif()

if( GEOS_ENABLE_SIMPLEPDE )
  add_subdirectory( simplePDE )
  list( APPEND dependencyLists simplePDESolvers )
endif()

if( GEOS_ENABLE_SOLIDMECHANICS )
  add_subdirectory( solidMechanics )
  list( APPEND dependencyList solidMechanicsSolvers )
  list( APPEND kernelTemplateFileList 
        solidMechanics/kernels/SolidMechanicsKernels.cpp.template
        solidMechanics/kernels/SolidMechanicsFixedStressThermoPoromechanicsKernels.cpp.template )
endif()

if( GEOS_ENABLE_SURFACEGENERATION )
  add_subdirectory( surfaceGeneration )
  list( APPEND dependencyList surfaceGeneration )
endif()

if( GEOS_ENABLE_WAVEPROPAGATION )
  add_subdirectory( wavePropagation )
  list( APPEND dependencyList wavePropagationSolvers )
endif()

if( GEOS_ENABLE_MULTIPHYSICS )
  add_subdirectory( multiphysics )
  list( APPEND dependencyList multiPhysicsSolvers )
  list( APPEND kernelTemplateFileList 
        multiphysics/poromechanicsKernels/PoromechanicsKernels.cpp.template
        multiphysics/poromechanicsKernels/PoromechanicsEFEMKernels.cpp.template
        multiphysics/poromechanicsKernels/ThermoPoromechanicsKernels.cpp.template
        multiphysics/poromechanicsKernels/PoromechanicsDamageKernels.cpp.template )

endif()

if(NOT DEFINED GEOS_KERNEL_SPEC)
  set(GEOS_KERNEL_SPEC "${CMAKE_CURRENT_LIST_DIR}/kernels-all.json" )
endif()
  
file( READ ${GEOS_KERNEL_SPEC} kernelSpecs )

 
foreach( kernelTemplateFile ${kernelTemplateFileList} )
  get_filename_component( jsonKey ${kernelTemplateFile} NAME_WE )
  generateKernels( TEMPLATE ${kernelTemplateFile}
                   JSON kernelSpecs
                   KEY ${jsonKey}
                   HEADERS headerFiles
                   SOURCES sourceFiles )

  list(APPEND physicsSolvers_headers ${headerFiles})
  list(APPEND physicsSolvers_sources ${sourceFiles})
endforeach()
 
# TODO: The two kernels below have non-matching file names and JSON keys.
#       Either fix them to follow pattern, or come up with another mechanism.

generateKernels( TEMPLATE "fluidFlow/CompositionalMultiphaseFVMKernel_isothermal.cpp.template"
                JSON kernelSpecs
                KEY "IsothermalCompositionalMultiphaseFVM_DirichletFluxComputeKernels"
                HEADERS headerFiles
                SOURCES sourceFiles )
list(APPEND physicsSolvers_headers ${headerFiles})
list(APPEND physicsSolvers_sources ${sourceFiles})

generateKernels( TEMPLATE "fluidFlow/CompositionalMultiphaseFVMKernel_thermal.cpp.template"
                JSON kernelSpecs
                KEY "ThermalCompositionalMultiphaseFVM_DirichletFluxComputeKernels"
                HEADERS headerFiles
                SOURCES sourceFiles )
list(APPEND physicsSolvers_headers ${headerFiles})
list(APPEND physicsSolvers_sources ${sourceFiles})


geos_decorate_link_dependencies( LIST decoratedDependencies
                                 DEPENDENCIES ${dependencyList} )

blt_add_library( NAME       physicsSolvers
                 SOURCES    ${physicsSolvers_sources}
                 HEADERS    ${physicsSolvers_headers}
                 DEPENDS_ON ${decoratedDependencies} ${externalComponentDeps}
                 OBJECT     ${GEOS_BUILD_OBJ_LIBS}
                 SHARED     ${GEOS_BUILD_SHARED_LIBS}
               )

target_include_directories( physicsSolvers PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  target_link_options(physicsSolvers PUBLIC -Wl,--no-as-needed)
endif()
              
install( TARGETS physicsSolvers LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib )