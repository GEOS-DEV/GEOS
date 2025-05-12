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
Package: multiPhysicsSolvers
Contains:
#]]


# Specify solver headers
set( multiPhysicsSolvers_headers
    CompositionalMultiphaseReservoirAndWells.hpp
    CoupledReservoirAndWellsBase.hpp
    CoupledSolver.hpp
    PoromechanicsSolver.hpp
    FlowProppantTransportSolver.hpp
    HydrofractureSolver.hpp
    HydrofractureSolverKernels.hpp
    MultiphasePoromechanics.hpp
    OneWayCoupledFractureFlowContactMechanics.hpp
    MultiphasePoromechanicsConformingFractures.hpp
    PhaseFieldFractureSolver.hpp
    PhaseFieldPoromechanicsSolver.hpp
    PoromechanicsInitialization.hpp
    PoromechanicsFields.hpp
    PoromechanicsInitialization.hpp
    poromechanicsKernels/MultiphasePoromechanics.hpp
    poromechanicsKernels/MultiphasePoromechanics_impl.hpp
    poromechanicsKernels/MultiphasePoromechanicsConformingFractures.hpp
    poromechanicsKernels/PoromechanicsBase.hpp
    poromechanicsKernels/SinglePhasePoromechanics.hpp
    poromechanicsKernels/SinglePhasePoromechanics_impl.hpp
    poromechanicsKernels/SinglePhasePoromechanicsConformingFractures.hpp
    poromechanicsKernels/SinglePhasePoromechanicsDamage.hpp
    poromechanicsKernels/SinglePhasePoromechanicsDamage_impl.hpp
    poromechanicsKernels/SinglePhasePoromechanicsEFEM.hpp
    poromechanicsKernels/SinglePhasePoromechanicsEFEM_impl.hpp
    poromechanicsKernels/SinglePhasePoromechanicsFractures.hpp
    poromechanicsKernels/SinglePhasePoromechanicsEmbeddedFractures.hpp
    poromechanicsKernels/ThermalMultiphasePoromechanics.hpp
    poromechanicsKernels/ThermalMultiphasePoromechanics_impl.hpp
    poromechanicsKernels/ThermalSinglePhasePoromechanics.hpp
    poromechanicsKernels/ThermalSinglePhasePoromechanics_impl.hpp
    poromechanicsKernels/ThermalSinglePhasePoromechanicsEFEM.hpp
    poromechanicsKernels/ThermalSinglePhasePoromechanicsEFEM_impl.hpp
    poromechanicsKernels/ThermalSinglePhasePoromechanicsConformingFractures.hpp
    poromechanicsKernels/ThermalSinglePhasePoromechanicsEmbeddedFractures.hpp
    SinglePhasePoromechanics.hpp
    SinglePhasePoromechanicsEmbeddedFractures.hpp
    SinglePhasePoromechanicsConformingFractures.hpp
    SinglePhaseReservoirAndWells.hpp )

# Specify solver sources
set( multiPhysicsSolvers_sources
    CompositionalMultiphaseReservoirAndWells.cpp
    CoupledReservoirAndWellsBase.cpp
    FlowProppantTransportSolver.cpp
    HydrofractureSolver.cpp
    MultiphasePoromechanics.cpp
    OneWayCoupledFractureFlowContactMechanics.cpp
    MultiphasePoromechanicsConformingFractures.cpp
    PhaseFieldFractureSolver.cpp
    PhaseFieldPoromechanicsSolver.cpp
    PoromechanicsInitialization.cpp
    SinglePhasePoromechanics.cpp
    SinglePhasePoromechanicsEmbeddedFractures.cpp
    SinglePhasePoromechanicsConformingFractures.cpp
    SinglePhaseReservoirAndWells.cpp )

#include(poromechanicsKernels/PoromechanicsKernels.cmake)

set( dependencyList ${parallelDeps} fluidFlowSolvers simplePDESolvers solidMechanicsSolvers )

geos_decorate_link_dependencies( LIST decoratedDependencies DEPENDENCIES ${dependencyList} )


file( READ "${CMAKE_CURRENT_LIST_DIR}/kernelSpecs.json" kernelSpecs )
set( kernelTemplateFileList "" )

  list( APPEND kernelTemplateFileList 
        poromechanicsKernels/PoromechanicsKernels.cpp.template
        poromechanicsKernels/PoromechanicsEFEMKernels.cpp.template
        poromechanicsKernels/ThermoPoromechanicsKernels.cpp.template
        poromechanicsKernels/PoromechanicsDamageKernels.cpp.template )


foreach( kernelTemplateFile ${kernelTemplateFileList} )
  get_filename_component( jsonKey ${kernelTemplateFile} NAME_WE )
  generateKernels( TEMPLATE ${kernelTemplateFile}
                   JSON kernelSpecs
                   KEY ${jsonKey}
                   HEADERS headerFiles
                   SOURCES sourceFiles )

list(APPEND multiPhysicsSolvers_headers ${headerFiles})
list(APPEND multiPhysicsSolvers_sources ${sourceFiles})
endforeach()

blt_add_library( NAME       multiPhysicsSolvers
                 SOURCES    ${multiPhysicsSolvers_sources}
                 HEADERS    ${multiPhysicsSolvers_headers}
                 DEPENDS_ON ${decoratedDependencies} ${externalComponentDeps}
                 OBJECT     ${GEOS_BUILD_OBJ_LIBS}
                 SHARED     ${GEOS_BUILD_SHARED_LIBS}
               )

target_include_directories( multiPhysicsSolvers PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )

install( TARGETS multiPhysicsSolvers LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib )

if( externalComponentDeps )
  target_include_directories( multiPhysicsSolvers PUBLIC ${CMAKE_SOURCE_DIR}/externalComponents )
endif()
