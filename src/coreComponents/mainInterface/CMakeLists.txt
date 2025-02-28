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
Package: mainInterface

Contains:
 - GEOS `main` entry point.
 - the class managing the operation flow of the problem being ran in GEOS.
 - basic initialization and environment setup routines.
#]]

#
# Specify all headers
#
set( mainInterface_headers
     GeosxState.hpp
     ProblemManager.hpp
     initialization.hpp
     version.hpp
   )

#
# Specify all sources
#
set( mainInterface_sources
     GeosxState.cpp
     ProblemManager.cpp
     initialization.cpp
     version.cpp
   )

set( dependencyList ${parallelDeps} constitutiveDrivers optionparser )

if( GEOS_ENABLE_FLUIDFLOW )
  list( APPEND dependencyList fluidFlowSolvers )
endif()

if( GEOS_ENABLE_INDUCEDSEISMICITY )
  list( APPEND dependencyList inducedSeismicitySolvers )
endif()

if( GEOS_ENABLE_SIMPLEPDE )
  list( APPEND dependencyLists simplePDESolvers )
endif()

if( GEOS_ENABLE_SOLIDMECHANICS )
  list( APPEND dependencyList solidMechanicsSolvers )
endif()

if( GEOS_ENABLE_SURFACEGENERATION )
  list( APPEND dependencyList surfaceGeneration )
endif()

if( GEOS_ENABLE_WAVEPROPAGATION )
  list( APPEND dependencyList wavePropagationSolvers )
endif()

if( GEOS_ENABLE_MULTIPHYSICS )
  list( APPEND dependencyList multiPhysicsSolvers )
endif()

geos_decorate_link_dependencies( LIST decoratedDependencies
                                 DEPENDENCIES ${dependencyList} )

blt_add_library( NAME       mainInterface
                 SOURCES    ${mainInterface_sources}
                 HEADERS    ${mainInterface_headers}
                 DEPENDS_ON ${decoratedDependencies}
                 OBJECT     ${GEOS_BUILD_OBJ_LIBS}
                 SHARED     ${GEOS_BUILD_SHARED_LIBS}
               )

add_dependencies( mainInterface generate_version )

target_include_directories( mainInterface PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )

install( TARGETS mainInterface LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib )
