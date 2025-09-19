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
Package: gravitySolvers
Contains:
#]]

# Specify solver headers
set( gravitySolvers_headers
     GravityLogLevelsInfo.hpp
     GravityFields.hpp
     GravitySolverBase.hpp
     GravityFE.hpp
     GravityFEKernel.hpp )

# Specify solver sources
set( gravitySolvers_sources
     GravitySolverBase.cpp
     GravityFE.cpp )

set( dependencyList ${parallelDeps} physicsSolversBase )

# Conditionally add compositional multiphase support if fluid flow is enabled
if( GEOS_ENABLE_FLUIDFLOW AND TARGET fluidFlowSolvers )
  list( APPEND gravitySolvers_headers GravityFE_CompositionalMultiphaseBase.hpp )
  list( APPEND gravitySolvers_sources GravityFE_CompositionalMultiphaseBase.cpp )
  list( APPEND dependencyList fluidFlowSolvers )
  
  # Define a preprocessor macro so the code knows fluid flow is available
  add_definitions( -DGEOS_GRAVITY_WITH_FLUIDFLOW )
endif()


geos_decorate_link_dependencies( LIST decoratedDependencies DEPENDENCIES ${dependencyList} )

blt_add_library( NAME       gravitySolvers
                 SOURCES    ${gravitySolvers_sources}
                 HEADERS    ${gravitySolvers_headers}
                 DEPENDS_ON ${decoratedDependencies} ${externalComponentDeps}
                 OBJECT     ${GEOS_BUILD_OBJ_LIBS}
                 SHARED     ${GEOS_BUILD_SHARED_LIBS}
               )

target_include_directories( gravitySolvers PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )

install( TARGETS gravitySolvers LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib )

if( externalComponentDeps )
  target_include_directories( gravitySolvers PUBLIC ${CMAKE_SOURCE_DIR}/externalComponents )
endif()    
