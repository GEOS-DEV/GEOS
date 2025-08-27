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
     GravitySolverBase.hpp
     GravityFE.hpp
     GravityFE_CompositionalMultiphaseFVM.hpp
     GravityFEKernel.hpp )

# Specify solver sources
set( gravitySolvers_sources
     GravitySolverBase.cpp
     GravityFE_CompositionalMultiphaseFVM.cpp
     GravityFE.cpp )

set( dependencyList ${parallelDeps} physicsSolversBase )

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
