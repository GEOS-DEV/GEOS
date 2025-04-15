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
Package: surfaceGeneration
Contains:
#]]

# Specify solver headers
set( surfaceGeneration_headers
     EmbeddedSurfaceGenerator.hpp
     EmbeddedSurfacesParallelSynchronization.hpp
     LogLevelsInfo.hpp
     ParallelTopologyChange.hpp
     SurfaceGenerator.hpp
     SurfaceGeneratorFields.hpp
     LogLevelsInfo.hpp
     kernels/surfaceGenerationKernels.hpp
     kernels/surfaceGenerationKernelsHelpers.hpp )

# Specify solver sources
set( surfaceGeneration_sources
     EmbeddedSurfaceGenerator.cpp
     EmbeddedSurfacesParallelSynchronization.cpp
     ParallelTopologyChange.cpp
     ParallelTopologyChangeNoFixup.cpp
     SurfaceGenerator.cpp )

set( dependencyList ${parallelDeps} physicsSolversBase )

geos_decorate_link_dependencies( LIST decoratedDependencies DEPENDENCIES ${dependencyList} )

blt_add_library( NAME       surfaceGeneration
                 SOURCES    ${surfaceGeneration_sources}
                 HEADERS    ${surfaceGeneration_headers}
                 DEPENDS_ON ${decoratedDependencies} ${externalComponentDeps}
                 OBJECT     ${GEOS_BUILD_OBJ_LIBS}
                 SHARED     ${GEOS_BUILD_SHARED_LIBS}
               )

target_include_directories( surfaceGeneration PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )

install( TARGETS surfaceGeneration LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib )

if( externalComponentDeps )
  target_include_directories( surfaceGeneration PUBLIC ${CMAKE_SOURCE_DIR}/externalComponents )
endif()