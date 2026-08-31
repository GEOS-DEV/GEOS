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
Package: mixedVEM

Lowest-order mixed virtual element discretization of linear elasticity on polyhedral
meshes, based on the Hellinger-Reissner principle. The package provides the element
level operators B_E and K_E and the helpers that scatter them into the global
saddle-point system.
#]]

#
# Specify all headers
#
set( mixedVEM_headers
     HybridMixedVEM.hpp
     MixedVEMAssembly.hpp
     MixedVEMCellOutput.hpp
     MixedVEMCompliance.hpp
     MixedVEMDiscretization.hpp
     MixedVEMElementOperators.hpp
     MixedVEMFaceBasis.hpp
     MixedVEMFields.hpp
     MixedVEMGeometry.hpp
     MixedVEMManager.hpp
     MixedVEMTypes.hpp
   )

#
# Specify all sources
#
set( mixedVEM_sources
     MixedVEMCompliance.cpp
     MixedVEMDiscretization.cpp
     MixedVEMManager.cpp
   )

set( dependencyList ${parallelDeps} mesh denseLinearAlgebra )

geos_decorate_link_dependencies( LIST decoratedDependencies
                                 DEPENDENCIES ${dependencyList} )

blt_add_library( NAME       mixedVEM
                 SOURCES    ${mixedVEM_sources}
                 HEADERS    ${mixedVEM_headers}
                 DEPENDS_ON ${decoratedDependencies}
                 OBJECT     ${GEOS_BUILD_OBJ_LIBS}
                 SHARED     ${GEOS_BUILD_SHARED_LIBS}
               )

target_include_directories( mixedVEM PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )

install( TARGETS mixedVEM LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib )

if( GEOS_ENABLE_TESTS )
  add_subdirectory( unitTests )
endif()
