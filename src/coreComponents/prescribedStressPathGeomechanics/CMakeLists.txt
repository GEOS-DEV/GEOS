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
Package : prescribedStressPathGeomechanics

Contains classes for computing the effective normal stress over the fracture plane analytically.
#]]

#
# Specify all headers
#
set( prescribedStressPathGeomechanics_headers
     OedometricStressPath.hpp
   )

#
# Specify all sources
#
set( prescribedStressPathGeomechanics_sources
     OedometricStressPath.cpp
   )

set( dependencyList ${parallelDeps} common dataRepository )

geos_decorate_link_dependencies( LIST decoratedDependencies
                                 DEPENDENCIES ${dependencyList} )

blt_add_library( NAME       prescribedStressPathGeomechanics
                 SOURCES    ${prescribedStressPathGeomechanics_sources}
                 HEADERS    ${prescribedStressPathGeomechanics_headers}
                 DEPENDS_ON ${decoratedDependencies}
                 OBJECT     ${GEOS_BUILD_OBJ_LIBS}
                 SHARED     ${GEOS_BUILD_SHARED_LIBS}
               )

target_include_directories( prescribedStressPathGeomechanics PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )

install( TARGETS prescribedStressPathGeomechanics LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib )
