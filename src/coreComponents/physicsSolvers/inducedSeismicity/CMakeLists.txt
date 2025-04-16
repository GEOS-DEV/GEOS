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
Package: inducedSeismicitySolvers
Contains:
#]]

# Specify solver headers
set( inducedSeismicitySolvers_headers
     ExplicitQDRateAndState.hpp
     ImplicitQDRateAndState.hpp
     inducedSeismicityFields.hpp
     QDRateAndStateBase.hpp
     QuasiDynamicEarthQuake.hpp
     rateAndStateFields.hpp
     SeismicityRate.hpp
     SpringSlider.hpp
     kernels/EmbeddedRungeKuttaKernels.hpp
     kernels/ExplicitRateAndStateKernels.hpp
     kernels/ImplicitRateAndStateKernels.hpp
     kernels/RateAndStateKernelsBase.hpp
     kernels/SeismicityRateKernels.hpp )

# Specify solver sources
set( inducedSeismicitySolvers_sources
     ExplicitQDRateAndState.cpp 
     ImplicitQDRateAndState.cpp 
     QDRateAndStateBase.cpp
     QuasiDynamicEarthQuake.cpp
     SeismicityRate.cpp
     SpringSlider.cpp )

set( dependencyList ${parallelDeps} physicsSolversBase )

geos_decorate_link_dependencies( LIST decoratedDependencies DEPENDENCIES ${dependencyList} )

blt_add_library( NAME       inducedSeismicitySolvers
                 SOURCES    ${inducedSeismicitySolvers_sources}
                 HEADERS    ${inducedSeismicitySolvers_headers}
                 DEPENDS_ON ${decoratedDependencies} ${externalComponentDeps}
                 OBJECT     ${GEOS_BUILD_OBJ_LIBS}
                 SHARED     ${GEOS_BUILD_SHARED_LIBS}
               )

target_include_directories( inducedSeismicitySolvers PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )

install( TARGETS inducedSeismicitySolvers LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib )

if( externalComponentDeps )
  target_include_directories( inducedSeismicitySolvers PUBLIC ${CMAKE_SOURCE_DIR}/externalComponents )
endif()    
