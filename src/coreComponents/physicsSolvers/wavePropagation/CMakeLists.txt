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
Package: wavePropagationSolvers
Contains:
#]]

# Specify solver headers
set( wavePropagationSolvers_headers
     LogLevelsInfo.hpp
     shared/WaveSolverBase.hpp
     shared/WaveSolverUtils.hpp
     shared/PrecomputeSourcesAndReceiversKernel.hpp
     sem/acoustic/shared/AcousticFields.hpp
     sem/acoustic/secondOrderEqn/isotropic/AcousticWaveEquationSEM.hpp
     sem/acoustic/secondOrderEqn/isotropic/AcousticWaveEquationSEMKernel.hpp
     sem/acoustic/shared/AcousticTimeSchemeSEMKernel.hpp
     sem/acoustic/shared/AcousticMatricesSEMKernel.hpp
     sem/acoustic/secondOrderEqn/isotropic/AcousticPMLSEMKernel.hpp
     sem/elastic/shared/ElasticFields.hpp
     sem/elastic/secondOrderEqn/isotropic/ElasticWaveEquationSEM.hpp
     sem/elastic/secondOrderEqn/isotropic/ElasticWaveEquationSEMKernel.hpp
     sem/elastic/shared/ElasticTimeSchemeSEMKernel.hpp
     sem/elastic/firstOrderEqn/isotropic/ElasticFirstOrderWaveEquationSEM.hpp
     sem/elastic/firstOrderEqn/isotropic/ElasticFirstOrderWaveEquationSEMKernel.hpp
     sem/acoustic/firstOrderEqn/isotropic/AcousticFirstOrderWaveEquationSEM.hpp
     sem/acoustic/firstOrderEqn/isotropic/AcousticFirstOrderWaveEquationSEMKernel.hpp
     sem/elastic/shared/ElasticMatricesSEMKernel.hpp
     sem/elastic/secondOrderEqn/anisotropic/ElasticVTIFields.hpp
     sem/elastic/secondOrderEqn/anisotropic/ElasticVTIWaveEquationSEMKernel.hpp
     sem/elastic/secondOrderEqn/anisotropic/ElasticTTIFields.hpp
     sem/elastic/secondOrderEqn/anisotropic/ElasticTTIWaveEquationSEMKernel.hpp
     sem/acoustic/secondOrderEqn/anisotropic/AcousticVTIFields.hpp
     sem/acoustic/secondOrderEqn/anisotropic/AcousticVTIWaveEquationSEM.hpp
     sem/acoustic/secondOrderEqn/anisotropic/AcousticVTIWaveEquationSEMKernel.hpp
     sem/acoustoelastic/secondOrderEqn/isotropic/AcoustoElasticFields.hpp
     sem/acoustoelastic/secondOrderEqn/isotropic/AcoustoElasticTimeSchemeSEMKernel.hpp
     sem/acoustoelastic/secondOrderEqn/isotropic/AcousticElasticWaveEquationSEM.hpp
     sem/acoustoelastic/secondOrderEqn/isotropic/AcousticElasticWaveEquationSEMKernel.hpp )

# Specify solver sources
set( wavePropagationSolvers_sources
     shared/WaveSolverBase.cpp
     sem/acoustic/secondOrderEqn/isotropic/AcousticWaveEquationSEM.cpp
     sem/elastic/secondOrderEqn/isotropic/ElasticWaveEquationSEM.cpp
     sem/elastic/firstOrderEqn/isotropic/ElasticFirstOrderWaveEquationSEM.cpp
     sem/acoustic/firstOrderEqn/isotropic/AcousticFirstOrderWaveEquationSEM.cpp
     sem/acoustic/secondOrderEqn/anisotropic/AcousticVTIWaveEquationSEM.cpp
     sem/acoustoelastic/secondOrderEqn/isotropic/AcousticElasticWaveEquationSEM.cpp )

set( dependencyList ${parallelDeps} physicsSolversBase )

geos_decorate_link_dependencies( LIST decoratedDependencies DEPENDENCIES ${dependencyList} )

blt_add_library( NAME       wavePropagationSolvers
                 SOURCES    ${wavePropagationSolvers_sources}
                 HEADERS    ${wavePropagationSolvers_headers}
                 DEPENDS_ON ${decoratedDependencies} ${externalComponentDeps}
                 OBJECT     ${GEOS_BUILD_OBJ_LIBS}
                 SHARED     ${GEOS_BUILD_SHARED_LIBS}
               )

target_include_directories( wavePropagationSolvers PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )

install( TARGETS wavePropagationSolvers LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib )

if( externalComponentDeps )
  target_include_directories( wavePropagationSolvers PUBLIC ${CMAKE_SOURCE_DIR}/externalComponents )
endif()
