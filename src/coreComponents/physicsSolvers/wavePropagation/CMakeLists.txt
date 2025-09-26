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
     shared/TaperKernel.hpp
     sem/acoustic/shared/AcousticFields.hpp
     sem/acoustic/secondOrderEqn/anisotropic/AcousticSecondOrderAnisotropicWaveEquationSEM.hpp
     sem/acoustic/secondOrderEqn/anisotropic/AcousticSecondOrderAnisotropicVTISEM.hpp
     sem/acoustic/secondOrderEqn/anisotropic/AcousticSecondOrderAnisotropicTTISEM.hpp
     sem/acoustic/secondOrderEqn/anisotropic/AcousticSecondOrderVTIFletcher.hpp
     sem/acoustic/secondOrderEqn/anisotropic/AcousticSecondOrderTTIFletcher.hpp
     sem/acoustic/secondOrderEqn/anisotropic/AcousticSecondOrderVTIZhang.hpp
     sem/acoustic/secondOrderEqn/anisotropic/AcousticSecondOrderTTIZhang.hpp
     sem/acoustic/secondOrderEqn/anisotropic/DampingMatrixComputers.hpp
     sem/acoustic/secondOrderEqn/isotropic/AcousticWaveEquationSEM.hpp
     sem/acoustic/secondOrderEqn/isotropic/AcousticWaveEquationSEMKernel.hpp
     sem/acoustic/secondOrderEqn/anisotropic/AcousticTTIFields.hpp
     sem/acoustic/secondOrderEqn/anisotropic/AcousticTTIFletcherWaveEquationSEM.hpp
     sem/acoustic/secondOrderEqn/anisotropic/AcousticTTIFletcherWaveEquationSEMKernel.hpp
     sem/acoustic/secondOrderEqn/anisotropic/AcousticTTIZhangWaveEquationSEM.hpp
     sem/acoustic/secondOrderEqn/anisotropic/AcousticTTIZhangWaveEquationSEMKernel.hpp
     sem/acoustic/secondOrderEqn/anisotropic/AcousticVTIFields.hpp
     sem/acoustic/secondOrderEqn/anisotropic/AcousticVTIFletcherAdjointWaveEquationSEMKernel.hpp
     sem/acoustic/secondOrderEqn/anisotropic/AcousticVTIFletcherWaveEquationSEM.hpp
     sem/acoustic/secondOrderEqn/anisotropic/AcousticVTIFletcherWaveEquationSEMKernel.hpp
     sem/acoustic/secondOrderEqn/anisotropic/AcousticVTIZhangAdjointWaveEquationSEMKernel.hpp
     sem/acoustic/secondOrderEqn/anisotropic/AcousticVTIZhangWaveEquationSEM.hpp
     sem/acoustic/secondOrderEqn/anisotropic/AcousticVTIZhangWaveEquationSEMKernel.hpp     
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
     sem/acoustoelastic/secondOrderEqn/isotropic/AcoustoElasticFields.hpp
     sem/acoustoelastic/secondOrderEqn/isotropic/AcoustoElasticTimeSchemeSEMKernel.hpp
     sem/acoustoelastic/secondOrderEqn/isotropic/AcousticElasticWaveEquationSEM.hpp
     sem/acoustoelastic/secondOrderEqn/isotropic/AcousticElasticWaveEquationSEMKernel.hpp
     dg/acoustic/shared/AcousticFieldsDG.hpp
     dg/acoustic/secondOrderEqn/isotropic/AcousticWaveEquationDG.hpp 
     dg/acoustic/secondOrderEqn/isotropic/AcousticWaveEquationDGKernel.hpp)

# Specify solver sources
set( wavePropagationSolvers_sources
     shared/WaveSolverBase.cpp
     sem/acoustic/secondOrderEqn/isotropic/AcousticWaveEquationSEM.cpp
     sem/acoustic/secondOrderEqn/anisotropic/AcousticSecondOrderAnisotropicWaveEquationSEM.cpp
     sem/acoustic/secondOrderEqn/anisotropic/AcousticSecondOrderVTIFletcher.cpp
     sem/acoustic/secondOrderEqn/anisotropic/AcousticSecondOrderTTIFletcher.cpp     
     sem/acoustic/secondOrderEqn/anisotropic/AcousticSecondOrderVTIZhang.cpp
     sem/acoustic/secondOrderEqn/anisotropic/AcousticSecondOrderTTIZhang.cpp
     sem/acoustic/secondOrderEqn/anisotropic/AcousticSecondOrderAnisotropicTTISEM.cpp
     sem/acoustic/secondOrderEqn/anisotropic/AcousticTTIFletcherWaveEquationSEM.cpp
     sem/acoustic/secondOrderEqn/anisotropic/AcousticTTIZhangWaveEquationSEM.cpp
     sem/acoustic/secondOrderEqn/anisotropic/AcousticVTIFletcherWaveEquationSEM.cpp
     sem/acoustic/secondOrderEqn/anisotropic/AcousticVTIZhangWaveEquationSEM.cpp     
     sem/elastic/secondOrderEqn/isotropic/ElasticWaveEquationSEM.cpp
     sem/elastic/firstOrderEqn/isotropic/ElasticFirstOrderWaveEquationSEM.cpp
     sem/acoustic/firstOrderEqn/isotropic/AcousticFirstOrderWaveEquationSEM.cpp
     sem/acoustoelastic/secondOrderEqn/isotropic/AcousticElasticWaveEquationSEM.cpp
     dg/acoustic/secondOrderEqn/isotropic/AcousticWaveEquationDG.cpp )

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
