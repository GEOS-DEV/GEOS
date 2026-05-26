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
Package: events

Contains GEOS event types.
Manages the events and tasks.
#]]

#
# Specify all headers
#
set( events_headers
     EventBase.hpp
     EventManager.hpp
     HaltEvent.hpp
     LogLevelsInfo.hpp
     PeriodicEvent.hpp
     SoloEvent.hpp
     tasks/TaskBase.hpp
     tasks/TasksManager.hpp
     mpmEvents/MPMEventManager.hpp
     mpmEvents/MPMEventBase.hpp
     mpmEvents/MaterialSwapMPMEvent.hpp
     mpmEvents/AnnealMPMEvent.hpp
     mpmEvents/InitializeStressMPMEvent.hpp
     mpmEvents/HealMPMEvent.hpp
     mpmEvents/PolymerHealMPMEvent.hpp
     mpmEvents/CrystalHealMPMEvent.hpp
     mpmEvents/InsertPeriodicContactSurfacesMPMEvent.hpp
     mpmEvents/MachineSampleMPMEvent.hpp
     mpmEvents/FrictionCoefficientSwapMPMEvent.hpp
     mpmEvents/BodyForceUpdateMPMEvent.hpp
     mpmEvents/DeformationUpdateMPMEvent.hpp
     mpmEvents/BoreholePressureMPMEvent.hpp
     mpmEvents/ConfiningPressureMPMEvent.hpp
     mpmEvents/ReferenceCohesiveZonesMPMEvent.hpp
     mpmEvents/ResetDeformationGradientMPMEvent.hpp
     mpmEvents/TemperatureProfileMPMEvent.hpp
     mpmEvents/TransformParticlesMPMEvent.hpp
     mpmEvents/UpdateSurfacesMPMEvent.hpp
   )

#
# Specify all sources
#
set( events_sources
     EventBase.cpp
     EventManager.cpp
     HaltEvent.cpp
     PeriodicEvent.cpp
     SoloEvent.cpp
     tasks/TaskBase.cpp
     tasks/TasksManager.cpp
     mpmEvents/MPMEventManager.cpp
     mpmEvents/MPMEventBase.cpp
     mpmEvents/MaterialSwapMPMEvent.cpp
     mpmEvents/AnnealMPMEvent.cpp
     mpmEvents/InitializeStressMPMEvent.cpp
     mpmEvents/HealMPMEvent.cpp
     mpmEvents/PolymerHealMPMEvent.cpp
     mpmEvents/CrystalHealMPMEvent.cpp
     mpmEvents/InsertPeriodicContactSurfacesMPMEvent.cpp
     mpmEvents/MachineSampleMPMEvent.cpp
     mpmEvents/FrictionCoefficientSwapMPMEvent.cpp
     mpmEvents/BodyForceUpdateMPMEvent.cpp
     mpmEvents/DeformationUpdateMPMEvent.cpp
     mpmEvents/BoreholePressureMPMEvent.cpp
     mpmEvents/ConfiningPressureMPMEvent.cpp
     mpmEvents/ReferenceCohesiveZonesMPMEvent.cpp
     mpmEvents/BoreholePressureMPMEvent.cpp
     mpmEvents/ResetDeformationGradientMPMEvent.cpp
     mpmEvents/TemperatureProfileMPMEvent.cpp
     mpmEvents/TransformParticlesMPMEvent.cpp
     mpmEvents/UpdateSurfacesMPMEvent.cpp
   )

set( dependencyList ${parallelDeps} functions )

geos_decorate_link_dependencies( LIST decoratedDependencies
                                 DEPENDENCIES ${dependencyList} )

blt_add_library( NAME       events
                 SOURCES    ${events_sources}
                 HEADERS    ${events_headers}
                 DEPENDS_ON ${decoratedDependencies}
                 OBJECT     ${GEOS_BUILD_OBJ_LIBS}
                 SHARED     ${GEOS_BUILD_SHARED_LIBS}
               )

target_include_directories( events PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )

install( TARGETS events LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib )
