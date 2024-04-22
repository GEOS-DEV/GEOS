/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file PoromechanicsInitialization.cpp
 */

#include "PoromechanicsInitialization.hpp"

#include "physicsSolvers/PhysicsSolverManager.hpp"
#include "physicsSolvers/multiphysics/MultiphasePoromechanics.hpp"
#include "physicsSolvers/multiphysics/SinglePhasePoromechanics.hpp"
#include "physicsSolvers/multiphysics/SinglePhasePoromechanicsConformingFractures.hpp"
#include "physicsSolvers/multiphysics/SinglePhasePoromechanicsEmbeddedFractures.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBase.hpp"
#include "physicsSolvers/multiphysics/SinglePhaseReservoirAndWells.hpp"
#include "physicsSolvers/multiphysics/CompositionalMultiphaseReservoirAndWells.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsStatistics.hpp"
#include "events/tasks/TasksManager.hpp"

namespace geos
{

using namespace dataRepository;

template< typename POROMECHANICS_SOLVER >
PoromechanicsInitialization< POROMECHANICS_SOLVER >::
PoromechanicsInitialization( const string & name,
                             Group * const parent ):
  TaskBase( name, parent ),
  m_poromechanicsSolverName(),
  m_solidMechanicsStatistics(),
  m_solidMechanicsStateResetTask( name, parent )
{
  enableLogLevelInput();

  registerWrapper( viewKeyStruct::poromechanicsSolverNameString(), &m_poromechanicsSolverName ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of the poromechanics solver" );

  registerWrapper( viewKeyStruct::solidMechanicsStatisticsNameString(), &m_solidMechanicsStatisticsName ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( "" ).
    setDescription( "Name of the solid mechanics statistics" );
}

template< typename POROMECHANICS_SOLVER >
PoromechanicsInitialization< POROMECHANICS_SOLVER >::~PoromechanicsInitialization() {}

template< typename POROMECHANICS_SOLVER >
void
PoromechanicsInitialization< POROMECHANICS_SOLVER >::
postProcessInput()
{
  ProblemManager & problemManager = this->getGroupByPath< ProblemManager >( "/Problem" );
  PhysicsSolverManager & physicsSolverManager = problemManager.getPhysicsSolverManager();

  GEOS_THROW_IF( !physicsSolverManager.hasGroup( m_poromechanicsSolverName ),
                 GEOS_FMT( "{}: physics solver named {} not found",
                           getWrapperDataContext( viewKeyStruct::poromechanicsSolverNameString() ),
                           m_poromechanicsSolverName ),
                 InputError );

  m_poromechanicsSolver = &physicsSolverManager.getGroup< POROMECHANICS_SOLVER >( m_poromechanicsSolverName );

  if( !m_solidMechanicsStatisticsName.empty())
  {
    TasksManager & tasksManager = problemManager.getTasksManager();

    GEOS_THROW_IF( !tasksManager.hasGroup( m_solidMechanicsStatisticsName ),
                   GEOS_FMT( "{}: statistics task named {} not found",
                             getWrapperDataContext( viewKeyStruct::solidMechanicsStatisticsNameString() ),
                             m_solidMechanicsStatisticsName ),
                   InputError );

    m_solidMechanicsStatistics = &tasksManager.getGroup< SolidMechanicsStatistics >( m_solidMechanicsStatisticsName );
  }

  m_solidMechanicsStateResetTask.setLogLevel( getLogLevel());
  m_solidMechanicsStateResetTask.m_solidSolverName = m_poromechanicsSolver->solidMechanicsSolver()->getName();
  m_solidMechanicsStateResetTask.postProcessInput();
}

template< typename POROMECHANICS_SOLVER >
bool
PoromechanicsInitialization< POROMECHANICS_SOLVER >::
execute( real64 const time_n,
         real64 const dt,
         integer const cycleNumber,
         integer const eventCounter,
         real64 const eventProgress,
         DomainPartition & domain )
{
  GEOS_LOG_LEVEL_RANK_0( 1, GEOS_FMT( "Task `{}`: at time {}s, physics solver `{}` is set to perform stress initialization during the next time step(s)",
                                      getName(), time_n, m_poromechanicsSolverName ) );
  m_poromechanicsSolver->setStressInitialization( true );

  m_solidMechanicsStateResetTask.execute( time_n, dt, cycleNumber, eventCounter, eventProgress, domain );

  m_poromechanicsSolver->execute( time_n, dt, cycleNumber, eventCounter, eventProgress, domain );

  GEOS_LOG_LEVEL_RANK_0( 1, GEOS_FMT( "Task `{}`: at time {}s, physics solver `{}` has completed stress initialization",
                                      getName(), time_n + dt, m_poromechanicsSolverName ) );
  m_poromechanicsSolver->setStressInitialization( false );

  if( m_solidMechanicsStatistics != nullptr )
  {
    m_solidMechanicsStatistics->execute( time_n, dt, cycleNumber, eventCounter, eventProgress, domain );
  }

  m_solidMechanicsStateResetTask.execute( time_n, dt, cycleNumber, eventCounter, eventProgress, domain );

  // always returns false because we don't want early return (see EventManager.cpp)
  return false;
}

namespace
{
typedef PoromechanicsInitialization< MultiphasePoromechanics<> > MultiphasePoromechanicsInitialization;
typedef PoromechanicsInitialization< MultiphasePoromechanics< CompositionalMultiphaseReservoirAndWells<> > > MultiphaseReservoirPoromechanicsInitialization;
typedef PoromechanicsInitialization< SinglePhasePoromechanics<> > SinglePhasePoromechanicsInitialization;
typedef PoromechanicsInitialization< SinglePhasePoromechanicsConformingFractures<> > SinglePhasePoromechanicsConformingFracturesInitialization;
typedef PoromechanicsInitialization< SinglePhasePoromechanicsEmbeddedFractures > SinglePhasePoromechanicsEmbeddedFracturesInitialization;
typedef PoromechanicsInitialization< SinglePhasePoromechanics< SinglePhaseReservoirAndWells<> > > SinglePhaseReservoirPoromechanicsInitialization;
REGISTER_CATALOG_ENTRY( TaskBase, MultiphasePoromechanicsInitialization, string const &, Group * const )
REGISTER_CATALOG_ENTRY( TaskBase, MultiphaseReservoirPoromechanicsInitialization, string const &, Group * const )
REGISTER_CATALOG_ENTRY( TaskBase, SinglePhasePoromechanicsInitialization, string const &, Group * const )
REGISTER_CATALOG_ENTRY( TaskBase, SinglePhasePoromechanicsConformingFracturesInitialization, string const &, Group * const )
REGISTER_CATALOG_ENTRY( TaskBase, SinglePhasePoromechanicsEmbeddedFracturesInitialization, string const &, Group * const )
REGISTER_CATALOG_ENTRY( TaskBase, SinglePhaseReservoirPoromechanicsInitialization, string const &, Group * const )
}

} /* namespace geos */
