/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file SolidMechanicsInitialization.cpp
 */

#include "SolidMechanicsInitialization.hpp"

#include "events/tasks/TasksManager.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsStatistics.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"
#include "physicsSolvers/solidMechanics/contact/SolidMechanicsLagrangeContact.hpp"
#include "physicsSolvers/solidMechanics/contact/SolidMechanicsEmbeddedFractures.hpp"
#include "physicsSolvers/solidMechanics/contact/SolidMechanicsAugmentedLagrangianContact.hpp"
#include "physicsSolvers/LogLevelsInfo.hpp"

namespace geos
{

using namespace dataRepository;

template< typename SOLID_SOLVER >
SolidMechanicsInitialization< SOLID_SOLVER >::SolidMechanicsInitialization( const string & name,
                                                                            Group * const parent ):
  TaskBase( name, parent ),
  m_solidSolverName(),
  m_solidMechanicsStatisticsName(),
  m_solidSolver( nullptr ),
  m_solidMechanicsStatistics( nullptr ),
  m_solidMechanicsStateResetTask( name, parent )
{
  registerWrapper( viewKeyStruct::solidSolverNameString(), &m_solidSolverName ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of the solid mechanics solver" );

  registerWrapper( viewKeyStruct::solidMechanicsStatisticsNameString(), &m_solidMechanicsStatisticsName ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( "" ).
    setDescription( "Name of the solid mechanics statistics task" );

  registerWrapper( viewKeyStruct::resetDisplacementsString(), &m_resetDisplacements ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 1 ).
    setDescription( "Flag to reset displacements (1: reset, 0: do not reset)" );

  addLogLevel< logInfo::SolverInitialization >();
}

template< typename SOLID_SOLVER >
SolidMechanicsInitialization< SOLID_SOLVER >::~SolidMechanicsInitialization() = default;

template< typename SOLID_SOLVER >
void SolidMechanicsInitialization< SOLID_SOLVER >::postInputInitialization()
{
  ProblemManagerBase & problemManager = getProblemManagerBase( *this );
  PhysicsSolverManager & physicsSolverManager = problemManager.getPhysicsSolverManager();

  GEOS_THROW_IF( !physicsSolverManager.hasGroup( m_solidSolverName ),
                 GEOS_FMT( "{}: {} solver named {} not found",
                           getWrapperDataContext( viewKeyStruct::solidSolverNameString() ),
                           SOLID_SOLVER::catalogName(),
                           m_solidSolverName ),
                 InputError, getWrapperDataContext( viewKeyStruct::solidSolverNameString() ) );

  m_solidSolver = &physicsSolverManager.getGroup< SOLID_SOLVER >( m_solidSolverName );

  if( !m_solidMechanicsStatisticsName.empty() )
  {
    TasksManager & tasksManager = problemManager.getTasksManager();

    GEOS_THROW_IF( !tasksManager.hasGroup( m_solidMechanicsStatisticsName ),
                   GEOS_FMT( "{}: {} task named {} not found",
                             getWrapperDataContext( viewKeyStruct::solidMechanicsStatisticsNameString() ),
                             SolidMechanicsStatistics::catalogName(),
                             m_solidMechanicsStatisticsName ),
                   InputError, getWrapperDataContext( viewKeyStruct::solidMechanicsStatisticsNameString() ) );

    m_solidMechanicsStatistics = &tasksManager.getGroup< SolidMechanicsStatistics >( m_solidMechanicsStatisticsName );
  }

  m_solidMechanicsStateResetTask.setLogLevel( getLogLevel() );
  m_solidMechanicsStateResetTask.m_solidSolverName = m_solidSolverName;
  m_solidMechanicsStateResetTask.postInputInitialization();
}

template< typename SOLID_SOLVER >
bool SolidMechanicsInitialization< SOLID_SOLVER >::execute( real64 const time_n,
                                                            real64 const dt,
                                                            integer const cycleNumber,
                                                            integer const eventCounter,
                                                            real64 const eventProgress,
                                                            DomainPartition & domain )
{
  GEOS_LOG_LEVEL_RANK_0( logInfo::SolverInitialization,
                         GEOS_FMT( "Task `{}`: at time {}s, solid solver `{}` is set to perform stress initialization during the next time step(s)",
                                   getName(), time_n, m_solidSolverName ) );

  m_solidSolver->setStressInitialization( true );
  if( m_resetDisplacements )
  {
    m_solidMechanicsStateResetTask.execute( time_n, dt, cycleNumber, eventCounter, eventProgress, domain );
  }
  m_solidSolver->execute( time_n, dt, cycleNumber, eventCounter, eventProgress, domain );

  GEOS_LOG_LEVEL_RANK_0( logInfo::SolverInitialization,
                         GEOS_FMT( "Task `{}`: at time {}s, solid solver `{}` has completed stress initialization",
                                   getName(), time_n + dt, m_solidSolverName ) );
  m_solidSolver->setStressInitialization( false );

  if( m_solidMechanicsStatistics != nullptr )
  {
    m_solidMechanicsStatistics->execute( time_n, dt, cycleNumber, eventCounter, eventProgress, domain );
  }

  if( m_resetDisplacements )
  {
    m_solidMechanicsStateResetTask.execute( time_n, dt, cycleNumber, eventCounter, eventProgress, domain );
  }

  return false;
}

namespace
{
typedef SolidMechanicsInitialization< SolidMechanicsLagrangianFEM > SolidMechanicsLagrangianFEMInitialization;
typedef SolidMechanicsInitialization< SolidMechanicsLagrangeContact > SolidMechanicsLagrangeContactInitialization;
typedef SolidMechanicsInitialization< SolidMechanicsAugmentedLagrangianContact > SolidMechanicsAugmentedLagrangianContactInitialization;
typedef SolidMechanicsInitialization< SolidMechanicsEmbeddedFractures > SolidMechanicsEmbeddedFracturesInitialization;
REGISTER_CATALOG_ENTRY( TaskBase, SolidMechanicsLagrangianFEMInitialization, string const &, Group * const )
REGISTER_CATALOG_ENTRY( TaskBase, SolidMechanicsLagrangeContactInitialization, string const &, Group * const )
REGISTER_CATALOG_ENTRY( TaskBase, SolidMechanicsAugmentedLagrangianContactInitialization, string const &, Group * const )
REGISTER_CATALOG_ENTRY( TaskBase, SolidMechanicsEmbeddedFracturesInitialization, string const &, Group * const )
}

} /* namespace geos */
