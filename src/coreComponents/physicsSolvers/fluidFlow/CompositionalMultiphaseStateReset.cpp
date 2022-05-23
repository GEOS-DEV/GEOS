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
 * @file CompositionalMultiphaseStateReset.cpp
 */

#include "CompositionalMultiphaseStateReset.hpp"

#include "physicsSolvers/PhysicsSolverManager.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBase.hpp"
#include "mainInterface/ProblemManager.hpp"

namespace geosx
{

using namespace constitutive;
using namespace dataRepository;

CompositionalMultiphaseStateReset::CompositionalMultiphaseStateReset( const string & name,
                                                                      Group * const parent ):
  TaskBase( name, parent ),
  m_flowSolverName()
{
  enableLogLevelInput();

  registerWrapper( viewKeyStruct::flowSolverNameString(), &m_flowSolverName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of the flow solver" );

  registerWrapper( viewKeyStruct::freezeFlowVariablesDuringStepString(), &m_freezeFlowVariablesDuringStep ).
    setApplyDefaultValue( true ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Flag to freeze the flow variables during a time step" );
}

CompositionalMultiphaseStateReset::~CompositionalMultiphaseStateReset()
{}

void CompositionalMultiphaseStateReset::postProcessInput()
{
  ProblemManager & problemManager = this->getGroupByPath< ProblemManager >( "/Problem" );
  PhysicsSolverManager & physicsSolverManager = problemManager.getPhysicsSolverManager();

  GEOSX_THROW_IF( !physicsSolverManager.hasGroup( m_flowSolverName ),
                  GEOSX_FMT( "Task {}: physics solver named {} not found",
                             getName(), m_flowSolverName ),
                  InputError );

  m_flowSolver = &physicsSolverManager.getGroup< CompositionalMultiphaseBase >( m_flowSolverName );
}

bool CompositionalMultiphaseStateReset::execute( real64 const time_n,
                                                 real64 const GEOSX_UNUSED_PARAM( dt ),
                                                 integer const GEOSX_UNUSED_PARAM( cycleNumber ),
                                                 integer const GEOSX_UNUSED_PARAM( eventCounter ),
                                                 real64 const GEOSX_UNUSED_PARAM( eventProgress ),
                                                 DomainPartition & GEOSX_UNUSED_PARAM( domain ) )
{
  if( m_freezeFlowVariablesDuringStep )
  {
    GEOSX_LOG_LEVEL_RANK_0( 1, GEOSX_FMT( "Task `{}`: at time {}s, physics solver `{}` is set to freeze the flow variables during the next time step(s)",
                                          getName(), time_n, m_flowSolverName ) );
    m_flowSolver->freezeFlowVariablesDuringStep( true );
  }
  else
  {
    GEOSX_LOG_LEVEL_RANK_0( 1, GEOSX_FMT( "Task `{}`: at time {}s, physics solver `{}` is set to unfreeze the flow variables during the next time step(s)",
                                          getName(), time_n, m_flowSolverName ) );
    m_flowSolver->freezeFlowVariablesDuringStep( false );
  }

  return false;
}

REGISTER_CATALOG_ENTRY( TaskBase,
                        CompositionalMultiphaseStateReset,
                        string const &, dataRepository::Group * const )

} /* namespace geosx */
