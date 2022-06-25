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
 * @file SinglePhasePoromechanicsStateReset.cpp
 */

#include "SinglePhasePoromechanicsStateReset.hpp"

#include "physicsSolvers/PhysicsSolverManager.hpp"
#include "physicsSolvers/multiphysics/SinglePhasePoromechanicsSolver.hpp"
#include "mainInterface/ProblemManager.hpp"

namespace geosx
{

using namespace constitutive;
using namespace dataRepository;

SinglePhasePoromechanicsStateReset::SinglePhasePoromechanicsStateReset( const string & name,
                                                                        Group * const parent ):
  TaskBase( name, parent ),
  m_singlePhasePoromechanicsSolverName()
{
  enableLogLevelInput();

  registerWrapper( viewKeyStruct::singlePhasePoromechanicsSolverNameString(), &m_singlePhasePoromechanicsSolverName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of the singlePhase poromechanics solver" );

  registerWrapper( viewKeyStruct::useInitializationSolverConfigurationString(), &m_useInitializationSolverConfiguration ).
    setApplyDefaultValue( true ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Flag to reset the singlePhase poromechanics solver configuration" );
}

SinglePhasePoromechanicsStateReset::~SinglePhasePoromechanicsStateReset()
{}

void SinglePhasePoromechanicsStateReset::postProcessInput()
{
  ProblemManager & problemManager = this->getGroupByPath< ProblemManager >( "/Problem" );
  PhysicsSolverManager & physicsSolverManager = problemManager.getPhysicsSolverManager();

  GEOSX_THROW_IF( !physicsSolverManager.hasGroup( m_singlePhasePoromechanicsSolverName ),
                  GEOSX_FMT( "Task {}: physics solver named {} not found",
                             getName(), m_singlePhasePoromechanicsSolverName ),
                  InputError );

  m_singlePhasePoromechanicsSolver =
    &physicsSolverManager.getGroup< SinglePhasePoromechanicsSolver >( m_singlePhasePoromechanicsSolverName );
}

bool SinglePhasePoromechanicsStateReset::execute( real64 const time_n,
                                                  real64 const GEOSX_UNUSED_PARAM( dt ),
                                                  integer const GEOSX_UNUSED_PARAM( cycleNumber ),
                                                  integer const GEOSX_UNUSED_PARAM( eventCounter ),
                                                  real64 const GEOSX_UNUSED_PARAM( eventProgress ),
                                                  DomainPartition & GEOSX_UNUSED_PARAM( domain ) )
{
  if( m_useInitializationSolverConfiguration )
  {
    GEOSX_LOG_LEVEL_RANK_0( 1, GEOSX_FMT( "Task `{}`: at time {}s, physics solver `{}` is set to use a solver configuration designed specifically for initialization",
                                          getName(), time_n, m_singlePhasePoromechanicsSolverName ) );
    m_singlePhasePoromechanicsSolver->useInitializationSolverConfiguration( true );
  }
  else
  {
    GEOSX_LOG_LEVEL_RANK_0( 1, GEOSX_FMT( "Task `{}`: at time {}s, physics solver `{}` is set to use a standard solver configuration",
                                          getName(), time_n, m_singlePhasePoromechanicsSolverName ) );
    m_singlePhasePoromechanicsSolver->useInitializationSolverConfiguration( false );
  }

  return false;
}

REGISTER_CATALOG_ENTRY( TaskBase,
                        SinglePhasePoromechanicsStateReset,
                        string const &, dataRepository::Group * const )

} /* namespace geosx */
