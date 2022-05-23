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
 * @file MultiphasePoromechanicsStateReset.cpp
 */

#include "MultiphasePoromechanicsStateReset.hpp"

#include "physicsSolvers/PhysicsSolverManager.hpp"
#include "physicsSolvers/multiphysics/MultiphasePoromechanicsSolver.hpp"
#include "mainInterface/ProblemManager.hpp"

namespace geosx
{

using namespace constitutive;
using namespace dataRepository;

MultiphasePoromechanicsStateReset::MultiphasePoromechanicsStateReset( const string & name,
                                                                      Group * const parent ):
  TaskBase( name, parent ),
  m_multiphasePoromechanicsSolverName()
{
  enableLogLevelInput();

  registerWrapper( viewKeyStruct::multiphasePoromechanicsSolverNameString(), &m_multiphasePoromechanicsSolverName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of the multiphase poromechanics solver" );

  registerWrapper( viewKeyStruct::useInitializationSolverConfigurationString(), &m_useInitializationSolverConfiguration ).
    setApplyDefaultValue( true ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Flag to reset the multiphase poromechanics solver configuration" );
}

MultiphasePoromechanicsStateReset::~MultiphasePoromechanicsStateReset()
{}

void MultiphasePoromechanicsStateReset::postProcessInput()
{
  ProblemManager & problemManager = this->getGroupByPath< ProblemManager >( "/Problem" );
  PhysicsSolverManager & physicsSolverManager = problemManager.getPhysicsSolverManager();

  GEOSX_THROW_IF( !physicsSolverManager.hasGroup( m_multiphasePoromechanicsSolverName ),
                  GEOSX_FMT( "Task {}: physics solver named {} not found",
                             getName(), m_multiphasePoromechanicsSolverName ),
                  InputError );

  m_multiphasePoromechanicsSolver =
    &physicsSolverManager.getGroup< MultiphasePoromechanicsSolver >( m_multiphasePoromechanicsSolverName );
}

bool MultiphasePoromechanicsStateReset::execute( real64 const time_n,
                                                 real64 const GEOSX_UNUSED_PARAM( dt ),
                                                 integer const GEOSX_UNUSED_PARAM( cycleNumber ),
                                                 integer const GEOSX_UNUSED_PARAM( eventCounter ),
                                                 real64 const GEOSX_UNUSED_PARAM( eventProgress ),
                                                 DomainPartition & GEOSX_UNUSED_PARAM( domain ) )
{
  if( m_useInitializationSolverConfiguration )
  {
    GEOSX_LOG_LEVEL_RANK_0( 1, GEOSX_FMT( "Task `{}`: at time {}s, physics solver `{}` is set to use a solver configuration designed specifically for initialization",
                                          getName(), time_n, m_multiphasePoromechanicsSolverName ) );
    m_multiphasePoromechanicsSolver->useInitializationSolverConfiguration( true );
  }
  else
  {
    GEOSX_LOG_LEVEL_RANK_0( 1, GEOSX_FMT( "Task `{}`: at time {}s, physics solver `{}` is set to use a standard solver configuration",
                                          getName(), time_n, m_multiphasePoromechanicsSolverName ) );
    m_multiphasePoromechanicsSolver->useInitializationSolverConfiguration( false );
  }

  return false;
}

REGISTER_CATALOG_ENTRY( TaskBase,
                        MultiphasePoromechanicsStateReset,
                        string const &, dataRepository::Group * const )

} /* namespace geosx */
