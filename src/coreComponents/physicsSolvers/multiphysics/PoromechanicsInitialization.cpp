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
#include "physicsSolvers/multiphysics/MultiphasePoromechanicsSolver.hpp"
#include "physicsSolvers/multiphysics/SinglePhasePoromechanicsSolver.hpp"
#include "mainInterface/ProblemManager.hpp"

namespace geosx
{

using namespace dataRepository;

namespace
{

// This is meant to be specialized to work, see below
template< typename POROMECHANICS_SOLVER > class
  PoromechanicsCatalogNames {};

// Class specialization for a POROMECHANICS_SOLVER set to SinglePhasePoromechanicsSolver
template<> class PoromechanicsCatalogNames< SinglePhasePoromechanicsSolver >
{
public:
  static string name() { return "SinglePhasePoromechanicsInitialization"; }
};
// Class specialization for a POROMECHANICS_SOLVER set to MultiphasePoromechanics
template<> class PoromechanicsCatalogNames< MultiphasePoromechanicsSolver >
{
public:
  static string name() { return "MultiphasePoromechanicsInitialization"; }
};

}

// provide a definition for catalogName()
template< typename POROMECHANICS_SOLVER >
string
PoromechanicsInitialization< POROMECHANICS_SOLVER >::
catalogName()
{
  return PoromechanicsCatalogNames< POROMECHANICS_SOLVER >::name();
}

template< typename POROMECHANICS_SOLVER >
PoromechanicsInitialization< POROMECHANICS_SOLVER >::
PoromechanicsInitialization( const string & name,
                             Group * const parent ):
  TaskBase( name, parent ),
  m_poromechanicsSolverName()
{
  enableLogLevelInput();

  registerWrapper( viewKeyStruct::poromechanicsSolverNameString(), &m_poromechanicsSolverName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of the poromechanics solver" );

  registerWrapper( viewKeyStruct::performStressInitializationString(), &m_performStressInitialization ).
    setApplyDefaultValue( true ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Flag to indicate that the solver is going to perform stress initialization" );
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

  GEOSX_THROW_IF( !physicsSolverManager.hasGroup( m_poromechanicsSolverName ),
                  GEOSX_FMT( "Task {}: physics solver named {} not found",
                             getName(), m_poromechanicsSolverName ),
                  InputError );

  m_poromechanicsSolver = &physicsSolverManager.getGroup< POROMECHANICS_SOLVER >( m_poromechanicsSolverName );
}

template< typename POROMECHANICS_SOLVER >
bool
PoromechanicsInitialization< POROMECHANICS_SOLVER >::
execute( real64 const time_n,
         real64 const GEOSX_UNUSED_PARAM( dt ),
         integer const GEOSX_UNUSED_PARAM( cycleNumber ),
         integer const GEOSX_UNUSED_PARAM( eventCounter ),
         real64 const GEOSX_UNUSED_PARAM( eventProgress ),
         DomainPartition & GEOSX_UNUSED_PARAM( domain ) )
{
  if( m_performStressInitialization )
  {
    GEOSX_LOG_LEVEL_RANK_0( 1, GEOSX_FMT( "Task `{}`: at time {}s, physics solver `{}` is set to perform stress initialization during the next time step(s)",
                                          getName(), time_n, m_poromechanicsSolverName ) );
    m_poromechanicsSolver->performStressInitialization( true );
  }
  else
  {
    GEOSX_LOG_LEVEL_RANK_0( 1, GEOSX_FMT( "Task `{}`: at time {}s, physics solver `{}` has completed stress initialization",
                                          getName(), time_n, m_poromechanicsSolverName ) );
    m_poromechanicsSolver->performStressInitialization( false );
  }

  return false;
}

namespace
{
typedef PoromechanicsInitialization< MultiphasePoromechanicsSolver > MultiphasePoromechanicsInitialization;
typedef PoromechanicsInitialization< SinglePhasePoromechanicsSolver > SinglePhasePoromechanicsInitialization;
REGISTER_CATALOG_ENTRY( TaskBase, MultiphasePoromechanicsInitialization, string const &, Group * const )
REGISTER_CATALOG_ENTRY( TaskBase, SinglePhasePoromechanicsInitialization, string const &, Group * const )
}

} /* namespace geosx */
