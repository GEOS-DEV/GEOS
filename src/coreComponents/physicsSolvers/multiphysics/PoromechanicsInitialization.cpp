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
#include "mainInterface/ProblemManager.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBase.hpp"
#include "physicsSolvers/multiphysics/SinglePhaseReservoirAndWells.hpp"
#include "physicsSolvers/multiphysics/CompositionalMultiphaseReservoirAndWells.hpp"

namespace geos
{

using namespace dataRepository;

namespace
{

// This is meant to be specialized to work, see below
template< typename POROMECHANICS_SOLVER > class
  PoromechanicsCatalogNames {};

// Class specializations for a POROMECHANICS_SOLVER set to SinglePhasePoromechanics
template<> class PoromechanicsCatalogNames< SinglePhasePoromechanics< SinglePhaseBase > >
{
public:
  static string name() { return SinglePhasePoromechanics< SinglePhaseBase >::catalogName() + "Initialization"; }
};
template<> class PoromechanicsCatalogNames< SinglePhasePoromechanics< SinglePhaseReservoirAndWells< SinglePhaseBase > > >
{
public:
  static string name() { return SinglePhasePoromechanics< SinglePhaseReservoirAndWells< SinglePhaseBase > >::catalogName() + "Initialization"; }
};
// Class specializations for a POROMECHANICS_SOLVER set to MultiphasePoromechanics
template<> class PoromechanicsCatalogNames< MultiphasePoromechanics< CompositionalMultiphaseBase > >
{
public:
  static string name() { return MultiphasePoromechanics< CompositionalMultiphaseBase >::catalogName() + "Initialization"; }
};
template<> class PoromechanicsCatalogNames< MultiphasePoromechanics< CompositionalMultiphaseReservoirAndWells< CompositionalMultiphaseBase > > >
{
public:
  static string name() { return MultiphasePoromechanics< CompositionalMultiphaseReservoirAndWells< CompositionalMultiphaseBase > >::catalogName() + "Initialization"; }
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
  m_poromechanicsSolverName(),
  m_solidMechanicsStateResetTask( name, parent )
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

  GEOS_THROW_IF( !physicsSolverManager.hasGroup( m_poromechanicsSolverName ),
                 GEOS_FMT( "{}: physics solver named {} not found",
                           getWrapperDataContext( viewKeyStruct::poromechanicsSolverNameString() ),
                           m_poromechanicsSolverName ),
                 InputError );

  m_poromechanicsSolver = &physicsSolverManager.getGroup< POROMECHANICS_SOLVER >( m_poromechanicsSolverName );

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
  if( m_performStressInitialization )
  {
    GEOS_LOG_LEVEL_RANK_0( 1, GEOS_FMT( "Task `{}`: at time {}s, physics solver `{}` is set to perform stress initialization during the next time step(s)",
                                        getName(), time_n, m_poromechanicsSolverName ) );
    m_poromechanicsSolver->setStressInitialization( true );
  }
  else
  {
    GEOS_LOG_LEVEL_RANK_0( 1, GEOS_FMT( "Task `{}`: at time {}s, physics solver `{}` has completed stress initialization",
                                        getName(), time_n, m_poromechanicsSolverName ) );
    m_poromechanicsSolver->setStressInitialization( false );
  }

  m_solidMechanicsStateResetTask.execute( time_n, dt, cycleNumber, eventCounter, eventProgress, domain );

  // always returns false because we don't want early return (see EventManager.cpp)
  return false;
}

namespace
{
typedef PoromechanicsInitialization< MultiphasePoromechanics< CompositionalMultiphaseBase > > MultiphasePoromechanicsInitialization;
typedef PoromechanicsInitialization< MultiphasePoromechanics< CompositionalMultiphaseReservoirAndWells< CompositionalMultiphaseBase > > > MultiphaseReservoirPoromechanicsInitialization;
typedef PoromechanicsInitialization< SinglePhasePoromechanics< SinglePhaseBase > > SinglePhasePoromechanicsInitialization;
typedef PoromechanicsInitialization< SinglePhasePoromechanics< SinglePhaseReservoirAndWells< SinglePhaseBase > > > SinglePhaseReservoirPoromechanicsInitialization;
REGISTER_CATALOG_ENTRY( TaskBase, MultiphasePoromechanicsInitialization, string const &, Group * const )
REGISTER_CATALOG_ENTRY( TaskBase, MultiphaseReservoirPoromechanicsInitialization, string const &, Group * const )
REGISTER_CATALOG_ENTRY( TaskBase, SinglePhasePoromechanicsInitialization, string const &, Group * const )
REGISTER_CATALOG_ENTRY( TaskBase, SinglePhaseReservoirPoromechanicsInitialization, string const &, Group * const )
}

} /* namespace geos */
