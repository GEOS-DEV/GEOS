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
 * @file ThermalSinglePhasePoromechanicsSolver.cpp
 *
 */

#include "ThermalSinglePhasePoromechanicsSolver.hpp"

#include "constitutive/ConstitutiveManager.hpp"
#include "discretizationMethods/NumericalMethodsManager.hpp"
#include "finiteElement/Kinematics.h"
#include "mesh/DomainPartition.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "mesh/utilities/ComputationalGeometry.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBase.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;

ThermalSinglePhasePoromechanicsSolver::ThermalSinglePhasePoromechanicsSolver( const string & name,
                                                                              Group * const parent ):
  SolverBase( name, parent ),
  m_solidSolverName(),
  m_fluidSolverName(),
  m_couplingTypeOption( CouplingTypeOption::FixedStress )

{
  registerWrapper( viewKeyStruct::solidSolverNameString(), &m_solidSolverName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription(
    "Name of the solid mechanics solver to use in the ThermalSinglePhasePoromechanics solver" );

  registerWrapper( viewKeyStruct::fluidSolverNameString(), &m_fluidSolverName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription(
    "Name of the fluid solver to use in the ThermalSinglePhasePoromechanics solver" );

  registerWrapper( viewKeyStruct::couplingTypeOptionString(), &m_couplingTypeOption ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Coupling option. Valid options:\n* " + EnumStrings< CouplingTypeOption >::concat( "\n* " ) );

  registerWrapper( viewKeyStruct::subcyclingOptionString(), &m_subcyclingOption ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "turn on subcycling on each load step" );

}

void ThermalSinglePhasePoromechanicsSolver::registerDataOnMesh( Group & GEOSX_UNUSED_PARAM( meshBodies ))
{}

void ThermalSinglePhasePoromechanicsSolver::implicitStepComplete( real64 const & GEOSX_UNUSED_PARAM( time_n ),
                                                                  real64 const & GEOSX_UNUSED_PARAM( dt ),
                                                                  DomainPartition & GEOSX_UNUSED_PARAM( domain ) )
{}

void ThermalSinglePhasePoromechanicsSolver::postProcessInput()
{
  if( m_couplingTypeOption == CouplingTypeOption::FixedStress )
  {
    // For this coupled solver the minimum number of Newton Iter should be 0 for both flow and solid solver otherwise it
    // will never converge.
    SolidMechanicsLagrangianFEM &
    solidSolver = this->getParent().getGroup< SolidMechanicsLagrangianFEM >( m_solidSolverName );
    integer & minNewtonIterSolid = solidSolver.getNonlinearSolverParameters().m_minIterNewton;

    SinglePhaseBase &
    fluidSolver = this->getParent().getGroup< SinglePhaseBase >( m_fluidSolverName );
    integer & minNewtonIterFluid = fluidSolver.getNonlinearSolverParameters().m_minIterNewton;

    minNewtonIterSolid = 0;
    minNewtonIterFluid = 0;
  }
}

void ThermalSinglePhasePoromechanicsSolver::initializePostInitialConditionsPreSubGroups()
{}

ThermalSinglePhasePoromechanicsSolver::~ThermalSinglePhasePoromechanicsSolver()
{
  // TODO Auto-generated destructor stub
}

void ThermalSinglePhasePoromechanicsSolver::resetStateToBeginningOfStep( DomainPartition & GEOSX_UNUSED_PARAM( domain ) )
{}

real64 ThermalSinglePhasePoromechanicsSolver::solverStep( real64 const & time_n,
                                                          real64 const & dt,
                                                          int const cycleNumber,
                                                          DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;
  real64 dtReturn = dt;
  if( m_couplingTypeOption == CouplingTypeOption::FixedStress )
  {
    dtReturn = splitOperatorStep( time_n, dt, cycleNumber, domain );
  }
  else if( m_couplingTypeOption == CouplingTypeOption::TightlyCoupled )
  {
    GEOSX_ERROR( "CouplingTypeOption::FullyImplicit not yet implemented" );
  }
  return dtReturn;
}

real64 ThermalSinglePhasePoromechanicsSolver::splitOperatorStep( real64 const & time_n,
                                                                 real64 const & dt,
                                                                 integer const cycleNumber,
                                                                 DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;
  real64 dtReturn = dt;
  real64 dtReturnTemporary;

  SolidMechanicsLagrangianFEM &
  solidSolver = this->getParent().getGroup< SolidMechanicsLagrangianFEM >( m_solidSolverName );

  SinglePhaseBase &
  fluidSolver = this->getParent().getGroup< SinglePhaseBase >( m_fluidSolverName );

  fluidSolver.setupSystem( domain,
                           fluidSolver.getDofManager(),
                           fluidSolver.getLocalMatrix(),
                           fluidSolver.getSystemRhs(),
                           fluidSolver.getSystemSolution(),
                           true );

  solidSolver.setupSystem( domain,
                           solidSolver.getDofManager(),
                           solidSolver.getLocalMatrix(),
                           solidSolver.getSystemRhs(),
                           solidSolver.getSystemSolution() );

  fluidSolver.implicitStepSetup( time_n, dt, domain );

  solidSolver.implicitStepSetup( time_n, dt, domain );

  // this->implicitStepSetup( time_n, dt, domain );

  NonlinearSolverParameters & solverParams = getNonlinearSolverParameters();
  integer & iter = solverParams.m_numNewtonIterations;
  iter = 0;
  bool isConverged = false;
  while( iter < solverParams.m_maxIterNewton )
  {
    if( iter == 0 )
    {
      // reset the states of all slave solvers if any of them has been reset
      fluidSolver.resetStateToBeginningOfStep( domain );
      solidSolver.resetStateToBeginningOfStep( domain );
      resetStateToBeginningOfStep( domain );
    }

    GEOSX_LOG_LEVEL_RANK_0( 1, "\tIteration: " << iter+1 << ", MechanicsSolver: " );

    dtReturnTemporary = solidSolver.nonlinearImplicitStep( time_n,
                                                           dtReturn,
                                                           cycleNumber,
                                                           domain );

    if( dtReturnTemporary < dtReturn )
    {
      iter = 0;
      dtReturn = dtReturnTemporary;
      continue;
    }

    if( solidSolver.getNonlinearSolverParameters().m_numNewtonIterations == 0 && iter > 0 )
    {
      GEOSX_LOG_LEVEL_RANK_0( 1, "***** The iterative coupling has converged in " << iter << " iterations! *****\n" );
      isConverged = true;
      break;
    }
    else if( m_subcyclingOption == 0 && iter > 0 )
    {
      GEOSX_LOG_LEVEL_RANK_0( 1, "***** Single Pass solver, no subcycling *****\n" );
      isConverged = true;
      break;
    }

    GEOSX_LOG_LEVEL_RANK_0( 1, "\tIteration: " << iter+1 << ", FlowSolver: " );

    dtReturnTemporary = fluidSolver.nonlinearImplicitStep( time_n,
                                                           dtReturn,
                                                           cycleNumber,
                                                           domain );

    if( dtReturnTemporary < dtReturn )
    {
      iter = 0;
      dtReturn = dtReturnTemporary;
      continue;
    }

    if( m_subcyclingOption == 0 )
    {
      GEOSX_LOG_LEVEL_RANK_0( 1, "***** Single Pass solver, no subcycling *****\n" );
      isConverged = true;
      break;
    }

    ++iter;
  }

  GEOSX_ERROR_IF( !isConverged, "ThermalSinglePhasePoromechanicsSolver::SplitOperatorStep() did not converge" );

  fluidSolver.implicitStepComplete( time_n, dt, domain );
  solidSolver.implicitStepComplete( time_n, dt, domain );
  this->implicitStepComplete( time_n, dt, domain );

  return dtReturn;
}

REGISTER_CATALOG_ENTRY( SolverBase, ThermalSinglePhasePoromechanicsSolver, string const &, Group * const )

} /* namespace geosx */
