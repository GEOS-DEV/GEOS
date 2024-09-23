/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file FlowProppantTransportSolver.cpp
 */

#include "FlowProppantTransportSolver.hpp"

#include "mesh/DomainPartition.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBase.hpp"
#include "physicsSolvers/fluidFlow/proppantTransport/ProppantTransport.hpp"

namespace geos
{

using namespace dataRepository;
using namespace constitutive;

FlowProppantTransportSolver::FlowProppantTransportSolver( const string & name,
                                                          Group * const parent ):
  Base( name, parent )
{}

void FlowProppantTransportSolver::preStepUpdate( real64 const & time_n,
                                                 real64 const & dt,
                                                 DomainPartition & domain )
{
  if( time_n <= 0.0 )
  {
    // We need resize composition array in fractures after they are generated by SurfaceGenerator solver
    forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                                 MeshLevel & mesh,
                                                                 arrayView1d< string const > const & regionNames )
    {
      proppantTransportSolver()->resizeFractureFields( mesh, regionNames );
      // We need re-apply initial conditions to fractures after they are generated
      FieldSpecificationManager const & boundaryConditionManager = FieldSpecificationManager::getInstance();
      boundaryConditionManager.applyInitialConditions( mesh );
    } );
  }

  flowSolver()->setupSystem( domain,
                             flowSolver()->getDofManager(),
                             flowSolver()->getLocalMatrix(),
                             flowSolver()->getSystemRhs(),
                             flowSolver()->getSystemSolution() );


  flowSolver()->implicitStepSetup( time_n, dt, domain );

  proppantTransportSolver()->setupSystem( domain,
                                          proppantTransportSolver()->getDofManager(),
                                          proppantTransportSolver()->getLocalMatrix(),
                                          proppantTransportSolver()->getSystemRhs(),
                                          proppantTransportSolver()->getSystemSolution() );


  proppantTransportSolver()->implicitStepSetup( time_n, dt, domain );

  proppantTransportSolver()->preStepUpdate( time_n, dt, domain );
}

void FlowProppantTransportSolver::postStepUpdate( real64 const & time_n,
                                                  real64 const & dt,
                                                  DomainPartition & domain )
{
  flowSolver()->implicitStepComplete( time_n, dt, domain );
  proppantTransportSolver()->implicitStepComplete( time_n, dt, domain );
  proppantTransportSolver()->postStepUpdate( time_n, dt, domain );
}

real64 FlowProppantTransportSolver::sequentiallyCoupledSolverStep( real64 const & time_n,
                                                                   real64 const & dt,
                                                                   int const cycleNumber,
                                                                   DomainPartition & domain )
{
  real64 dtReturn = dt;
  real64 dtReturnTemporary;

  preStepUpdate( time_n, dt, domain );

  int iter = 0;
  while( iter < this->m_nonlinearSolverParameters.m_maxIterNewton )
  {
    if( iter == 0 )
    {
      // reset the states of all sub-solvers if any of them has been reset
      resetStateToBeginningOfStep( domain );
    }

    GEOS_LOG_LEVEL_RANK_0( 1, "  Iteration: " << iter+1  << ", FlowSolver: " );

    dtReturnTemporary = flowSolver()->nonlinearImplicitStep( time_n, dtReturn, cycleNumber, domain );

    if( dtReturnTemporary < dtReturn )
    {
      iter = 0;
      dtReturn = dtReturnTemporary;
      continue;
    }

    NonlinearSolverParameters const & fluidNonLinearParams = flowSolver()->getNonlinearSolverParameters();
    if( fluidNonLinearParams.m_numNewtonIterations <= this->m_nonlinearSolverParameters.m_minIterNewton && iter > 0 )
    {
      m_solverStatistics.logNonlinearIteration();
      GEOS_LOG_LEVEL_RANK_0( 1, GEOS_FMT( "***** The iterative coupling has converged in {} iterations *****", iter ) );
      break;
    }

    GEOS_LOG_LEVEL_RANK_0( 1, GEOS_FMT( "  Iteration: {}, Proppant Solver: ", iter+1 ) );

    dtReturnTemporary = proppantTransportSolver()->nonlinearImplicitStep( time_n, dtReturn, cycleNumber, domain );

    if( dtReturnTemporary < dtReturn )
    {
      iter = 0;
      dtReturn = dtReturnTemporary;
      continue;
    }

    ++iter;
  }

  postStepUpdate( time_n, dtReturn, domain );

  return dtReturn;
}

REGISTER_CATALOG_ENTRY( PhysicsSolverBase, FlowProppantTransportSolver, string const &, Group * const )

} /* namespace geos */
