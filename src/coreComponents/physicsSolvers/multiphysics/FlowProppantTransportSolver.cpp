/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file FlowProppantTransportSolver.cpp
 *
 */


#include "FlowProppantTransportSolver.hpp"

#include "constitutive/ConstitutiveManager.hpp"
#include "finiteElement/Kinematics.h"
#include "managers/NumericalMethodsManager.hpp"
#include "managers/DomainPartition.hpp"
#include "managers/FieldSpecification/FieldSpecificationManager.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseFVM.hpp"
#include "physicsSolvers/fluidFlow/ProppantTransport.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;

FlowProppantTransportSolver::FlowProppantTransportSolver( const std::string & name, Group * const parent ):
  SolverBase( name, parent ),
  m_proppantSolverName(),
  m_flowSolverName()
{
  registerWrapper( viewKeyStruct::proppantSolverNameString, &m_proppantSolverName )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Name of the proppant transport solver to use in the flowProppantTransport solver" );

  registerWrapper( viewKeyStruct::flowSolverNameString, &m_flowSolverName )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Name of the flow solver to use in the flowProppantTransport solver" );

}

void FlowProppantTransportSolver::RegisterDataOnMesh( dataRepository::Group * const )
{}

void FlowProppantTransportSolver::ImplicitStepSetup( real64 const & GEOSX_UNUSED_PARAM( time_n ),
                                                     real64 const & GEOSX_UNUSED_PARAM( dt ),
                                                     DomainPartition * const GEOSX_UNUSED_PARAM( domain ),
                                                     DofManager & GEOSX_UNUSED_PARAM( dofManager ),
                                                     ParallelMatrix & GEOSX_UNUSED_PARAM( matrix ),
                                                     ParallelVector & GEOSX_UNUSED_PARAM( rhs ),
                                                     ParallelVector & GEOSX_UNUSED_PARAM( solution ) )
{}

void FlowProppantTransportSolver::ImplicitStepComplete( real64 const & GEOSX_UNUSED_PARAM( time_n ),
                                                        real64 const & GEOSX_UNUSED_PARAM( dt ),
                                                        DomainPartition * const GEOSX_UNUSED_PARAM( domain ) )
{}

void FlowProppantTransportSolver::PostProcessInput()
{}

void FlowProppantTransportSolver::InitializePostInitialConditions_PreSubGroups( Group * const GEOSX_UNUSED_PARAM( problemManager ))
{}

FlowProppantTransportSolver::~FlowProppantTransportSolver()
{
  // TODO Auto-generated destructor stub
}

void FlowProppantTransportSolver::ResetStateToBeginningOfStep( DomainPartition * const GEOSX_UNUSED_PARAM( domain ) )
{}

real64 FlowProppantTransportSolver::SolverStep( real64 const & time_n,
                                                real64 const & dt,
                                                int const cycleNumber,
                                                DomainPartition * const domain )
{

  real64 dtReturn = dt;
  real64 dtReturnTemporary;

  ProppantTransport &
  proppantSolver = *(this->getParent()->GetGroup( m_proppantSolverName )->group_cast< ProppantTransport * >());

  SinglePhaseFVM< SinglePhaseProppantBase > &
  flowSolver = *(this->getParent()->GetGroup( m_flowSolverName )->group_cast< SinglePhaseFVM< SinglePhaseProppantBase > * >());

  proppantSolver.ResizeFractureFields( time_n, dt, domain );

  if( cycleNumber == 0 )
  {

    FieldSpecificationManager const & boundaryConditionManager = FieldSpecificationManager::get();

    boundaryConditionManager.ApplyInitialConditions( domain );

  }

  flowSolver.SetupSystem( domain,
                          flowSolver.getDofManager(),
                          flowSolver.getSystemMatrix(),
                          flowSolver.getSystemRhs(),
                          flowSolver.getSystemSolution() );


  flowSolver.ImplicitStepSetup( time_n, dt, domain,
                                flowSolver.getDofManager(),
                                flowSolver.getSystemMatrix(),
                                flowSolver.getSystemRhs(),
                                flowSolver.getSystemSolution() );

  proppantSolver.SetupSystem( domain,
                              proppantSolver.getDofManager(),
                              proppantSolver.getSystemMatrix(),
                              proppantSolver.getSystemRhs(),
                              proppantSolver.getSystemSolution() );


  proppantSolver.ImplicitStepSetup( time_n, dt, domain,
                                    proppantSolver.getDofManager(),
                                    proppantSolver.getSystemMatrix(),
                                    proppantSolver.getSystemRhs(),
                                    proppantSolver.getSystemSolution() );


  proppantSolver.PreStepUpdate( time_n, dt, cycleNumber, domain );

  this->ImplicitStepSetup( time_n, dt, domain, m_dofManager, m_matrix, m_rhs, m_solution );

  int iter = 0;
  while( iter <  this->m_nonlinearSolverParameters.m_maxIterNewton )
  {
    if( iter == 0 )
    {
      // reset the states of all slave solvers if any of them has been reset
      flowSolver.ResetStateToBeginningOfStep( domain );
      proppantSolver.ResetStateToBeginningOfStep( domain );
      ResetStateToBeginningOfStep( domain );
    }
    if( getLogLevel() >= 1 )
    {
      GEOSX_LOG_RANK_0( "\tIteration: " << iter+1  << ", FlowSolver: " );
    }

    dtReturnTemporary = flowSolver.NonlinearImplicitStep( time_n,
                                                          dtReturn,
                                                          cycleNumber,
                                                          domain,
                                                          flowSolver.getDofManager(),
                                                          flowSolver.getSystemMatrix(),
                                                          flowSolver.getSystemRhs(),
                                                          flowSolver.getSystemSolution() );


    if( dtReturnTemporary < dtReturn )
    {
      iter = 0;
      dtReturn = dtReturnTemporary;
      continue;
    }

    NonlinearSolverParameters const & fluidNonLinearParams = flowSolver.getNonlinearSolverParameters();
    if( fluidNonLinearParams.m_numNewtonIterations<=this->m_nonlinearSolverParameters.m_minIterNewton && iter > 0 &&
        getLogLevel() >= 1 )
    {
      GEOSX_LOG_RANK_0( "***** The iterative coupling has converged in " << iter  << " iterations! *****\n" );
      break;
    }

    if( getLogLevel()  >= 1 )
    {
      GEOSX_LOG_RANK_0( "\tIteration: " << iter+1  << ", Proppant Solver: " );
    }

    dtReturnTemporary = proppantSolver.NonlinearImplicitStep( time_n,
                                                              dtReturn,
                                                              cycleNumber,
                                                              domain,
                                                              proppantSolver.getDofManager(),
                                                              proppantSolver.getSystemMatrix(),
                                                              proppantSolver.getSystemRhs(),
                                                              proppantSolver.getSystemSolution() );

    if( dtReturnTemporary < dtReturn )
    {
      iter = 0;
      dtReturn = dtReturnTemporary;
      continue;
    }

    ++iter;
  }


  flowSolver.ImplicitStepComplete( time_n, dt, domain );
  proppantSolver.ImplicitStepComplete( time_n, dt, domain );

  proppantSolver.PostStepUpdate( time_n, dtReturn, cycleNumber, domain );
  this->ImplicitStepComplete( time_n, dt, domain );

  return dtReturn;
}

REGISTER_CATALOG_ENTRY( SolverBase, FlowProppantTransportSolver, std::string const &, Group * const )

} /* namespace geosx */
