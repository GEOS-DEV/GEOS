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

/*
 * ContactSolverBase.cpp
 */

#include "ContactSolverBase.hpp"

#include "SolidMechanicsEFEMKernels.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/contact/ContactSelector.hpp"
#include "constitutive/solid/ElasticIsotropic.hpp"
#include "finiteElement/elementFormulations/FiniteElementBase.hpp"
#include "linearAlgebra/utilities/LAIHelperFunctions.hpp"
#include "mesh/DomainPartition.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "mesh/NodeManager.hpp"
#include "mesh/SurfaceElementRegion.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"
#include "common/GEOS_RAJA_Interface.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;

ContactSolverBase::ContactSolverBase( const string & name,
                                                                  Group * const parent ):
  SolverBase( name, parent ),
  m_solidSolverName(),
  m_fractureRegionName(),
  m_solidSolver( nullptr )
{
  registerWrapper( viewKeyStruct::solidSolverNameString(), &m_solidSolverName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of the solid mechanics solver in the rock matrix" );

  registerWrapper( viewKeyStruct::fractureRegionNameString(), &m_fractureRegionName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of the fracture region." );

  registerWrapper( viewKeyStruct::contactRelationNameString(), &m_contactRelationName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of contact relation to enforce constraints on fracture boundary." );

  this->getWrapper< string >( viewKeyStruct::discretizationString() ).
    setInputFlag( InputFlags::FALSE );

}

ContactSolverBase::~ContactSolverBase()
{
  // TODO Auto-generated destructor stub
}

void ContactSolverBase::postProcessInput()
{
  m_solidSolver = &this->getParent().getGroup< SolidMechanicsLagrangianFEM >( m_solidSolverName );
  SolverBase::postProcessInput();
}

void ContactSolverBase::registerDataOnMesh( dataRepository::Group & meshBodies )
{
  meshBodies.forSubGroups< MeshBody >( [&] ( MeshBody & meshBody )
  {
    MeshLevel & meshLevel = meshBody.getMeshLevel( 0 );

    ElementRegionManager & elemManager = meshLevel.getElemManager();
    {
      elemManager.forElementRegions< SurfaceElementRegion >( [&] ( SurfaceElementRegion & region )
      {
        region.forElementSubRegions< EmbeddedSurfaceSubRegion >( [&]( EmbeddedSurfaceSubRegion & subRegion )
        {

          subRegion.registerWrapper< array2d< real64 > >( viewKeyStruct::dispJumpString() ).
            setPlotLevel( PlotLevel::LEVEL_0 ).
            reference().resizeDimension< 1 >( 3 );

          subRegion.registerWrapper< array2d< real64 > >( viewKeyStruct::deltaDispJumpString() ).
            reference().resizeDimension< 1 >( 3 );

          subRegion.registerWrapper< array2d< real64 > >( viewKeyStruct::oldDispJumpString() ).
            reference().resizeDimension< 1 >( 3 );

          subRegion.registerWrapper< array2d< real64 > >( viewKeyStruct::fractureTractionString() ).
            reference().resizeDimension< 1 >( 3 );

          subRegion.registerWrapper< array3d< real64 > >( viewKeyStruct::dTraction_dJumpString() ).
            reference().resizeDimension< 1, 2 >( 3, 3 );
        } );
      } );
    }
  } );
}

real64 ContactSolverBase::nonlinearImplicitStep( real64 const & time_n,
                                                 real64 const & dt,
                                                 integer const cycleNumber,
                                                 DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;
  // dt may be cut during the course of this step, so we are keeping a local
  // value to track the achieved dt for this step.
  real64 stepDt = dt;

  integer const maxNewtonIter = m_nonlinearSolverParameters.m_maxIterNewton;
  integer const minNewtonIter = m_nonlinearSolverParameters.m_minIterNewton;
  real64 const newtonTol = m_nonlinearSolverParameters.m_newtonTol;

  integer const maxNumberDtCuts = m_nonlinearSolverParameters.m_maxTimeStepCuts;
  real64 const dtCutFactor = m_nonlinearSolverParameters.m_timeStepCutFactor;

  bool const allowNonConverged = m_nonlinearSolverParameters.m_allowNonConverged > 0;

  integer & dtAttempt = m_nonlinearSolverParameters.m_numdtAttempts;

  // a flag to denote whether we have converged
  bool isNewtonConverged = false;

  bool isActiveSetConverged = false;

  // outer loop attempts to apply full timestep, and manages timestep cut if required.
  for( dtAttempt = 0; dtAttempt < maxNumberDtCuts; ++dtAttempt )
  {
    // reset the solver state, since we are restarting the time step
    if( dtAttempt > 0 )
    {
      resetStateToBeginningOfStep( domain );
      globalIndex numStick, numSlip, numOpen;
      computeFractureStateStatistics( domain, numStick, numSlip, numOpen, true );
    }

    bool useElasticStep = !isFractureAllInStickCondition( domain );

    integer & activeSetIter = m_activeSetIter;
    for( activeSetIter = 0; activeSetIter < m_activeSetMaxIter; ++activeSetIter )
    {
      // *******************************
      // Newton loop: begin
      // *******************************
      isNewtonConverged = false;
      // keep residual from previous iteration in case we need to do a line search
      real64 lastResidual = 1e99;
      integer & newtonIter = m_nonlinearSolverParameters.m_numNewtonIterations;
      real64 scaleFactor = 1.0;

      // main Newton loop
      bool computeResidual = true;
      for( newtonIter = 0; newtonIter < maxNewtonIter; ++newtonIter )
      {
        GEOSX_LOG_LEVEL_RANK_0( 1, GEOSX_FMT( "    Attempt: {:2}, ActiveSetIter: {:2} ; NewtonIter: {:2} ; ", dtAttempt, activeSetIter, newtonIter ) );

        // zero out matrix/rhs before assembly
        m_localMatrix.zero();
        m_localRhs.zero();

        // call assemble to fill the matrix and the rhs
        assembleSystem( time_n,
                        stepDt,
                        domain,
                        m_dofManager,
                        m_localMatrix.toViewConstSizes(),
                        m_localRhs );

        // apply boundary conditions to system
        applyBoundaryConditions( time_n,
                                 stepDt,
                                 domain,
                                 m_dofManager,
                                 m_localMatrix.toViewConstSizes(),
                                 m_localRhs );

        // TODO: maybe add scale function here?
        // Scale()

        real64 residualNorm;
        // get residual norm
        if( computeResidual )
        {
          residualNorm = calculateResidualNorm( domain, m_dofManager, m_localRhs );
        }
        else
        {
          residualNorm = lastResidual;
        }

        if( getLogLevel() >= 1 && logger::internal::rank==0 )
        {
          if( newtonIter!=0 )
          {
            std::cout << GEOSX_FMT( "Last LinSolve(iter,tol) = ({:4}, {:4.2e}) ; ",
                                    m_linearSolverResult.numIterations,
                                    m_linearSolverResult.residualReduction );
          }
          std::cout << std::endl;
        }

        // if the residual norm is less than the Newton tolerance we denote that we have
        // converged and break from the Newton loop immediately.
        if( residualNorm < newtonTol && newtonIter >= minNewtonIter )
        {
          isNewtonConverged = true;
          break;
        }

        // if using adaptive Krylov tolerance scheme, update tolerance.
        LinearSolverParameters::Krylov & krylovParams = m_linearSolverParameters.get().krylov;
        if( krylovParams.useAdaptiveTol )
        {
          krylovParams.relTolerance = eisenstatWalker( residualNorm, lastResidual, krylovParams.weakestTol );
        }

        // Compose parallel LA matrix/rhs out of local LA matrix/rhs
        m_matrix.create( m_localMatrix.toViewConst(), MPI_COMM_GEOSX );
        m_rhs.create( m_localRhs, MPI_COMM_GEOSX );
        m_solution.createWithLocalSize( m_matrix.numLocalCols(), MPI_COMM_GEOSX );

        // Output the linear system matrix/rhs for debugging purposes
        debugOutputSystem( time_n, cycleNumber, newtonIter, m_matrix, m_rhs );

        // Solve the linear system
        solveSystem( m_dofManager, m_matrix, m_rhs, m_solution );

        // Output the linear system solution for debugging purposes
        debugOutputSolution( time_n, cycleNumber, newtonIter, m_solution );

        // Copy solution from parallel vector back to local
        // TODO: This step will not be needed when we teach LA vectors to wrap our pointers
        m_solution.extract( m_localSolution );

        scaleFactor = scalingForSystemSolution( domain, m_dofManager, m_localSolution );

        // do line search in case residual has increased
        if( m_nonlinearSolverParameters.m_lineSearchAction != NonlinearSolverParameters::LineSearchAction::None && newtonIter > 0 )
        {
          bool lineSearchSuccess = lineSearch( time_n,
                                               stepDt,
                                               cycleNumber,
                                               domain,
                                               m_dofManager,
                                               m_localMatrix.toViewConstSizes(),
                                               m_localRhs,
                                               m_localSolution,
                                               scaleFactor,
                                               residualNorm );

          if( !lineSearchSuccess )
          {
            if( m_nonlinearSolverParameters.m_lineSearchAction == NonlinearSolverParameters::LineSearchAction::Attempt )
            {
              GEOSX_LOG_LEVEL_RANK_0( 1, "        Line search failed to produce reduced residual. Accepting iteration." );
            }
            else if( m_nonlinearSolverParameters.m_lineSearchAction == NonlinearSolverParameters::LineSearchAction::Require )
            {
              // if line search failed, then break out of the main Newton loop. Timestep will be cut.
              GEOSX_LOG_LEVEL_RANK_0( 1, "        Line search failed to produce reduced residual. Exiting Newton Loop." );
              break;
            }
          }
          // Residual norm already computed in line search and stored in "residualNorm"
          computeResidual = false;
        }
        else
        {
          // apply the system solution to the fields/variables
          applySystemSolution( m_dofManager, m_localSolution, scaleFactor, domain );
          // Need to compute the residual norm
          computeResidual = true;
        }

        if( !checkSystemSolution( domain, m_dofManager, m_localSolution, scaleFactor ) )
        {
          // TODO try chopping (similar to line search)
          GEOSX_LOG_RANK_0( "    Solution check failed. Newton loop terminated." );
          break;
        }

        lastResidual = residualNorm;
      }
      // *******************************
      // Newton loop: end
      // *******************************

      // *******************************
      // Active set check: begin
      // *******************************
      bool const isPreviousFractureStateValid = updateFractureState( domain );
      GEOSX_LOG_LEVEL_RANK_0( 1, "active set flag: " << std::boolalpha << isPreviousFractureStateValid );

      if( getLogLevel() >= 1 )
      {
        globalIndex numStick, numSlip, numOpen;
        computeFractureStateStatistics( domain, numStick, numSlip, numOpen, true );
      }
      // *******************************
      // Active set check: end
      // *******************************

      GEOSX_LOG_LEVEL_RANK_0( 1, "isPreviousFractureStateValid: " << std::boolalpha << isPreviousFractureStateValid <<
                              " | isNewtonConverged: " << isNewtonConverged << " | useElasticStep: " << useElasticStep );
      if( isNewtonConverged )
      {
        isActiveSetConverged = isPreviousFractureStateValid;
        if( isActiveSetConverged )
        {
          break;
        }
      }
      else if( useElasticStep )
      {
        GEOSX_LOG_LEVEL_RANK_0( 1, "Trying with an elastic step" );
        useElasticStep = false;
        resetStateToBeginningOfStep( domain );
        setFractureStateForElasticStep( domain );
      }
      else
      {
        GEOSX_LOG_LEVEL_RANK_0( 1, "Newton did not converge in active set loop" );
        break;
      }
    }
    if( !isNewtonConverged )
    {
      // cut timestep, go back to beginning of step and restart the Newton loop
      stepDt *= dtCutFactor;
      GEOSX_LOG_LEVEL_RANK_0 ( 1, "New dt = " <<  stepDt );
    }
    if( isActiveSetConverged )
    {
      break;
    }
  }

  if( !isNewtonConverged )
  {
    GEOSX_LOG_RANK_0( "Convergence not achieved." );

    if( allowNonConverged )
    {
      GEOSX_LOG_RANK_0( "The accepted solution may be inaccurate." );
    }
    else
    {
      GEOSX_ERROR( "Nonconverged solutions not allowed. Terminating..." );
    }
  }

  if( !isActiveSetConverged )
  {
    GEOSX_ERROR( "Active set did not reach a solution. Terminating..." );
  }
  else
  {
    GEOSX_LOG_RANK_0( "Number of active set iterations: " << m_activeSetIter );
  }

  // return the achieved timestep
  return stepDt;
}


} /* namespace geosx */
