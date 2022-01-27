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
 * @file LagrangianContactFlowSolver.cpp
 *
 */

#include "LagrangianContactFlowSolver.hpp"

#include "common/TimingMacros.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/fluid/SingleFluidBase.hpp"
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "mesh/DomainPartition.hpp"
#include "discretizationMethods/NumericalMethodsManager.hpp"
#include "mesh/SurfaceElementRegion.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "mesh/mpiCommunications/NeighborCommunicator.hpp"
#include "physicsSolvers/multiphysics/LagrangianContactSolver.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBase.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"
#include "physicsSolvers/surfaceGeneration/SurfaceGenerator.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "linearAlgebra/utilities/LAIHelperFunctions.hpp"
#include "linearAlgebra/solvers/PreconditionerJacobi.hpp"
#include "linearAlgebra/solvers/PreconditionerBlockJacobi.hpp"
#include "linearAlgebra/solvers/BlockPreconditionerGeneral.hpp"
#include "linearAlgebra/solvers/SeparateComponentPreconditioner.hpp"
#include "math/interpolation/Interpolation.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;
using namespace interpolation;

LagrangianContactFlowSolver::LagrangianContactFlowSolver( const std::string & name,
                                                          Group * const parent ):
  SinglePhasePoromechanicsSolver( name, parent ),
  m_contactSolverName(),
  m_flowSolverName(),
  m_contactSolver( nullptr ),
  m_flowSolver( nullptr ),
  m_stabilizationName( "" )
{
  registerWrapper( viewKeyStruct::contactSolverNameString(), &m_contactSolverName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of the contact solid mechanics solver to use in the contact with flow solver" );

  registerWrapper( viewKeyStruct::flowSolverNameString(), &m_flowSolverName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of the flow mechanics solver to use in the contact with flow solver" );

  registerWrapper( viewKeyStruct::stabilizationNameString(), &m_stabilizationName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of the stabilization to use in the lagrangian contact with flow solver" );

  this->getWrapper< string >( viewKeyStruct::discretizationString() ).
    setInputFlag( InputFlags::FALSE );

  /*
  m_linearSolverParameters.get().mgr.strategy = "LagrangianContactMechanicsFlow";
  m_linearSolverParameters.get().mgr.separateComponents = true;
  m_linearSolverParameters.get().mgr.displacementFieldName = keys::TotalDisplacement;
  m_linearSolverParameters.get().dofsPerNode = 3;
  */
}

void LagrangianContactFlowSolver::registerDataOnMesh( Group & )
{}

void LagrangianContactFlowSolver::initializePreSubGroups()
{
  SinglePhasePoromechanicsSolver::initializePreSubGroups();
}

void LagrangianContactFlowSolver::implicitStepSetup( real64 const & time_n,
                                                     real64 const & dt,
                                                     DomainPartition & domain )
{
  m_contactSolver->implicitStepSetup( time_n, dt, domain );
  m_flowSolver->implicitStepSetup( time_n, dt, domain );
}

void LagrangianContactFlowSolver::setupSystem( DomainPartition & domain,
                                               DofManager & dofManager,
                                               CRSMatrix< real64, globalIndex > & localMatrix,
                                               ParallelVector & rhs,
                                               ParallelVector & solution,
                                               bool const setSparsity )
{

  GEOSX_UNUSED_VAR( setSparsity );

  if( m_precond )
  {
    m_precond->clear();
  }

  // Laura
  // Uncomment these lines to run the lagrangiancontactflow solver alone
  // Comment them to run the fully coupled porolagrangian solver

//  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );
//  m_flowSolver->resetViews( mesh );
//
//  dofManager.setMesh( mesh );
//
//  setupDofs( domain, dofManager );
//  dofManager.reorderByRank();
  // Laura
  
  localIndex const numLocalRows = dofManager.numLocalDofs();

  SparsityPattern< globalIndex > patternOriginal;
  dofManager.setSparsityPattern( patternOriginal );

  // Get the original row lengths (diagonal blocks only)
  array1d< localIndex > rowLengths( patternOriginal.numRows() );
  for( localIndex localRow = 0; localRow < patternOriginal.numRows(); ++localRow )
  {
    rowLengths[localRow] = patternOriginal.numNonZeros( localRow );
  }

  // Add the number of nonzeros induced by coupling
  addTransmissibilityCouplingNNZ( domain, dofManager, rowLengths.toView() );

  // Create a new pattern with enough capacity for coupled matrix
  SparsityPattern< globalIndex > pattern;
  pattern.resizeFromRowCapacities< parallelHostPolicy >( patternOriginal.numRows(),
                                                         patternOriginal.numColumns(),
                                                         rowLengths.data() );

  // Copy the original nonzeros
  for( localIndex localRow = 0; localRow < patternOriginal.numRows(); ++localRow )
  {
    globalIndex const * cols = patternOriginal.getColumns( localRow ).dataIfContiguous();
    pattern.insertNonZeros( localRow, cols, cols + patternOriginal.numNonZeros( localRow ) );
  }

  // Add the nonzeros from coupling
  addTransmissibilityCouplingPattern( domain, dofManager, pattern.toView() );

  localMatrix.assimilate< parallelDevicePolicy<> >( std::move( pattern ) );
  localMatrix.setName( this->getName() + "/matrix" );

  rhs.setName( this->getName() + "/rhs" );
  rhs.create( numLocalRows, MPI_COMM_GEOSX );

  solution.setName( this->getName() + "/solution" );
  solution.create( numLocalRows, MPI_COMM_GEOSX );

  setUpDflux_dApertureMatrix( domain, dofManager, localMatrix );

  if( !m_precond && m_linearSolverParameters.get().solverType != LinearSolverParameters::SolverType::direct )
  {
    createPreconditioner( domain );
  }
}

void LagrangianContactFlowSolver::implicitStepComplete( real64 const & time_n,
                                                        real64 const & dt,
                                                        DomainPartition & domain )
{
  m_contactSolver->implicitStepComplete( time_n, dt, domain );
  m_flowSolver->implicitStepComplete( time_n, dt, domain );
}

void LagrangianContactFlowSolver::postProcessInput()
{
  m_contactSolver = &this->getParent().getGroup< LagrangianContactSolver >( m_contactSolverName );
  m_flowSolver = &this->getParent().getGroup< SinglePhaseBase >( m_flowSolverName );

  SinglePhasePoromechanicsSolver::postProcessInput();
}

void LagrangianContactFlowSolver::initializePostInitialConditionsPreSubGroups()
{}

LagrangianContactFlowSolver::~LagrangianContactFlowSolver()
{
  // TODO Auto-generated destructor stub
}

void LagrangianContactFlowSolver::resetStateToBeginningOfStep( DomainPartition & domain )
{
  m_contactSolver->resetStateToBeginningOfStep( domain );
  m_flowSolver->resetStateToBeginningOfStep( domain );
  updateOpeningForFlow( domain );
}

real64 LagrangianContactFlowSolver::solverStep( real64 const & time_n,
                                                real64 const & dt,
                                                int const cycleNumber,
                                                DomainPartition & domain )
{
  real64 dtReturn = dt;

  implicitStepSetup( time_n,
                     dt,
                     domain );

  // Need this strange call to allow the two physics solvers to initialize internal data structures
  // (like derivativeFluxResidual_dAperture in m_flowSolver) but I need to use my dofManager
  DofManager localDofManager( "localDofManager" );
  m_contactSolver->setupSystem( domain, localDofManager, m_localMatrix, m_rhs, m_solution );
  m_flowSolver->setupSystem( domain, localDofManager, m_localMatrix, m_rhs, m_solution );
  localDofManager.clear();

  setupSystem( domain,
               m_dofManager,
               m_localMatrix,
               m_rhs,
               m_solution );

  // currently the only method is implicit time integration
  dtReturn = this->nonlinearImplicitStep( time_n, dt, cycleNumber, domain );

  // final step for completion of timestep. typically secondary variable updates and cleanup.
  implicitStepComplete( time_n, dtReturn, domain );

  return dtReturn;
}

void LagrangianContactFlowSolver::updateOpeningForFlow( DomainPartition & domain )
{
  MeshLevel & meshLevel = domain.getMeshBody( 0 ).getMeshLevel( 0 );
  ElementRegionManager & elemManager = meshLevel.getElemManager();

  elemManager.forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion & subRegion )
  {
    if( subRegion.hasWrapper( m_pressureKey ) )
    {
      arrayView1d< integer const > const & ghostRank = subRegion.ghostRank();
      arrayView2d< real64 const > const & dispJump = subRegion.getReference< array2d< real64 > >( m_dispJumpKey );
      arrayView1d< real64 const > const & area = subRegion.getElementArea().toViewConst();
      arrayView1d< real64 const > const & volume = subRegion.getElementVolume();
      arrayView1d< real64 > const &
      aperture = subRegion.getReference< array1d< real64 > >( extrinsicMeshData::flow::hydraulicAperture::key() );
      arrayView1d< real64 > const &
      deltaVolume = subRegion.getReference< array1d< real64 > >( extrinsicMeshData::flow::deltaVolume::key() );

      forAll< serialPolicy >( subRegion.size(), [&]( localIndex const kfe )
      {
        if( ghostRank[kfe] < 0 )
        {
          if( m_contactSolver->isElementInOpenState( subRegion, kfe ) )
          {
            aperture[kfe] = dispJump[kfe][0];
            deltaVolume[kfe] = aperture[kfe] * area[kfe] - volume[kfe];
          }
          else
          {
            aperture[kfe] = 0.0;
            deltaVolume[kfe] = 0.0;
          }
        }
      } );
    }
  } );

  return;
}

real64 LagrangianContactFlowSolver::splitOperatorStep( real64 const & GEOSX_UNUSED_PARAM( time_n ),
                                                       real64 const & dt,
                                                       integer const GEOSX_UNUSED_PARAM( cycleNumber ),
                                                       DomainPartition & GEOSX_UNUSED_PARAM( domain ) )
{
  real64 dtReturn = dt;
  return dtReturn;
}

real64 LagrangianContactFlowSolver::explicitStep( real64 const & GEOSX_UNUSED_PARAM( time_n ),
                                                  real64 const & dt,
                                                  const int GEOSX_UNUSED_PARAM( cycleNumber ),
                                                  DomainPartition & GEOSX_UNUSED_PARAM( domain ) )
{
  GEOSX_MARK_FUNCTION;
  GEOSX_ERROR( "ExplicitStep non available for LagrangianContactFlowSolver!" );
  return dt;
}

real64 LagrangianContactFlowSolver::nonlinearImplicitStep( real64 const & time_n,
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

  bool useElasticStep = !m_contactSolver->isFractureAllInStickCondition( domain );

  // outer loop attempts to apply full timestep, and managed the cutting of the timestep if
  // required.
  for( dtAttempt = 0; dtAttempt < maxNumberDtCuts; ++dtAttempt )
  {
    // reset the solver state, since we are restarting the time step
    if( dtAttempt > 0 )
    {
      resetStateToBeginningOfStep( domain );
      globalIndex numStick, numSlip, numOpen;
      m_contactSolver->computeFractureStateStatistics( domain, numStick, numSlip, numOpen, true );
    }

    integer & activeSetIter = m_activeSetIter;
    for( activeSetIter = 0; activeSetIter < m_contactSolver->getActiveSetMaxIter(); ++activeSetIter )
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
        if( getLogLevel() >= 1 && logger::internal::rank==0 )
        {
          char output[55] = {0};
          sprintf( output, "    Attempt: %2d, ActiveSetIter: %2d ; NewtonIter: %2d ; ",
                   dtAttempt, activeSetIter, newtonIter );
          std::cout<<output<<std::endl;
        }

        // zero out matrix/rhs before assembly
        m_localMatrix.zero();
        m_rhs.zero();

        {
          arrayView1d< real64 > const localRhs = m_rhs.open();

          // call assemble to fill the matrix and the rhs
          assembleSystem( time_n,
                          stepDt,
                          domain,
                          m_dofManager,
                          m_localMatrix.toViewConstSizes(),
                          localRhs );

          // apply boundary conditions to system
          applyBoundaryConditions( time_n,
                                   stepDt,
                                   domain,
                                   m_dofManager,
                                   m_localMatrix.toViewConstSizes(),
                                   localRhs );

          m_rhs.close();
        }
        // TODO: maybe add scale function here?
        // Scale()

        real64 residualNorm;
        // get residual norm
        if( computeResidual )
        {
          residualNorm = calculateResidualNorm( domain, m_dofManager, m_rhs.values() );
        }
        else
        {
          residualNorm = lastResidual;
        }

        if( getLogLevel() >= 1 && logger::internal::rank==0 )
        {
          if( newtonIter!=0 )
          {
            char output[46] = {0};
            sprintf( output,
                     "Last LinSolve(iter,tol) = (%4d, %4.2e) ; ",
                     m_linearSolverResult.numIterations,
                     m_linearSolverResult.residualReduction );
            std::cout<<output;
          }
          std::cout<<std::endl;
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
        m_matrix.create( m_localMatrix.toViewConst(), m_dofManager.numLocalDofs(), MPI_COMM_GEOSX );

        // Output the linear system matrix/rhs for debugging purposes
        debugOutputSystem( time_n, cycleNumber, 10*activeSetIter+newtonIter, m_matrix, m_rhs );

        // Solve the linear system
        solveSystem( m_dofManager, m_matrix, m_rhs, m_solution );

        // Output the linear system solution for debugging purposes
        debugOutputSolution( time_n, cycleNumber, newtonIter, m_solution );

        scaleFactor = scalingForSystemSolution( domain, m_dofManager, m_solution.values() );

        // do line search in case residual has increased
        if( m_nonlinearSolverParameters.m_lineSearchAction != NonlinearSolverParameters::LineSearchAction::None && newtonIter > 0 )
        {
          bool lineSearchSuccess = lineSearch( time_n,
                                               stepDt,
                                               cycleNumber,
                                               domain,
                                               m_dofManager,
                                               m_localMatrix.toViewConstSizes(),
                                               m_rhs,
                                               m_solution,
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
          applySystemSolution( m_dofManager, m_solution.values(), scaleFactor, domain );
          // Need to compute the residual norm
          computeResidual = true;
        }

        if( !checkSystemSolution( domain, m_dofManager, m_solution.values(), scaleFactor ) )
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
      bool const isPreviousFractureStateValid = m_contactSolver->updateFractureState( domain );
      GEOSX_LOG_LEVEL_RANK_0( 1, "active set flag: " << std::boolalpha << isPreviousFractureStateValid );

      if( getLogLevel() > 2 )
      {
        globalIndex numStick, numSlip, numOpen;
        m_contactSolver->computeFractureStateStatistics( domain, numStick, numSlip, numOpen, true );
      }
      // *******************************
      // Active set check: end
      // *******************************

      GEOSX_LOG_LEVEL_RANK_0( 2, "isPreviousFractureStateValid: " << std::boolalpha << isPreviousFractureStateValid <<
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
        GEOSX_LOG_LEVEL_RANK_0( 2, "Trying with an elastic step" );
        useElasticStep = false;
        resetStateToBeginningOfStep( domain );
        m_contactSolver->setFractureStateForElasticStep( domain );
        updateOpeningForFlow( domain );
      }
      else
      {
        GEOSX_LOG_LEVEL_RANK_0( 2, "Newton did not converge in active set loop" );
        break;
      }
    }
    if( !isNewtonConverged )
    {
      // cut timestep, go back to beginning of step and restart the Newton loop
      // if (stepDt >= 0.25)
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
    //GEOSX_ERROR( "Active set did not reached a solution. Terminating..." );
  }
  else
  {
    GEOSX_LOG_RANK_0( "Number of active set iterations: " << m_activeSetIter );
  }

  // return the achieved timestep
  return stepDt;
}

bool LagrangianContactFlowSolver::lineSearch( real64 const & time_n,
                                              real64 const & dt,
                                              integer const GEOSX_UNUSED_PARAM( cycleNumber ),
                                              DomainPartition & domain,
                                              DofManager const & dofManager,
                                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                              ParallelVector & rhs,
                                              ParallelVector & solution,
                                              real64 const scaleFactor,
                                              real64 & lastResidual )
{
  bool lineSearchSuccess = true;

  integer const maxNumberLineSearchCuts = m_nonlinearSolverParameters.m_lineSearchMaxCuts;

  real64 const sigma1 = 0.5;
  real64 const alpha = 1.e-4;

  real64 localScaleFactor = scaleFactor;
  real64 lamm = scaleFactor;
  real64 lamc = localScaleFactor;
  integer lineSearchIteration = 0;

  // get residual norm
  real64 residualNorm0 = lastResidual;

  applySystemSolution( dofManager, solution.values(), scaleFactor, domain );

  // re-assemble system
  localMatrix.zero();
  rhs.zero();

  {
    arrayView1d< real64 > const localRhs = rhs.open();
    assembleSystem( time_n, dt, domain, dofManager, localMatrix, localRhs );
    applyBoundaryConditions( time_n, dt, domain, dofManager, localMatrix, localRhs );
    rhs.close();
  }

  // get residual norm
  real64 residualNormT = calculateResidualNorm( domain, dofManager, rhs.values() );

  real64 ff0 = residualNorm0*residualNorm0;
  real64 ffT = residualNormT*residualNormT;
  real64 ffm = ffT;
  real64 cumulativeScale = scaleFactor;

  while( residualNormT >= (1.0 - alpha*localScaleFactor)*residualNorm0 )
  {
    real64 const previousLocalScaleFactor = localScaleFactor;
    // Apply the three point parabolic model
    if( lineSearchIteration == 0 )
    {
      localScaleFactor *= sigma1;
    }
    else
    {
      localScaleFactor = parabolicInterpolationThreePoints( lamc, lamm, ff0, ffT, ffm );
    }

    // Update x; keep the books on lambda
    real64 const deltaLocalScaleFactor = ( localScaleFactor - previousLocalScaleFactor );
    cumulativeScale += deltaLocalScaleFactor;

    if( !checkSystemSolution( domain, dofManager, solution.values(), deltaLocalScaleFactor ) )
    {
      GEOSX_LOG_LEVEL_RANK_0( 1, "        Line search " << lineSearchIteration << ", solution check failed" );
      continue;
    }

    applySystemSolution( dofManager, solution.values(), deltaLocalScaleFactor, domain );
    lamm = lamc;
    lamc = localScaleFactor;

    // Keep the books on the function norms
    // re-assemble system
    // TODO: add a flag to avoid a completely useless Jacobian computation: rhs is enough
    localMatrix.zero();
    rhs.zero();

    {
      arrayView1d< real64 > const localRhs = rhs.open();
      assembleSystem( time_n, dt, domain, dofManager, localMatrix, localRhs );
      applyBoundaryConditions( time_n, dt, domain, dofManager, localMatrix, localRhs );
      rhs.close();
    }

    if( getLogLevel() >= 1 && logger::internal::rank==0 )
    {
      char output[100];
      sprintf( output, "        Line search @ %0.3f:      ", cumulativeScale );
      std::cout<<output;
    }

    // get residual norm
    residualNormT = calculateResidualNorm( domain, dofManager, rhs.values() );
    ffm = ffT;
    ffT = residualNormT*residualNormT;
    lineSearchIteration += 1;

    if( lineSearchIteration > maxNumberLineSearchCuts )
    {
      lineSearchSuccess = false;
      break;
    }
  }

  lastResidual = residualNormT;

  return lineSearchSuccess;
}

void LagrangianContactFlowSolver::setupDofs( DomainPartition const & domain,
                                             DofManager & dofManager ) const
{
  GEOSX_MARK_FUNCTION;
  m_contactSolver->setupDofs( domain, dofManager );
  m_flowSolver->setupDofs( domain, dofManager );

  // restrict coupling to fracture regions only
  ElementRegionManager const & elemManager = domain.getMeshBody( 0 ).getMeshLevel( 0 ).getElemManager();
  string_array fractureRegions;
  elemManager.forElementRegions< SurfaceElementRegion >( [&]( SurfaceElementRegion const & elementRegion )
  {
    fractureRegions.emplace_back( elementRegion.getName() );
  } );

  dofManager.addCoupling( keys::TotalDisplacement,
                          m_pressureKey,
                          DofManager::Connector::Elem,
                          fractureRegions );
}

void LagrangianContactFlowSolver::assembleSystem( real64 const time,
                                                  real64 const dt,
                                                  DomainPartition & domain,
                                                  DofManager const & dofManager,
                                                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                  arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  // Need to synchronize the two iteration counters
  m_contactSolver->getNonlinearSolverParameters().m_numNewtonIterations = m_nonlinearSolverParameters.m_numNewtonIterations;

  // SynchronizeFractureState is called in AssembleSystem and it is needed by:
  // - UpdateOpeningForFlow
  // - AssembleFluidMassResidualDerivativeWrtDisplacement
  m_contactSolver->assembleSystem( time,
                                   dt,
                                   domain,
                                   dofManager,
                                   localMatrix,
                                   localRhs );

  updateOpeningForFlow( domain );

  //m_flowSolver->assembleSystem( time,
  //                              dt,
  //                              domain,
  //                              dofManager,
  //                              localMatrix,
  //                              localRhs );

  //m_flowSolver->resetViews( domain.getMeshBody( 0 ).getMeshLevel( 0 ) );

  m_flowSolver->assembleAccumulationTerms< parallelDevicePolicy<> >( domain,
                                                                     dofManager,
                                                                     localMatrix,
                                                                     localRhs );
  m_flowSolver->assembleHydrofracFluxTerms( time,
                                            dt,
                                            domain,
                                            dofManager,
                                            localMatrix,
                                            localRhs,
                                            getDerivativeFluxResidual_dAperture() );

  assembleForceResidualDerivativeWrtPressure( domain, dofManager, localMatrix, localRhs );
  assembleFluidMassResidualDerivativeWrtDisplacement( domain, dofManager, localMatrix, localRhs );
  assembleStabilization( domain, dofManager, localMatrix, localRhs );

  getRefDerivativeFluxResidual_dAperture()->zero();
}

void LagrangianContactFlowSolver::applyBoundaryConditions( real64 const time,
                                                           real64 const dt,
                                                           DomainPartition & domain,
                                                           DofManager const & dofManager,
                                                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                           arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;
  m_contactSolver->applyBoundaryConditions( time,
                                            dt,
                                            domain,
                                            dofManager,
                                            localMatrix,
                                            localRhs );
  m_flowSolver->applyBoundaryConditions( time,
                                         dt,
                                         domain,
                                         dofManager,
                                         localMatrix,
                                         localRhs );
}

real64 LagrangianContactFlowSolver::calculateResidualNorm( DomainPartition const & domain,
                                                           DofManager const & dofManager,
                                                           arrayView1d< real64 const > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel const & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  string const & tracDofKey = dofManager.getKey( LagrangianContactSolver::viewKeyStruct::tractionString() );
  string const & presDofKey = dofManager.getKey( m_pressureKey );
  globalIndex const rankOffset = dofManager.rankOffset();

  NodeManager const & nodeManager = mesh.getNodeManager();

  arrayView1d< globalIndex const > const & dispDofNumber =
    nodeManager.getReference< array1d< globalIndex > >( dofManager.getKey( keys::TotalDisplacement ) );

  arrayView1d< integer const > const & elemGhostRank = nodeManager.ghostRank();

  RAJA::ReduceSum< parallelDeviceReduce, real64 > localSum0( 0.0 );
  forAll< parallelDevicePolicy<> >( nodeManager.size(),
                                    [localRhs, localSum0, dispDofNumber, rankOffset, elemGhostRank] GEOSX_HOST_DEVICE ( localIndex const k )
  {
    if( elemGhostRank[k] < 0 )
    {
      localIndex const localRow = LvArray::integerConversion< localIndex >( dispDofNumber[k] - rankOffset );
      for( localIndex dim = 0; dim < 3; ++dim )
      {
        localSum0 += localRhs[localRow + dim] * localRhs[localRow + dim];
      }
    }
  } );
  real64 const momentumR2 = localSum0.get();

  real64 contactR2 = 0.0;
  real64 pressureR2 = 0.0;

  forTargetSubRegions< FaceElementSubRegion >( mesh, [&]( localIndex const, FaceElementSubRegion const & subRegion )
  {
    arrayView1d< globalIndex const > const & tracDofNumber = subRegion.getReference< array1d< globalIndex > >( tracDofKey );
    arrayView1d< globalIndex const > const & presDofNumber = subRegion.getReference< array1d< globalIndex > >( presDofKey );
    arrayView1d< integer const > const & ghostRank = subRegion.ghostRank();

    RAJA::ReduceSum< parallelHostReduce, real64 > localSum1( 0.0 );
    RAJA::ReduceSum< parallelHostReduce, real64 > localSum2( 0.0 );
    forAll< parallelHostPolicy >( subRegion.size(), [=] ( localIndex const k )
    {
      if( ghostRank[k] < 0 )
      {
        localIndex const tracLocalRow = LvArray::integerConversion< localIndex >( tracDofNumber[k] - rankOffset );
        for( localIndex dim = 0; dim < 3; ++dim )
        {
          localSum1 += localRhs[tracLocalRow + dim] * localRhs[tracLocalRow + dim];
        }
        localIndex const presLocalRow = LvArray::integerConversion< localIndex >( presDofNumber[k] - rankOffset );
        localSum2 += localRhs[presLocalRow] * localRhs[presLocalRow];
      }
    } );
    contactR2 += localSum1.get();
    pressureR2 += localSum2.get();
  } );

  real64 localR2[3] = { momentumR2, contactR2, pressureR2 };
  real64 globalResidualNorm[4]{};

  int const rank = MpiWrapper::commRank( MPI_COMM_GEOSX );
  int const size = MpiWrapper::commSize( MPI_COMM_GEOSX );
  array1d< real64 > globalR2( 3 * size );
  globalR2.setValues< serialPolicy >( 0 );

  // Everything is done on rank 0
  MpiWrapper::gather( localR2,
                      3,
                      globalR2.data(),
                      3,
                      0,
                      MPI_COMM_GEOSX );

  if( rank==0 )
  {
    globalResidualNorm[0] = 0.0;
    globalResidualNorm[1] = 0.0;
    globalResidualNorm[2] = 0.0;
    for( int r=0; r<size; ++r )
    {
      // sum across all ranks
      globalResidualNorm[0] += globalR2[3 * r + 0];
      globalResidualNorm[1] += globalR2[3 * r + 1];
      globalResidualNorm[2] += globalR2[3 * r + 2];
    }
    globalResidualNorm[3] = globalResidualNorm[0] + globalResidualNorm[1] + globalResidualNorm[2];
    globalResidualNorm[0] = sqrt( globalResidualNorm[0] );
    globalResidualNorm[1] = sqrt( globalResidualNorm[1] );
    globalResidualNorm[2] = sqrt( globalResidualNorm[2] );
    globalResidualNorm[3] = sqrt( globalResidualNorm[3] );
  }

  MpiWrapper::bcast( globalResidualNorm, 4, 0, MPI_COMM_GEOSX );

  if( m_nonlinearSolverParameters.m_numNewtonIterations == 0 )
  {
    m_initialResidual[0] = globalResidualNorm[0];
    m_initialResidual[1] = globalResidualNorm[1];
    m_initialResidual[2] = globalResidualNorm[2];
    m_initialResidual[3] = globalResidualNorm[3];
    globalResidualNorm[0] = 1.0;
    globalResidualNorm[1] = 1.0;
    globalResidualNorm[2] = 1.0;
    globalResidualNorm[3] = 1.0;
  }
  else
  {
    globalResidualNorm[0] /= (m_initialResidual[0]+1.0);
    globalResidualNorm[1] /= (m_initialResidual[1]+1.0);
    globalResidualNorm[2] /= (m_initialResidual[2]+1.0);
    globalResidualNorm[3] /= (m_initialResidual[3]+0.0);
  }

  char output[122] = {0};
  sprintf( output,
           "( Rdisplacement, Rtraction, Rpressure, Rtotal ) = ( %15.6e, %15.6e, %15.6e, %15.6e );",
           globalResidualNorm[0],
           globalResidualNorm[1],
           globalResidualNorm[2],
           globalResidualNorm[3] );
  GEOSX_LOG_RANK_0( output );

  return globalResidualNorm[3];
}

void LagrangianContactFlowSolver::createPreconditioner( DomainPartition const & domain )
{
  if( m_linearSolverParameters.get().preconditionerType == LinearSolverParameters::PreconditionerType::block )
  {
    // TODO: move among inputs (xml)
    string const leadingBlockApproximation = "blockJacobi";
    //string const leadingBlockApproximation = "jacobi";

    LinearSolverParameters mechParams = m_contactSolver->getSolidSolver()->getLinearSolverParameters();
    // Because of boundary conditions
    mechParams.isSymmetric = false;

    std::vector< SchurComplementOption > schurOptions( 2 );
    std::unique_ptr< BlockPreconditionerGeneral< LAInterface > > precond;
    std::unique_ptr< PreconditionerBase< LAInterface > > tracPrecond;
    std::unique_ptr< PreconditionerBase< LAInterface > > flowPrecond;

    if( leadingBlockApproximation == "jacobi" )
    {
      schurOptions[0] = SchurComplementOption::Diagonal;

      // Using LAI implementation of Jacobi preconditioner
      LinearSolverParameters tracParams;
      tracParams.preconditionerType = LinearSolverParameters::PreconditionerType::jacobi;
      tracPrecond = LAInterface::createPreconditioner( tracParams );
    }
    else if( leadingBlockApproximation == "blockJacobi" )
    {
      schurOptions[0] = SchurComplementOption::UserDefined;

      tracPrecond = std::make_unique< PreconditionerBlockJacobi< LAInterface > >( mechParams.dofsPerNode );
    }
    else
    {
      GEOSX_ERROR( "LagrangianContactSolver::CreatePreconditioner leadingBlockApproximation option " << leadingBlockApproximation << " not supported" );
    }

    // Flow + Jacobi: using LAI implementation of Jacobi preconditioner
    schurOptions[1] = SchurComplementOption::Diagonal;

    LinearSolverParameters flowParams = m_flowSolver->getLinearSolverParameters();
    flowPrecond = LAInterface::createPreconditioner( flowParams );

    precond = std::make_unique< BlockPreconditionerGeneral< LAInterface > >( 3,
                                                                             BlockShapeOption::LowerUpperTriangular,
                                                                             schurOptions );

    precond->setupBlock( 0,
                         { { m_tractionKey, { 3, 0, 3 } } },
                         std::move( tracPrecond ) );

    precond->setupBlock( 1,
                         { { m_pressureKey, { 1, 0, 1 } } },
                         std::move( flowPrecond ) );

    if( mechParams.amg.nullSpaceType == LinearSolverParameters::AMG::NullSpaceType::rigidBodyModes )
    {
      if( m_contactSolver->getSolidSolver()->getRigidBodyModes().empty() )
      {
        MeshLevel const & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );
        LAIHelperFunctions::computeRigidBodyModes( mesh,
                                                   m_dofManager,
                                                   { keys::TotalDisplacement },
                                                   m_contactSolver->getSolidSolver()->getRigidBodyModes() );
      }
    }

    // Preconditioner for the Schur complement: mechPrecond
    std::unique_ptr< PreconditionerBase< LAInterface > > mechPrecond = LAInterface::createPreconditioner( mechParams, m_contactSolver->getSolidSolver()->getRigidBodyModes() );
    precond->setupBlock( 2,
                         { { keys::TotalDisplacement, { 3, 0, 3 } } },
                         std::move( mechPrecond ) );

    m_precond = std::move( precond );
  }
  else if( m_linearSolverParameters.get().preconditionerType == LinearSolverParameters::PreconditionerType::block_fs )
  {
    // TODO: move among inputs (xml)
    string const leadingBlockApproximation = "blockJacobi";
    //string const leadingBlockApproximation = "jacobi";

    LinearSolverParameters mechParams = m_contactSolver->getSolidSolver()->getLinearSolverParameters();
    // Because of boundary conditions
    mechParams.isSymmetric = false;

    std::vector< SchurComplementOption > schurOptions( 2 );
    std::unique_ptr< BlockPreconditionerGeneral< LAInterface > > precond;
    std::unique_ptr< PreconditionerBase< LAInterface > > tracPrecond;
    std::unique_ptr< PreconditionerBase< LAInterface > > flowPrecond;

    if( leadingBlockApproximation == "jacobi" )
    {
      schurOptions[0] = SchurComplementOption::Diagonal;

      // Using LAI implementation of Jacobi preconditioner
      LinearSolverParameters tracParams;
      tracParams.preconditionerType = LinearSolverParameters::PreconditionerType::jacobi;
      tracPrecond = LAInterface::createPreconditioner( tracParams );
    }
    else if( leadingBlockApproximation == "blockJacobi" )
    {
      schurOptions[0] = SchurComplementOption::UserDefined;

      tracPrecond = std::make_unique< PreconditionerBlockJacobi< LAInterface > >( mechParams.dofsPerNode );
    }
    else
    {
      GEOSX_ERROR( "LagrangianContactSolver::CreatePreconditioner leadingBlockApproximation option " << leadingBlockApproximation << " not supported" );
    }

    // Mechanics + Jacobi: using LAI implementation of Jacobi preconditioner
    schurOptions[1] = SchurComplementOption::Diagonal;

    LinearSolverParameters flowParams = m_flowSolver->getLinearSolverParameters();
    flowPrecond = LAInterface::createPreconditioner( flowParams );

    precond = std::make_unique< BlockPreconditionerGeneral< LAInterface > >( 3,
                                                                             BlockShapeOption::LowerUpperTriangular,
                                                                             schurOptions );

    precond->setupBlock( 0,
                         { { m_tractionKey, { 3, 0, 3 } } },
                         std::move( tracPrecond ) );

    if( mechParams.amg.nullSpaceType == LinearSolverParameters::AMG::NullSpaceType::rigidBodyModes )
    {
      if( m_contactSolver->getSolidSolver()->getRigidBodyModes().empty() )
      {
        MeshLevel const & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );
        LAIHelperFunctions::computeRigidBodyModes( mesh,
                                                   m_dofManager,
                                                   { keys::TotalDisplacement },
                                                   m_contactSolver->getSolidSolver()->getRigidBodyModes() );
      }
    }

    // Preconditioner for the Schur complement: mechPrecond
    std::unique_ptr< PreconditionerBase< LAInterface > > mechPrecond = LAInterface::createPreconditioner( mechParams, m_contactSolver->getSolidSolver()->getRigidBodyModes() );
    precond->setupBlock( 1,
                         { { keys::TotalDisplacement, { 3, 0, 3 } } },
                         std::move( mechPrecond ) );

    precond->setupBlock( 2,
                         { { m_pressureKey, { 1, 0, 1 } } },
                         std::move( flowPrecond ) );

    m_precond = std::move( precond );
  }
  else
  {
    //TODO: Revisit this part such that is coherent across physics solver
    //m_precond = LAInterface::createPreconditioner( m_linearSolverParameters.get() );
  }
}

void LagrangianContactFlowSolver::
  addTransmissibilityCouplingNNZ( DomainPartition const & domain,
                                  DofManager const & dofManager,
                                  arrayView1d< localIndex > const & rowLengths ) const
{
  GEOSX_MARK_FUNCTION;

  std::cout << "in LagrangianContactFlowSolver::addTransmissibilityCouplingNNZ\n";

  MeshLevel const & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  ElementRegionManager const & elemManager = mesh.getElemManager();

  string const presDofKey = dofManager.getKey( m_pressureKey );

  globalIndex const rankOffset = dofManager.rankOffset();

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & stabilizationMethod = fvManager.getFluxApproximation( m_stabilizationName );

  stabilizationMethod.forStencils< SurfaceElementStencil >( mesh, [&]( SurfaceElementStencil const & stencil )
  {
    for( localIndex iconn=0; iconn<stencil.size(); ++iconn )
    {
      localIndex const numFluxElems = stencil.stencilSize( iconn );
      typename SurfaceElementStencil::IndexContainerViewConstType const & seri = stencil.getElementRegionIndices();
      typename SurfaceElementStencil::IndexContainerViewConstType const & sesri = stencil.getElementSubRegionIndices();
      typename SurfaceElementStencil::IndexContainerViewConstType const & sei = stencil.getElementIndices();

      FaceElementSubRegion const & elementSubRegion =
        elemManager.getRegion( seri[iconn][0] ).getSubRegion< FaceElementSubRegion >( sesri[iconn][0] );

      ArrayOfArraysView< localIndex const > const elemsToNodes = elementSubRegion.nodeList().toViewConst();

      arrayView1d< globalIndex const > const faceElementDofNumber =
        elementSubRegion.getReference< array1d< globalIndex > >( presDofKey );

      for( localIndex k0=0; k0<numFluxElems; ++k0 )
      {
        globalIndex const activeFlowDOF = faceElementDofNumber[sei[iconn][k0]];
        globalIndex const rowNumber = activeFlowDOF - rankOffset;

        if( rowNumber >= 0 && rowNumber < rowLengths.size() )
        {
          for( localIndex k1=0; k1<numFluxElems; ++k1 )
          {
            // The coupling with the nodal displacements of the cell itself has already been added by the dofManager
            // so we only add the coupling with the nodal displacements of the neighbours.
            if( k1 != k0 )
            {
              localIndex const numNodesPerElement = elemsToNodes[sei[iconn][k1]].size();
              rowLengths[rowNumber] += 3*numNodesPerElement;
            }
          }
        }
      }
    }
  } );
}

void LagrangianContactFlowSolver::
  addTransmissibilityCouplingPattern( DomainPartition const & domain,
                                      DofManager const & dofManager,
                                      SparsityPatternView< globalIndex > const & pattern ) const
{
  GEOSX_MARK_FUNCTION;

  MeshLevel const & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  FaceManager const & faceManager = mesh.getFaceManager();
  NodeManager const & nodeManager = mesh.getNodeManager();
  ElementRegionManager const & elemManager = mesh.getElemManager();

  string const dispDofKey = dofManager.getKey( keys::TotalDisplacement );
  string const presDofKey = dofManager.getKey( m_pressureKey );

  arrayView1d< globalIndex const > const &
  dispDofNumber = nodeManager.getReference< globalIndex_array >( dispDofKey );
  ArrayOfArraysView< localIndex const > const & faceToNodeMap = faceManager.nodeList().toViewConst();

  // Get the finite volume method used to compute the stabilization
  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & stabilizationMethod = fvManager.getFluxApproximation( m_stabilizationName );

  // Form the SurfaceGenerator, get the fracture name and use it to retrieve the faceMap (from fracture element to face)
  SurfaceGenerator const &
  surfaceGenerator = this->getParent().getGroup< SurfaceGenerator >( "SurfaceGen" );
  SurfaceElementRegion const & fractureRegion = elemManager.getRegion< SurfaceElementRegion >( surfaceGenerator.getFractureRegionName() );
  FaceElementSubRegion const & fractureSubRegion = fractureRegion.getSubRegion< FaceElementSubRegion >( "faceElementSubRegion" );
  GEOSX_ERROR_IF( !fractureSubRegion.hasWrapper( m_pressureKey ), "The fracture subregion must contain pressure field." );
  arrayView2d< localIndex const > const faceMap = fractureSubRegion.faceList();
  GEOSX_ERROR_IF( faceMap.size( 1 ) != 2, "A fracture face has to be shared by two cells." );

  arrayView1d< globalIndex const > const &
  presDofNumber = fractureSubRegion.getReference< globalIndex_array >( presDofKey );

  arrayView2d< localIndex const > const & elemsToFaces = fractureSubRegion.faceList();

  stabilizationMethod.forStencils< SurfaceElementStencil >( mesh, [&]( SurfaceElementStencil const & stencil )
  {
    forAll< serialPolicy >( stencil.size(), [=] ( localIndex const iconn )
    {
      localIndex const numFluxElems = stencil.stencilSize( iconn );

      // A fracture connector has to be an edge shared by two faces
      if( numFluxElems == 2 )
      {
        typename SurfaceElementStencil::IndexContainerViewConstType const & sei = stencil.getElementIndices();

        // First index: face element. Second index: node
        for( localIndex kf = 0; kf < 2; ++kf )
        {
          // Set row DOF index
          globalIndex const rowIndex = presDofNumber[sei[iconn][1-kf]];

          // Get fracture, face and region/subregion/element indices (for elements on both sides)
          localIndex fractureIndex = sei[iconn][kf];

          // Get the number of nodes
          localIndex const numNodesPerFace = faceToNodeMap.sizeOfArray( elemsToFaces[fractureIndex][0] );

          // Loop over the two sides of each fracture element
          for( localIndex kf1 = 0; kf1 < 2; ++kf1 )
          {
            localIndex faceIndex = faceMap[fractureIndex][kf1];

            // Save the list of DOF associated with nodes
            for( localIndex a=0; a<numNodesPerFace; ++a )
            {
              for( localIndex i = 0; i < 3; ++i )
              {
                globalIndex const colIndex = dispDofNumber[faceToNodeMap( faceIndex, a )] + LvArray::integerConversion< globalIndex >( i );
                pattern.insertNonZero( rowIndex, colIndex );
              }
            }
          }
        }
      }
    } );
  } );
}

void LagrangianContactFlowSolver::
  assembleForceResidualDerivativeWrtPressure( DomainPartition & domain,
                                              DofManager const & dofManager,
                                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                              arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;
  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  FaceManager const & faceManager = mesh.getFaceManager();
  NodeManager & nodeManager = mesh.getNodeManager();
  ElementRegionManager const & elemManager = mesh.getElemManager();

  ArrayOfArraysView< localIndex const > const & faceToNodeMap = faceManager.nodeList().toViewConst();
  arrayView2d< real64 const > const & faceNormal = faceManager.faceNormal();

  string const & dispDofKey = dofManager.getKey( keys::TotalDisplacement );
  string const & presDofKey = dofManager.getKey( m_pressureKey );

  arrayView1d< globalIndex > const &
  dispDofNumber = nodeManager.getReference< globalIndex_array >( dispDofKey );
  globalIndex const rankOffset = dofManager.rankOffset();

  // Get the coordinates for all nodes
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition = nodeManager.referencePosition();

  elemManager.forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion const & subRegion )
  {
    if( subRegion.hasWrapper( m_pressureKey ) )
    {
      arrayView1d< globalIndex const > const &
      presDofNumber = subRegion.getReference< globalIndex_array >( presDofKey );
      arrayView1d< real64 const > const & pressure = subRegion.getReference< array1d< real64 > >( m_pressureKey );
      arrayView1d< real64 const > const & deltaPressure = subRegion.getReference< array1d< real64 > >( m_deltaPressureKey );
      arrayView2d< localIndex const > const & elemsToFaces = subRegion.faceList();

      forAll< serialPolicy >( subRegion.size(), [=]( localIndex const kfe )
      {
        localIndex const kf0 = elemsToFaces[kfe][0];
        localIndex const numNodesPerFace = faceToNodeMap.sizeOfArray( kf0 );

        array1d< real64 > Nbar( 3 );
        Nbar[ 0 ] = faceNormal[elemsToFaces[kfe][0]][0] - faceNormal[elemsToFaces[kfe][1]][0];
        Nbar[ 1 ] = faceNormal[elemsToFaces[kfe][0]][1] - faceNormal[elemsToFaces[kfe][1]][1];
        Nbar[ 2 ] = faceNormal[elemsToFaces[kfe][0]][2] - faceNormal[elemsToFaces[kfe][1]][2];
        LvArray::tensorOps::normalize< 3 >( Nbar );

        globalIndex rowDOF[12];
        real64 nodeRHS[12];
        stackArray1d< real64, 12 > dRdP( 3*numNodesPerFace );
        globalIndex colDOF[1];
        colDOF[0] = presDofNumber[kfe];

        for( localIndex kf=0; kf<2; ++kf )
        {
          localIndex const faceIndex = elemsToFaces[kfe][kf];

          for( localIndex a=0; a<numNodesPerFace; ++a )
          {
            // Compute local area contribution for each node
            array1d< real64 > nodalArea;
            m_contactSolver->computeFaceNodalArea( nodePosition, faceToNodeMap, elemsToFaces[kfe][kf], nodalArea );

            real64 const nodalForceMag = -( pressure[kfe] + deltaPressure[kfe] ) * nodalArea[a];
            array1d< real64 > globalNodalForce( 3 );
            LvArray::tensorOps::scaledCopy< 3 >( globalNodalForce, Nbar, nodalForceMag );

            for( localIndex i=0; i<3; ++i )
            {
              rowDOF[3*a+i] = dispDofNumber[faceToNodeMap( faceIndex, a )] + LvArray::integerConversion< globalIndex >( i );
              // Opposite sign w.r.t. theory because of minus sign in stiffness matrix definition (K < 0)
              nodeRHS[3*a+i] = +globalNodalForce[i] * pow( -1, kf );

              // Opposite sign w.r.t. theory because of minus sign in stiffness matrix definition (K < 0)
              dRdP( 3*a+i ) = -nodalArea[a] * Nbar[i] * pow( -1, kf );
            }
          }

          for( localIndex idof = 0; idof < numNodesPerFace * 3; ++idof )
          {
            localIndex const localRow = LvArray::integerConversion< localIndex >( rowDOF[idof] - rankOffset );

            if( localRow >= 0 && localRow < localMatrix.numRows() )
            {
              localMatrix.addToRow< parallelHostAtomic >( localRow,
                                                          colDOF,
                                                          &dRdP[idof],
                                                          1 );
              RAJA::atomicAdd( parallelHostAtomic{}, &localRhs[localRow], nodeRHS[idof] );
            }
          }
        }
      } );
    }
  } );
}

void LagrangianContactFlowSolver::
  assembleFluidMassResidualDerivativeWrtDisplacement( DomainPartition const & domain,
                                                      DofManager const & dofManager,
                                                      CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                      arrayView1d< real64 > const & GEOSX_UNUSED_PARAM( localRhs ) )
{
  GEOSX_MARK_FUNCTION;
  MeshLevel const & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  FaceManager const & faceManager = mesh.getFaceManager();
  NodeManager const & nodeManager = mesh.getNodeManager();

  arrayView2d< real64 const > const & faceNormal = faceManager.faceNormal();
  ArrayOfArraysView< localIndex const > const & faceToNodeMap = faceManager.nodeList().toViewConst();

  CRSMatrixView< real64 const, localIndex const > const &
  dFluxResidual_dAperture = getDerivativeFluxResidual_dAperture().toViewConst();

  string const & dispDofKey = dofManager.getKey( keys::TotalDisplacement );
  string const & presDofKey = dofManager.getKey( m_pressureKey );

  arrayView1d< globalIndex const > const &
  dispDofNumber = nodeManager.getReference< globalIndex_array >( dispDofKey );
  globalIndex const rankOffset = dofManager.rankOffset();

  // Get the coordinates for all nodes
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition = nodeManager.referencePosition();

  forTargetSubRegionsComplete< FaceElementSubRegion >( mesh,
                                                       [&]( localIndex const,
                                                            localIndex const,
                                                            localIndex const,
                                                            ElementRegionBase const & region,
                                                            FaceElementSubRegion const & subRegion )
  {
    if( subRegion.hasWrapper( m_pressureKey ) )
    {
      string const &
      fluidName = m_flowSolver->fluidModelNames()[m_flowSolver->targetRegionIndex( region.getName() )];
      SingleFluidBase const & fluid = getConstitutiveModel< SingleFluidBase >( subRegion, fluidName );
      arrayView2d< real64 const > const & density = fluid.density();

      arrayView1d< globalIndex const > const &
      presDofNumber = subRegion.getReference< array1d< globalIndex > >( presDofKey );
      arrayView2d< localIndex const > const & elemsToFaces = subRegion.faceList();
      arrayView1d< real64 const > const & area = subRegion.getElementArea().toViewConst();

      forAll< serialPolicy >( subRegion.size(), [&]( localIndex const kfe )
      {
        localIndex const numNodesPerFace = faceToNodeMap.sizeOfArray( elemsToFaces[kfe][0] );
        globalIndex nodeDOF[24];
        globalIndex elemDOF[1];
        elemDOF[0] = presDofNumber[kfe];

        array1d< real64 > Nbar( 3 );
        Nbar[ 0 ] = faceNormal[elemsToFaces[kfe][0]][0] - faceNormal[elemsToFaces[kfe][1]][0];
        Nbar[ 1 ] = faceNormal[elemsToFaces[kfe][0]][1] - faceNormal[elemsToFaces[kfe][1]][1];
        Nbar[ 2 ] = faceNormal[elemsToFaces[kfe][0]][2] - faceNormal[elemsToFaces[kfe][1]][2];
        LvArray::tensorOps::normalize< 3 >( Nbar );

        stackArray1d< real64, 2*3*4 > dRdU( 2*3*numNodesPerFace );

        // Accumulation derivative
        if( m_contactSolver->isElementInOpenState( subRegion, kfe ) )
        {
          for( localIndex kf=0; kf<2; ++kf )
          {
            // Compute local area contribution for each node
            array1d< real64 > nodalArea;
            m_contactSolver->computeFaceNodalArea( nodePosition, faceToNodeMap, elemsToFaces[kfe][kf], nodalArea );

            for( localIndex a=0; a<numNodesPerFace; ++a )
            {
              real64 const dAccumulationResidualdAperture = density[kfe][0] * nodalArea[a];
              for( localIndex i=0; i<3; ++i )
              {
                nodeDOF[ kf*3*numNodesPerFace + 3*a+i ] = dispDofNumber[faceToNodeMap( elemsToFaces[kfe][kf], a )]
                                                          + LvArray::integerConversion< globalIndex >( i );
                real64 const dAper_dU = -pow( -1, kf ) * Nbar[i];
                dRdU( kf*3*numNodesPerFace + 3*a+i ) = dAccumulationResidualdAperture * dAper_dU;
              }
            }
          }

          localIndex const localRow = LvArray::integerConversion< localIndex >( elemDOF[0] - rankOffset );

          localMatrix.addToRowBinarySearchUnsorted< serialAtomic >( localRow,
                                                                    nodeDOF,
                                                                    dRdU.data(),
                                                                    2 * 3 * numNodesPerFace );
        }

        // flux derivative
        bool skipAssembly = true;
        localIndex const numColumns = dFluxResidual_dAperture.numNonZeros( kfe );
        arraySlice1d< localIndex const > const & columns = dFluxResidual_dAperture.getColumns( kfe );
        arraySlice1d< real64 const > const & values = dFluxResidual_dAperture.getEntries( kfe );

        skipAssembly &= !( m_contactSolver->isElementInOpenState( subRegion, kfe ) );

        for( localIndex kfe1=0; kfe1<numColumns; ++kfe1 )
        {
          real64 const dR_dAper = values[kfe1];
          localIndex const kfe2 = columns[kfe1];

          skipAssembly &= !( m_contactSolver->isElementInOpenState( subRegion, kfe2 ) );

          for( localIndex kf=0; kf<2; ++kf )
          {
            // Compute local area contribution for each node
            array1d< real64 > nodalArea;
            m_contactSolver->computeFaceNodalArea( nodePosition, faceToNodeMap, elemsToFaces[kfe2][kf], nodalArea );

            for( localIndex a=0; a<numNodesPerFace; ++a )
            {
              for( localIndex i=0; i<3; ++i )
              {
                nodeDOF[ kf*3*numNodesPerFace + 3*a+i ] = dispDofNumber[faceToNodeMap( elemsToFaces[kfe2][kf], a )]
                                                          + LvArray::integerConversion< globalIndex >( i );
                real64 const dAper_dU = -pow( -1, kf ) * Nbar[i] * ( nodalArea[a] / area[kfe2] );
                dRdU( kf*3*numNodesPerFace + 3*a+i ) = dR_dAper * dAper_dU;
              }
            }
          }

          if( !skipAssembly )
          {
            localIndex const localRow = LvArray::integerConversion< localIndex >( elemDOF[0] - rankOffset );

            localMatrix.addToRowBinarySearchUnsorted< serialAtomic >( localRow,
                                                                      nodeDOF,
                                                                      dRdU.data(),
                                                                      2 * 3 * numNodesPerFace );
          }
        }
      } );
    }
  } );
}

void LagrangianContactFlowSolver::assembleStabilization( DomainPartition const & domain,
                                                         DofManager const & dofManager,
                                                         CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                         arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel const & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  FaceManager const & faceManager = mesh.getFaceManager();
  NodeManager const & nodeManager = mesh.getNodeManager();
  ElementRegionManager const & elemManager = mesh.getElemManager();

  string const & presDofKey = dofManager.getKey( m_pressureKey );
  globalIndex const rankOffset = dofManager.rankOffset();

  // Get the finite volume method used to compute the stabilization
  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & stabilizationMethod = fvManager.getFluxApproximation( m_stabilizationName );

  // Get the "face to element" map (valid for the entire mesh)
  FaceManager::ElemMapType const & faceToElem = faceManager.toElementRelation();
  arrayView2d< localIndex const > const & faceToElemRegion = faceToElem.m_toElementRegion.toViewConst();
  arrayView2d< localIndex const > const & faceToElemSubRegion = faceToElem.m_toElementSubRegion.toViewConst();
  arrayView2d< localIndex const > const & faceToElemIndex = faceToElem.m_toElementIndex.toViewConst();

  // Form the SurfaceGenerator, get the fracture name and use it to retrieve the faceMap (from fracture element to face)
  SurfaceGenerator const &
  surfaceGenerator = this->getParent().getGroup< SurfaceGenerator >( "SurfaceGen" );
  SurfaceElementRegion const & fractureRegion = elemManager.getRegion< SurfaceElementRegion >( surfaceGenerator.getFractureRegionName() );
  FaceElementSubRegion const & fractureSubRegion = fractureRegion.getSubRegion< FaceElementSubRegion >( "faceElementSubRegion" );
  GEOSX_ERROR_IF( !fractureSubRegion.hasWrapper( m_pressureKey ), "The fracture subregion must contain pressure field." );
  arrayView2d< localIndex const > const faceMap = fractureSubRegion.faceList();
  GEOSX_ERROR_IF( faceMap.size( 1 ) != 2, "A fracture face has to be shared by two cells." );

  // Get the pressures
  arrayView1d< real64 const > const &
  pressure = fractureSubRegion.getReference< array1d< real64 > >( m_pressureKey );
  arrayView1d< real64 const > const &
  deltaPressure = fractureSubRegion.getReference< array1d< real64 > >( m_deltaPressureKey );

  string const &
  fluidName = m_flowSolver->fluidModelNames()[m_flowSolver->targetRegionIndex( fractureRegion.getName() )];
  SingleFluidBase const & fluid = getConstitutiveModel< SingleFluidBase >( fractureSubRegion, fluidName );
  arrayView2d< real64 const > const & density = fluid.density();

  // Get the volume for all elements
  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > const elemVolume =
    elemManager.constructViewAccessor< real64_array, arrayView1d< real64 const > >( ElementSubRegionBase::viewKeyStruct::elementVolumeString() );

  // Get the coordinates for all nodes
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition = nodeManager.referencePosition();

  // Get area and rotation matrix for all faces
  ArrayOfArraysView< localIndex const > const & faceToNodeMap = faceManager.nodeList().toViewConst();
  arrayView1d< real64 const > const & faceArea = faceManager.faceArea();
  arrayView2d< real64 const > const & faceNormal = faceManager.faceNormal();

  // Bulk modulus accessor
  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > const bulkModulus =
    elemManager.constructMaterialViewAccessor< array1d< real64 >, arrayView1d< real64 const > >( ElasticIsotropic::viewKeyStruct::bulkModulusString(),
                                                                                                 m_contactSolver->getSolidSolver()->targetRegionNames(),
                                                                                                 m_contactSolver->getSolidSolver()->solidMaterialNames() );
  // Shear modulus accessor
  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > const shearModulus =
    elemManager.constructMaterialViewAccessor< array1d< real64 >, arrayView1d< real64 const > >( ElasticIsotropic::viewKeyStruct::shearModulusString(),
                                                                                                 m_contactSolver->getSolidSolver()->targetRegionNames(),
                                                                                                 m_contactSolver->getSolidSolver()->solidMaterialNames() );

  using NodeMapViewType = arrayView2d< localIndex const, cells::NODE_MAP_USD >;
  ElementRegionManager::ElementViewAccessor< NodeMapViewType > const elemToNode =
    elemManager.constructViewAccessor< CellBlock::NodeMapType, NodeMapViewType >( ElementSubRegionBase::viewKeyStruct::nodeListString() );
  ElementRegionManager::ElementViewConst< NodeMapViewType > const elemToNodeView = elemToNode.toNestedViewConst();

  arrayView1d< globalIndex const > const & presDofNumber = fractureSubRegion.getReference< globalIndex_array >( presDofKey );

  stabilizationMethod.forStencils< SurfaceElementStencil >( mesh, [&]( SurfaceElementStencil const & stencil )
  {
    typename SurfaceElementStencil::IndexContainerViewConstType const & sei = stencil.getElementIndices();

    forAll< serialPolicy >( stencil.size(), [=] ( localIndex const iconn )
    {
      localIndex const numFluxElems = sei.sizeOfArray( iconn );

      // A fracture connector has to be an edge shared by two faces
      if( numFluxElems == 2 )
      {
        // Find shared edge (pair of nodes)
        array1d< real64 > Nbar0( 3 ), Nbar1( 3 );
        Nbar0[ 0 ] = faceNormal[ faceMap[sei[iconn][0]][0] ][0] - faceNormal[ faceMap[sei[iconn][0]][1] ][0];
        Nbar0[ 1 ] = faceNormal[ faceMap[sei[iconn][0]][0] ][1] - faceNormal[ faceMap[sei[iconn][0]][1] ][1];
        Nbar0[ 2 ] = faceNormal[ faceMap[sei[iconn][0]][0] ][2] - faceNormal[ faceMap[sei[iconn][0]][1] ][2];
        LvArray::tensorOps::normalize< 3 >( Nbar0 );
        Nbar1[ 0 ] = faceNormal[ faceMap[sei[iconn][1]][0] ][0] - faceNormal[ faceMap[sei[iconn][1]][1] ][0];
        Nbar1[ 1 ] = faceNormal[ faceMap[sei[iconn][1]][0] ][1] - faceNormal[ faceMap[sei[iconn][1]][1] ][1];
        Nbar1[ 2 ] = faceNormal[ faceMap[sei[iconn][1]][0] ][2] - faceNormal[ faceMap[sei[iconn][1]][1] ][2];
        LvArray::tensorOps::normalize< 3 >( Nbar1 );

        real64 normalProduct = LvArray::tensorOps::AiBi< 3 >( Nbar0, Nbar1 );

        localIndex const id1 = ( normalProduct > 0.0 ) ? 0 : 1;

        localIndex const numNodesPerFace0 = faceToNodeMap.sizeOfArray( faceMap[sei[iconn][0]][0] );
        array1d< localIndex > nodes0( numNodesPerFace0 );
        for( localIndex i = 0; i < numNodesPerFace0; ++i )
        {
          nodes0[i] = faceToNodeMap( faceMap[sei[iconn][0]][0], i );
        }
        localIndex const numNodesPerFace1 = faceToNodeMap.sizeOfArray( faceMap[sei[iconn][1]][0] );
        array1d< localIndex > nodes1( numNodesPerFace1 );
        for( localIndex i = 0; i < numNodesPerFace1; ++i )
        {
          nodes1[i] = faceToNodeMap( faceMap[sei[iconn][1]][id1], i );
        }
        std::sort( nodes0.begin(), nodes0.end() );
        std::sort( nodes1.begin(), nodes1.end() );
        array1d< localIndex > edge( std::max( numNodesPerFace0, numNodesPerFace1 ) );
        edge.setValues< serialPolicy >( -1 );
        std::set_intersection( nodes0.begin(), nodes0.end(), nodes1.begin(), nodes1.end(), edge.begin() );
        localIndex realNodes = 0;
        for( localIndex i = 0; i < edge.size(); ++i )
        {
          if( edge[i] > -1 )
          {
            realNodes++;
          }
        }
        GEOSX_ERROR_IF( realNodes != 2, "An edge shared by two fracture elements must have 2 nodes." );
        edge.resize( realNodes );

        // Compute nodal area factor
        localIndex node0index0 = -1;
        localIndex node1index0 = -1;
        for( localIndex i = 0; i < numNodesPerFace0; ++i )
        {
          if( edge[0] == faceToNodeMap( faceMap[sei[iconn][0]][0], i ) )
          {
            node0index0 = i;
          }
          if( edge[1] == faceToNodeMap( faceMap[sei[iconn][0]][0], i ) )
          {
            node1index0 = i;
          }
        }
        localIndex node0index1 = -1;
        localIndex node1index1 = -1;
        for( localIndex i = 0; i < numNodesPerFace1; ++i )
        {
          if( edge[0] == faceToNodeMap( faceMap[sei[iconn][1]][id1], i ) )
          {
            node0index1 = i;
          }
          if( edge[1] == faceToNodeMap( faceMap[sei[iconn][1]][id1], i ) )
          {
            node1index1 = i;
          }
        }
        array1d< real64 > nodalArea0, nodalArea1;
        m_contactSolver->computeFaceNodalArea( nodePosition, faceToNodeMap, faceMap[sei[iconn][0]][0], nodalArea0 );
        m_contactSolver->computeFaceNodalArea( nodePosition, faceToNodeMap, faceMap[sei[iconn][1]][id1], nodalArea1 );
        real64 const areafac = nodalArea0[node0index0] * nodalArea1[node0index1] + nodalArea0[node1index0] * nodalArea1[node1index1];

        // first index: face, second index: element (T/B), third index: dof (x, y, z)
        real64 stiffDiagApprox[ 2 ][ 2 ][ 3 ];
        for( localIndex kf = 0; kf < 2; ++kf )
        {
          // Get fracture, face and region/subregion/element indices (for elements on both sides)
          localIndex const fractureIndex = sei[iconn][kf];

          for( localIndex i = 0; i < 2; ++i )
          {
            localIndex const faceIndex = ( kf == 0 || id1 == 0 ) ? faceMap[fractureIndex][i] : faceMap[fractureIndex][1-i];
            localIndex const ke = faceToElemIndex[faceIndex][0] >= 0 ? 0 : 1;

            localIndex const er  = faceToElemRegion[faceIndex][ke];
            localIndex const esr = faceToElemSubRegion[faceIndex][ke];
            localIndex const ei  = faceToElemIndex[faceIndex][ke];

            real64 const volume = elemVolume[er][esr][ei];

            // Get the "element to node" map for the specific region/subregion
            NodeMapViewType const & cellElemsToNodes = elemToNodeView[er][esr];
            localIndex const numNodesPerElem = cellElemsToNodes.size( 1 );

            // Compute the box size
            real64 maxSize[3];
            real64 minSize[3];
            for( localIndex j = 0; j < 3; ++j )
            {
              maxSize[j] = nodePosition[cellElemsToNodes[ei][0]][j];
              minSize[j] = nodePosition[cellElemsToNodes[ei][0]][j];
            }
            for( localIndex a = 1; a < numNodesPerElem; ++a )
            {
              for( localIndex j = 0; j < 3; ++j )
              {
                maxSize[j] = fmax( maxSize[j], nodePosition[cellElemsToNodes[ei][a]][j] );
                minSize[j] = fmin( minSize[j], nodePosition[cellElemsToNodes[ei][a]][j] );
              }
            }
            real64 boxSize[3];
            for( localIndex j = 0; j < 3; ++j )
            {
              boxSize[j] = maxSize[j] - minSize[j];
            }

            // Get linear elastic isotropic constitutive parameters for the element
            real64 const K = bulkModulus[er][esr][ei];
            real64 const G = shearModulus[er][esr][ei];
            real64 const E = 9.0 * K * G / ( 3.0 * K + G );
            real64 const nu = ( 3.0 * K - 2.0 * G ) / ( 2.0 * ( 3.0 * K + G ) );

            for( localIndex j = 0; j < 3; ++j )
            {
              stiffDiagApprox[ kf ][ i ][ j ] = E / ( ( 1.0 + nu )*( 1.0 - 2.0*nu ) ) * 2.0 / 9.0 * ( 2.0 - 3.0 * nu ) * volume / ( boxSize[j]*boxSize[j] );
            }
          }
        }
        real64 invTotStiffApprox[ 3 ][ 3 ] = { { 0 } };
        for( localIndex i = 0; i < 3; ++i )
        {
          // K(i,i)^-1 = Ka(i,i)^-1 + Kb(i,i)^-1
          // T -> top (index 0), B -> bottom (index 1)
          // Ka(i,i) = KT(i,i) + KB(i,i)
          // Kb(i,i) = KT(i,i) + KB(i,i)
          invTotStiffApprox[ i ][ i ] = 1.0 / ( stiffDiagApprox[ 0 ][ 0 ][ i ] + stiffDiagApprox[ 1 ][ 0 ][ i ] )
                                        + 1.0 / ( stiffDiagApprox[ 0 ][ 1 ][ i ] + stiffDiagApprox[ 1 ][ 1 ][ i ] );
        }

        array1d< real64 > avgNbar( 3 );

        // To be able to compute an average rotation matrix, normal has to point in the same direction.
        if( normalProduct < 0.0 )
        {
          LvArray::tensorOps::scale< 3 >( Nbar1, -1.0 );
          normalProduct *= -1.0;
        }
        // If the surfaces are co-planar, then use the first rotation matrix
        if( std::abs( normalProduct - 1.0 ) < 1.e+2*machinePrecision )
        {
          LvArray::tensorOps::copy< 3 >( avgNbar, Nbar0 );
        }
        // otherwise, compute the average rotation matrix
        else
        {
          avgNbar[ 0 ] = faceArea[faceMap[ sei[iconn][0] ][0]] * Nbar0[0] + faceArea[faceMap[ sei[iconn][1] ][0]] * Nbar1[0];
          avgNbar[ 1 ] = faceArea[faceMap[ sei[iconn][0] ][0]] * Nbar0[1] + faceArea[faceMap[ sei[iconn][1] ][0]] * Nbar1[1];
          avgNbar[ 2 ] = faceArea[faceMap[ sei[iconn][0] ][0]] * Nbar0[2] + faceArea[faceMap[ sei[iconn][1] ][0]] * Nbar1[2];
          LvArray::tensorOps::normalize< 3 >( avgNbar );
        }

        // Compute n^T * (invK) * n
        real64 temp[ 3 ];
        LvArray::tensorOps::Ri_eq_AijBj< 3, 3 >( temp, invTotStiffApprox, avgNbar );
        real64 const rotatedInvStiffApprox = LvArray::tensorOps::AiBi< 3 >( temp, avgNbar );

        // Add nodal area contribution
        stackArray1d< real64, 1 > totalInvStiffApproxDiag( 1 );
        totalInvStiffApproxDiag( 0 ) = rotatedInvStiffApprox * areafac;

        // Get DOF numbering
        localIndex fractureIndex[2];
        localIndex nDof[2];
        globalIndex elemDOF[2];
        for( localIndex kf = 0; kf < 2; ++kf )
        {
          fractureIndex[kf] = sei[iconn][kf];
          elemDOF[kf] = presDofNumber[fractureIndex[kf]];
          nDof[kf] = 1;
        }

        // Add mean density contribution
        totalInvStiffApproxDiag( 0 ) *= 0.5 * ( density[fractureIndex[0]][0] + density[fractureIndex[1]][0] );

        stackArray1d< real64, 1 > totalInvStiffApproxOffDiag( 1 );
        totalInvStiffApproxOffDiag( 0 ) = -totalInvStiffApproxDiag( 0 );

        // Compute rhs
        real64 rhs0 = 0.0;
        if( nDof[0] > 0 )
        {
          rhs0 -= totalInvStiffApproxDiag( 0 ) * ( -( pressure[fractureIndex[0]] + deltaPressure[fractureIndex[0]] ) );
        }
        if( nDof[1] > 0 )
        {
          rhs0 += totalInvStiffApproxDiag( 0 ) * ( -( pressure[fractureIndex[1]] + deltaPressure[fractureIndex[1]] ) );
        }
        real64 rhs1 = -rhs0;

        // Global matrix and rhs assembly
        for( localIndex kf = 0; kf < 2; ++kf )
        {
          localIndex const localRow = LvArray::integerConversion< localIndex >( elemDOF[kf] - rankOffset );

          real64 const & rhs = ( kf == 0 ) ? rhs0 : rhs1;

          // Only assemble contribution if "row" fracture element is local
          // TODO: use parallel atomics
          if( localRow >= 0 && localRow < localMatrix.numRows() )
          {
            // (i,i)-block
            localMatrix.addToRowBinarySearchUnsorted< parallelHostAtomic >( localRow,
                                                                            &elemDOF[kf],
                                                                            totalInvStiffApproxDiag.data(),
                                                                            nDof[kf] );
            // (i,j)-block
            if( nDof[1-kf] > 0 )
            {
              localMatrix.addToRowBinarySearchUnsorted< parallelHostAtomic >( localRow,
                                                                              &elemDOF[1 - kf],
                                                                              totalInvStiffApproxOffDiag.data(),
                                                                              nDof[1 - kf] );
            }

            // residual
            RAJA::atomicAdd( parallelHostAtomic{}, &localRhs[localRow], rhs );
          }
        }
      }
    } );
  } );
}

void LagrangianContactFlowSolver::applySystemSolution( DofManager const & dofManager,
                                                       arrayView1d< real64 const > const & localSolution,
                                                       real64 const scalingFactor,
                                                       DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  m_contactSolver->applySystemSolution( dofManager,
                                        localSolution,
                                        scalingFactor,
                                        domain );

  m_flowSolver->applySystemSolution( dofManager,
                                     localSolution,
                                     -scalingFactor,
                                     domain );

  updateOpeningForFlow( domain );
}

void LagrangianContactFlowSolver::solveSystem( DofManager const & dofManager,
                                               ParallelMatrix & matrix,
                                               ParallelVector & rhs,
                                               ParallelVector & solution )
{
  GEOSX_MARK_FUNCTION;

  //rhs.write( "rhs.petsc", LAIOutputFormat::MATRIX_MARKET );
  //matrix.write( "matrix.petsc", LAIOutputFormat::NATIVE_BINARY );

  if( getLogLevel() > 3 )
  {
    matrix.write( "matrix.mtx", LAIOutputFormat::MATRIX_MARKET );
    rhs.write( "rhs.mtx", LAIOutputFormat::MATRIX_MARKET );
  }

  SolverBase::solveSystem( dofManager, matrix, rhs, solution );

  //solution.write( "sol.petsc", LAIOutputFormat::MATRIX_MARKET );

  if( getLogLevel() > 3 )
  {
    solution.write( "sol.mtx", LAIOutputFormat::MATRIX_MARKET );
  }

  // int rank = MpiWrapper::commRank( MPI_COMM_GEOSX );
  // if( rank == 0 )
  // {
  //   string str;
  //   std::getline( std::cin, str );
  //   if( str.length() > 0 )
  //   {
  //     GEOSX_ERROR( "STOP" );
  //   }
  // }
  // MpiWrapper::barrier( MPI_COMM_GEOSX );
}

void LagrangianContactFlowSolver::setNextDt( real64 const & currentDt,
                                             real64 & nextDt )
{
  nextDt = currentDt;
}


void LagrangianContactFlowSolver::setUpDflux_dApertureMatrix( DomainPartition & domain,
                                                              DofManager const & dofManager,
                                                              CRSMatrix< real64, globalIndex > & localMatrix )
{
  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  std::unique_ptr< CRSMatrix< real64, localIndex > > &
  derivativeFluxResidual_dAperture = this->getRefDerivativeFluxResidual_dAperture();

  {
    localIndex numRows = 0;
    this->template forTargetSubRegions< FaceElementSubRegion >( mesh, [&]( localIndex const,
                                                                           auto const & elementSubRegion )
    {
      numRows += elementSubRegion.size();
    } );

    derivativeFluxResidual_dAperture = std::make_unique< CRSMatrix< real64, localIndex > >( numRows, numRows );
    derivativeFluxResidual_dAperture->setName( this->getName() + "/derivativeFluxResidual_dAperture" );

    derivativeFluxResidual_dAperture->reserveNonZeros( localMatrix.numNonZeros() );
    localIndex maxRowSize = -1;
    for( localIndex row = 0; row < localMatrix.numRows(); ++row )
    {
      localIndex const rowSize = localMatrix.numNonZeros( row );
      maxRowSize = maxRowSize > rowSize ? maxRowSize : rowSize;
    }
    // TODO This is way too much. The With the full system rowSize is not a good estimate for this.
    for( localIndex row = 0; row < numRows; ++row )
    {
      derivativeFluxResidual_dAperture->reserveNonZeros( row, maxRowSize );
    }
  }

  string const presDofKey = dofManager.getKey( extrinsicMeshData::flow::pressure::key() );

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_flowSolver->getDiscretizationName() );

  fluxApprox.forStencils< SurfaceElementStencil >( mesh, [&]( SurfaceElementStencil const & stencil )
  {
    for( localIndex iconn = 0; iconn < stencil.size(); ++iconn )
    {
      localIndex const numFluxElems = stencil.stencilSize( iconn );
      typename SurfaceElementStencil::IndexContainerViewConstType const & sei = stencil.getElementIndices();

      for( localIndex k0 = 0; k0 < numFluxElems; ++k0 )
      {
        for( localIndex k1 = 0; k1 < numFluxElems; ++k1 )
        {
          derivativeFluxResidual_dAperture->insertNonZero( sei[iconn][k0], sei[iconn][k1], 0.0 );
        }
      }
    }
  } );
}

REGISTER_CATALOG_ENTRY( SolverBase, LagrangianContactFlowSolver, string const &, Group * const )

} /* namespace geosx */
