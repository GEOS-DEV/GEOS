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
#include "constitutive/contact/ContactRelationBase.hpp"
#include "constitutive/fluid/SingleFluidBase.hpp"
#include "finiteElement/Kinematics.h"
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "managers/DomainPartition.hpp"
#include "managers/FieldSpecification/FieldSpecificationManager.hpp"
#include "managers/NumericalMethodsManager.hpp"
#include "mesh/FaceElementRegion.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "meshUtilities/ComputationalGeometry.hpp"
#include "mpiCommunications/NeighborCommunicator.hpp"
#include "physicsSolvers/multiphysics/LagrangianContactSolver.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBase.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"
#include "physicsSolvers/surfaceGeneration/SurfaceGenerator.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"
#include "linearAlgebra/utilities/LAIHelperFunctions.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;

LagrangianContactFlowSolver::LagrangianContactFlowSolver( const std::string & name,
                                                          Group * const parent ):
  SolverBase( name, parent ),
  m_contactSolverName(),
  m_flowSolverName(),
  m_contactSolver( nullptr ),
  m_flowSolver( nullptr ),
  m_stabilizationName( "" ),
  m_defaultConductivity( -1.0 )
{
  registerWrapper( viewKeyStruct::contactSolverNameString, &m_contactSolverName, false )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Name of the contact solid mechanics solver to use in the contact with flow solver" );

  registerWrapper( viewKeyStruct::flowSolverNameString, &m_flowSolverName, false )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Name of the flow mechanics solver to use in the contact with flow solver" );

  registerWrapper( viewKeyStruct::stabilizationNameString, &m_stabilizationName, false )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Name of the stabilization to use in the lagrangian contact with flow solver" );

  registerWrapper( viewKeyStruct::defaultConductivityString, &m_defaultConductivity, false )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Value of the default conductivity C_{f,0} for the fracture" );
}

void LagrangianContactFlowSolver::RegisterDataOnMesh( dataRepository::Group * const MeshBodies )
{
  m_contactSolver->RegisterDataOnMesh( MeshBodies );
  m_flowSolver->RegisterDataOnMesh( MeshBodies );
}

void LagrangianContactFlowSolver::InitializePreSubGroups( Group * const rootGroup )
{
  SolverBase::InitializePreSubGroups( rootGroup );
}

void LagrangianContactFlowSolver::ImplicitStepSetup( real64 const & time_n,
                                                     real64 const & dt,
                                                     DomainPartition * const domain,
                                                     DofManager & dofManager,
                                                     ParallelMatrix & matrix,
                                                     ParallelVector & rhs,
                                                     ParallelVector & solution )
{
  m_contactSolver->ImplicitStepSetup( time_n, dt, domain,
                                      dofManager,
                                      matrix,
                                      rhs,
                                      solution );

  m_flowSolver->ImplicitStepSetup( time_n, dt, domain,
                                   dofManager,
                                   matrix,
                                   rhs,
                                   solution );
}

void LagrangianContactFlowSolver::ImplicitStepComplete( real64 const & time_n,
                                                        real64 const & dt,
                                                        DomainPartition * const domain )
{
  m_contactSolver->ImplicitStepComplete( time_n, dt, domain );
  m_flowSolver->ImplicitStepComplete( time_n, dt, domain );
}

void LagrangianContactFlowSolver::PostProcessInput()
{
  m_contactSolver = this->getParent()->GetGroup< LagrangianContactSolver >( m_contactSolverName );
  GEOSX_ERROR_IF( m_contactSolver == nullptr, this->getName() << ": invalid contact solver name: " << m_contactSolverName );

  m_flowSolver = this->getParent()->GetGroup< FlowSolverBase >( m_flowSolverName );
  GEOSX_ERROR_IF( m_flowSolver == nullptr, this->getName() << ": invalid flow solver name: " << m_flowSolverName );

  SolverBase::PostProcessInput();
}

void LagrangianContactFlowSolver::InitializePostInitialConditions_PreSubGroups( Group * const GEOSX_UNUSED_PARAM( problemManager ) )
{}

LagrangianContactFlowSolver::~LagrangianContactFlowSolver()
{
  // TODO Auto-generated destructor stub
}

void LagrangianContactFlowSolver::ResetStateToBeginningOfStep( DomainPartition * const domain )
{
  m_contactSolver->ResetStateToBeginningOfStep( domain );
  m_flowSolver->ResetStateToBeginningOfStep( domain );
  UpdateOpeningForFlow( domain );
}

real64 LagrangianContactFlowSolver::SolverStep( real64 const & time_n,
                                                real64 const & dt,
                                                int const cycleNumber,
                                                DomainPartition * const domain )
{
  {
    // FIXME: to put somewhere else ...
    MeshLevel * const meshLevel = domain->getMeshBody( 0 )->getMeshLevel( 0 );
    ElementRegionManager * const elemManager = meshLevel->getElemManager();
    elemManager->forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion & subRegion )
    {
      if( subRegion.hasWrapper( m_pressureKey ) )
      {
        subRegion.getElementDefaultConductivity() = m_defaultConductivity;
      }
    } );
  }

  real64 dtReturn = dt;

  ImplicitStepSetup( time_n,
                     dt,
                     domain,
                     m_dofManager,
                     m_matrix,
                     m_rhs,
                     m_solution );

  SetupSystem( domain,
               m_dofManager,
               m_matrix,
               m_rhs,
               m_solution );

  // currently the only method is implicit time integration
  dtReturn = this->NonlinearImplicitStep( time_n,
                                          dt,
                                          cycleNumber,
                                          domain,
                                          m_dofManager,
                                          m_matrix,
                                          m_rhs,
                                          m_solution );

  m_contactSolver->getSolidSolver()->updateStress( domain );

  // final step for completion of timestep. typically secondary variable updates and cleanup.
  ImplicitStepComplete( time_n, dtReturn, domain );

  return dtReturn;
}

void LagrangianContactFlowSolver::UpdateOpeningForFlow( DomainPartition * const domain )
{
  MeshLevel * const meshLevel = domain->getMeshBody( 0 )->getMeshLevel( 0 );
  ElementRegionManager * const elemManager = meshLevel->getElemManager();

  elemManager->forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion & subRegion )
  {
    if( subRegion.hasWrapper( m_pressureKey ) )
    {
      arrayView1d< integer const > const & ghostRank = subRegion.ghostRank();
      arrayView2d< real64 const > const & localJump = subRegion.getReference< array2d< real64 > >( m_localJumpKey );
      arrayView1d< real64 const > const & area = subRegion.getElementArea();
      arrayView1d< real64 const > const & volume = subRegion.getElementVolume();
      arrayView1d< real64 > const &
      aperture = subRegion.getReference< array1d< real64 > >( FlowSolverBase::viewKeyStruct::effectiveApertureString );
      arrayView1d< real64 > const &
      deltaVolume = subRegion.getReference< array1d< real64 > >( FlowSolverBase::viewKeyStruct::deltaVolumeString );

      forAll< serialPolicy >( subRegion.size(), [&]( localIndex const kfe )
      {
        if( ghostRank[kfe] < 0 )
        {
          if( m_contactSolver->IsElementInOpenState( subRegion, kfe ) )
          {
            aperture[kfe] = localJump[kfe][0];
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

real64 LagrangianContactFlowSolver::SplitOperatorStep( real64 const & GEOSX_UNUSED_PARAM( time_n ),
                                                       real64 const & dt,
                                                       integer const GEOSX_UNUSED_PARAM( cycleNumber ),
                                                       DomainPartition * const GEOSX_UNUSED_PARAM( domain ) )
{
  real64 dtReturn = dt;
  return dtReturn;
}

real64 LagrangianContactFlowSolver::ExplicitStep( real64 const & GEOSX_UNUSED_PARAM( time_n ),
                                                  real64 const & dt,
                                                  const int GEOSX_UNUSED_PARAM( cycleNumber ),
                                                  DomainPartition * const GEOSX_UNUSED_PARAM( domain ) )
{
  GEOSX_MARK_FUNCTION;
  GEOSX_ERROR( "ExplicitStep non available for LagrangianContactFlowSolver!" );
  return dt;
}

real64 LagrangianContactFlowSolver::NonlinearImplicitStep( real64 const & time_n,
                                                           real64 const & dt,
                                                           integer const cycleNumber,
                                                           DomainPartition * const domain,
                                                           DofManager const & dofManager,
                                                           ParallelMatrix & matrix,
                                                           ParallelVector & rhs,
                                                           ParallelVector & solution )
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

  bool useElasticStep = !m_contactSolver->IsFractureAllInStickCondition( domain );

  // outer loop attempts to apply full timestep, and managed the cutting of the timestep if
  // required.
  for( dtAttempt = 0; dtAttempt < maxNumberDtCuts; ++dtAttempt )
  {
    // reset the solver state, since we are restarting the time step
    if( dtAttempt > 0 )
    {
      ResetStateToBeginningOfStep( domain );
      globalIndex numStick, numSlip, numOpen;
      m_contactSolver->ComputeFractureStateStatistics( domain, numStick, numSlip, numOpen, true );
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

        // call assemble to fill the matrix and the rhs
        matrix.zero();
        rhs.zero();
        AssembleSystem( time_n, stepDt, domain, dofManager, matrix, rhs );

        // apply boundary conditions to system
        ApplyBoundaryConditions( time_n, stepDt, domain, dofManager, matrix, rhs );

        // TODO: maybe add scale function here?
        // Scale()

        real64 residualNorm;
        // get residual norm
        if( computeResidual )
        {
          residualNorm = CalculateResidualNorm( domain, dofManager, rhs );
        }
        else
        {
          residualNorm = lastResidual;
        }

        if( m_linearSolverParameters.solverType != "direct" && getLogLevel() >= 1 && logger::internal::rank==0 )
        {
          if( newtonIter!=0 )
          {
            char output[46] = {0};
            sprintf( output,
                     "Last LinSolve(iter,tol) = (%4d, %4.2e) ; ",
                     m_systemSolverParameters.m_numKrylovIter,
                     m_systemSolverParameters.m_krylovTol );
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
        // TODO: need to combine overlapping usage on LinearSolverParameters and SystemSolverParamters
        if( m_systemSolverParameters.useAdaptiveKrylovTol())
        {
          m_systemSolverParameters.m_krylovTol = LinearSolverParameters::eisenstatWalker( residualNorm, lastResidual );
        }

        // call the default linear solver on the system
        SolveSystem( dofManager, matrix, rhs, solution );

        scaleFactor = ScalingForSystemSolution( domain, dofManager, solution );

        // do line search in case residual has increased
        if( m_nonlinearSolverParameters.m_lineSearchAction>0 && newtonIter > 0 )
        {
          bool lineSearchSuccess = LineSearch( time_n, stepDt, cycleNumber, domain, dofManager,
                                               matrix, rhs, solution, scaleFactor, residualNorm );

          if( !lineSearchSuccess )
          {
            if( m_nonlinearSolverParameters.m_lineSearchAction==1 )
            {
              GEOSX_LOG_LEVEL_RANK_0( 1, "        Line search failed to produce reduced residual. Accepting iteration." );
            }
            else if( m_nonlinearSolverParameters.m_lineSearchAction==2 )
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
          ApplySystemSolution( dofManager, solution, scaleFactor, domain );
          // Need to compute the residual norm
          computeResidual = true;
        }

        if( !CheckSystemSolution( domain, dofManager, solution, scaleFactor ) )
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
      bool const isPreviousFractureStateValid = m_contactSolver->UpdateFractureState( domain );
      GEOSX_LOG_LEVEL_RANK_0( 1, "active set flag: " << std::boolalpha << isPreviousFractureStateValid );

      if( getLogLevel() > 2 )
      {
        globalIndex numStick, numSlip, numOpen;
        m_contactSolver->ComputeFractureStateStatistics( domain, numStick, numSlip, numOpen, true );
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
        ResetStateToBeginningOfStep( domain );
        m_contactSolver->SetFractureStateForElasticStep( domain );
        UpdateOpeningForFlow( domain );
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
    GEOSX_ERROR( "Active set did not reached a solution. Terminating..." );
  }
  else
  {
    GEOSX_LOG_RANK_0( "Number of active set iterations: " << m_activeSetIter );
  }

  // return the achieved timestep
  return stepDt;
}

bool LagrangianContactFlowSolver::LineSearch( real64 const & time_n,
                                              real64 const & dt,
                                              integer const GEOSX_UNUSED_PARAM( cycleNumber ),
                                              DomainPartition * const domain,
                                              DofManager const & dofManager,
                                              ParallelMatrix & matrix,
                                              ParallelVector & rhs,
                                              ParallelVector const & solution,
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

  ApplySystemSolution( dofManager, solution, scaleFactor, domain );

  // re-assemble system
  matrix.zero();
  rhs.zero();
  AssembleSystem( time_n, dt, domain, dofManager, matrix, rhs );

  // apply boundary conditions to system
  ApplyBoundaryConditions( time_n, dt, domain, dofManager, matrix, rhs );

  // get residual norm
  real64 residualNormT = CalculateResidualNorm( domain, dofManager, rhs );

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
      localScaleFactor = m_contactSolver->ParabolicInterpolationThreePoints( lamc, lamm, ff0, ffT, ffm );
    }

    // Update x; keep the books on lambda
    real64 const deltaLocalScaleFactor = ( localScaleFactor - previousLocalScaleFactor );
    cumulativeScale += deltaLocalScaleFactor;

    if( !CheckSystemSolution( domain, dofManager, solution, deltaLocalScaleFactor ) )
    {
      GEOSX_LOG_LEVEL_RANK_0( 1, "        Line search " << lineSearchIteration << ", solution check failed" );
      continue;
    }

    ApplySystemSolution( dofManager, solution, deltaLocalScaleFactor, domain );
    lamm = lamc;
    lamc = localScaleFactor;

    // Keep the books on the function norms
    // re-assemble system
    // TODO: add a flag to avoid a completely useless Jacobian computation: rhs is enough
    matrix.zero();
    rhs.zero();
    AssembleSystem( time_n, dt, domain, dofManager, matrix, rhs );

    // apply boundary conditions to system
    ApplyBoundaryConditions( time_n, dt, domain, dofManager, matrix, rhs );

    if( getLogLevel() >= 1 && logger::internal::rank==0 )
    {
      char output[100];
      sprintf( output, "        Line search @ %0.3f:      ", cumulativeScale );
      std::cout<<output;
    }

    // get residual norm
    residualNormT = CalculateResidualNorm( domain, dofManager, rhs );
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

void LagrangianContactFlowSolver::SetupDofs( DomainPartition const * const domain,
                                             DofManager & dofManager ) const
{
  GEOSX_MARK_FUNCTION;
  m_contactSolver->SetupDofs( domain, dofManager );
  m_flowSolver->SetupDofs( domain, dofManager );

  // restrict coupling to fracture regions only
  ElementRegionManager const * const elemManager = domain->getMeshBody( 0 )->getMeshLevel( 0 )->getElemManager();
  string_array fractureRegions;
  elemManager->forElementRegions< FaceElementRegion >( [&]( FaceElementRegion const & elementRegion )
  {
    fractureRegions.push_back( elementRegion.getName() );
  } );

  dofManager.addCoupling( keys::TotalDisplacement,
                          m_pressureKey,
                          DofManager::Connector::Elem,
                          fractureRegions );
}

void LagrangianContactFlowSolver::SetupSystem( DomainPartition * const domain,
                                               DofManager & dofManager,
                                               ParallelMatrix & matrix,
                                               ParallelVector & rhs,
                                               ParallelVector & solution )
{
  GEOSX_MARK_FUNCTION;
  m_flowSolver->ResetViews( domain );

  // setup DofManager
  dofManager.setMesh( domain, 0, 0 );

  // add traction and coupling
  SetupDofs( domain, dofManager );
  dofManager.reorderByRank();

  localIndex const numDisplacementDofs = dofManager.numLocalDofs( keys::TotalDisplacement );
  localIndex const numTractionDofs = dofManager.numLocalDofs( m_tractionKey );
  localIndex const numPressureDofs = dofManager.numLocalDofs( m_pressureKey );
  GEOSX_LOG_RANK( numDisplacementDofs << " " << numTractionDofs << " " << numPressureDofs );

  matrix.createWithLocalSize( numDisplacementDofs + numTractionDofs + numPressureDofs,
                              numDisplacementDofs + numTractionDofs + numPressureDofs,
                              5*(3*27+3*12+12),
                              MPI_COMM_GEOSX );
  rhs.createWithLocalSize( numDisplacementDofs + numTractionDofs + numPressureDofs,
                           MPI_COMM_GEOSX );
  solution.createWithLocalSize( numDisplacementDofs + numTractionDofs + numPressureDofs,
                                MPI_COMM_GEOSX );

  dofManager.setSparsityPattern( matrix, false );

  AddTransmissibilityDerivativePattern( domain,
                                        dofManager,
                                        &matrix );
  matrix.close();
}

void LagrangianContactFlowSolver::AssembleSystem( real64 const time,
                                                  real64 const dt,
                                                  DomainPartition * const domain,
                                                  DofManager const & dofManager,
                                                  ParallelMatrix & matrix,
                                                  ParallelVector & rhs )
{
  GEOSX_MARK_FUNCTION;

  // Need to synchronize the two iteration counters
  m_contactSolver->getNonlinearSolverParameters().m_numNewtonIterations = m_nonlinearSolverParameters.m_numNewtonIterations;
  m_contactSolver->AssembleSystem( time,
                                   dt,
                                   domain,
                                   dofManager,
                                   matrix,
                                   rhs );

  UpdateOpeningForFlow( domain );
  m_flowSolver->AssembleSystem( time,
                                dt,
                                domain,
                                dofManager,
                                matrix,
                                rhs );

  AssembleForceResidualDerivativeWrtPressure( domain, dofManager, &matrix, &rhs );
  AssembleFluidMassResidualDerivativeWrtDisplacement( domain, dofManager, &matrix, &rhs );
  AssembleStabiliziation( domain, dofManager, &matrix, &rhs );
}

void LagrangianContactFlowSolver::ApplyBoundaryConditions( real64 const time,
                                                           real64 const dt,
                                                           DomainPartition * const domain,
                                                           DofManager const & dofManager,
                                                           ParallelMatrix & matrix,
                                                           ParallelVector & rhs )
{
  GEOSX_MARK_FUNCTION;
  m_contactSolver->ApplyBoundaryConditions( time,
                                            dt,
                                            domain,
                                            dofManager,
                                            matrix,
                                            rhs );
  m_flowSolver->ApplyBoundaryConditions( time,
                                         dt,
                                         domain,
                                         dofManager,
                                         matrix,
                                         rhs );
}

real64 LagrangianContactFlowSolver::CalculateResidualNorm( DomainPartition const * const GEOSX_UNUSED_PARAM( domain ),
                                                           DofManager const & dofManager,
                                                           ParallelVector const & rhs )
{
  GEOSX_MARK_FUNCTION;

  localIndex numDispDofs = dofManager.numLocalDofs( keys::TotalDisplacement );
  localIndex numTracDofs = dofManager.numLocalDofs( m_tractionKey );
  localIndex numPresDofs = dofManager.numLocalDofs( m_pressureKey );
  real64 const * localResidual = rhs.extractLocalVector();
  real64 localResidualNorm[4] = {0.0, 0.0, 0.0, 0.0};
  for( localIndex i=0; i<numDispDofs; ++i )
  {
    localResidualNorm[0] += localResidual[i] * localResidual[i];
  }
  for( localIndex i=numDispDofs; i<numDispDofs+numTracDofs; ++i )
  {
    localResidualNorm[1] += localResidual[i] * localResidual[i];
  }
  for( localIndex i=numDispDofs+numTracDofs; i<numDispDofs+numTracDofs+numPresDofs; ++i )
  {
    localResidualNorm[2] += localResidual[i] * localResidual[i];
  }
  localResidualNorm[3] = localResidualNorm[0] + localResidualNorm[1] + localResidualNorm[2];

  real64 globalResidualNorm[4] = {0.0, 0.0, 0.0, 0.0};
  int const rank = MpiWrapper::Comm_rank( MPI_COMM_GEOSX );
  int const size = MpiWrapper::Comm_size( MPI_COMM_GEOSX );
  real64_array globalValues( 4*size );
  globalValues = 0;

  // Everything is done on rank 0
  MpiWrapper::gather( localResidualNorm,
                      4,
                      globalValues.data(),
                      4,
                      0,
                      MPI_COMM_GEOSX );

  if( rank==0 )
  {
    for( int r=0; r<size; ++r )
    {
      // sum across all ranks
      globalResidualNorm[0] += globalValues[4*r];
      globalResidualNorm[1] += globalValues[4*r+1];
      globalResidualNorm[2] += globalValues[4*r+2];
      globalResidualNorm[3] += globalValues[4*r+3];
    }
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
    // Add 0 just to match Matlab code results
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

void LagrangianContactFlowSolver::AddTransmissibilityDerivativePattern( DomainPartition * const domain,
                                                                        DofManager const & dofManager,
                                                                        ParallelMatrix * const matrix )
{
  MeshLevel const * const mesh = domain->getMeshBody( 0 )->getMeshLevel( 0 );

  FaceManager const * const faceManager = mesh->getFaceManager();
  NodeManager const * const nodeManager = mesh->getNodeManager();
  ElementRegionManager const * const elemManager = mesh->getElemManager();

  string const dispDofKey = dofManager.getKey( keys::TotalDisplacement );
  string const presDofKey = dofManager.getKey( m_pressureKey );

  arrayView1d< globalIndex const > const &
  dispDofNumber = nodeManager->getReference< globalIndex_array >( dispDofKey );
  ArrayOfArraysView< localIndex const > const & faceToNodeMap = faceManager->nodeList();

  // Get the finite volume method used to compute the stabilization
  NumericalMethodsManager const * const numericalMethodManager = domain->getParent()->GetGroup< NumericalMethodsManager >( keys::numericalMethodsManager );
  FiniteVolumeManager const * const fvManager = numericalMethodManager->GetGroup< FiniteVolumeManager >( keys::finiteVolumeManager );
  FluxApproximationBase const * const fluxApprox = fvManager->getFluxApproximation( m_stabilizationName );

  // Form the SurfaceGenerator, get the fracture name and use it to retrieve the faceMap (from fracture element to face)
  SurfaceGenerator const * const
  surfaceGenerator = this->getParent()->GetGroup< SolverBase >( "SurfaceGen" )->group_cast< SurfaceGenerator const * >();
  FaceElementRegion const * const fractureRegion = elemManager->GetRegion< FaceElementRegion >( surfaceGenerator->getFractureRegionName() );
  FaceElementSubRegion const * const fractureSubRegion = fractureRegion->GetSubRegion< FaceElementSubRegion >( "default" );
  GEOSX_ERROR_IF( !fractureSubRegion->hasWrapper( m_pressureKey ), "The fracture subregion must contain pressure field." );
  FaceElementSubRegion::FaceMapType const & faceMap = fractureSubRegion->faceList();
  GEOSX_ERROR_IF( faceMap.size( 1 ) != 2, "A fracture face has to be shared by two cells." );

  arrayView1d< globalIndex const > const &
  presDofNumber = fractureSubRegion->getReference< globalIndex_array >( presDofKey );

  arrayView2d< localIndex const > const & elemsToFaces = fractureSubRegion->faceList();

  fluxApprox->forStencils< FaceElementStencil >( [&]( FaceElementStencil const & stencil )
  {
    // Get ghost rank
    ArrayOfArraysView< integer const > const & isGhostConnectors = stencil.getIsGhostConnectors();

    for( localIndex iconn=0; iconn<stencil.size(); ++iconn )
    {
      localIndex const numFluxElems = stencil.stencilSize( iconn );

      // A fracture connector has to be an edge shared by two faces
      if( numFluxElems == 2 && isGhostConnectors[iconn][0] < 0 )
      {
        typename FaceElementStencil::IndexContainerViewConstType const & sei = stencil.getElementIndices();

        // First index: face element. Second index: node
        globalIndex rowDOF[1];
        globalIndex colDOF[24];
        real64_array ones( 24 );
        ones = 1.0;
        for( localIndex kf = 0; kf < 2; ++kf )
        {
          // Set row DOF index
          rowDOF[0] = presDofNumber[sei[iconn][1-kf]];

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
                colDOF[kf1*3*numNodesPerFace + 3*a+i] = dispDofNumber[faceToNodeMap( faceIndex, a )] + integer_conversion< globalIndex >( i );
              }
            }
          }

          matrix->insert( rowDOF,
                          colDOF,
                          ones.data(),
                          1,
                          2*3*numNodesPerFace );
        }
      }
    }
  } );
}

void LagrangianContactFlowSolver::AssembleForceResidualDerivativeWrtPressure( DomainPartition * const domain,
                                                                              DofManager const & dofManager,
                                                                              ParallelMatrix * const matrix,
                                                                              ParallelVector * const rhs )
{
  GEOSX_MARK_FUNCTION;
  MeshLevel * const mesh = domain->getMeshBody( 0 )->getMeshLevel( 0 );

  FaceManager const * const faceManager = mesh->getFaceManager();
  NodeManager * const nodeManager = mesh->getNodeManager();
  ElementRegionManager * const elemManager = mesh->getElemManager();

  arrayView1d< R1Tensor const > const & faceNormal = faceManager->faceNormal();
  ArrayOfArraysView< localIndex const > const & faceToNodeMap = faceManager->nodeList();

  arrayView1d< R1Tensor > const &
  fext = nodeManager->getReference< array1d< R1Tensor > >( SolidMechanicsLagrangianFEM::viewKeyStruct::forceExternal );

  string const dispDofKey = dofManager.getKey( keys::TotalDisplacement );
  string const presDofKey = dofManager.getKey( m_pressureKey );

  arrayView1d< globalIndex > const &
  dispDofNumber = nodeManager->getReference< globalIndex_array >( dispDofKey );

  matrix->open();
  rhs->open();

  elemManager->forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion const & subRegion )
  {
    if( subRegion.hasWrapper( m_pressureKey ) )
    {
      arrayView1d< globalIndex const > const &
      presDofNumber = subRegion.getReference< globalIndex_array >( presDofKey );
      arrayView1d< real64 const > const & pressure = subRegion.getReference< array1d< real64 > >( m_pressureKey );
      arrayView1d< real64 const > const & deltaPressure = subRegion.getReference< array1d< real64 > >( m_deltaPressureKey );
      arrayView1d< integer const > const & ghostRank = subRegion.ghostRank();
      arrayView1d< real64 const > const & area = subRegion.getElementArea();
      arrayView2d< localIndex const > const & elemsToFaces = subRegion.faceList();

      forAll< serialPolicy >( subRegion.size(), [=]( localIndex const kfe )
      {
        if( ghostRank[kfe] < 0 )
        {
          localIndex const kf0 = elemsToFaces[kfe][0];
          localIndex const numNodesPerFace = faceToNodeMap.sizeOfArray( kf0 );

          R1Tensor Nbar = faceNormal[elemsToFaces[kfe][0]];
          Nbar -= faceNormal[elemsToFaces[kfe][1]];
          Nbar.Normalize();

          globalIndex rowDOF[12];
          real64 nodeRHS[12];
          stackArray1d< real64, 12 > dRdP( 3*numNodesPerFace );
          dRdP = 0.0;
          globalIndex colDOF[1];
          colDOF[0] = presDofNumber[kfe];

          real64 const Ja = area[kfe];
          real64 const nodalArea = Ja / static_cast< real64 >( numNodesPerFace );
          real64 nodalForceMag = ( pressure[kfe] + deltaPressure[kfe] ) * nodalArea;
          R1Tensor globalNodalForce( Nbar );
          globalNodalForce *= nodalForceMag;

          for( localIndex kf=0; kf<2; ++kf )
          {
            localIndex const faceIndex = elemsToFaces[kfe][kf];

            for( localIndex a=0; a<numNodesPerFace; ++a )
            {
              for( localIndex i=0; i<3; ++i )
              {
                rowDOF[3*a+i] = dispDofNumber[faceToNodeMap( faceIndex, a )] + integer_conversion< globalIndex >( i );
                // Opposite sign w.r.t. theory because of minus sign in stiffness matrix definition (K < 0)
                nodeRHS[3*a+i] = +globalNodalForce[i] * pow( -1, kf );
                fext[faceToNodeMap( faceIndex, a )][i] += +globalNodalForce[i] * pow( -1, kf );

                // Opposite sign w.r.t. theory because of minus sign in stiffness matrix definition (K < 0)
                dRdP( 3*a+i ) = -nodalArea * Nbar( i ) * pow( -1, kf );
              }
            }

            rhs->add( rowDOF,
                      nodeRHS,
                      3 * numNodesPerFace );

            matrix->add( rowDOF,
                         colDOF,
                         dRdP.data(),
                         3 * numNodesPerFace,
                         1 );
          }
        }
      } );
    }
  } );

  matrix->close();
  rhs->close();
}

void LagrangianContactFlowSolver::AssembleFluidMassResidualDerivativeWrtDisplacement( DomainPartition const * const domain,
                                                                                      DofManager const & dofManager,
                                                                                      ParallelMatrix * const matrix,
                                                                                      ParallelVector * const rhs )
{
  GEOSX_MARK_FUNCTION;
  MeshLevel const * const mesh = domain->getMeshBody( 0 )->getMeshLevel( 0 );

  FaceManager const * const faceManager = mesh->getFaceManager();
  NodeManager const * const nodeManager = mesh->getNodeManager();
  ElementRegionManager const * const elemManager = mesh->getElemManager();
  ConstitutiveManager const * const constitutiveManager = domain->getConstitutiveManager();

  arrayView1d< R1Tensor const > const & faceNormal = faceManager->faceNormal();
  ArrayOfArraysView< localIndex const > const & faceToNodeMap = faceManager->nodeList();

  CRSMatrixView< real64 const, localIndex const > const &
  dFluxResidual_dAperture = m_flowSolver->getDerivativeFluxResidual_dAperture();

  string const dispDofKey = dofManager.getKey( keys::TotalDisplacement );
  string const presDofKey = dofManager.getKey( m_pressureKey );
  string const constitutiveName = constitutiveManager->GetGroup( m_flowSolver->fluidIndex())->getName();

  arrayView1d< globalIndex const > const &
  dispDofNumber = nodeManager->getReference< globalIndex_array >( dispDofKey );

  matrix->open();
  rhs->open();

  elemManager->forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion const & subRegion )
  {
    if( subRegion.hasWrapper( m_pressureKey ) )
    {
      arrayView1d< globalIndex const > const &
      presDofNumber = subRegion.getReference< globalIndex_array >( presDofKey );
      arrayView1d< integer const > const & ghostRank = subRegion.ghostRank();
      arrayView1d< real64 const > const & area = subRegion.getElementArea();
      arrayView2d< localIndex const > const & elemsToFaces = subRegion.faceList();

      dataRepository::Group const * const constitutiveGroup = subRegion.GetConstitutiveModels();
      dataRepository::Group const * const constitutiveRelation = constitutiveGroup->GetGroup( constitutiveName );
      arrayView2d< real64 const > const &
      density = constitutiveRelation->getReference< array2d< real64 > >( SingleFluidBase::viewKeyStruct::densityString );

      forAll< serialPolicy >( subRegion.size(), [&]( localIndex const kfe )
      {
        if( ghostRank[kfe] < 0 )
        {
          localIndex const numNodesPerFace = faceToNodeMap.sizeOfArray( elemsToFaces[kfe][0] );
          globalIndex nodeDOF[24];
          globalIndex elemDOF[1];
          elemDOF[0] = presDofNumber[kfe];

          real64 const Ja = area[kfe];
          R1Tensor Nbar = faceNormal[elemsToFaces[kfe][0]];
          Nbar -= faceNormal[elemsToFaces[kfe][1]];
          Nbar.Normalize();

          real64 const dAccumulationResidualdAperture = density[kfe][0] * Ja;

          stackArray1d< real64, 2*3*4 > dRdU( 2*3*numNodesPerFace );
          dRdU = 0.0;

          // Accumulation derivative
          if( m_contactSolver->IsElementInOpenState( subRegion, kfe ) )
          {
            for( localIndex kf=0; kf<2; ++kf )
            {
              for( localIndex a=0; a<numNodesPerFace; ++a )
              {
                for( int i=0; i<3; ++i )
                {
                  nodeDOF[ kf*3*numNodesPerFace + 3*a+i]  = dispDofNumber[faceToNodeMap( elemsToFaces[kfe][kf], a )] + integer_conversion< globalIndex >( i );
                  real64 const dAper_dU = -pow( -1, kf ) * Nbar[i] / numNodesPerFace;
                  dRdU( kf*3*numNodesPerFace + 3*a+i ) = dAccumulationResidualdAperture * dAper_dU;
                }
              }
            }
            matrix->add( elemDOF,
                         nodeDOF,
                         dRdU.data(),
                         1,
                         2*3*numNodesPerFace );
          }

          dRdU = 0.0;
          // flux derivative
          bool skipAssembly = true;
          localIndex const numColumns = dFluxResidual_dAperture.numNonZeros( kfe );
          arraySlice1d< localIndex const > const & columns = dFluxResidual_dAperture.getColumns( kfe );
          arraySlice1d< real64 const > const & values = dFluxResidual_dAperture.getEntries( kfe );

          skipAssembly &= !( m_contactSolver->IsElementInOpenState( subRegion, kfe ) );

          for( localIndex kfe1=0; kfe1<numColumns; ++kfe1 )
          {
            real64 const dR_dAper = values[kfe1];
            localIndex const kfe2 = columns[kfe1];

            skipAssembly &= !( m_contactSolver->IsElementInOpenState( subRegion, kfe2 ) );

            for( localIndex kf=0; kf<2; ++kf )
            {
              for( localIndex a=0; a<numNodesPerFace; ++a )
              {
                for( int i=0; i<3; ++i )
                {
                  nodeDOF[ kf*3*numNodesPerFace + 3*a+i ] = dispDofNumber[faceToNodeMap( elemsToFaces[kfe2][kf], a )]
                                                            + integer_conversion< globalIndex >( i );
                  real64 const dAper_dU = -pow( -1, kf ) * Nbar[i] / numNodesPerFace;
                  dRdU( kf*3*numNodesPerFace + 3*a+i ) = dR_dAper * dAper_dU;
                }
              }
            }

            if( !skipAssembly )
            {
              matrix->add( elemDOF,
                           nodeDOF,
                           dRdU.data(),
                           1,
                           2*3*numNodesPerFace );
            }
          }
        }
      } );
    }
  } );

  matrix->close();
  rhs->close();
}

void LagrangianContactFlowSolver::AssembleStabiliziation( DomainPartition const * const domain,
                                                          DofManager const & dofManager,
                                                          ParallelMatrix * const matrix,
                                                          ParallelVector * const rhs )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel const * const mesh = domain->getMeshBody( 0 )->getMeshLevel( 0 );

  FaceManager const * const faceManager = mesh->getFaceManager();
  NodeManager const * const nodeManager = mesh->getNodeManager();
  ElementRegionManager const * const elemManager = mesh->getElemManager();
  ConstitutiveManager const * const constitutiveManager = domain->getConstitutiveManager();

  string const presDofKey = dofManager.getKey( m_pressureKey );
  string const constitutiveName = constitutiveManager->GetGroup( m_flowSolver->fluidIndex())->getName();

  // Get the finite volume method used to compute the stabilization
  NumericalMethodsManager const * const numericalMethodManager = domain->getParent()->GetGroup< NumericalMethodsManager >( keys::numericalMethodsManager );
  FiniteVolumeManager const * const fvManager = numericalMethodManager->GetGroup< FiniteVolumeManager >( keys::finiteVolumeManager );
  FluxApproximationBase const * const stabilizationMethod = fvManager->getFluxApproximation( m_stabilizationName );

  // Get the "face to element" map (valid for the entire mesh)
  FaceManager::ElemMapType const & faceToElem = faceManager->toElementRelation();

  // Form the SurfaceGenerator, get the fracture name and use it to retrieve the faceMap (from fracture element to face)
  SurfaceGenerator const * const
  surfaceGenerator = this->getParent()->GetGroup< SolverBase >( "SurfaceGen" )->group_cast< SurfaceGenerator const * >();
  FaceElementRegion const * const fractureRegion = elemManager->GetRegion< FaceElementRegion >( surfaceGenerator->getFractureRegionName() );
  FaceElementSubRegion const * const fractureSubRegion = fractureRegion->GetSubRegion< FaceElementSubRegion >( "default" );
  GEOSX_ERROR_IF( !fractureSubRegion->hasWrapper( m_pressureKey ), "The fracture subregion must contain pressure field." );
  FaceElementSubRegion::FaceMapType const & faceMap = fractureSubRegion->faceList();
  GEOSX_ERROR_IF( faceMap.size( 1 ) != 2, "A fracture face has to be shared by two cells." );

  // Get the pressures
  arrayView1d< real64 const > const &
  pressure = fractureSubRegion->getReference< array1d< real64 > >( m_pressureKey );
  arrayView1d< real64 const > const &
  deltaPressure = fractureSubRegion->getReference< array1d< real64 > >( m_deltaPressureKey );

  // Get the density
  dataRepository::Group const * const constitutiveGroup = fractureSubRegion->GetConstitutiveModels();
  dataRepository::Group const * const constitutiveRelation = constitutiveGroup->GetGroup( constitutiveName );
  arrayView2d< real64 const > const &
  density = constitutiveRelation->getReference< array2d< real64 > >( SingleFluidBase::viewKeyStruct::densityString );

  // Get the volume for all elements
  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > const elemVolume =
    elemManager->ConstructViewAccessor< real64_array, arrayView1d< real64 const > >( ElementSubRegionBase::viewKeyStruct::elementVolumeString );

  // Get the coordinates for all nodes
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition = nodeManager->referencePosition();

  // Get area and rotation matrix for all faces
  arrayView1d< real64 const > const & faceArea = faceManager->faceArea();
  arrayView1d< R1Tensor const > const & faceNormal = faceManager->faceNormal();

  ElementRegionManager::ConstitutiveRelationAccessor< ConstitutiveBase const > const
  constitutiveRelations = elemManager->ConstructFullConstitutiveAccessor< ConstitutiveBase const >( constitutiveManager );

  arrayView1d< globalIndex const > const &
  presDofNumber = fractureSubRegion->getReference< globalIndex_array >( presDofKey );

  matrix->open();
  rhs->open();

  stabilizationMethod->forStencils< FaceElementStencil >( [&]( FaceElementStencil const & stencil )
  {
    // Get ghost rank
    ArrayOfArraysView< integer const > const & isGhostConnectors = stencil.getIsGhostConnectors();

    for( localIndex iconn=0; iconn<stencil.size(); ++iconn )
    {
      localIndex const numFluxElems = stencil.stencilSize( iconn );

      // A fracture connector has to be an edge shared by two faces
      if( numFluxElems == 2 && isGhostConnectors[iconn][0] < 0 )
      {
        typename FaceElementStencil::IndexContainerViewConstType const & sei = stencil.getElementIndices();

        // First index: face element. Second index: node
        real64_array2d nodalArea( 2, 2 );
        real64_array rotatedInvStiffApprox( 2 );
        for( localIndex kf = 0; kf < 2; ++kf )
        {
          // Get fracture, face and region/subregion/element indices (for elements on both sides)
          localIndex fractureIndex = sei[iconn][kf];

          localIndex faceIndexRef = faceMap[fractureIndex][0];
          real64 const area = faceArea[faceIndexRef];
          R1Tensor const & Nbar = faceNormal[faceIndexRef];
          // TODO: use higher order integration scheme
          nodalArea[kf][0] = area / 4.0;
          nodalArea[kf][1] = area / 4.0;

          real64_array2d invStiffApprox( 2, 3 );
          for( localIndex i = 0; i < 2; ++i )
          {
            localIndex faceIndex = faceMap[fractureIndex][i];
            localIndex er = faceToElem.m_toElementRegion[faceIndex][0];
            localIndex esr = faceToElem.m_toElementSubRegion[faceIndex][0];
            localIndex ei = faceToElem.m_toElementIndex[faceIndex][0];

            real64 const volume = elemVolume[er][esr][ei];

            // Get the "element to node" map for the specific region/subregion
            CellElementSubRegion const * const
            cellElementSubRegion = elemManager->GetRegion( er )->GetSubRegion< CellElementSubRegion >( esr );
            arrayView2d< localIndex const, cells::NODE_MAP_USD > const & cellElemsToNodes = cellElementSubRegion->nodeList();
            localIndex numNodesPerElem = cellElementSubRegion->numNodesPerElement();

            // Compute the box size
            real64_array maxSize( 3 ), minSize( 3 );
            for( localIndex j = 0; j < 3; ++j )
            {
              maxSize[j] = nodePosition[cellElemsToNodes[ei][0]][j];
              minSize[j] = nodePosition[cellElemsToNodes[ei][0]][j];
            }
            for( localIndex a=1; a<numNodesPerElem; ++a )
            {
              for( localIndex j = 0; j < 3; ++j )
              {
                maxSize[j] = std::max( maxSize[j], nodePosition[cellElemsToNodes[ei][a]][j] );
                minSize[j] = std::min( minSize[j], nodePosition[cellElemsToNodes[ei][a]][j] );
              }
            }
            real64_array boxSize( 3 );
            for( localIndex j = 0; j < 3; ++j )
            {
              boxSize[j] = maxSize[j] - minSize[j];
            }

            // Get linear elastic isotropic constitutive parameters for the element
            LinearElasticIsotropic const * const constitutiveRelation0 =
              dynamic_cast< LinearElasticIsotropic const * >( constitutiveRelations[er][esr][m_contactSolver->getSolidSolver()->getSolidMaterialFullIndex()] );
            real64 const K = constitutiveRelation0->bulkModulus()[ei];
            real64 const G = constitutiveRelation0->shearModulus()[ei];
            real64 const E = 9.0 * K * G / ( 3.0 * K + G );
            real64 const nu = ( 3.0 * K - 2.0 * G ) / ( 2.0 * ( 3.0 * K + G ) );

            for( localIndex j = 0; j < 3; ++j )
            {
              invStiffApprox[i][j] = 1.0 / ( E / ( ( 1.0 + nu )*( 1.0 - 2.0*nu ) ) * 4.0 / 9.0 * ( 2.0 - 3.0 * nu ) * volume / ( boxSize[j]*boxSize[j] ) );
            }
          }

          R2Tensor invStiffApproxTotal;
          invStiffApproxTotal = 0.0;
          for( localIndex i = 0; i < 2; ++i )
          {
            for( localIndex j = 0; j < 3; ++j )
            {
              invStiffApproxTotal( j, j ) += invStiffApprox[i][j];
            }
          }
          // Compute n^T * (invK) * n
          R1Tensor tmpTensor;
          tmpTensor.AijBj( invStiffApproxTotal, Nbar );
          rotatedInvStiffApprox[kf] = Dot( Nbar, tmpTensor);
        }

        // Compose local nodal-based local stiffness matrices
        stackArray1d< real64, 1 > totalInvStiffApprox( 1 );
        for( localIndex kf = 0; kf < 2; ++kf )
        {
          rotatedInvStiffApprox[kf] *= nodalArea[0][kf] * nodalArea[1][kf];
        }
        // Local assembly
        totalInvStiffApprox( 0 ) = ( rotatedInvStiffApprox[0] + rotatedInvStiffApprox[1] );

        // Get DOF numbering
        localIndex fractureIndex[2];
        localIndex nDof[2];
        globalIndex elemDOF[2][1];
        for( localIndex kf = 0; kf < 2; ++kf )
        {
          fractureIndex[kf] = sei[iconn][kf];
          elemDOF[kf][0] = -1;
          nDof[kf] = 0;
          if( !m_contactSolver->IsElementInOpenState( *fractureSubRegion, fractureIndex[kf] ) )
          {
            elemDOF[kf][0] = presDofNumber[fractureIndex[kf]];
            nDof[kf] = 1;
          }
        }

        // Add mean density contribution
        totalInvStiffApprox( 0 ) *= 0.5 * ( density[fractureIndex[0]][0] + density[fractureIndex[1]][0] );

        // Compute rhs
        real64 rhs0 = 0.0;
        if( nDof[0] > 0)
        {
          rhs0 -= totalInvStiffApprox( 0 ) * ( pressure[fractureIndex[0]] + deltaPressure[fractureIndex[0]] );
        }
        if( nDof[1] > 0)
        {
          rhs0 += totalInvStiffApprox( 0 ) * ( pressure[fractureIndex[1]] + deltaPressure[fractureIndex[1]] );
        }
        real64 rhs1 = -rhs0;

        // Global matrix and rhs assembly
        if( std::max( nDof[0], nDof[1] ) > 0 )
        {
          matrix->add( elemDOF[0],
                       elemDOF[0],
                       totalInvStiffApprox.data(),
                       nDof[0],
                       nDof[0] );

          matrix->add( elemDOF[1],
                       elemDOF[1],
                       totalInvStiffApprox.data(),
                       nDof[1],
                       nDof[1] );

          // Change sign
          totalInvStiffApprox( 0 ) *= -1.0;

          matrix->add( elemDOF[0],
                       elemDOF[1],
                       totalInvStiffApprox.data(),
                       nDof[0],
                       nDof[1] );

          matrix->add( elemDOF[1],
                       elemDOF[0],
                       totalInvStiffApprox.data(),
                       nDof[1],
                       nDof[0] );

          rhs->add( elemDOF[0],
                    &rhs0,
                    nDof[0] );
          rhs->add( elemDOF[1],
                    &rhs1,
                    nDof[1] );
        }
      }
    }
  } );

  matrix->close();
  rhs->close();
}

void LagrangianContactFlowSolver::ApplySystemSolution( DofManager const & dofManager,
                                                       ParallelVector const & solution,
                                                       real64 const scalingFactor,
                                                       DomainPartition * const domain )
{
  GEOSX_MARK_FUNCTION;

  m_contactSolver->ApplySystemSolution( dofManager,
                                        solution,
                                        scalingFactor,
                                        domain );

  m_flowSolver->ApplySystemSolution( dofManager,
                                     solution,
                                     scalingFactor,
                                     domain );

  UpdateOpeningForFlow( domain );
}

void LagrangianContactFlowSolver::SolveSystem( DofManager const & dofManager,
                                               ParallelMatrix & matrix,
                                               ParallelVector & rhs,
                                               ParallelVector & solution )
{
  GEOSX_MARK_FUNCTION;

  if( getLogLevel() > 3 )
  {
    matrix.write( "matrix.mtx", LAIOutputFormat::MATRIX_MARKET );
    rhs.write( "rhs.mtx", LAIOutputFormat::MATRIX_MARKET );
  }

  SolverBase::SolveSystem( dofManager, matrix, rhs, solution );

  if( getLogLevel() > 3 )
  {
    solution.write( "sol.mtx", LAIOutputFormat::MATRIX_MARKET );
  }

  int rank = MpiWrapper::Comm_rank( MPI_COMM_GEOSX );
  if( rank == 0 )
  {
    string str;
    std::getline( std::cin, str );
    if( str.length() > 0 )
    {
      GEOSX_ERROR( "STOP" );
    }
  }
  MpiWrapper::Barrier( MPI_COMM_GEOSX );
}

void LagrangianContactFlowSolver::SetNextDt( real64 const & currentDt,
                                             real64 & nextDt )
{
  nextDt = currentDt;
}

REGISTER_CATALOG_ENTRY( SolverBase, LagrangianContactFlowSolver, std::string const &, Group * const )
} /* namespace geosx */
