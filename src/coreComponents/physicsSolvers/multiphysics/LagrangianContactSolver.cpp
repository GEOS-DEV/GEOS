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
 * @file LagrangianContactSolver.cpp
 *
 */

#include "LagrangianContactSolver.hpp"

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
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"
#include "physicsSolvers/surfaceGeneration/SurfaceGenerator.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"
#include "linearAlgebra/utilities/LAIHelperFunctions.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;

LagrangianContactSolver::LagrangianContactSolver( const std::string & name,
                                                  Group * const parent ):
  SolverBase( name, parent ),
  m_solidSolverName(),
  m_solidSolver( nullptr ),
  m_stabilizationName()
{
  registerWrapper( viewKeyStruct::solidSolverNameString, &m_solidSolverName, false )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Name of the solid mechanics solver to use in the lagrangian contact solver" );

  registerWrapper( viewKeyStruct::stabilizationNameString, &m_stabilizationName, false )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Name of the stabilization to use in the lagrangian contact solver" );

  registerWrapper( viewKeyStruct::activeSetMaxIterString, &m_activeSetMaxIter, false )->
    setApplyDefaultValue( 10 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Maximum number of iteration for the active set strategy in the lagrangian contact solver" );
}

void LagrangianContactSolver::RegisterDataOnMesh( dataRepository::Group * const MeshBodies )
{
  GEOSX_LOG_RANK( "In RegisterDataOnMesh!" );

  for( auto & mesh : MeshBodies->GetSubGroups() )
  {
    MeshLevel * meshLevel = Group::group_cast< MeshBody * >( mesh.second )->getMeshLevel( 0 );

    ElementRegionManager * const elemManager = meshLevel->getElemManager();
    elemManager->forElementRegions< FaceElementRegion >( [&] ( FaceElementRegion * const region )
    {
      region->forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion * const subRegion )
      {
        subRegion->registerWrapper< array2d< real64 > >( viewKeyStruct::tractionString )->
          setApplyDefaultValue( 0.0 )->
          setPlotLevel( PlotLevel::LEVEL_0 )->
          setRegisteringObjects( this->getName())->
          setDescription( "An array that holds the tractions on the fracture." );
        subRegion->getReference< array2d< real64 > >( viewKeyStruct::tractionString ).resizeDimension< 1 >( 3 );

        subRegion->registerWrapper< array2d< real64 > >( viewKeyStruct::deltaTractionString )->
          setApplyDefaultValue( 0.0 )->
          setPlotLevel( PlotLevel::NOPLOT )->
          setRegisteringObjects( this->getName())->
          setDescription( "An array that holds the traction increments on the fracture." );
        subRegion->getReference< array2d< real64 > >( viewKeyStruct::deltaTractionString ).resizeDimension< 1 >( 3 );

        subRegion->registerWrapper< array1d< FractureState > >( viewKeyStruct::fractureStateString )->
          setPlotLevel( PlotLevel::LEVEL_1 )->
          setRegisteringObjects( this->getName())->
          setDescription( "An array that holds the fracture state." );
        InitializeFractureState( meshLevel );

        subRegion->registerWrapper< array2d< real64 > >( viewKeyStruct::localJumpString )->
          setApplyDefaultValue( 0.0 )->
          setPlotLevel( PlotLevel::LEVEL_0 )->
          setRegisteringObjects( this->getName())->
          setDescription( "An array that holds the local jump on the fracture at the current time step." );
        subRegion->getReference< array2d< real64 > >( viewKeyStruct::localJumpString ).resizeDimension< 1 >( 3 );

        subRegion->registerWrapper< array2d< real64 > >( viewKeyStruct::previousLocalJumpString )->
          setApplyDefaultValue( 0.0 )->
          setPlotLevel( PlotLevel::NOPLOT )->
          setRegisteringObjects( this->getName())->
          setDescription( "An array that holds the local jump on the fracture at the previous time step." );
        subRegion->getReference< array2d< real64 > >( viewKeyStruct::previousLocalJumpString ).resizeDimension< 1 >( 3 );

        subRegion->registerWrapper< array2d< real64 > >( viewKeyStruct::localJumpCorrectionString )->
          setApplyDefaultValue( 0.0 )->
          setPlotLevel( PlotLevel::NOPLOT )->
          setRegisteringObjects( this->getName())->
          setDescription( "An array that holds the jump correction due to the stabilization contribution at the current time step." );
        subRegion->getReference< array2d< real64 > >( viewKeyStruct::localJumpCorrectionString ).resizeDimension< 1 >( 3 );

        subRegion->registerWrapper< array2d< real64 > >( viewKeyStruct::previousLocalJumpCorrectionString )->
          setApplyDefaultValue( 0.0 )->
          setPlotLevel( PlotLevel::NOPLOT )->
          setRegisteringObjects( this->getName())->
          setDescription( "An array that holds the jump correction due to the stabilization contribution at the previous time step." );
        subRegion->getReference< array2d< real64 > >( viewKeyStruct::previousLocalJumpCorrectionString ).resizeDimension< 1 >( 3 );

        // Needed just because SurfaceGenerator initialize the field "pressure" (NEEDED!!!)
        subRegion->registerWrapper< array1d< real64 > >( "pressure" )->
          setPlotLevel( PlotLevel::NOPLOT )->
          setRegisteringObjects( this->getName());
      } );
    } );
  }
}

void LagrangianContactSolver::ImplicitStepSetup( real64 const & time_n,
                                                 real64 const & dt,
                                                 DomainPartition * const domain,
                                                 DofManager & dofManager,
                                                 ParallelMatrix & matrix,
                                                 ParallelVector & rhs,
                                                 ParallelVector & solution )
{
  this->UpdateDeformationForCoupling( domain );

  m_solidSolver->ImplicitStepSetup( time_n, dt, domain,
                                    dofManager,
                                    matrix,
                                    rhs,
                                    solution );
}

void LagrangianContactSolver::ImplicitStepComplete( real64 const & time_n,
                                                    real64 const & dt,
                                                    DomainPartition * const domain )
{
  GEOSX_LOG_RANK( "In ImplicitStepComplete!" );

  m_solidSolver->ImplicitStepComplete( time_n, dt, domain );

  MeshLevel * const meshLevel = domain->getMeshBody( 0 )->getMeshLevel( 0 );
  ElementRegionManager * const elemManager = meshLevel->getElemManager();

  elemManager->forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion * const subRegion )->void
  {
    if( subRegion->hasWrapper( m_tractionKey ) )
    {
      arrayView2d< real64 const > const &
      localJump = subRegion->getReference< array2d< real64 > >( viewKeyStruct::localJumpString );
      arrayView2d< real64 > const &
      previousLocalJump = subRegion->getReference< array2d< real64 > >( viewKeyStruct::previousLocalJumpString );

      forall_in_range< serialPolicy >( 0,
                                       subRegion->size(),
                                       GEOSX_LAMBDA ( localIndex const kfe )
        {
          for( localIndex i = 0 ; i < 3 ; ++i )
          {
            previousLocalJump[kfe][i] = localJump[kfe][i];
          }
        } );
    }
  } );
  GEOSX_LOG_RANK( "here ..." );
}

void LagrangianContactSolver::PostProcessInput()
{
  m_solidSolver = this->getParent()->GetGroup< SolidMechanicsLagrangianFEM >( m_solidSolverName );
  GEOSX_ERROR_IF( m_solidSolver == nullptr, this->getName() << ": invalid solid solver name: " << m_solidSolverName );

  SolverBase::PostProcessInput();
}

void LagrangianContactSolver::InitializePostInitialConditions_PreSubGroups( Group * const GEOSX_UNUSED_PARAM( problemManager ) )
{}

LagrangianContactSolver::~LagrangianContactSolver()
{
  // TODO Auto-generated destructor stub
}

void LagrangianContactSolver::ResetStateToBeginningOfStep( DomainPartition * const domain )
{
  GEOSX_LOG_RANK( "In ResetStateToBeginningOfStep!" );

  m_solidSolver->ResetStateToBeginningOfStep( domain );

  MeshLevel * const meshLevel = domain->getMeshBody( 0 )->getMeshLevel( 0 );
  ElementRegionManager * const elemManager = meshLevel->getElemManager();

  elemManager->forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion * const subRegion )->void
  {
    if( subRegion->hasWrapper( m_tractionKey ) )
    {
      arrayView2d< real64 > const & traction = subRegion->getReference< array2d< real64 > >( viewKeyStruct::tractionString );
      arrayView2d< real64 > const & deltaTraction = subRegion->getReference< array2d< real64 > >( viewKeyStruct::deltaTractionString );
      arrayView2d< real64 > const & localJump = subRegion->getReference< array2d< real64 > >( viewKeyStruct::localJumpString );
      arrayView2d< real64 const > const & previousLocalJump = subRegion->getReference< array2d< real64 > >( viewKeyStruct::previousLocalJumpString );
      arrayView2d< real64 > const & jumpCorrection = subRegion->getReference< array2d< real64 > >( viewKeyStruct::localJumpCorrectionString );
      arrayView2d< real64 > const & previousLocalJumpCorrection = subRegion->getReference< array2d< real64 > >(
        viewKeyStruct::previousLocalJumpCorrectionString );

      forall_in_range< serialPolicy >( 0,
                                       subRegion->size(),
                                       GEOSX_LAMBDA ( localIndex const kfe )
        {
          for( localIndex i = 0 ; i < 3 ; ++i )
          {
            traction[kfe][i] -= deltaTraction[kfe][i];
            deltaTraction[kfe][i] = 0.0;

            localJump[kfe][i] = previousLocalJump[kfe][i];
            jumpCorrection[kfe][i] = previousLocalJumpCorrection[kfe][i];
            previousLocalJumpCorrection[kfe][i] = 0.0;
          }
        } );
    }
  } );

}

real64 LagrangianContactSolver::SolverStep( real64 const & time_n,
                                            real64 const & dt,
                                            int const cycleNumber,
                                            DomainPartition * const domain )
{
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

  m_solidSolver->updateStress( domain );

  // final step for completion of timestep. typically secondary variable updates and cleanup.
  ImplicitStepComplete( time_n, dtReturn, domain );

  return dtReturn;
}

void LagrangianContactSolver::UpdateDeformationForCoupling( DomainPartition * const domain )
{
  MeshLevel * const meshLevel = domain->getMeshBody( 0 )->getMeshLevel( 0 );

  NodeManager const * const nodeManager = meshLevel->getNodeManager();
  FaceManager * const faceManager = meshLevel->getFaceManager();
  ElementRegionManager * const elemManager = meshLevel->getElemManager();

  ArrayOfArraysView< localIndex const > const & faceToNodeMap = faceManager->nodeList();

  arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const & u = nodeManager->totalDisplacement();

  elemManager->forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion * const subRegion )->void
  {
    if( subRegion->hasWrapper( m_tractionKey ) )
    {
      arrayView1d< integer const > const & ghostRank = subRegion->GhostRank();
      arrayView1d< R2Tensor const > const & rotationMatrix = subRegion->getElementRotationMatrix();
      arrayView2d< localIndex const > const & elemsToFaces = subRegion->faceList();
      arrayView2d< real64 > const & localJump = subRegion->getReference< array2d< real64 > >( viewKeyStruct::localJumpString );

      forall_in_range< serialPolicy >( 0,
                                       subRegion->size(),
                                       GEOSX_LAMBDA ( localIndex const kfe )
        {
          {
            localIndex const numNodesPerFace = faceToNodeMap.sizeOfArray( elemsToFaces[kfe][0] );

            // Contact constraints
            if( ghostRank[kfe] < 0 )
            {
              R1Tensor globalJump( 0.0, 0.0, 0.0 );
              for( localIndex a=0 ; a<numNodesPerFace ; ++a )
              {
                for( localIndex i=0 ; i<3 ; ++i )
                {
                  globalJump( i ) +=
                    ( -u[faceToNodeMap( elemsToFaces[kfe][0], a )][i]
                      + u[faceToNodeMap( elemsToFaces[kfe][1], a )][i] ) / static_cast< real64 >(numNodesPerFace);
                }
              }
              R1Tensor localJumpTensor;
              localJumpTensor.AijBi( rotationMatrix[kfe], globalJump );
              localJump[kfe][0] = localJumpTensor( 0 );
              localJump[kfe][1] = localJumpTensor( 1 );
              localJump[kfe][2] = localJumpTensor( 2 );

//              GEOSX_LOG_RANK( "Element: " << kfe << " localJump: " << localJumpTensor );
            }
          }
        } );
    }
  } );

  return;
}

real64 LagrangianContactSolver::SplitOperatorStep( real64 const & GEOSX_UNUSED_PARAM( time_n ),
                                                   real64 const & dt,
                                                   integer const GEOSX_UNUSED_PARAM( cycleNumber ),
                                                   DomainPartition * const GEOSX_UNUSED_PARAM( domain ) )
{
  real64 dtReturn = dt;
  return dtReturn;
}

real64 LagrangianContactSolver::ExplicitStep( real64 const & GEOSX_UNUSED_PARAM( time_n ),
                                              real64 const & dt,
                                              const int GEOSX_UNUSED_PARAM( cycleNumber ),
                                              DomainPartition * const GEOSX_UNUSED_PARAM( domain ) )
{
  GEOSX_MARK_FUNCTION;
  GEOSX_ERROR( "ExplicitStep non available for LagrangianContactSolver!" );
  return dt;
}

real64 LagrangianContactSolver::NonlinearImplicitStep( real64 const & time_n,
                                                       real64 const & dt,
                                                       integer const cycleNumber,
                                                       DomainPartition * const domain,
                                                       DofManager const & dofManager,
                                                       ParallelMatrix & matrix,
                                                       ParallelVector & rhs,
                                                       ParallelVector & solution )
{
  real64 stepDt = dt;

  bool converged = false;
  integer & activeSetIter = m_activeSetIter;
  for( activeSetIter = 0 ; activeSetIter < m_activeSetMaxIter ; ++activeSetIter )
  {
    stepDt = SolverBase::NonlinearImplicitStep( time_n,
                                                dt,
                                                cycleNumber,
                                                domain,
                                                dofManager,
                                                matrix,
                                                rhs,
                                                solution );

    bool const isPreviousFractureStateValid = UpdateFractureState( domain );
    GEOSX_LOG_RANK_0( "active set flag ---> " << isPreviousFractureStateValid );
    if( isPreviousFractureStateValid )
    {
      converged = true;
      break;
    }
    else
    {
      GEOSX_LOG_RANK_0( "**** NonlinearImplicitStep: done! ****" );
    }
  }
  if( !converged )
  {
    GEOSX_ERROR( "Active set did not reached a solution. Terminating..." );
  }
  else
  {
    GEOSX_LOG_RANK_0( "Number of active set iterations: " << activeSetIter );
  }

  return stepDt;
}

void LagrangianContactSolver::SetupDofs( DomainPartition const * const domain,
                                         DofManager & dofManager ) const
{
  GEOSX_MARK_FUNCTION;
  m_solidSolver->SetupDofs( domain, dofManager );

  // restrict coupling to fracture regions only
  ElementRegionManager const * const elemManager = domain->getMeshBody( 0 )->getMeshLevel( 0 )->getElemManager();
  string_array fractureRegions;
  elemManager->forElementRegions< FaceElementRegion >( [&]( FaceElementRegion const * const elementRegion )
  {
    fractureRegions.push_back( elementRegion->getName() );
  } );

  dofManager.addField( viewKeyStruct::tractionString,
                       DofManager::Location::Elem,
                       3,
                       fractureRegions );
  dofManager.addCoupling( viewKeyStruct::tractionString,
                          viewKeyStruct::tractionString,
                          DofManager::Connectivity::Face,
                          fractureRegions );
  dofManager.addCoupling( keys::TotalDisplacement,
                          viewKeyStruct::tractionString,
                          DofManager::Connectivity::Elem,
                          fractureRegions );
}

void LagrangianContactSolver::SetupSystem( DomainPartition * const domain,
                                           DofManager & dofManager,
                                           ParallelMatrix & matrix,
                                           ParallelVector & rhs,
                                           ParallelVector & solution )
{
  GEOSX_MARK_FUNCTION;

  // setup DofManager
  dofManager.setMesh( domain, 0, 0 );

  // add traction and coupling
  SetupDofs( domain, dofManager );
  dofManager.reorderByRank();

  localIndex const numDisplacementDofs = dofManager.numLocalDofs( keys::TotalDisplacement );
  localIndex const numTractionDofs = dofManager.numLocalDofs( viewKeyStruct::tractionString );
  GEOSX_LOG_RANK( numDisplacementDofs << " " << numTractionDofs );

  matrix.createWithLocalSize( numDisplacementDofs + numTractionDofs,
                              numDisplacementDofs + numTractionDofs,
                              MPI_COMM_GEOSX );
  rhs.createWithLocalSize( numDisplacementDofs + numTractionDofs,
                           MPI_COMM_GEOSX );
  solution.createWithLocalSize( numDisplacementDofs + numTractionDofs,
                                MPI_COMM_GEOSX );

  dofManager.setSparsityPattern( matrix );
}

void LagrangianContactSolver::AssembleSystem( real64 const time,
                                              real64 const dt,
                                              DomainPartition * const domain,
                                              DofManager const & dofManager,
                                              ParallelMatrix & matrix,
                                              ParallelVector & rhs )
{
  GEOSX_MARK_FUNCTION;

  GEOSX_LOG_RANK( time << " " << dt );

  matrix.open();
  matrix.zero();
  rhs.open();
  rhs.zero();

  m_solidSolver->AssembleSystem( time,
                                 dt,
                                 domain,
                                 dofManager,
                                 matrix,
                                 rhs );

  AssembleForceResidualDerivativeWrtTraction( domain, dofManager, &matrix, &rhs );
  AssembleTractionResidualDerivativeWrtDisplacementAndTraction( domain, dofManager, &matrix, &rhs );
  AssembleStabliziation( domain, dofManager, &matrix, &rhs );

  matrix.close();
  rhs.close();
}

void LagrangianContactSolver::ApplyBoundaryConditions( real64 const time,
                                                       real64 const dt,
                                                       DomainPartition * const domain,
                                                       DofManager const & dofManager,
                                                       ParallelMatrix & matrix,
                                                       ParallelVector & rhs )
{
  GEOSX_MARK_FUNCTION;
  m_solidSolver->ApplyBoundaryConditions( time,
                                          dt,
                                          domain,
                                          dofManager,
                                          matrix,
                                          rhs );

}

real64 LagrangianContactSolver::CalculateResidualNorm( DomainPartition const * const GEOSX_UNUSED_PARAM( domain ),
                                                       DofManager const & dofManager,
                                                       ParallelVector const & rhs )
{
  GEOSX_MARK_FUNCTION;
  GEOSX_LOG_RANK( "In CalculateResidualNorm!" );

  string const tracDofKey = dofManager.getKey( viewKeyStruct::tractionString );
  string const dispDofKey = dofManager.getKey( keys::TotalDisplacement );

  localIndex numDispDofs = dofManager.numLocalDofs( keys::TotalDisplacement );
  localIndex numTracDofs = dofManager.numLocalDofs( viewKeyStruct::tractionString );
  real64 const * localResidual = rhs.extractLocalVector();
  real64 localResidualNorm[2] = {0.0, 0.0};
  for( localIndex i=0 ; i<numDispDofs ; ++i )
  {
    localResidualNorm[0] += localResidual[i] * localResidual[i];
  }
  for( localIndex i=numDispDofs ; i<numDispDofs+numTracDofs ; ++i )
  {
    localResidualNorm[1] += localResidual[i] * localResidual[i];
  }

  real64 globalResidualNorm[2] = {0.0, 0.0};
  int const rank = MpiWrapper::Comm_rank( MPI_COMM_GEOSX );
  int const size = MpiWrapper::Comm_size( MPI_COMM_GEOSX );
  array1d< real64 > globalValues( 2*size );
  globalValues = 0;

  // Everything is done on rank 0
  MpiWrapper::gather( localResidualNorm,
                      2,
                      globalValues.data(),
                      2,
                      0,
                      MPI_COMM_GEOSX );

  if( rank==0 )
  {
    for( int r=0 ; r<size ; ++r )
    {
      // sum across all ranks
      globalResidualNorm[0] += globalValues[2*r];
      globalResidualNorm[1] += globalValues[2*r+1];
    }
    globalResidualNorm[0] = sqrt( globalResidualNorm[0] );
    globalResidualNorm[1] = sqrt( globalResidualNorm[1] );
  }

  MpiWrapper::bcast( globalResidualNorm, 2, 0, MPI_COMM_GEOSX );

  if( m_nonlinearSolverParameters.m_numNewtonIterations == 0 )
  {
    m_initialResidual[0] = globalResidualNorm[0];
    m_initialResidual[1] = globalResidualNorm[1];
    globalResidualNorm[0] = 1.0;
    globalResidualNorm[1] = 1.0;
  }
  else
  {
    globalResidualNorm[0] /= (m_initialResidual[0]+1.0);
    globalResidualNorm[1] /= (m_initialResidual[1]+1.0);
  }

  char output[69] = {0};
  sprintf( output,
           "( Rdisplacement, Rtraction ) = ( %15.6e, %15.6e );",
           globalResidualNorm[0],
           globalResidualNorm[1] );
  GEOSX_LOG_RANK_0( output );

  return ( globalResidualNorm[0] + globalResidualNorm[1] );
}

void LagrangianContactSolver::AssembleForceResidualDerivativeWrtTraction( DomainPartition * const domain,
                                                                          DofManager const & dofManager,
                                                                          ParallelMatrix * const matrix,
                                                                          ParallelVector * const rhs )
{
  GEOSX_MARK_FUNCTION;
  MeshLevel * const mesh = domain->getMeshBody( 0 )->getMeshLevel( 0 );

  FaceManager const * const faceManager = mesh->getFaceManager();
  NodeManager * const nodeManager = mesh->getNodeManager();
  ElementRegionManager * const elemManager = mesh->getElemManager();

  ArrayOfArraysView< localIndex const > const & faceToNodeMap = faceManager->nodeList();

  arrayView1d< R1Tensor > const &
  fext = nodeManager->getReference< array1d< R1Tensor > >( SolidMechanicsLagrangianFEM::viewKeyStruct::forceExternal );
  fext = {0, 0, 0};

  string const tracDofKey = dofManager.getKey( viewKeyStruct::tractionString );
  string const dispDofKey = dofManager.getKey( keys::TotalDisplacement );

  arrayView1d< globalIndex > const &
  dispDofNumber = nodeManager->getReference< globalIndex_array >( dispDofKey );

  elemManager->forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion const * const subRegion )->void
  {
    if( subRegion->hasWrapper( m_tractionKey ) )
    {
      arrayView1d< globalIndex const > const &
      tracDofNumber = subRegion->getReference< array1d< globalIndex > >( tracDofKey );
      arrayView2d< real64 const > const & traction = subRegion->getReference< array2d< real64 > >( viewKeyStruct::tractionString );
      arrayView1d< integer const > const & ghostRank = subRegion->GhostRank();
      arrayView1d< real64 const > const & area = subRegion->getElementArea();
      arrayView1d< R2Tensor const > const & rotationMatrix = subRegion->getElementRotationMatrix();
      arrayView2d< localIndex const > const & elemsToFaces = subRegion->faceList();

      forall_in_range< serialPolicy >( 0,
                                       subRegion->size(),
                                       GEOSX_LAMBDA( localIndex const kfe )
        {
          localIndex const kf0 = elemsToFaces[kfe][0];
          localIndex const numNodesPerFace = faceToNodeMap.sizeOfArray( kf0 );

          globalIndex rowDOF[12];
          real64 nodeRHS[12];
          stackArray2d< real64, 3*4*3 > dRdT( 3*numNodesPerFace, 3 );
          dRdT = 0.0;
          globalIndex colDOF[3];
          for( localIndex i=0 ; i<3 ; ++i )
          {
            colDOF[i] = tracDofNumber[kfe] + integer_conversion< globalIndex >( i );
          }

          real64 const nodalArea = area[kfe] / static_cast< real64 >( numNodesPerFace );
          real64_array nodalForceVec( 3 );
          nodalForceVec[0] = ( traction[kfe][0] ) * nodalArea;
          nodalForceVec[1] = ( traction[kfe][1] ) * nodalArea;
          nodalForceVec[2] = ( traction[kfe][2] ) * nodalArea;
          R1Tensor localNodalForce( nodalForceVec[0], nodalForceVec[1], nodalForceVec[2] );
          R1Tensor globalNodalForce;
          globalNodalForce.AijBj( rotationMatrix[kfe], localNodalForce );

          for( localIndex kf=0 ; kf<2 ; ++kf )
          {
            localIndex const faceIndex = elemsToFaces[kfe][kf];

            for( localIndex a=0 ; a<numNodesPerFace ; ++a )
            {
              for( localIndex i=0 ; i<3 ; ++i )
              {
                rowDOF[3*a+i] = dispDofNumber[faceToNodeMap( faceIndex, a )] + integer_conversion< globalIndex >( i );
                // Opposite sign w.r.t. theory because of minus sign in stiffness matrix definition (K < 0)
                nodeRHS[3*a+i] = +globalNodalForce[i] * pow( -1, kf );
                fext[faceToNodeMap( faceIndex, a )][i] += -globalNodalForce[i] * pow( -1, kf );

                // Opposite sign w.r.t. theory because of minus sign in stiffness matrix definition (K < 0)
                dRdT( 3*a+i, 0 ) = +nodalArea * rotationMatrix[kfe]( i, 0 ) * pow( -1, kf );
                dRdT( 3*a+i, 1 ) = +nodalArea * rotationMatrix[kfe]( i, 1 ) * pow( -1, kf );
                dRdT( 3*a+i, 2 ) = +nodalArea * rotationMatrix[kfe]( i, 2 ) * pow( -1, kf );
              }
            }

            if( ghostRank[kfe] < 0 )
            {
              rhs->add( rowDOF,
                        nodeRHS,
                        3 * numNodesPerFace );

              matrix->add( rowDOF,
                           colDOF,
                           dRdT.data(),
                           3 * numNodesPerFace,
                           3 );
            }
          }
        } );
    }
  } );
}

void LagrangianContactSolver::AssembleTractionResidualDerivativeWrtDisplacementAndTraction( DomainPartition const * const domain,
                                                                                            DofManager const & dofManager,
                                                                                            ParallelMatrix * const matrix,
                                                                                            ParallelVector * const rhs )
{
  GEOSX_MARK_FUNCTION;
  MeshLevel const * const mesh = domain->getMeshBody( 0 )->getMeshLevel( 0 );

  FaceManager const * const faceManager = mesh->getFaceManager();
  NodeManager const * const nodeManager = mesh->getNodeManager();
  ElementRegionManager const * const elemManager = mesh->getElemManager();

  ArrayOfArraysView< localIndex const > const & faceToNodeMap = faceManager->nodeList();

  string const tracDofKey = dofManager.getKey( viewKeyStruct::tractionString );
  string const dispDofKey = dofManager.getKey( keys::TotalDisplacement );

  arrayView1d< globalIndex const > const &
  dispDofNumber = nodeManager->getReference< array1d< globalIndex > >( dispDofKey );

  elemManager->forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion const * const subRegion )->void
  {
    if( subRegion->hasWrapper( m_tractionKey ) )
    {
      arrayView1d< globalIndex const > const &
      tracDofNumber = subRegion->getReference< array1d< globalIndex > >( tracDofKey );
      arrayView1d< integer const > const & ghostRank = subRegion->GhostRank();
      arrayView1d< real64 const > const & area = subRegion->getElementArea();
      arrayView1d< R2Tensor const > const & rotationMatrix = subRegion->getElementRotationMatrix();
      arrayView2d< localIndex const > const & elemsToFaces = subRegion->faceList();
      arrayView2d< real64 const > const &
      traction = subRegion->getReference< array2d< real64 > >( viewKeyStruct::tractionString );
      arrayView1d< FractureState const > const &
      fractureState = subRegion->getReference< array1d< FractureState > >( viewKeyStruct::fractureStateString );
      arrayView2d< real64 const > const &
      localJump = subRegion->getReference< array2d< real64 > >( viewKeyStruct::localJumpString );
      arrayView2d< real64 const > const &
      previousLocalJump = subRegion->getReference< array2d< real64 > >( viewKeyStruct::previousLocalJumpString );
      arrayView2d< real64 const > const &
      previousLocalJumpCorrection = subRegion->getReference< array2d< real64 > >( viewKeyStruct::previousLocalJumpCorrectionString );

      forall_in_range< serialPolicy >( 0,
                                       subRegion->size(),
                                       GEOSX_LAMBDA ( localIndex const kfe )
        {
          {
            localIndex const numNodesPerFace = faceToNodeMap.sizeOfArray( elemsToFaces[kfe][0] );
            globalIndex nodeDOF[24];
            globalIndex elemDOF[3];
            for( localIndex i=0 ; i<3 ; ++i )
            {
              elemDOF[i] = tracDofNumber[kfe] + integer_conversion< globalIndex >( i );
            }

            real64 elemRHS[3];
            real64 const Ja = area[kfe];
            real64 const nodalArea = Ja / static_cast< real64 >( numNodesPerFace );

            stackArray2d< real64, 2*3*4*3 > dRdU( 3, 2*3*numNodesPerFace );
            stackArray2d< real64, 3*3 > dRdT( 3, 3 );
            dRdU = 0.0;
            dRdT = 0.0;

            // Contact constraints
            if( ghostRank[kfe] < 0 )
            {
//              GEOSX_LOG_RANK( "localJump[" << kfe << "]: " << localJump[kfe] << " | previousLocalJump[" << kfe << "]:
// " << previousLocalJump[kfe] );
              switch( fractureState[kfe] )
              {
                case FractureState::STICK:
                  {
                    for( localIndex i=0 ; i<3 ; ++i )
                    {
                      if( i == 0 )
                      {
                        elemRHS[i] = +Ja * localJump[kfe][i];
                      }
                      else
                      {
                        elemRHS[i] = +Ja * ( localJump[kfe][i] - previousLocalJump[kfe][i] );
                      }
                    }

                    for( localIndex kf=0 ; kf<2 ; ++kf )
                    {
                      for( localIndex a=0 ; a<numNodesPerFace ; ++a )
                      {
                        for( localIndex i=0 ; i<3 ; ++i )
                        {
                          nodeDOF[ kf*3*numNodesPerFace + 3*a+
                                   i ] = dispDofNumber[faceToNodeMap( elemsToFaces[kfe][kf], a )] + integer_conversion< globalIndex >( i );

                          dRdU( 0, kf*3*numNodesPerFace + 3*a+i ) = -nodalArea * rotationMatrix[kfe]( 0, i ) * pow( -1, kf );
                          dRdU( 1, kf*3*numNodesPerFace + 3*a+i ) = -nodalArea * rotationMatrix[kfe]( 1, i ) * pow( -1, kf );
                          dRdU( 2, kf*3*numNodesPerFace + 3*a+i ) = -nodalArea * rotationMatrix[kfe]( 2, i ) * pow( -1, kf );
                        }
                      }
                    }
                    break;
                  }
                case FractureState::SLIP:
                case FractureState::NEW_SLIP:
                  {
                    elemRHS[0] = +Ja * localJump[kfe][0];

                    for( localIndex kf=0 ; kf<2 ; ++kf )
                    {
                      for( localIndex a=0 ; a<numNodesPerFace ; ++a )
                      {
                        for( localIndex i=0 ; i<3 ; ++i )
                        {
                          nodeDOF[ kf*3*numNodesPerFace + 3*a+
                                   i ] = dispDofNumber[faceToNodeMap( elemsToFaces[kfe][kf], a )] + integer_conversion< globalIndex >( i );
                          dRdU( 0, kf*3*numNodesPerFace + 3*a+i ) = -nodalArea * rotationMatrix[kfe]( 0, i ) * pow( -1, kf );
                        }
                      }
                    }

                    real64 limitTau = m_cohesion - traction[kfe][0] * std::tan( m_frictionAngle );

                    R1TensorT< 2 > sliding( localJump[kfe][1] - previousLocalJump[kfe][1], localJump[kfe][2] - previousLocalJump[kfe][2] );
                    if( fractureState[kfe] == FractureState::NEW_SLIP )
                    {
                      GEOSX_LOG_RANK( "new slip " << previousLocalJump[kfe][1] << " " << previousLocalJumpCorrection[kfe][1] );
                      sliding( 0 ) -= previousLocalJumpCorrection[kfe][1];
                      sliding( 1 ) -= previousLocalJumpCorrection[kfe][2];
                    }
                    real64 slidingNorm = sqrt( sliding( 0 )*sliding( 0 ) + sliding( 1 )*sliding( 1 ) );

//                  if( !( ( m_nonlinearSolverParameters.m_numNewtonIterations == 0 ) && ( fractureState[kfe] ==
// FractureState::NEW_SLIP ) )
//                      && slidingNorm > m_slidingTolerance )
                    if( slidingNorm > m_slidingTolerance )
                    {
                      for( localIndex i=1 ; i<3 ; ++i )
                      {
                        elemRHS[i] = +Ja * ( traction[kfe][i] - limitTau * sliding( i-1 ) / slidingNorm );
                      }

                      R2TensorT< 2 > dUdgT;
                      dUdgT.dyadic_aa( sliding );
                      dUdgT( 0, 0 ) = (slidingNorm*slidingNorm - dUdgT( 0, 0 )) * limitTau / std::pow( slidingNorm, 3 );
                      dUdgT( 0, 1 ) *= -limitTau / std::pow( slidingNorm, 3 );
                      dUdgT( 1, 0 ) *= -limitTau / std::pow( slidingNorm, 3 );
                      dUdgT( 1, 1 ) = (slidingNorm*slidingNorm - dUdgT( 1, 1 )) * limitTau / std::pow( slidingNorm, 3 );

                      for( localIndex kf=0 ; kf<2 ; ++kf )
                      {
                        for( localIndex a=0 ; a<numNodesPerFace ; ++a )
                        {
                          for( localIndex i=0 ; i<3 ; ++i )
                          {
                            R1TensorT< 2 > localRowB( rotationMatrix[kfe]( 1, i ), rotationMatrix[kfe]( 2, i ) );
                            R1TensorT< 2 > localRowE;
                            localRowE.AijBj( dUdgT, localRowB );

                            dRdU( 1, kf*3*numNodesPerFace + 3*a+i ) = nodalArea * localRowE( 0 ) * pow( -1, kf );
                            dRdU( 2, kf*3*numNodesPerFace + 3*a+i ) = nodalArea * localRowE( 1 ) * pow( -1, kf );
                          }
                        }
                      }
                      for( localIndex i=1 ; i<3 ; ++i )
                      {
                        dRdT( i, 0 ) = Ja * std::tan( m_frictionAngle ) * sliding( i-1 ) / slidingNorm;
                        dRdT( i, i ) = Ja;
                      }
                    }
                    else
                    {
                      R1TensorT< 2 > vaux( traction[kfe][1], traction[kfe][2] );
                      real64 vauxNorm = sqrt( vaux( 0 )*vaux( 0 ) + vaux( 1 )*vaux( 1 ) );
                      if( vauxNorm > 0.0 )
                      {
                        for( localIndex i=1 ; i<3 ; ++i )
                        {
                          elemRHS[i] = +Ja * ( traction[kfe][i] - limitTau * vaux( i-1 ) / vauxNorm );
                        }
                        for( localIndex i=1 ; i<3 ; ++i )
                        {
                          dRdT( i, i ) = Ja;
                        }
                      }
                      else
                      {
                        for( localIndex i=1 ; i<3 ; ++i )
                        {
                          elemRHS[i] = 0.0;
                        }
                        for( localIndex i=1 ; i<3 ; ++i )
                        {
                          dRdT( i, i ) = Ja;
                        }
                      }
                    }
                    break;
                  }
                case FractureState::OPEN:
                  {
                    for( localIndex i=0 ; i<3 ; ++i )
                    {
                      elemRHS[i] = +Ja * traction[kfe][i];
                    }

                    for( localIndex i=0 ; i<3 ; ++i )
                    {
                      dRdT( i, i ) = Ja;
                    }
                    break;
                  }
              }

//              for(int i=0;i<3;++i) GEOSX_LOG_RANK( "rhs[" << i << "]: "<<elemRHS[i]);

              rhs->add( elemDOF,
                        elemRHS,
                        3 );

              if( fractureState[kfe] != FractureState::OPEN )
              {
                matrix->add( elemDOF,
                             nodeDOF,
                             dRdU.data(),
                             3,
                             2*3*numNodesPerFace );
              }

              if( fractureState[kfe] != FractureState::STICK )
              {
                matrix->add( elemDOF,
                             elemDOF,
                             dRdT.data(),
                             3,
                             3 );
              }
            }
          }
        } );
    }
  } );
}

void LagrangianContactSolver::AssembleStabliziation( DomainPartition * const domain,
                                                     DofManager const & dofManager,
                                                     ParallelMatrix * const matrix,
                                                     ParallelVector * const rhs )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel * const mesh = domain->getMeshBody( 0 )->getMeshLevel( 0 );

  FaceManager const * const faceManager = mesh->getFaceManager();
  NodeManager const * const nodeManager = mesh->getNodeManager();
  ElementRegionManager * const elemManager = mesh->getElemManager();

  string const tracDofKey = dofManager.getKey( viewKeyStruct::tractionString );

  // Get the finite volume method used to compute the stabilization
  NumericalMethodsManager const * numericalMethodManager = domain->getParent()->GetGroup< NumericalMethodsManager >( keys::numericalMethodsManager );
  FiniteVolumeManager const * fvManager = numericalMethodManager->GetGroup< FiniteVolumeManager >( keys::finiteVolumeManager );
  FluxApproximationBase const * stabilizationMethod = fvManager->getFluxApproximation( m_stabilizationName );

  // Get the "face to element" map (valid for the entire mesh)
  FaceManager::ElemMapType const & faceToElem = faceManager->toElementRelation();

  // Form the SurfaceGenerator, get the fracture name and use it to retrieve the faceMap (from fracture element to face)
  SurfaceGenerator const * const
  surfaceGenerator = this->getParent()->GetGroup< SolverBase >( "SurfaceGen" )->group_cast< SurfaceGenerator const * >();
  FaceElementRegion * const fractureRegion = elemManager->GetRegion< FaceElementRegion >( surfaceGenerator->getFractureRegionName() );
  FaceElementSubRegion * const fractureSubRegion = fractureRegion->GetSubRegion< FaceElementSubRegion >( "default" );
  GEOSX_ERROR_IF( !fractureSubRegion->hasWrapper( m_tractionKey ), "The fracture subregion must contain traction field." );
  FaceElementSubRegion::FaceMapType const & faceMap = fractureSubRegion->faceList();
  GEOSX_ERROR_IF( faceMap.size( 1 ) != 2, "A fracture face has to be shared by two cells." );

  // Get the state of fracture elements
  arrayView1d< FractureState const > const &
  fractureState = fractureSubRegion->getReference< array1d< FractureState > >( viewKeyStruct::fractureStateString );

  // Get the tractions and stabilization contribution to the local jump
  arrayView2d< real64 const > const &
  traction = fractureSubRegion->getReference< array2d< real64 > >( viewKeyStruct::tractionString );
  arrayView2d< real64 > const &
  localJumpCorrection = fractureSubRegion->getReference< array2d< real64 > >( viewKeyStruct::localJumpCorrectionString );
  arrayView2d< real64 > const &
  previousLocalJumpCorrection = fractureSubRegion->getReference< array2d< real64 > >( viewKeyStruct::previousLocalJumpCorrectionString );

  // Get the volume for all elements
  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > const elemVolume =
    elemManager->ConstructViewAccessor< array1d< real64 >, arrayView1d< real64 const > >( ElementSubRegionBase::viewKeyStruct::elementVolumeString );

  // Get the coordinates for all nodes
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition = nodeManager->referencePosition();

  // Get area and rotation matrix for all faces
  arrayView1d< real64 const > const & faceArea = faceManager->faceArea();
  arrayView1d< R2Tensor const > const & faceRotationMatrix = faceManager->faceRotationMatrix();

  ConstitutiveManager const * const constitutiveManager = domain->GetGroup< ConstitutiveManager >( keys::ConstitutiveManager );
  ElementRegionManager::ConstitutiveRelationAccessor< ConstitutiveBase const > const
  constitutiveRelations = elemManager->ConstructFullConstitutiveAccessor< ConstitutiveBase const >( constitutiveManager );

  arrayView1d< globalIndex const > const &
  tracDofNumber = fractureSubRegion->getReference< array1d< globalIndex > >( tracDofKey );

  stabilizationMethod->forStencils< FaceElementStencil >( [&]( FaceElementStencil const & stencil )
  {
    for( localIndex iconn=0 ; iconn<stencil.size() ; ++iconn )
    {
      localIndex const numFluxElems = stencil.stencilSize( iconn );
      GEOSX_ERROR_IF( numFluxElems != 2, "A fracture connector has to be an edge shared by two faces." );
      typename FaceElementStencil::IndexContainerViewConstType const & sei = stencil.getElementIndices();

      // First index: face element. Second index: node
      real64_array2d nodalArea( 2, 2 );
      array1d< R2Tensor > rotatedInvStiffApprox( 2 );
      for( localIndex kf = 0 ; kf < 2 ; ++kf )
      {
        // Get fracture, face and region/subregion/element indices (for elements on both sides)
        localIndex fractureIndex = sei[iconn][kf];

        localIndex faceIndexRef = faceMap[fractureIndex][0];
        real64 const area = faceArea[faceIndexRef];
        R2Tensor const & rotationMatrix = faceRotationMatrix[faceIndexRef];
        // TODO: use higher order integration scheme
        nodalArea[kf][0] = area / 4.0;
        nodalArea[kf][1] = area / 4.0;

        real64_array2d invStiffApprox( 2, 3 );
        for( localIndex i = 0 ; i < 2 ; ++i )
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
          for( localIndex j = 0 ; j < 3 ; ++j )
          {
            maxSize[j] = nodePosition[cellElemsToNodes[ei][0]][j];
            minSize[j] = nodePosition[cellElemsToNodes[ei][0]][j];
          }
          for( localIndex a=1 ; a<numNodesPerElem ; ++a )
          {
            for( localIndex j = 0 ; j < 3 ; ++j )
            {
              maxSize[j] = std::max( maxSize[j], nodePosition[cellElemsToNodes[ei][a]][j] );
              minSize[j] = std::min( minSize[j], nodePosition[cellElemsToNodes[ei][a]][j] );
            }
          }
          real64_array boxSize( 3 );
          for( localIndex j = 0 ; j < 3 ; ++j )
          {
            boxSize[j] = maxSize[j] - minSize[j];
          }

          // Get linear elastic isotropic constitutive parameters for the element
          LinearElasticIsotropic const * const constitutiveRelation0 =
            dynamic_cast< LinearElasticIsotropic const * >( constitutiveRelations[er][esr][m_solidSolver->getSolidMaterialFullIndex()] );
          real64 const K = constitutiveRelation0->bulkModulus()[ei];
          real64 const G = constitutiveRelation0->shearModulus()[ei];
          real64 const E = 9.0 * K * G / ( 3.0 * K + G );
          real64 const nu = ( 3.0 * K - 2.0 * G ) / ( 2.0 * ( 3.0 * K + G ) );

          for( localIndex j = 0 ; j < 3 ; ++j )
          {
            invStiffApprox[i][j] = 1.0 / ( E / ( ( 1.0 + nu )*( 1.0 - 2.0*nu ) ) * 4.0 / 9.0 * ( 2.0 - 3.0 * nu ) * volume / ( boxSize[j]*boxSize[j] ) );
          }
        }

        R2Tensor invStiffApproxTotal;
        invStiffApproxTotal = 0.0;
        for( localIndex i = 0 ; i < 2 ; ++i )
        {
          for( localIndex j = 0 ; j < 3 ; ++j )
          {
            invStiffApproxTotal( j, j ) += invStiffApprox[i][j];
          }
        }
        // Compute R^T * (invK) * R
        R2Tensor tmpTensor;
        tmpTensor.AjiBjk( rotationMatrix, invStiffApproxTotal );
        rotatedInvStiffApprox[kf].AijBjk( tmpTensor, rotationMatrix );
      }

      // Compose local nodal-based local stiffness matrices
      stackArray2d< real64, 3*3 > totalInvStiffApprox( 3, 3 );
      for( localIndex kf = 0 ; kf < 2 ; ++kf )
      {
        for( localIndex i = 0 ; i < 3 ; ++i )
        {
          for( localIndex j = 0 ; j < 3 ; ++j )
          {
            rotatedInvStiffApprox[kf]( i, j ) *= nodalArea[0][kf] * nodalArea[1][kf];
          }
        }
      }
      // Local assembly
      for( localIndex i = 0 ; i < 3 ; ++i )
      {
        for( localIndex j = 0 ; j < 3 ; ++j )
        {
          totalInvStiffApprox( i, j ) = -( rotatedInvStiffApprox[0]( i, j ) + rotatedInvStiffApprox[1]( i, j ) );
        }
      }

      // Get DOF numbering
      localIndex fractureIndex[2];
      localIndex nDof[2];
      globalIndex elemDOF[2][3];
      for( localIndex kf = 0 ; kf < 2 ; ++kf )
      {
        fractureIndex[kf] = sei[iconn][kf];
        for( localIndex i=0 ; i<3 ; ++i )
        {
          elemDOF[kf][i] = -1;
        }
        nDof[kf] = 0;
        switch( fractureState[fractureIndex[kf]] )
        {
          case ( FractureState::STICK ):
            {
              nDof[kf] = 3;
              for( localIndex i=0 ; i<3 ; ++i )
              {
                elemDOF[kf][i] = tracDofNumber[fractureIndex[kf]] + integer_conversion< globalIndex >( i );
              }
              break;
            }
          case ( FractureState::NEW_SLIP ):
          case ( FractureState::SLIP ):
            {
              nDof[kf] = 1;
              elemDOF[kf][0] = tracDofNumber[fractureIndex[kf]] + integer_conversion< globalIndex >( 0 );
              break;
            }
          case ( FractureState::OPEN ):
            {
              nDof[kf] = 0;
              break;
            }
        }
      }

      // No interaction between adjacent fracture elements
      if( std::min( nDof[0], nDof[1] ) == 0 )
      {
        nDof[0] = 0;
        nDof[1] = 0;
      }

      // Global matrix assembly
      stackArray2d< real64, 3*3 > totalInvStiffApprox00( 3, 3 );
      totalInvStiffApprox00 = totalInvStiffApprox;
      totalInvStiffApprox00.resizeDimension< 0 >( nDof[0] );
      totalInvStiffApprox00.resizeDimension< 1 >( nDof[0] );
      matrix->add( elemDOF[0],
                   elemDOF[0],
                   totalInvStiffApprox00.data(),
                   nDof[0],
                   nDof[0] );

      stackArray2d< real64, 3*3 > totalInvStiffApprox11( 3, 3 );
      totalInvStiffApprox11 = totalInvStiffApprox;
      totalInvStiffApprox11.resizeDimension< 0 >( nDof[1] );
      totalInvStiffApprox11.resizeDimension< 1 >( nDof[1] );
      matrix->add( elemDOF[1],
                   elemDOF[1],
                   totalInvStiffApprox11.data(),
                   nDof[1],
                   nDof[1] );

      // Change sign
      for( localIndex i = 0 ; i < 3 ; ++i )
      {
        for( localIndex j = 0 ; j < 3 ; ++j )
        {
          totalInvStiffApprox( i, j ) *= -1.0;
        }
      }

      stackArray2d< real64, 3*3 > totalInvStiffApprox01( 3, 3 );
      totalInvStiffApprox01 = totalInvStiffApprox;
      totalInvStiffApprox01.resizeDimension< 0 >( nDof[0] );
      totalInvStiffApprox01.resizeDimension< 1 >( nDof[1] );
      matrix->add( elemDOF[0],
                   elemDOF[1],
                   totalInvStiffApprox01.data(),
                   nDof[0],
                   nDof[1] );

      stackArray2d< real64, 3*3 > totalInvStiffApprox10( 3, 3 );
      totalInvStiffApprox10 = totalInvStiffApprox;
      totalInvStiffApprox10.resizeDimension< 0 >( nDof[1] );
      totalInvStiffApprox10.resizeDimension< 1 >( nDof[0] );
      matrix->add( elemDOF[1],
                   elemDOF[0],
                   totalInvStiffApprox10.data(),
                   nDof[1],
                   nDof[0] );

      // Compute rhs
      real64_array rhs0( 3 );
      rhs0 = 0.0;
      for( localIndex i = 0 ; i < nDof[0] ; ++i )
      {
        for( localIndex j = 0 ; j < nDof[0] ; ++j )
        {
          rhs0( i ) += totalInvStiffApprox00( i, j ) * ( traction[fractureIndex[0]][j] );
        }
        for( localIndex j = 0 ; j < nDof[1] ; ++j )
        {
          rhs0( i ) += totalInvStiffApprox01( i, j ) * ( traction[fractureIndex[1]][j] );
        }
      }

      real64_array rhs1( 3 );
      rhs1 = 0.0;
      for( localIndex i = 0 ; i < nDof[1] ; ++i )
      {
        for( localIndex j = 0 ; j < nDof[0] ; ++j )
        {
          rhs1( i ) += totalInvStiffApprox10( i, j ) * ( traction[fractureIndex[0]][j] );
        }
        for( localIndex j = 0 ; j < nDof[1] ; ++j )
        {
          rhs1( i ) += totalInvStiffApprox11( i, j ) * ( traction[fractureIndex[1]][j] );
        }
      }

      // Save rhs in previousLocalJumpCorretion array
      for( localIndex i = 0 ; i < nDof[0] ; ++i )
      {
        localJumpCorrection[fractureIndex[0]][i] = rhs0( i );
      }
      for( localIndex i = 0 ; i < nDof[1] ; ++i )
      {
        localJumpCorrection[fractureIndex[1]][i] = rhs1( i );
      }

      // Compute the difference of rhs
      if( fractureState[fractureIndex[0]] == FractureState::NEW_SLIP || fractureState[fractureIndex[0]] == FractureState::SLIP )
      {
        for( localIndex i = 0 ; i < nDof[0] ; ++i )
        {
          if( i > 0 )
          {
            rhs0( i ) -= previousLocalJumpCorrection[fractureIndex[0]][i];
          }
        }
      }
      if( fractureState[fractureIndex[1]] == FractureState::NEW_SLIP || fractureState[fractureIndex[1]] == FractureState::SLIP )
      {
        for( localIndex i = 0 ; i < nDof[1] ; ++i )
        {
          if( i > 0 )
          {
            rhs1( i ) -= previousLocalJumpCorrection[fractureIndex[1]][i];
          }
        }
      }

//      GEOSX_LOG_RANK( "pair: " << tracDofNumber[fractureIndex[0]] << " and " << tracDofNumber[fractureIndex[1]] );
//      GEOSX_LOG_RANK( traction[fractureIndex[0]] << " " << traction[fractureIndex[1]] );
//      GEOSX_LOG_RANK( elemDOF[0][0] << " " << elemDOF[0][1] << " " << elemDOF[0][2] );
//      GEOSX_LOG_RANK( elemDOF[1][0] << " " << elemDOF[1][1] << " " << elemDOF[1][2] );
//      GEOSX_LOG_RANK( rhs0 << " " << rhs1 );

      // Global rhs assembly
      rhs->add( elemDOF[0],
                rhs0.data(),
                nDof[0] );
      rhs->add( elemDOF[1],
                rhs1.data(),
                nDof[1] );
    }
  } );
}

void LagrangianContactSolver::ApplySystemSolution( DofManager const & dofManager,
                                                   ParallelVector const & solution,
                                                   real64 const scalingFactor,
                                                   DomainPartition * const domain )
{
  GEOSX_MARK_FUNCTION;
  GEOSX_LOG_RANK( "In ApplySystemSolution!" );

  m_solidSolver->ApplySystemSolution( dofManager,
                                      solution,
                                      scalingFactor,
                                      domain );

  dofManager.addVectorToField( solution, viewKeyStruct::tractionString, viewKeyStruct::deltaTractionString, -scalingFactor );
  dofManager.addVectorToField( solution, viewKeyStruct::tractionString, viewKeyStruct::tractionString, -scalingFactor );

  std::map< string, string_array > fieldNames;
  fieldNames["elems"].push_back( viewKeyStruct::tractionString );
  fieldNames["elems"].push_back( viewKeyStruct::deltaTractionString );
  fieldNames["elems"].push_back( viewKeyStruct::fractureStateString );
  fieldNames["elems"].push_back( viewKeyStruct::localJumpString );
  //  fieldNames["elems"].push_back( viewKeyStruct::previousLocalJumpString );
//  fieldNames["elems"].push_back( viewKeyStruct::localJumpCorrectionString );
//  fieldNames["elems"].push_back( viewKeyStruct::previousLocalJumpCorrectionString );

  CommunicationTools::SynchronizeFields( fieldNames,
                                         domain->getMeshBody( 0 )->getMeshLevel( 0 ),
                                         domain->getNeighbors() );

  this->UpdateDeformationForCoupling( domain );
}

void LagrangianContactSolver::InitializeFractureState( MeshLevel * const mesh )
{
  GEOSX_MARK_FUNCTION;
  ElementRegionManager * const elemManager = mesh->getElemManager();

  elemManager->forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion * const subRegion )->void
  {
    if( subRegion->hasWrapper( m_tractionKey ) )
    {
      arrayView1d< FractureState > const & fractureState = subRegion->getReference< array1d< FractureState > >( viewKeyStruct::fractureStateString );
      fractureState = FractureState::STICK;
    }
  } );
}

bool LagrangianContactSolver::UpdateFractureState( DomainPartition * const domain )
{
  GEOSX_MARK_FUNCTION;
  GEOSX_LOG_RANK( "In UpdateFractureState!" );

  MeshLevel * const mesh = domain->getMeshBody( 0 )->getMeshLevel( 0 );
  ElementRegionManager * const elemManager = mesh->getElemManager();

  bool checkActiveSet = true;

  elemManager->forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion * const subRegion )->void
  {
    if( subRegion->hasWrapper( m_tractionKey ) )
    {
      arrayView1d< integer const > const & ghostRank = subRegion->GhostRank();
      arrayView2d< real64 const > const & traction = subRegion->getReference< array2d< real64 > >( viewKeyStruct::tractionString );
      arrayView2d< real64 const > const & localJump = subRegion->getReference< array2d< real64 > >( viewKeyStruct::localJumpString );
      arrayView1d< FractureState > const & fractureState = subRegion->getReference< array1d< FractureState > >( viewKeyStruct::fractureStateString );

      forall_in_range< serialPolicy >( 0,
                                       subRegion->size(),
                                       [&]( localIndex const kfe )
      {
        {
          if( ghostRank[kfe] < 0 )
          {
            FractureState const originalFractureState = fractureState[kfe];
            if( originalFractureState == FractureState::OPEN )
            {
              if( localJump[kfe][0] > -m_normalDisplacementTolerance )
              {
                fractureState[kfe] = FractureState::OPEN;
              }
              else
              {
                fractureState[kfe] = FractureState::STICK;
              }
            }
            else if( traction[kfe][0] > m_normalTractionTolerance )
            {
              fractureState[kfe] = FractureState::OPEN;
            }
            else
            {
              real64 currentTau = sqrt( traction[kfe][1]*traction[kfe][1] + traction[kfe][2]*traction[kfe][2] );
              real64 limitTau = m_cohesion - traction[kfe][0] * std::tan( m_frictionAngle );
              if( originalFractureState == FractureState::STICK && currentTau >= limitTau )
              {
                currentTau *= (1.0 - m_alpha);
              }
              else if( originalFractureState != FractureState::STICK && currentTau <= limitTau )
              {
                currentTau *= (1.0 + m_alpha);
              }
              if( currentTau > limitTau )
              {
                if( originalFractureState == FractureState::STICK )
                {
                  fractureState[kfe] = FractureState::NEW_SLIP;
                }
                else
                {
                  fractureState[kfe] = FractureState::SLIP;
                }
              }
              else
              {
                fractureState[kfe] = FractureState::STICK;
              }
            }
            GEOSX_LOG_RANK( "element " << kfe << " traction: " << traction[kfe]
                                       << " previous state <"
                                       << FractureStateToString( originalFractureState )
                                       << "> current state <"
                                       << FractureStateToString( fractureState[kfe] )
                                       << ">" );
            checkActiveSet &= CompareFractureStates( originalFractureState, fractureState[kfe] );
          }
        }
      } );
    }
  } );

  bool globalCheckActiveSet;
  MpiWrapper::allReduce( &checkActiveSet,
                         &globalCheckActiveSet,
                         1,
                         MPI_LAND,
                         MPI_COMM_GEOSX );

  return globalCheckActiveSet;
}

void LagrangianContactSolver::SolveSystem( DofManager const & dofManager,
                                           ParallelMatrix & matrix,
                                           ParallelVector & rhs,
                                           ParallelVector & solution )
{
  GEOSX_MARK_FUNCTION;
  GEOSX_LOG_RANK( "In SolveSystem!" );

  if( getLogLevel() >= 2 )
  {
    matrix.write( "matrix" );
    rhs.write( "rhs" );
  }

  SolverBase::SolveSystem( dofManager, matrix, rhs, solution );

  if( getLogLevel() >= 2 )
  {
    solution.write( "sol" );
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

void LagrangianContactSolver::SetNextDt( real64 const & currentDt,
                                         real64 & nextDt )
{
  nextDt = currentDt;
}

REGISTER_CATALOG_ENTRY( SolverBase, LagrangianContactSolver, std::string const &, Group * const )
} /* namespace geosx */
