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
#include "physicsSolvers/fluidFlow/FlowSolverBase.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"
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
  m_maxNumResolves( 10 )
{
  registerWrapper( viewKeyStruct::solidSolverNameString, &m_solidSolverName, 0 )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Name of the solid mechanics solver to use in the poroelastic solver" );

  m_numResolves[0] = 0;
}

void LagrangianContactSolver::RegisterDataOnMesh( dataRepository::Group * const MeshBodies )
{
  std::cout << "In RegisterDataOnMesh!" << std::endl;
  for( auto & mesh : MeshBodies->GetSubGroups() )
  {
    MeshLevel * meshLevel = Group::group_cast< MeshBody * >( mesh.second )->getMeshLevel( 0 );

    ElementRegionManager * const elemManager = meshLevel->getElemManager();
    elemManager->forElementRegions< FaceElementRegion >( [&] ( FaceElementRegion * const region )
    {
      region->forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion * const subRegion )
      {
        subRegion->registerWrapper< array2d< real64 > >( viewKeyStruct::tractionString )->
          setApplyDefaultValue(0.0)->
          setPlotLevel( PlotLevel::LEVEL_0 )->
          setRegisteringObjects( this->getName())->
          setDescription( "An array that holds the tractions on the fracture." );
        subRegion->getReference< array2d< real64 > >( viewKeyStruct::tractionString ).resizeDimension< 1 >( 3 );

        subRegion->registerWrapper< array2d< real64 > >( viewKeyStruct::deltaTractionString )->
          setApplyDefaultValue(0.0)->
          setPlotLevel( PlotLevel::LEVEL_0 )->
          setRegisteringObjects( this->getName())->
          setDescription( "An array that holds the traction increments on the fracture." );
        subRegion->getReference< array2d< real64 > >( viewKeyStruct::deltaTractionString ).resizeDimension< 1 >( 3 );

        subRegion->registerWrapper< array1d< FractureState > >( viewKeyStruct::fractureStateString )->
          setPlotLevel( PlotLevel::LEVEL_1 )->
          setRegisteringObjects( this->getName())->
          setDescription( "An array that holds the fracture state." );
        InitializeFractureState( meshLevel );

//        // TODO: maybe not needed???
//        subRegion->registerWrapper< array1d< R1Tensor > >( FaceElementSubRegion::viewKeyStruct::elementLocalJumpString )->
//          setPlotLevel( PlotLevel::LEVEL_1 )->
//          setRegisteringObjects( this->getName())->
//          setDescription( "An array that holds the local jump on the fracture." );
      } );
    } );
  }
}

void LagrangianContactSolver::ImplicitStepSetup( real64 const & time_n,
                                                 real64 const & dt,
                                                 DomainPartition * const domain,
                                                 DofManager & GEOSX_UNUSED_PARAM( dofManager ),
                                                 ParallelMatrix & GEOSX_UNUSED_PARAM( matrix ),
                                                 ParallelVector & GEOSX_UNUSED_PARAM( rhs ),
                                                 ParallelVector & GEOSX_UNUSED_PARAM( solution ) )
{
  this->UpdateDeformationForCoupling( domain );

  m_solidSolver->ImplicitStepSetup( time_n, dt, domain,
                                    m_solidSolver->getDofManager(),
                                    m_solidSolver->getSystemMatrix(),
                                    m_solidSolver->getSystemRhs(),
                                    m_solidSolver->getSystemSolution() );

#ifdef GEOSX_USE_SEPARATION_COEFFICIENT
  domain->getMeshBody( 0 )->getMeshLevel( 0 )->getElemManager()->forElementRegions< FaceElementRegion >( [&]( FaceElementRegion * const faceElemRegion )
  {
    faceElemRegion->forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion * const subRegion )
    {
      arrayView1d< real64 > const &
      separationCoeff0 = subRegion->getReference< array1d< real64 > >( viewKeyStruct::separationCoeff0String );
      arrayView1d< real64 const > const &
      separationCoeff = subRegion->getSeparationCoefficient();
      for( localIndex k=0 ; k<separationCoeff0.size() ; ++k )
      {
        separationCoeff0[k] = separationCoeff[k];
      }
    } );
  } );
#endif

}

void LagrangianContactSolver::ImplicitStepComplete( real64 const & time_n,
                                                    real64 const & dt,
                                                    DomainPartition * const domain )
{
  m_solidSolver->ImplicitStepComplete( time_n, dt, domain );
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
  std::cout << "In ResetStateToBeginningOfStep!" << std::endl;

  m_solidSolver->ResetStateToBeginningOfStep( domain );

  string const tractionKey = viewKeyStruct::tractionString;
  MeshLevel * const meshLevel = domain->getMeshBody( 0 )->getMeshLevel( 0 );
  ElementRegionManager * const elemManager = meshLevel->getElemManager();

  elemManager->forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion * const subRegion )->void
  {
    if( subRegion->hasWrapper( tractionKey ) )
    {
      arrayView2d< real64 > const & traction = subRegion->getReference< array2d< real64 > >( viewKeyStruct::tractionString );
      arrayView2d< real64 > const & deltaTraction = subRegion->getReference< array2d< real64 > >( viewKeyStruct::deltaTractionString );

      forall_in_range< serialPolicy >( 0,
                                       subRegion->size(),
                                       GEOSX_LAMBDA ( localIndex const kfe )
        {
          for( int i = 0 ; i < 3 ; ++i )
          {
            traction[kfe]( i ) -= deltaTraction[kfe]( i );
            deltaTraction[kfe]( i ) = 0.0;
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

  m_numResolves[1] = m_numResolves[0];

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
  m_numResolves[1] += 1;

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

  string const tractionKey = viewKeyStruct::tractionString;

  elemManager->forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion * const subRegion )->void
  {
    if( subRegion->hasWrapper( tractionKey ) )
    {
      arrayView1d< integer const > const & ghostRank = subRegion->GhostRank();
      arrayView2d< real64 > const & localJump = subRegion->getElementLocalJump();
      arrayView1d< R2Tensor const > const & rotationMatrix = subRegion->getElementRotationMatrix();
      arrayView2d< localIndex const > const & elemsToFaces = subRegion->faceList();

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
                for( int i=0 ; i<3 ; ++i )
                {
                  globalJump( i ) +=
                    ( - u[faceToNodeMap( elemsToFaces[kfe][0], a )][i]
                      + u[faceToNodeMap( elemsToFaces[kfe][1], a )][i] ) / static_cast< real64 >(numNodesPerFace);
                }
              }
              R1Tensor localJumpTensor;
              localJumpTensor.AijBi( rotationMatrix[kfe], globalJump );
              localJump[kfe][0] = localJumpTensor[0];
              localJump[kfe][1] = localJumpTensor[1];
              localJump[kfe][2] = localJumpTensor[2];

              std::cout << "Element: " << kfe << " localJump: " << localJumpTensor << std::endl;
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

real64 LagrangianContactSolver::ExplicitStep( real64 const & time_n,
                                              real64 const & dt,
                                              const int cycleNumber,
                                              DomainPartition * const domain )
{
  GEOSX_MARK_FUNCTION;
  m_solidSolver->ExplicitStep( time_n, dt, cycleNumber, domain );
//  m_flowSolver->SolverStep( time_n, dt, cycleNumber, domain );

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
  real64 stepDt;
  bool continueLoop = true;
  while( continueLoop )
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
    std::cout << isPreviousFractureStateValid << std::endl;
    if( isPreviousFractureStateValid )
    {
      continueLoop = false;
    }
    else
    {
      std::cout << "*******************************************" << std::endl;
      std::cout << "*******************************************" << std::endl;
      std::cout << "*******************************************" << std::endl;
    }
  }

  return stepDt;
}

void LagrangianContactSolver::SetupDofs( DomainPartition const * const domain,
                                         DofManager & dofManager ) const
{
  GEOSX_MARK_FUNCTION;
  m_solidSolver->SetupDofs( domain, dofManager );

  // restrict coupling to fracture regions only (as done originally in SetupSystem)
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
                          DofManager::Connectivity::Face );

  dofManager.addCoupling( keys::TotalDisplacement,
                          viewKeyStruct::tractionString,
                          DofManager::Connectivity::Elem,
                          fractureRegions );
}

void LagrangianContactSolver::SetupSystem( DomainPartition * const domain,
                                           DofManager & dofManager,
                                           ParallelMatrix & matrix,
                                           ParallelVector & rhs,
                                           ParallelVector &  solution )
{
  GEOSX_MARK_FUNCTION;

  m_solidSolver->SetupSystem( domain,
                              m_solidSolver->getDofManager(),
                              m_solidSolver->getSystemMatrix(),
                              m_solidSolver->getSystemRhs(),
                              m_solidSolver->getSystemSolution() );
  // setup DofManager
  dofManager.setMesh( domain, 0, 0 );

  // add traction and coupling
  SetupDofs( domain, dofManager );

  // By not calling dofManager.reorderByRank(), we keep separate dof numbering for each field,
  // which allows constructing separate sparsity patterns for off-diagonal blocks of the matrix.
  // Once the solver moves to monolithic matrix, we can remove this method and just use SolverBase::SetupSystem.
  localIndex const numDisplacementDofs = dofManager.numLocalDofs( keys::TotalDisplacement );
  localIndex const numTractionDofs = dofManager.numLocalDofs( viewKeyStruct::tractionString );
  m_matrix01.createWithLocalSize( numDisplacementDofs,
                                  numTractionDofs,
                                  12,
                                  MPI_COMM_GEOSX );
  m_matrix10.createWithLocalSize( numTractionDofs,
                                  numDisplacementDofs,
                                  24,
                                  MPI_COMM_GEOSX );
  m_matrix11.createWithLocalSize( numTractionDofs,
                                  numTractionDofs,
                                  15,
                                  MPI_COMM_GEOSX );

  m_rhs0.createWithLocalSize( numDisplacementDofs,
                              MPI_COMM_GEOSX );
  m_rhs1.createWithLocalSize( numTractionDofs,
                              MPI_COMM_GEOSX );

  dofManager.setSparsityPattern( m_matrix01, keys::TotalDisplacement, viewKeyStruct::tractionString );
  dofManager.setSparsityPattern( m_matrix10, viewKeyStruct::tractionString, keys::TotalDisplacement );
  dofManager.setSparsityPattern( m_matrix11, viewKeyStruct::tractionString, viewKeyStruct::tractionString );

  matrix.createWithLocalSize( m_solidSolver->getSystemMatrix().localRows() + m_matrix11.localRows(),
                              m_solidSolver->getSystemMatrix().localCols() + m_matrix11.localCols(),
                              MPI_COMM_GEOSX );

  globalIndex numDispDofs = m_solidSolver->getSystemRhs().globalSize();
  matrix.open();
  for( globalIndex i=m_solidSolver->getSystemMatrix().ilower() ; i<m_solidSolver->getSystemMatrix().iupper() ; ++i )
  {
    globalIndex_array colIndices;
    real64_array values;
    m_solidSolver->getSystemMatrix().getRowCopy( i, colIndices, values );
    matrix.insert( i, colIndices, values );

    m_matrix01.getRowCopy( i, colIndices, values );
    for( localIndex j = 0 ; j < colIndices.size() ; ++j )
    {
      colIndices[j] += numDispDofs;
    }
    matrix.insert( i, colIndices, values );
  }
  for( globalIndex i=m_matrix11.ilower() ; i<m_matrix11.iupper() ; ++i )
  {
    globalIndex_array colIndices;
    real64_array values;
    m_matrix10.getRowCopy( i, colIndices, values );
    matrix.insert( numDispDofs+i, colIndices, values );

    m_matrix11.getRowCopy( i, colIndices, values );
    for( localIndex j = 0 ; j < colIndices.size() ; ++j )
    {
      colIndices[j] += numDispDofs;
    }
    matrix.insert( numDispDofs+i, colIndices, values );
  }
  matrix.close();

  rhs.createWithLocalSize( m_solidSolver->getSystemRhs().localSize() + m_rhs1.localSize(),
                           MPI_COMM_GEOSX );

  solution.createWithLocalSize( m_solidSolver->getSystemSolution().localSize() + m_matrix11.localCols(),
                                MPI_COMM_GEOSX );
}

void LagrangianContactSolver::AssembleSystem( real64 const time,
                                              real64 const dt,
                                              DomainPartition * const domain,
                                              DofManager const & GEOSX_UNUSED_PARAM( dofManager ),
                                              ParallelMatrix & GEOSX_UNUSED_PARAM( matrix ),
                                              ParallelVector & GEOSX_UNUSED_PARAM( rhs ) )
{
  GEOSX_MARK_FUNCTION;

  m_solidSolver->AssembleSystem( time,
                                 dt,
                                 domain,
                                 m_solidSolver->getDofManager(),
                                 m_solidSolver->getSystemMatrix(),
                                 m_solidSolver->getSystemRhs() );

  AssembleForceResidualDerivativeWrtTraction( domain, &m_matrix01, &m_rhs0 );
  AssembleTractionResidualDerivativeWrtDisplacementAndTraction( domain, &m_matrix10, &m_matrix11, &m_rhs1 );
}

void LagrangianContactSolver::ApplyBoundaryConditions( real64 const time,
                                                       real64 const dt,
                                                       DomainPartition * const domain,
                                                       DofManager const & GEOSX_UNUSED_PARAM( dofManager ),
                                                       ParallelMatrix & GEOSX_UNUSED_PARAM( matrix ),
                                                       ParallelVector & GEOSX_UNUSED_PARAM( rhs ) )
{
  GEOSX_MARK_FUNCTION;
  m_solidSolver->ApplyBoundaryConditions( time,
                                          dt,
                                          domain,
                                          m_solidSolver->getDofManager(),
                                          m_solidSolver->getSystemMatrix(),
                                          m_solidSolver->getSystemRhs() );

  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup< MeshBody >( 0 )->getMeshLevel( 0 );

  FieldSpecificationManager const & fsManager = FieldSpecificationManager::get();
  string const dispDofKey = m_solidSolver->getDofManager().getKey( keys::TotalDisplacement );
  NodeManager const * const nodeManager = mesh->getNodeManager();
  arrayView1d< globalIndex const > const & dispDofNumber = nodeManager->getReference< globalIndex_array >( dispDofKey );
  arrayView1d< integer const > const & nodeGhostRank = nodeManager->GhostRank();

  fsManager.Apply( time + dt,
                   domain,
                   "nodeManager",
                   keys::TotalDisplacement,
                   [&]( FieldSpecificationBase const * const bc,
                        string const &,
                        SortedArrayView< localIndex const > const & targetSet,
                        Group * const,
                        string const )
  {
    SortedArray< localIndex > localSet;
    for( auto const & a : targetSet )
    {
      if( nodeGhostRank[a]<0 )
      {
        localSet.insert( a );
      }
    }
    bc->ZeroSystemRowsForBoundaryCondition< LAInterface >( localSet,
                                                           dispDofNumber,
                                                           m_matrix01 );
  } );
}

real64
LagrangianContactSolver::
  CalculateResidualNorm( DomainPartition const * const GEOSX_UNUSED_PARAM( domain ),
                         DofManager const & GEOSX_UNUSED_PARAM( dofManager ),
                         ParallelVector const & GEOSX_UNUSED_PARAM( rhs ) )
{
  GEOSX_MARK_FUNCTION;

  std::cout << "In CalculateResidualNorm!" << std::endl;
  ParallelVector const & solidResidual = m_solidSolver->getSystemRhs();
  real64 const tractionResidualNorm = m_rhs1.norm2();

  ParallelVector rhs0;
  rhs0.create( solidResidual );
  rhs0.copy( solidResidual );
  rhs0.axpy( -1.0, m_rhs0 );
  real64 const solidResidualNorm = rhs0.norm2();

  char output[60] = {0};
  sprintf( output,
           "( Rsolid, Rtraction ) = (%15.6e, %15.6e);",
           solidResidualNorm,
           tractionResidualNorm );
  GEOSX_LOG_RANK_0( output );

  return solidResidualNorm + tractionResidualNorm;
}

void
LagrangianContactSolver::
  AssembleForceResidualDerivativeWrtTraction( DomainPartition * const domain,
                                              ParallelMatrix * const matrix01,
                                              ParallelVector * const rhs0 )
{
  GEOSX_MARK_FUNCTION;
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup< MeshBody >( 0 )->getMeshLevel( 0 );

  FaceManager const * const faceManager = mesh->getFaceManager();
  NodeManager * const nodeManager = mesh->getNodeManager();
  ElementRegionManager * const elemManager = mesh->getElemManager();

  ArrayOfArraysView< localIndex const > const & faceToNodeMap = faceManager->nodeList();

  arrayView1d< R1Tensor > const &
  fext = nodeManager->getReference< array1d< R1Tensor > >( SolidMechanicsLagrangianFEM::viewKeyStruct::forceExternal );
  fext = {0, 0, 0};

  string const tractionKey = viewKeyStruct::tractionString;
  string const tracDofKey = getDofManager().getKey( viewKeyStruct::tractionString );
  string const dispDofKey = m_solidSolver->getDofManager().getKey( keys::TotalDisplacement );

  arrayView1d< globalIndex > const &
  dispDofNumber = nodeManager->getReference< globalIndex_array >( dispDofKey );

  matrix01->open();
  matrix01->zero();
  rhs0->open();
  rhs0->zero();

  elemManager->forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion const * const subRegion )->void
  {
    if( subRegion->hasWrapper( tractionKey ) )
    {
      arrayView1d< globalIndex const > const &
      tracDofNumber = subRegion->getReference< array1d< globalIndex > >( tracDofKey );
      arrayView2d< real64 const > const & traction = subRegion->getReference< array2d< real64 > >( viewKeyStruct::tractionString );
      arrayView1d< integer const > const & ghostRank = subRegion->GhostRank();
      arrayView1d< real64 const > const & area = subRegion->getElementArea();
      arrayView1d< R2Tensor const > const & rotationMatrix = subRegion->getElementRotationMatrix();
      arrayView2d< localIndex const > const & elemsToFaces = subRegion->faceList();

      std::cout << traction << std::endl;

      forall_in_range< serialPolicy >( 0,
                                       subRegion->size(),
                                       GEOSX_LAMBDA ( localIndex const kfe )
        {
          localIndex const kf0 = elemsToFaces[kfe][0];
          localIndex const numNodesPerFace = faceToNodeMap.sizeOfArray( kf0 );

          globalIndex rowDOF[12];
          real64 nodeRHS[12];
          stackArray2d< real64, 3*4*3 > dRdT( 3*numNodesPerFace, 3 );
          dRdT = 0.0;
          globalIndex colDOF[3];
          for( int i=0 ; i<3 ; ++i )
          {
            colDOF[i] = tracDofNumber[kfe] + i;
          }

          real64 const nodalArea = area[kfe] / static_cast< real64 >( numNodesPerFace );
          real64_array nodalForceVec( 3 );
          nodalForceVec( 0 ) = ( traction[kfe]( 0 ) ) * nodalArea;
          nodalForceVec( 1 ) = ( traction[kfe]( 1 ) ) * nodalArea;
          nodalForceVec( 2 ) = ( traction[kfe]( 2 ) ) * nodalArea;
          R1Tensor localNodalForce( nodalForceVec( 0 ), nodalForceVec( 1 ), nodalForceVec( 2 ));
          R1Tensor globalNodalForce;
          globalNodalForce.AijBj( rotationMatrix[kfe], localNodalForce );

          for( localIndex kf=0 ; kf<2 ; ++kf )
          {
            localIndex const faceIndex = elemsToFaces[kfe][kf];

            for( localIndex a=0 ; a<numNodesPerFace ; ++a )
            {
              for( int i=0 ; i<3 ; ++i )
              {
                rowDOF[3*a+i] = dispDofNumber[faceToNodeMap( faceIndex, a )] + i;
                nodeRHS[3*a+i] = -globalNodalForce[i] * pow( -1, kf );
                fext[faceToNodeMap( faceIndex, a )][i] += -globalNodalForce[i] * pow( -1, kf );

                dRdT( 3*a+i, 0 ) = - nodalArea * rotationMatrix[kfe]( i, 0 ) * pow( -1, kf );
                dRdT( 3*a+i, 1 ) = - nodalArea * rotationMatrix[kfe]( i, 1 ) * pow( -1, kf );
                dRdT( 3*a+i, 2 ) = - nodalArea * rotationMatrix[kfe]( i, 2 ) * pow( -1, kf );
              }
            }

            if( ghostRank[kfe] < 0 )
            {
              rhs0->add( rowDOF,
                         nodeRHS,
                         3 * numNodesPerFace );

              matrix01->add( rowDOF,
                             colDOF,
                             dRdT.data(),
                             3 * numNodesPerFace,
                             3 );
            }
          }
        } );
    }
  } );

  rhs0->close();
  matrix01->close();
}

void
LagrangianContactSolver::
  AssembleTractionResidualDerivativeWrtDisplacementAndTraction( DomainPartition const * const domain,
                                                                ParallelMatrix * const matrix10,
                                                                ParallelMatrix * const matrix11,
                                                                ParallelVector * const rhs1 )
{
  GEOSX_MARK_FUNCTION;
  MeshLevel const * const mesh = domain->getMeshBodies()->GetGroup< MeshBody >( 0 )->getMeshLevel( 0 );

  FaceManager const * const faceManager = mesh->getFaceManager();
  NodeManager const * const nodeManager = mesh->getNodeManager();
  ElementRegionManager const * const elemManager = mesh->getElemManager();

  ArrayOfArraysView< localIndex const > const & faceToNodeMap = faceManager->nodeList();

  string const tractionKey = viewKeyStruct::tractionString;
  string const tracDofKey = getDofManager().getKey( viewKeyStruct::tractionString );
  string const dispDofKey = m_solidSolver->getDofManager().getKey( keys::TotalDisplacement );

  arrayView1d< globalIndex const > const &
  dispDofNumber = nodeManager->getReference< array1d< globalIndex > >( dispDofKey );

  matrix10->open();
  matrix10->zero();
  matrix11->open();
  matrix11->zero();
  rhs1->open();
  rhs1->zero();

  elemManager->forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion const * const subRegion )->void
  {
    if( subRegion->hasWrapper( tractionKey ) )
    {
      arrayView1d< globalIndex const > const &
      tracDofNumber = subRegion->getReference< array1d< globalIndex > >( tracDofKey );
      arrayView1d< integer const > const & ghostRank = subRegion->GhostRank();
      arrayView1d< real64 const > const & area = subRegion->getElementArea();
      arrayView1d< R2Tensor const > const & rotationMatrix = subRegion->getElementRotationMatrix();
      arrayView2d< localIndex const > const & elemsToFaces = subRegion->faceList();
      arrayView2d< real64 const > const & localJump = subRegion->getElementLocalJump();
      arrayView2d< real64 const > const &
      traction = subRegion->getReference< array2d< real64 > >( viewKeyStruct::tractionString );
      arrayView1d< FractureState const > const &
      fractureState = subRegion->getReference< array1d< FractureState > >( viewKeyStruct::fractureStateString );

      forall_in_range< serialPolicy >( 0,
                                       subRegion->size(),
                                       GEOSX_LAMBDA ( localIndex const kfe )
        {
          {
            localIndex const numNodesPerFace = faceToNodeMap.sizeOfArray( elemsToFaces[kfe][0] );
            globalIndex nodeDOF[24];
            globalIndex elemDOF[3];
            for( int i=0 ; i<3 ; ++i )
            {
              elemDOF[i] = tracDofNumber[kfe] + i;
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
              switch( fractureState[kfe] )
              {
                case FractureState::STICK:
                {
                  std::cout << kfe << " STICK" << std::endl;
                  for( int i=0 ; i<3 ; ++i )
                  {
                    elemRHS[i] = + Ja * localJump[kfe][i];
                  }

                  for( localIndex kf=0 ; kf<2 ; ++kf )
                  {
                    for( localIndex a=0 ; a<numNodesPerFace ; ++a )
                    {
                      for( int i=0 ; i<3 ; ++i )
                      {
                        nodeDOF[ kf*3*numNodesPerFace + 3*a+i ] = dispDofNumber[faceToNodeMap( elemsToFaces[kfe][kf], a )] + i;

                        dRdU( 0, kf*3*numNodesPerFace + 3*a+i ) = - nodalArea * rotationMatrix[kfe]( 0, i ) * pow( -1, kf );
                        dRdU( 1, kf*3*numNodesPerFace + 3*a+i ) = - nodalArea * rotationMatrix[kfe]( 1, i ) * pow( -1, kf );
                        dRdU( 2, kf*3*numNodesPerFace + 3*a+i ) = - nodalArea * rotationMatrix[kfe]( 2, i ) * pow( -1, kf );
                      }
                    }
                  }
                  break;
                }
                case FractureState::SLIP:
                case FractureState::NEW_SLIP:
                {
                  std::cout << kfe << " SLIP" << std::endl;
                  elemRHS[0] = + Ja * localJump[kfe](0);

                  for( localIndex kf=0 ; kf<2 ; ++kf )
                  {
                    for( localIndex a=0 ; a<numNodesPerFace ; ++a )
                    {
                      for( int i=0 ; i<3 ; ++i )
                      {
                        nodeDOF[ kf*3*numNodesPerFace + 3*a+i ] = dispDofNumber[faceToNodeMap( elemsToFaces[kfe][kf], a )] + i;
                        dRdU( 0, kf*3*numNodesPerFace + 3*a+i ) = - nodalArea * rotationMatrix[kfe]( 0, i ) * pow( -1, kf );
                      }
                    }
                  }

                  real64 limitTau = m_cohesion - traction[kfe](0) * std::tan( m_frictionAngle );
                  R1TensorT<2> sliding( localJump[kfe](1), localJump[kfe](2) );
                  real64 slidingNorm = sqrt( sliding(0)*sliding(0) + sliding(1)*sliding(1) );

                  if( !( ( m_nonlinearSolverParameters.m_numNewtonIterations == 0 ) && ( fractureState[kfe] == FractureState::NEW_SLIP ) )
                      && slidingNorm > m_slidingTolerance )
                  {
                    for( int i=1 ; i<3 ; ++i )
                    {
                      elemRHS[i] = + Ja * ( traction[kfe](i) - limitTau * sliding(i-1) / slidingNorm );
                    }

                    R2TensorT<2> dUdgT;
                    dUdgT.dyadic_aa( sliding );
                    dUdgT(0,0) = (slidingNorm*slidingNorm - dUdgT(0,0)) * limitTau / std::pow(slidingNorm, 3);
                    dUdgT(0,1) *= - limitTau / std::pow(slidingNorm, 3);
                    dUdgT(1,0) *= - limitTau / std::pow(slidingNorm, 3);
                    dUdgT(1,1) = (slidingNorm*slidingNorm - dUdgT(1,1)) * limitTau / std::pow(slidingNorm, 3);

                    for( localIndex kf=0 ; kf<2 ; ++kf )
                    {
                      for( localIndex a=0 ; a<numNodesPerFace ; ++a )
                      {
                        for( int i=0 ; i<3 ; ++i )
                        {
                          R1TensorT<2> localRowB( rotationMatrix[kfe]( 1, i ), rotationMatrix[kfe]( 2, i ) );
                          R1TensorT<2> localRowE;
                          localRowE.AijBj( dUdgT, localRowB );

                          dRdU( 1, kf*3*numNodesPerFace + 3*a+i ) = - nodalArea * localRowE( 0 ) * pow( -1, kf );
                          dRdU( 2, kf*3*numNodesPerFace + 3*a+i ) = - nodalArea * localRowE( 1 ) * pow( -1, kf );
                        }
                      }
                    }
                    for( int i=1 ; i<3 ; ++i )
                    {
                      dRdT( i, 0 ) = Ja * std::tan( m_frictionAngle ) * sliding(i-1) / slidingNorm;
                      dRdT( i, i ) = Ja;
                    }
                    std::cout << "dRdT" << dRdT << std::endl;
                  }
                  else
                  {
                    R1TensorT<2> vaux( traction[kfe](1), traction[kfe](2) );
                    real64 vauxNorm = sqrt( vaux(0)*vaux(0) + vaux(1)*vaux(1) );
                    if( vauxNorm > 0.0 )
                    {
                      for( int i=1 ; i<3 ; ++i )
                      {
                        elemRHS[i] = + Ja * ( traction[kfe](i) - limitTau * vaux(i-1) / vauxNorm );
                      }
                      for( int i=1 ; i<3 ; ++i )
                      {
                        dRdT( i, i ) = Ja;
                      }
                    }
                    else
                    {
                      for( int i=1 ; i<3 ; ++i )
                      {
                        elemRHS[i] = 0.0;
                      }
                      for( int i=1 ; i<3 ; ++i )
                      {
                        dRdT( i, i ) = Ja;
                      }
                    }
                  }
                  for( int i=0 ; i<3 ; ++i )
                  {
                    std::cout << i << " " << elemRHS[i] << std::endl;
                  }
                  break;
                }
                case FractureState::OPEN:
                {
                  std::cout << kfe << " OPEN" << std::endl;
                  for( int i=0 ; i<3 ; ++i )
                  {
                    elemRHS[i] = + Ja * traction[kfe][i];
                  }

                  for( int i=0 ; i<3 ; ++i )
                  {
                    dRdT( i, i ) = Ja;
                  }
                  break;
                }
              }

              rhs1->add( elemDOF,
                         elemRHS,
                         3 );

              if( fractureState[kfe] != FractureState::OPEN )
              {
                matrix10->add( elemDOF,
                               nodeDOF,
                               dRdU.data(),
                               3,
                               2*3*numNodesPerFace );
              }

              if( fractureState[kfe] != FractureState::STICK )
              {
                matrix11->add( elemDOF,
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

  matrix10->close();
  matrix11->close();
  rhs1->close();
}

void
LagrangianContactSolver::
  ApplySystemSolution( DofManager const & GEOSX_UNUSED_PARAM( dofManager ),
                       ParallelVector const & solution,
                       real64 const scalingFactor,
                       DomainPartition * const domain )
{
  GEOSX_MARK_FUNCTION;

  std::cout << "In ApplySystemSolution!" << std::endl;

  globalIndex numDispDofs = m_solidSolver->getSystemRhs().globalSize();

  ParallelVector & sol0 = m_solidSolver->getSystemSolution();
  sol0.open();
  sol0.zero();
  for( globalIndex i=sol0.ilower() ; i<sol0.iupper() ; ++i )
  {
    if( solution.getLocalRowID( i ) >= 0 )
    {
      real64 const value = solution.get( i );
      sol0.add( i, value );
    }
  }
  sol0.close();

  m_solidSolver->ApplySystemSolution( m_solidSolver->getDofManager(),
                                      m_solidSolver->getSystemSolution(),
                                      scalingFactor,
                                      domain );

  ParallelVector sol1;
  sol1.createWithLocalSize( m_rhs1.localSize(),
                            MPI_COMM_GEOSX );
  sol1.open();
  sol1.zero();
  for( globalIndex i=sol1.ilower() ; i<sol1.iupper() ; ++i )
  {
    if( solution.getLocalRowID( i ) >= 0 )
    {
      real64 const value = solution.get( numDispDofs+i );
      sol1.add( i, value );
    }
  }
  sol1.close();

  getDofManager().addVectorToField( sol1, viewKeyStruct::tractionString, viewKeyStruct::deltaTractionString, -scalingFactor );
  getDofManager().addVectorToField( sol1, viewKeyStruct::tractionString, viewKeyStruct::tractionString, -scalingFactor );

  std::map< string, string_array > fieldNames;
  fieldNames["elem"].push_back( viewKeyStruct::tractionString );
  fieldNames["elem"].push_back( viewKeyStruct::deltaTractionString );
  fieldNames["elem"].push_back( viewKeyStruct::fractureStateString );
  fieldNames["elem"].push_back( FaceElementSubRegion::viewKeyStruct::elementLocalJumpString );

  CommunicationTools::SynchronizeFields( fieldNames,
                                         domain->getMeshBody( 0 )->getMeshLevel( 0 ),
                                         domain->getReference< array1d< NeighborCommunicator > >( domain->viewKeys.neighbors ) );

  this->UpdateDeformationForCoupling( domain );
}

void LagrangianContactSolver::InitializeFractureState( MeshLevel * const mesh )
{
  GEOSX_MARK_FUNCTION;
  ElementRegionManager * const elemManager = mesh->getElemManager();

  string const tractionKey = viewKeyStruct::tractionString;

  elemManager->forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion * const subRegion )->void
  {
    if( subRegion->hasWrapper( tractionKey ) )
    {
      arrayView1d< FractureState > const & fractureState = subRegion->getReference< array1d< FractureState > >( viewKeyStruct::fractureStateString );
      fractureState = FractureState::STICK;
    }
  });
}

bool LagrangianContactSolver::UpdateFractureState( DomainPartition * const domain )
{
  GEOSX_MARK_FUNCTION;
  std::cout << "In UpdateFractureState!" << std::endl;

  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup< MeshBody >( 0 )->getMeshLevel( 0 );
  ElementRegionManager * const elemManager = mesh->getElemManager();

  bool checkActiveSet = true;

  elemManager->forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion * const subRegion )->void
  {
    if( subRegion->hasWrapper( m_tractionKey ) )
    {
      arrayView2d< real64 const > const & traction = subRegion->getReference< array2d< real64 > >( viewKeyStruct::tractionString );
      arrayView2d< real64 const > const & localJump = subRegion->getReference< array2d< real64 > >( FaceElementSubRegion::viewKeyStruct::elementLocalJumpString );
      arrayView1d< FractureState > const & fractureState = subRegion->getReference< array1d< FractureState > >( viewKeyStruct::fractureStateString );
      array1d< bool > localCheck( subRegion->size() );

      forall_in_range< serialPolicy >( 0,
                                       subRegion->size(),
                                       [&]( localIndex const kfe )
        {
          {
            std::cout << traction[kfe] << std::endl;

            FractureState const originalFractureState = fractureState[kfe];
            if( originalFractureState == FractureState::OPEN )
            {
              if( localJump[kfe](0) > - m_normalDisplacementTolerance )
              {
                fractureState[kfe] = FractureState::OPEN;
              }
              else
              {
                fractureState[kfe] = FractureState::STICK;
              }
            }
            else if( traction[kfe](0) > m_normalTractionTolerance )
            {
              std::cout << "need to be open!!!\n";
              fractureState[kfe] = FractureState::OPEN;
            }
            else
            {
              real64 currentTau = sqrt( traction[kfe](1)*traction[kfe](1) + traction[kfe](2)*traction[kfe](2) );
              real64 limitTau = m_cohesion - traction[kfe](0) * std::tan( m_frictionAngle );
              if( originalFractureState != FractureState::SLIP && currentTau >= limitTau )
              {
                currentTau *= (1.0 - m_alpha);
              }
              else if( originalFractureState == FractureState::SLIP && currentTau <= limitTau )
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
            std::cout << kfe << " " << FractureStateToString( originalFractureState ) << std::endl;
            std::cout << kfe << " " << FractureStateToString( fractureState[kfe] ) << std::endl;
            localCheck[kfe] = ( originalFractureState == fractureState[kfe] );
            std::cout << kfe << " *** " << localCheck[kfe] << std::endl;
          }
        });

      std::cout << localCheck << std::endl;
      for( localIndex kfe = 0; kfe < subRegion->size(); ++kfe )
      {
        std::cout << kfe << " " << checkActiveSet << " " << localCheck[kfe] << std::endl;
        checkActiveSet &= localCheck[kfe];
      }
    }
  });

  std::cout << "output of UpdateFractureState -> " << checkActiveSet << std::endl;
  return checkActiveSet;
}

void LagrangianContactSolver::SolveSystem( DofManager const & dofManager,
                                           ParallelMatrix & matrix,
                                           ParallelVector & rhs,
                                           ParallelVector & solution )
{
  GEOSX_MARK_FUNCTION;

  std::cout << "In SolveSystem!" << std::endl;
  globalIndex numDispDofs = m_solidSolver->getSystemRhs().globalSize();
  globalIndex numTracDofs = m_matrix11.globalRows();
  GEOSX_LOG_RANK_0( "size = " << numDispDofs << " + " << numTracDofs );

  matrix.open();
  matrix.zero();
  rhs.open();
  rhs.zero();
  for( globalIndex i=m_solidSolver->getSystemMatrix().ilower() ; i<m_solidSolver->getSystemMatrix().iupper() ; ++i )
  {
    globalIndex_array colIndices;
    real64_array values;
    m_solidSolver->getSystemMatrix().getRowCopy( i, colIndices, values );
    // Attention: change sign to K (m_solidSolver provide K < 0, I'd prefer K > 0)
    for( localIndex j = 0 ; j < values.size() ; ++j )
    {
      values[j] *= -1.0;
    }
    matrix.set( i, colIndices, values );

    real64 value = m_solidSolver->getSystemRhs().get( i );
    // Attention: change sign to rhs (m_solidSolver provide K < 0, I'd prefer K > 0)
    rhs.add( i, -value );

    m_matrix01.getRowCopy( i, colIndices, values );
    for( localIndex j = 0 ; j < colIndices.size() ; ++j )
    {
      colIndices[j] += numDispDofs;
    }
    matrix.set( i, colIndices, values );

    value = m_rhs0.get( i );
    rhs.add( i, value );
  }
  for( globalIndex i=m_matrix11.ilower() ; i<m_matrix11.iupper() ; ++i )
  {
    globalIndex_array colIndices;
    real64_array values;
    m_matrix10.getRowCopy( i, colIndices, values );
    matrix.set( numDispDofs+i, colIndices, values );

    m_matrix11.getRowCopy( i, colIndices, values );
    for( localIndex j = 0 ; j < colIndices.size() ; ++j )
    {
      colIndices[j] += numDispDofs;
    }
    matrix.set( numDispDofs+i, colIndices, values );

    real64 const value = m_rhs1.get( i );
    rhs.add( numDispDofs+i, value );
  }
  matrix.close();
  rhs.close();

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
  std::cin.get();
}

real64
LagrangianContactSolver::ScalingForSystemSolution( DomainPartition const * const domain,
                                                   DofManager const & GEOSX_UNUSED_PARAM( dofManager ),
                                                   ParallelVector const & GEOSX_UNUSED_PARAM( solution ) )
{
  // TODO: add something for tractions!!!
  return m_solidSolver->ScalingForSystemSolution( domain,
                                                  m_solidSolver->getDofManager(),
                                                  m_solidSolver->getSystemSolution() );
}

void LagrangianContactSolver::SetNextDt( real64 const & currentDt,
                                         real64 & nextDt )
{

  if( m_numResolves[0] == 0 && m_numResolves[1] == 0 )
  {
    this->SetNextDtBasedOnNewtonIter( currentDt, nextDt );
  }
  else
  {
    std::cout << "call to surfaceGenerator" << std::endl;
    SolverBase * const surfaceGenerator =  this->getParent()->GetGroup< SolverBase >( "SurfaceGen" );
    nextDt = surfaceGenerator->GetTimestepRequest() < 1e99 ? surfaceGenerator->GetTimestepRequest() : currentDt;
  }
  GEOSX_LOG_LEVEL_RANK_0( 3, this->getName() << ": nextDt request is "  << nextDt );
}

void LagrangianContactSolver::initializeNewFaceElements( DomainPartition const & )
{
//  m_flowSolver->
}

REGISTER_CATALOG_ENTRY( SolverBase, LagrangianContactSolver, std::string const &, Group * const )
} /* namespace geosx */
