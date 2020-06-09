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
#include "math/interpolation/Interpolation.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;
using namespace interpolation;

LagrangianContactSolver::LagrangianContactSolver( const std::string & name,
                                                  Group * const parent ):
  SolverBase( name, parent ),
  m_solidSolverName(),
  m_solidSolver( nullptr ),
  m_stabilizationName(),
  m_contactRelationName(),
  m_activeSetMaxIter()
{
  registerWrapper( viewKeyStruct::solidSolverNameString, &m_solidSolverName )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Name of the solid mechanics solver to use in the lagrangian contact solver" );

  registerWrapper( viewKeyStruct::stabilizationNameString, &m_stabilizationName )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Name of the stabilization to use in the lagrangian contact solver" );

  registerWrapper( viewKeyStruct::contactRelationNameString, &m_contactRelationName )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Name of the constitutive law used for fracture elements" );

  registerWrapper( viewKeyStruct::activeSetMaxIterString, &m_activeSetMaxIter )->
    setApplyDefaultValue( 10 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Maximum number of iteration for the active set strategy in the lagrangian contact solver" );
}

void LagrangianContactSolver::RegisterDataOnMesh( dataRepository::Group * const MeshBodies )
{
  for( auto & mesh : MeshBodies->GetSubGroups() )
  {
    MeshLevel * meshLevel = Group::group_cast< MeshBody * >( mesh.second )->getMeshLevel( 0 );

    ElementRegionManager * const elemManager = meshLevel->getElemManager();
    elemManager->forElementRegions< FaceElementRegion >( [&] ( FaceElementRegion & region )
    {
      region.forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion & subRegion )
      {
        subRegion.registerWrapper< array2d< real64 > >( viewKeyStruct::tractionString )->
          setApplyDefaultValue( 0.0 )->
          setPlotLevel( PlotLevel::LEVEL_0 )->
          setRegisteringObjects( this->getName())->
          setDescription( "An array that holds the tractions on the fracture." )->
          reference().resizeDimension< 1 >( 3 );

        subRegion.registerWrapper< array2d< real64 > >( viewKeyStruct::deltaTractionString )->
          setApplyDefaultValue( 0.0 )->
          setPlotLevel( PlotLevel::NOPLOT )->
          setRegisteringObjects( this->getName())->
          setDescription( "An array that holds the traction increments on the fracture." )->
          reference().resizeDimension< 1 >( 3 );

        subRegion.registerWrapper< array1d< integer > >( viewKeyStruct::fractureStateString )->
          setPlotLevel( PlotLevel::LEVEL_0 )->
          setRegisteringObjects( this->getName())->
          setDescription( "An array that holds the fracture state." );
        InitializeFractureState( meshLevel, viewKeyStruct::fractureStateString );

        subRegion.registerWrapper< array1d< integer > >( viewKeyStruct::previousFractureStateString )->
          setPlotLevel( PlotLevel::NOPLOT )->
          setRegisteringObjects( this->getName())->
          setDescription( "An array that holds the fracture state." );
        InitializeFractureState( meshLevel, viewKeyStruct::previousFractureStateString );

        subRegion.registerWrapper< array2d< real64 > >( viewKeyStruct::localJumpString )->
          setApplyDefaultValue( 0.0 )->
          setPlotLevel( PlotLevel::LEVEL_0 )->
          setRegisteringObjects( this->getName())->
          setDescription( "An array that holds the local jump on the fracture at the current time step." )->
          reference().resizeDimension< 1 >( 3 );

        subRegion.registerWrapper< array2d< real64 > >( viewKeyStruct::previousLocalJumpString )->
          setApplyDefaultValue( 0.0 )->
          setPlotLevel( PlotLevel::NOPLOT )->
          setRegisteringObjects( this->getName())->
          setDescription( "An array that holds the local jump on the fracture at the previous time step." )->
          reference().resizeDimension< 1 >( 3 );

        subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::normalTractionToleranceString )->
          setPlotLevel( PlotLevel::NOPLOT )->
          setRegisteringObjects( this->getName())->
          setDescription( "An array that holds the normal traction tolerance." );

        subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::normalDisplacementToleranceString )->
          setPlotLevel( PlotLevel::NOPLOT )->
          setRegisteringObjects( this->getName())->
          setDescription( "An array that holds the normal displacement tolerance." );

        subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::slidingToleranceString )->
          setPlotLevel( PlotLevel::NOPLOT )->
          setRegisteringObjects( this->getName())->
          setDescription( "An array that holds the sliding tolerance." );

        // Needed just because SurfaceGenerator initialize the field "pressure" (NEEDED!!!)
        // It is used in "TwoPointFluxApproximation.cpp", called by "SurfaceGenerator.cpp"
        subRegion.registerWrapper< real64_array >( "pressure" )->
          setPlotLevel( PlotLevel::NOPLOT )->
          setRegisteringObjects( this->getName());
      } );
    } );
  }
}

void LagrangianContactSolver::InitializePreSubGroups( Group * const rootGroup )
{
  SolverBase::InitializePreSubGroups( rootGroup );

  DomainPartition * domain = rootGroup->GetGroup< DomainPartition >( keys::domain );
  ConstitutiveManager const * const cm = domain->getConstitutiveManager();

  ConstitutiveBase const * const contactRelation  = cm->GetConstitutiveRelation< ConstitutiveBase >( m_contactRelationName );
  GEOSX_ERROR_IF( contactRelation == nullptr, "fracture constitutive model " + m_contactRelationName + " not found" );
  m_contactRelationFullIndex = contactRelation->getIndexInParent();
}

void LagrangianContactSolver::ImplicitStepSetup( real64 const & time_n,
                                                 real64 const & dt,
                                                 DomainPartition * const domain,
                                                 DofManager & dofManager,
                                                 ParallelMatrix & matrix,
                                                 ParallelVector & rhs,
                                                 ParallelVector & solution )
{
  ComputeTolerances( domain );

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
  m_solidSolver->ImplicitStepComplete( time_n, dt, domain );

  MeshLevel * const meshLevel = domain->getMeshBody( 0 )->getMeshLevel( 0 );
  ElementRegionManager * const elemManager = meshLevel->getElemManager();

  elemManager->forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion & subRegion )
  {
    if( subRegion.hasWrapper( m_tractionKey ) )
    {
      arrayView1d< integer const > const & ghostRank = subRegion.ghostRank();
      arrayView2d< real64 > const &
      deltaTraction = subRegion.getReference< array2d< real64 > >( viewKeyStruct::deltaTractionString );
      arrayView2d< real64 const > const &
      localJump = subRegion.getReference< array2d< real64 > >( viewKeyStruct::localJumpString );
      arrayView2d< real64 > const &
      previousLocalJump = subRegion.getReference< array2d< real64 > >( viewKeyStruct::previousLocalJumpString );
      arrayView1d< integer const > const &
      fractureState = subRegion.getReference< array1d< integer > >( viewKeyStruct::fractureStateString );
      arrayView1d< integer > const &
      previousFractureState = subRegion.getReference< array1d< integer > >( viewKeyStruct::previousFractureStateString );

      forAll< serialPolicy >( subRegion.size(), [=]( localIndex const kfe )
      {
        if( ghostRank[kfe] < 0 )
        {
          for( localIndex i = 0; i < 3; ++i )
          {
            deltaTraction[kfe][i] = 0.0;
            previousLocalJump[kfe][i] = localJump[kfe][i];
          }
          previousFractureState[kfe] = fractureState[kfe];
        }
      } );
    }
  } );

  // Need a synchronization of deltaTraction as will be used in AssembleStabilization
  std::map< string, string_array > fieldNames;
  fieldNames["elems"].push_back( viewKeyStruct::deltaTractionString );
  CommunicationTools::SynchronizeFields( fieldNames,
                                         domain->getMeshBody( 0 )->getMeshLevel( 0 ),
                                         domain->getNeighbors() );

  GEOSX_LOG_LEVEL_RANK_0( 1, " ***** ImplicitStepComplete *****" );
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

void LagrangianContactSolver::ComputeTolerances( DomainPartition * const domain ) const
{
  GEOSX_MARK_FUNCTION;

  MeshLevel * const mesh = domain->getMeshBody( 0 )->getMeshLevel( 0 );

  FaceManager const * const faceManager = mesh->getFaceManager();
  NodeManager const * const nodeManager = mesh->getNodeManager();
  ElementRegionManager * const elemManager = mesh->getElemManager();

  // Get the "face to element" map (valid for the entire mesh)
  FaceManager::ElemMapType const & faceToElem = faceManager->toElementRelation();

  // Get the volume for all elements
  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > const elemVolume =
    elemManager->ConstructViewAccessor< real64_array, arrayView1d< real64 const > >( ElementSubRegionBase::viewKeyStruct::elementVolumeString );

  // Get the coordinates for all nodes
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition = nodeManager->referencePosition();

  // Bulk modulus accessor
  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > const bulkModulus =
    elemManager->ConstructMaterialViewAccessor< array1d< real64 >, arrayView1d< real64 const > >( LinearElasticIsotropic::viewKeyStruct::bulkModulusString,
                                                                                                  m_solidSolver->targetRegionNames(),
                                                                                                  m_solidSolver->solidMaterialNames() );
  // Shear modulus accessor
  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > const shearModulus =
    elemManager->ConstructMaterialViewAccessor< array1d< real64 >, arrayView1d< real64 const > >( LinearElasticIsotropic::viewKeyStruct::shearModulusString,
                                                                                                  m_solidSolver->targetRegionNames(),
                                                                                                  m_solidSolver->solidMaterialNames() );
  elemManager->forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion & subRegion )
  {
    if( subRegion.hasWrapper( m_tractionKey ) )
    {
      arrayView1d< integer const > const & ghostRank = subRegion.ghostRank();
      arrayView1d< real64 const > const & faceArea = subRegion.getElementArea();
      arrayView3d< real64 const > const & faceRotationMatrix = subRegion.getElementRotationMatrix();
      arrayView2d< localIndex const > const & elemsToFaces = subRegion.faceList();

      arrayView1d< real64 > const &
      normalTractionTolerance = subRegion.getReference< array1d< real64 > >( viewKeyStruct::normalTractionToleranceString );
      arrayView1d< real64 > const &
      normalDisplacementTolerance = subRegion.getReference< array1d< real64 > >( viewKeyStruct::normalDisplacementToleranceString );
      arrayView1d< real64 > const &
      slidingTolerance = subRegion.getReference< array1d< real64 > >( viewKeyStruct::slidingToleranceString );

      forAll< serialPolicy >( subRegion.size(), [=]( localIndex const kfe )
      {
        if( ghostRank[kfe] < 0 )
        {
          real64 const area = faceArea[kfe];
          R1Tensor stiffApprox[ 2 ];
          real64 averageYoungModulus = 0.0;
          real64 averageConstrainedModulus = 0.0;
          real64 averageBoxSize0 = 0.0;
          for( localIndex i = 0; i < 2; ++i )
          {
            localIndex faceIndex = elemsToFaces[kfe][i];
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

            real64 boxSize[ 3 ];
            for( localIndex j = 0; j < 3; ++j )
            {
              boxSize[j] = maxSize[j] - minSize[j];
            }

            // Get linear elastic isotropic constitutive parameters for the element
            real64 const K = bulkModulus[er][esr][ei];
            real64 const G = shearModulus[er][esr][ei];
            real64 const E = 9.0 * K * G / ( 3.0 * K + G );
            real64 const nu = ( 3.0 * K - 2.0 * G ) / ( 2.0 * ( 3.0 * K + G ) );
            real64 const M = K + 4.0 / 3.0 * G;

            for( localIndex j = 0; j < 3; ++j )
            {
              stiffApprox[i]( j ) = E / ( ( 1.0 + nu )*( 1.0 - 2.0*nu ) ) * 8.0 / 9.0 * ( 2.0 - 3.0 * nu ) * volume / ( boxSize[j]*boxSize[j] );
            }

            averageYoungModulus += 0.5*E;
            averageConstrainedModulus += 0.5*M;
            averageBoxSize0 += 0.5*boxSize[0];
          }

          real64 invStiffApprox[ 3 ][ 3 ] = { { 0 } };
          for( localIndex j = 0; j < 3; ++j )
          {
            invStiffApprox[ j ][ j ] = ( stiffApprox[0]( j ) + stiffApprox[1]( j ) ) / ( stiffApprox[0]( j ) * stiffApprox[1]( j ) );
          }

          // Compute R^T * (invK) * R
          real64 temp[ 3 ][ 3 ];
          LvArray::tensorOps::AkiBkj< 3, 3, 3 >( temp, faceRotationMatrix[ kfe ], invStiffApprox );
          real64 rotatedInvStiffApprox[ 3 ][ 3 ];
          LvArray::tensorOps::AikBkj< 3, 3, 3 >( rotatedInvStiffApprox, temp, faceRotationMatrix[ kfe ] );
          LvArray::tensorOps::scale< 3, 3 >( rotatedInvStiffApprox, area );

          normalDisplacementTolerance[kfe] = rotatedInvStiffApprox[ 0 ][ 0 ] * averageYoungModulus / 2.e+7;
          slidingTolerance[kfe] = std::sqrt( rotatedInvStiffApprox[ 1 ][ 1 ] * rotatedInvStiffApprox[ 1 ][ 1 ] +
                                             rotatedInvStiffApprox[ 2 ][ 2 ] * rotatedInvStiffApprox[ 2 ][ 2 ] ) * averageYoungModulus / 2.e+7;

          normalTractionTolerance[kfe] = 1.0 / 2.0 * averageConstrainedModulus / averageBoxSize0 * normalDisplacementTolerance[kfe];
        }
      } );
    }
  } );
}

void LagrangianContactSolver::ResetStateToBeginningOfStep( DomainPartition * const domain )
{
  m_solidSolver->ResetStateToBeginningOfStep( domain );

  MeshLevel * const meshLevel = domain->getMeshBody( 0 )->getMeshLevel( 0 );
  ElementRegionManager * const elemManager = meshLevel->getElemManager();

  elemManager->forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion & subRegion )
  {
    if( subRegion.hasWrapper( m_tractionKey ) )
    {
      arrayView1d< integer const > const & ghostRank = subRegion.ghostRank();
      arrayView2d< real64 > const & traction = subRegion.getReference< array2d< real64 > >( viewKeyStruct::tractionString );
      arrayView2d< real64 > const & deltaTraction = subRegion.getReference< array2d< real64 > >( viewKeyStruct::deltaTractionString );
      arrayView2d< real64 > const & localJump = subRegion.getReference< array2d< real64 > >( viewKeyStruct::localJumpString );
      arrayView2d< real64 const > const & previousLocalJump = subRegion.getReference< array2d< real64 > >( viewKeyStruct::previousLocalJumpString );
      arrayView1d< integer > const & fractureState = subRegion.getReference< array1d< integer > >( viewKeyStruct::fractureStateString );
      arrayView1d< integer const > const & previousFractureState =
        subRegion.getReference< array1d< integer > >( viewKeyStruct::previousFractureStateString );

      forAll< serialPolicy >( subRegion.size(), [=]( localIndex const kfe )
      {
        if( ghostRank[kfe] < 0 )
        {
          for( localIndex i = 0; i < 3; ++i )
          {
            traction[kfe][i] -= deltaTraction[kfe][i];
            deltaTraction[kfe][i] = 0.0;

            localJump[kfe][i] = previousLocalJump[kfe][i];
          }
          fractureState[kfe] = previousFractureState[kfe];
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

  FaceManager::NodeMapType const & faceToNodeMap = faceManager->nodeList();

  arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const & u = nodeManager->totalDisplacement();

  elemManager->forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion & subRegion )
  {
    if( subRegion.hasWrapper( m_tractionKey ) )
    {
      arrayView1d< integer const > const & ghostRank = subRegion.ghostRank();
      arrayView3d< real64 const > const & rotationMatrix = subRegion.getElementRotationMatrix();
      arrayView2d< localIndex const > const & elemsToFaces = subRegion.faceList();
      arrayView2d< real64 > const & localJump = subRegion.getReference< array2d< real64 > >( viewKeyStruct::localJumpString );

      forAll< serialPolicy >( subRegion.size(), [=]( localIndex const kfe )
      {
        // Contact constraints
        if( ghostRank[kfe] < 0 )
        {
          localIndex const numNodesPerFace = faceToNodeMap.sizeOfArray( elemsToFaces[kfe][0] );
          real64 globalJumpTemp[ 3 ] = { 0 };
          for( localIndex a=0; a<numNodesPerFace; ++a )
          {
            for( localIndex i=0; i<3; ++i )
            {
              globalJumpTemp[ i ] +=
                ( -u[faceToNodeMap( elemsToFaces[kfe][0], a )][i]
                  + u[faceToNodeMap( elemsToFaces[kfe][1], a )][i] ) / numNodesPerFace;
            }
          }

          real64 localJumpTemp[ 3 ];
          LvArray::tensorOps::AjiBj< 3, 3 >( localJumpTemp, rotationMatrix[ kfe ], globalJumpTemp );
          LvArray::tensorOps::copy< 3 >( localJump[ kfe ], localJumpTemp );

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

  bool useElasticStep = !IsFractureAllInStickCondition( domain );

  // outer loop attempts to apply full timestep, and managed the cutting of the timestep if
  // required.
  for( dtAttempt = 0; dtAttempt < maxNumberDtCuts; ++dtAttempt )
  {
    // reset the solver state, since we are restarting the time step
    if( dtAttempt > 0 )
    {
      ResetStateToBeginningOfStep( domain );
      globalIndex numStick, numSlip, numOpen;
      ComputeFractureStateStatistics( domain, numStick, numSlip, numOpen, true );
    }

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
          krylovParams.relTolerance = EisenstatWalker( residualNorm, lastResidual, krylovParams.weakestTol );
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
      bool const isPreviousFractureStateValid = UpdateFractureState( domain );
      GEOSX_LOG_LEVEL_RANK_0( 1, "active set flag: " << std::boolalpha << isPreviousFractureStateValid );

      if( getLogLevel() >= 1 )
      {
        globalIndex numStick, numSlip, numOpen;
        ComputeFractureStateStatistics( domain, numStick, numSlip, numOpen, true );
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
        ResetStateToBeginningOfStep( domain );
        SetFractureStateForElasticStep( domain );
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
    GEOSX_ERROR( "Active set did not reached a solution. Terminating..." );
  }
  else
  {
    GEOSX_LOG_RANK_0( "Number of active set iterations: " << m_activeSetIter );
  }

  // return the achieved timestep
  return stepDt;
}

bool LagrangianContactSolver::LineSearch( real64 const & time_n,
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
      localScaleFactor = ParabolicInterpolationThreePoints( lamc, lamm, ff0, ffT, ffm );
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

void LagrangianContactSolver::SetupDofs( DomainPartition const * const domain,
                                         DofManager & dofManager ) const
{
  GEOSX_MARK_FUNCTION;
  m_solidSolver->SetupDofs( domain, dofManager );

  // restrict coupling to fracture regions only
  ElementRegionManager const * const elemManager = domain->getMeshBody( 0 )->getMeshLevel( 0 )->getElemManager();
  string_array fractureRegions;
  elemManager->forElementRegions< FaceElementRegion >( [&]( FaceElementRegion const & elementRegion )
  {
    fractureRegions.push_back( elementRegion.getName() );
  } );

  dofManager.addField( viewKeyStruct::tractionString,
                       DofManager::Location::Elem,
                       3,
                       fractureRegions );
  dofManager.addCoupling( viewKeyStruct::tractionString,
                          viewKeyStruct::tractionString,
                          DofManager::Connector::Face,
                          fractureRegions );
  dofManager.addCoupling( keys::TotalDisplacement,
                          viewKeyStruct::tractionString,
                          DofManager::Connector::Elem,
                          fractureRegions );
}

void LagrangianContactSolver::SetupSystem( DomainPartition * const domain,
                                           DofManager & dofManager,
                                           ParallelMatrix & matrix,
                                           ParallelVector & rhs,
                                           ParallelVector & solution,
                                           bool const GEOSX_UNUSED_PARAM( setSparsity ) )
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
                              5*(3*27+3*12),
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

  SynchronizeFractureState( domain );

  m_solidSolver->AssembleSystem( time,
                                 dt,
                                 domain,
                                 dofManager,
                                 matrix,
                                 rhs );

  AssembleForceResidualDerivativeWrtTraction( domain, dofManager, &matrix, &rhs );
  AssembleTractionResidualDerivativeWrtDisplacementAndTraction( domain, dofManager, &matrix, &rhs );
  AssembleStabilization( domain, dofManager, &matrix, &rhs );
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

  localIndex numDispDofs = dofManager.numLocalDofs( keys::TotalDisplacement );
  localIndex numTracDofs = dofManager.numLocalDofs( viewKeyStruct::tractionString );
  real64 const * localResidual = rhs.extractLocalVector();
  real64 localResidualNorm[3] = {0.0, 0.0, 0.0};
  for( localIndex i=0; i<numDispDofs; ++i )
  {
    localResidualNorm[0] += localResidual[i] * localResidual[i];
  }
  for( localIndex i=numDispDofs; i<numDispDofs+numTracDofs; ++i )
  {
    localResidualNorm[1] += localResidual[i] * localResidual[i];
  }
  localResidualNorm[2] = localResidualNorm[0] + localResidualNorm[1];

  real64 globalResidualNorm[3] = {0.0, 0.0, 0.0};
  int const rank = MpiWrapper::Comm_rank( MPI_COMM_GEOSX );
  int const size = MpiWrapper::Comm_size( MPI_COMM_GEOSX );
  real64_array globalValues( 3*size );
  globalValues = 0;

  // Everything is done on rank 0
  MpiWrapper::gather( localResidualNorm,
                      3,
                      globalValues.data(),
                      3,
                      0,
                      MPI_COMM_GEOSX );

  if( rank==0 )
  {
    for( int r=0; r<size; ++r )
    {
      // sum across all ranks
      globalResidualNorm[0] += globalValues[3*r];
      globalResidualNorm[1] += globalValues[3*r+1];
      globalResidualNorm[2] += globalValues[3*r+2];
    }
    globalResidualNorm[0] = sqrt( globalResidualNorm[0] );
    globalResidualNorm[1] = sqrt( globalResidualNorm[1] );
    globalResidualNorm[2] = sqrt( globalResidualNorm[2] );
  }

  MpiWrapper::bcast( globalResidualNorm, 3, 0, MPI_COMM_GEOSX );

  if( m_nonlinearSolverParameters.m_numNewtonIterations == 0 )
  {
    m_initialResidual[0] = globalResidualNorm[0];
    m_initialResidual[1] = globalResidualNorm[1];
    m_initialResidual[2] = globalResidualNorm[2];
    globalResidualNorm[0] = 1.0;
    globalResidualNorm[1] = 1.0;
    globalResidualNorm[2] = 1.0;
  }
  else
  {
    globalResidualNorm[0] /= (m_initialResidual[0]+1.0);
    globalResidualNorm[1] /= (m_initialResidual[1]+1.0);
    // Add 0 just to match Matlab code results
    globalResidualNorm[2] /= (m_initialResidual[2]+0.0);
  }

//  char output[69] = {0};
//  sprintf( output,
//           "( Rdisplacement, Rtraction ) = ( %15.6e, %15.6e );",
//           globalResidualNorm[0],
//           globalResidualNorm[1] );
//  GEOSX_LOG_RANK_0( output );

//  return ( globalResidualNorm[0] + globalResidualNorm[1] );

  char output[94] = {0};
  sprintf( output,
           "( Rdisplacement, Rtraction, Rtotal ) = ( %15.6e, %15.6e, %15.6e );",
           globalResidualNorm[0],
           globalResidualNorm[1],
           globalResidualNorm[2] );
  GEOSX_LOG_RANK_0( output );

  return globalResidualNorm[2];
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

  FaceManager::NodeMapType const & faceToNodeMap = faceManager->nodeList();

  arrayView1d< R1Tensor > const &
  fext = nodeManager->getReference< array1d< R1Tensor > >( SolidMechanicsLagrangianFEM::viewKeyStruct::forceExternal );
  fext = {0, 0, 0};

  string const tracDofKey = dofManager.getKey( viewKeyStruct::tractionString );
  string const dispDofKey = dofManager.getKey( keys::TotalDisplacement );

  arrayView1d< globalIndex > const &
  dispDofNumber = nodeManager->getReference< globalIndex_array >( dispDofKey );

  matrix->open();
  rhs->open();

  elemManager->forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion const & subRegion )
  {
    if( subRegion.hasWrapper( m_tractionKey ) )
    {
      arrayView1d< globalIndex const > const &
      tracDofNumber = subRegion.getReference< globalIndex_array >( tracDofKey );
      arrayView2d< real64 const > const & traction = subRegion.getReference< array2d< real64 > >( viewKeyStruct::tractionString );
      arrayView1d< integer const > const & ghostRank = subRegion.ghostRank();
      arrayView1d< real64 const > const & area = subRegion.getElementArea();
      arrayView3d< real64 const > const & rotationMatrix = subRegion.getElementRotationMatrix();
      arrayView2d< localIndex const > const & elemsToFaces = subRegion.faceList();

      forAll< serialPolicy >( subRegion.size(), [=]( localIndex const kfe )
      {
        if( ghostRank[kfe] < 0 )
        {
          localIndex const kf0 = elemsToFaces[kfe][0];
          localIndex const numNodesPerFace = faceToNodeMap.sizeOfArray( kf0 );

          globalIndex rowDOF[12];
          real64 nodeRHS[12];
          stackArray2d< real64, 3*4*3 > dRdT( 3*numNodesPerFace, 3 );
          dRdT = 0.0;
          globalIndex colDOF[3];
          for( localIndex i=0; i<3; ++i )
          {
            colDOF[i] = tracDofNumber[kfe] + LvArray::integerConversion< globalIndex >( i );
          }

          real64 const nodalArea = area[kfe] / static_cast< real64 >( numNodesPerFace );
          real64 const localNodalForce[ 3 ] = { traction( kfe, 0 ) * nodalArea, traction( kfe, 1 ) * nodalArea, traction( kfe, 2 ) * nodalArea };

          real64 globalNodalForce[ 3 ];
          LvArray::tensorOps::AijBj< 3, 3 >( globalNodalForce, rotationMatrix[ kfe ], localNodalForce );

          for( localIndex kf=0; kf<2; ++kf )
          {
            localIndex const faceIndex = elemsToFaces[kfe][kf];

            for( localIndex a=0; a<numNodesPerFace; ++a )
            {
              for( localIndex i=0; i<3; ++i )
              {
                rowDOF[3*a+i] = dispDofNumber[faceToNodeMap( faceIndex, a )] + LvArray::integerConversion< globalIndex >( i );
                // Opposite sign w.r.t. theory because of minus sign in stiffness matrix definition (K < 0)
                nodeRHS[3*a+i] = +globalNodalForce[i] * pow( -1, kf );
                fext[faceToNodeMap( faceIndex, a )][i] += +globalNodalForce[i] * pow( -1, kf );

                // Opposite sign w.r.t. theory because of minus sign in stiffness matrix definition (K < 0)
                dRdT( 3*a+i, 0 ) = +nodalArea * rotationMatrix( kfe, i, 0 ) * pow( -1, kf );
                dRdT( 3*a+i, 1 ) = +nodalArea * rotationMatrix( kfe, i, 1 ) * pow( -1, kf );
                dRdT( 3*a+i, 2 ) = +nodalArea * rotationMatrix( kfe, i, 2 ) * pow( -1, kf );
              }
            }

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

  matrix->close();
  rhs->close();
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

  ConstitutiveManager const * const
  constitutiveManager = domain->GetGroup< ConstitutiveManager >( keys::ConstitutiveManager );
  ContactRelationBase const * const
  contactRelation = constitutiveManager->GetGroup< ContactRelationBase const >( m_contactRelationName );

  FaceManager::NodeMapType const & faceToNodeMap = faceManager->nodeList();

  string const tracDofKey = dofManager.getKey( viewKeyStruct::tractionString );
  string const dispDofKey = dofManager.getKey( keys::TotalDisplacement );

  arrayView1d< globalIndex const > const &
  dispDofNumber = nodeManager->getReference< globalIndex_array >( dispDofKey );

  matrix->open();
  rhs->open();

  elemManager->forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion const & subRegion )
  {
    if( subRegion.hasWrapper( m_tractionKey ) )
    {
      arrayView1d< globalIndex const > const &
      tracDofNumber = subRegion.getReference< globalIndex_array >( tracDofKey );
      arrayView1d< integer const > const & ghostRank = subRegion.ghostRank();
      arrayView1d< real64 const > const & area = subRegion.getElementArea();
      arrayView3d< real64 const > const & rotationMatrix = subRegion.getElementRotationMatrix();
      arrayView2d< localIndex const > const & elemsToFaces = subRegion.faceList();
      arrayView2d< real64 const > const &
      traction = subRegion.getReference< array2d< real64 > >( viewKeyStruct::tractionString );
      arrayView1d< integer const > const &
      fractureState = subRegion.getReference< array1d< integer > >( viewKeyStruct::fractureStateString );
      arrayView2d< real64 const > const &
      localJump = subRegion.getReference< array2d< real64 > >( viewKeyStruct::localJumpString );
      arrayView2d< real64 const > const &
      previousLocalJump = subRegion.getReference< array2d< real64 > >( viewKeyStruct::previousLocalJumpString );
      arrayView1d< real64 const > const &
      slidingTolerance = subRegion.getReference< array1d< real64 > >( viewKeyStruct::slidingToleranceString );

      forAll< serialPolicy >( subRegion.size(), [=]( localIndex const kfe )
      {
        if( ghostRank[kfe] < 0 )
        {
          localIndex const numNodesPerFace = faceToNodeMap.sizeOfArray( elemsToFaces[kfe][0] );
          globalIndex nodeDOF[24];
          globalIndex elemDOF[3];
          for( localIndex i=0; i<3; ++i )
          {
            elemDOF[i] = tracDofNumber[kfe] + LvArray::integerConversion< globalIndex >( i );
          }

          real64 elemRHS[3];
          real64 const Ja = area[kfe];
          real64 const nodalArea = Ja / static_cast< real64 >( numNodesPerFace );

          stackArray2d< real64, 2*3*4*3 > dRdU( 3, 2*3*numNodesPerFace );
          stackArray2d< real64, 3*3 > dRdT( 3, 3 );
          dRdU = 0.0;
          dRdT = 0.0;

          switch( fractureState[kfe] )
          {
            case FractureState::STICK:
              {
                for( localIndex i=0; i<3; ++i )
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

                for( localIndex kf=0; kf<2; ++kf )
                {
                  for( localIndex a=0; a<numNodesPerFace; ++a )
                  {
                    for( localIndex i=0; i<3; ++i )
                    {
                      nodeDOF[ kf*3*numNodesPerFace + 3*a+i ] = dispDofNumber[faceToNodeMap( elemsToFaces[kfe][kf], a )] +
                                                                LvArray::integerConversion< globalIndex >( i );

                      dRdU( 0, kf*3*numNodesPerFace + 3*a+i ) = -nodalArea * rotationMatrix( kfe, i, 0 ) * pow( -1, kf );
                      dRdU( 1, kf*3*numNodesPerFace + 3*a+i ) = -nodalArea * rotationMatrix( kfe, i, 1 ) * pow( -1, kf );
                      dRdU( 2, kf*3*numNodesPerFace + 3*a+i ) = -nodalArea * rotationMatrix( kfe, i, 2 ) * pow( -1, kf );
                    }
                  }
                }
                break;
              }
            case FractureState::SLIP:
            case FractureState::NEW_SLIP:
              {
                elemRHS[0] = +Ja * localJump[kfe][0];

                for( localIndex kf=0; kf<2; ++kf )
                {
                  for( localIndex a=0; a<numNodesPerFace; ++a )
                  {
                    for( localIndex i=0; i<3; ++i )
                    {
                      nodeDOF[ kf*3*numNodesPerFace + 3*a+i ] = dispDofNumber[faceToNodeMap( elemsToFaces[kfe][kf], a )] +
                                                                LvArray::integerConversion< globalIndex >( i );
                      dRdU( 0, kf*3*numNodesPerFace + 3*a+i ) = -nodalArea * rotationMatrix( kfe, i, 0 ) * pow( -1, kf );
                    }
                  }
                }

                real64 const limitTau = contactRelation->limitTangentialTractionNorm( traction[kfe][0] );
                R1TensorT< 2 > sliding( localJump[kfe][1] - previousLocalJump[kfe][1], localJump[kfe][2] - previousLocalJump[kfe][2] );
                real64 slidingNorm = sqrt( sliding( 0 )*sliding( 0 ) + sliding( 1 )*sliding( 1 ) );

                GEOSX_LOG_LEVEL_BY_RANK( 3, "element: " << kfe << " sliding: " << sliding );

                if( !( ( m_nonlinearSolverParameters.m_numNewtonIterations == 0 ) && ( fractureState[kfe] == FractureState::NEW_SLIP ) )
                    && slidingNorm > slidingTolerance[kfe] )
                {
                  for( localIndex i=1; i<3; ++i )
                  {
                    elemRHS[i] = +Ja * ( traction[kfe][i] - limitTau * sliding( i-1 ) / slidingNorm );
                  }

                  // A symmetric 2x2 matrix.
                  real64 dUdgT[ 3 ];
                  dUdgT[ 0 ] = (slidingNorm * slidingNorm - sliding[ 0 ] * sliding[ 0 ]) * limitTau / std::pow( slidingNorm, 3 );
                  dUdgT[ 1 ] = (slidingNorm * slidingNorm - sliding[ 1 ] * sliding[ 1 ]) * limitTau / std::pow( slidingNorm, 3 );
                  dUdgT[ 2 ] = -sliding[ 0 ] * sliding[ 1 ] * limitTau / std::pow( slidingNorm, 3 );

                  for( localIndex kf=0; kf<2; ++kf )
                  {
                    for( localIndex a=0; a<numNodesPerFace; ++a )
                    {
                      for( localIndex i=0; i<3; ++i )
                      {
                        real64 const localRowB[ 2 ] = { rotationMatrix( kfe, i, 1 ), rotationMatrix( kfe, i, 2 ) };
                        real64 localRowE[ 2 ];
                        LvArray::tensorOps::symAijBj< 2 >( localRowE, dUdgT, localRowB );

                        dRdU( 1, kf * 3 * numNodesPerFace + 3 * a + i ) = nodalArea * localRowE[ 0 ] * pow( -1, kf );
                        dRdU( 2, kf * 3 * numNodesPerFace + 3 * a + i ) = nodalArea * localRowE[ 1 ] * pow( -1, kf );
                      }
                    }
                  }
                  for( localIndex i=1; i<3; ++i )
                  {
                    dRdT( i, 0 ) = Ja * contactRelation->dLimitTangentialTractionNorm_dNormalTraction( traction[kfe][0] ) * sliding( i-1 ) / slidingNorm;
                    dRdT( i, i ) = Ja;
                  }
                }
                else
                {
                  R1TensorT< 2 > vaux( traction[kfe][1], traction[kfe][2] );
                  real64 vauxNorm = sqrt( vaux( 0 )*vaux( 0 ) + vaux( 1 )*vaux( 1 ) );
                  if( vauxNorm > 0.0 )
                  {
                    for( localIndex i=1; i<3; ++i )
                    {
                      elemRHS[i] = +Ja * ( traction[kfe][i] - limitTau * vaux( i-1 ) / vauxNorm );
                    }
                    for( localIndex i=1; i<3; ++i )
                    {
                      dRdT( i, i ) = Ja;
                    }
                  }
                  else
                  {
                    for( localIndex i=1; i<3; ++i )
                    {
                      elemRHS[i] = 0.0;
                    }
                    for( localIndex i=1; i<3; ++i )
                    {
                      dRdT( i, i ) = Ja;
                    }
                  }
                }
                break;
              }
            case FractureState::OPEN:
              {
                GEOSX_LOG_LEVEL_BY_RANK( 3, "element: " << kfe << " opening: " << localJump[kfe][0] );

                for( localIndex i=0; i<3; ++i )
                {
                  elemRHS[i] = +Ja * traction[kfe][i];
                }

                for( localIndex i=0; i<3; ++i )
                {
                  dRdT( i, i ) = Ja;
                }
                break;
              }
          }

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
      } );
    }
  } );

  matrix->close();
  rhs->close();
}

void LagrangianContactSolver::AssembleStabilization( DomainPartition const * const domain,
                                                     DofManager const & dofManager,
                                                     ParallelMatrix * const matrix,
                                                     ParallelVector * const rhs )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel const * const mesh = domain->getMeshBody( 0 )->getMeshLevel( 0 );

  FaceManager const * const faceManager = mesh->getFaceManager();
  NodeManager const * const nodeManager = mesh->getNodeManager();
  ElementRegionManager const * const elemManager = mesh->getElemManager();

  string const tracDofKey = dofManager.getKey( viewKeyStruct::tractionString );

  // Get the finite volume method used to compute the stabilization
  NumericalMethodsManager const & numericalMethodManager = domain->getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const * const stabilizationMethod = fvManager.getFluxApproximation( m_stabilizationName );

  // Get the "face to element" map (valid for the entire mesh)
  FaceManager::ElemMapType const & faceToElem = faceManager->toElementRelation();

  // Form the SurfaceGenerator, get the fracture name and use it to retrieve the faceMap (from fracture element to face)
  SurfaceGenerator const * const
  surfaceGenerator = this->getParent()->GetGroup< SolverBase >( "SurfaceGen" )->group_cast< SurfaceGenerator const * >();
  FaceElementRegion const * const fractureRegion = elemManager->GetRegion< FaceElementRegion >( surfaceGenerator->getFractureRegionName() );
  FaceElementSubRegion const * const fractureSubRegion = fractureRegion->GetSubRegion< FaceElementSubRegion >( "default" );
  GEOSX_ERROR_IF( !fractureSubRegion->hasWrapper( m_tractionKey ), "The fracture subregion must contain traction field." );
  FaceElementSubRegion::FaceMapType const & faceMap = fractureSubRegion->faceList();
  GEOSX_ERROR_IF( faceMap.size( 1 ) != 2, "A fracture face has to be shared by two cells." );

  // Get the state of fracture elements
  arrayView1d< integer const > const &
  fractureState = fractureSubRegion->getReference< array1d< integer > >( viewKeyStruct::fractureStateString );

  // Get the tractions and stabilization contribution to the local jump
  arrayView2d< real64 const > const &
  traction = fractureSubRegion->getReference< array2d< real64 > >( viewKeyStruct::tractionString );
  arrayView2d< real64 const > const &
  deltaTraction = fractureSubRegion->getReference< array2d< real64 > >( viewKeyStruct::deltaTractionString );

  // Get the volume for all elements
  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > const elemVolume =
    elemManager->ConstructViewAccessor< real64_array, arrayView1d< real64 const > >( ElementSubRegionBase::viewKeyStruct::elementVolumeString );

  // Get the coordinates for all nodes
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition = nodeManager->referencePosition();

  // Get area and rotation matrix for all faces
  arrayView1d< real64 const > const & faceArea = faceManager->faceArea();
  arrayView3d< real64 const > const & faceRotationMatrix = faceManager->faceRotationMatrix();

  // Bulk modulus accessor
  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > const bulkModulus =
    elemManager->ConstructMaterialViewAccessor< array1d< real64 >, arrayView1d< real64 const > >( LinearElasticIsotropic::viewKeyStruct::bulkModulusString,
                                                                                                  m_solidSolver->targetRegionNames(),
                                                                                                  m_solidSolver->solidMaterialNames() );
  // Shear modulus accessor
  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > const shearModulus =
    elemManager->ConstructMaterialViewAccessor< array1d< real64 >, arrayView1d< real64 const > >( LinearElasticIsotropic::viewKeyStruct::shearModulusString,
                                                                                                  m_solidSolver->targetRegionNames(),
                                                                                                  m_solidSolver->solidMaterialNames() );
  arrayView1d< globalIndex const > const &
  tracDofNumber = fractureSubRegion->getReference< globalIndex_array >( tracDofKey );

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
        real64 rotatedInvStiffApprox[ 2 ][ 3 ][ 3 ];
        for( localIndex kf = 0; kf < 2; ++kf )
        {
          // Get fracture, face and region/subregion/element indices (for elements on both sides)
          localIndex fractureIndex = sei[iconn][kf];

          localIndex faceIndexRef = faceMap[fractureIndex][0];
          real64 const area = faceArea[faceIndexRef];
          // TODO: use higher order integration scheme
          nodalArea[kf][0] = area / 4.0;
          nodalArea[kf][1] = area / 4.0;

          real64 invStiffApprox[ 2 ][ 3 ];
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
            real64 const K = bulkModulus[er][esr][ei];
            real64 const G = shearModulus[er][esr][ei];
            real64 const E = 9.0 * K * G / ( 3.0 * K + G );
            real64 const nu = ( 3.0 * K - 2.0 * G ) / ( 2.0 * ( 3.0 * K + G ) );

            // The factor is 8/9 / 4 (number of nodes) = 2/9
            for( localIndex j = 0; j < 3; ++j )
            {
              invStiffApprox[ i ][ j ] = 1.0 / ( E / ( ( 1.0 + nu )*( 1.0 - 2.0*nu ) ) * 2.0 / 9.0 * ( 2.0 - 3.0 * nu ) * volume / ( boxSize[j]*boxSize[j] ) );
            }
          }

          real64 invStiffApproxTotal[ 3 ][ 3 ] = { { 0 } };
          for( localIndex i = 0; i < 2; ++i )
          {
            for( localIndex j = 0; j < 3; ++j )
            {
              invStiffApproxTotal[ j ][ j ] += invStiffApprox[ i ][ j ];
            }
          }

          // Compute R^T * (invK) * R
          real64 temp[ 3 ][ 3 ];
          LvArray::tensorOps::AkiBkj< 3, 3, 3 >( temp, faceRotationMatrix[ faceIndexRef ], invStiffApproxTotal );
          LvArray::tensorOps::AikBkj< 3, 3, 3 >( rotatedInvStiffApprox[ kf ], temp, faceRotationMatrix[ faceIndexRef ] );
        }

        // Compose local nodal-based local stiffness matrices
        stackArray2d< real64, 3*3 > totalInvStiffApprox( 3, 3 );
        for( localIndex kf = 0; kf < 2; ++kf )
        {
          for( localIndex i = 0; i < 3; ++i )
          {
            for( localIndex j = 0; j < 3; ++j )
            {
              rotatedInvStiffApprox[ kf ][ i ][ j ] *= nodalArea[0][kf] * nodalArea[1][kf];
            }
          }
        }
        // Local assembly
        for( localIndex i = 0; i < 3; ++i )
        {
          for( localIndex j = 0; j < 3; ++j )
          {
            totalInvStiffApprox( i, j ) = -( rotatedInvStiffApprox[ 0 ][ i ][ j ] + rotatedInvStiffApprox[ 1 ][ i ][ j ] );
          }
        }

        // Get DOF numbering
        localIndex fractureIndex[2];
        localIndex nDof[2];
        globalIndex elemDOF[2][3];
        for( localIndex kf = 0; kf < 2; ++kf )
        {
          fractureIndex[kf] = sei[iconn][kf];
          for( localIndex i=0; i<3; ++i )
          {
            elemDOF[kf][i] = tracDofNumber[fractureIndex[kf]] + LvArray::integerConversion< globalIndex >( i );
          }
          nDof[kf] = 0;
          switch( fractureState[fractureIndex[kf]] )
          {
            case ( FractureState::STICK ):
              {
                nDof[kf] = 3;
                break;
              }
            case ( FractureState::NEW_SLIP ):
            case ( FractureState::SLIP ):
              {
                nDof[kf] = 1;
                break;
              }
            case ( FractureState::OPEN ):
              {
                nDof[kf] = 0;
                break;
              }
          }
        }

        // Define local "transmissibility" matrices
        stackArray2d< real64, 3*3 > totalInvStiffApprox00( nDof[0], nDof[0] );
        stackArray2d< real64, 3*3 > totalInvStiffApprox01( nDof[0], nDof[1] );
        stackArray2d< real64, 3*3 > totalInvStiffApprox10( nDof[1], nDof[0] );
        stackArray2d< real64, 3*3 > totalInvStiffApprox11( nDof[1], nDof[1] );
        for( localIndex i = 0; i < nDof[0]; ++i )
        {
          for( localIndex j = 0; j < nDof[0]; ++j )
          {
            totalInvStiffApprox00( i, j ) = totalInvStiffApprox( i, j );
          }
          for( localIndex j = 0; j < nDof[1]; ++j )
          {
            totalInvStiffApprox01( i, j ) = -totalInvStiffApprox( i, j );
          }
        }

        for( localIndex i = 0; i < nDof[1]; ++i )
        {
          for( localIndex j = 0; j < nDof[0]; ++j )
          {
            totalInvStiffApprox10( i, j ) = -totalInvStiffApprox( i, j );
          }
          for( localIndex j = 0; j < nDof[1]; ++j )
          {
            totalInvStiffApprox11( i, j ) = totalInvStiffApprox( i, j );
          }
        }

        // Compute rhs
        real64_array rhs0( 3 );
        rhs0 = 0.0;
        if( nDof[0] > 0 )
        {
          for( localIndex j = 0; j < nDof[0]; ++j )
          {
            rhs0( 0 ) += totalInvStiffApprox00( 0, j ) * ( traction[fractureIndex[0]][j] );
          }
          for( localIndex j = 0; j < nDof[1]; ++j )
          {
            rhs0( 0 ) += totalInvStiffApprox01( 0, j ) * ( traction[fractureIndex[1]][j] );
          }
          for( localIndex i = 1; i < nDof[0]; ++i )
          {
            for( localIndex j = 0; j < nDof[0]; ++j )
            {
              rhs0( i ) += totalInvStiffApprox00( i, j ) * ( deltaTraction[fractureIndex[0]][j] );
            }
            for( localIndex j = 0; j < nDof[1]; ++j )
            {
              rhs0( i ) += totalInvStiffApprox01( i, j ) * ( deltaTraction[fractureIndex[1]][j] );
            }
          }
        }

        real64_array rhs1( 3 );
        rhs1 = 0.0;
        if( nDof[1] > 0 )
        {
          for( localIndex j = 0; j < nDof[0]; ++j )
          {
            rhs1( 0 ) += totalInvStiffApprox10( 0, j ) * ( traction[fractureIndex[0]][j] );
          }
          for( localIndex j = 0; j < nDof[1]; ++j )
          {
            rhs1( 0 ) += totalInvStiffApprox11( 0, j ) * ( traction[fractureIndex[1]][j] );
          }
          for( localIndex i = 1; i < nDof[1]; ++i )
          {
            for( localIndex j = 0; j < nDof[0]; ++j )
            {
              rhs1( i ) += totalInvStiffApprox10( i, j ) * ( deltaTraction[fractureIndex[0]][j] );
            }
            for( localIndex j = 0; j < nDof[1]; ++j )
            {
              rhs1( i ) += totalInvStiffApprox11( i, j ) * ( deltaTraction[fractureIndex[1]][j] );
            }
          }
        }

        // Global matrix and rhs assembly
        if( std::max( nDof[0], nDof[1] ) > 0 )
        {
          matrix->add( elemDOF[0],
                       elemDOF[0],
                       totalInvStiffApprox00.data(),
                       nDof[0],
                       nDof[0] );

          matrix->add( elemDOF[1],
                       elemDOF[1],
                       totalInvStiffApprox11.data(),
                       nDof[1],
                       nDof[1] );

          matrix->add( elemDOF[0],
                       elemDOF[1],
                       totalInvStiffApprox01.data(),
                       nDof[0],
                       nDof[1] );

          matrix->add( elemDOF[1],
                       elemDOF[0],
                       totalInvStiffApprox10.data(),
                       nDof[1],
                       nDof[0] );

          rhs->add( elemDOF[0],
                    rhs0.data(),
                    nDof[0] );
          rhs->add( elemDOF[1],
                    rhs1.data(),
                    nDof[1] );
        }
      }
    }
  } );

  matrix->close();
  rhs->close();
}

void LagrangianContactSolver::ApplySystemSolution( DofManager const & dofManager,
                                                   ParallelVector const & solution,
                                                   real64 const scalingFactor,
                                                   DomainPartition * const domain )
{
  GEOSX_MARK_FUNCTION;

  m_solidSolver->ApplySystemSolution( dofManager,
                                      solution,
                                      scalingFactor,
                                      domain );

  dofManager.addVectorToField( solution, viewKeyStruct::tractionString, viewKeyStruct::deltaTractionString, -scalingFactor );
  dofManager.addVectorToField( solution, viewKeyStruct::tractionString, viewKeyStruct::tractionString, -scalingFactor );

  std::map< string, string_array > fieldNames;
  fieldNames["elems"].push_back( viewKeyStruct::tractionString );
  fieldNames["elems"].push_back( viewKeyStruct::deltaTractionString );
  // This is used locally only, synchronized just for output reasons
  fieldNames["elems"].push_back( viewKeyStruct::localJumpString );
  // fractureStateString is synchronized in UpdateFractureState
  // previousFractureStateString and previousLocalJumpString used locally only

  CommunicationTools::SynchronizeFields( fieldNames,
                                         domain->getMeshBody( 0 )->getMeshLevel( 0 ),
                                         domain->getNeighbors() );

  this->UpdateDeformationForCoupling( domain );
}

void LagrangianContactSolver::InitializeFractureState( MeshLevel * const mesh,
                                                       string const fieldName ) const
{
  GEOSX_MARK_FUNCTION;
  ElementRegionManager * const elemManager = mesh->getElemManager();

  elemManager->forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion & subRegion )
  {
    if( subRegion.hasWrapper( m_tractionKey ) )
    {
      arrayView1d< integer > const & fractureState = subRegion.getReference< array1d< integer > >( fieldName );
      forAll< serialPolicy >( subRegion.size(), [&]( localIndex const kfe )
      {
        fractureState[kfe] = FractureState::STICK;
      } );
    }
  } );
}

void LagrangianContactSolver::SetFractureStateForElasticStep( DomainPartition * const domain ) const
{
  GEOSX_MARK_FUNCTION;

  MeshLevel * const mesh = domain->getMeshBody( 0 )->getMeshLevel( 0 );
  ElementRegionManager * const elemManager = mesh->getElemManager();

  elemManager->forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion & subRegion )
  {
    if( subRegion.hasWrapper( m_tractionKey ) )
    {
      arrayView1d< integer const > const & ghostRank = subRegion.ghostRank();
      arrayView1d< integer > const & fractureState = subRegion.getReference< array1d< integer > >( viewKeyStruct::fractureStateString );
      forAll< serialPolicy >( subRegion.size(), [&]( localIndex const kfe )
      {
        if( ghostRank[kfe] < 0 )
        {
          if( fractureState[kfe] != FractureState::OPEN )
          {
            fractureState[kfe] = FractureState::STICK;
          }
        }
      } );
    }
  } );
}

bool LagrangianContactSolver::UpdateFractureState( DomainPartition * const domain ) const
{
  GEOSX_MARK_FUNCTION;

  MeshLevel * const mesh = domain->getMeshBody( 0 )->getMeshLevel( 0 );
  ElementRegionManager * const elemManager = mesh->getElemManager();

  ConstitutiveManager const * const
  constitutiveManager = domain->GetGroup< ConstitutiveManager >( keys::ConstitutiveManager );
  ContactRelationBase const * const
  contactRelation = constitutiveManager->GetGroup< ContactRelationBase const >( m_contactRelationName );

  bool checkActiveSet = true;

  elemManager->forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion & subRegion )
  {
    if( subRegion.hasWrapper( m_tractionKey ) )
    {
      arrayView1d< integer const > const & ghostRank = subRegion.ghostRank();
      arrayView2d< real64 const > const & traction = subRegion.getReference< array2d< real64 > >( viewKeyStruct::tractionString );
      arrayView2d< real64 const > const & localJump = subRegion.getReference< array2d< real64 > >( viewKeyStruct::localJumpString );
      arrayView1d< integer > const & fractureState = subRegion.getReference< array1d< integer > >( viewKeyStruct::fractureStateString );
      arrayView1d< real64 const > const &
      normalTractionTolerance = subRegion.getReference< array1d< real64 > >( viewKeyStruct::normalTractionToleranceString );
      arrayView1d< real64 const > const &
      normalDisplacementTolerance = subRegion.getReference< array1d< real64 > >( viewKeyStruct::normalDisplacementToleranceString );

      forAll< serialPolicy >( subRegion.size(), [&]( localIndex const kfe )
      {
        if( ghostRank[kfe] < 0 )
        {
          integer const originalFractureState = fractureState[kfe];
          if( originalFractureState == FractureState::OPEN )
          {
            if( localJump[kfe][0] > -normalDisplacementTolerance[kfe] )
            {
              fractureState[kfe] = FractureState::OPEN;
            }
            else
            {
              fractureState[kfe] = FractureState::STICK;
            }
          }
          else if( traction[kfe][0] > normalTractionTolerance[kfe] )
          {
            fractureState[kfe] = FractureState::OPEN;
          }
          else
          {
            real64 currentTau = sqrt( traction[kfe][1]*traction[kfe][1] + traction[kfe][2]*traction[kfe][2] );
            real64 limitTau = contactRelation->limitTangentialTractionNorm( traction[kfe][0] );
            if( originalFractureState == FractureState::STICK && currentTau >= limitTau )
            {
              currentTau *= (1.0 - m_slidingCheckTolerance);
            }
            else if( originalFractureState != FractureState::STICK && currentTau <= limitTau )
            {
              currentTau *= (1.0 + m_slidingCheckTolerance);
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

          if( originalFractureState != fractureState[kfe] )
          {
            GEOSX_LOG_LEVEL_BY_RANK( 3, "element " << kfe << " traction: " << traction[kfe]
                                                   << " previous state <"
                                                   << FractureStateToString( originalFractureState )
                                                   << "> current state <"
                                                   << FractureStateToString( fractureState[kfe] )
                                                   << ">" );
          }
          checkActiveSet &= CompareFractureStates( originalFractureState, fractureState[kfe] );
        }
      } );
    }
  } );

  // Need to synchronize the fracture state due to the use will be made of in AssemblyStabilization
  SynchronizeFractureState( domain );

  // Compute if globally the fracture state has changed
  bool globalCheckActiveSet;
  MpiWrapper::allReduce( &checkActiveSet,
                         &globalCheckActiveSet,
                         1,
                         MPI_LAND,
                         MPI_COMM_GEOSX );

  return globalCheckActiveSet;
}

void LagrangianContactSolver::SynchronizeFractureState( DomainPartition * const domain ) const
{
  std::map< string, string_array > fieldNames;
  fieldNames["elems"].push_back( viewKeyStruct::fractureStateString );

  CommunicationTools::SynchronizeFields( fieldNames,
                                         domain->getMeshBody( 0 )->getMeshLevel( 0 ),
                                         domain->getNeighbors() );
}

bool LagrangianContactSolver::IsFractureAllInStickCondition( DomainPartition const * const domain ) const
{
  globalIndex numStick, numSlip, numOpen;
  ComputeFractureStateStatistics( domain, numStick, numSlip, numOpen, false );
  return ( ( numSlip + numOpen ) == 0 );
}

void LagrangianContactSolver::ComputeFractureStateStatistics( DomainPartition const * const domain,
                                                              globalIndex & numStick,
                                                              globalIndex & numSlip,
                                                              globalIndex & numOpen,
                                                              bool printAll ) const
{
  MeshLevel const * const mesh = domain->getMeshBody( 0 )->getMeshLevel( 0 );
  ElementRegionManager const * const elemManager = mesh->getElemManager();

  globalIndex_array localCounter( 3 );
  localCounter = 0;

  elemManager->forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion const & subRegion )
  {
    if( subRegion.hasWrapper( m_tractionKey ) )
    {
      arrayView1d< integer const > const & ghostRank = subRegion.ghostRank();
      arrayView1d< integer const > const & fractureState = subRegion.getReference< array1d< integer > >( viewKeyStruct::fractureStateString );
      arrayView2d< real64 const > const & traction = subRegion.getReference< array2d< real64 > >( viewKeyStruct::tractionString );

      forAll< serialPolicy >( subRegion.size(), [&]( localIndex const kfe )
      {
        if( ghostRank[kfe] < 0 )
        {
          switch( fractureState[kfe] )
          {
            case FractureState::STICK:
              {
                localCounter[0] += 1;
                break;
              }
            case FractureState::NEW_SLIP:
            case FractureState::SLIP:
              {
                localCounter[1] += 1;
                break;
              }
            case FractureState::OPEN:
              {
                localCounter[2] += 1;
                break;
              }
          }
          if( printAll )
          {
            GEOSX_LOG_LEVEL_BY_RANK( 3, "element " << kfe << " traction: " << traction[kfe]
                                                   << " state <"
                                                   << FractureStateToString( fractureState[kfe] )
                                                   << ">" );
          }
        }
      } );
    }
  } );

  globalIndex_array totalCounter( 3 );
  totalCounter = 0;
  int const rank = MpiWrapper::Comm_rank( MPI_COMM_GEOSX );
  int const size = MpiWrapper::Comm_size( MPI_COMM_GEOSX );
  globalIndex_array globalCounter( 3*size );
  globalCounter = 0;

  // Everything is done on rank 0
  MpiWrapper::gather( localCounter.data(),
                      3,
                      globalCounter.data(),
                      3,
                      0,
                      MPI_COMM_GEOSX );

  if( rank==0 )
  {
    for( int r=0; r<size; ++r )
    {
      // sum across all ranks
      totalCounter[0] += globalCounter[3*r];
      totalCounter[1] += globalCounter[3*r+1];
      totalCounter[2] += globalCounter[3*r+2];
    }
  }

  MpiWrapper::bcast( totalCounter.data(), 3, 0, MPI_COMM_GEOSX );

  numStick = totalCounter[0];
  numSlip  = totalCounter[1];
  numOpen  = totalCounter[2];

  char output[108] = {0};
  sprintf( output,
           " Number of element for each fracture state:"
           " stick: %12lli | slip:  %12lli | open:  %12lli",
           numStick,
           numSlip,
           numOpen );
  GEOSX_LOG_RANK_0( output );
}

void LagrangianContactSolver::SolveSystem( DofManager const & dofManager,
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

//  int rank = MpiWrapper::Comm_rank( MPI_COMM_GEOSX );
//  if( rank == 0 )
//  {
//    string str;
//    std::getline( std::cin, str );
//    if( str.length() > 0 )
//    {
//      GEOSX_ERROR( "STOP" );
//    }
//  }
//  MpiWrapper::Barrier( MPI_COMM_GEOSX );
}

void LagrangianContactSolver::SetNextDt( real64 const & currentDt,
                                         real64 & nextDt )
{
  nextDt = currentDt;
}

REGISTER_CATALOG_ENTRY( SolverBase, LagrangianContactSolver, std::string const &, Group * const )
} /* namespace geosx */
