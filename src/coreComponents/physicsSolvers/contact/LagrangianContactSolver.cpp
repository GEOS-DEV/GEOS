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
 * @file LagrangianContactSolver.cpp
 *
 */

#include "LagrangianContactSolver.hpp"

#include "common/TimingMacros.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/contact/ContactSelector.hpp"
#include "constitutive/fluid/singlefluid/SingleFluidBase.hpp"
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "mesh/DomainPartition.hpp"
#include "discretizationMethods/NumericalMethodsManager.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mesh/SurfaceElementRegion.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "mesh/mpiCommunications/NeighborCommunicator.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp" // needed to register pressure(_n)
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"
#include "physicsSolvers/surfaceGeneration/SurfaceGenerator.hpp"
#include "physicsSolvers/contact/ContactFields.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "linearAlgebra/utilities/LAIHelperFunctions.hpp"
#include "linearAlgebra/solvers/PreconditionerJacobi.hpp"
#include "linearAlgebra/solvers/PreconditionerBlockJacobi.hpp"
#include "linearAlgebra/solvers/BlockPreconditioner.hpp"
#include "linearAlgebra/solvers/SeparateComponentPreconditioner.hpp"
#include "finiteElement/elementFormulations/H1_TriangleFace_Lagrange1_Gauss1.hpp"
#include "finiteElement/elementFormulations/H1_QuadrilateralFace_Lagrange1_GaussLegendre2.hpp"

#if defined( __INTEL_COMPILER )
#pragma GCC optimize "O0"
#endif

namespace geos
{

using namespace constitutive;
using namespace dataRepository;
using namespace fields;
using namespace finiteElement;

LagrangianContactSolver::LagrangianContactSolver( const string & name,
                                                  Group * const parent ):
  ContactSolverBase( name, parent )
{
  registerWrapper( viewKeyStruct::stabilizationNameString(), &m_stabilizationName ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of the stabilization to use in the lagrangian contact solver" );

  LinearSolverParameters & linSolParams = m_linearSolverParameters.get();
  linSolParams.mgr.strategy = LinearSolverParameters::MGR::StrategyType::lagrangianContactMechanics;
  linSolParams.mgr.separateComponents = true;
  linSolParams.mgr.displacementFieldName = solidMechanics::totalDisplacement::key();
  linSolParams.dofsPerNode = 3;
}

void LagrangianContactSolver::registerDataOnMesh( Group & meshBodies )
{
  ContactSolverBase::registerDataOnMesh( meshBodies );

  forFractureRegionOnMeshTargets( meshBodies, [&] ( SurfaceElementRegion & fractureRegion )
  {
    fractureRegion.forElementSubRegions< SurfaceElementSubRegion >( [&]( SurfaceElementSubRegion & subRegion )
    {
      subRegion.registerWrapper< array3d< real64 > >( viewKeyStruct::rotationMatrixString() ).
        setPlotLevel( PlotLevel::NOPLOT ).
        setRegisteringObjects( this->getName()).
        setDescription( "An array that holds the rotation matrices on the fracture." ).
        reference().resizeDimension< 1, 2 >( 3, 3 );

      subRegion.registerField< fields::contact::deltaTraction >( getName() ).
        reference().resizeDimension< 1 >( 3 );

      subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::normalTractionToleranceString() ).
        setPlotLevel( PlotLevel::NOPLOT ).
        setRegisteringObjects( this->getName()).
        setDescription( "An array that holds the normal traction tolerance." );

      subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::normalDisplacementToleranceString() ).
        setPlotLevel( PlotLevel::NOPLOT ).
        setRegisteringObjects( this->getName()).
        setDescription( "An array that holds the normal displacement tolerance." );

      subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::slidingToleranceString() ).
        setPlotLevel( PlotLevel::NOPLOT ).
        setRegisteringObjects( this->getName()).
        setDescription( "An array that holds the sliding tolerance." );

      // Needed just because SurfaceGenerator initialize the field "pressure" (NEEDED!!!)
      // It is used in "TwoPointFluxApproximation.cpp", called by "SurfaceGenerator.cpp"
      subRegion.registerField< flow::pressure >( getName() ).
        setPlotLevel( PlotLevel::NOPLOT ).
        setRegisteringObjects( this->getName());
      subRegion.registerField< flow::pressure_n >( getName() ).
        setPlotLevel( PlotLevel::NOPLOT ).
        setRegisteringObjects( this->getName());

    } );

    forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                      MeshLevel & mesh,
                                                      arrayView1d< string const > const & )
    {
      FaceManager & faceManager = mesh.getFaceManager();

      faceManager.registerWrapper< array1d< real64 > >( viewKeyStruct::transMultiplierString() ).
        setApplyDefaultValue( 1.0 ).
        setPlotLevel( PlotLevel::LEVEL_0 ).
        setRegisteringObjects( this->getName() ).
        setDescription( "An array that holds the permeability transmissibility multipliers" );
    } );

  } );
}

void LagrangianContactSolver::initializePreSubGroups()
{
  ContactSolverBase::initializePreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
  //ConstitutiveManager const & cm = domain.getConstitutiveManager();

  // fill stencil targetRegions
  NumericalMethodsManager & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager & fvManager = numericalMethodManager.getFiniteVolumeManager();

  if( fvManager.hasGroup< FluxApproximationBase >( m_stabilizationName ) )
  {

    FluxApproximationBase & fluxApprox = fvManager.getFluxApproximation( m_stabilizationName );
    fluxApprox.addFieldName( contact::traction::key() );
    fluxApprox.setCoeffName( "penaltyStiffness" );

    forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const & meshBodyName,
                                                                  MeshLevel &,
                                                                  arrayView1d< string const > const & regionNames )
    {
      array1d< string > & stencilTargetRegions = fluxApprox.targetRegions( meshBodyName );
      std::set< string > stencilTargetRegionsSet( stencilTargetRegions.begin(), stencilTargetRegions.end() );
      stencilTargetRegionsSet.insert( regionNames.begin(), regionNames.end() );
      stencilTargetRegions.clear();
      for( auto const & targetRegion: stencilTargetRegionsSet )
      {
        stencilTargetRegions.emplace_back( targetRegion );
      }
    } );
  }

}

void LagrangianContactSolver::setupSystem( DomainPartition & domain,
                                           DofManager & dofManager,
                                           CRSMatrix< real64, globalIndex > & localMatrix,
                                           ParallelVector & rhs,
                                           ParallelVector & solution,
                                           bool const GEOS_UNUSED_PARAM( setSparsity ) )
{
  if( m_precond )
  {
    m_precond->clear();
  }

  // setup monolithic coupled system
  SolverBase::setupSystem( domain, dofManager, localMatrix, rhs, solution, true );

  if( !m_precond && m_linearSolverParameters.get().solverType != LinearSolverParameters::SolverType::direct )
  {
    createPreconditioner( domain );
  }
}

void LagrangianContactSolver::implicitStepSetup( real64 const & time_n,
                                                 real64 const & dt,
                                                 DomainPartition & domain )
{
  computeRotationMatrices( domain );
  computeTolerances( domain );
  computeFaceDisplacementJump( domain );

  SolidMechanicsLagrangianFEM::implicitStepSetup( time_n, dt, domain );
}

void LagrangianContactSolver::implicitStepComplete( real64 const & time_n,
                                                    real64 const & dt,
                                                    DomainPartition & domain )
{
  //if( m_setupSolidSolverDofs )
  {
    SolidMechanicsLagrangianFEM::implicitStepComplete( time_n, dt, domain );
  }

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & )
  {
    mesh.getElemManager().forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion & subRegion )
    {
      arrayView2d< real64 > const & deltaTraction = subRegion.getField< contact::deltaTraction >();
      arrayView2d< real64 const > const & dispJump = subRegion.getField< contact::dispJump >();
      arrayView2d< real64 > const & oldDispJump = subRegion.getField< contact::oldDispJump >();
      arrayView1d< integer const > const & fractureState = subRegion.getField< contact::fractureState >();
      arrayView1d< integer > const & oldFractureState = subRegion.getField< contact::oldFractureState >();

      forAll< parallelHostPolicy >( subRegion.size(), [=] ( localIndex const kfe )
      {
        for( localIndex i = 0; i < 3; ++i )
        {
          deltaTraction[kfe][i] = 0.0;
          oldDispJump[kfe][i] = dispJump[kfe][i];
        }
        oldFractureState[kfe] = fractureState[kfe];
      } );
    } );

    // Need a synchronization of deltaTraction as will be used in AssembleStabilization
    FieldIdentifiers fieldsToBeSync;
    fieldsToBeSync.addElementFields( { contact::deltaTraction::key() },
                                     { getFractureRegionName() } );

    CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync,
                                                         mesh,
                                                         domain.getNeighbors(),
                                                         true );

  } );
}

LagrangianContactSolver::~LagrangianContactSolver()
{
  // TODO Auto-generated destructor stub
}

void LagrangianContactSolver::computeTolerances( DomainPartition & domain ) const
{
  GEOS_MARK_FUNCTION;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & )
  {
    FaceManager const & faceManager = mesh.getFaceManager();
    NodeManager const & nodeManager = mesh.getNodeManager();
    ElementRegionManager & elemManager = mesh.getElemManager();

    // Get the "face to element" map (valid for the entire mesh)
    FaceManager::ElemMapType const & faceToElem = faceManager.toElementRelation();
    arrayView2d< localIndex const > const & faceToElemRegion = faceToElem.m_toElementRegion;
    arrayView2d< localIndex const > const & faceToElemSubRegion = faceToElem.m_toElementSubRegion;
    arrayView2d< localIndex const > const & faceToElemIndex = faceToElem.m_toElementIndex;

    // Get the volume for all elements
    ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > const elemVolume =
      elemManager.constructViewAccessor< array1d< real64 >, arrayView1d< real64 const > >( ElementSubRegionBase::viewKeyStruct::elementVolumeString() );

    // Get the coordinates for all nodes
    arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition = nodeManager.referencePosition();

    // Bulk modulus accessor
    ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > const bulkModulus =
      elemManager.constructMaterialViewAccessor< ElasticIsotropic, array1d< real64 >, arrayView1d< real64 const > >( ElasticIsotropic::viewKeyStruct::bulkModulusString() );
    // Shear modulus accessor
    ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > const shearModulus =
      elemManager.constructMaterialViewAccessor< ElasticIsotropic, array1d< real64 >, arrayView1d< real64 const > >( ElasticIsotropic::viewKeyStruct::shearModulusString() );

    using NodeMapViewType = arrayView2d< localIndex const, cells::NODE_MAP_USD >;
    ElementRegionManager::ElementViewAccessor< NodeMapViewType > const elemToNode =
      elemManager.constructViewAccessor< CellElementSubRegion::NodeMapType, NodeMapViewType >( ElementSubRegionBase::viewKeyStruct::nodeListString() );
    ElementRegionManager::ElementViewConst< NodeMapViewType > const elemToNodeView = elemToNode.toNestedViewConst();

    elemManager.forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion & subRegion )
    {
      if( subRegion.hasField< contact::traction >() )
      {
        arrayView1d< integer const > const & ghostRank = subRegion.ghostRank();
        arrayView1d< real64 const > const & faceArea = subRegion.getElementArea().toViewConst();
        arrayView3d< real64 const > const & faceRotationMatrix = subRegion.getReference< array3d< real64 > >( viewKeyStruct::rotationMatrixString() );
        ArrayOfArraysView< localIndex const > const & elemsToFaces = subRegion.faceList().toViewConst();

        arrayView1d< real64 > const & normalTractionTolerance =
          subRegion.getReference< array1d< real64 > >( viewKeyStruct::normalTractionToleranceString() );
        arrayView1d< real64 > const & normalDisplacementTolerance =
          subRegion.getReference< array1d< real64 > >( viewKeyStruct::normalDisplacementToleranceString() );
        arrayView1d< real64 > const & slidingTolerance =
          subRegion.getReference< array1d< real64 > >( viewKeyStruct::slidingToleranceString() );

        forAll< parallelHostPolicy >( subRegion.size(), [=] ( localIndex const kfe )
        {
          if( ghostRank[kfe] < 0 )
          {
            real64 const area = faceArea[kfe];
            // approximation of the stiffness along coordinate directions
            // ( first, second ) index -> ( element index, direction )
            // 1. T -> top (index 0), B -> bottom (index 1)
            // 2. the coordinate direction (x, y, z)
            real64 stiffDiagApprox[ 2 ][ 3 ];
            real64 averageYoungModulus = 0.0;
            real64 averageConstrainedModulus = 0.0;
            real64 averageBoxSize0 = 0.0;

            for( localIndex i = 0; i < elemsToFaces.sizeOfArray( kfe ); ++i )
            {
              localIndex const faceIndex = elemsToFaces[kfe][i];
              localIndex const er = faceToElemRegion[faceIndex][0];
              localIndex const esr = faceToElemSubRegion[faceIndex][0];
              localIndex const ei = faceToElemIndex[faceIndex][0];

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
              real64 const M = K + 4.0 / 3.0 * G;

              // Combine E and nu to obtain a stiffness approximation (like it was an hexahedron)
              for( localIndex j = 0; j < 3; ++j )
              {
                stiffDiagApprox[ i ][ j ] = E / ( ( 1.0 + nu )*( 1.0 - 2.0*nu ) ) * 4.0 / 9.0 * ( 2.0 - 3.0 * nu ) * volume / ( boxSize[j]*boxSize[j] );
              }

              averageYoungModulus += 0.5*E;
              averageConstrainedModulus += 0.5*M;
              averageBoxSize0 += 0.5*boxSize[0];
            }

            // Average the stiffness and compute the inverse
            real64 invStiffApprox[ 3 ][ 3 ] = { { 0 } };
            for( localIndex j = 0; j < 3; ++j )
            {
              invStiffApprox[ j ][ j ] = ( stiffDiagApprox[ 0 ][ j ] + stiffDiagApprox[ 1 ][ j ] ) / ( stiffDiagApprox[ 0 ][ j ] * stiffDiagApprox[ 1 ][ j ] );
            }

            // Rotate in the local reference system, computing R^T * (invK) * R
            real64 temp[ 3 ][ 3 ];
            LvArray::tensorOps::Rij_eq_AkiBkj< 3, 3, 3 >( temp, faceRotationMatrix[ kfe ], invStiffApprox );
            real64 rotatedInvStiffApprox[ 3 ][ 3 ];
            LvArray::tensorOps::Rij_eq_AikBkj< 3, 3, 3 >( rotatedInvStiffApprox, temp, faceRotationMatrix[ kfe ] );
            LvArray::tensorOps::scale< 3, 3 >( rotatedInvStiffApprox, area );

            // Finally, compute tolerances for the given fracture element
            normalDisplacementTolerance[kfe] = rotatedInvStiffApprox[ 0 ][ 0 ] * averageYoungModulus / 2.e+7;
            slidingTolerance[kfe] = sqrt( rotatedInvStiffApprox[ 1 ][ 1 ] * rotatedInvStiffApprox[ 1 ][ 1 ] +
                                          rotatedInvStiffApprox[ 2 ][ 2 ] * rotatedInvStiffApprox[ 2 ][ 2 ] ) * averageYoungModulus / 2.e+7;
            normalTractionTolerance[kfe] = 1.0 / 2.0 * averageConstrainedModulus / averageBoxSize0 * normalDisplacementTolerance[kfe];
          }
        } );
      }
    } );
  } );
}

void LagrangianContactSolver::resetStateToBeginningOfStep( DomainPartition & domain )
{
  SolidMechanicsLagrangianFEM::resetStateToBeginningOfStep( domain );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & )
  {
    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion & subRegion )
    {
      arrayView2d< real64 > const & traction = subRegion.getField< contact::traction >();
      arrayView2d< real64 > const & deltaTraction = subRegion.getField< contact::deltaTraction >();
      arrayView2d< real64 > const & dispJump = subRegion.getField< contact::dispJump >();
      arrayView2d< real64 const > const & oldDispJump = subRegion.getField< contact::oldDispJump >();

      arrayView1d< integer > const & fractureState = subRegion.getField< contact::fractureState >();
      arrayView1d< integer const > const & oldFractureState = subRegion.getField< contact::oldFractureState >();

      forAll< parallelHostPolicy >( subRegion.size(), [=] ( localIndex const kfe )
      {
        for( localIndex i = 0; i < 3; ++i )
        {
          traction[kfe][i] -= deltaTraction[kfe][i];
          deltaTraction[kfe][i] = 0.0;

          dispJump[kfe][i] = oldDispJump[kfe][i];
        }
        fractureState[kfe] = oldFractureState[kfe];
      } );
    } );
  } );
}

void LagrangianContactSolver::computeFaceDisplacementJump( DomainPartition & domain )
{
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    NodeManager const & nodeManager = mesh.getNodeManager();
    FaceManager & faceManager = mesh.getFaceManager();
    ElementRegionManager & elemManager = mesh.getElemManager();

    // Get the coordinates for all nodes
    arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition = nodeManager.referencePosition();

    ArrayOfArraysView< localIndex const > const faceToNodeMap = faceManager.nodeList().toViewConst();

    arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const & u =
      nodeManager.getField< solidMechanics::totalDisplacement >();

    elemManager.forElementSubRegions< FaceElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   FaceElementSubRegion & subRegion )
    {
      if( subRegion.hasField< contact::traction >() )
      {
        arrayView3d< real64 > const &
        rotationMatrix = subRegion.getReference< array3d< real64 > >( viewKeyStruct::rotationMatrixString() );
        ArrayOfArraysView< localIndex const > const & elemsToFaces = subRegion.faceList().toViewConst();
        arrayView2d< real64 > const & dispJump = subRegion.getField< contact::dispJump >();
        arrayView1d< real64 const > const & area = subRegion.getElementArea().toViewConst();

        forAll< parallelHostPolicy >( subRegion.size(), [=] ( localIndex const kfe )
        {
          if( elemsToFaces.sizeOfArray( kfe ) != 2 )
          {
            return;
          }

          // Contact constraints
          localIndex const numNodesPerFace = faceToNodeMap.sizeOfArray( elemsToFaces[kfe][0] );

          array1d< real64 > nodalArea0, nodalArea1;
          computeFaceNodalArea( nodePosition, faceToNodeMap, elemsToFaces[kfe][0], nodalArea0 );
          computeFaceNodalArea( nodePosition, faceToNodeMap, elemsToFaces[kfe][1], nodalArea1 );

          real64 globalJumpTemp[ 3 ] = { 0 };
          for( localIndex a = 0; a < numNodesPerFace; ++a )
          {
            for( localIndex i = 0; i < 3; ++i )
            {
              globalJumpTemp[ i ] +=
                ( -u[faceToNodeMap( elemsToFaces[kfe][0], a )][i] * nodalArea0[a]
                  +u[faceToNodeMap( elemsToFaces[kfe][1], a )][i] * nodalArea1[a] ) / area[kfe];
            }
          }

          real64 dispJumpTemp[ 3 ];
          LvArray::tensorOps::Ri_eq_AjiBj< 3, 3 >( dispJumpTemp, rotationMatrix[ kfe ], globalJumpTemp );
          LvArray::tensorOps::copy< 3 >( dispJump[ kfe ], dispJumpTemp );
        } );
      }
    } );
  } );
}

void LagrangianContactSolver::setupDofs( DomainPartition const & domain,
                                         DofManager & dofManager ) const
{
  GEOS_MARK_FUNCTION;
  //if( m_setupSolidSolverDofs )
  {
    SolidMechanicsLagrangianFEM::setupDofs( domain, dofManager );
  }
  // restrict coupling to fracture regions only
  map< std::pair< string, string >, array1d< string > > meshTargets;
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const & meshBodyName,
                                                                MeshLevel const & meshLevel,
                                                                arrayView1d< string const > const & regionNames )
  {
    array1d< string > regions;
    ElementRegionManager const & elementRegionManager = meshLevel.getElemManager();
    elementRegionManager.forElementRegions< SurfaceElementRegion >( regionNames,
                                                                    [&]( localIndex const,
                                                                         SurfaceElementRegion const & region )
    {
      regions.emplace_back( region.getName() );
    } );
    meshTargets[std::make_pair( meshBodyName, meshLevel.getName())] = std::move( regions );
  } );

  dofManager.addField( contact::traction::key(),
                       FieldLocation::Elem,
                       3,
                       meshTargets );

  dofManager.addCoupling( contact::traction::key(),
                          contact::traction::key(),
                          DofManager::Connector::Face,
                          meshTargets );

  dofManager.addCoupling( solidMechanics::totalDisplacement::key(),
                          contact::traction::key(),
                          DofManager::Connector::Elem,
                          meshTargets );
}

void LagrangianContactSolver::assembleSystem( real64 const time,
                                              real64 const dt,
                                              DomainPartition & domain,
                                              DofManager const & dofManager,
                                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                              arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  synchronizeFractureState( domain );

  SolidMechanicsLagrangianFEM::assembleSystem( time,
                                               dt,
                                               domain,
                                               dofManager,
                                               localMatrix,
                                               localRhs );

  assembleContact( domain, dofManager, localMatrix, localRhs );
}

void LagrangianContactSolver::assembleContact( DomainPartition & domain,
                                               DofManager const & dofManager,
                                               CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                               arrayView1d< real64 > const & localRhs )
{
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    /// assemble Kut
    assembleForceResidualDerivativeWrtTraction( mesh, regionNames, dofManager, localMatrix, localRhs );
    /// assemble Ktu, Ktt blocks.
    assembleTractionResidualDerivativeWrtDisplacementAndTraction( mesh, regionNames, dofManager, localMatrix, localRhs );
    /// assemble stabilization
    assembleStabilization( mesh, domain.getNumericalMethodManager(), dofManager, localMatrix, localRhs );
  } );
}


real64 LagrangianContactSolver::calculateResidualNorm( real64 const & GEOS_UNUSED_PARAM( time ),
                                                       real64 const & GEOS_UNUSED_PARAM( dt ),
                                                       DomainPartition const & domain,
                                                       DofManager const & dofManager,
                                                       arrayView1d< real64 const > const & localRhs )
{
  GEOS_MARK_FUNCTION;
  real64 momentumR2 = 0.0;
  real64 contactR2 = 0.0;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel const & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    NodeManager const & nodeManager = mesh.getNodeManager();
    arrayView1d< globalIndex const > const & dispDofNumber =
      nodeManager.getReference< array1d< globalIndex > >( dofManager.getKey( solidMechanics::totalDisplacement::key() ) );

    string const & dofKey = dofManager.getKey( contact::traction::key() );
    globalIndex const rankOffset = dofManager.rankOffset();

    arrayView1d< integer const > const & elemGhostRank = nodeManager.ghostRank();

    RAJA::ReduceSum< parallelDeviceReduce, real64 > localSum0( 0.0 );
    forAll< parallelDevicePolicy<> >( nodeManager.size(),
                                      [localRhs, localSum0, dispDofNumber, rankOffset, elemGhostRank] GEOS_HOST_DEVICE ( localIndex const k )
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
    momentumR2 += localSum0.get();

    mesh.getElemManager().forElementSubRegions< FaceElementSubRegion >( regionNames,
                                                                        [&]( localIndex const, FaceElementSubRegion const & subRegion )
    {
      arrayView1d< globalIndex const > const & dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );
      arrayView1d< integer const > const & ghostRank = subRegion.ghostRank();

      RAJA::ReduceSum< parallelHostReduce, real64 > localSum( 0.0 );
      forAll< parallelHostPolicy >( subRegion.size(), [=] ( localIndex const k )
      {
        if( ghostRank[k] < 0 )
        {
          localIndex const localRow = LvArray::integerConversion< localIndex >( dofNumber[k] - rankOffset );
          for( localIndex dim = 0; dim < 3; ++dim )
          {
            localSum += localRhs[localRow + dim] * localRhs[localRow + dim];
          }
        }
      } );
      contactR2 += localSum.get();
    } );
  } );
  real64 localR2[2] = { momentumR2, contactR2 };
  real64 globalResidualNorm[3]{};

  int const rank = MpiWrapper::commRank( MPI_COMM_GEOSX );
  int const size = MpiWrapper::commSize( MPI_COMM_GEOSX );
  array1d< real64 > globalR2( 2 * size );
  globalR2.zero();

  // Everything is done on rank 0
  MpiWrapper::gather( localR2,
                      2,
                      globalR2.data(),
                      2,
                      0,
                      MPI_COMM_GEOSX );

  if( rank==0 )
  {
    globalResidualNorm[0] = 0.0;
    globalResidualNorm[1] = 0.0;
    for( int r=0; r<size; ++r )
    {
      // sum across all ranks
      globalResidualNorm[0] += globalR2[2 * r + 0];
      globalResidualNorm[1] += globalR2[2 * r + 1];
    }
    globalResidualNorm[2] = globalResidualNorm[0] + globalResidualNorm[1];
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
    globalResidualNorm[2] /= (m_initialResidual[2]+1.0);
  }
  if( getLogLevel() >= 1 && logger::internal::rank == 0 )
  {
    std::cout<< GEOS_FMT(
      "        ( Rdisplacement, Rtraction, Rtotal ) = ( {:15.6e}, {:15.6e}, {:15.6e} )",
      globalResidualNorm[0],
      globalResidualNorm[1],
      globalResidualNorm[2] );
  }
  return globalResidualNorm[2];
}

void LagrangianContactSolver::createPreconditioner( DomainPartition const & domain )
{
  if( m_linearSolverParameters.get().preconditionerType == LinearSolverParameters::PreconditionerType::block )
  {
    // TODO: move among inputs (xml)
    string const leadingBlockApproximation = "blockJacobi";

    LinearSolverParameters mechParams = getLinearSolverParameters();
    // Because of boundary conditions
    mechParams.isSymmetric = false;

    std::unique_ptr< BlockPreconditioner< LAInterface > > precond;
    std::unique_ptr< PreconditionerBase< LAInterface > > tracPrecond;

    if( leadingBlockApproximation == "jacobi" )
    {
      precond = std::make_unique< BlockPreconditioner< LAInterface > >( BlockShapeOption::LowerUpperTriangular,
                                                                        SchurComplementOption::FirstBlockDiagonal,
                                                                        BlockScalingOption::UserProvided );
      // Using GEOSX implementation of Jacobi preconditioner
      // tracPrecond = std::make_unique< PreconditionerJacobi< LAInterface > >();

      // Using LAI implementation of Jacobi preconditioner
      LinearSolverParameters tracParams;
      tracParams.preconditionerType = LinearSolverParameters::PreconditionerType::jacobi;
      tracPrecond = LAInterface::createPreconditioner( tracParams );
    }
    else if( leadingBlockApproximation == "blockJacobi" )
    {
      precond = std::make_unique< BlockPreconditioner< LAInterface > >( BlockShapeOption::LowerUpperTriangular,
                                                                        SchurComplementOption::FirstBlockUserDefined,
                                                                        BlockScalingOption::UserProvided );
      tracPrecond = std::make_unique< PreconditionerBlockJacobi< LAInterface > >( mechParams.dofsPerNode );
    }
    else
    {
      GEOS_ERROR( "LagrangianContactSolver::CreatePreconditioner leadingBlockApproximation option " << leadingBlockApproximation << " not supported" );
    }

    // Preconditioner for the leading block: tracPrecond
    precond->setupBlock( 0,
                         { { contact::traction::key(), { 3, true } } },
                         std::move( tracPrecond ) );

    if( mechParams.amg.nullSpaceType == LinearSolverParameters::AMG::NullSpaceType::rigidBodyModes )
    {
      if( getRigidBodyModes().empty() )
      {
        MeshLevel const & mesh = domain.getMeshBody( 0 ).getBaseDiscretization();
        LAIHelperFunctions::computeRigidBodyModes( mesh,
                                                   m_dofManager,
                                                   { solidMechanics::totalDisplacement::key() },
                                                   getRigidBodyModes() );
      }
    }

    // Preconditioner for the Schur complement: mechPrecond
    std::unique_ptr< PreconditionerBase< LAInterface > > mechPrecond = LAInterface::createPreconditioner( mechParams, getRigidBodyModes() );
    precond->setupBlock( 1,
                         { { solidMechanics::totalDisplacement::key(), { 3, true } } },
                         std::move( mechPrecond ) );

    m_precond = std::move( precond );
  }
  else
  {
    //TODO: Revisit this part such that is coherent across physics solver
    //m_precond = LAInterface::createPreconditioner( m_linearSolverParameters.get() );
  }
}

void LagrangianContactSolver::computeRotationMatrices( DomainPartition & domain ) const
{
  GEOS_MARK_FUNCTION;
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    FaceManager const & faceManager = mesh.getFaceManager();
    ElementRegionManager & elemManager = mesh.getElemManager();

    arrayView2d< real64 const > const & faceNormal = faceManager.faceNormal();

    elemManager.forElementSubRegions< FaceElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   FaceElementSubRegion & subRegion )
    {
      ArrayOfArraysView< localIndex const > const & elemsToFaces = subRegion.faceList().toViewConst();

      arrayView3d< real64 > const &
      rotationMatrix = subRegion.getReference< array3d< real64 > >( viewKeyStruct::rotationMatrixString() );

      forAll< parallelHostPolicy >( subRegion.size(), [=]( localIndex const kfe )
      {
        if( elemsToFaces.sizeOfArray( kfe ) != 2 )
        {
          return;
        }

        stackArray1d< real64, 3 > Nbar( 3 );
        localIndex const & f0 = elemsToFaces[kfe][0], f1 = elemsToFaces[kfe][1];
        Nbar[ 0 ] = faceNormal[f0][0] - faceNormal[f1][0];
        Nbar[ 1 ] = faceNormal[f0][1] - faceNormal[f1][1];
        Nbar[ 2 ] = faceNormal[f0][2] - faceNormal[f1][2];
        LvArray::tensorOps::normalize< 3 >( Nbar );

        computationalGeometry::RotationMatrix_3D( Nbar.toSliceConst(), rotationMatrix[kfe] );
      } );
    } );
  } );
}

void LagrangianContactSolver::computeFaceNodalArea( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition,
                                                    ArrayOfArraysView< localIndex const > const & faceToNodeMap,
                                                    localIndex const kf0,
                                                    array1d< real64 > & nodalArea ) const
{
  // I've tried to access the finiteElement::dispatch3D with
  // finiteElement::FiniteElementBase const &
  // fe = fractureSubRegion->getReference< finiteElement::FiniteElementBase >( surfaceGenerator->getDiscretizationName() );
  // but it's either empty (unknown discretization) or for 3D only (e.g., hexahedra)
  GEOS_MARK_FUNCTION;

  localIndex const TriangularPermutation[3] = { 0, 1, 2 };
  localIndex const QuadrilateralPermutation[4] = { 0, 1, 3, 2 };

  localIndex const numNodesPerFace = faceToNodeMap.sizeOfArray( kf0 );

  nodalArea.resize( numNodesPerFace );
  for( localIndex a = 0; a < numNodesPerFace; ++a )
  {
    nodalArea[a] = 0.0;
  }
  localIndex const * const permutation = ( numNodesPerFace == 3 ) ? TriangularPermutation : QuadrilateralPermutation;
  if( numNodesPerFace == 3 )
  {
    real64 xLocal[3][3];
    for( localIndex a = 0; a < numNodesPerFace; ++a )
    {
      for( localIndex j = 0; j < 3; ++j )
      {
        xLocal[a][j] = nodePosition[faceToNodeMap( kf0, permutation[a] )][j];
      }
    }
    real64 N[3];
    for( localIndex q=0; q<H1_TriangleFace_Lagrange1_Gauss1::numQuadraturePoints; ++q )
    {
      real64 const detJ = H1_TriangleFace_Lagrange1_Gauss1::transformedQuadratureWeight( q, xLocal );
      H1_TriangleFace_Lagrange1_Gauss1::calcN( q, N );
      for( localIndex a = 0; a < numNodesPerFace; ++a )
      {
        nodalArea[a] += detJ * N[permutation[a]];
      }
    }
  }
  else if( numNodesPerFace == 4 )
  {
    real64 xLocal[4][3];
    for( localIndex a = 0; a < numNodesPerFace; ++a )
    {
      for( localIndex j = 0; j < 3; ++j )
      {
        xLocal[a][j] = nodePosition[faceToNodeMap( kf0, permutation[a] )][j];
      }
    }
    real64 N[4];
    for( localIndex q=0; q<H1_QuadrilateralFace_Lagrange1_GaussLegendre2::numQuadraturePoints; ++q )
    {
      real64 const detJ = H1_QuadrilateralFace_Lagrange1_GaussLegendre2::transformedQuadratureWeight( q, xLocal );
      H1_QuadrilateralFace_Lagrange1_GaussLegendre2::calcN( q, N );
      for( localIndex a = 0; a < numNodesPerFace; ++a )
      {
        nodalArea[a] += detJ * N[permutation[a]];
      }
    }
  }
  else
  {
    GEOS_ERROR( "LagrangianContactSolver " << getDataContext() << ": face with " << numNodesPerFace <<
                " nodes. Only triangles and quadrilaterals are supported." );
  }
}

void LagrangianContactSolver::
  assembleForceResidualDerivativeWrtTraction( MeshLevel const & mesh,
                                              arrayView1d< string const > const & regionNames,
                                              DofManager const & dofManager,
                                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                              arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  FaceManager const & faceManager = mesh.getFaceManager();
  NodeManager const & nodeManager = mesh.getNodeManager();
  ElementRegionManager const & elemManager = mesh.getElemManager();

  ArrayOfArraysView< localIndex const > const faceToNodeMap = faceManager.nodeList().toViewConst();

  string const & tracDofKey = dofManager.getKey( contact::traction::key() );
  string const & dispDofKey = dofManager.getKey( solidMechanics::totalDisplacement::key() );

  arrayView1d< globalIndex const > const & dispDofNumber = nodeManager.getReference< globalIndex_array >( dispDofKey );
  globalIndex const rankOffset = dofManager.rankOffset();

  // Get the coordinates for all nodes
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition = nodeManager.referencePosition();

  elemManager.forElementSubRegions< FaceElementSubRegion >( regionNames,
                                                            [&]( localIndex const,
                                                                 FaceElementSubRegion const & subRegion )
  {
    arrayView1d< globalIndex const > const & tracDofNumber = subRegion.getReference< globalIndex_array >( tracDofKey );
    arrayView2d< real64 const > const & traction = subRegion.getReference< array2d< real64 > >( contact::traction::key() );
    arrayView3d< real64 const > const & rotationMatrix = subRegion.getReference< array3d< real64 > >( viewKeyStruct::rotationMatrixString() );
    ArrayOfArraysView< localIndex const > const & elemsToFaces = subRegion.faceList().toViewConst();

    constexpr localIndex TriangularPermutation[3] = { 0, 1, 2 };
    constexpr localIndex QuadrilateralPermutation[4] = { 0, 1, 3, 2 };

    forAll< parallelHostPolicy >( subRegion.size(), [=] ( localIndex const kfe )
    {
      if( elemsToFaces.sizeOfArray( kfe ) != 2 )
      {
        return;
      }

      localIndex const numNodesPerFace = faceToNodeMap.sizeOfArray( elemsToFaces[kfe][0] );
      localIndex const numQuadraturePointsPerElem = numNodesPerFace==3 ? 1 : 4;

      globalIndex rowDOF[12];
      real64 nodeRHS[12];
      stackArray2d< real64, 3*4*3 > dRdT( 3*numNodesPerFace, 3 );
      globalIndex colDOF[3];
      for( localIndex i = 0; i < 3; ++i )
      {
        colDOF[i] = tracDofNumber[kfe] + i;
      }

      localIndex const * const permutation = ( numNodesPerFace == 3 ) ? TriangularPermutation : QuadrilateralPermutation;
      real64 xLocal[2][4][3];
      for( localIndex kf = 0; kf < 2; ++kf )
      {
        localIndex const faceIndex = elemsToFaces[kfe][kf];
        for( localIndex a = 0; a < numNodesPerFace; ++a )
        {
          for( localIndex j = 0; j < 3; ++j )
          {
            xLocal[kf][a][j] = nodePosition[ faceToNodeMap( faceIndex, permutation[a] ) ][j];
          }
        }
      }

      real64 N[4];

      for( localIndex q=0; q<numQuadraturePointsPerElem; ++q )
      {
        if( numNodesPerFace==3 )
        {
          using NT = real64[3];
          H1_TriangleFace_Lagrange1_Gauss1::calcN( q, reinterpret_cast< NT & >(N) );
        }
        else if( numNodesPerFace==4 )
        {
          H1_QuadrilateralFace_Lagrange1_GaussLegendre2::calcN( q, N );
        }

        constexpr int normalSign[2] = { 1, -1 };
        for( localIndex kf = 0; kf < 2; ++kf )
        {
          localIndex const faceIndex = elemsToFaces[kfe][kf];
          using xLocalTriangle = real64[3][3];
          real64 const detJxW = numNodesPerFace==3 ?
                                H1_TriangleFace_Lagrange1_Gauss1::transformedQuadratureWeight( q, reinterpret_cast< xLocalTriangle & >( xLocal[kf] ) ) :
                                H1_QuadrilateralFace_Lagrange1_GaussLegendre2::transformedQuadratureWeight( q, xLocal[kf] );

          for( localIndex a = 0; a < numNodesPerFace; ++a )
          {
            real64 const NaDetJxQ = N[permutation[a]] * detJxW;
            real64 const localNodalForce[ 3 ] = { traction( kfe, 0 ) * NaDetJxQ,
                                                  traction( kfe, 1 ) * NaDetJxQ,
                                                  traction( kfe, 2 ) * NaDetJxQ };
            real64 globalNodalForce[ 3 ];
            LvArray::tensorOps::Ri_eq_AijBj< 3, 3 >( globalNodalForce, rotationMatrix[ kfe ], localNodalForce );

            for( localIndex i = 0; i < 3; ++i )
            {
              rowDOF[3*a+i] = dispDofNumber[faceToNodeMap( faceIndex, a )] + i;
              // Opposite sign w.r.t. to formulation presented in
              // Algebraically Stabilized Lagrange Multiplier Method for Frictional Contact Mechanics with
              // Hydraulically Active Fractures
              // Franceschini, A., Castelletto, N., White, J. A., Tchelepi, H. A.
              // Computer Methods in Applied Mechanics and Engineering (2020) 368, 113161
              // doi: 10.1016/j.cma.2020.113161
              nodeRHS[3*a+i] = +globalNodalForce[i] * normalSign[ kf ];

              // Opposite sign w.r.t. to the same formulation as above
              dRdT( 3*a+i, 0 ) = rotationMatrix( kfe, i, 0 ) * normalSign[ kf ] * NaDetJxQ;
              dRdT( 3*a+i, 1 ) = rotationMatrix( kfe, i, 1 ) * normalSign[ kf ] * NaDetJxQ;
              dRdT( 3*a+i, 2 ) = rotationMatrix( kfe, i, 2 ) * normalSign[ kf ] * NaDetJxQ;
            }
          }

          for( localIndex idof = 0; idof < numNodesPerFace * 3; ++idof )
          {
            localIndex const localRow = LvArray::integerConversion< localIndex >( rowDOF[idof] - rankOffset );

            if( localRow >= 0 && localRow < localMatrix.numRows() )
            {
              localMatrix.addToRow< parallelHostAtomic >( localRow,
                                                          colDOF,
                                                          dRdT[idof].dataIfContiguous(),
                                                          3 );
              RAJA::atomicAdd( parallelHostAtomic{}, &localRhs[localRow], nodeRHS[idof] );
            }
          }
        }
      }
    } );
  } );
}

void LagrangianContactSolver::
  assembleTractionResidualDerivativeWrtDisplacementAndTraction( MeshLevel const & mesh,
                                                                arrayView1d< string const > const & regionNames,
                                                                DofManager const & dofManager,
                                                                CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                                arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;
  FaceManager const & faceManager = mesh.getFaceManager();
  NodeManager const & nodeManager = mesh.getNodeManager();
  ElementRegionManager const & elemManager = mesh.getElemManager();

  ArrayOfArraysView< localIndex const > const faceToNodeMap = faceManager.nodeList().toViewConst();

  string const & tracDofKey = dofManager.getKey( contact::traction::key() );
  string const & dispDofKey = dofManager.getKey( solidMechanics::totalDisplacement::key() );

  arrayView1d< globalIndex const > const & dispDofNumber = nodeManager.getReference< globalIndex_array >( dispDofKey );
  globalIndex const rankOffset = dofManager.rankOffset();

  // Get the coordinates for all nodes
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition = nodeManager.referencePosition();

  elemManager.forElementSubRegions< FaceElementSubRegion >( regionNames,
                                                            [&]( localIndex const,
                                                                 FaceElementSubRegion const & subRegion )
  {
    string const & contactRelationName = subRegion.template getReference< string >( viewKeyStruct::contactRelationNameString() );
    ContactBase const & contact = getConstitutiveModel< ContactBase >( subRegion, contactRelationName );

    arrayView1d< globalIndex const > const & tracDofNumber = subRegion.getReference< globalIndex_array >( tracDofKey );
    arrayView1d< integer const > const & ghostRank = subRegion.ghostRank();
    arrayView1d< real64 const > const & area = subRegion.getElementArea();
    arrayView3d< real64 const > const &
    rotationMatrix = subRegion.getReference< array3d< real64 > >( viewKeyStruct::rotationMatrixString() );
    ArrayOfArraysView< localIndex const > const & elemsToFaces = subRegion.faceList().toViewConst();
    arrayView2d< real64 const > const & traction = subRegion.getField< contact::traction >();
    arrayView1d< integer const > const & fractureState = subRegion.getField< contact::fractureState >();
    arrayView2d< real64 const > const & dispJump = subRegion.getField< contact::dispJump >();
    arrayView2d< real64 const > const & previousDispJump = subRegion.getField< contact::oldDispJump >();
    arrayView1d< real64 const > const & slidingTolerance = subRegion.getReference< array1d< real64 > >( viewKeyStruct::slidingToleranceString() );

    constitutiveUpdatePassThru( contact, [&] ( auto & castedContact )
    {
      using ContactType = TYPEOFREF( castedContact );
      typename ContactType::KernelWrapper contactWrapper = castedContact.createKernelWrapper();

      forAll< parallelHostPolicy >( subRegion.size(), [=] ( localIndex const kfe )
      {
        if( elemsToFaces.sizeOfArray( kfe ) != 2 )
        {
          return;
        }

        if( ghostRank[kfe] < 0 )
        {
          localIndex const numNodesPerFace = faceToNodeMap.sizeOfArray( elemsToFaces[kfe][0] );
          globalIndex nodeDOF[24];
          globalIndex elemDOF[3];
          for( localIndex i = 0; i < 3; ++i )
          {
            elemDOF[i] = tracDofNumber[kfe] + i;
          }

          real64 elemRHS[3] = {0.0, 0.0, 0.0};
          real64 const Ja = area[kfe];

          stackArray2d< real64, 2 * 3 * 4 * 3 > dRdU( 3, 2 * 3 * numNodesPerFace );
          stackArray2d< real64, 3 * 3 > dRdT( 3, 3 );

          switch( fractureState[kfe] )
          {
            case contact::FractureState::Stick:
              {
                for( localIndex i = 0; i < 3; ++i )
                {
                  if( i == 0 )
                  {
                    elemRHS[i] = +Ja * dispJump[kfe][i];
                  }
                  else
                  {
                    elemRHS[i] = +Ja * ( dispJump[kfe][i] - previousDispJump[kfe][i] );
                  }
                }

                for( localIndex kf = 0; kf < 2; ++kf )
                {
                  // Compute local area contribution for each node
                  array1d< real64 > nodalArea;
                  computeFaceNodalArea( nodePosition, faceToNodeMap, elemsToFaces[kfe][kf], nodalArea );

                  for( localIndex a = 0; a < numNodesPerFace; ++a )
                  {
                    for( localIndex i = 0; i < 3; ++i )
                    {
                      nodeDOF[kf * 3 * numNodesPerFace + 3 * a + i] = dispDofNumber[faceToNodeMap( elemsToFaces[kfe][kf], a )] + i;

                      for( localIndex j = 0; j < 3; ++j )
                      {
                        dRdU( j, kf * 3 * numNodesPerFace + 3 * a + i ) = -nodalArea[a] * rotationMatrix( kfe, i, j ) * pow( -1, kf );
                      }
                    }
                  }
                }
                break;
              }
            case contact::FractureState::Slip:
            case contact::FractureState::NewSlip:
              {
                elemRHS[0] = +Ja * dispJump[kfe][0];

                for( localIndex kf = 0; kf < 2; ++kf )
                {
                  // Compute local area contribution for each node
                  array1d< real64 > nodalArea;
                  computeFaceNodalArea( nodePosition, faceToNodeMap, elemsToFaces[kfe][kf], nodalArea );

                  for( localIndex a = 0; a < numNodesPerFace; ++a )
                  {
                    for( localIndex i = 0; i < 3; ++i )
                    {
                      nodeDOF[kf * 3 * numNodesPerFace + 3 * a + i] = dispDofNumber[faceToNodeMap( elemsToFaces[kfe][kf], a )] +
                                                                      LvArray::integerConversion< globalIndex >( i );
                      dRdU( 0, kf * 3 * numNodesPerFace + 3 * a + i ) = -nodalArea[a] * rotationMatrix( kfe, i, 0 ) * pow( -1, kf );
                    }
                  }
                }

                real64 dLimitTau_dNormalTraction = 0;
                real64 const limitTau = contactWrapper.computeLimitTangentialTractionNorm( traction[kfe][0],
                                                                                           dLimitTau_dNormalTraction );

                real64 sliding[ 2 ] = { dispJump[kfe][1] - previousDispJump[kfe][1], dispJump[kfe][2] - previousDispJump[kfe][2] };
                real64 slidingNorm = sqrt( sliding[ 0 ]*sliding[ 0 ] + sliding[ 1 ]*sliding[ 1 ] );

                if( !( ( m_nonlinearSolverParameters.m_numNewtonIterations == 0 ) && ( fractureState[kfe] == contact::FractureState::NewSlip ) )
                    && slidingNorm > slidingTolerance[kfe] )
                {
                  for( localIndex i = 1; i < 3; ++i )
                  {
                    elemRHS[i] = +Ja * ( traction[kfe][i] - limitTau * sliding[ i-1 ] / slidingNorm );
                  }

                  // A symmetric 2x2 matrix.
                  real64 dUdgT[ 3 ];
                  dUdgT[ 0 ] = (slidingNorm * slidingNorm - sliding[ 0 ] * sliding[ 0 ]) * limitTau / std::pow( slidingNorm, 3 );
                  dUdgT[ 1 ] = (slidingNorm * slidingNorm - sliding[ 1 ] * sliding[ 1 ]) * limitTau / std::pow( slidingNorm, 3 );
                  dUdgT[ 2 ] = -sliding[ 0 ] * sliding[ 1 ] * limitTau / std::pow( slidingNorm, 3 );

                  for( localIndex kf = 0; kf < 2; ++kf )
                  {
                    // Compute local area contribution for each node
                    array1d< real64 > nodalArea;
                    computeFaceNodalArea( nodePosition, faceToNodeMap, elemsToFaces[kfe][kf], nodalArea );

                    for( localIndex a = 0; a < numNodesPerFace; ++a )
                    {
                      for( localIndex i = 0; i < 3; ++i )
                      {
                        real64 const localRowB[ 2 ] = { rotationMatrix( kfe, i, 1 ), rotationMatrix( kfe, i, 2 ) };
                        real64 localRowE[ 2 ];
                        LvArray::tensorOps::Ri_eq_symAijBj< 2 >( localRowE, dUdgT, localRowB );

                        dRdU( 1, kf * 3 * numNodesPerFace + 3 * a + i ) = nodalArea[a] * localRowE[ 0 ] * pow( -1, kf );
                        dRdU( 2, kf * 3 * numNodesPerFace + 3 * a + i ) = nodalArea[a] * localRowE[ 1 ] * pow( -1, kf );
                      }
                    }
                  }
                  for( localIndex i = 1; i < 3; ++i )
                  {
                    dRdT( i, 0 ) = Ja * dLimitTau_dNormalTraction * sliding[ i-1 ] / slidingNorm;
                    dRdT( i, i ) = Ja;
                  }
                }
                else
                {
                  real64 vaux[ 2 ] = { traction[kfe][1], traction[kfe][2] };
                  real64 vauxNorm = sqrt( vaux[ 0 ]*vaux[ 0 ] + vaux[ 1 ]*vaux[ 1 ] );
                  if( vauxNorm > 0.0 )
                  {
                    for( localIndex i = 1; i < 3; ++i )
                    {
                      elemRHS[i] = +Ja * ( traction[kfe][i] - limitTau * vaux[ i-1 ] / vauxNorm );
                    }
                    for( localIndex i = 1; i < 3; ++i )
                    {
                      dRdT( i, i ) = Ja;
                    }
                  }
                  else
                  {
                    for( localIndex i = 1; i < 3; ++i )
                    {
                      elemRHS[i] = 0.0;
                    }
                    for( localIndex i = 1; i < 3; ++i )
                    {
                      dRdT( i, i ) = Ja;
                    }
                  }
                }
                break;
              }
            case contact::FractureState::Open:
              {
                for( localIndex i = 0; i < 3; ++i )
                {
                  elemRHS[i] = +Ja * traction[kfe][i];
                }

                for( localIndex i = 0; i < 3; ++i )
                {
                  dRdT( i, i ) = Ja;
                }
                break;
              }
          }

          localIndex const localRow = LvArray::integerConversion< localIndex >( elemDOF[0] - rankOffset );


          for( localIndex idof = 0; idof < 3; ++idof )
          {
            localRhs[localRow + idof] += elemRHS[idof];

            if( fractureState[kfe] != contact::FractureState::Open )
            {
              localMatrix.addToRowBinarySearchUnsorted< serialAtomic >( localRow + idof,
                                                                        nodeDOF,
                                                                        dRdU[idof].dataIfContiguous(),
                                                                        2 * 3 * numNodesPerFace );
            }

            if( fractureState[kfe] != contact::FractureState::Stick )
            {
              localMatrix.addToRow< serialAtomic >( localRow + idof,
                                                    elemDOF,
                                                    dRdT[idof].dataIfContiguous(),
                                                    3 );
            }
          }
        }
      } );
    } );
  } );
}

void LagrangianContactSolver::assembleStabilization( MeshLevel const & mesh,
                                                     NumericalMethodsManager const & numericalMethodManager,
                                                     DofManager const & dofManager,
                                                     CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                     arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  FaceManager const & faceManager = mesh.getFaceManager();
  NodeManager const & nodeManager = mesh.getNodeManager();
  ElementRegionManager const & elemManager = mesh.getElemManager();

  string const & tracDofKey = dofManager.getKey( contact::traction::key() );
  globalIndex const rankOffset = dofManager.rankOffset();

  // Get the finite volume method used to compute the stabilization
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & stabilizationMethod = fvManager.getFluxApproximation( m_stabilizationName );

  // Get the "face to element" map (valid for the entire mesh)
  FaceManager::ElemMapType const & faceToElem = faceManager.toElementRelation();
  arrayView2d< localIndex const > const & faceToElemRegion = faceToElem.m_toElementRegion.toViewConst();
  arrayView2d< localIndex const > const & faceToElemSubRegion = faceToElem.m_toElementSubRegion.toViewConst();
  arrayView2d< localIndex const > const & faceToElemIndex = faceToElem.m_toElementIndex.toViewConst();

  SurfaceElementRegion const & fractureRegion = elemManager.getRegion< SurfaceElementRegion >( getFractureRegionName() );
  FaceElementSubRegion const & fractureSubRegion = fractureRegion.getUniqueSubRegion< FaceElementSubRegion >();

  GEOS_ERROR_IF( !fractureSubRegion.hasField< contact::traction >(),
                 getDataContext() << ": The fracture subregion must contain traction field." );
  ArrayOfArraysView< localIndex const > const elem2dToFaces = fractureSubRegion.faceList().toViewConst();

  // Get the state of fracture elements
  arrayView1d< integer const > const & fractureState =
    fractureSubRegion.getReference< array1d< integer > >( viewKeyStruct::fractureStateString() );

  // Get the tractions and stabilization contribution to the local jump
  arrayView2d< real64 const > const & traction = fractureSubRegion.getField< contact::traction >();
  arrayView2d< real64 const > const & deltaTraction = fractureSubRegion.getField< contact::deltaTraction >();

  // Get the volume for all elements
  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > const elemVolume =
    elemManager.constructViewAccessor< real64_array, arrayView1d< real64 const > >( ElementSubRegionBase::viewKeyStruct::elementVolumeString() );

  // Get the coordinates for all nodes
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition = nodeManager.referencePosition();

  // Get area and rotation matrix for all faces
  ArrayOfArraysView< localIndex const > const & faceToNodeMap = faceManager.nodeList().toViewConst();
  arrayView1d< real64 const > const & faceArea = faceManager.faceArea();
  arrayView3d< real64 const > const &
  faceRotationMatrix = fractureSubRegion.getReference< array3d< real64 > >( viewKeyStruct::rotationMatrixString() );

  // Bulk modulus accessor
  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > const bulkModulus =
    elemManager.constructMaterialViewAccessor< ElasticIsotropic, array1d< real64 >, arrayView1d< real64 const > >( ElasticIsotropic::viewKeyStruct::bulkModulusString() );
  // Shear modulus accessor
  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > const shearModulus =
    elemManager.constructMaterialViewAccessor< ElasticIsotropic, array1d< real64 >, arrayView1d< real64 const > >( ElasticIsotropic::viewKeyStruct::shearModulusString() );

  using NodeMapViewType = arrayView2d< localIndex const, cells::NODE_MAP_USD >;
  ElementRegionManager::ElementViewAccessor< NodeMapViewType > const elemToNode =
    elemManager.constructViewAccessor< CellElementSubRegion::NodeMapType, NodeMapViewType >( ElementSubRegionBase::viewKeyStruct::nodeListString() );
  ElementRegionManager::ElementViewConst< NodeMapViewType > const elemToNodeView = elemToNode.toNestedViewConst();

  arrayView1d< globalIndex const > const & tracDofNumber = fractureSubRegion.getReference< globalIndex_array >( tracDofKey );

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
        Nbar0[ 0 ] = faceRotationMatrix[ sei[iconn][0] ][0][0];
        Nbar0[ 1 ] = faceRotationMatrix[ sei[iconn][0] ][1][0];
        Nbar0[ 2 ] = faceRotationMatrix[ sei[iconn][0] ][2][0];
        Nbar1[ 0 ] = faceRotationMatrix[ sei[iconn][1] ][0][0];
        Nbar1[ 1 ] = faceRotationMatrix[ sei[iconn][1] ][1][0];
        Nbar1[ 2 ] = faceRotationMatrix[ sei[iconn][1] ][2][0];

        real64 normalProduct = LvArray::tensorOps::AiBi< 3 >( Nbar0, Nbar1 );

        localIndex const id1 = ( normalProduct > 0.0 ) ? 0 : 1;

        localIndex const numNodesPerFace0 = faceToNodeMap.sizeOfArray( elem2dToFaces[sei[iconn][0]][0] );
        array1d< localIndex > nodes0( numNodesPerFace0 );
        for( localIndex i = 0; i < numNodesPerFace0; ++i )
        {
          nodes0[i] = faceToNodeMap( elem2dToFaces[sei[iconn][0]][0], i );
        }
        localIndex const numNodesPerFace1 = faceToNodeMap.sizeOfArray( elem2dToFaces[sei[iconn][1]][0] );
        array1d< localIndex > nodes1( numNodesPerFace1 );
        for( localIndex i = 0; i < numNodesPerFace1; ++i )
        {
          nodes1[i] = faceToNodeMap( elem2dToFaces[sei[iconn][1]][id1], i );
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
        GEOS_ERROR_IF( realNodes != 2,
                       getDataContext() << ": An edge shared by two fracture elements must have 2 nodes." );
        edge.resize( realNodes );

        // Compute nodal area factor
        localIndex node0index0 = -1;
        localIndex node1index0 = -1;
        for( localIndex i = 0; i < numNodesPerFace0; ++i )
        {
          if( edge[0] == faceToNodeMap( elem2dToFaces[sei[iconn][0]][0], i ) )
          {
            node0index0 = i;
          }
          if( edge[1] == faceToNodeMap( elem2dToFaces[sei[iconn][0]][0], i ) )
          {
            node1index0 = i;
          }
        }
        localIndex node0index1 = -1;
        localIndex node1index1 = -1;
        for( localIndex i = 0; i < numNodesPerFace1; ++i )
        {
          if( edge[0] == faceToNodeMap( elem2dToFaces[sei[iconn][1]][id1], i ) )
          {
            node0index1 = i;
          }
          if( edge[1] == faceToNodeMap( elem2dToFaces[sei[iconn][1]][id1], i ) )
          {
            node1index1 = i;
          }
        }
        array1d< real64 > nodalArea0, nodalArea1;
        computeFaceNodalArea( nodePosition, faceToNodeMap, elem2dToFaces[sei[iconn][0]][0], nodalArea0 );
        computeFaceNodalArea( nodePosition, faceToNodeMap, elem2dToFaces[sei[iconn][1]][id1], nodalArea1 );
        real64 const areafac = nodalArea0[node0index0] * nodalArea1[node0index1] + nodalArea0[node1index0] * nodalArea1[node1index1];

        // first index: face, second index: element (T/B), third index: dof (x, y, z)
        real64 stiffDiagApprox[ 2 ][ 2 ][ 3 ];
        for( localIndex kf = 0; kf < 2; ++kf )
        {
          // Get fracture, face and region/subregion/element indices (for elements on both sides)
          localIndex const fractureIndex = sei[iconn][kf];

          for( localIndex i = 0; i < 2; ++i )
          {
            localIndex const faceIndex = ( kf == 0 || id1 == 0 ) ? elem2dToFaces[fractureIndex][i] : elem2dToFaces[fractureIndex][1 - i];
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

            // Combine E and nu to obtain a stiffness approximation (like it was an hexahedron)
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

        array2d< real64 > avgRotationMatrix( 3, 3 );

        // To be able to compute an average rotation matrix, normal has to point in the same direction.
        if( normalProduct < 0.0 )
        {
          LvArray::tensorOps::scale< 3 >( Nbar1, -1.0 );
          normalProduct *= -1.0;
        }
        // If the surfaces are co-planar, then use the first rotation matrix
        if( std::abs( normalProduct - 1.0 ) < 1.e+2*machinePrecision )
        {
          LvArray::tensorOps::copy< 3, 3 >( avgRotationMatrix, faceRotationMatrix[ sei[iconn][0] ] );
        }
        // otherwise, compute the average rotation matrix
        else
        {
          array1d< real64 > avgNbar( 3 );
          avgNbar[ 0 ] = faceArea[elem2dToFaces[ sei[iconn][0] ][0]] * Nbar0[0] + faceArea[elem2dToFaces[ sei[iconn][1] ][0]] * Nbar1[0];
          avgNbar[ 1 ] = faceArea[elem2dToFaces[ sei[iconn][0] ][0]] * Nbar0[1] + faceArea[elem2dToFaces[ sei[iconn][1] ][0]] * Nbar1[1];
          avgNbar[ 2 ] = faceArea[elem2dToFaces[ sei[iconn][0] ][0]] * Nbar0[2] + faceArea[elem2dToFaces[ sei[iconn][1] ][0]] * Nbar1[2];
          LvArray::tensorOps::normalize< 3 >( avgNbar );

          computationalGeometry::RotationMatrix_3D( avgNbar.toSliceConst(), avgRotationMatrix );
        }

        // Compute R^T * (invK) * R
        real64 temp[ 3 ][ 3 ];
        real64 rotatedInvStiffApprox[ 3 ][ 3 ];
        LvArray::tensorOps::Rij_eq_AkiBkj< 3, 3, 3 >( temp, avgRotationMatrix, invTotStiffApprox );
        LvArray::tensorOps::Rij_eq_AikBkj< 3, 3, 3 >( rotatedInvStiffApprox, temp, avgRotationMatrix );

        // Add nodal area contribution
        stackArray2d< real64, 3*3 > totalInvStiffApprox( 3, 3 );
        for( localIndex i = 0; i < 3; ++i )
        {
          for( localIndex j = 0; j < 3; ++j )
          {
            totalInvStiffApprox( i, j ) = -rotatedInvStiffApprox[ i ][ j ] * areafac;
          }
        }

        // Get DOF numbering
        localIndex fractureIndex[2];
        localIndex nDof[2];
        globalIndex elemDOF[2][3];
        for( localIndex kf = 0; kf < 2; ++kf )
        {
          fractureIndex[kf] = sei[iconn][kf];
          for( localIndex i = 0; i < 3; ++i )
          {
            elemDOF[kf][i] = tracDofNumber[fractureIndex[kf]] + i;
          }
          nDof[kf] = 0;
          switch( fractureState[fractureIndex[kf]] )
          {
            case ( contact::FractureState::Stick ):
              {
                nDof[kf] = 3;
                break;
              }
            case ( contact::FractureState::NewSlip ):
            case ( contact::FractureState::Slip ):
              {
                nDof[kf] = 1;
                break;
              }
            case ( contact::FractureState::Open ):
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
        stackArray1d< real64, 3 > rhs0( 3 );
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

        stackArray1d< real64, 3 > rhs1( 3 );
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
        for( localIndex kf = 0; kf < 2; ++kf )
        {
          localIndex const localRow = LvArray::integerConversion< localIndex >( elemDOF[kf][0] - rankOffset );

          stackArray2d< real64, 3*3 > const & totalInvStiffApproxDiag = ( kf == 0 ) ? totalInvStiffApprox00 : totalInvStiffApprox11;
          stackArray2d< real64, 3*3 > const & totalInvStiffApproxOffDiag = ( kf == 0 ) ? totalInvStiffApprox01 : totalInvStiffApprox10;
          stackArray1d< real64, 3 > const & rhs = ( kf == 0 ) ? rhs0 : rhs1;

          // Only assemble contribution if "row" fracture element is local
          // TODO: use parallel atomics
          if( localRow >= 0 && localRow < localMatrix.numRows() )
          {
            for( localIndex idof = 0; idof < nDof[kf]; ++idof )
            {
              // (i,i)-block
              localMatrix.addToRowBinarySearchUnsorted< parallelHostAtomic >( localRow + idof,
                                                                              elemDOF[kf],
                                                                              totalInvStiffApproxDiag[idof].dataIfContiguous(),
                                                                              nDof[kf] );
              // (i,j)-block
              if( nDof[1-kf] > 0 )
              {
                localMatrix.addToRowBinarySearchUnsorted< parallelHostAtomic >( localRow + idof,
                                                                                elemDOF[1 - kf],
                                                                                totalInvStiffApproxOffDiag[idof].dataIfContiguous(),
                                                                                nDof[1 - kf] );
              }

              // residual
              RAJA::atomicAdd( parallelHostAtomic{}, &localRhs[localRow + idof], rhs[idof] );
            }
          }
        }
      }
    } );
  } );
}

void LagrangianContactSolver::applySystemSolution( DofManager const & dofManager,
                                                   arrayView1d< real64 const > const & localSolution,
                                                   real64 const scalingFactor,
                                                   real64 const dt,
                                                   DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  //if( m_setupSolidSolverDofs )
  {
    SolidMechanicsLagrangianFEM::applySystemSolution( dofManager, localSolution, scalingFactor, dt, domain );
  }

  dofManager.addVectorToField( localSolution,
                               contact::traction::key(),
                               contact::deltaTraction::key(),
                               scalingFactor );

  dofManager.addVectorToField( localSolution,
                               contact::traction::key(),
                               contact::traction::key(),
                               scalingFactor );

  // fractureStateString is synchronized in UpdateFractureState
  // oldFractureStateString and oldDispJumpString used locally only

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & )
  {
    FieldIdentifiers fieldsToBeSync;

    fieldsToBeSync.addElementFields( { contact::traction::key(),
                                       contact::deltaTraction::key(),
                                       contact::dispJump::key() },
                                     { getFractureRegionName() } );

    CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync,
                                                         mesh,
                                                         domain.getNeighbors(),
                                                         true );
  } );
}

void LagrangianContactSolver::updateState( DomainPartition & domain )
{
  computeFaceDisplacementJump( domain );
}

bool LagrangianContactSolver::resetConfigurationToDefault( DomainPartition & domain ) const
{
  GEOS_MARK_FUNCTION;

  using namespace fields::contact;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< FaceElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                FaceElementSubRegion & subRegion )
    {
      if( subRegion.hasField< contact::traction >() )
      {
        arrayView1d< integer > const & fractureState = subRegion.getField< contact::fractureState >();
        forAll< parallelHostPolicy >( subRegion.size(), [=] ( localIndex const kfe )
        {
          if( fractureState[kfe] != FractureState::Open )
          {
            fractureState[kfe] = FractureState::Stick;
          }
        } );
      }
    } );
  } );
  return false;
}

bool LagrangianContactSolver::updateConfiguration( DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  using namespace fields::contact;

  int hasConfigurationConverged = true;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< FaceElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                FaceElementSubRegion & subRegion )
    {
      string const & contactRelationName = subRegion.template getReference< string >( viewKeyStruct::contactRelationNameString() );
      ContactBase const & contact = getConstitutiveModel< ContactBase >( subRegion, contactRelationName );

      arrayView1d< integer const > const & ghostRank = subRegion.ghostRank();
      arrayView2d< real64 const > const & traction = subRegion.getField< contact::traction >();
      arrayView2d< real64 const > const & dispJump = subRegion.getField< contact::dispJump >();
      arrayView1d< integer > const & fractureState = subRegion.getField< contact::fractureState >();

      arrayView1d< real64 const > const & normalTractionTolerance =
        subRegion.getReference< array1d< real64 > >( viewKeyStruct::normalTractionToleranceString() );
      arrayView1d< real64 const > const & normalDisplacementTolerance =
        subRegion.getReference< array1d< real64 > >( viewKeyStruct::normalDisplacementToleranceString() );

      RAJA::ReduceMin< parallelHostReduce, integer > checkActiveSetSub( 1 );

      constitutiveUpdatePassThru( contact, [&] ( auto & castedContact )
      {
        using ContactType = TYPEOFREF( castedContact );
        typename ContactType::KernelWrapper contactWrapper = castedContact.createKernelWrapper();

        forAll< parallelHostPolicy >( subRegion.size(), [=] ( localIndex const kfe )
        {
          if( ghostRank[kfe] < 0 )
          {
            integer const originalFractureState = fractureState[kfe];
            if( originalFractureState == contact::FractureState::Open )
            {
              if( dispJump[kfe][0] > -normalDisplacementTolerance[kfe] )
              {
                fractureState[kfe] = contact::FractureState::Open;
              }
              else
              {
                fractureState[kfe] = contact::FractureState::Stick;
              }
            }
            else if( traction[kfe][0] > normalTractionTolerance[kfe] )
            {
              fractureState[kfe] = contact::FractureState::Open;
            }
            else
            {
              real64 currentTau = sqrt( traction[kfe][1]*traction[kfe][1] + traction[kfe][2]*traction[kfe][2] );

              real64 dLimitTangentialTractionNorm_dTraction = 0.0;
              real64 const limitTau =
                contactWrapper.computeLimitTangentialTractionNorm( traction[kfe][0],
                                                                   dLimitTangentialTractionNorm_dTraction );

              if( originalFractureState == contact::FractureState::Stick && currentTau >= limitTau )
              {
                currentTau *= (1.0 - m_slidingCheckTolerance);
              }
              else if( originalFractureState != contact::FractureState::Stick && currentTau <= limitTau )
              {
                currentTau *= (1.0 + m_slidingCheckTolerance);
              }
              if( currentTau > limitTau )
              {
                if( originalFractureState == contact::FractureState::Stick )
                {
                  fractureState[kfe] = contact::FractureState::NewSlip;
                }
                else
                {
                  fractureState[kfe] = contact::FractureState::Slip;
                }
              }
              else
              {
                fractureState[kfe] = contact::FractureState::Stick;
              }
            }
            checkActiveSetSub.min( compareFractureStates( originalFractureState, fractureState[kfe] ) );
          }
        } );
      } );

      hasConfigurationConverged &= checkActiveSetSub.get();
    } );
  } );
  // Need to synchronize the fracture state due to the use will be made of in AssemblyStabilization
  synchronizeFractureState( domain );

  // Compute if globally the fracture state has changed
  int hasConfigurationConvergedGlobally;
  MpiWrapper::allReduce( &hasConfigurationConverged,
                         &hasConfigurationConvergedGlobally,
                         1,
                         MPI_LAND,
                         MPI_COMM_GEOSX );

  return hasConfigurationConvergedGlobally;
}

bool LagrangianContactSolver::isFractureAllInStickCondition( DomainPartition const & domain ) const
{
  globalIndex numStick, numSlip, numOpen;
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel const & mesh,
                                                                arrayView1d< string const > const & )
  {
    computeFractureStateStatistics( mesh, numStick, numSlip, numOpen );
  } );

  return ( ( numSlip + numOpen ) == 0 );
}

real64 LagrangianContactSolver::setNextDt( real64 const & currentDt,
                                           DomainPartition & domain )
{
  GEOS_UNUSED_VAR( domain );
  return currentDt;
}

REGISTER_CATALOG_ENTRY( SolverBase, LagrangianContactSolver, string const &, Group * const )

} /* namespace geos */
