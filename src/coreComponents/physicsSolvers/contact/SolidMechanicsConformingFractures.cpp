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
 * @file SolidMechanicsConformingFractures.cpp
 *
 */

#include "SolidMechanicsConformingFractures.hpp"

#include "common/TimingMacros.hpp"

#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/contact/ContactSelector.hpp"

#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"

namespace geosx
{

using namespace constitutive;
using namespace dataRepository;
using namespace fields;
using namespace finiteElement;

SolidMechanicsConformingFractures::SolidMechanicsConformingFractures( const string & name,
                                                                      Group * const parent ):
  ContactSolverBase( name, parent ),
  m_contactEnforcementMethod( ContactEnforcementMethod::Penalty )
{
  registerWrapper( viewKeyStruct::contactEnforcementMethodString(), &m_contactEnforcementMethod ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue(ContactEnforcementMethod::Penalty).
    setDescription( " add description here" );
}

void SolidMechanicsConformingFractures::registerDataOnMesh( Group & meshBodies )
{
  // Matteo: Do we need element-based fields added inside?
  ContactSolverBase::registerDataOnMesh( meshBodies );

  using namespace fields::contact;

  forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                    MeshLevel & meshLevel,
                                                    arrayView1d< string const > const & regionNames )
  {
    NodeManager & nodes = meshLevel.getNodeManager();

    nodes.registerField< contact::traction > ( getName() ).
      reference().resizeDimension< 1 >( 3 );
    
    nodes.registerField< contact::deltaTraction > ( getName() ).
      reference().resizeDimension< 1 >( 3 );

    nodes.registerField< contact::dispJump > ( getName() ).
      reference().resizeDimension< 1 >( 3 );

    nodes.registerField< contact::oldDispJump > ( getName() ).
      reference().resizeDimension< 1 >( 3 );

    nodes.registerField< contact::dTraction_dJump >( getName() ).
        reference().resizeDimension< 1, 2 >( 3, 3 );

    nodes.registerField< contact::fractureState >( getName() );

    nodes.registerField< contact::oldFractureState >( getName() );

    ElementRegionManager & elemManager = meshLevel.getElemManager();
    elemManager.forElementSubRegions< FaceElementSubRegion >( regionNames, [&] ( localIndex const,
                                                                                 SurfaceElementSubRegion & subRegion )
    {

      /// ALEKS: you may not need this.
      subRegion.registerWrapper< array3d< real64 > >( viewKeyStruct::rotationMatrixString() ).
        setPlotLevel( PlotLevel::NOPLOT ).
        setRegisteringObjects( this->getName()).
        setDescription( "An array that holds the rotation matrices on the fracture." ).
        reference().resizeDimension< 1, 2 >( 3, 3 );

    } );

  } );
}

void SolidMechanicsConformingFractures::setConstitutiveNames( ElementSubRegionBase & subRegion ) const
{
  subRegion.registerWrapper< string >( viewKeyStruct::contactRelationNameString() ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setSizedFromParent( 0 );

  string & contactRelationName = subRegion.getReference< string >( viewKeyStruct::contactRelationNameString() );
  contactRelationName = this->m_contactRelationName;
  GEOSX_ERROR_IF( contactRelationName.empty(), GEOSX_FMT( "Solid model not found on subregion {}", subRegion.getName() ) );
}


void SolidMechanicsConformingFractures::initializePreSubGroups()
{
  SolverBase::initializePreSubGroups();

}

void SolidMechanicsConformingFractures::setupSystem( DomainPartition & domain,
                                                     DofManager & dofManager,
                                                     CRSMatrix< real64, globalIndex > & localMatrix,
                                                     ParallelVector & rhs,
                                                     ParallelVector & solution,
                                                     bool const setSparsity )
{
  GEOSX_MARK_FUNCTION;

  if( m_precond )
  {
    m_precond->clear();
  }

  // setup monolithic coupled system
  SolverBase::setupSystem( domain, dofManager, localMatrix, rhs, solution, setSparsity );
}

void SolidMechanicsConformingFractures::implicitStepSetup( real64 const & time_n,
                                                           real64 const & dt,
                                                           DomainPartition & domain )
{
  computeRotationMatrices( domain );

  m_solidSolver->implicitStepSetup( time_n, dt, domain );
}

void SolidMechanicsConformingFractures::implicitStepComplete( real64 const & time_n,
                                                              real64 const & dt,
                                                              DomainPartition & domain )
{
  if( m_setupSolidSolverDofs )
  {
    m_solidSolver->implicitStepComplete( time_n, dt, domain );
  }

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & meshLevel,
                                                                arrayView1d< string const > const & )
  {
    meshLevel.getElemManager().forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion & subRegion )
    {
      arrayView2d< real64 const > const & dispJump = subRegion.getField< contact::dispJump >();
      arrayView2d< real64 > const & oldDispJump = subRegion.getField< contact::oldDispJump >();
      arrayView1d< integer const > const & fractureState = subRegion.getField< contact::fractureState >();
      arrayView1d< integer > const & oldFractureState = subRegion.getField< contact::oldFractureState >();

      forAll< parallelHostPolicy >( subRegion.size(), [=] ( localIndex const kfe )
      {
        for( localIndex i = 0; i < 3; ++i )
        {
          oldDispJump[kfe][i] = dispJump[kfe][i];
        }
        oldFractureState[kfe] = fractureState[kfe];
      } );
    } );


    /// ALEKS: you may want to sync some fields so I left this here as an example.

    // Need a synchronization of deltaTraction as will be used in AssembleStabilization
    // FieldIdentifiers fieldsToBeSync;
    // fieldsToBeSync.addElementFields( { contact::deltaTraction::key() },
    //                                  { getFractureRegionName() } );

    // CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync,
    //                                                      mesh,
    //                                                      domain.getNeighbors(),
    //                                                      true );

  } );

  GEOSX_LOG_LEVEL_RANK_0( 1, " ***** ImplicitStepComplete *****" );
}

void SolidMechanicsConformingFractures::postProcessInput()
{
  m_solidSolver = &this->getParent().getGroup< SolidMechanicsLagrangianFEM >( m_solidSolverName );
  SolverBase::postProcessInput();
}

SolidMechanicsConformingFractures::~SolidMechanicsConformingFractures()
{
  // TODO Auto-generated destructor stub
}

void SolidMechanicsConformingFractures::resetStateToBeginningOfStep( DomainPartition & domain )
{
  m_solidSolver->resetStateToBeginningOfStep( domain );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & meshLevel,
                                                                arrayView1d< string const > const & )
  {
    NodeManager & nodes = meshLevel.getNodeManager();

    arrayView2d< real64 > const & traction = nodes.getField< contact::traction >();
    arrayView2d< real64 > const & deltaTraction = nodes.getField< contact::deltaTraction >();
    arrayView2d< real64 > const & dispJump = nodes.getField< contact::dispJump >();
    arrayView2d< real64 const > const & oldDispJump = nodes.getField< contact::oldDispJump >();
    arrayView1d< integer > const & fractureState = nodes.getField< contact::fractureState >();
    arrayView1d< integer const > const & oldFractureState = nodes.getField< contact::oldFractureState >();

    forAll< parallelHostPolicy >( nodes.size(), [=] ( localIndex const kn )
    {
        for( localIndex i = 0; i < 3; ++i )
        {
          traction[kn][i] -= deltaTraction[kn][i];
          deltaTraction[kn][i] = 0.0;
          dispJump[kn][i] = oldDispJump[kn][i];
        }
        fractureState[kn] = oldFractureState[kn];
    } );
  } );
}

void SolidMechanicsConformingFractures::computeFaceDisplacementJump( DomainPartition & domain ) const
{
  // TODO: rewrite for node-based dispJump
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & meshLevel,
                                                                arrayView1d< string const > const & regionNames )
  {
    NodeManager const & nodeManager = meshLevel.getNodeManager();
    FaceManager & faceManager = meshLevel.getFaceManager();
    ElementRegionManager & elemManager = meshLevel.getElemManager();

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
        arrayView2d< localIndex const > const & elemsToFaces = subRegion.faceList();
        arrayView2d< real64 > const & dispJump = subRegion.getField< contact::dispJump >();
        arrayView1d< real64 const > const & area = subRegion.getElementArea().toViewConst();

        forAll< parallelHostPolicy >( subRegion.size(), [=] ( localIndex const kfe )
        {
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

void SolidMechanicsConformingFractures::setupDofs( DomainPartition const & domain,
                                                   DofManager & dofManager ) const
{
  GEOSX_MARK_FUNCTION;
  if( m_setupSolidSolverDofs )
  {
    m_solidSolver->setupDofs( domain, dofManager );
  }

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

  if (m_contactEnforcementMethod == ContactEnforcementMethod::Penalty)
  {
    dofManager.addField( solidMechanics::totalDisplacement::key(),
                          FieldLocation::Node,
                          3,
                          meshTargets );

    dofManager.addCoupling( solidMechanics::totalDisplacement::key(),
                            solidMechanics::totalDisplacement::key(),
                            DofManager::Connector::Elem,
                            meshTargets );

    // dofManager.addField(  contact::dispJump::key(),
    //                       FieldLocation::Node,
    //                       3,
    //                       meshTargets );

    // dofManager.addCoupling( contact::dispJump::key(),
    //                         contact::dispJump::key(),
    //                         DofManager::Connector::Elem,
    //                         meshTargets );
  }
}

void SolidMechanicsConformingFractures::assembleSystem( real64 const time,
                                                        real64 const dt,
                                                        DomainPartition & domain,
                                                        DofManager const & dofManager,
                                                        CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                        arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  synchronizeFractureState( domain );

  m_solidSolver->assembleSystem( time,
                                 dt,
                                 domain,
                                 dofManager,
                                 localMatrix,
                                 localRhs );

  /// ALEKS: I would call here the functions to assemble contact constraints
  if (m_contactEnforcementMethod == ContactEnforcementMethod::Penalty)
  {
    assemblePenalizedContact(time, dt, domain, dofManager, localMatrix, localRhs);
  }
  /*else if (m_contactEnforcementMethod == ContactEnforcementMethod::NodalLagrangeMultiplier)
  {
    assembleNodalLagrangeMultiplierContact();
  }
  else if (m_contactEnforcementMethod == ContactEnforcementMethod::FaceLagrangeMultiplier)
  {
    assembleFaceLagrangeMultiplierContact();
  }*/

}

void SolidMechanicsConformingFractures::assemblePenalizedContact( real64 const time,
                                                                  real64 const dt,
                                                                  DomainPartition & domain,
                                                                  DofManager const & dofManager,
                                                                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                                  arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;


  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & meshLevel,
                                                                arrayView1d< string const > const & regionNames )
  {
    //FaceManager const & faceManager = meshLevel.getFaceManager();
    //NodeManager const & nodeManager = meshLevel.getNodeManager();
    //ElementRegionManager const & elemManager = meshLevel.getElemManager();

    //ArrayOfArraysView< localIndex const > const faceToNodeMap = faceManager.nodeList().toViewConst();

    //string const & tracDofKey = dofManager.getKey( contact::traction::key() );
    //string const & dispDofKey = dofManager.getKey( solidMechanics::totalDisplacement::key() );

    //arrayView1d< globalIndex const > const & dispDofNumber = nodeManager.getReference< globalIndex_array >( dispDofKey );
    //globalIndex const rankOffset = dofManager.rankOffset();

    // Get the coordinates for all nodes
    //arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition = nodeManager.referencePosition();

  });
}

real64 SolidMechanicsConformingFractures::calculateResidualNorm( real64 const & GEOSX_UNUSED_PARAM( time ),
                                                                 real64 const & GEOSX_UNUSED_PARAM( dt ),
                                                                 DomainPartition const & domain,
                                                                 DofManager const & dofManager,
                                                                 arrayView1d< real64 const > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  return 0.0;
}

void SolidMechanicsConformingFractures::computeRotationMatrices( DomainPartition & domain ) const
{
  /// ALEKS: I guess you may need this function in case you want, for example, to compute the jump.
  // It should probably be somewhere else though. Plus I am not sure we need to carry around the matrix tbh.

  GEOSX_MARK_FUNCTION;
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & meshLevel,
                                                                arrayView1d< string const > const & regionNames )
  {
    FaceManager const & faceManager = meshLevel.getFaceManager();
    ElementRegionManager & elemManager = meshLevel.getElemManager();

    arrayView2d< real64 const > const & faceNormal = faceManager.faceNormal();

    elemManager.forElementSubRegions< FaceElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   FaceElementSubRegion & subRegion )
    {
      arrayView2d< localIndex const > const & elemsToFaces = subRegion.faceList();

      arrayView3d< real64 > const &
      rotationMatrix = subRegion.getReference< array3d< real64 > >( viewKeyStruct::rotationMatrixString() );

      forAll< parallelHostPolicy >( subRegion.size(), [=]( localIndex const kfe )
      {
        stackArray1d< real64, 3 > Nbar( 3 );
        Nbar[ 0 ] = faceNormal[elemsToFaces[kfe][0]][0] - faceNormal[elemsToFaces[kfe][1]][0];
        Nbar[ 1 ] = faceNormal[elemsToFaces[kfe][0]][1] - faceNormal[elemsToFaces[kfe][1]][1];
        Nbar[ 2 ] = faceNormal[elemsToFaces[kfe][0]][2] - faceNormal[elemsToFaces[kfe][1]][2];
        LvArray::tensorOps::normalize< 3 >( Nbar );

        computationalGeometry::RotationMatrix_3D( Nbar.toSliceConst(), rotationMatrix[kfe] );
      } );
    } );
  } );
}

void SolidMechanicsConformingFractures::computeFaceNodalArea( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition,
                                                                  ArrayOfArraysView< localIndex const > const & faceToNodeMap,
                                                                  localIndex const kf0,
                                                                  array1d< real64 > & nodalArea ) const
{
  // I've tried to access the finiteElement::dispatch3D with
  // finiteElement::FiniteElementBase const &
  // fe = fractureSubRegion->getReference< finiteElement::FiniteElementBase >( surfaceGenerator->getDiscretizationName() );
  // but it's either empty (unknown discretization) or for 3D only (e.g., hexahedra)
  GEOSX_MARK_FUNCTION;

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
    GEOSX_ERROR( "SolidMechanicsConformingFractures: face with " << numNodesPerFace << " nodes. Only triangles and quadrilaterals are supported." );
  }
}

void SolidMechanicsConformingFractures::applySystemSolution( DofManager const & dofManager,
                                                             arrayView1d< real64 const > const & localSolution,
                                                             real64 const scalingFactor,
                                                             DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  if( m_setupSolidSolverDofs )
  {
    m_solidSolver->applySystemSolution( dofManager, localSolution, scalingFactor, domain );
  }

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & meshLevel,
                                                                arrayView1d< string const > const & )
  {
    FieldIdentifiers fieldsToBeSync;

    fieldsToBeSync.addElementFields( { contact::dispJump::key() },
                                     { getFractureRegionName() } );

    CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync,
                                                         meshLevel,
                                                         domain.getNeighbors(),
                                                         true );
  } );
}

void SolidMechanicsConformingFractures::updateState( DomainPartition & domain )
{
  //computeFaceDisplacementJump( domain );
}

bool SolidMechanicsConformingFractures::resetConfigurationToDefault( DomainPartition & domain ) const
{
  GEOSX_MARK_FUNCTION;

  using namespace fields::contact;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & meshLevel,
                                                                arrayView1d< string const > const & regionNames )
  {
    NodeManager & nodes = meshLevel.getNodeManager();
    arrayView1d< integer > const & fractureState = nodes.getField< contact::fractureState >();

    forAll< parallelHostPolicy >( nodes.size(), [=] ( localIndex const kn )
    {
      if( fractureState[kn] != FractureState::Open )
      {
        fractureState[kn] = FractureState::Stick;
      }
    } );
  } );
  return false;
}

bool SolidMechanicsConformingFractures::updateConfiguration( DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  using namespace fields::contact;

  int hasConfigurationConverged = true;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & meshLevel,
                                                                arrayView1d< string const > const & regionNames )
  {} );
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

real64 SolidMechanicsConformingFractures::setNextDt( real64 const & currentDt,
                                                     DomainPartition & domain )
{
  GEOSX_UNUSED_VAR( domain );
  return currentDt;
}

REGISTER_CATALOG_ENTRY( SolverBase, SolidMechanicsConformingFractures, string const &, Group * const )

} /* namespace geosx */
