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
  ContactSolverBase::registerDataOnMesh( meshBodies );

  using namespace fields::contact;

  forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                    MeshLevel & mesh,
                                                    arrayView1d< string const > const & regionNames )
  {
    ElementRegionManager & elemManager = mesh.getElemManager();
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
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & )
  {
    mesh.getElemManager().forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion & subRegion )
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

void SolidMechanicsConformingFractures::computeFaceDisplacementJump( DomainPartition & domain ) const
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
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & )
  {
    FieldIdentifiers fieldsToBeSync;

    fieldsToBeSync.addElementFields( { contact::dispJump::key() },
                                     { getFractureRegionName() } );

    CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync,
                                                         mesh,
                                                         domain.getNeighbors(),
                                                         true );
  } );
}

void SolidMechanicsConformingFractures::updateState( DomainPartition & domain )
{
  computeFaceDisplacementJump( domain );
}

bool SolidMechanicsConformingFractures::resetConfigurationToDefault( DomainPartition & domain ) const
{
  GEOSX_MARK_FUNCTION;

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

bool SolidMechanicsConformingFractures::updateConfiguration( DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  using namespace fields::contact;

  int hasConfigurationConverged = true;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
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
