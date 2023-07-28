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

namespace geos
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
  // ContactSolverBase::registerDataOnMesh( meshBodies );
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

    // The tractionDOFs node-based field is used to store the traction DOFs for the contact enforcement without duplicates
    nodes.registerField< contact::tractionDOFs >( getName() );
    //nodes.registerField< contact::dualGridNodalArea >( getName() );

    nodes.registerWrapper< array1d< real64 > >( viewKeyStruct::normalTractionToleranceString() ).
      setPlotLevel( PlotLevel::NOPLOT ).
      setRegisteringObjects( this->getName()).
      setDescription( "An array that holds the normal traction tolerance." );

    nodes.registerWrapper< array1d< real64 > >( viewKeyStruct::normalDisplacementToleranceString() ).
      setPlotLevel( PlotLevel::NOPLOT ).
      setRegisteringObjects( this->getName()).
      setDescription( "An array that holds the normal displacement tolerance." );

    nodes.registerWrapper< array1d< real64 > >( viewKeyStruct::slidingToleranceString() ).
      setPlotLevel( PlotLevel::NOPLOT ).
      setRegisteringObjects( this->getName()).
      setDescription( "An array that holds the sliding tolerance." );

    ElementRegionManager & elemManager = meshLevel.getElemManager();
    elemManager.forElementSubRegions< FaceElementSubRegion >( regionNames, [&] ( localIndex const,
                                                                                 SurfaceElementSubRegion & subRegion )
    {
      // assuming rotation matrices are constant over one fracture element. This could be problematic for nonplanar fractures
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
  GEOS_ERROR_IF( contactRelationName.empty(), GEOS_FMT( "Solid model not found on subregion {}", subRegion.getName() ) );
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
  GEOS_MARK_FUNCTION;

  if( m_precond )
  {
    m_precond->clear();
  }

  GEOS_UNUSED_VAR( setSparsity );
  dofManager.setDomain( domain );
  setupDofs( domain, dofManager );
  dofManager.reorderByRank();

  // Set the sparsity pattern of stiffness matrix blocks Kuu without Kut and Ktt blocks.
  SparsityPattern< globalIndex > patternDiag;
  dofManager.setSparsityPattern( patternDiag );

  computeNodalTractionDOFs( domain, dofManager ); 
  // Get the original row lengths (diagonal blocks only)
  array1d< localIndex > rowLengths( patternDiag.numRows() );
  for( localIndex localRow = 0; localRow < patternDiag.numRows(); ++localRow )
  {
    rowLengths[localRow] = patternDiag.numNonZeros( localRow );
  }

  // Add the number of nonzeros induced by coupling
  // it should modify the entry of rowLengths for all rows corresponding to a displacement 
  // dof that talks // to a lagrange multipllier (B matrix). 
  //addCouplingNumNonzeros( domain, dofManager, rowLengths.toView() );

  // Create a new pattern with enough capacity for coupled matrix
  SparsityPattern< globalIndex > pattern;
  pattern.resizeFromRowCapacities< parallelHostPolicy >( patternDiag.numRows(), patternDiag.numColumns(), rowLengths.data() );

  // Copy the original nonzeros
  for( localIndex localRow = 0; localRow < patternDiag.numRows(); ++localRow )
  {
    globalIndex const * cols = patternDiag.getColumns( localRow ).dataIfContiguous();
    pattern.insertNonZeros( localRow, cols, cols + patternDiag.numNonZeros( localRow ) );
  }

  // Add the nonzeros from coupling
  // addCouplingSparsityPattern( domain, dofManager, pattern.toView() );

  // Finally, steal the pattern into a CRS matrix
  localMatrix.assimilate< parallelDevicePolicy<> >( std::move( pattern ) );
  localMatrix.setName( this->getName() + "/localMatrix" );

  rhs.setName( this->getName() + "/rhs" );
  rhs.create( dofManager.numLocalDofs(), MPI_COMM_GEOSX );

  solution.setName( this->getName() + "/solution" );
  solution.create( dofManager.numLocalDofs(), MPI_COMM_GEOSX );

}

// void SolidMechanicsEmbeddedFractures::addCouplingNumNonzeros( DomainPartition & domain,
//                                                               DofManager & dofManager,
//                                                               arrayView1d< localIndex > const & rowLengths ) const
// {
//   // Find all rows related to the displacement DOFs associated with the fracture
//   forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
//                                                                 MeshLevel const & mesh,
//                                                                 arrayView1d< string const > const & regionNames )
//   {
//     NodeManager const & nodeManager          = mesh.getNodeManager();
//     ElementRegionManager const & elemManager = mesh.getElemManager();

//     string const jumpDofKey = dofManager.getKey( contact::traction::key() );
//     string const dispDofKey = dofManager.getKey( solidMechanics::totalDisplacement::key() );

//     arrayView1d< globalIndex const > const &
//     dispDofNumber =  nodeManager.getReference< globalIndex_array >( dispDofKey );

//     globalIndex const rankOffset = dofManager.rankOffset();

//     SurfaceElementRegion const & fractureRegion = elemManager.getRegion< SurfaceElementRegion >( getFractureRegionName() );

//     EmbeddedSurfaceSubRegion const & embeddedSurfaceSubRegion = fractureRegion.getSubRegion< EmbeddedSurfaceSubRegion >( 0 );

//     arrayView1d< globalIndex const > const &
//     jumpDofNumber = embeddedSurfaceSubRegion.getReference< array1d< globalIndex > >( jumpDofKey );

//     elemManager.forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const, CellElementSubRegion const & cellElementSubRegion )
//     {

//       SortedArrayView< localIndex const > const fracturedElements = cellElementSubRegion.fracturedElementsList();

//       ArrayOfArraysView< localIndex const > const cellsToEmbeddedSurfaces = cellElementSubRegion.embeddedSurfacesList().toViewConst();

//       localIndex const numDispDof = 3*cellElementSubRegion.numNodesPerElement();

//       for( localIndex ei=0; ei<fracturedElements.size(); ++ei )
//       {
//         localIndex const cellIndex = fracturedElements[ei];

//         localIndex k = cellsToEmbeddedSurfaces[cellIndex][0];
//         localIndex const localRow = LvArray::integerConversion< localIndex >( jumpDofNumber[k] - rankOffset );
//         if( localRow >= 0 && localRow < rowLengths.size() )
//         {
//           for( localIndex i=0; i<3; ++i )
//           {
//             rowLengths[localRow + i] += numDispDof;
//           }
//         }

//         for( localIndex a=0; a<cellElementSubRegion.numNodesPerElement(); ++a )
//         {
//           const localIndex & node = cellElementSubRegion.nodeList( cellIndex, a );
//           localIndex const localDispRow = LvArray::integerConversion< localIndex >( dispDofNumber[node] - rankOffset );

//           if( localDispRow >= 0 && localDispRow < rowLengths.size() )
//           {
//             for( int d=0; d<3; ++d )
//             {
//               rowLengths[localDispRow + d] += 3;
//             }
//           }
//         }
//       }
//     } );
//   } );
// }


// void SolidMechanicsEmbeddedFractures::addCouplingSparsityPattern( DomainPartition const & domain,
//                                                                   DofManager const & dofManager,
//                                                                   SparsityPatternView< globalIndex > const & pattern ) const
// {
//   forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
//                                                                 MeshLevel const & mesh,
//                                                                 arrayView1d< string const > const & regionNames )
//   {
//     NodeManager const & nodeManager          = mesh.getNodeManager();
//     ElementRegionManager const & elemManager = mesh.getElemManager();

//     string const jumpDofKey = dofManager.getKey( contact::dispJump::key() );
//     string const dispDofKey = dofManager.getKey( solidMechanics::totalDisplacement::key() );

//     arrayView1d< globalIndex const > const &
//     dispDofNumber =  nodeManager.getReference< globalIndex_array >( dispDofKey );

//     globalIndex const rankOffset = dofManager.rankOffset();

//     SurfaceElementRegion const & fractureRegion = elemManager.getRegion< SurfaceElementRegion >( getFractureRegionName() );

//     EmbeddedSurfaceSubRegion const & embeddedSurfaceSubRegion = fractureRegion.getSubRegion< EmbeddedSurfaceSubRegion >( 0 );

//     arrayView1d< globalIndex const > const &
//     jumpDofNumber = embeddedSurfaceSubRegion.getReference< array1d< globalIndex > >( jumpDofKey );

//     static constexpr int maxNumDispDof = 3 * 8;

//     elemManager.forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const,
//                                                                                 CellElementSubRegion const & cellElementSubRegion )
//     {

//       SortedArrayView< localIndex const > const fracturedElements = cellElementSubRegion.fracturedElementsList();

//       ArrayOfArraysView< localIndex const > const cellsToEmbeddedSurfaces = cellElementSubRegion.embeddedSurfacesList().toViewConst();

//       localIndex const numDispDof = 3*cellElementSubRegion.numNodesPerElement();

//       for( localIndex ei=0; ei<fracturedElements.size(); ++ei )
//       {
//         localIndex const cellIndex = fracturedElements[ei];
//         localIndex const k = cellsToEmbeddedSurfaces[cellIndex][0];

//         // working arrays
//         stackArray1d< globalIndex, maxNumDispDof > eqnRowIndicesDisp ( numDispDof );
//         stackArray1d< globalIndex, 3 > eqnRowIndicesJump( 3 );
//         stackArray1d< globalIndex, maxNumDispDof > dofColIndicesDisp ( numDispDof );
//         stackArray1d< globalIndex, 3 > dofColIndicesJump( 3 );

//         for( localIndex idof = 0; idof < 3; ++idof )
//         {
//           eqnRowIndicesJump[idof] = jumpDofNumber[k] + idof - rankOffset;
//           dofColIndicesJump[idof] = jumpDofNumber[k] + idof;
//         }

//         for( localIndex a=0; a<cellElementSubRegion.numNodesPerElement(); ++a )
//         {
//           const localIndex & node = cellElementSubRegion.nodeList( cellIndex, a );
//           for( localIndex idof = 0; idof < 3; ++idof )
//           {
//             eqnRowIndicesDisp[3*a + idof] = dispDofNumber[node] + idof - rankOffset;
//             dofColIndicesDisp[3*a + idof] = dispDofNumber[node] + idof;
//           }
//         }

//         for( localIndex i = 0; i < eqnRowIndicesDisp.size(); ++i )
//         {
//           if( eqnRowIndicesDisp[i] >= 0 && eqnRowIndicesDisp[i] < pattern.numRows() )
//           {
//             for( localIndex j = 0; j < dofColIndicesJump.size(); ++j )
//             {
//               pattern.insertNonZero( eqnRowIndicesDisp[i], dofColIndicesJump[j] );
//             }
//           }
//         }

//         for( localIndex i = 0; i < eqnRowIndicesJump.size(); ++i )
//         {
//           if( eqnRowIndicesJump[i] >= 0 && eqnRowIndicesJump[i] < pattern.numRows() )
//           {
//             for( localIndex j=0; j < dofColIndicesDisp.size(); ++j )
//             {
//               pattern.insertNonZero( eqnRowIndicesJump[i], dofColIndicesDisp[j] );
//             }
//           }
//         }
//       }

//     } );
//   } );

// }

void SolidMechanicsConformingFractures::implicitStepSetup( real64 const & time_n,
                                                           real64 const & dt,
                                                           DomainPartition & domain )
{
  computeRotationMatrices( domain );
  computeNodalDisplacementJump( domain );

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
                                                                arrayView1d< string const > const & regionNames )
  {
    NodeManager & nodeManager = meshLevel.getNodeManager();
    FaceManager const & faceManager = meshLevel.getFaceManager();
    ElementRegionManager & elemManager = meshLevel.getElemManager();
    // node-based fields
    arrayView2d< real64 const > const & dispJump = nodeManager.getField< contact::dispJump >();
    arrayView2d< real64 > const & oldDispJump = nodeManager.getField< contact::oldDispJump >();
    arrayView1d< integer const > const & fractureState = nodeManager.getField< contact::fractureState >();
    arrayView1d< integer > const & oldFractureState = nodeManager.getField< contact::oldFractureState >();

    ArrayOfArraysView< localIndex const > const faceToNodeMap = faceManager.nodeList().toViewConst();
    elemManager.forElementSubRegions< FaceElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   FaceElementSubRegion & subRegion )
    {
      arrayView2d< localIndex const > const & elemsToFaces = subRegion.faceList();
      forAll< parallelHostPolicy >( subRegion.size(), [=] ( localIndex const kfe )
      {
        localIndex const numNodesPerFace = faceToNodeMap.sizeOfArray( elemsToFaces[kfe][0] );
        for ( localIndex lagIndex = 0; lagIndex < numNodesPerFace; ++lagIndex )
        {
          for( localIndex faceSide = 0; faceSide < 2; ++faceSide )
          {
            localIndex const nodeIndex = faceToNodeMap( elemsToFaces[kfe][faceSide], lagIndex ); // global index of the node
            for( localIndex i = 0; i < 3; ++i )
            {
              oldDispJump[nodeIndex][i] = dispJump[nodeIndex][i];
            }
            oldFractureState[nodeIndex] = fractureState[nodeIndex];
          }
        }
      } );
    } );
    // Need a synchronization of deltaTraction as will be used in AssembleStabilization
    FieldIdentifiers fieldsToBeSync;
    CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync,
                                                          meshLevel,
                                                          domain.getNeighbors(),
                                                          true );
  } );

  GEOS_LOG_LEVEL_RANK_0( 1, " ***** ImplicitStepComplete *****" );
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

void SolidMechanicsConformingFractures::computeNodalDisplacementJump( DomainPartition & domain ) const
{
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & meshLevel,
                                                                arrayView1d< string const > const & regionNames )
  {
    NodeManager & nodeManager = meshLevel.getNodeManager();
    FaceManager & faceManager = meshLevel.getFaceManager();
    ElementRegionManager & elemManager = meshLevel.getElemManager();

    // Get the coordinates for all nodes
    arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition = nodeManager.referencePosition();
    ArrayOfArraysView< localIndex const > const faceToNodeMap = faceManager.nodeList().toViewConst();
    arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const & u =
      nodeManager.getField< solidMechanics::totalDisplacement >();

    // Get node-based fields
    arrayView2d< real64 > const & dispJump = nodeManager.getField< contact::dispJump >();

    elemManager.forElementSubRegions< FaceElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   FaceElementSubRegion & subRegion )
    {
      arrayView3d< real64 > const &
      rotationMatrix = subRegion.getReference< array3d< real64 > >( viewKeyStruct::rotationMatrixString() );
      arrayView2d< localIndex const > const & elemsToFaces = subRegion.faceList();
      arrayView1d< real64 const > const & area = subRegion.getElementArea().toViewConst();
      arrayView1d< integer const > const & ghostRank = subRegion.ghostRank();

      forAll< parallelHostPolicy >( subRegion.size(), [=] ( localIndex const kfe )
      {
        // Contact constraints
        if (ghostRank[ kfe ] < 0)
        {
          localIndex const numNodesPerFace = faceToNodeMap.sizeOfArray( elemsToFaces[kfe][0] );
          array1d< real64 > nodalArea0, nodalArea1;
          computeFaceNodalArea( nodePosition, faceToNodeMap, elemsToFaces[kfe][0], nodalArea0 );
          computeFaceNodalArea( nodePosition, faceToNodeMap, elemsToFaces[kfe][1], nodalArea1 );

          real64 globalJumpTemp[ 3 ] = { 0 };

          for( localIndex lagIndex = 0; lagIndex < numNodesPerFace; ++lagIndex )
          {
            for( localIndex faceSide = 0; faceSide < 2; ++faceSide )
            {
              // Get the node index for the current node at the both side of the face
              localIndex const nodeIndex = faceToNodeMap( elemsToFaces[kfe][faceSide], lagIndex ); // global index of the node
              for( localIndex i = 0; i < 3; ++i )
              {
                globalJumpTemp[ i ] +=
                    -u[faceToNodeMap( elemsToFaces[kfe][0], lagIndex )][i] * nodalArea0[lagIndex]
                    +u[faceToNodeMap( elemsToFaces[kfe][1], lagIndex )][i] * nodalArea1[lagIndex];
              }
              
              real64 dispJumpTemp[ 3 ];
              LvArray::tensorOps::Ri_eq_AjiBj< 3, 3 >( dispJumpTemp, rotationMatrix[ kfe ], globalJumpTemp );
              LvArray::tensorOps::copy< 3 >( dispJump[ nodeIndex ], dispJumpTemp );
            }
          }
        }
      } );
    } );
  } );
}

void SolidMechanicsConformingFractures::computeNodalTractionDOFs( DomainPartition & domain,
                                                   DofManager & dofManager ) const
{
  // compute the nodal traction DOFs
  GEOS_MARK_FUNCTION;
  std::set< localIndex > uniqueNodeSet;
  localIndex const numComponents = 3;
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & meshLevel,
                                                                arrayView1d< string const > const & regionNames ) 
  {
    NodeManager & nodeManager = meshLevel.getNodeManager();
    FaceManager const & faceManager = meshLevel.getFaceManager();
    // we need to recompute the dofs for tractions. For each pair of nodes, we use the smaller value of the two dofs as the 
    // dof for the traction.

    arrayView1d< localIndex > const & tracDofNumber = nodeManager.getField< contact::tractionDOFs >();
    ArrayOfArraysView< localIndex const > const faceToNodeMap = faceManager.nodeList().toViewConst();

    meshLevel.getElemManager().forElementSubRegions< FaceElementSubRegion >( regionNames,
                                                              [&]( localIndex,
                                                                   FaceElementSubRegion & subRegion )
    {
      arrayView2d< localIndex const > const & elemsToFaces = subRegion.faceList();
      forAll< parallelHostPolicy >( subRegion.size(), [&] ( localIndex const kfe )
      {
        localIndex const numNodesPerFace = faceToNodeMap.sizeOfArray( elemsToFaces[kfe][0] );

        for ( localIndex lagIndex = 0; lagIndex < numNodesPerFace; ++lagIndex )
        {
          // Get the smaller value of the two dofs as the dof for the traction
          // Only get the index of nodes on the first face
          localIndex const firstTractionIndex = faceToNodeMap( elemsToFaces[kfe][0], lagIndex );
          uniqueNodeSet.insert( firstTractionIndex );
          tracDofNumber[firstTractionIndex] = dofManager.numLocalDofs() + uniqueNodeSet.size()*numComponents;
        }
      } );
    } );
  } );
  
}

void SolidMechanicsConformingFractures::setupDofs( DomainPartition const & domain,
                                                   DofManager & dofManager ) const
{
  GEOS_MARK_FUNCTION;
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
  }
  else if (m_contactEnforcementMethod == ContactEnforcementMethod::NodalLagrangeMultiplier)
  {

    // dofManager.addField( contact::traction::key(),
    //                    FieldLocation::Node,
    //                    3,
    //                    meshTargets );

    // dofManager.addCoupling( contact::traction::key(),
    //                         contact::traction::key(),
    //                         DofManager::Connector::Elem,
    //                         meshTargets );
                         
  }
  else
  {
    GEOS_ERROR( "Unknown contact enforcement method." );
  }
}

void SolidMechanicsConformingFractures::assembleSystem( real64 const time,
                                                        real64 const dt,
                                                        DomainPartition & domain,
                                                        DofManager const & dofManager,
                                                        CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                        arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  synchronizeFractureState( domain );

  m_solidSolver->assembleSystem( time,
                                 dt,
                                 domain,
                                 dofManager,
                                 localMatrix,
                                 localRhs );

  if (m_contactEnforcementMethod == ContactEnforcementMethod::Penalty)
  {
    // TODO: INTENTIONALLY LEFT BLANK FOR PENALIZED CONTACT
  }
  else if (m_contactEnforcementMethod == ContactEnforcementMethod::NodalLagrangeMultiplier)
  {
    assembleNodalLagrangeMultiplierContact( domain, dofManager, localMatrix, localRhs);
  }
  
  else if (m_contactEnforcementMethod == ContactEnforcementMethod::FaceLagrangeMultiplier)
  {
    // TODO: CALL assembleFaceLagrangeMultiplierContact()
    //assembleFaceLagrangeMultiplierContact();
  }

}

void SolidMechanicsConformingFractures::assembleNodalLagrangeMultiplierContact( DomainPartition & domain,
                                                                                DofManager const & dofManager,
                                                                                CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                                                arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    assembleForceResidualDerivativeWrtTraction( mesh, regionNames, dofManager, localMatrix, localRhs );
    assembleTractionResidualDerivativeWrtDisplacementAndTraction( mesh, regionNames, dofManager, localMatrix, localRhs );
  } );

}

void SolidMechanicsConformingFractures::
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

  // nodal-based fields
  arrayView1d< globalIndex const > const & dispDofNumber = nodeManager.getReference< globalIndex_array >( dispDofKey );
  arrayView1d< globalIndex const > const & tracDofNumber = nodeManager.getReference< globalIndex_array >( tracDofKey );

  arrayView2d< real64 const > const & dispJump = nodeManager.getField< contact::dispJump >();
  arrayView2d< real64 const > const & previousDispJump = nodeManager.getField< contact::oldDispJump >();
  arrayView1d< real64 const > const & slidingTolerance = nodeManager.getReference< array1d< real64 > >( viewKeyStruct::slidingToleranceString() );  
  //arrayView1d< real64 const > const & nodalArea = nodeManager.getField< contact::dualGridNodalArea >();

  arrayView2d< real64 const > const & traction = nodeManager.getField< contact::traction >();
  arrayView1d< integer const > const & fractureState = nodeManager.getField< contact::fractureState >();
  globalIndex const rankOffset = dofManager.rankOffset();

  // Get the coordinates for all nodes
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition = nodeManager.referencePosition();

  elemManager.forElementSubRegions< FaceElementSubRegion >( regionNames,
                                                            [&]( localIndex const,
                                                                 FaceElementSubRegion const & subRegion )
  {
    // Loop through all subregions of the face elements
    // face-based fields
    ContactBase const & contact = getConstitutiveModel< ContactBase >( subRegion, m_contactRelationName );
    arrayView1d< integer const > const & ghostRank = subRegion.ghostRank();
    arrayView3d< real64 const > const &
    rotationMatrix = subRegion.getReference< array3d< real64 > >( viewKeyStruct::rotationMatrixString() );
    arrayView2d< localIndex const > const & elemsToFaces = subRegion.faceList();
    
    constitutiveUpdatePassThru( contact, [&] ( auto & castedContact )
    {
      using ContactType = TYPEOFREF( castedContact );
      typename ContactType::KernelWrapper contactWrapper = castedContact.createKernelWrapper();

      forAll< parallelHostPolicy >( subRegion.size(), [=] ( localIndex const kfe )
      {
        // Loop through all face elements in the subregion
        if (ghostRank[kfe] < 0)
        {
          // Filter out ghost elements
          localIndex const numNodesPerFace = faceToNodeMap.sizeOfArray( elemsToFaces[kfe][0] );
          array1d< real64 > nodalArea;

          for ( localIndex lagIndex = 0; lagIndex < numNodesPerFace; ++lagIndex )
          {
            // loop through all nodes on the face element; operations for each node are duplicated for nodes living on the other side of the face element
            //globalIndex elemDOF[3] = { tracDofNumber[nodeIndex], tracDofNumber[nodeIndex] + 1, tracDofNumber[nodeIndex] + 2 }; // traction DOF for the node
            for( localIndex faceSide = 0; faceSide < 2; ++faceSide )
            {

              globalIndex nodeDOF[3]; // displacement DOF for the node 
              globalIndex elemDOF[3]; // displacement DOF for the traction node

              stackArray1d< real64, 3 > elemRHS( 3 );
              stackArray2d< real64, 3 * 3 > dRdU( 3 , 3 );
              stackArray2d< real64, 3 * 3 > dRdT( 3 , 3 );
              
              localIndex const nodeIndex = faceToNodeMap( elemsToFaces[kfe][faceSide], lagIndex ); // global index of the node
              computeFaceNodalArea( nodePosition, faceToNodeMap, elemsToFaces[kfe][faceSide], nodalArea ); // compute nodal area (dual mesh size)
                            
              switch( fractureState[nodeIndex] )
              {
                case contact::FractureState::Stick:
                {
                  // Node is in stick state
                  for( localIndex i = 0; i < 3; ++i )
                  {
                    if (faceSide == 0)
                    {
                      elemDOF[i] = tracDofNumber[nodeIndex] + i;
                      elemRHS[i] = +nodalArea[lagIndex] * (dispJump[nodeIndex][i] - (i == 0 ? 0 : previousDispJump[nodeIndex][i]));
                    }

                    nodeDOF[i] = dispDofNumber[nodeIndex] + i;
                    dRdU( 0, i ) = -nodalArea[lagIndex] * rotationMatrix( kfe, i, 0 ) * pow( -1, faceSide ); 
                    dRdU( 1, i ) = -nodalArea[lagIndex] * rotationMatrix( kfe, i, 1 ) * pow( -1, faceSide );
                    dRdU( 2, i ) = -nodalArea[lagIndex] * rotationMatrix( kfe, i, 2 ) * pow( -1, faceSide );
                  }
                  break;                
                }
                case contact::FractureState::Slip:
                case contact::FractureState::NewSlip:
                {
                  // Node is in slip state
                  if (faceSide == 0)
                  {
                    elemRHS[0] = +nodalArea[lagIndex] * dispJump[nodeIndex][0];
                  }
                  for( localIndex i = 0; i < 3; ++i )
                  {
                    nodeDOF[i] = dispDofNumber[nodeIndex] + LvArray::integerConversion< globalIndex >( i );
                    dRdU( 0, i ) = -nodalArea[lagIndex] * rotationMatrix( kfe, i, 0 ) * pow( -1, faceSide ); 
                  } 

                  real64 dLimitTau_dNormalTraction = 0.0;
                  real64 const limitTau = contactWrapper.computeLimitTangentialTractionNorm( traction[nodeIndex][0],
                                                                                            dLimitTau_dNormalTraction );
                  real64 sliding[ 2 ] = { dispJump[nodeIndex][1] - previousDispJump[nodeIndex][1], dispJump[nodeIndex][2] - previousDispJump[nodeIndex][2] };
                  real64 slidingNorm = sqrt( sliding[ 0 ]*sliding[ 0 ] + sliding[ 1 ]*sliding[ 1 ] );

                  if(!( ( m_nonlinearSolverParameters.m_numNewtonIterations == 0 ) && ( fractureState[nodeIndex] == contact::FractureState::NewSlip ) )
                    && slidingNorm > slidingTolerance[nodeIndex] )
                  {
                    real64 dUdgT[ 3 ];
                    if (faceSide == 0)
                    {
                      for( localIndex i = 1; i < 3; ++i )
                      {
                        elemRHS[i] = +nodalArea[lagIndex] * ( traction[nodeIndex][i] - limitTau * sliding[ i-1 ] / slidingNorm );
                      }
                      // A symmetric 2x2 matrix.
                      dUdgT[ 0 ] = (slidingNorm * slidingNorm - sliding[ 0 ] * sliding[ 0 ]) * limitTau / std::pow( slidingNorm, 3 );
                      dUdgT[ 1 ] = (slidingNorm * slidingNorm - sliding[ 1 ] * sliding[ 1 ]) * limitTau / std::pow( slidingNorm, 3 );
                      dUdgT[ 2 ] = -sliding[ 0 ] * sliding[ 1 ] * limitTau / std::pow( slidingNorm, 3 );

                      for( localIndex i = 1; i < 3; ++i )
                      {
                        dRdT( i, 0 ) = nodalArea[lagIndex] * dLimitTau_dNormalTraction * sliding[ i-1 ] / slidingNorm;
                        dRdT( i, i ) = nodalArea[lagIndex];
                      }
                    }

                    {
                      for( localIndex a = 0; a < numNodesPerFace; ++a )
                      {
                        for( localIndex i = 0; i < 3; ++i )
                        {
                          real64 const localRowB[ 2 ] = { rotationMatrix( kfe, i, 1 ), rotationMatrix( kfe, i, 2 ) };
                          real64 localRowE[ 2 ];
                          LvArray::tensorOps::Ri_eq_symAijBj< 2 >( localRowE, dUdgT, localRowB );

                          dRdU( 1, i ) = nodalArea[a] * localRowE[ 0 ] * pow( -1, faceSide );
                          dRdU( 2, i ) = nodalArea[a] * localRowE[ 1 ] * pow( -1, faceSide );
                        }
                      }
                    }
                    for( localIndex i = 1; i < 3; ++i )
                    {
                      dRdT( i, 0 ) = nodalArea[lagIndex] * dLimitTau_dNormalTraction * sliding[ i-1 ] / slidingNorm;
                      dRdT( i, i ) = nodalArea[lagIndex];
                    }
                  }
                  else
                  {
                    // directly change from the stick mode to the slip mode
                    real64 vaux[ 2 ] = { traction[nodeIndex][1], traction[nodeIndex][2] };
                    real64 vauxNorm = sqrt( vaux[ 0 ]*vaux[ 0 ] + vaux[ 1 ]*vaux[ 1 ] );
                    if( vauxNorm > 0.0 )
                    {
                      for( localIndex i = 1; i < 3; ++i )
                      {
                        elemRHS[i] = +nodalArea[lagIndex] * ( traction[nodeIndex][i] - limitTau * vaux[ i-1 ] / vauxNorm );
                      }
                      for( localIndex i = 1; i < 3; ++i )
                      {
                        dRdT( i, i ) = nodalArea[lagIndex];
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
                        dRdT( i, i ) = nodalArea[lagIndex];
                      }
                    }
                  }
                  break;
                }
                case contact::FractureState::Open:
                {
                  // Node is in open state
                  for( localIndex i = 0; i < 3; ++i )
                  {
                    if (faceSide == 0)
                    {
                      elemRHS[i] = +nodalArea[lagIndex] * traction[nodeIndex][i];
                      dRdT( i, i ) = nodalArea[lagIndex];
                    }
                  }
                }
                break;
              }

              localIndex const localRow = LvArray::integerConversion< localIndex >( elemDOF[0] - rankOffset );
              // assemble the residual and Jacobian coefficients associated with the node
              for( localIndex idof = 0; idof < 3; ++idof )
              {
                localRhs[localRow + idof] += elemRHS[idof];

                if( fractureState[nodeIndex] != contact::FractureState::Open )
                {
                  localMatrix.addToRowBinarySearchUnsorted< serialAtomic >( localRow + idof,
                                                                            nodeDOF,
                                                                            dRdU[idof].dataIfContiguous(),
                                                                            3 );
                }

                if( fractureState[nodeIndex] != contact::FractureState::Stick )
                {
                  localMatrix.addToRow< serialAtomic >( localRow + idof,
                                                        elemDOF,
                                                        dRdT[idof].dataIfContiguous(),
                                                        3 );
                }
              }
            }
          }
        }
      });
    });
  });
}

void SolidMechanicsConformingFractures::
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
  globalIndex const rankOffset = dofManager.rankOffset();

  // Get the coordinates for all nodes
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition = nodeManager.referencePosition();
  // nodal-based fields
  arrayView1d< globalIndex const > const & dispDofNumber = nodeManager.getReference< globalIndex_array >( dispDofKey );
  arrayView1d< globalIndex const > const & tracDofNumber = nodeManager.getReference< globalIndex_array >( tracDofKey );
  arrayView2d< real64 const > const & traction = nodeManager.getReference< array2d< real64 > >( contact::traction::key() );
  
  elemManager.forElementSubRegions< FaceElementSubRegion >( regionNames,
                                                            [&]( localIndex const,
                                                                 FaceElementSubRegion const & subRegion )
  {
    arrayView3d< real64 const > const &
    rotationMatrix = subRegion.getReference< array3d< real64 > >( viewKeyStruct::rotationMatrixString() );
    arrayView2d< localIndex const > const & elemsToFaces = subRegion.faceList();
    arrayView1d< integer const > const & ghostRank = subRegion.ghostRank();

    constexpr localIndex TriangularPermutation[3] = { 0, 1, 2 };
    constexpr localIndex QuadrilateralPermutation[4] = { 0, 1, 3, 2 };

    forAll< parallelHostPolicy >( subRegion.size(), [=] ( localIndex const kfe )
    {
      if ( ghostRank[kfe] < 0 )
      {
        localIndex const numNodesPerFace = faceToNodeMap.sizeOfArray( elemsToFaces[kfe][0] );
        localIndex const numQuadraturePointsPerElem = numNodesPerFace==3 ? 1 : 4;

        std::vector<globalIndex> rowDOF(3 * numNodesPerFace);
        std::vector<real64> nodeRHS(3 * numNodesPerFace);

        stackArray2d< real64, 3*4*3 > dRdT( 3*4, 3 );
        globalIndex colDOF[3];
      
        localIndex const * const permutation = ( numNodesPerFace == 3 ) ? TriangularPermutation : QuadrilateralPermutation;
        real64 xLocal[2][4][3]; // TODO: this is a hack to get around the fact that the code below assumes 4 or less quadrature points
                                //      this should be fixed when dealing with arbitrary number of quadrature points

        for ( localIndex lagIndex = 0; lagIndex < numNodesPerFace; ++lagIndex )
        {
          // loop through all nodes on the face element; operations for each node are duplicated for nodes living on the other side of the face element
          // globalIndex elemDOF[3] = { tracDofNumber[nodeIndex], tracDofNumber[nodeIndex] + 1, tracDofNumber[nodeIndex] + 2 }; // traction DOF for the node
      
          for( localIndex faceSide = 0; faceSide < 2; ++faceSide )
          {
            localIndex const nodeIndex = faceToNodeMap( elemsToFaces[kfe][faceSide], lagIndex ); // global index of the node
            for( localIndex i = 0; i < 3; ++i )
            {
              if (faceSide == 0)
              {
                colDOF[i] = tracDofNumber[nodeIndex] + i;
              }
              xLocal[faceSide][lagIndex][i] = nodePosition[ nodeIndex ][i];
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
          for( localIndex faceSide = 0; faceSide < 2; ++faceSide )
          {
            using xLocalTriangle = real64[3][3];
            real64 const detJxW = numNodesPerFace==3 ?
                                  H1_TriangleFace_Lagrange1_Gauss1::transformedQuadratureWeight( q, reinterpret_cast< xLocalTriangle & >( xLocal[faceSide] ) ) :
                                  H1_QuadrilateralFace_Lagrange1_GaussLegendre2::transformedQuadratureWeight( q, xLocal[faceSide] );

            for( localIndex lagIndex = 0; lagIndex < numNodesPerFace; ++lagIndex )
            {
              localIndex const nodeIndex = faceToNodeMap( elemsToFaces[kfe][faceSide], lagIndex ); // global index of the node
              real64 const NaDetJxQ = N[permutation[lagIndex]] * detJxW;
              real64 const localNodalForce[ 3 ] = { traction( nodeIndex, 0 ) * NaDetJxQ,
                                                    traction( nodeIndex, 1 ) * NaDetJxQ,
                                                    traction( nodeIndex, 2 ) * NaDetJxQ };
              real64 globalNodalForce[ 3 ];
              LvArray::tensorOps::Ri_eq_AijBj< 3, 3 >( globalNodalForce, rotationMatrix[ kfe ], localNodalForce );

              for( localIndex i = 0; i < 3; ++i )
              {
                rowDOF[3*lagIndex+i] = dispDofNumber[ nodeIndex ] + i;
                nodeRHS[3*lagIndex+i] = +globalNodalForce[i] * normalSign[ faceSide ];

                // Opposite sign w.r.t. to the same formulation as above
                dRdT( 3*lagIndex+i, 0 ) = rotationMatrix( kfe, i, 0 ) * normalSign[ faceSide ] * NaDetJxQ;
                dRdT( 3*lagIndex+i, 1 ) = rotationMatrix( kfe, i, 1 ) * normalSign[ faceSide ] * NaDetJxQ;
                dRdT( 3*lagIndex+i, 2 ) = rotationMatrix( kfe, i, 2 ) * normalSign[ faceSide ] * NaDetJxQ;
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
      }
    });
  });
}

real64 SolidMechanicsConformingFractures::calculateResidualNorm( real64 const & GEOS_UNUSED_PARAM( time ),
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
    FaceManager const & faceManager = mesh.getFaceManager();
    arrayView1d< globalIndex const > const & dispDofNumber =
      nodeManager.getReference< array1d< globalIndex > >( dofManager.getKey( solidMechanics::totalDisplacement::key() ) );
    ArrayOfArraysView< localIndex const > const faceToNodeMap = faceManager.nodeList().toViewConst();

    string const & dofKey = dofManager.getKey( contact::traction::key() );
    globalIndex const rankOffset = dofManager.rankOffset();

    arrayView1d< integer const > const & elemGhostRank = nodeManager.ghostRank();
    RAJA::ReduceSum< parallelDeviceReduce, real64 > localSum0( 0.0 );
    // 1. Compute the norm of the residuals associated with displacement
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

    string const & tracDofKey = dofManager.getKey( contact::traction::key() );
    arrayView1d< globalIndex const > const & tracDofNumber = nodeManager.getReference< globalIndex_array >( tracDofKey );
    // 2. Compute the norm of the residuals associated with traction
    mesh.getElemManager().forElementSubRegions< FaceElementSubRegion >( regionNames,
                                                                        [&]( localIndex const, FaceElementSubRegion const & subRegion )
    {
      arrayView2d< localIndex const > const & elemsToFaces = subRegion.faceList();
      arrayView1d< integer const > const & ghostRank = subRegion.ghostRank();
      RAJA::ReduceSum< parallelHostReduce, real64 > localSum( 0.0 );
      forAll< parallelHostPolicy >( subRegion.size(), [=] ( localIndex const kfe )
      {
        if( ghostRank[kfe] < 0 )
        {
          localIndex const numNodesPerFace = faceToNodeMap.sizeOfArray( elemsToFaces[kfe][0] );
          for ( localIndex lagIndex = 0; lagIndex < numNodesPerFace; ++lagIndex )
          {
            for( localIndex faceSide = 0; faceSide < 2; ++faceSide )
            {
              localIndex const nodeIndex = faceToNodeMap( elemsToFaces[kfe][faceSide], lagIndex );
              localIndex const localRow = LvArray::integerConversion< localIndex >( tracDofNumber[kfe] - rankOffset );
              for( localIndex dim = 0; dim < 3; ++dim )
              {
                localSum += localRhs[localRow + dim] * localRhs[localRow + dim];
              }
            }
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
  GEOS_LOG_LEVEL_RANK_0( 1, GEOS_FMT( "    ( Rdisplacement, Rtraction, Rtotal ) = ( {:15.6e}, {:15.6e}, {:15.6e} );",
                                      globalResidualNorm[0],
                                      globalResidualNorm[1],
                                      globalResidualNorm[2] ) );
  return globalResidualNorm[2];
}

void SolidMechanicsConformingFractures::computeRotationMatrices( DomainPartition & domain ) const
{
  GEOS_MARK_FUNCTION;
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
  // TODO: This function assumes only triangular and quadrilateral faces. It should be generalized to arbitrary polygons.
  //       For polygons with more than 4 nodes, the area should be computed by using the barycenter of the polygon. 
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
    GEOS_ERROR( "SolidMechanicsConformingFractures: face with " << numNodesPerFace << " nodes. Only triangles and quadrilaterals are supported." );
  }
}

void SolidMechanicsConformingFractures::applySystemSolution( DofManager const & dofManager,
                                                             arrayView1d< real64 const > const & localSolution,
                                                             real64 const scalingFactor,
                                                             DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

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
  computeNodalDisplacementJump( domain );
}

bool SolidMechanicsConformingFractures::resetConfigurationToDefault( DomainPartition & domain ) const
{
  GEOS_MARK_FUNCTION;

  using namespace fields::contact;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & meshLevel,
                                                                arrayView1d< string const > const & regionNames )
  {
    GEOS_UNUSED_VAR( regionNames );
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
  GEOS_MARK_FUNCTION;
  using namespace fields::contact;
  int hasConfigurationConverged = true;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    // node-based fields
    NodeManager & nodeManager = mesh.getNodeManager();
    FaceManager const & faceManager = mesh.getFaceManager();
    ElementRegionManager & elemManager = mesh.getElemManager(); 

    ArrayOfArraysView< localIndex const > const faceToNodeMap = faceManager.nodeList().toViewConst();
    arrayView2d< real64 const > const & traction = nodeManager.getField< contact::traction >();
    arrayView2d< real64 const > const & dispJump = nodeManager.getField< contact::dispJump >();
    arrayView1d< integer > const & fractureState = nodeManager.getField< contact::fractureState >();

    arrayView1d< real64 const > const & normalTractionTolerance =
    nodeManager.getReference< array1d< real64 > >( viewKeyStruct::normalTractionToleranceString() ); 
    arrayView1d< real64 const > const & normalDisplacementTolerance =
    nodeManager.getReference< array1d< real64 > >( viewKeyStruct::normalDisplacementToleranceString() );

    elemManager.forElementSubRegions< FaceElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                FaceElementSubRegion & subRegion )
    {
      // face-based fields
      ContactBase const & contact = getConstitutiveModel< ContactBase >( subRegion, m_contactRelationName );
      arrayView2d< localIndex const > const & elemsToFaces = subRegion.faceList();
      arrayView1d< integer const > const & ghostRank = subRegion.ghostRank();
      RAJA::ReduceMin< parallelHostReduce, integer > checkActiveSetSub( 1 );
      constitutiveUpdatePassThru( contact, [&] ( auto & castedContact )
      {
        using ContactType = TYPEOFREF( castedContact );
        typename ContactType::KernelWrapper contactWrapper = castedContact.createKernelWrapper();

        forAll< parallelHostPolicy >( subRegion.size(), [=] ( localIndex const kfe )
        {
          if( ghostRank[kfe] < 0 )
          {
            localIndex const numNodesPerFace = faceToNodeMap.sizeOfArray( elemsToFaces[kfe][0] );
            array1d< real64 > nodalArea;

            for ( localIndex lagIndex = 0; lagIndex < numNodesPerFace; ++lagIndex )
            {
              for( localIndex faceSide = 0; faceSide < 2; ++faceSide )
              {
                localIndex const nodeIndex = faceToNodeMap( elemsToFaces[kfe][faceSide], lagIndex ); // global index of the node
                integer const originalFractureState = fractureState[nodeIndex];
                if( originalFractureState == contact::FractureState::Open )
                {
                  if( dispJump[nodeIndex][0] > -normalDisplacementTolerance[nodeIndex] )
                  {
                    fractureState[nodeIndex] = contact::FractureState::Open;
                  }
                  else
                  {
                    fractureState[nodeIndex] = contact::FractureState::Stick;
                  }
                }
                else if( traction[nodeIndex][0] > normalTractionTolerance[nodeIndex] )
                {
                  fractureState[nodeIndex] = contact::FractureState::Open;
                }
                else
                {
                  real64 currentTau = sqrt( traction[nodeIndex][1]*traction[nodeIndex][1] + traction[nodeIndex][2]*traction[nodeIndex][2] );

                  real64 dLimitTangentialTractionNorm_dTraction = 0.0;
                  real64 const limitTau =
                    contactWrapper.computeLimitTangentialTractionNorm( traction[nodeIndex][0],
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
                      fractureState[nodeIndex] = contact::FractureState::NewSlip;
                    }
                    else
                    {
                      fractureState[nodeIndex] = contact::FractureState::Slip;
                    }
                  }
                  else
                  {
                    fractureState[nodeIndex] = contact::FractureState::Stick;
                  }
                }
                checkActiveSetSub.min( compareFractureStates( originalFractureState, fractureState[kfe] ) );
              }
            }
          }
        } );
      } );
    } );
  } );

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
  GEOS_UNUSED_VAR( domain );
  return currentDt;
}

REGISTER_CATALOG_ENTRY( SolverBase, SolidMechanicsConformingFractures, string const &, Group * const )

} /* namespace geos */