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
 * @file MeshLevel.cpp
 */

#include "MeshLevel.hpp"

#include "EdgeManager.hpp"
#include "ElementRegionManager.hpp"
#include "NodeManager.hpp"
#include "FaceManager.hpp"

namespace geosx
{
using namespace dataRepository;

MeshLevel::MeshLevel( string const & name,
                      Group * const parent ):
  Group( name, parent ),
  m_nodeManager( new NodeManager( groupStructKeys::nodeManagerString(), this ) ),
  m_edgeManager( new EdgeManager( groupStructKeys::edgeManagerString(), this ) ),
  m_faceManager( new FaceManager( groupStructKeys::faceManagerString(), this ) ),
  m_elementManager( new ElementRegionManager( groupStructKeys::elemManagerString(), this ) ),
  m_embSurfNodeManager( new EmbeddedSurfaceNodeManager( groupStructKeys::embSurfNodeManagerString, this ) ),
  m_embSurfEdgeManager( new EdgeManager( groupStructKeys::embSurfEdgeManagerString, this ) ),
  m_modificationTimestamp( 0 ),
  m_isShallowCopy( false ),
  m_shallowParent( nullptr )
{

  registerGroup( groupStructKeys::nodeManagerString(), m_nodeManager );

  registerGroup( groupStructKeys::edgeManagerString(), m_edgeManager );


  registerGroup< FaceManager >( groupStructKeys::faceManagerString(), m_faceManager );
  m_faceManager->nodeList().setRelatedObject( *m_nodeManager );


  registerGroup< ElementRegionManager >( groupStructKeys::elemManagerString(), m_elementManager );

  registerGroup< EdgeManager >( groupStructKeys::embSurfEdgeManagerString, m_embSurfEdgeManager );

  registerGroup< EmbeddedSurfaceNodeManager >( groupStructKeys::embSurfNodeManagerString, m_embSurfNodeManager );

  registerWrapper< integer >( viewKeys.meshLevel );

  // increment the modification timestamp at mesh level creation
  // this is to make sure that the actions that depend on this timestamp (such as system setup) are performed at the beginning of the
  // simulations
  modified();
}


MeshLevel::MeshLevel( string const & name,
                      Group * const parent,
                      MeshLevel & source ):
  Group( name, parent ),
  m_nodeManager( source.m_nodeManager ),
  m_edgeManager( source.m_edgeManager ),
  m_faceManager( source.m_faceManager ),
  m_elementManager( source.m_elementManager ),
  m_embSurfNodeManager( source.m_embSurfNodeManager ),
  m_embSurfEdgeManager( source.m_embSurfEdgeManager ),
  m_modificationTimestamp( 0 ),
  m_isShallowCopy( true ),
  m_shallowParent( &source )
{
  this->setRestartFlags( RestartFlags::NO_WRITE );

  registerGroup( groupStructKeys::nodeManagerString(), m_nodeManager );

  registerGroup( groupStructKeys::edgeManagerString(), m_edgeManager );


  registerGroup< FaceManager >( groupStructKeys::faceManagerString(), m_faceManager );
  m_faceManager->nodeList().setRelatedObject( *m_nodeManager );


  registerGroup< ElementRegionManager >( groupStructKeys::elemManagerString(), m_elementManager );

  registerGroup< EdgeManager >( groupStructKeys::embSurfEdgeManagerString, m_embSurfEdgeManager );

  registerGroup< EmbeddedSurfaceNodeManager >( groupStructKeys::embSurfNodeManagerString, m_embSurfNodeManager );

  registerWrapper< integer >( viewKeys.meshLevel );

  // increment the modification timestamp at mesh level creation
  // this is to make sure that the actions that depend on this timestamp (such as system setup) are performed at the beginning of the
  // simulations
  modified();
}


MeshLevel::MeshLevel( string const & name,
                      Group * const parent,
                      MeshLevel const & source,
                      CellBlockManagerABC & cellBlockManager,
                      int const order ):
  MeshLevel( name, parent )
{

  GEOSX_MARK_FUNCTION;
  localIndex const numBasisSupportPoints = order+1;


  // Edges
  localIndex const numNodesPerEdge = numBasisSupportPoints;
  localIndex const numNonVertexNodesPerEdge = numNodesPerEdge - 2;

  m_edgeManager->resize( source.m_edgeManager->size() );
  localIndex const numInternalEdgeNodes = m_edgeManager->size() * numNonVertexNodesPerEdge;

  m_faceManager->resize( source.m_faceManager->size() );
  m_faceManager->edgeList() = source.m_faceManager->edgeList();
  m_faceManager->faceCenter() = source.m_faceManager->faceCenter();
  m_faceManager->faceNormal() = source.m_faceManager->faceNormal();
  m_faceManager->getDomainBoundaryIndicator() = source.m_faceManager->getDomainBoundaryIndicator();
  m_faceManager->elementList() = source.m_faceManager->elementList();

  // Faces
  ArrayOfArrays< localIndex > & faceToNodeMapNew = m_faceManager->nodeList();
  ArrayOfArraysView< localIndex const > const & facesToEdges = m_faceManager->edgeList().toViewConst();
  localIndex const numNodesPerFace = pow( order+1, 2 );

  // number of elements in each row of the map as capacity
  array1d< localIndex > counts( faceToNodeMapNew.size());
  counts.setValues< parallelHostPolicy >( numNodesPerFace );

  //  reconstructs the faceToNodeMap with the provided capacity in counts
  faceToNodeMapNew.resizeFromCapacities< parallelHostPolicy >( faceToNodeMapNew.size(), counts.data() );

  // setup initial values of the faceToNodeMap using emplaceBack
  forAll< parallelHostPolicy >( faceToNodeMapNew.size(),
                                [faceToNodeMapNew = faceToNodeMapNew.toView()]
                                  ( localIndex const faceIndex )
  {
    for( localIndex i = 0; i < faceToNodeMapNew.capacityOfArray( faceIndex ); ++i )
    {
      faceToNodeMapNew.emplaceBack( faceIndex, -1 );
    }
  } );

  // add the number of non-edge face nodes
  localIndex numInternalFaceNodes = 0;
  localIndex const numNonEdgeNodesPerFace = pow( order-1, 2 );
  for( localIndex kf=0; kf<m_faceManager->size(); ++kf )
  {
    localIndex const numEdgesPerFace = facesToEdges.sizeOfArray( kf );
    if( numEdgesPerFace==4 )
      numInternalFaceNodes += numNonEdgeNodesPerFace;
    else
    {
      GEOSX_ERROR( "need more support for face geometry" );
    }
  }

  // add the number of non-face element nodes
  localIndex numInternalElementNodes = 0;
  source.m_elementManager->forElementRegions< CellElementRegion >( [&]( CellElementRegion const & sourceRegion )
  {
    sourceRegion.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion const & sourceSubRegion )
    {
      if( sourceSubRegion.getElementType() == ElementType::Hexahedron )
      {
        numInternalElementNodes += sourceSubRegion.size() * pow( order-1, 3 );
      }

    } );
  } );



  localIndex const numNodes = source.m_nodeManager->size()
                              + numInternalEdgeNodes
                              + numInternalFaceNodes
                              + numInternalElementNodes;

  m_nodeManager->resize( numNodes );

  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const refPosSource = source.m_nodeManager->referencePosition();
  arrayView2d< real64, nodes::REFERENCE_POSITION_USD > const refPosNew = m_nodeManager->referencePosition().toView();

  {
    Group & nodeSets = m_nodeManager->sets();
    SortedArray< localIndex > & allNodes  = nodeSets.registerWrapper< SortedArray< localIndex > >( string( "all" ) ).reference();
    allNodes.reserve( m_nodeManager->size() );

    for( localIndex a=0; a<m_nodeManager->size(); ++a )
    {
      allNodes.insert( a );
    }

  }


  // ---- to calculate the node local to global map ----
  cellBlockManager.setNumNodes( numNodes );

  arrayView1d< globalIndex > nodelocalToGlobal = cellBlockManager.getNodeLocalToGlobal();

  CellBlock & currentcellblock = cellBlockManager.getCellBlocks().getGroup< CellBlock >( 0 );

  array1d< localIndex > partitionInfo = currentcellblock.getPartitionInformation();

  //GEOSX_LOG_RANK("!!!! INFO !!!! MeshLevel::MeshLevel partitionInfo="<<partitionInfo);

  // the number of dimensions for the current elements
  localIndex m_dim = 3;
  // the number of nodes in each direction for the current partition
  array1d< localIndex > numNodesInDir{3}; //= { 1, 1, 1 };
  array1d< localIndex > numElemsInDir{3}; // = { 1, 1, 1 };

  localIndex localNodeIndex= 0;
  // the partitionInfo are passing from InternalMeshGenerator.cpp:
  // partitionInfo[0-2] are firstElemIndexInPartition
  // partitionInfo[3-5] are lastElemIndexInPartition
  for( localIndex i = 0; i < m_dim; ++i )
  {
    numNodesInDir[i] = (partitionInfo[i+3] - partitionInfo[i]) * 3 + 4;
    numElemsInDir[i] = (partitionInfo[i+3] - partitionInfo[i]) + 1;
  }

  // numNodesInDir[0] is the number of nodes on the x direction of the current partition
  for( localIndex i = 0; i < numNodesInDir[0]; ++i )
  {
    // numNodesInDir[0] is the number of nodes on the y direction of the current partition
    for( localIndex j = 0; j < numNodesInDir[1]; ++j )
    {
      // numNodesInDir[0] is the number of nodes on the z direction of the current partition
      for( localIndex k = 0; k < numNodesInDir[2]; ++k )
      {
        // globalnodeIJK is initilized as the global node index of the current partition
        // plus the global information of the partition which is firstElemIndexInPartition
        localIndex globalnodeIJK[3] = { i + partitionInfo[0] * 3, j + partitionInfo[1] * 3, k + partitionInfo[2] *3 };

        //               15__________________ 63
        //               /|                  /|
        //              / |                 / |
        //             /  |                /  |
        //           3/___|_______________/51 |
        //            |   |               |   |
        //            |   |               |   |
        //            |   |               |   |
        //            2   |               |   |
        //            |   |               |   |
        //            |   12-----28-----44|---/16       z
        //            1  /                |  /          |   y
        //            | /                 | /           |  /
        //            |/__________________|/            | /
        //            0      16     32    48            |/____ x

        // assigning the localNodeIndex and globalnodeIJK for nodes
        nodelocalToGlobal[localNodeIndex] = globalnodeIJK[0]*(partitionInfo[7]*3+1)*(partitionInfo[8]*3+1)
                                            + globalnodeIJK[1]*(partitionInfo[8]*3+1)
                                            + globalnodeIJK[2];

        // local node index is simply incremental
        ++localNodeIndex;
      }
    }
  }


  ArrayOfArraysView< localIndex const > const & faceToNodeMapSource = source.m_faceManager->nodeList().toViewConst();

  FaceManager::ElemMapType const & faceToElem = source.m_faceManager->toElementRelation();

  arrayView2d< localIndex const > const & faceToElemIndex = faceToElem.m_toElementIndex.toViewConst();


  source.m_elementManager->forElementRegions< CellElementRegion >( [&]( CellElementRegion const & sourceRegion )
  {
    CellElementRegion & region = *(dynamic_cast< CellElementRegion * >( m_elementManager->createChild( sourceRegion.getCatalogName(),
                                                                                                       sourceRegion.getName() ) ) );
    region.addCellBlockNames( sourceRegion.getCellBlockNames() );

    sourceRegion.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion const & sourceSubRegion )
    {

      localIndex const numNodesPerElem = pow( order+1, 3 );

      CellElementSubRegion & newSubRegion = region.getSubRegions().registerGroup< CellElementSubRegion >( sourceSubRegion.getName() );

      // assigning the element to nodes map through cell block
      arrayView2d< localIndex, cells::NODE_MAP_USD > elemsToNodes = currentcellblock.getElemToNode();

      localIndex localElemIndex= 0;

      // numElemsInDir[0] is the number of elements on the x direction of the current partition
      for( localIndex i = 0; i < numElemsInDir [0]; ++i )
      {
        // numElemsInDir[1] is the number of elements on the y direction of the current partition
        for( localIndex j = 0; j < numElemsInDir[1]; ++j )
        {
          // numElemsInDir[2] is the number of elements on the z direction of the current partition
          for( localIndex k = 0; k < numElemsInDir[2]; ++k )
          {
            // get the first node index in the current partition
            localIndex firstNodeIndex = i * numNodesInDir[2] * numNodesInDir[1] * 3 + j * numNodesInDir[2] * 3 + k * 3;
            //GEOSX_LOG_RANK("!!!! INFO !!!! CellBlockManager::setO3Information firstNodeIndex=" << firstNodeIndex);

            // assign the element to node index with in the current partition
            for( localIndex zNodeIndex = 0; zNodeIndex < 4; ++zNodeIndex )
            {
              for( localIndex yNodeIndex = 0; yNodeIndex < 4; ++yNodeIndex )
              {
                for( localIndex xNodeIndex = 0; xNodeIndex < 4; ++xNodeIndex )
                {

                  //               60__________________ 63
                  //               /                   /|
                  //              /                   / |
                  //             /                   /  |
                  //          48/__________________ /51 |
                  //            |                   |   |
                  //            |                   |   |
                  //            |                   |   |
                  //            |                   |   |
                  //            |                   |   |
                  //            |                   |   /15       z
                  //            |                   |  /          |   y
                  //            |                   | /           |  /
                  //            |___________________|/            | /
                  //            0     1      2      3             |/____ x

                  elemsToNodes[localElemIndex][xNodeIndex + yNodeIndex * 4 + zNodeIndex * 16 ] =
                    // note that the indexing of elemsToNodes are in the x-y-z ordering
                    // and it needs to consider the postion of each element in the current partition
                    numNodesInDir[1] * numNodesInDir[2] * xNodeIndex + numNodesInDir[2] * yNodeIndex + zNodeIndex + firstNodeIndex;
                }
              }
            }

            // local node index is simply incremental
            ++localElemIndex;
          }
        }
      }

      newSubRegion.setElementType( sourceSubRegion.getElementType() );

      newSubRegion.resizePerElementValues( numNodesPerElem,
                                           sourceSubRegion.numEdgesPerElement(),
                                           sourceSubRegion.numFacesPerElement() );

      newSubRegion.resize( sourceSubRegion.size() );

      arrayView2d< localIndex const > const & elemToFaces = sourceSubRegion.faceList();

      arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodesSource = sourceSubRegion.nodeList().toViewConst();
      array2d< localIndex, cells::NODE_MAP_PERMUTATION > & elemsToNodesNew = newSubRegion.nodeList();
      elemsToNodesNew =  elemsToNodes;

      array2d< localIndex > & elemToFacesNew = newSubRegion.faceList();
      elemToFacesNew = elemToFaces;

      // Copy a new elemCenter map from the old one
      array2d< real64 > & elemCenterNew = newSubRegion.getElementCenter();
      arrayView2d< real64 const > const elemCenterOld = sourceSubRegion.getElementCenter();
      elemCenterNew = elemCenterOld;
      for( localIndex elem = 0; elem < elemsToNodesNew.size( 0 ); ++elem )
      {
        for( localIndex a = 0; a < 3; ++a )
        {
          elemCenterNew[elem][a] = elemCenterOld[elem][a];
        }
      }

      // Fill a temporary table which knowing the global number of a degree of freedom and a face, gives you the local number of this degree
      // of freedom on the considering face
      array2d< localIndex > localElemToLocalFace( 6, numNodesPerElem );

      //Init arrays
      for( localIndex i = 0; i < 6; ++i )
      {
        for( localIndex j = 0; j < numNodesPerElem; ++j )
        {
          localElemToLocalFace[i][j]=-1;
        }

      }



      for( localIndex i = 0; i < faceToNodeMapSource.size(); ++i )
      {
        for( localIndex j = 0; j < numNodesPerFace; ++j )
        {
          faceToNodeMapNew[i][j] = -1;
        }

      }

      //Face 0
      for( localIndex i = 0; i < order+1; i++ )
      {
        for( localIndex j = 0; j <order+1; j++ )
        {
          localElemToLocalFace[0][i + numNodesPerFace * j] = i + (order+1)*j;
        }

      }

      //Face 1
      for( localIndex i = 0; i < order+1; i++ )
      {
        for( localIndex k = 0; k < order+1; k++ )
        {
          localElemToLocalFace[1][k+(order+1)*i] = k + (order+1)*i;
        }

      }

      //Face 2
      for( localIndex k = 0; k < order+1; k++ )
      {
        for( localIndex j = 0; j < order+1; j++ )
        {
          localElemToLocalFace[2][k * numNodesPerFace +j*(order+1)] = j + (order+1)*k;
        }

      }

      //Face 3
      for( localIndex j = 0; j < order+1; j++ )
      {
        for( localIndex k = 0; k < order+1; k++ )
        {
          localElemToLocalFace[3][order +k*(order+1)+j * numNodesPerFace] = k + (order+1)*j;
        }

      }

      //Face 4
      for( localIndex j = 0; j < order+1; j++ )
      {
        for( localIndex i = 0; i < order+1; i++ )
        {
          localElemToLocalFace[4][order*(order+1)+i+j * numNodesPerFace] = i + (order+1)*j;
        }

      }

      //Face 5
      for( localIndex k = 0; k < order+1; k++ )
      {
        for( localIndex i = 0; i < order+1; i++ )
        {
          localElemToLocalFace[5][order * numNodesPerFace +i+k*(order+1)] = i + (order+1)*k;
        }

      }

      //Initialisation of elemToNodes
      /*
         for( localIndex e = 0; e < elemsToNodesNew.size( 0 ); ++e )
         {
         for( localIndex i = 0; i < numNodesPerElem; i++ )
         {
          elemsToNodesNew[e][i] =-1;
         }

         }
       */

      //localIndex count=0;

      for( localIndex elem = 0; elem < elemsToNodesNew.size( 0 ); ++elem )
      {
        for( localIndex k = 0; k < order+1; ++k )
        {
          for( localIndex j = 0; j< order+1; ++j )
          {
            for( localIndex i = 0; i <order+1; ++i )
            {

              localIndex face = 0;
              localIndex foundFace = 0;
              localIndex nodeindex = i + (order+1)*j + (order+1)*(order+1)*k;

              while( face<6 && foundFace<1 )
              {
                localIndex m = localElemToLocalFace[face][nodeindex];

                if( m != -1 )
                {
                  if( faceToNodeMapNew[elemToFaces[elem][face]][m] != -1 )
                  {
                    foundFace = 1;
                    for( localIndex l = 0; l < 2; ++l )
                    {
                      localIndex elemNeighbour = faceToElemIndex[elemToFaces[elem][face]][l];
                      if( elemNeighbour != elem && elemNeighbour != -1 )
                      {
                        for( localIndex node = 0; node < numNodesPerElem; ++node )
                        {
                          if( elemsToNodesNew[elemNeighbour][node] == faceToNodeMapNew[elemToFaces[elem][face]][m] )
                          {
                            elemsToNodesNew[elem][nodeindex] = elemsToNodesNew[elemNeighbour][node];
                            break;
                          }
                        }
                        break;
                      }
                    }
                    for( localIndex face2 = 0; face2 < 6; ++face2 )
                    {
                      localIndex m2 = localElemToLocalFace[face2][nodeindex];

                      if( m2 != -1 && face2!=face )
                      {
                        faceToNodeMapNew[elemToFaces[elem][face2]][m2] = faceToNodeMapNew[elemToFaces[elem][face]][m];
                      }
                    }
                  }
                  else
                  {
                    faceToNodeMapNew[elemToFaces[elem][face]][m] = elemsToNodesNew[elem][nodeindex];
                    //faceToNodeMapNew[elemToFaces[elem][face]][m] = count;
                  }
                }
                face++;
              }
              /*
                 if( face > 5 && foundFace < 1 )
                 {
                 elemsToNodesNew[elem][nodeindex] = count;
                 count++;
                 }
               */
            }
          }
        }
      }

      //Fill a temporary array which contains the Gauss-Lobatto points depending on the order
      array1d< real64 > GaussLobattoPts( order+1 );

      switch( order )
      {
        case 1:
          GaussLobattoPts[0] = -1.0;
          GaussLobattoPts[1] = 1.0;
          break;
        case 2:
          GaussLobattoPts[0] = -1.0;
          GaussLobattoPts[1] = 0.0;
          GaussLobattoPts[2] = 1.0;
          break;
        case 3:
          static constexpr real64 sqrt5 = 2.2360679774997897;
          GaussLobattoPts[0] = -1.0;
          GaussLobattoPts[1] = -1./sqrt5;
          GaussLobattoPts[2] = 1./sqrt5;
          GaussLobattoPts[3] = 1.;
          break;
        case 4:
          static constexpr real64 sqrt3_7 = 0.6546536707079771;
          GaussLobattoPts[0] = -1.0;
          GaussLobattoPts[1] = -sqrt3_7;
          GaussLobattoPts[2] = 0.0;
          GaussLobattoPts[3] = sqrt3_7;
          GaussLobattoPts[4] = 1.0;
          break;
        case 5:
          static constexpr real64 sqrt__7_plus_2sqrt7__ = 3.50592393273573196;
          static constexpr real64 sqrt__7_mins_2sqrt7__ = 1.30709501485960033;
          static constexpr real64 sqrt_inv21 = 0.218217890235992381;
          GaussLobattoPts[0] = -1.0;
          GaussLobattoPts[1] = -sqrt_inv21*sqrt__7_plus_2sqrt7__;
          GaussLobattoPts[2] = -sqrt_inv21*sqrt__7_mins_2sqrt7__;
          GaussLobattoPts[3] = sqrt_inv21*sqrt__7_mins_2sqrt7__;
          GaussLobattoPts[4] = sqrt_inv21*sqrt__7_plus_2sqrt7__;
          GaussLobattoPts[5] = 1.0;
          break;
      }

      //Three 1D arrays to contains the GL points in the new coordinates knowing the mesh nodes
      array1d< real64 > x( order+1 );
      array1d< real64 > y( order+1 );
      array1d< real64 > z( order+1 );

      for( localIndex e = 0; e < elemsToNodesNew.size( 0 ); e++ )
      {
        //Fill the three 1D array
        for( localIndex i = 0; i < order+1; i++ )
        {
          x[i] = refPosSource[elemsToNodesSource[e][0]][0] + ((refPosSource[elemsToNodesSource[e][1]][0]-refPosSource[elemsToNodesSource[e][0]][0])/2.)*GaussLobattoPts[i]
                 + (refPosSource[elemsToNodesSource[e][1]][0]-refPosSource[elemsToNodesSource[e][0]][0])/2;
          y[i] = refPosSource[elemsToNodesSource[e][0]][1] + ((refPosSource[elemsToNodesSource[e][2]][1]-refPosSource[elemsToNodesSource[e][0]][1])/2.)*GaussLobattoPts[i]
                 + (refPosSource[elemsToNodesSource[e][2]][1]-refPosSource[elemsToNodesSource[e][0]][1])/2;
          z[i] = refPosSource[elemsToNodesSource[e][0]][2] + ((refPosSource[elemsToNodesSource[e][4]][2]-refPosSource[elemsToNodesSource[e][0]][2])/2.)*GaussLobattoPts[i]
                 + (refPosSource[elemsToNodesSource[e][4]][2]-refPosSource[elemsToNodesSource[e][0]][2])/2;
        }


        //Fill new refPos array
        for( localIndex k = 0; k< order+1; k++ )
        {
          for( localIndex j = 0; j < order+1; j++ )
          {
            for( localIndex i = 0; i < order+1; i++ )
            {
              localIndex const nodeIndex = elemsToNodesNew( e, i+j*(order+1)+k* numNodesPerFace );

              refPosNew( nodeIndex, 0 ) = x[i];
              refPosNew( nodeIndex, 1 ) = y[j];
              refPosNew( nodeIndex, 2 ) = z[k];

            }

          }

        }

      }
    } );
  } );
}



MeshLevel::~MeshLevel()
{
  if( !m_isShallowCopy )
  {
    delete m_nodeManager;
    delete m_edgeManager;
    delete m_faceManager;
    delete m_elementManager;
    delete m_embSurfNodeManager;
    delete m_embSurfEdgeManager;
  }

}

void MeshLevel::initializePostInitialConditionsPostSubGroups()
{
  m_elementManager->forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion & subRegion )
  {
    subRegion.calculateElementGeometricQuantities( *m_nodeManager, *m_faceManager );
  } );
}

void MeshLevel::generateAdjacencyLists( arrayView1d< localIndex const > const & seedNodeList,
                                        localIndex_array & nodeAdjacencyList,
                                        localIndex_array & edgeAdjacencyList,
                                        localIndex_array & faceAdjacencyList,
                                        ElementRegionManager::ElementViewAccessor< ReferenceWrapper< localIndex_array > > & elementAdjacencyList,
                                        integer const depth )
{
  NodeManager const & nodeManager = getNodeManager();

  ArrayOfArraysView< localIndex const > const & nodeToElementRegionList = nodeManager.elementRegionList().toViewConst();
  ArrayOfArraysView< localIndex const > const & nodeToElementSubRegionList = nodeManager.elementSubRegionList().toViewConst();
  ArrayOfArraysView< localIndex const > const & nodeToElementList = nodeManager.elementList().toViewConst();

  FaceManager const & faceManager = this->getFaceManager();
  ArrayOfArraysView< localIndex const > const & faceToEdges = faceManager.edgeList().toViewConst();

  ElementRegionManager const & elemManager = this->getElemManager();

  std::set< localIndex > nodeAdjacencySet;
  std::set< localIndex > edgeAdjacencySet;
  std::set< localIndex > faceAdjacencySet;
  std::vector< std::vector< std::set< localIndex > > > elementAdjacencySet( elemManager.numRegions() );

  for( localIndex a=0; a<elemManager.numRegions(); ++a )
  {
    elementAdjacencySet[a].resize( elemManager.getRegion( a ).numSubRegions() );
  }

  nodeAdjacencySet.insert( seedNodeList.begin(), seedNodeList.end() );

  for( integer d=0; d<depth; ++d )
  {
    for( localIndex const nodeIndex : nodeAdjacencySet )
    {
      for( localIndex b=0; b<nodeToElementRegionList.sizeOfArray( nodeIndex ); ++b )
      {
        localIndex const regionIndex = nodeToElementRegionList[nodeIndex][b];
        localIndex const subRegionIndex = nodeToElementSubRegionList[nodeIndex][b];
        localIndex const elementIndex = nodeToElementList[nodeIndex][b];
        elementAdjacencySet[regionIndex][subRegionIndex].insert( elementIndex );
      }
    }

    for( typename dataRepository::indexType kReg=0; kReg<elemManager.numRegions(); ++kReg )
    {
      ElementRegionBase const & elemRegion = elemManager.getRegion( kReg );

      elemRegion.forElementSubRegionsIndex< CellElementSubRegion,
                                            WellElementSubRegion >( [&]( localIndex const kSubReg,
                                                                         auto const & subRegion )
      {
        arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodes = subRegion.nodeList();
        arrayView2d< localIndex const > const elemsToFaces = subRegion.faceList();
        for( auto const elementIndex : elementAdjacencySet[kReg][kSubReg] )
        {
          for( localIndex a=0; a<elemsToNodes.size( 1 ); ++a )
          {
            nodeAdjacencySet.insert( elemsToNodes[elementIndex][a] );
          }

          for( localIndex a=0; a<elemsToFaces.size( 1 ); ++a )
          {
            faceAdjacencySet.insert( elemsToFaces[elementIndex][a] );

            localIndex const faceIndex = elemsToFaces[elementIndex][a];
            localIndex const numEdges = faceToEdges.sizeOfArray( faceIndex );
            for( localIndex b=0; b<numEdges; ++b )
            {
              edgeAdjacencySet.insert( faceToEdges( faceIndex, b ));
            }

          }

        }
      } );
    }
  }

  nodeAdjacencyList.resize( LvArray::integerConversion< localIndex >( nodeAdjacencySet.size()));
  std::copy( nodeAdjacencySet.begin(), nodeAdjacencySet.end(), nodeAdjacencyList.begin() );

  edgeAdjacencyList.resize( LvArray::integerConversion< localIndex >( edgeAdjacencySet.size()));
  std::copy( edgeAdjacencySet.begin(), edgeAdjacencySet.end(), edgeAdjacencyList.begin() );

  faceAdjacencyList.resize( LvArray::integerConversion< localIndex >( faceAdjacencySet.size()));
  std::copy( faceAdjacencySet.begin(), faceAdjacencySet.end(), faceAdjacencyList.begin() );

  for( localIndex kReg=0; kReg<elemManager.numRegions(); ++kReg )
  {
    ElementRegionBase const & elemRegion = elemManager.getRegion( kReg );

    for( localIndex kSubReg = 0; kSubReg < elemRegion.numSubRegions(); ++kSubReg )
    {
      elementAdjacencyList[kReg][kSubReg].get().clear();
      elementAdjacencyList[kReg][kSubReg].get().resize( LvArray::integerConversion< localIndex >( elementAdjacencySet[kReg][kSubReg].size()) );
      std::copy( elementAdjacencySet[kReg][kSubReg].begin(),
                 elementAdjacencySet[kReg][kSubReg].end(),
                 elementAdjacencyList[kReg][kSubReg].get().begin() );

    }
  }

}

void MeshLevel::generateSets()
{
  GEOSX_MARK_FUNCTION;

  NodeManager const & nodeManager = *m_nodeManager;

  dataRepository::Group const & nodeSets = nodeManager.sets();

  map< string, array1d< bool > > nodeInSet; // map to contain indicator of whether a node is in a set.
  string_array setNames; // just a holder for the names of the sets

  // loop over all wrappers and fill the nodeIndSet arrays for each set
  for( auto & wrapper : nodeSets.wrappers() )
  {
    string const & name = wrapper.second->getName();
    nodeInSet[name].resize( nodeManager.size() );
    nodeInSet[name].setValues< serialPolicy >( false );

    if( nodeSets.hasWrapper( name ) )
    {
      setNames.emplace_back( name );
      SortedArrayView< localIndex const > const & set = nodeSets.getReference< SortedArray< localIndex > >( name );
      for( localIndex const a : set )
      {
        nodeInSet[name][a] = true;
      }
    }
  }


  ElementRegionManager & elementRegionManager = *m_elementManager;
  elementRegionManager.forElementSubRegions( [&]( auto & subRegion )
  {
    dataRepository::Group & elementSets = subRegion.sets();

    auto const & elemToNodeMap = subRegion.nodeList();

    for( string const & setName : setNames )
    {
      arrayView1d< bool const > const nodeInCurSet = nodeInSet[setName];

      SortedArray< localIndex > & targetSet = elementSets.registerWrapper< SortedArray< localIndex > >( setName ).reference();
      for( localIndex k = 0; k < subRegion.size(); ++k )
      {
        localIndex const numNodes = subRegion.numNodesPerElement( k );

        localIndex elementInSet = true;
        for( localIndex i = 0; i < numNodes; ++i )
        {
          if( !nodeInCurSet( elemToNodeMap[ k ][ i ] ) )
          {
            elementInSet = false;
            break;
          }
        }

        if( elementInSet )
        {
          targetSet.insert( k );
        }
      }
    }
  } );
}

MeshLevel const & MeshLevel::getShallowParent() const
{
  return ( m_shallowParent!=nullptr ? *m_shallowParent : *this);
}

MeshLevel & MeshLevel::getShallowParent()
{
  return ( m_shallowParent!=nullptr ? *m_shallowParent : *this);
}


bool MeshLevel::isShallowCopyOf( MeshLevel const & comparisonLevel ) const
{
  return ( m_nodeManager == comparisonLevel.m_nodeManager ) &&
         ( m_edgeManager == comparisonLevel.m_edgeManager ) &&
         ( m_faceManager == comparisonLevel.m_faceManager ) &&
         ( m_elementManager == comparisonLevel.m_elementManager ) &&
         ( m_embSurfNodeManager == comparisonLevel.m_embSurfNodeManager) &&
         ( m_embSurfEdgeManager == comparisonLevel.m_embSurfEdgeManager ) &&
         isShallowCopy();
}


} /* namespace geosx */
