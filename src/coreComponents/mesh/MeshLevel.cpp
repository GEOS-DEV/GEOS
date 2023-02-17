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

  ////////////////////////////////
  // Get the new number of nodes 
  ////////////////////////////////

  // the total number of nodes: to add the number of non-vertex edge nodes
  localIndex const numEdgesPartition = source.m_edgeManager->size();
  localIndex const numBasisSupportPoints = order+1;
  localIndex const numNodesPerEdge = numBasisSupportPoints;
  localIndex const numNonVertexNodesPerEdge = numNodesPerEdge - 2; 
  localIndex const numInternalEdgeNodes = numEdgesPartition * numNonVertexNodesPerEdge;

  // the total number of nodes: to add the number of non-edge face nodes
  localIndex const numFacesPartition = source.m_faceManager->size();
  localIndex const numNodesPerFace = pow( order+1, 2 ); 
  localIndex const numNonEdgeNodesPerFace = pow( order-1, 2 ); 
  localIndex numInternalFaceNodes = 0; 
  localIndex numEdgesPerFace;

  ArrayOfArraysView< localIndex const > const & facesToEdges = source.m_faceManager->edgeList().toViewConst();
  for( localIndex kf=0; kf<numFacesPartition; ++kf )
  {
    numEdgesPerFace = facesToEdges.sizeOfArray( kf );
    if( numEdgesPerFace==4 )
      numInternalFaceNodes += numNonEdgeNodesPerFace;
    else 
    {    
      GEOSX_ERROR( "need more support for face geometry" );
    }    
  }

  // the total number of nodes: to add the number of non-face element nodes
  localIndex const numInternalNodesPerElement = ( order - 1 ) * ( order - 1 ) * ( order - 1 ); 
  localIndex numInternalElementNodes = 0; 

  source.m_elementManager->forElementRegions< CellElementRegion >( [&]( CellElementRegion const & sourceRegion )
  {
    sourceRegion.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion const & sourceSubRegion )
    {    
      if( sourceSubRegion.getElementType() == ElementType::Hexahedron )
      {    
        numInternalElementNodes += sourceSubRegion.size() * numInternalNodesPerElement;
      }    
    } ); 
  } ); 

  // calculate the total number of nodes for the current partition
  localIndex const numNodesPartitionSource = source.m_nodeManager->size();
  localIndex const numNodesPartitionNew = numNodesPartitionSource
                              + numInternalEdgeNodes
                              + numInternalFaceNodes
                              + numInternalElementNodes;


  //
  // ---- to initialize the node local to global map ----
  //   

  cellBlockManager.setNumNodes( numNodesPartitionNew, order );

  arrayView1d< globalIndex const > const nodeLocalToGlobalSource = source.m_nodeManager->localToGlobalMap();
  arrayView1d< globalIndex > nodeLocalToGlobalNew = cellBlockManager.getNodeLocalToGlobal();



  //////////////////////////
  // Edges
  //////////////////////////

  // the total number of nodes: to add the number of non-vertex edge nodes
  m_edgeManager->resize( numEdgesPartition );
  m_edgeManager->getDomainBoundaryIndicator() = source.m_edgeManager->getDomainBoundaryIndicator();

  arrayView1d< globalIndex const > const edgeLocalToGlobal = m_edgeManager->localToGlobalMap();
  for (localIndex iter_localToGlobalsize = 0; iter_localToGlobalsize < numEdgesPartition; iter_localToGlobalsize++)
  {
    m_edgeManager->localToGlobalMap()[iter_localToGlobalsize] = source.m_edgeManager->localToGlobalMap()[iter_localToGlobalsize];
  }
  m_edgeManager->constructGlobalToLocalMap();

  //
  // ---- initialize edge-to-node map ----
  //

  // get information from the source (base mesh-level) edge-to-node map
  arrayView2d< localIndex const > const & edgeToNodeMapSource = source.m_edgeManager->nodeList();
  arrayView2d< localIndex > edgeToNodeMapNew = cellBlockManager.getEdgeToNodes();
  m_edgeManager->nodeList().resize( numEdgesPartition, numNodesPerEdge );


  /////////////////////////
  // Faces
  //////////////////////////

  m_faceManager->resize( numFacesPartition );
  m_faceManager->faceCenter() = source.m_faceManager->faceCenter();
  m_faceManager->faceNormal() = source.m_faceManager->faceNormal();

  // copy the faces-to-edgs map from source
  m_faceManager->edgeList() = source.m_faceManager->edgeList();
  // copy the faces-to-elements map from source
  m_faceManager->elementList() = source.m_faceManager->elementList();
  // copy the faces-boundaryIndicator from source
  m_faceManager->getDomainBoundaryIndicator() = source.m_faceManager->getDomainBoundaryIndicator();

  arrayView1d< globalIndex const > const faceLocalToGlobal = m_faceManager->localToGlobalMap();
  for (localIndex iter_localToGlobalsize = 0; iter_localToGlobalsize < numFacesPartition; iter_localToGlobalsize++)
  {
    m_faceManager->localToGlobalMap()[iter_localToGlobalsize] = source.m_faceManager->localToGlobalMap()[iter_localToGlobalsize];
  }
  m_faceManager->constructGlobalToLocalMap();

  ArrayOfArraysView< localIndex const > const & faceToNodeMapSource = source.m_faceManager->nodeList().toViewConst();
  ArrayOfArrays< localIndex > & faceToNodeMapNew = cellBlockManager.getFaceToNodes();

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


  //////////////////////////
  // Nodes
  //////////////////////////

  // the number of nodes per elements depends on the order number
  localIndex const numNodesPerElem = pow( order+1, 3 );
  m_nodeManager->resize( numNodesPartitionNew );

  Group & nodeSets = m_nodeManager->sets();
  SortedArray< localIndex > & allNodes  = nodeSets.registerWrapper< SortedArray< localIndex > >( string( "all" ) ).reference();
  allNodes.reserve( numNodesPartitionNew );

  for( localIndex iter_nodes=0; iter_nodes< numNodesPartitionNew; ++iter_nodes )
  {
    allNodes.insert( iter_nodes );
  }

  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const refPosSource = source.m_nodeManager->referencePosition();
  arrayView2d< real64, nodes::REFERENCE_POSITION_USD > refPosNew = cellBlockManager.getNodePositions();
  refPosNew.setValues< parallelHostPolicy >( -1.0 );


  /////////////////////////
  // Elements 
  //////////////////////////

  source.m_elementManager->forElementRegions< CellElementRegion >( [&]( CellElementRegion const & sourceRegion )
  {
    // create element region with the same name as source element region "Region"
    CellElementRegion & region = *(dynamic_cast< CellElementRegion * >( m_elementManager->createChild( sourceRegion.getCatalogName(),
                                                                                                       sourceRegion.getName() ) ) );
    // add cell block to the new element region with the same name as cell block name from source element region
    region.addCellBlockNames( sourceRegion.getCellBlockNames() );

    sourceRegion.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion const & sourceSubRegion )
    {
      // create element sub region with the same name as source element sub region "cb"
      CellElementSubRegion & newSubRegion = region.getSubRegions().registerGroup< CellElementSubRegion >( sourceSubRegion.getName() );
      newSubRegion.setElementType( sourceSubRegion.getElementType() );

      // resize per elements value for the new sub region with the new number of nodes per element
      newSubRegion.resizePerElementValues( numNodesPerElem,
                                           sourceSubRegion.numEdgesPerElement(),
                                           sourceSubRegion.numFacesPerElement() );

      newSubRegion.resize( sourceSubRegion.size() );

      // copy new elemCenter map from source
      array2d< real64 > & elemCenterNew = newSubRegion.getElementCenter();
      arrayView2d< real64 const > const elemCenterOld = sourceSubRegion.getElementCenter();
      elemCenterNew = elemCenterOld;

      // copy the elements-to-faces map from source
      arrayView2d< localIndex const > const & elemsToFacesSource = sourceSubRegion.faceList();
      array2d< localIndex > & elemsToFacesNew = newSubRegion.faceList();
      elemsToFacesNew = elemsToFacesSource;

      // copy the elements-to-edgs map from source
      arrayView2d< localIndex const > const & elemsToEdgesSource = sourceSubRegion.edgeList();
      array2d< localIndex > & elemsToEdgesNew = newSubRegion.edgeList();
      elemsToEdgesNew = elemsToEdgesSource;


      // initialize the elements-to-nodes map 
      arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodesSource = sourceSubRegion.nodeList().toViewConst();
      arrayView2d< localIndex, cells::NODE_MAP_USD > elemsToNodesNew = cellBlockManager.getElemToNodes(newSubRegion.getName());

      arrayView1d< globalIndex const > const elementLocalToGlobal = sourceSubRegion.localToGlobalMap();
      for (localIndex iter_localToGlobalsize = 0; iter_localToGlobalsize < sourceSubRegion.localToGlobalMap().size(); iter_localToGlobalsize++)
      {
        newSubRegion.localToGlobalMap()[iter_localToGlobalsize] = sourceSubRegion.localToGlobalMap()[iter_localToGlobalsize];
      }
      newSubRegion.constructGlobalToLocalMap();

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

      array1d< localIndex > foundNodeLocalToGlobal ( numNodesPartitionNew );
      foundNodeLocalToGlobal.setValues< parallelHostPolicy > ( -1 );

      globalIndex const firstGlobalNodeIDforEdges = source.m_nodeManager->maxGlobalIndex() + 1;
      globalIndex const firstGlobalNodeIDforFaces = firstGlobalNodeIDforEdges + numNonVertexNodesPerEdge * ( source.m_edgeManager->maxGlobalIndex() + 1 );
      globalIndex const firstGlobalNodeIDforElems = firstGlobalNodeIDforFaces + numNonEdgeNodesPerFace * ( source.m_faceManager->maxGlobalIndex() + 1 );

      localIndex const numOfElements = elemsToNodesSource.size( 0 );
      localIndex const numOfNodesPerElementForBase = elemsToNodesSource.size(1);
      localIndex const numOfEdgesPerElement = elemsToEdgesSource.size(1);
      localIndex const numOfFacesPerElement = elemsToFacesSource.size(1);

      //Three 1D arrays to contains the GL points in the new coordinates knowing the mesh nodes
      array1d< real64 > corrdX( numBasisSupportPoints );
      array1d< real64 > corrdY( numBasisSupportPoints );
      array1d< real64 > corrdZ( numBasisSupportPoints );

      real64 tempPositionX;
      real64 tempPositionY;
      real64 tempPositionZ;

      localIndex localElemNodeID = 0;
      localIndex iter_sourcenodes = 0;
      localIndex iter_nodes = 0;

      for( localIndex iter_elem = 0; iter_elem < numOfElements; ++iter_elem )
      {
        //Fill the three 1D array
        for( localIndex iter_nodeDim = 0; iter_nodeDim < numBasisSupportPoints; iter_nodeDim++ )
        {
          corrdX[iter_nodeDim] = refPosSource[elemsToNodesSource[iter_elem][0]][0] +
                 ((refPosSource[elemsToNodesSource[iter_elem][1]][0]-refPosSource[elemsToNodesSource[iter_elem][0]][0])/2.)*GaussLobattoPts[iter_nodeDim]
                 + (refPosSource[elemsToNodesSource[iter_elem][1]][0]-refPosSource[elemsToNodesSource[iter_elem][0]][0])/2;

          corrdY[iter_nodeDim] = refPosSource[elemsToNodesSource[iter_elem][0]][1] +
                 ((refPosSource[elemsToNodesSource[iter_elem][2]][1]-refPosSource[elemsToNodesSource[iter_elem][0]][1])/2.)*GaussLobattoPts[iter_nodeDim]
                 + (refPosSource[elemsToNodesSource[iter_elem][2]][1]-refPosSource[elemsToNodesSource[iter_elem][0]][1])/2;

          corrdZ[iter_nodeDim] = refPosSource[elemsToNodesSource[iter_elem][0]][2] +
                 ((refPosSource[elemsToNodesSource[iter_elem][4]][2]-refPosSource[elemsToNodesSource[iter_elem][0]][2])/2.)*GaussLobattoPts[iter_nodeDim]
                 + (refPosSource[elemsToNodesSource[iter_elem][4]][2]-refPosSource[elemsToNodesSource[iter_elem][0]][2])/2;
        }

        for ( localIndex iter_elemNodeZ = 0; iter_elemNodeZ < numBasisSupportPoints;  ++iter_elemNodeZ )
        {
          for ( localIndex iter_elemNodeY = 0; iter_elemNodeY < numBasisSupportPoints;  ++iter_elemNodeY )
          {
            for ( localIndex iter_elemNodeX = 0; iter_elemNodeX < numBasisSupportPoints;  ++iter_elemNodeX )
            {
              localIndex iter_elemNode = iter_elemNodeX + iter_elemNodeY * numBasisSupportPoints
                                       + iter_elemNodeZ * numBasisSupportPoints * numBasisSupportPoints;

              tempPositionX = corrdX[iter_elemNodeX];
              tempPositionY = corrdY[iter_elemNodeY];
              tempPositionZ = corrdZ[iter_elemNodeZ];
              iter_nodes = 0;
              while (  iter_nodes < numNodesPartitionNew )
              {
                if ( (  std::fabs ( refPosNew[iter_nodes][0] - tempPositionX ) > 1e-5 ) ||
                     (  std::fabs ( refPosNew[iter_nodes][1] - tempPositionY ) > 1e-5 ) ||
                     (  std::fabs ( refPosNew[iter_nodes][2] - tempPositionZ ) > 1e-5 ) )
                {
                  iter_nodes++;
                }
                else
                {
                  elemsToNodesNew[ iter_elem ][ iter_elemNode ] = iter_nodes;
                  break;
                }
              }
              if ( iter_nodes == numNodesPartitionNew )
              {
                elemsToNodesNew[ iter_elem ][ iter_elemNode ] = localElemNodeID;
                refPosNew[localElemNodeID][0] = tempPositionX;
                refPosNew[localElemNodeID][1] = tempPositionY;
                refPosNew[localElemNodeID][2] = tempPositionZ;

                iter_sourcenodes = 0;
                while ( (  iter_sourcenodes < numNodesPartitionSource ) && ( foundNodeLocalToGlobal[ localElemNodeID ] < 0 ) )
                {
                  if ( (  std::fabs ( refPosSource[iter_sourcenodes][0] - tempPositionX ) < 1e-5 ) &&
                       (  std::fabs ( refPosSource[iter_sourcenodes][1] - tempPositionY ) < 1e-5 ) &&
                       (  std::fabs ( refPosSource[iter_sourcenodes][2] - tempPositionZ ) < 1e-5 ) )
                  {
                    nodeLocalToGlobalNew[ localElemNodeID ] = nodeLocalToGlobalSource[ iter_sourcenodes ];
                    foundNodeLocalToGlobal[ localElemNodeID ] = 1;
                    break;
                  }
                  else
                    iter_sourcenodes++;
                }
                localElemNodeID++;
              }
            }
          }
        }
      }

      for( localIndex iter_edge = 0; iter_edge < numEdgesPartition; ++iter_edge )
      {
        globalIndex edgeGlobalID = edgeLocalToGlobal[iter_edge];
        localIndex counterEdgeGlobalID = 0;

        for (localIndex iter_edgeNode = 0; iter_edgeNode < numBasisSupportPoints;  ++iter_edgeNode)
        {
          tempPositionX = refPosSource[ edgeToNodeMapSource[iter_edge][0]][0]
                                  + ((refPosSource[ edgeToNodeMapSource[iter_edge][1]][0]-refPosSource[ edgeToNodeMapSource[iter_edge][0]][0])/2.)
                                    * GaussLobattoPts[iter_edgeNode]
                                  + (refPosSource[ edgeToNodeMapSource[iter_edge][1]][0]-refPosSource[ edgeToNodeMapSource[iter_edge][0]][0])/2;

          tempPositionY = refPosSource[ edgeToNodeMapSource[iter_edge][0]][1]
                                  + ((refPosSource[ edgeToNodeMapSource[iter_edge][1]][1]-refPosSource[ edgeToNodeMapSource[iter_edge][0]][1])/2.)
                                    * GaussLobattoPts[iter_edgeNode]
                                  + (refPosSource[ edgeToNodeMapSource[iter_edge][1]][1]-refPosSource[ edgeToNodeMapSource[iter_edge][0]][1])/2;

          tempPositionZ = refPosSource[ edgeToNodeMapSource[iter_edge][0]][2]
                                  + ((refPosSource[ edgeToNodeMapSource[iter_edge][1]][2]-refPosSource[ edgeToNodeMapSource[iter_edge][0]][2])/2.)
                                    * GaussLobattoPts[iter_edgeNode]
                                  + (refPosSource[ edgeToNodeMapSource[iter_edge][1]][2]-refPosSource[ edgeToNodeMapSource[iter_edge][0]][2])/2;

          iter_nodes = 0;
          while (  iter_nodes < numNodesPartitionNew )
          {
            if ( (  std::fabs ( refPosNew[iter_nodes][0] - tempPositionX ) > 1e-5 ) ||
                 (  std::fabs ( refPosNew[iter_nodes][1] - tempPositionY ) > 1e-5 ) ||
                 (  std::fabs ( refPosNew[iter_nodes][2] - tempPositionZ ) > 1e-5 ) )
             {
               iter_nodes++;
             }
             else break;
          }
          if ( iter_nodes == numNodesPartitionNew )
          {
            GEOSX_LOG_RANK ("!!!! ERROR !!!! CREATE edgeToNodeMapNew: CANNOT FIND NODES");
          }
          else
          {
            edgeToNodeMapNew[iter_edge][ iter_edgeNode ] = iter_nodes;
            if ( foundNodeLocalToGlobal[ iter_nodes ] < 0 )
            {
               nodeLocalToGlobalNew[ iter_nodes ] = firstGlobalNodeIDforEdges + edgeGlobalID * numNonVertexNodesPerEdge +  counterEdgeGlobalID;
               foundNodeLocalToGlobal[ iter_nodes ] = 1;
               counterEdgeGlobalID++;
            }
          }
        }
      }

      for( localIndex iter_face = 0; iter_face < numFacesPartition; ++iter_face )
      {
        globalIndex faceGlobalID = faceLocalToGlobal[iter_face];
        localIndex counterFaceGlobalID = 0;

        real64 refPosXmax = std::max( std::max( refPosSource[faceToNodeMapSource[iter_face][0]][0], refPosSource[faceToNodeMapSource[iter_face][1]][0] ),
                                refPosSource[faceToNodeMapSource[iter_face][2]][0] );
        real64 refPosYmax = std::max( std::max( refPosSource[faceToNodeMapSource[iter_face][0]][1], refPosSource[faceToNodeMapSource[iter_face][1]][1] ),
                                refPosSource[faceToNodeMapSource[iter_face][2]][1] );
        real64 refPosZmax = std::max( std::max( refPosSource[faceToNodeMapSource[iter_face][0]][2], refPosSource[faceToNodeMapSource[iter_face][1]][2] ),
                                refPosSource[faceToNodeMapSource[iter_face][2]][2] );


        for (localIndex iter_faceNodeA = 0; iter_faceNodeA < numEdgesPerFace;  ++iter_faceNodeA)
        {
          for (localIndex iter_faceNodeB = 0; iter_faceNodeB < numEdgesPerFace;  ++iter_faceNodeB)
          {
            localIndex iter_faceNode = iter_faceNodeA + iter_faceNodeB * numEdgesPerFace;

            if ( std::fabs ( refPosXmax - refPosSource[faceToNodeMapSource[iter_face][0]][0] ) < 1e-5 )
            {
              tempPositionX = refPosSource[faceToNodeMapSource[iter_face][0]][0];
              tempPositionY = refPosSource[faceToNodeMapSource[iter_face][0]][1]
                              + ((refPosYmax-refPosSource[faceToNodeMapSource[iter_face][0]][1])/2.)*GaussLobattoPts[ iter_faceNodeA ]
                              + (refPosYmax-refPosSource[faceToNodeMapSource[iter_face][0]][1])/2;
              tempPositionZ = refPosSource[faceToNodeMapSource[iter_face][0]][2]
                              + ((refPosZmax-refPosSource[faceToNodeMapSource[iter_face][0]][2])/2.)*GaussLobattoPts[ iter_faceNodeB ]
                              + (refPosZmax-refPosSource[faceToNodeMapSource[iter_face][0]][2])/2;
            }

            else if ( std::fabs ( refPosYmax - refPosSource[faceToNodeMapSource[iter_face][0]][1] ) < 1e-5 )
            {
              tempPositionX = refPosSource[faceToNodeMapSource[iter_face][0]][0]
                              + ((refPosXmax-refPosSource[faceToNodeMapSource[iter_face][0]][0])/2.)*GaussLobattoPts[ iter_faceNodeA ]
                              + (refPosXmax-refPosSource[faceToNodeMapSource[iter_face][0]][0])/2;
              tempPositionY = refPosSource[faceToNodeMapSource[iter_face][0]][1];
              tempPositionZ = refPosSource[faceToNodeMapSource[iter_face][0]][2]
                              + ((refPosZmax-refPosSource[faceToNodeMapSource[iter_face][0]][2])/2.)*GaussLobattoPts[ iter_faceNodeB ]
                              + (refPosZmax-refPosSource[faceToNodeMapSource[iter_face][0]][2])/2;
            }
            else
          {
              tempPositionX = refPosSource[faceToNodeMapSource[iter_face][0]][0]
                              + ((refPosXmax-refPosSource[faceToNodeMapSource[iter_face][0]][0])/2.)*GaussLobattoPts[ iter_faceNodeA ]
                              + (refPosXmax-refPosSource[faceToNodeMapSource[iter_face][0]][0])/2;
              tempPositionY = refPosSource[faceToNodeMapSource[iter_face][0]][1]
                              + ((refPosYmax-refPosSource[faceToNodeMapSource[iter_face][0]][1])/2.)*GaussLobattoPts[ iter_faceNodeB ]
                              + (refPosYmax-refPosSource[faceToNodeMapSource[iter_face][0]][1])/2;
              tempPositionZ = refPosSource[faceToNodeMapSource[iter_face][0]][2];
            }

            iter_nodes = 0;
            while (  iter_nodes < numNodesPartitionNew )
            {
              if ( (  std::fabs ( refPosNew[iter_nodes][0] - tempPositionX ) > 1e-5 ) ||
                   (  std::fabs ( refPosNew[iter_nodes][1] - tempPositionY ) > 1e-5 ) ||
                   (  std::fabs ( refPosNew[iter_nodes][2] - tempPositionZ ) > 1e-5 ) )
                 iter_nodes++;
               else break;
            }
            if ( iter_nodes == numNodesPartitionNew )
            {
              GEOSX_LOG_RANK ("!!!! ERROR !!!! CREATE faceToNodeMapNew: CANNOT FIND NODES");
            }
            else
            {
              faceToNodeMapNew[iter_face][ iter_faceNode ] = iter_nodes;
              if ( foundNodeLocalToGlobal[ iter_nodes ] < 0 )
              {
                nodeLocalToGlobalNew[ iter_nodes ] = firstGlobalNodeIDforFaces + faceGlobalID * numNonEdgeNodesPerFace +  counterFaceGlobalID;
                foundNodeLocalToGlobal[ iter_nodes ] = 1;
                counterFaceGlobalID++;
              }
            }
          }
        }
      }
      for( localIndex iter_elem = 0; iter_elem < numOfElements; ++iter_elem )
      {
        globalIndex elemGlobalID = elementLocalToGlobal[iter_elem];
        localIndex counterElemGlobalID = 0;

        for( iter_nodes = 0; iter_nodes < numNodesPartitionNew; ++iter_nodes )
        {
          if ( foundNodeLocalToGlobal[ iter_nodes ] < 0 )
          {
            nodeLocalToGlobalNew[ iter_nodes ] = firstGlobalNodeIDforElems + elemGlobalID * numInternalNodesPerElement + counterElemGlobalID;
            foundNodeLocalToGlobal[ iter_nodes ] = 1;
            counterElemGlobalID++;
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
