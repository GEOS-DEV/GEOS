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
#include "common/MpiWrapper.hpp"
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
  MeshLevel(name, parent, 1)
{}

MeshLevel::MeshLevel( string const & name,
                      Group * const parent,
                      int order ):
  Group( name, parent ),
  m_nodeManager( new NodeManager( groupStructKeys::nodeManagerString(), this ) ),
  m_edgeManager( new EdgeManager( groupStructKeys::edgeManagerString(), this ) ),
  m_faceManager( new FaceManager( groupStructKeys::faceManagerString(), this ) ),
  m_elementManager( new ElementRegionManager( groupStructKeys::elemManagerString(), this ) ),
  m_embSurfNodeManager( new EmbeddedSurfaceNodeManager( groupStructKeys::embSurfNodeManagerString, this ) ),
  m_embSurfEdgeManager( new EdgeManager( groupStructKeys::embSurfEdgeManagerString, this ) ),
  m_modificationTimestamp( 0 ),
  m_isShallowCopy( false ),
  m_shallowParent( nullptr ),
  m_order( order )
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
  m_shallowParent( &source ),
  m_order( source.getOrder() )
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

template< typename T >
struct NodeKeyHasher {
  std::size_t operator()(const std::array< T, 6 >& arr) const {
    std::size_t hash = 0;
    // use a boost-style hash function
    for (auto v : arr) {
        hash ^= std::hash<T>{}( v )  + 0x9e3779b9 + ( hash << 6 ) + ( hash >> 2 );
    }
    return hash;
  }
};

template< typename T >
struct NodeKeyEqual {
  bool operator()(const std::array< T, 6 > & lhs, const std::array< T, 6 > & rhs) const
  {
     for( int i=0; i< 6; i++ )
     {
       if( lhs[ i ] != rhs[ i ] )
       {

         return false;
       }
     }
     return true;
  }
};

template< typename T >
static std::array< T, 6 > createNodeKey( T v )
{
  return std::array< T, 6 > { v, -1, -1 ,-1 ,-1, -1 };
}

template< typename T >
static std::array< T, 6 > createNodeKey( T v1, T v2, int a, int order )
{
  if( a == 0) return createNodeKey( v1 );
  if( a == order ) return createNodeKey( v2 );
  if( v1 < v2 )
  {
    return std::array< T, 6 > { v1, v2, -1 ,-1 ,a , -1 };
  }
  else
  {
    return std::array< T, 6 > { v2, v1, -1 ,-1 ,order - a , -1 };
  }
}

template< typename T >
static std::array< T, 6 > createNodeKey( T v1, T v2, T v3, T v4, int a, int b, int order )
{
  if( a == 0 ) return createNodeKey( v1, v3, b, order );
  if( a == order ) return createNodeKey( v2, v4, b, order );
  if( b == 0 ) return createNodeKey( v1, v2, a, order );
  if( b == order ) return createNodeKey( v3, v4, a, order );
  // arrange the vertices of the face such that v1 is the lowest value, and v2 is lower than v3
  // this ensures a coherent orientation of all face nodes
  while( v1 > v2 || v1 > v3 || v1 > v4 || v2 > v3 )
  {
    if( v1 > v2 )
    {
      std::swap( v1, v2 );
      std::swap( v3, v4 );
      a = order - a;
    }
    if( v1 > v3 )
    {
      std::swap( v1, v3 );
      std::swap( v2, v4 );
      b = order - b;
    }
    if( v1 > v4 )
    {
      std::swap( v1, v4 );
      std::swap( a, b );
      a = order - a;
      b = order - b;
    }
    if( v2 > v3 )
    {
      std::swap( v2, v3 );
      std::swap( a, b );
    }
  }
  return std::array< T, 6 > { v1, v2, v3, v4, a, b };
}

template< typename T >
static std::array< T, 6 > createNodeKey( T const (&elemNodes)[ 8 ], int q1, int q2, int q3, int order )
{
  bool extremal1 = q1 == 0 || q1 == order;
  bool extremal2 = q2 == 0 || q2 == order;
  bool extremal3 = q3 == 0 || q3 == order;
  int v1 = q1/order;
  int v2 = q2/order;
  int v3 = q3/order;
  if( extremal1 && extremal2 && extremal3 )
  {
    // vertex node
    return createNodeKey( elemNodes[ v1 + 2*v2 + 4*v3 ] );
  }
  else if( extremal1 && extremal2 )
  {
    // edge node on v1, v2
    return createNodeKey( elemNodes[ v1 + 2*v2 ], elemNodes[ v1 + 2*v2 + 4 ], q3, order );
  }
  else if( extremal1 && extremal3 )
  {
    // edge node on v1, v3
    return createNodeKey( elemNodes[ v1 + 4*v3 ], elemNodes[ v1 + 2 + 4*v3 ], q2, order );
  }
  else if( extremal2 && extremal3 )
  {
    // edge node on v2, v3
    return createNodeKey( elemNodes[ 2*v2 + 4*v3 ], elemNodes[ 1 + 2*v2 + 4*v3 ], q1, order );
  }
  else if( extremal1 )
  {
    // face node on the face of type 1
    return createNodeKey( elemNodes[ v1 ], elemNodes[ v1 + 2 ], elemNodes[ v1 + 4 ], elemNodes[ v1 + 2 + 4 ], q2, q3, order );
  }
  else if( extremal2 )
  {
    // face node on the face of type 2
    return createNodeKey( elemNodes[ 2*v2 ], elemNodes[ 1 + 2*v2 ], elemNodes[ 2*v2 + 4 ], elemNodes[ 1 + 2*v2 + 4 ], q1, q3, order );
  }
  else if( extremal3 )
  {
    // face node on the face of type 3
    return createNodeKey( elemNodes[ 4*v3 ], elemNodes[ 1 + 4*v3 ], elemNodes[ 2 + 4*v3 ], elemNodes[ 1 + 2 + 4*v3 ], q1, q2, order );
  }
  else
  {
   // node internal to the cell -- no need for key, it will be created
   return createNodeKey( -1 );
  }
}

static array1d< real64 > gaussLobattoPoints( int order )
{
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
  return GaussLobattoPts;
}

static void trilinearInterp( real64 const alpha,
                             real64 const beta,
                             real64 const gamma,
                             real64 const (&X)[8][3],
                             real64 (&coords)[3] )
{
  for( int i=0; i<3; i++ )
  {
    coords[i] = X[0][i]*( 1.0-alpha )*( 1.0-beta )*( 1.0-gamma )+
                X[1][i]*    alpha    *( 1.0-beta )*( 1.0-gamma )+
                X[2][i]*( 1.0-alpha )*    beta    *( 1.0-gamma )+
                X[3][i]*    alpha    *    beta    *( 1.0-gamma )+
                X[4][i]*( 1.0-alpha )*( 1.0-beta )*  gamma+
                X[5][i]*    alpha    *( 1.0-beta )*  gamma+
                X[6][i]*( 1.0-alpha )*    beta    *  gamma+
                X[7][i]*    alpha    *    beta    *  gamma;
  }
}

MeshLevel::MeshLevel( string const & name,
                      Group * const parent,
                      MeshLevel const & source,
                      CellBlockManagerABC & cellBlockManager,
                      int const order ):
  MeshLevel( name, parent, order )
{
  GEOSX_MARK_FUNCTION;

  // check that all elements are hexahedra
  source.m_elementManager->forElementRegions< CellElementRegion >( [&]( CellElementRegion const & sourceRegion )
  {
    sourceRegion.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion const & sourceSubRegion )
    {
      if( sourceSubRegion.getElementType() != ElementType::Hexahedron )
      {
        GEOSX_ERROR( "Orders higher than one are only available for hexahedral meshes" );
      }
    } );
  } );

  // constants for hex mesh
  localIndex const numVerticesPerCell = 8;
  localIndex const numNodesPerEdge = ( order+1 );
  localIndex const numNodesPerFace = ( order+1 )*( order+1 );
  localIndex const numNodesPerCell = ( order+1 )*( order+1 )*( order+1 );
  localIndex const numInternalNodesPerEdge = ( order-1 );
  localIndex const numInternalNodesPerFace = ( order-1 )*( order-1 );
  localIndex const numInternalNodesPerCell = ( order-1 )*( order-1 )*( order-1 );
  localIndex const numLocalVertices = source.m_nodeManager->size();
  localIndex const numLocalEdges = source.m_edgeManager->size();
  localIndex const numLocalFaces = source.m_faceManager->size();
  localIndex const maxVertexGlobalID = source.m_nodeManager->maxGlobalIndex() + 1;
  localIndex const maxEdgeGlobalID = source.m_edgeManager->maxGlobalIndex() + 1;
  localIndex const maxFaceGlobalID = source.m_faceManager->maxGlobalIndex() + 1;
  localIndex n1 = 0;
  source.m_elementManager->forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion const & sourceSubRegion )
  {
   n1 += sourceSubRegion.size();
   localIndex maxRegionIndex = sourceSubRegion.maxGlobalIndex() + 1;
  } );
  localIndex const numLocalCells = n1;
  ////////////////////////////////
  // Get the new number of nodes
  ////////////////////////////////
  localIndex numLocalNodes = numLocalVertices
                             + numLocalEdges * numInternalNodesPerEdge
                             + numLocalFaces * numInternalNodesPerFace
                             + numLocalCells * numInternalNodesPerCell;
  // initialize the node local to global map

  cellBlockManager.setNumNodes( numLocalNodes, order );

  arrayView1d< globalIndex const > const nodeLocalToGlobalSource = source.m_nodeManager->localToGlobalMap();
  arrayView1d< globalIndex > nodeLocalToGlobalNew = cellBlockManager.getNodeLocalToGlobal();

  // hash map that contains unique vertices and their indices for shared nodes.
  // a vertex is identified by 6 integers [i1 i2 i3 i4 a b], as follows
  // - nodes on a vertex v are identified by the vector [v -1 -1 -1 -1 -1]
  // - edge nodes are given by a linear interpolation between vertices v1 and v2, 'a' steps away from v1.
  //   We assume that v1 < v2 and identify these nodes with [v1 v2 -1 -1 a -1].
  // - face nodes are given by a bilinear interpolation between edges v1-v2 and v3-v4 (v1-v4 and v2-v3 are the diagonals), with interpolation parameters 'a' and 'b'.
  //   We assume that v1 is the smallest, and that v2 < v3. Then these nodes are identified with [v1 v2 v3 v4 a b]
  // - cell nodes are encountered only once, and thus do not need to be put in the hash map
  std::unordered_map< std::array< localIndex, 6 >, localIndex, NodeKeyHasher< localIndex >, NodeKeyEqual< localIndex > > nodeIDs;

  // Create new nodes, with local and global IDs
  localIndex localNodeID = 0;
  for( localIndex iter_vertex=0; iter_vertex < numLocalVertices; iter_vertex++)
  {
    nodeLocalToGlobalNew[ localNodeID ] = nodeLocalToGlobalSource[ iter_vertex ];
    nodeIDs[ createNodeKey( iter_vertex ) ] = localNodeID;
    localNodeID++;
  }

  //////////////////////////
  // Edges
  //////////////////////////

  // the total number of nodes: to add the number of non-vertex edge nodes
  m_edgeManager->resize( numLocalEdges );
  m_edgeManager->getDomainBoundaryIndicator() = source.m_edgeManager->getDomainBoundaryIndicator();

  arrayView1d< globalIndex const > const edgeLocalToGlobal = m_edgeManager->localToGlobalMap();
  for (localIndex iter_localToGlobalsize = 0; iter_localToGlobalsize < numLocalEdges; iter_localToGlobalsize++)
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
  m_edgeManager->nodeList().resize( numLocalEdges, numNodesPerEdge );
  // create / retrieve nodes on edges
  localIndex offset = maxVertexGlobalID;
  for( localIndex iter_edge = 0; iter_edge < numLocalEdges; iter_edge++ )
  {
    localIndex newEdgeNodes = 0;
    localIndex v1 = source.m_edgeManager->nodeList()[ iter_edge ][ 0 ];
    localIndex v2 = source.m_edgeManager->nodeList()[ iter_edge ][ 1 ];
    for( int q=0;q<numNodesPerEdge; q++ )
    {
      localIndex nodeID;
      std::array< localIndex, 6 > nodeKey = createNodeKey( v1, v2, q, order );
      if( nodeIDs.count( nodeKey ) == 0 )
      {
        // this is an internal edge node: create it
        nodeID = localNodeID;
        nodeIDs[ nodeKey ] = nodeID;
        nodeLocalToGlobalNew[ nodeID ] = offset + edgeLocalToGlobal[ iter_edge ] * numInternalNodesPerEdge + newEdgeNodes;
        localNodeID++;
        newEdgeNodes++;
      }
      else
      {
        nodeID = nodeIDs[ nodeKey ];
      }
      edgeToNodeMapNew[ iter_edge ][ q ] = nodeID;
    }
  }

  /////////////////////////
  // Faces
  //////////////////////////

  m_faceManager->resize( numLocalFaces );
  m_faceManager->faceCenter() = source.m_faceManager->faceCenter();
  m_faceManager->faceNormal() = source.m_faceManager->faceNormal();

  // copy the faces-to-edgs map from source
  m_faceManager->edgeList() = source.m_faceManager->edgeList();
  // copy the faces-to-elements map from source
  m_faceManager->elementList() = source.m_faceManager->elementList();
  m_faceManager->elementRegionList() = source.m_faceManager->elementRegionList();
  m_faceManager->elementSubRegionList() = source.m_faceManager->elementSubRegionList();
  // copy the faces-boundaryIndicator from source
  m_faceManager->getDomainBoundaryIndicator() = source.m_faceManager->getDomainBoundaryIndicator();

  arrayView1d< globalIndex const > const faceLocalToGlobal = m_faceManager->localToGlobalMap();
  for (localIndex iter_localToGlobalsize = 0; iter_localToGlobalsize < numLocalFaces; iter_localToGlobalsize++)
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


  // create / retrieve nodes on faces
  // by default, faces are oriented so that the normal has an outward orientation for the current rank.
  // For this reason :
  // - The 3rd and 4th node need to be swapped, to be consistent with the GL ordering
  // - The global IDs of the internal nodes must be referred to a global "reference" orientation (using the createNodeKey method)
  offset = maxVertexGlobalID + maxEdgeGlobalID * numInternalNodesPerEdge;
  for( localIndex iter_face = 0; iter_face < numLocalFaces; iter_face++ )
  {
    localIndex v1 = source.m_faceManager->nodeList()[ iter_face ][ 0 ];
    localIndex v2 = source.m_faceManager->nodeList()[ iter_face ][ 1 ];
    localIndex v3 = source.m_faceManager->nodeList()[ iter_face ][ 2 ];
    localIndex v4 = source.m_faceManager->nodeList()[ iter_face ][ 3 ];
    std::swap(v3, v4);
    globalIndex gv1 = nodeLocalToGlobalSource[v1];
    globalIndex gv2 = nodeLocalToGlobalSource[v2];
    globalIndex gv3 = nodeLocalToGlobalSource[v3];
    globalIndex gv4 = nodeLocalToGlobalSource[v4];
    for( int q1=0;q1<numNodesPerEdge; q1++ )
    {
      for( int q2=0;q2<numNodesPerEdge; q2++ )
      {
        localIndex nodeID;
        std::array< localIndex, 6 > nodeKey = createNodeKey( v1, v2, v3, v4, q1, q2, order );
        if( nodeIDs.count( nodeKey ) == 0 )
        {
          // this is an internal face node: create it
          nodeID = localNodeID;
          nodeIDs[ nodeKey ] = nodeID;
          std::array< globalIndex, 6 > referenceOrientation = createNodeKey( gv1, gv2, gv3, gv4, q1, q2, order );
          int gq1 = referenceOrientation[4];
          int gq2 = referenceOrientation[5];
          nodeLocalToGlobalNew[ nodeID ] = offset + faceLocalToGlobal[ iter_face ] * numInternalNodesPerFace + gq1* numInternalNodesPerEdge + gq2; 
          localNodeID++;
        }
        else
        {
          nodeID = nodeIDs[ nodeKey ];
        }
        faceToNodeMapNew[ iter_face ][ q2 + q1*numNodesPerEdge ] = nodeID;
      }
    }
  }

  // add all nodes to the target set "all"
  m_nodeManager->resize( numLocalNodes );
  
  SortedArray< localIndex > & allNodesSet = cellBlockManager.getNodeSets()[ "all" ];
  allNodesSet.reserve( numLocalNodes );

  for( localIndex iter_nodes=0; iter_nodes< numLocalNodes; ++iter_nodes )
  {
    allNodesSet.insert( iter_nodes );
  }

  /////////////////////////
  // Elements
  //////////////////////////
  // also assign node coordinates using trilinear interpolation in th elements
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const refPosSource = source.m_nodeManager->referencePosition();
  arrayView2d< real64, nodes::REFERENCE_POSITION_USD > refPosNew = cellBlockManager.getNodePositions();
  refPosNew.setValues< parallelHostPolicy >( -1.0 );

  real64 Xmesh[ numVerticesPerCell ][ 3 ] = { { } };
  real64 X[ 3 ] = { { } };
  array1d< real64 > glCoords = gaussLobattoPoints( order );
  localIndex elemMeshVertices[ numVerticesPerCell ] = { };
  offset = maxVertexGlobalID + maxEdgeGlobalID * numInternalNodesPerEdge + maxFaceGlobalID * numInternalNodesPerFace;
  std::array< localIndex, 6 > const nullKey = std::array< localIndex, 6 >{ -1, -1 ,-1, -1, -1, -1 };
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
      newSubRegion.resizePerElementValues( numNodesPerCell,
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

      // copy the elements-to-edges map from source
      arrayView2d< localIndex const > const & elemsToEdgesSource = sourceSubRegion.edgeList();
      array2d< localIndex > & elemsToEdgesNew = newSubRegion.edgeList();
      elemsToEdgesNew = elemsToEdgesSource;


      // initialize the elements-to-nodes map
      arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodesSource = sourceSubRegion.nodeList().toViewConst();
      arrayView2d< localIndex, cells::NODE_MAP_USD > elemsToNodesNew = cellBlockManager.getElemToNodes(newSubRegion.getName(), numNodesPerCell);

      arrayView1d< globalIndex const > const elementLocalToGlobal = sourceSubRegion.localToGlobalMap();
      for (localIndex iter_localToGlobalsize = 0; iter_localToGlobalsize < sourceSubRegion.localToGlobalMap().size(); iter_localToGlobalsize++)
      {
        newSubRegion.localToGlobalMap()[iter_localToGlobalsize] = sourceSubRegion.localToGlobalMap()[iter_localToGlobalsize];
      }
      newSubRegion.constructGlobalToLocalMap();

      // then loop through all the elements and assign the globalID according to the globalID of the Element
      // and insert the new local to global ID ( for the internal nodes of elements ) into the nodeLocalToGlobal
      // retrieve finite element type
      for( localIndex iter_elem = 0; iter_elem < numLocalCells; ++iter_elem )
      {
        localIndex newCellNodes = 0;
        for( localIndex iter_vertex = 0; iter_vertex < numVerticesPerCell; iter_vertex++ )
        {
          elemMeshVertices[ iter_vertex ] = elemsToNodesSource[ iter_elem ][ iter_vertex ];
          for( int i =0; i < 3; i++ )
          {
            Xmesh[ iter_vertex ][ i ] = refPosSource[ elemMeshVertices[ iter_vertex ] ][ i ];
          }
        }

        for( int q = 0; q < numNodesPerCell; q++ )
        {
          localIndex nodeID;
          int dof = q;
          int q1 = dof % numNodesPerEdge;
          dof /= ( numNodesPerEdge );
          int q2 = dof % ( numNodesPerEdge );
          dof /= ( numNodesPerEdge );
          int q3 = dof % ( numNodesPerEdge );
          // compute node coords
          real64 alpha = ( glCoords[ q1 ] + 1.0 ) / 2.0;
          real64 beta = ( glCoords[ q2 ] + 1.0 ) / 2.0;
          real64 gamma = ( glCoords[ q3 ] + 1.0 ) / 2.0;
          trilinearInterp( alpha, beta, gamma, Xmesh, X );
          // find node ID
          std::array< localIndex, 6 > nodeKey = createNodeKey( elemMeshVertices, q1, q2, q3, order );
          if( nodeKey == nullKey )
          {
            // the node is internal to a cell -- create it
            nodeID = localNodeID;
            nodeLocalToGlobalNew[ nodeID ] = offset + elementLocalToGlobal[ iter_elem ] * numInternalNodesPerCell + newCellNodes;
            localNodeID++;
            newCellNodes++;
          }
          else
          {
            nodeID = nodeIDs[ nodeKey ];
          }
          for( int i=0; i<3; i++ )
          {
            refPosNew( nodeID, i ) = X[ i ];
          }
          elemsToNodesNew[ iter_elem ][ q ] = nodeID;
        }
      }
    } );
  } );
  this->generateSets();
}

MeshLevel::MeshLevel( string const & name,
                      Group * const parent,
                      MeshLevel const & source,
                      int const order ):
  MeshLevel( name, parent, order )
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
  m_faceManager->elementRegionList() = source.m_faceManager->elementRegionList();
  m_faceManager->elementSubRegionList() = source.m_faceManager->elementSubRegionList();

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
      newSubRegion.setElementType( sourceSubRegion.getElementType() );

      newSubRegion.resizePerElementValues( numNodesPerElem,
                                           sourceSubRegion.numEdgesPerElement(),
                                           sourceSubRegion.numFacesPerElement() );

      newSubRegion.resize( sourceSubRegion.size() );

      arrayView2d< localIndex const > const & elemToFaces = sourceSubRegion.faceList();

      arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodesSource = sourceSubRegion.nodeList().toViewConst();
      array2d< localIndex, cells::NODE_MAP_PERMUTATION > & elemsToNodesNew = newSubRegion.nodeList();

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
      for( localIndex e = 0; e < elemsToNodesNew.size( 0 ); ++e )
      {
        for( localIndex i = 0; i < numNodesPerElem; i++ )
        {
          elemsToNodesNew[e][i] =-1;
        }

      }

      localIndex count=0;

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
                    faceToNodeMapNew[elemToFaces[elem][face]][m] = count;
                  }
                }
                face++;
              }
              if( face > 5 && foundFace < 1 )
              {
                elemsToNodesNew[elem][nodeindex] = count;
                count++;
              }
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
  if( getOrder() > 1 )
  {
    m_elementManager->forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion & subRegion )
    {
      subRegion.calculateElementGeometricQuantities( *m_nodeManager, *m_faceManager );
    } );
  }
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

int MeshLevel::getOrder() const
{
  return m_order;
}

} /* namespace geosx */
