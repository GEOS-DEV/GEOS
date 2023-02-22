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

struct NodeKeyHasher {
  std::size_t operator()(const std::array< localIndex, 6 >& arr) const {
    std::size_t hash = 0;
    // use a boost-style hash function
    for (auto v : arr) {
        hash ^= std::hash<localIndex>{}( v )  + 0x9e3779b9 + ( hash << 6 ) + ( hash >> 2 ); 
    }
    GEOSX_LOG_RANK ("!!!! INFO !!!! creating hash for array "<< arr[0] << " " << arr[1] << " " << arr[2] << " " << arr[3] << " " << arr[4] << " " << arr[5] <<" : " << hash );
    return hash;
  }   
};

struct NodeKeyEqual {
  bool operator()(const std::array< localIndex, 6 > & lhs, const std::array< localIndex, 6 > & rhs) const
  {
    GEOSX_LOG_RANK ("!!!! INFO !!!! comparing array "<< lhs[0] << " " << lhs[1] << " " << lhs[2] << " " << lhs[3] << " " << lhs[4] << " " << lhs[5] << " with " << rhs[0] << " " << rhs[1] << " " << rhs[2] << " " << rhs[3] << " " << rhs[4] << " " << rhs[5] );
     for( int i=0; i< 6; i++ )
     { 
       if( lhs[ i ] != rhs[ i ] )
       {
    GEOSX_LOG_RANK ("!!!! INFO !!!! they are different." );
       
         return false;
       }
     }
    GEOSX_LOG_RANK ("!!!! INFO !!!! they are equal." );
     return true;
  }
};

static std::array< localIndex, 6 > createNodeKey( localIndex v )
{
  GEOSX_LOG_RANK ("!!!! INFO !!!! creating key for vertex "<< v );
  return std::array< localIndex, 6 > { v, -1, -1 ,-1 ,-1, -1 };
}

static std::array< localIndex, 6 > createNodeKey( localIndex v1, localIndex v2, int a, int order )
{
  if( a == 0) return createNodeKey( v1 );
  if( a == order ) return createNodeKey( v2 );
  if( v1 < v2 )
  {
  GEOSX_LOG_RANK ("!!!! INFO !!!! creating key for edge "<< v1 << " " << v2 << ", location " << a );
    return std::array< localIndex, 6 > { v1, v2, -1 ,-1 ,a , -1 };
  }
  else
  {
  GEOSX_LOG_RANK ("!!!! INFO !!!! creatingkey for edge "<< v2 << " " << v1 << ", location "<< (order - a) );
    return std::array< localIndex, 6 > { v2, v1, -1 ,-1 ,order - a , -1 };
  }
}


static std::array< localIndex, 6 > createNodeKey( localIndex v1, localIndex v2, localIndex v3, localIndex v4, int a, int b, int order )
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
  GEOSX_LOG_RANK ("!!!! INFO !!!! creating key for face "<< v1 << " " << v2 << " "<< v3 << " "<< v4 << ", location " << a << " " << b );
  return std::array< localIndex, 6 > { v1, v2, v3, v4, a, b };
}

static std::array< localIndex, 6 > createNodeKey( localIndex const (&elemNodes)[ 8 ], int q1, int q2, int q3, int order )
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
                      int const order,
                      int const strategy ):
  MeshLevel( name, parent )
{

if ( strategy == 1 )
{

  GEOSX_MARK_FUNCTION;

  ////////////////////////////////
  // Get the new number of nodes 
  ////////////////////////////////

  // the total number of nodes: to add the number of non-vertex edge nodes
  localIndex const numBasisSupportPoints = order+1;
  localIndex const numNodesPerEdge = numBasisSupportPoints;
  localIndex const numNonVertexNodesPerEdge = numNodesPerEdge - 2; 
  localIndex const numInternalEdgeNodes = source.m_edgeManager->size() * numNonVertexNodesPerEdge;

  // the total number of nodes: to add the number of non-edge face nodes
  localIndex const numNodesPerFace = pow( order+1, 2 ); 
  localIndex const numNonEdgeNodesPerFace = pow( order-1, 2 ); 
  localIndex numInternalFaceNodes = 0; 

  ArrayOfArraysView< localIndex const > const & facesToEdges = source.m_faceManager->edgeList().toViewConst();
  for( localIndex kf=0; kf<source.m_faceManager->size(); ++kf )
  {
    localIndex const numEdgesPerFace = facesToEdges.sizeOfArray( kf );
    if( numEdgesPerFace==4 )
      numInternalFaceNodes += numNonEdgeNodesPerFace;
    else 
    {    
      GEOSX_ERROR( "need more support for face geometry" );
    }    
  }

  // the total number of nodes: to add the number of non-face element nodes
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

  // calculate the total number of nodes
  localIndex const numNodes = source.m_nodeManager->size()
                              + numInternalEdgeNodes
                              + numInternalFaceNodes
                              + numInternalElementNodes;

  //
  // ---- to initialize the node local to global map ----
  //   

  cellBlockManager.setNumNodes( numNodes, order );

  arrayView1d< globalIndex const > const sourceNodeLocalToGlobal = source.m_nodeManager->localToGlobalMap();
  arrayView1d< globalIndex > nodeLocalToGlobal = cellBlockManager.getNodeLocalToGlobal();

  //GEOSX_LOG_RANK ("!!!! INFO !!!! nodeLocalToGlobal = "<< nodeLocalToGlobal );
  //GEOSX_LOG_RANK ("!!!! INFO !!!! sourceNodeLocalToGlobal = "<< sourceNodeLocalToGlobal);



  //////////////////////////
  // Edges
  //////////////////////////

  // the total number of nodes: to add the number of non-vertex edge nodes
  m_edgeManager->resize( source.m_edgeManager->size() );
  m_edgeManager->getDomainBoundaryIndicator() = source.m_edgeManager->getDomainBoundaryIndicator();

  arrayView1d< globalIndex const > const edgeLocalToGlobal = m_edgeManager->localToGlobalMap();
  for (localIndex iter_localToGlobalsize = 0; iter_localToGlobalsize < source.m_edgeManager->localToGlobalMap().size(); iter_localToGlobalsize++)
  {
    m_edgeManager->localToGlobalMap()[iter_localToGlobalsize] = source.m_edgeManager->localToGlobalMap()[iter_localToGlobalsize];
  }
  m_edgeManager->constructGlobalToLocalMap();
  //GEOSX_LOG_RANK ("!!!! INFO !!!! edgeLocalToGlobal = "<< edgeLocalToGlobal <<"; numInternalEdgeNodes="<<numInternalEdgeNodes);
  //GEOSX_LOG_RANK_0 ("!!!! INFO !!!! test="<<test);

  //
  // ---- initialize edge-to-node map ----
  //

  // get information from the source (base mesh-level) edge-to-node map
  arrayView2d< localIndex const > const & edgeToNodeMapSource = source.m_edgeManager->nodeList();
  //array2d< localIndex > & edgeToNodeMapNew = m_edgeManager->nodeList();
  arrayView2d< localIndex > edgeToNodeMapNew = cellBlockManager.getEdgeToNodes();
  m_edgeManager->nodeList().resize( source.m_edgeManager->size(), numNodesPerEdge );


  /////////////////////////
  // Faces
  //////////////////////////

  m_faceManager->resize( source.m_faceManager->size() );
  m_faceManager->faceCenter() = source.m_faceManager->faceCenter();
  m_faceManager->faceNormal() = source.m_faceManager->faceNormal();

  // copy the faces-to-edgs map from source
  m_faceManager->edgeList() = source.m_faceManager->edgeList();
  // copy the faces-to-elements map from source
  m_faceManager->elementList() = source.m_faceManager->elementList();
  // copy the faces-boundaryIndicator from source
  m_faceManager->getDomainBoundaryIndicator() = source.m_faceManager->getDomainBoundaryIndicator();

  arrayView1d< globalIndex const > const faceLocalToGlobal = m_faceManager->localToGlobalMap();
  for (localIndex iter_localToGlobalsize = 0; iter_localToGlobalsize < source.m_faceManager->localToGlobalMap().size(); iter_localToGlobalsize++)
  {
    m_faceManager->localToGlobalMap()[iter_localToGlobalsize] = source.m_faceManager->localToGlobalMap()[iter_localToGlobalsize];
  }
  m_faceManager->constructGlobalToLocalMap();
  //GEOSX_LOG_RANK ("!!!! INFO !!!! faceLocalToGlobal = "<< faceLocalToGlobal <<"; numInternalFaceNodes="<<numInternalFaceNodes);

  ArrayOfArraysView< localIndex const > const & faceToNodeMapSource = source.m_faceManager->nodeList().toViewConst();
  ArrayOfArrays< localIndex > & faceToNodesRelation = m_faceManager->nodeList();
  ArrayOfArrays< localIndex > & faceToNodeMapNew = cellBlockManager.getFaceToNodes();

  // number of elements in each row of the map as capacity
  array1d< localIndex > counts( faceToNodeMapNew.size());
  counts.setValues< parallelHostPolicy >( numNodesPerFace );

  //  reconstructs the faceToNodeMap with the provided capacity in counts
  faceToNodeMapNew.resizeFromCapacities< parallelHostPolicy >( faceToNodeMapNew.size(), counts.data() );
  faceToNodesRelation.resizeFromCapacities< parallelHostPolicy >( faceToNodesRelation.size(), counts.data() );

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
  
  //GEOSX_LOG_RANK ("!!!! INFO !!!! faceToNodeMapSource = "<< faceToNodeMapSource ); 
  //GEOSX_LOG_RANK ("!!!! INFO !!!! faceToNodeMapNew = "<< faceToNodeMapNew ); 
  //GEOSX_LOG_RANK ("!!!! INFO !!!! faceToNodesRelation = "<< faceToNodesRelation ); 


  //////////////////////////
  // Nodes
  //////////////////////////

  m_nodeManager->resize( numNodes );

  Group & nodeSets = m_nodeManager->sets();
  SortedArray< localIndex > & allNodes  = nodeSets.registerWrapper< SortedArray< localIndex > >( string( "all" ) ).reference();
  allNodes.reserve( m_nodeManager->size() );

  for( localIndex a=0; a<m_nodeManager->size(); ++a )
  {
    allNodes.insert( a );
  }

  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const refPosSource = source.m_nodeManager->referencePosition();
  //arrayView2d< real64, nodes::REFERENCE_POSITION_USD > const refPosNew = m_nodeManager->referencePosition().toView();
  arrayView2d< real64, nodes::REFERENCE_POSITION_USD > refPosNew = cellBlockManager.getNodePositions();


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

      // the number of nodes per elements depends on the order number
      localIndex const numNodesPerElem = pow( order+1, 3 );

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
      arrayView2d< localIndex const > const & elemsToFaces = sourceSubRegion.faceList();
      array2d< localIndex > & elemsToFacesNew = newSubRegion.faceList();
      elemsToFacesNew = elemsToFaces;

      // copy the elements-to-edgs map from source
      arrayView2d< localIndex const > const & elemsToEdges = sourceSubRegion.edgeList();
      array2d< localIndex > & elemsToEdgesNew = newSubRegion.edgeList();
      elemsToEdgesNew = elemsToEdges;


      // initialize the elements-to-nodes map 
      arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodesSource = sourceSubRegion.nodeList().toViewConst();
      //array2d< localIndex, cells::NODE_MAP_PERMUTATION > & elemsToNodesNew = newSubRegion.nodeList();
      arrayView2d< localIndex, cells::NODE_MAP_USD > elemsToNodesNew = cellBlockManager.getElemToNodes(newSubRegion.getName(), numNodesPerElem);
      //GEOSX_LOG_RANK ("!!!! INFO !!!! elemsToNodesNew = "<< elemsToNodesNew.size(1));

      arrayView1d< globalIndex const > const elementLocalToGlobal = sourceSubRegion.localToGlobalMap();

      for (localIndex iter_localToGlobalsize = 0; iter_localToGlobalsize < sourceSubRegion.localToGlobalMap().size(); iter_localToGlobalsize++)
      {
        newSubRegion.localToGlobalMap()[iter_localToGlobalsize] = sourceSubRegion.localToGlobalMap()[iter_localToGlobalsize];
      }
      newSubRegion.constructGlobalToLocalMap();

      //GEOSX_LOG_RANK ("!!!! INFO !!!! elementLocalToGlobal="<<elementLocalToGlobal); 

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

      // initiate the localNodeID to be the totalnumber of local nodes from the current rank
      localIndex localNodeID = source.m_nodeManager->size();

      // firstly loop through all the nodes from source ( base-mesh level )
      // and copy the nodeLocalToGlobal map for the existing local nodes
      for( localIndex iter_node = 0; iter_node < localNodeID; ++iter_node )
      {
        nodeLocalToGlobal[iter_node] =  sourceNodeLocalToGlobal[iter_node];
      }

      // secondly loop through all the edges and assign the globalID according to the globalID of the Edge
      // and insert the new local to global ID ( for the internal edge nodes ) into the nodeLocalToGlobal
      globalIndex const firstGlobalNodeIDforEdges = source.m_nodeManager->maxGlobalIndex() + 1;
      localIndex const numOfNodesPerEdgeForBase = 2;

      for( localIndex iter_edge = 0; iter_edge < m_edgeManager->size(); ++iter_edge )
      {
        for ( localIndex iter_nodesPerEdgeBase  = 0; iter_nodesPerEdgeBase < numOfNodesPerEdgeForBase; iter_nodesPerEdgeBase++)
        {
          edgeToNodeMapNew[iter_edge][ iter_nodesPerEdgeBase ] = edgeToNodeMapSource[iter_edge][ iter_nodesPerEdgeBase ];
        }

        globalIndex edgeGlobalID = edgeLocalToGlobal[iter_edge];
        for (localIndex iter_internalnode = 0; iter_internalnode < numNonVertexNodesPerEdge;  ++iter_internalnode)
        {
          nodeLocalToGlobal[localNodeID] = firstGlobalNodeIDforEdges + edgeGlobalID * numNonVertexNodesPerEdge +  iter_internalnode;
          edgeToNodeMapNew[iter_edge][ numOfNodesPerEdgeForBase + iter_internalnode ] = localNodeID;

          refPosNew[localNodeID][0] = refPosSource[ edgeToNodeMapSource[iter_edge][0]][0]
                                    + ((refPosSource[ edgeToNodeMapSource[iter_edge][1]][0]-refPosSource[ edgeToNodeMapSource[iter_edge][0]][0])/2.)
                                      * GaussLobattoPts[iter_internalnode + 1]
                                    + (refPosSource[ edgeToNodeMapSource[iter_edge][1]][0]-refPosSource[ edgeToNodeMapSource[iter_edge][0]][0])/2;

          refPosNew[localNodeID][1] = refPosSource[ edgeToNodeMapSource[iter_edge][0]][1]
                                    + ((refPosSource[ edgeToNodeMapSource[iter_edge][1]][1]-refPosSource[ edgeToNodeMapSource[iter_edge][0]][1])/2.)
                                      * GaussLobattoPts[iter_internalnode + 1]
                                    + (refPosSource[ edgeToNodeMapSource[iter_edge][1]][1]-refPosSource[ edgeToNodeMapSource[iter_edge][0]][1])/2;

          refPosNew[localNodeID][2] = refPosSource[ edgeToNodeMapSource[iter_edge][0]][2]
                                    + ((refPosSource[ edgeToNodeMapSource[iter_edge][1]][2]-refPosSource[ edgeToNodeMapSource[iter_edge][0]][2])/2.)
                                      * GaussLobattoPts[iter_internalnode + 1]
                                    + (refPosSource[ edgeToNodeMapSource[iter_edge][1]][2]-refPosSource[ edgeToNodeMapSource[iter_edge][0]][2])/2;
          localNodeID++;
        }
      }

      //GEOSX_LOG_RANK ("!!!! INFO !!!! refPosNew = "<< refPosNew );
      //GEOSX_LOG_RANK ("!!!! INFO !!!! edgeToNodeMapSource = "<< edgeToNodeMapSource);
      //GEOSX_LOG_RANK ("!!!! INFO !!!! edgeToNodeMapNew = "<< edgeToNodeMapNew);

      // then loop through all the faces and assign the globalID according to the globalID of the Face 
      // and insert the new local to global ID ( for the internal nodes of faces ) into the nodeLocalToGlobal
      globalIndex const firstGlobalNodeIDforFaces = firstGlobalNodeIDforEdges + numInternalEdgeNodes * ( source.m_edgeManager->maxGlobalIndex() + 1 );
      localIndex const numOfNodesPerFaceForBase = 4;
      localIndex const numOfEdgesPerFace = 4;


      for( localIndex iter_face = 0; iter_face < m_faceManager->size(); ++iter_face )
      {
        for ( localIndex iter_nodesPerFaceBase  = 0; iter_nodesPerFaceBase < numOfNodesPerFaceForBase; iter_nodesPerFaceBase++)
        {
          faceToNodeMapNew[iter_face][ iter_nodesPerFaceBase ] = faceToNodeMapSource[iter_face][ iter_nodesPerFaceBase ];
        }

        globalIndex faceGlobalID = faceLocalToGlobal[iter_face];

        for( localIndex iter_edge = 0; iter_edge < numOfEdgesPerFace; ++iter_edge )
        {
          for (localIndex iter_internalnode = 0; iter_internalnode < numNonVertexNodesPerEdge;  ++iter_internalnode)
            faceToNodeMapNew[iter_face][ numOfNodesPerFaceForBase + iter_edge * numNonVertexNodesPerEdge + iter_internalnode ] =
              edgeToNodeMapNew[facesToEdges[iter_face][iter_edge]][numOfNodesPerEdgeForBase + iter_internalnode] ;
        }

        real64 refPosXmax = std::max( std::max( refPosSource[faceToNodeMapSource[iter_face][0]][0], refPosSource[faceToNodeMapSource[iter_face][1]][0] ),
                                refPosSource[faceToNodeMapSource[iter_face][2]][0] );
        real64 refPosYmax = std::max( std::max( refPosSource[faceToNodeMapSource[iter_face][0]][1], refPosSource[faceToNodeMapSource[iter_face][1]][1] ),
                                refPosSource[faceToNodeMapSource[iter_face][2]][1] );
        real64 refPosZmax = std::max( std::max( refPosSource[faceToNodeMapSource[iter_face][0]][2], refPosSource[faceToNodeMapSource[iter_face][1]][2] ),
                                refPosSource[faceToNodeMapSource[iter_face][2]][2] );

        for (localIndex iter_internalnodeA = 0; iter_internalnodeA < numNonVertexNodesPerEdge;  ++iter_internalnodeA)
        {
          for (localIndex iter_internalnodeB = 0; iter_internalnodeB < numNonVertexNodesPerEdge;  ++iter_internalnodeB)
          {
            localIndex iter_internalnode = iter_internalnodeA + iter_internalnodeB * numNonVertexNodesPerEdge;
            nodeLocalToGlobal[localNodeID]= firstGlobalNodeIDforFaces + faceGlobalID * numNonEdgeNodesPerFace +  iter_internalnode;
            faceToNodeMapNew[iter_face][ numOfNodesPerFaceForBase + numOfEdgesPerFace * numNonVertexNodesPerEdge + iter_internalnode] = localNodeID ;

            if ( std::fabs ( refPosXmax - refPosSource[faceToNodeMapSource[iter_face][0]][0] ) < 1e-5 )
            {
              refPosNew[localNodeID][0] = refPosSource[faceToNodeMapSource[iter_face][0]][0];
              refPosNew[localNodeID][1] = refPosSource[faceToNodeMapSource[iter_face][0]][1]
                                        + ((refPosYmax-refPosSource[faceToNodeMapSource[iter_face][0]][1])/2.)*GaussLobattoPts[iter_internalnodeA + 1]
                                        + (refPosYmax-refPosSource[faceToNodeMapSource[iter_face][0]][1])/2;
              refPosNew[localNodeID][2] = refPosSource[faceToNodeMapSource[iter_face][0]][2]
                                        + ((refPosZmax-refPosSource[faceToNodeMapSource[iter_face][0]][2])/2.)*GaussLobattoPts[iter_internalnodeB + 1]
                                        + (refPosZmax-refPosSource[faceToNodeMapSource[iter_face][0]][2])/2;
            }

            else if ( std::fabs ( refPosYmax - refPosSource[faceToNodeMapSource[iter_face][0]][1] ) < 1e-5 )
            {
              refPosNew[localNodeID][0] = refPosSource[faceToNodeMapSource[iter_face][0]][0]
                                        + ((refPosXmax-refPosSource[faceToNodeMapSource[iter_face][0]][0])/2.)*GaussLobattoPts[iter_internalnodeA + 1]
                                        + (refPosXmax-refPosSource[faceToNodeMapSource[iter_face][0]][0])/2;
              refPosNew[localNodeID][1] = refPosSource[faceToNodeMapSource[iter_face][0]][1];
              refPosNew[localNodeID][2] = refPosSource[faceToNodeMapSource[iter_face][0]][2]
                                        + ((refPosZmax-refPosSource[faceToNodeMapSource[iter_face][0]][2])/2.)*GaussLobattoPts[iter_internalnodeB + 1]
                                        + (refPosZmax-refPosSource[faceToNodeMapSource[iter_face][0]][2])/2;
            }
            else
            {
              refPosNew[localNodeID][0] = refPosSource[faceToNodeMapSource[iter_face][0]][0]
                                        + ((refPosXmax-refPosSource[faceToNodeMapSource[iter_face][0]][0])/2.)*GaussLobattoPts[iter_internalnodeA + 1]
                                        + (refPosXmax-refPosSource[faceToNodeMapSource[iter_face][0]][0])/2;
              refPosNew[localNodeID][1] = refPosSource[faceToNodeMapSource[iter_face][0]][1]
                                        + ((refPosYmax-refPosSource[faceToNodeMapSource[iter_face][0]][1])/2.)*GaussLobattoPts[iter_internalnodeB + 1]
                                        + (refPosYmax-refPosSource[faceToNodeMapSource[iter_face][0]][1])/2;
              refPosNew[localNodeID][2] = refPosSource[faceToNodeMapSource[iter_face][0]][2];
            }

            localNodeID++;
          }
        }
      }

      //GEOSX_LOG_RANK ("!!!! INFO !!!! refPosNew = "<< refPosNew );
      //GEOSX_LOG_RANK ("!!!! INFO !!!! faceToNodeMapNew = "<< faceToNodeMapNew);
      //GEOSX_LOG_RANK ("!!!! INFO !!!! faceToNodeMapSource = "<< faceToNodeMapSource);


      // then loop through all the elements and assign the globalID according to the globalID of the Element 
      // and insert the new local to global ID ( for the internal nodes of elements ) into the nodeLocalToGlobal
      globalIndex const firstGlobalNodeIDforElems = firstGlobalNodeIDforFaces + numInternalFaceNodes * ( source.m_faceManager->maxGlobalIndex() + 1 );
      localIndex const numOfNodesPerElementForBase = 8;
      localIndex const numOfEdgesPerElement = 12;
      localIndex const numOfFacesPerElement = 6;
      localIndex const numInternalNodesPerElement = ( order - 1 ) * ( order - 1 ) * ( order - 1 );


      for( localIndex iter_elem = 0; iter_elem < elemsToNodesSource.size( 0 ); ++iter_elem )
      {
        for ( localIndex iter_nodesPerElemBase  = 0; iter_nodesPerElemBase < numOfNodesPerElementForBase; iter_nodesPerElemBase++)
        {
          elemsToNodesNew[iter_elem][ iter_nodesPerElemBase ] = elemsToNodesSource[iter_elem][ iter_nodesPerElemBase ];
          refPosNew[elemsToNodesSource[iter_elem][iter_nodesPerElemBase]][0] = refPosSource[elemsToNodesSource[iter_elem][iter_nodesPerElemBase]][0];
          refPosNew[elemsToNodesSource[iter_elem][iter_nodesPerElemBase]][1] = refPosSource[elemsToNodesSource[iter_elem][iter_nodesPerElemBase]][1];
          refPosNew[elemsToNodesSource[iter_elem][iter_nodesPerElemBase]][2] = refPosSource[elemsToNodesSource[iter_elem][iter_nodesPerElemBase]][2];
        }

        globalIndex elemGlobalID = elementLocalToGlobal[iter_elem];

        //Three 1D arrays to contains the GL points in the new coordinates knowing the mesh nodes
        array1d< real64 > x( order+1 );
        array1d< real64 > y( order+1 );
        array1d< real64 > z( order+1 );

        //Fill the three 1D array
        for( localIndex i = 0; i < order+1; i++ )
        {
          x[i] = refPosSource[elemsToNodesSource[iter_elem][0]][0] +
                 ((refPosSource[elemsToNodesSource[iter_elem][1]][0]-refPosSource[elemsToNodesSource[iter_elem][0]][0])/2.)*GaussLobattoPts[i]
                 + (refPosSource[elemsToNodesSource[iter_elem][1]][0]-refPosSource[elemsToNodesSource[iter_elem][0]][0])/2;
          y[i] = refPosSource[elemsToNodesSource[iter_elem][0]][1] +
                 ((refPosSource[elemsToNodesSource[iter_elem][2]][1]-refPosSource[elemsToNodesSource[iter_elem][0]][1])/2.)*GaussLobattoPts[i]
                 + (refPosSource[elemsToNodesSource[iter_elem][2]][1]-refPosSource[elemsToNodesSource[iter_elem][0]][1])/2;
          z[i] = refPosSource[elemsToNodesSource[iter_elem][0]][2] +
                 ((refPosSource[elemsToNodesSource[iter_elem][4]][2]-refPosSource[elemsToNodesSource[iter_elem][0]][2])/2.)*GaussLobattoPts[i]
                 + (refPosSource[elemsToNodesSource[iter_elem][4]][2]-refPosSource[elemsToNodesSource[iter_elem][0]][2])/2;
        }

        for( localIndex iter_edge = 0; iter_edge < numOfEdgesPerElement; ++iter_edge )
        {
          for (localIndex iter_internalnode = 0; iter_internalnode < numNonVertexNodesPerEdge;  ++iter_internalnode)
            elemsToNodesNew[iter_elem][ numOfNodesPerElementForBase + iter_edge * numNonVertexNodesPerEdge + iter_internalnode ]
              = edgeToNodeMapNew[elemsToEdges[iter_elem][iter_edge]][numOfNodesPerEdgeForBase + iter_internalnode];
        }

        localIndex elemToNodeCounter = numOfNodesPerElementForBase + numOfEdgesPerElement * numNonVertexNodesPerEdge;
        for( localIndex iter_face = 0; iter_face < numOfFacesPerElement; ++iter_face )
        {
          for (localIndex iter_internalnode = 0; iter_internalnode < numNonEdgeNodesPerFace;  ++iter_internalnode)
            elemsToNodesNew[iter_elem][ elemToNodeCounter + iter_face * numNonEdgeNodesPerFace + iter_internalnode ]
              = faceToNodeMapNew[elemsToFaces[iter_elem][iter_face]][numOfNodesPerFaceForBase + numOfEdgesPerFace * numNonVertexNodesPerEdge + iter_internalnode];
        }

        elemToNodeCounter += numOfFacesPerElement * numNonEdgeNodesPerFace;
        for (localIndex iter_internalnodeZ = 0; iter_internalnodeZ < numNonVertexNodesPerEdge;  ++iter_internalnodeZ)
        {
          for (localIndex iter_internalnodeY = 0; iter_internalnodeY < numNonVertexNodesPerEdge;  ++iter_internalnodeY)
          {
            for (localIndex iter_internalnodeX = 0; iter_internalnodeX < numNonVertexNodesPerEdge;  ++iter_internalnodeX)
            {
              localIndex iter_internalnode = iter_internalnodeX + iter_internalnodeY * numNonVertexNodesPerEdge
                                           + iter_internalnodeZ * numNonVertexNodesPerEdge * numNonVertexNodesPerEdge;
              nodeLocalToGlobal[localNodeID] = firstGlobalNodeIDforElems + elemGlobalID * numInternalNodesPerElement +  iter_internalnode;
              elemsToNodesNew[ iter_elem ][ elemToNodeCounter + iter_internalnode ] = localNodeID;

              refPosNew[localNodeID][0] = x[iter_internalnodeX+1];
              refPosNew[localNodeID][1] = y[iter_internalnodeY+1];
              refPosNew[localNodeID][2] = z[iter_internalnodeZ+1];

              localNodeID++;
            }
          }
        }
      }

      //GEOSX_LOG_RANK ("!!!! INFO !!!! refPosNew = "<< refPosNew );
      //GEOSX_LOG_RANK ("!!!! INFO !!!! nodeLocalToGlobal = "<< nodeLocalToGlobal );
      //GEOSX_LOG_RANK ("!!!! INFO !!!! sourceNodeLocalToGlobal = "<< sourceNodeLocalToGlobal);

      //GEOSX_LOG_RANK ("!!!! INFO !!!! faceToNodeMapNew = "<< faceToNodeMapNew);
      //GEOSX_LOG_RANK ("!!!! INFO !!!! edgeToNodeMapNew = "<< edgeToNodeMapNew);
      //GEOSX_LOG_RANK ("!!!! INFO !!!! elemsToNodesNew = "<< elemsToNodesNew);
      //GEOSX_LOG_RANK ("!!!! INFO !!!! elemsToNodesSource = "<< elemsToNodesSource);
      //GEOSX_LOG_RANK ("!!!! INFO !!!! refPosNew = "<< refPosNew );


    } );
  } );
}
if ( strategy == 2 )
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
  // localIndex const numVerticesPerEdge = 2;
  // localIndex const numVerticesPerFace = 4;
  localIndex const numVerticesPerCell = 8;
  // localIndex const numEdgesPerFace = 4;
  // localIndex const numEdgesPerCell = 12;
  // localIndex const numFacesPerCell = 6;
  localIndex const numNodesPerEdge = ( order+1 );
  localIndex const numNodesPerFace = ( order+1 )*( order+1 );
  localIndex const numNodesPerCell = ( order+1 )*( order+1 )*( order+1 );
  localIndex const numInternalNodesPerEdge = ( order-1 );
  localIndex const numInternalNodesPerFace = ( order-1 )*( order-1 );
  localIndex const numInternalNodesPerCell = ( order-1 )*( order-1 )*( order-1 );
  localIndex const numLocalVertices = source.m_nodeManager->size();
  localIndex const numLocalEdges = source.m_edgeManager->size();
  localIndex const numLocalFaces = source.m_faceManager->size();
  localIndex const numGlobalVertices = source.m_nodeManager->maxGlobalIndex() + 1;
  localIndex const numGlobalEdges = source.m_edgeManager->maxGlobalIndex() + 1;
  localIndex const numGlobalFaces = source.m_faceManager->maxGlobalIndex() + 1;
  localIndex n1 = 0;
  // localIndex n2 = 0;
  source.m_elementManager->forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion const & sourceSubRegion )
  {
   n1 += sourceSubRegion.size();
   localIndex maxRegionIndex = sourceSubRegion.maxGlobalIndex() + 1; 
   // n2 = std::max( n2, maxRegionIndex );
  } );
  localIndex const numLocalCells = n1;
  // localIndex const numGlobalCells = n2;
  ////////////////////////////////
  // Get the new number of nodes 
  ////////////////////////////////
  localIndex numLocalNodes = numLocalVertices 
                             + numLocalEdges * numInternalNodesPerEdge 
                             + numLocalFaces * numInternalNodesPerFace 
                             + numLocalCells * numInternalNodesPerCell; 
  GEOSX_LOG_RANK ("!!!! INFO !!!! numLocalVertices " << numLocalVertices );
  GEOSX_LOG_RANK ("!!!! INFO !!!! numLocalEdgess " << numLocalEdges );
  GEOSX_LOG_RANK ("!!!! INFO !!!! numLocalFaces " << numLocalFaces );
  GEOSX_LOG_RANK ("!!!! INFO !!!! numLocalCells " << numLocalCells );
  GEOSX_LOG_RANK ("!!!! INFO !!!! numGlobalVertices " << numGlobalVertices );
  GEOSX_LOG_RANK ("!!!! INFO !!!! numGlobalEdgess " << numGlobalEdges );
  GEOSX_LOG_RANK ("!!!! INFO !!!! numGlobalFaces " << numGlobalFaces );
  //GEOSX_LOG_RANK ("!!!! INFO !!!! numGlobalCells " << numGlobalCells );


  // initialize the node local to global map

  cellBlockManager.setNumNodes( numLocalNodes, order );

  arrayView1d< globalIndex const > const nodeLocalToGlobalSource = source.m_nodeManager->localToGlobalMap();
  arrayView1d< globalIndex > nodeLocalToGlobalNew = cellBlockManager.getNodeLocalToGlobal();

  //GEOSX_LOG_RANK ("!!!! INFO !!!! nodeLocalToGlobalSource = "<< nodeLocalToGlobalSource );
  //GEOSX_LOG_RANK ("!!!! INFO !!!! nodeLocalToGlobalNew = "<< nodeLocalToGlobalNew );

  // hash map that contains unique vertices and their indices for shared nodes.
  // a vertex is identified by 6 integers [i1 i2 i3 i4 a b], as follows
  // - nodes on a vertex v are identified by the vector [v -1 -1 -1 -1 -1]
  // - edge nodes are given by a linear interpolation between vertices v1 and v2, 'a' steps away from v1. 
  //   We assume that v1 < v2 and identify these nodes with [v1 v2 -1 -1 a -1].
  // - face nodes are given by a bilinear interpolation between edges v1-v2 and v3-v4 (v1-v4 and v2-v3 are the diagonals), with interpolation parameters 'a' and 'b'. 
  //   We assume that v1 is the smallest, and that v2 < v3. Then these nodes are identified with [v1 v2 v3 v4 -1 -1 -1 -1 a b -1]
  // - cell nodes are given as a trilinear interpolation between nodes v1--v8 in the usual lexicographical order of hex cell corners. 
  //   these nodes are not inserted in the hash map, since they are only visited once.
  std::unordered_map< std::array< localIndex, 6 >, localIndex, NodeKeyHasher, NodeKeyEqual > nodeIDs;

  // Create new nodes, with local and global IDs
  localIndex localNodeID = 0;
  for( localIndex iter_vertex=0; iter_vertex < numLocalVertices; iter_vertex++)
  {
    nodeLocalToGlobalNew[ localNodeID ] = nodeLocalToGlobalSource[ iter_vertex ];
    nodeIDs[ createNodeKey( iter_vertex ) ] = localNodeID;
    GEOSX_LOG_RANK ("!!!! INFO !!!! created node for vertex, node ID= "<<localNodeID );
    localNodeID++;
  }
  GEOSX_LOG_RANK ("!!!! INFO !!!! finished vertex nodes (" << localNodeID << " so far )" );




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
  //GEOSX_LOG_RANK ("!!!! INFO !!!! edgeLocalToGlobal = "<< edgeLocalToGlobal <<"; numInternalEdgeNodes="<<numInternalEdgeNodes);
  //GEOSX_LOG_RANK_0 ("!!!! INFO !!!! test="<<test);

  //
  // ---- initialize edge-to-node map ----
  //

  // get information from the source (base mesh-level) edge-to-node map
  arrayView2d< localIndex const > const & edgeToNodeMapSource = source.m_edgeManager->nodeList();
  arrayView2d< localIndex > edgeToNodeMapNew = cellBlockManager.getEdgeToNodes();
  m_edgeManager->nodeList().resize( numLocalEdges, numNodesPerEdge );

  // create / retrieve nodes on edges
  localIndex offset = numLocalVertices;
  for( localIndex iter_edge = 0; iter_edge < numLocalEdges; iter_edge++ )
  {
    localIndex v1 = source.m_edgeManager->nodeList()[ iter_edge ][ 0 ];
    localIndex v2 = source.m_edgeManager->nodeList()[ iter_edge ][ 1 ];
    for( int q=0;q<numNodesPerEdge; q++ )
    {
      localIndex nodeID;
      std::array< localIndex, 6 > nodeKey = createNodeKey( v1, v2, q, order );
      if( !nodeIDs.count( nodeKey ) == 1 )
      {
        // this is an internal edge node: create it
        nodeID = localNodeID;
        nodeIDs[ nodeKey ] = localNodeID;
        GEOSX_LOG_RANK ("!!!! INFO !!!! created node for edge, node ID= "<< nodeID );
        nodeLocalToGlobalNew[ nodeID ] = offset + edgeLocalToGlobal[ iter_edge ] * numInternalNodesPerEdge + (q - 1);
        localNodeID++;                 
      }
      else
      {
        nodeID = nodeIDs[ nodeKey ];
        GEOSX_LOG_RANK ("!!!! INFO !!!! found previous node, node ID= "<< nodeID );
      }
      edgeToNodeMapNew[ iter_edge ][ q ] = nodeID;
    }
  }
  GEOSX_LOG_RANK ("!!!! INFO !!!! finished edge nodes (" << localNodeID << " so far )" );

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
  // copy the faces-boundaryIndicator from source
  m_faceManager->getDomainBoundaryIndicator() = source.m_faceManager->getDomainBoundaryIndicator();

  arrayView1d< globalIndex const > const faceLocalToGlobal = m_faceManager->localToGlobalMap();
  for (localIndex iter_localToGlobalsize = 0; iter_localToGlobalsize < numLocalFaces; iter_localToGlobalsize++)
  {
    m_faceManager->localToGlobalMap()[iter_localToGlobalsize] = source.m_faceManager->localToGlobalMap()[iter_localToGlobalsize];
  }
  m_faceManager->constructGlobalToLocalMap();
  //GEOSX_LOG_RANK ("!!!! INFO !!!! faceLocalToGlobal = "<< faceLocalToGlobal <<"; numInternalFaceNodes="<<numInternalFaceNodes);

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
  offset = numLocalVertices + numLocalEdges * numInternalNodesPerEdge;
  for( localIndex iter_face = 0; iter_face < numLocalFaces; iter_face++ )
  {
    localIndex v1 = source.m_faceManager->nodeList()[ iter_face ][ 0 ];
    localIndex v2 = source.m_faceManager->nodeList()[ iter_face ][ 1 ];
    localIndex v3 = source.m_faceManager->nodeList()[ iter_face ][ 2 ];
    localIndex v4 = source.m_faceManager->nodeList()[ iter_face ][ 3 ];
    std::swap(v3, v4);
    for( int q1=0;q1<numNodesPerEdge; q1++ )
    {
      for( int q2=0;q2<numNodesPerEdge; q2++ )
      {
        localIndex nodeID;
        std::array< localIndex, 6 > nodeKey = createNodeKey( v1, v2, v3, v4, q1, q2, order );
        if( !nodeIDs.count( nodeKey ) == 1 )
        {
          // this is an internal face node: create it
          nodeID = localNodeID;
          nodeIDs[ nodeKey ] = localNodeID;
          GEOSX_LOG_RANK ("!!!! INFO !!!! created node for face, node ID= "<< nodeID );
          nodeLocalToGlobalNew[ nodeID ] = offset 
                                              + faceLocalToGlobal[ iter_face ] * numInternalNodesPerFace 
                                              + numInternalNodesPerEdge * ( q1 - 1 )
                                              + ( q2 - 1 );
          localNodeID++;                 
        }
        else
        {
          nodeID = nodeIDs[ nodeKey ];
          GEOSX_LOG_RANK ("!!!! INFO !!!! found previous node, node ID= "<< nodeID );
        }
        faceToNodeMapNew[ iter_face ][ q2 + q1*numNodesPerEdge ] = nodeID;
      }
    }
  }
  GEOSX_LOG_RANK ("!!!! INFO !!!! finished face nodes (" << localNodeID << " so far )" );

  GEOSX_LOG_RANK ("!!!! INFO !!!! faceToNodeMapSource = "<< faceToNodeMapSource ); 
  GEOSX_LOG_RANK ("!!!! INFO !!!! faceToNodeMapNew = "<< faceToNodeMapNew ); 


  //////////////////////////
  // Nodes
  //////////////////////////

  // the number of nodes per elements depends on the order number
  m_nodeManager->resize( numLocalNodes );

  Group & nodeSets = m_nodeManager->sets();
  SortedArray< localIndex > & allNodes  = nodeSets.registerWrapper< SortedArray< localIndex > >( string( "all" ) ).reference();
  allNodes.reserve( numLocalNodes );

  for( localIndex iter_nodes=0; iter_nodes< numLocalNodes; ++iter_nodes )
  {
    allNodes.insert( iter_nodes );
  }
  GEOSX_LOG_RANK ("!!!! INFO !!!! finished allnodes" );

  //GEOSX_LOG_RANK ("!!!! INFO !!!! refPosNew = "<< refPosNew );
  //GEOSX_LOG_RANK ("!!!! INFO !!!! refPosSource = "<< refPosSource );


  /////////////////////////
  // Elements 
  //////////////////////////
  // also assign node coordinates using trilinear interpolation in th elements
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const refPosSource = source.m_nodeManager->referencePosition();
  arrayView2d< real64, nodes::REFERENCE_POSITION_USD > refPosNew = cellBlockManager.getNodePositions();
  refPosNew.setValues< parallelHostPolicy >( -1.0 );
  GEOSX_LOG_RANK ("!!!! INFO !!!! retrieved refposnew" );

  real64 Xmesh[ numVerticesPerCell ][ 3 ] = { { } };
  real64 X[ 3 ] = { { } };
  array1d< real64 > glCoords = gaussLobattoPoints( order ); 
  localIndex elemMeshVertices[ numVerticesPerCell ] = { };
  offset = numLocalVertices + numLocalEdges * numInternalNodesPerEdge + numLocalFaces * numInternalNodesPerFace;
  std::array< localIndex, 6 > const nullKey = std::array< localIndex, 6 >{ -1, -1 ,-1, -1, -1, -1 };
  GEOSX_LOG_RANK ("!!!! INFO !!!! looping on regions" );
  source.m_elementManager->forElementRegions< CellElementRegion >( [&]( CellElementRegion const & sourceRegion )
  {
  GEOSX_LOG_RANK ("!!!! INFO !!!! in region" );
    // create element region with the same name as source element region "Region"
    CellElementRegion & region = *(dynamic_cast< CellElementRegion * >( m_elementManager->createChild( sourceRegion.getCatalogName(),
                                                                                                       sourceRegion.getName() ) ) );
    // add cell block to the new element region with the same name as cell block name from source element region
    region.addCellBlockNames( sourceRegion.getCellBlockNames() );

    sourceRegion.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion const & sourceSubRegion )
    {
  GEOSX_LOG_RANK ("!!!! INFO !!!! in subregion" );
 
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

      // copy the elements-to-edgs map from source
      arrayView2d< localIndex const > const & elemsToEdgesSource = sourceSubRegion.edgeList();
      array2d< localIndex > & elemsToEdgesNew = newSubRegion.edgeList();
      elemsToEdgesNew = elemsToEdgesSource;


      // initialize the elements-to-nodes map 
      arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodesSource = sourceSubRegion.nodeList().toViewConst();
      arrayView2d< localIndex, cells::NODE_MAP_USD > elemsToNodesNew = cellBlockManager.getElemToNodes(newSubRegion.getName(), numNodesPerCell);
      //GEOSX_LOG_RANK ("!!!! INFO !!!! elemsToNodesSource="<<elemsToNodesSource); 

      arrayView1d< globalIndex const > const elementLocalToGlobal = sourceSubRegion.localToGlobalMap();
      for (localIndex iter_localToGlobalsize = 0; iter_localToGlobalsize < sourceSubRegion.localToGlobalMap().size(); iter_localToGlobalsize++)
      {
        newSubRegion.localToGlobalMap()[iter_localToGlobalsize] = sourceSubRegion.localToGlobalMap()[iter_localToGlobalsize];
      }
      newSubRegion.constructGlobalToLocalMap();
      //GEOSX_LOG_RANK ("!!!! INFO !!!! elementLocalToGlobal="<<elementLocalToGlobal); 


      // then loop through all the elements and assign the globalID according to the globalID of the Element 
      // and insert the new local to global ID ( for the internal nodes of elements ) into the nodeLocalToGlobal
      // retrieve finite element type 
      for( localIndex iter_elem = 0; iter_elem < numLocalCells; ++iter_elem )
      {
        for( localIndex iter_vertex = 0; iter_vertex < numVerticesPerCell; iter_vertex++ )
        {
          elemMeshVertices[ iter_vertex ] = elemsToNodesSource[ iter_elem ][ iter_vertex ];
          for( int i =0; i < 3; i++ )
          {
            Xmesh[ iter_vertex ][ i ] = refPosSource[ elemMeshVertices[ iter_vertex ] ][ i ];
          }
        }
        GEOSX_LOG_RANK ("!!!! INFO !!!! cell vertices: " << elemMeshVertices[ 0 ] << " "<< elemMeshVertices[ 1 ] << " " << elemMeshVertices[ 2 ] << " " << elemMeshVertices[ 3 ] << " " << elemMeshVertices[ 4 ] << " " << elemMeshVertices[ 5 ]  << " " << elemMeshVertices[ 6 ] << " " << elemMeshVertices[ 7 ]);

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
          GEOSX_LOG_RANK ("!!!! INFO !!!! created node for cell, node ID= "<< nodeID );
            nodeLocalToGlobalNew[ nodeID ] = offset
                                             + elementLocalToGlobal[ iter_elem ] * numInternalNodesPerCell
                                             + numInternalNodesPerFace * (q1 - 1)
                                             + numInternalNodesPerEdge * (q2 - 1)
                                             + (q3 - 1);
            localNodeID++;
          }
          else
          {
            nodeID = nodeIDs[ nodeKey ]; 
          GEOSX_LOG_RANK ("!!!! INFO !!!! found previous node, node ID= "<< nodeID );
          }
          for( int i=0; i<3; i++ )
          {
            refPosNew( nodeID, i ) = X[ i ];
          }
          elemsToNodesNew[ iter_elem ][ q ] = nodeID;
        }
      }
  GEOSX_LOG_RANK ("!!!! INFO !!!! finished cell nodes (" << localNodeID << " so far )" );
      
       GEOSX_LOG_RANK ("!!!! INFO !!!! faceToNodeMapNew = "<< faceToNodeMapNew);


       GEOSX_LOG_RANK ("!!!! INFO !!!! refPosNew = "<< refPosNew );
       GEOSX_LOG_RANK ("!!!! INFO !!!! nodeLocalToGlobalSource = "<< nodeLocalToGlobalSource );
       GEOSX_LOG_RANK ("!!!! INFO !!!! nodeLocalToGlobalNew = "<< nodeLocalToGlobalNew );
       GEOSX_LOG_RANK ("!!!! INFO !!!! faceToNodeMapSource = "<< faceToNodeMapSource);
       GEOSX_LOG_RANK ("!!!! INFO !!!! faceToNodeMapNew = "<< faceToNodeMapNew);
       GEOSX_LOG_RANK ("!!!! INFO !!!! edgeToNodeMapSource = "<< edgeToNodeMapSource);
       GEOSX_LOG_RANK ("!!!! INFO !!!! edgeToNodeMapNew = "<< edgeToNodeMapNew);
       GEOSX_LOG_RANK ("!!!! INFO !!!! elemsToNodesNew = "<< elemsToNodesNew);

       GEOSX_LOG_RANK ("!!!! INFO !!!! elemsToNodesSource = "<< elemsToNodesSource);
       GEOSX_LOG_RANK ("!!!! INFO !!!! refPosNew = "<< refPosNew );

      //GEOSX_LOG_RANK_0 ("!!!! INFO !!!! elemsToNodesNew = "<< elemsToNodesNew);
      //GEOSX_LOG_RANK_0 ("!!!! INFO !!!! faceToNodeMapNew = "<< faceToNodeMapNew);
    } );
  } );
}
}

MeshLevel::MeshLevel( string const & name,
                      Group * const parent,
                      MeshLevel const & source,
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
