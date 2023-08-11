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
 * @file FaceElementSubRegion.cpp
 */

#include "FaceElementSubRegion.hpp"
#include "common/GEOS_RAJA_Interface.hpp"

#include "NodeManager.hpp"
#include "MeshLevel.hpp"
#include "BufferOps.hpp"
#include "LvArray/src/tensorOps.hpp"
#include "common/MpiWrapper.hpp"

namespace geos
{
using namespace dataRepository;


FaceElementSubRegion::FaceElementSubRegion( string const & name,
                                            dataRepository::Group * const parent ):
  SurfaceElementSubRegion( name, parent ),
  m_unmappedGlobalIndicesInToEdges(),
  m_unmappedGlobalIndicesInToFaces(),
  m_newFaceElements(),
  m_toFacesRelation()
{
  m_elementType = ElementType::Hexahedron;

  registerWrapper( viewKeyStruct::dNdXString(), &m_dNdX ).setSizedFromParent( 1 ).reference().resizeDimension< 3 >( 3 );

  registerWrapper( viewKeyStruct::detJString(), &m_detJ ).setSizedFromParent( 1 ).reference();

  registerWrapper( viewKeyStruct::faceListString(), &m_toFacesRelation ).
    setDescription( "Map to the faces attached to each FaceElement." ).
    reference().resize( 0, 2 );

  registerWrapper( viewKeyStruct::edgesTofractureConnectorsEdgesString(), &m_edgesTo2dFaces ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "A map of edge local indices to the fracture connector local indices." ).
    setSizedFromParent( 0 );

  registerWrapper( viewKeyStruct::fractureConnectorEdgesToEdgesString(), &m_2dFaceToEdge ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "A map of fracture connector local indices to edge local indices." ).
    setSizedFromParent( 0 );

  registerWrapper( viewKeyStruct::fractureConnectorsEdgesToFaceElementsIndexString(), &m_2dFaceTo2dElems ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "A map of fracture connector local indices face element local indices" ).
    setSizedFromParent( 0 );

#ifdef GEOSX_USE_SEPARATION_COEFFICIENT
  registerWrapper( viewKeyStruct::separationCoeffString(), &m_separationCoefficient ).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( dataRepository::PlotLevel::LEVEL_1 ).
    setDescription( "Scalar indicator of level of separation for a fracturing face." );
#endif

  excludeWrappersFromPacking( { viewKeyStruct::faceListString() } );

  m_2dElemToElems.resize( 0, 2 );

  m_numNodesPerElement = 8;
}

void FaceElementSubRegion::copyFromCellBlock( FaceBlockABC const & faceBlock )
{
  localIndex const num2dElements = faceBlock.num2dElements();
  resize( faceBlock.num2dElements() );

  m_toNodesRelation.base() = faceBlock.get2dElemToNodes();
  m_toEdgesRelation.base() = faceBlock.get2dElemToEdges();

  // `FaceBlockABC` is designed to be heterogeneous.
  // `FaceElementSubRegion` inherits from `ElementSubRegionBase` which is meant to be homogeneous.
  // But `FaceElementSubRegion` is sometimes used as an heterogeneous sub region,
  // which emphasizes the need of a refactoring.
  // In the meantime, we try to fill the face block into the sub region and hope for the best...
  {
    auto const hack = []( integer allSizes ) -> ElementType
    {
      GEOS_LOG_RANK( "All sizes : " << allSizes );
      switch( allSizes )
      {
        case 3:
        case 6:
          return ElementType::Wedge;
        case 4:
        case 8:
        case 0:
          return ElementType::Hexahedron;
        default:
          GEOS_ERROR( "Unsupported type of elements during the face element sub region creation." );
          return {};
      }
    };

    // Checking if all the 2d elements are homogeneous.
    // We rely on the number of nodes for each element to find out.
    std::vector< integer > sizes( num2dElements );
    for( int i = 0; i < num2dElements; ++i )
    {
      sizes[i] = m_toNodesRelation[i].size();
    }
    std::set< integer > const s( sizes.cbegin(), sizes.cend() );

    if( s.size() > 1 )
    {
      // If we have found that the input face block contains 2d elements of different types,
      // we inform the used that the situation may be at risk.
      // (We're storing the face block in a homogeneous container while it's actually heterogeneous).
      GEOS_WARNING( "Heterogeneous face element sub region found and stored as homogeneous. Use at your own risk." );
    }

    auto const it = std::max_element( s.cbegin(), s.cend() );
    integer const maxSize = *it;
    m_elementType = hack( maxSize );
    m_numNodesPerElement = maxSize;
  }

  // The `m_2dElemToElems` mappings involves element, sub regions and regions indices.
  // We store the element indices that are correct.
  // But we only have access to the cell block indices, not the sub regions indices.
  // Temporarily, and also because they share the same dimensions,
  // we store the cell block mapping at the sub region mapping location.
  // It will later be transformed into a sub regions mapping.  // Last, we fill the regions mapping with dummy -1 values that should all be replaced eventually.
  auto const elem2dToElems = faceBlock.get2dElemToElems();
  m_2dElemToElems.resize( num2dElements, 2 );
  for( int i = 0; i < num2dElements; ++i )
  {
    for( localIndex const & j: elem2dToElems.toCellIndex[i] )
    {
      m_2dElemToElems.m_toElementIndex.emplaceBack( i, j );
    }
    for( localIndex const & j: elem2dToElems.toBlockIndex[i] )
    {
      m_2dElemToElems.m_toElementSubRegion.emplaceBack( i, j );
    }
  }

  m_toFacesRelation.base() = faceBlock.get2dElemToFaces();

  m_2dFaceToEdge = faceBlock.get2dFaceToEdge();
  m_2dFaceTo2dElems = faceBlock.get2dFaceTo2dElems();

  m_localToGlobalMap = faceBlock.localToGlobalMap();
  this->constructGlobalToLocalMap();

  for( int i = 0; i < faceBlock.num2dFaces(); ++i )
  {
    m_recalculateConnectionsFor2dFaces.insert( i );
  }

  for( localIndex i = 0; i < faceBlock.num2dElements(); ++i )
  {
    m_newFaceElements.insert( i );
  }

  m_collocatedNodes = faceBlock.getCollocatedNodes();
  m_collocatedNodesOf2dElems = faceBlock.getCollocatedNodesOf2dElems();

  // TODO We still need to be able to import fields on the FaceElementSubRegion.
}

void FaceElementSubRegion::setupRelatedObjectsInRelations( MeshLevel const & mesh )
{
  this->m_toNodesRelation.setRelatedObject( mesh.getNodeManager() );
  this->m_toEdgesRelation.setRelatedObject( mesh.getEdgeManager() );
  this->m_toFacesRelation.setRelatedObject( mesh.getFaceManager() );
}

void FaceElementSubRegion::calculateSingleElementGeometricQuantities( localIndex const k,
                                                                      arrayView1d< real64 const > const & faceArea )
{
  m_elementArea[k] = faceArea[ m_toFacesRelation[k][0] ];
  m_elementVolume[k] = m_elementAperture[k] * faceArea[m_toFacesRelation[k][0]];
}

void FaceElementSubRegion::calculateElementGeometricQuantities( NodeManager const & GEOS_UNUSED_PARAM( nodeManager ),
                                                                FaceManager const & faceManager )
{
  arrayView1d< real64 const > const & faceArea = faceManager.faceArea();

  forAll< parallelHostPolicy >( this->size(), [=] ( localIndex const k )
  {
    calculateSingleElementGeometricQuantities( k, faceArea );
  } );
}



localIndex FaceElementSubRegion::packUpDownMapsSize( arrayView1d< localIndex const > const & packList ) const
{
  buffer_unit_type * junk = nullptr;
  return packUpDownMapsImpl< false >( junk, packList );
}

localIndex FaceElementSubRegion::packUpDownMaps( buffer_unit_type * & buffer,
                                                 arrayView1d< localIndex const > const & packList ) const
{
  return packUpDownMapsImpl< true >( buffer, packList );
}


template< bool DO_PACKING >
localIndex FaceElementSubRegion::packUpDownMapsImpl( buffer_unit_type * & buffer,
                                                     arrayView1d< localIndex const > const & packList ) const
{
  arrayView1d< globalIndex const > const localToGlobal = this->localToGlobalMap();
  arrayView1d< globalIndex const > const nodeLocalToGlobal = m_toNodesRelation.relatedObjectLocalToGlobal();
  arrayView1d< globalIndex const > const edgeLocalToGlobal = m_toEdgesRelation.relatedObjectLocalToGlobal();
  arrayView1d< globalIndex const > const faceLocalToGlobal = m_toFacesRelation.relatedObjectLocalToGlobal();

  localIndex packedSize = bufferOps::Pack< DO_PACKING >( buffer, string( viewKeyStruct::nodeListString() ) );
  packedSize += bufferOps::Pack< DO_PACKING >( buffer,
                                               m_toNodesRelation.base().toViewConst(),
                                               m_unmappedGlobalIndicesInToNodes,
                                               packList,
                                               localToGlobal,
                                               nodeLocalToGlobal );

  packedSize += bufferOps::Pack< DO_PACKING >( buffer, string( viewKeyStruct::edgeListString() ) );
  packedSize += bufferOps::Pack< DO_PACKING >( buffer,
                                               m_toEdgesRelation.base().toViewConst(),
                                               m_unmappedGlobalIndicesInToEdges,
                                               packList,
                                               localToGlobal,
                                               edgeLocalToGlobal );

  packedSize += bufferOps::Pack< DO_PACKING >( buffer, string( viewKeyStruct::faceListString() ) );
  packedSize += bufferOps::Pack< DO_PACKING >( buffer,
                                               m_toFacesRelation.toViewConst(),
                                               m_unmappedGlobalIndicesInToFaces,
                                               packList,
                                               localToGlobal,
                                               faceLocalToGlobal );

  packedSize += bufferOps::Pack< DO_PACKING >( buffer, string( viewKeyStruct::surfaceElementsToCellRegionsString() ) );
  packedSize += bufferOps::Pack< DO_PACKING >( buffer,
                                               m_2dElemToElems,
                                               packList,
                                               m_2dElemToElems.getElementRegionManager() );

//  packedSize += bufferOps::Pack< DO_PACKING >( buffer, string( viewKeyStruct::collocatedNodesOf2dElemString() ) );
//  packedSize += bufferOps::Pack< DO_PACKING >( buffer, m_collocatedNodesOf2dElems );

  packedSize += bufferOps::Pack< DO_PACKING >( buffer, string( viewKeyStruct::collocatedNodesString() ) );
  packedSize += bufferOps::Pack< DO_PACKING >( buffer, m_collocatedNodes );

//  packedSize += bufferOps::Pack< DO_PACKING >( buffer, string( "missingNodes" ) );
//  packedSize += bufferOps::Pack< DO_PACKING >( buffer, m_missingNodes );

  return packedSize;
}



localIndex FaceElementSubRegion::unpackUpDownMaps( buffer_unit_type const * & buffer,
                                                   localIndex_array & packList,
                                                   bool const overwriteUpMaps,
                                                   bool const GEOS_UNUSED_PARAM( overwriteDownMaps ) )
{
  localIndex unPackedSize = 0;

  string nodeListString;
  unPackedSize += bufferOps::Unpack( buffer, nodeListString );
  GEOS_ERROR_IF_NE( nodeListString, viewKeyStruct::nodeListString() );

  unPackedSize += bufferOps::Unpack( buffer,
                                     m_toNodesRelation,
                                     packList,
                                     m_unmappedGlobalIndicesInToNodes,
                                     this->globalToLocalMap(),
                                     m_toNodesRelation.relatedObjectGlobalToLocal() );


  string edgeListString;
  unPackedSize += bufferOps::Unpack( buffer, edgeListString );
  GEOS_ERROR_IF_NE( edgeListString, viewKeyStruct::edgeListString() );

  unPackedSize += bufferOps::Unpack( buffer,
                                     m_toEdgesRelation,
                                     packList,
                                     m_unmappedGlobalIndicesInToEdges,
                                     this->globalToLocalMap(),
                                     m_toEdgesRelation.relatedObjectGlobalToLocal() );

  string faceListString;
  unPackedSize += bufferOps::Unpack( buffer, faceListString );
  GEOS_ERROR_IF_NE( faceListString, viewKeyStruct::faceListString() );

  unPackedSize += bufferOps::Unpack( buffer,
                                     m_toFacesRelation,
                                     packList,
                                     m_unmappedGlobalIndicesInToFaces,
                                     this->globalToLocalMap(),
                                     m_toFacesRelation.relatedObjectGlobalToLocal() );

  string elementListString;
  unPackedSize += bufferOps::Unpack( buffer, elementListString );
  GEOS_ERROR_IF_NE( elementListString, viewKeyStruct::surfaceElementsToCellRegionsString() );

  unPackedSize += bufferOps::Unpack( buffer,
                                     m_2dElemToElems,
                                     packList.toViewConst(),
                                     m_2dElemToElems.getElementRegionManager(),
                                     overwriteUpMaps );

  string collocatedNodesString;
  unPackedSize += bufferOps::Unpack( buffer, collocatedNodesString );
  GEOS_ERROR_IF_NE( collocatedNodesString, viewKeyStruct::collocatedNodesString() );
  m_otherCollocatedNodes.push_back( ArrayOfArrays< globalIndex >{} );
  unPackedSize += bufferOps::Unpack( buffer, m_otherCollocatedNodes.back() );

  return unPackedSize;
}

void FaceElementSubRegion::fixUpDownMaps( bool const clearIfUnmapped )
{
  ObjectManagerBase::fixUpDownMaps( m_toNodesRelation,
                                    m_unmappedGlobalIndicesInToNodes,
                                    clearIfUnmapped );

  ObjectManagerBase::fixUpDownMaps( m_toEdgesRelation,
                                    m_unmappedGlobalIndicesInToEdges,
                                    clearIfUnmapped );

  ObjectManagerBase::fixUpDownMaps( m_toFacesRelation,
                                    m_unmappedGlobalIndicesInToFaces,
                                    clearIfUnmapped );
}

void FaceElementSubRegion::fixSecondaryMappings( NodeManager const & nodeManager,
                                                 EdgeManager const & edgeManager,
                                                 FaceManager const & faceManager,
                                                 ElementRegionManager const & elemManager )
{
  // Here I can fix the other mappings which are not properly defined...
  localIndex const num2dElems = this->size();

  std::map< globalIndex, globalIndex > referenceCollocatedNodes;  // Take the minimum duplicated node for each bucket.
  {
    std::set< std::set< globalIndex > > mergedCollocatedNodes;
    m_otherCollocatedNodes.push_back( m_collocatedNodes );
    for( ArrayOfArrays< globalIndex > const & dns: m_otherCollocatedNodes )
    {
      for( int i = 0; i < dns.size(); ++i )
      {
        std::set< globalIndex > tmp( dns[i].begin(), dns[i].end() );
        mergedCollocatedNodes.insert( tmp );
      }
    }

    for( std::set< globalIndex > const & collocatedNodes: mergedCollocatedNodes )
    {
      globalIndex const & ref = *std::min_element( collocatedNodes.cbegin(), collocatedNodes.cend() );
      for( globalIndex const & n: collocatedNodes )
      {
        referenceCollocatedNodes[n] = ref;
      }
    }
  }

  // The concept of 2d face is deeply local to the `FaceElementSubRegion`.
  // Therefore, it's not exchanged during the ghosting process.
  // Furthermore, all the information is there to reconstruct the mappings related to the 2d faces.

  // To recreate the mappings, the first step is to associate a `2f face index` based on the global indices.
  // But let's be precise: some the global indices of the edges are different while the
  // Then we can resize all the mappings related to the `2d faces`.
  ArrayOfArraysView< localIndex const > const elem2dToEdges = m_toEdgesRelation.base().toViewConst();
  std::set< localIndex > edges;  // Will include twin edges.
  for( int i = 0; i < elem2dToEdges.size(); ++i )
  {
    for( localIndex const & j: elem2dToEdges[i] )
    {
      edges.insert( j );
    }
  }

  auto const & edgeToNodes = edgeManager.nodeList();
  std::map< std::pair< globalIndex, globalIndex >, localIndex > edgesIds;
  for( localIndex const & edge: edges )
  {
    auto const nodes = edgeToNodes[edge];
    GEOS_ASSERT_EQ( nodes.size(), 2 );
    std::pair< globalIndex, globalIndex > const pg{ nodeManager.localToGlobalMap()[nodes[0]], nodeManager.localToGlobalMap()[nodes[1]] };
    edgesIds[pg] = edge;
  }

  std::map< std::pair< globalIndex, globalIndex >, std::set< localIndex > > uniqueEdgeIds;
  for( auto const & p: edgesIds )
  {
    std::pair< globalIndex, globalIndex > const & nodes = p.first;
    localIndex const & edge = p.second;

    auto it0 = referenceCollocatedNodes.find( nodes.first );
    globalIndex const n0 = it0 != referenceCollocatedNodes.cend() ? it0->second : nodes.first;

    auto it1 = referenceCollocatedNodes.find( nodes.second );
    globalIndex const n1 = it1 != referenceCollocatedNodes.cend() ? it1->second : nodes.second;

    std::pair< globalIndex, globalIndex > const edgeHash = std::minmax( n0, n1 );
    uniqueEdgeIds[edgeHash].insert( edge );
//    try
//    {
//      std::pair< globalIndex, globalIndex > const edgeHash = std::minmax( referenceCollocatedNodes.at( nodes.first ), referenceCollocatedNodes.at( nodes.second ) );
//      uniqueEdgeIds[edgeHash].insert( edge );
//    }
//    catch( std::exception & e )
//    {
//      GEOS_LOG_RANK( "excp: " << e.what() );
//    }
  }

  std::map< localIndex, localIndex > duplicatedEdges;

  std::size_t const num2dFaces = uniqueEdgeIds.size();
  arrayView1d< integer const > edgeGhostRanks = edgeManager.ghostRank().toViewConst();
  m_2dFaceToEdge.clear();
  m_2dFaceToEdge.reserve( num2dFaces );
  auto cmp = [=]( int e,
                  int f ) -> bool
  {
    int const re = edgeGhostRanks[e] < 0 ? 0 : 1;
    int const rf = edgeGhostRanks[f] < 0 ? 0 : 1;
    return std::tie( re, e ) < std::tie( rf, f );
  };
  std::map< localIndex, localIndex > rejectedEdges;  // Mapping from the removed edge to its replacement.
  for( auto const & p: uniqueEdgeIds )
  {
    std::set< localIndex > const & input = p.second;
    localIndex const e = *std::min_element( input.cbegin(), input.cend(), cmp );
    std::set< localIndex > rejectedEdgeIds( input );
    rejectedEdgeIds.erase( e );
    m_2dFaceToEdge.emplace_back( e );
    for( localIndex const & re: rejectedEdgeIds )
    {
      rejectedEdges[re] = e;
    }
  }

  // map< localIndex, localIndex >  m_edgesTo2dFaces
  // Simple map inversion
  m_edgesTo2dFaces.clear();
  for( std::size_t i = 0; i < num2dFaces; ++i )
  {
    m_edgesTo2dFaces[m_2dFaceToEdge[i]] = i;
  }

  m_newFaceElements.clear();
  m_newFaceElements.reserve( num2dElems );
  for( localIndex i = 0; i < num2dElems; ++i )
  {
    m_newFaceElements.insert( i );
  }
  m_recalculateConnectionsFor2dFaces.clear();
  m_recalculateConnectionsFor2dFaces.reserve( num2dFaces );  // TODO variable...
  for( std::size_t i = 0; i < num2dFaces; ++i )
  {
    m_recalculateConnectionsFor2dFaces.insert( i );
  }

//  ArrayOfArrays< localIndex > tmp;
  std::vector< std::vector< localIndex > > tmp( num2dFaces );
  for( auto i = 0; i < num2dElems; ++i )
  {
    for( auto const & e: elem2dToEdges[i] )
    {
      auto const fIt = m_edgesTo2dFaces.find( e );
      if( fIt != m_edgesTo2dFaces.cend() )
      {
        localIndex const & f = fIt->second;
        tmp[f].push_back( i );
      }
      else
      {
        tmp[m_edgesTo2dFaces.at( rejectedEdges.at( e ) )].push_back( i );
      }
    }
  }
  {
    std::vector< localIndex > sizes;
    sizes.reserve( tmp.size() );
    for( std::vector< localIndex > const & t: tmp )
    {
      sizes.push_back( t.size() );
    }
    for( auto i = 0; i < m_2dFaceTo2dElems.size(); ++i )
    {
      m_2dFaceTo2dElems.clearArray( i );
    }
    m_2dFaceTo2dElems.resizeFromCapacities< serialPolicy >( sizes.size(), sizes.data() );

    for( std::size_t i = 0; i < tmp.size(); ++i )
    {
      for( std::size_t j = 0; j < tmp[i].size(); ++j )
      {
        m_2dFaceTo2dElems.emplaceBack( i, tmp[i][j] );
      }
    }
  }

  // When a fracture element has only one neighbor, let's try to find the other one.
  // Reconnecting the elements from the side of the fracture.
  struct ElemPath
  {
    localIndex er;
    localIndex esr;
    localIndex ei;
    localIndex face;
    std::vector< localIndex > nodes;

    bool operator<( ElemPath const & other ) const
    {
      return std::tie( er, esr, ei, face, nodes ) < std::tie( other.er, other.esr, other.ei, other.face, other.nodes );
    }

//    string toString() const
//    {
//      string n = "[" + stringutilities::join( nodes, ", " ) + "]";
//      std::vector< localIndex > tmp{
//        er, esr, ei, face
//      };
//      return "[" + stringutilities::join( tmp, ", " ) + ", " + n + "]";
//    }
  };

  std::map< std::set< globalIndex >, std::set< ElemPath > > faceNodesToElems;  // TODO so bad: only consider some candidate elements.
  auto const buildFaceNodesToElems = [&]( localIndex const er,
                                          localIndex const esr,
                                          ElementRegionBase const & region,
                                          CellElementSubRegion const & subRegion )
  {
    auto const & elemToFaces = subRegion.faceList().base();
    auto const & faceToNodes = faceManager.nodeList();
    auto const & isOnBoundary = subRegion.getDomainBoundaryIndicator();
    for( localIndex ei = 0; ei < elemToFaces.size( 0 ); ++ei )
    {
//      if( not isOnBoundary[ei] )
//      { continue; }

      for( auto const & face: elemToFaces[ei] )
      {
        std::set< globalIndex > nodesOfFace;  // Signature of the face nodes.
        for( localIndex const & n: faceToNodes[face] )
        {
          auto const it = referenceCollocatedNodes.find( nodeManager.localToGlobalMap()[n] );
          if( it != referenceCollocatedNodes.cend() )
          {
            nodesOfFace.insert( it->second );
          }
          // No else: not finding is legit. Maybe `break`?
        }
        if( nodesOfFace.size() == LvArray::integerConversion< std::size_t >( faceToNodes[face].size() ) )
        {
          std::vector< localIndex > const nodes( faceToNodes[face].begin(), faceToNodes[face].end() );
          ElemPath const path = ElemPath{ er, esr, ei, face, nodes };
          faceNodesToElems[nodesOfFace].insert( path );
        }
      }
    }
  };
  elemManager.forElementSubRegionsComplete< CellElementSubRegion >( buildFaceNodesToElems );

//  {
//    string msg;
//    for( auto const p: faceNodesToElems )
//    {
//      using stringutilities::join;
//      string key = join( p.first, ", " );
//      std::vector <string> elems;
//      for( ElemPath const & ep: p.second )
//      {
//        elems.push_back( ep.toString() );
//      }
//      msg += "{" + key + "} -> {" + stringutilities::join( elems, ", " ) + "}; ";
//    }
//    GEOS_LOG_RANK( msg );
//  }

  std::vector< localIndex > misMatches;
  for( int e2d = 0; e2d < num2dElems; ++e2d )
  {
    if( m_2dElemToElems.m_toElementIndex.sizeOfArray( e2d ) >= 2 )
    { continue; }
//    else
//    {
//      GEOS_LOG_RANK( "working on e2d: " << e2d );
//    }

    std::set< globalIndex > refNodes;  // TODO think about those ref node twice.
    if( m_toNodesRelation[e2d].size() != 0 )
    {
      for( localIndex const & n: m_toNodesRelation[e2d] )  // TODO Refatctor this branch into something neat...
      {
        globalIndex const & gn = nodeManager.localToGlobalMap()[n];
        auto const it = referenceCollocatedNodes.find( gn );
        if( it == referenceCollocatedNodes.cend() )
        {
          GEOS_LOG_RANK( "Could not find node " << gn << " in the duplicated nodes dict 1." );
        }
        else
        {
          refNodes.insert( it->second );
        }
      }
    }
    else if( m_ghostRank[e2d] < 0 )
    {
      for( globalIndex const & gn: m_collocatedNodesOf2dElems[e2d] )
      {
        auto const it = referenceCollocatedNodes.find( gn );
        if( it == referenceCollocatedNodes.cend() )
        {
          GEOS_LOG_RANK( "Could not find node " << gn << " in the duplicated nodes dict 2." );
        }
        else
        {
          refNodes.insert( it->second );
        }
      }
    }

    auto const match = faceNodesToElems.find( refNodes );
    if( match != faceNodesToElems.cend() )
    {
//      GEOS_LOG_RANK( "Found a match for " << e2d );
      bool found = false;
      for( ElemPath const & path: match->second )
      {
        if( m_2dElemToElems.m_toElementIndex.sizeOfArray( e2d ) == 0 || m_2dElemToElems.m_toElementIndex[e2d][0] != path.ei )
        {
          found = true;
//          GEOS_LOG_RANK( "Found a correction for " << e2d );
          m_2dElemToElems.m_toElementRegion.emplaceBack( e2d, path.er );
          m_2dElemToElems.m_toElementSubRegion.emplaceBack( e2d, path.esr );
          m_2dElemToElems.m_toElementIndex.emplaceBack( e2d, path.ei );
          m_toFacesRelation.emplaceBack( e2d, path.face );
          for( localIndex const & n: path.nodes )
          {
            m_toNodesRelation.emplaceBack( e2d, n );
          }
        }
      }
      if( !found && m_ghostRank[e2d] < 0 )
      {
        misMatches.emplace_back(e2d );
//        GEOS_LOG_RANK( "ERROR " << stringutilities::join( refNodes, ", " ) );
        GEOS_ERROR( "  -->> match not validated " << stringutilities::join( refNodes, ", " ) << ", ghostRank = " << m_ghostRank[e2d] );
      }
    }
  }

  GEOS_ERROR_IF( !misMatches.empty(),
                 "Fracture " << this->getName() << " has elements {" << stringutilities::join( misMatches, ", " ) << "} without two neighbors." );

  // Checking that each face has two neighboring elements.
  {
    std::vector< localIndex > isolatedFractureElements;
    for( int e2d = 0; e2d < num2dElems; ++e2d )
    {
      if( m_2dElemToElems.m_toElementIndex.sizeOfArray( e2d ) < 2 && m_ghostRank[e2d] < 0 )
      {
        isolatedFractureElements.push_back( e2d );
      }
    }
    GEOS_ERROR_IF( !isolatedFractureElements.empty(),
                   "Fracture " << this->getName() << " has elements {" << stringutilities::join( isolatedFractureElements, ", " ) << "} with less than two neighbors." );
  }

  // TODO clear all useless mappings!
}

void FaceElementSubRegion::inheritGhostRankFromParentFace( FaceManager const & faceManager,
                                                           std::set< localIndex > const & indices )
{
  arrayView1d< integer const > const & faceGhostRank = faceManager.ghostRank();
  for( localIndex const & index: indices )
  {
    m_ghostRank[index] = faceGhostRank[m_toFacesRelation[index][0]];
  }
}

std::set< globalIndex > FaceElementSubRegion::getMissingNodes( unordered_map< globalIndex, localIndex > const & g2l ) const
{
  std::set< globalIndex > missingNodes;
  for( int i = 0; i < m_collocatedNodes.size(); ++i )
  {
    for( globalIndex const & n:  m_collocatedNodes[i] )
    {
      auto const it = g2l.find( n );
      if( it == g2l.cend() )
      {
        missingNodes.insert( n );
      }
    }
  }

  return missingNodes;
}

} /* namespace geos */
