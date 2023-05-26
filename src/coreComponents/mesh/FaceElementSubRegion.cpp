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
    auto const f = []( integer allSizes ) -> ElementType
    {
      if( allSizes == 6 )
      {
        return ElementType::Wedge;
      }
      else if( allSizes == 8 )
      {
        return ElementType::Hexahedron;
      }
      else
      {
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
    m_elementType = f( maxSize );
    m_numNodesPerElement = maxSize;
  }

  // The `m_2dElemToElems` mappings involves element, sub regions and regions indices.
  // We store the element indices that are correct.
  // But we only have access to the cell block indices, not the sub regions indices.
  // Temporarily, and also because they share the same dimensions,
  // we store the cell block mapping at the sub region mapping location.
  // It will later be transformed into a sub regions mapping.  // Last, we fill the regions mapping with dummy -1 values that should all be replaced eventually.
  auto const elem2dToElems = faceBlock.get2dElemToElems();
  m_2dFaceTo2dElems.resize( num2dElements, 2 );
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

  m_duplicatedNodes = faceBlock.getDuplicatedNodes();

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


template< class T >
ArrayOfArrays< T > myConvert( array2d< T > const & input )
{
  ArrayOfArrays< T > result;
  auto const n = input.size( 0 );
  GEOS_ERROR_IF_NE_MSG( input.size( 1 ), 2, "Wrong size budy" );
  array1d< localIndex > sizes( n );
  for( auto i = 0; i < n; ++i )
  {
    sizes[i] = input( i, 1 ) == -1 ? 1 : 2;
  }
  result.template resizeFromCapacities< serialPolicy >( n, sizes.data() );
  for( auto i = 0; i < n; ++i )
  {
    result.emplaceBack(i, input( i, 0 ));
    if( input( i, 1 ) != -1 )
    { result.emplaceBack( i, input( i, 1 ) ); }
  }
  return result;
}


template< class T >
array2d< T > myConvert( ArrayOfArrays< T > const & input )
{
  auto const n = input.size();
  std::vector< int > sizes;
  sizes.reserve( n );
  for( int i = 0; i < n; ++i )
  {
    sizes.push_back( input.sizeOfArray( i ) );
  }
  int const m = *std::max_element( sizes.cbegin(), sizes.cend() );
  array2d< T > result( n, m );
  result.template setValues< serialPolicy >( -1 );
  for( int i = 0; i < n; ++i )
  {
    for( int j = 0; j < input.sizeOfArray( i ); ++j )
    {
      result( i, j ) = input( i, j );
    }
  }

  return result;
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
  if( MpiWrapper::commSize() == 1 )
  {
    return; // TODO I want a good reference for serial runs.
  }
  // Here I can fix the other mappings which are not properly defined...
  localIndex const num2dElems = this->size();

  std::map< globalIndex, globalIndex > const referenceDuplicatedNodes{
    { 5,   5 },
    { 132, 5 },
    { 11,  11 },
    { 138, 11 },
    { 17,  17 },
    { 144, 17 },
    { 23,  23 },
    { 150, 23 },
    { 29,  29 },
    { 156, 29 },
    { 35,  35 },
    { 162, 35 },
    { 204, 35 },
    { 41,  41 },
    { 210, 41 },
    { 47,  47 },
    { 216, 47 },
    { 53,  53 },
    { 222, 53 },
    { 59,  59 },
    { 228, 59 },
    { 65,  65 },
    { 234, 65 },
    { 71,  71 },
    { 168, 71 },
    { 77,  77 },
    { 174, 77 },
    { 83,  83 },
    { 180, 83 },
    { 89,  89 },
    { 186, 89 },
    { 95,  95 },
    { 192, 95 },
    { 101, 101 },
    { 198, 101 },
    { 240, 101 },
    { 107, 107 },
    { 246, 107 },
    { 113, 113 },
    { 252, 113 },
    { 119, 119 },
    { 258, 119 },
    { 125, 125 },
    { 264, 125 },
    { 131, 131 },
    { 270, 131 },
    { 163, 163 },
    { 205, 163 },
    { 164, 164 },
    { 206, 164 },
    { 165, 165 },
    { 207, 165 },
    { 166, 166 },
    { 208, 166 },
    { 167, 167 },
    { 209, 167 },
    { 199, 199 },
    { 241, 199 },
    { 200, 200 },
    { 242, 200 },
    { 201, 201 },
    { 243, 201 },
    { 202, 202 },
    { 244, 202 },
    { 203, 203 },
    { 245, 203 }
  };

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
//  std::map< std::pair< globalIndex, globalIndex >, std::set< localIndex > > edgesIds;
  for( localIndex const & edge: edges )
  {
    auto const nodes = edgeToNodes[edge];
    GEOS_ASSERT_EQ( nodes.size(), 2 );
    std::pair< globalIndex, globalIndex > const pg{ nodeManager.localToGlobalMap()[ nodes[0] ], nodeManager.localToGlobalMap()[ nodes[1] ] };
//    edgesIds[pg] = edge;
    edgesIds[pg] = edge;
  }

  std::map< std::pair< globalIndex, globalIndex >, std::set< localIndex > > uniqueEdgeIds;
  for( auto const & p: edgesIds )
  {
    std::pair< globalIndex, globalIndex > const & nodes = p.first;
    localIndex const & edge = p.second;
    std::pair< globalIndex, globalIndex > const edgeHash = std::minmax( referenceDuplicatedNodes.at( nodes.first ), referenceDuplicatedNodes.at( nodes.second ) );
    uniqueEdgeIds[edgeHash].insert( edge );
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
  std::map<localIndex, localIndex> rejectedEdges;  // Mapping from the removed edge to its replacement.
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
        m_2dFaceTo2dElems.emplaceBack(i, tmp[i][j]);
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

    bool operator<( ElemPath const & other ) const
    {
      return std::tie( er, esr, ei, face ) < std::tie( other.er, other.esr, other.ei, other.face );
    }
  };

  std::map< std::set< globalIndex >, std::set< ElemPath > > faceNodesToElems;  // TODO so bad: only consider some candidate elements.
  auto const ff = [&]( localIndex const er,
                       localIndex const esr,
                       ElementRegionBase const & region,
                       CellElementSubRegion const & subRegion )
  {
    auto const & elemToFaces = subRegion.faceList().base();
    auto const & faceToNodes = faceManager.nodeList();
    for( localIndex ei = 0; ei < elemToFaces.size( 0 ); ++ei )
    {
      for( auto const & face: elemToFaces[ei] )
      {
        std::set< globalIndex > nodesOfFace;
        for( localIndex const & n: faceToNodes[face] )
        {
          auto const it = referenceDuplicatedNodes.find( nodeManager.localToGlobalMap()[n] );
          if( it != referenceDuplicatedNodes.cend() )
          {
            nodesOfFace.insert( it->second );
          }
        }
        ElemPath const path = ElemPath{ er, esr, ei, face };
        faceNodesToElems[nodesOfFace].insert( path );
      }
    }
  };
  elemManager.forElementSubRegionsComplete< CellElementSubRegion >( ff );

  for( int e2d = 0; e2d < num2dElems; ++e2d )
  {
    if( m_2dElemToElems.m_toElementIndex.sizeOfArray( e2d ) >= 2 )
    { continue; }

    std::set< globalIndex > nodes;
    for( localIndex const & n: m_toNodesRelation[e2d] )
    {
      nodes.insert( referenceDuplicatedNodes.at( nodeManager.localToGlobalMap()[n] ) );
    }

    auto const match = faceNodesToElems.find( nodes );
    if( match != faceNodesToElems.cend() )
    {
      GEOS_LOG_RANK( "Found a match for " << e2d );
      for( ElemPath const & path: match->second )
      {
        if(  m_2dElemToElems.m_toElementIndex.sizeOfArray( e2d ) == 0 || m_2dElemToElems.m_toElementIndex[e2d][0] != path.ei )
        {
          GEOS_LOG_RANK( "Found a correction for " << e2d );
          m_2dElemToElems.m_toElementRegion.emplaceBack( e2d, path.er );
          m_2dElemToElems.m_toElementSubRegion.emplaceBack( e2d, path.esr );
          m_2dElemToElems.m_toElementIndex.emplaceBack( e2d, path.ei );
          m_toFacesRelation.emplaceBack( e2d, path.face );
          for( globalIndex const & gn: nodes )
          {
            m_toNodesRelation.emplaceBack( e2d, nodeManager.globalToLocalMap().at( gn ) );
          }
        }
      }
    }
  }
}

void FaceElementSubRegion::inheritGhostRankFromParentFace( FaceManager const & faceManager,
                                                           std::set< localIndex > const & indices )
{
  arrayView1d< integer const > const & faceGhostRank = faceManager.ghostRank();
  for( localIndex const & index : indices )
  {
    m_ghostRank[index] = faceGhostRank[ m_toFacesRelation[index][0] ];
  }
}

} /* namespace geos */
