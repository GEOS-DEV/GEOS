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

  registerWrapper( viewKeyStruct::elem2dToCollocatedNodesBucketsString(), &m_2dElemToCollocatedNodesBuckets ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "A map eventually containing all the collocated nodes." ).
    setSizedFromParent( 1 );

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
    auto const deduce3dElemType = [&]( integer maxSize ) -> ElementType
    {
      switch( maxSize )
      {
        case 3:
          return ElementType::Wedge;
        case 4:
          return ElementType::Hexahedron;
        case 5:
          return ElementType::Prism5;
        case 6:
          return ElementType::Prism6;
        case 7:
          return ElementType::Prism7;
        case 8:
          return ElementType::Prism8;
        case 9:
          return ElementType::Prism9;
        case 10:
          return ElementType::Prism10;
        case 11:
          return ElementType::Prism11;
        case 0:
          // In the case the fracture is empty (on this rank), then we default to hexahedron. Otherwise, there's something wrong
          GEOS_ERROR_IF_NE_MSG( num2dElements, 0, "Could not determine the element type of the fracture \"" << getName() << "\"." );
          return ElementType::Hexahedron;
        default:
          GEOS_ERROR( "Unsupported type of elements during the face element sub region creation." );
          return {};
      }
    };

    m_2dElemToCollocatedNodesBuckets = faceBlock.get2dElemsToCollocatedNodesBuckets();
    // Checking if all the 2d elements are homogeneous.
    // We rely on the number of nodes for each element to find out.
    std::vector< integer > numNodesPerElement( num2dElements );
    for( int i = 0; i < num2dElements; ++i )
    {
      numNodesPerElement[i] = m_2dElemToCollocatedNodesBuckets[i].size();
    }
    std::set< integer > const sizes( numNodesPerElement.cbegin(), numNodesPerElement.cend() );

    if( sizes.size() > 1 )
    {
      // If we have found that the input face block contains 2d elements of different types,
      // we inform the used that the situation may be at risk.
      // (We're storing the face block in a homogeneous container while it's actually heterogeneous).
      GEOS_WARNING( "Heterogeneous face element sub region found and stored as homogeneous. Use at your own risk." );
    }

    auto const it = std::max_element( sizes.cbegin(), sizes.cend() );
    integer const maxSize = it != sizes.cend() ? *it : 0;
    m_elementType = deduce3dElemType( maxSize );
    m_numNodesPerElement = maxSize;
  }

  // The `m_2dElemToElems` mappings involves element, sub regions and regions indices.
  // We store the element indices that are correct.
  // But we only have access to the cell block indices, not the sub regions indices.
  // Temporarily, and also because they share the same dimensions,
  // we store the cell block mapping at the sub region mapping location.
  // It will later be transformed into a sub regions mapping.
  // Last, we fill the regions mapping with dummy -1 values that should all be replaced eventually.
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

ElementType FaceElementSubRegion::getElementType( localIndex ei ) const
{
  // We try to get the element type information from the bucket.
  // If this information does not appear to be reliable,
  // let's fall back on the 2d elem to nodes.
  auto const sizeFromBuckets = m_2dElemToCollocatedNodesBuckets[ei].size();
  auto const size = sizeFromBuckets > 0 ? sizeFromBuckets : m_toNodesRelation[ei].size() / 2;
  switch( size )
  {
    case 3:
      return ElementType::Triangle;
    case 4:
      return ElementType::Quadrilateral;
    default:
      return ElementType::Polygon;
  }
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

  packedSize += bufferOps::Pack< DO_PACKING >( buffer, string( viewKeyStruct::elem2dToCollocatedNodesBucketsString() ) );
  packedSize += bufferOps::Pack< DO_PACKING >( buffer,
                                               m_2dElemToCollocatedNodesBuckets.toViewConst(),
                                               packList,
                                               localToGlobal );

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

  string elem2dToCollocatedNodesBucketsString;
  unPackedSize += bufferOps::Unpack( buffer, elem2dToCollocatedNodesBucketsString );
  GEOS_ERROR_IF_NE( elem2dToCollocatedNodesBucketsString, viewKeyStruct::elem2dToCollocatedNodesBucketsString() );
  unPackedSize += bufferOps::Unpack( buffer,
                                     m_2dElemToCollocatedNodesBuckets,
                                     packList,
                                     this->globalToLocalMap() );

  return unPackedSize;
}


/**
 * @brief The two mappings @p elem2dToElems3d and @p elem2dToFaces have to be consistent
 * in the sense that each @p face for a given index must "belong" to @p 3d @p element at the same index.
 * @param[in] fractureName The name of the fracture for which we're checking the mapping consistency.
 * @param[in] elem2dToElems3d A mapping.
 * @param[inout] elem2dToFaces This mapping will be corrected if needed to match @p elem2dToElems3d.
 */
void fixNeighborMappingsInconsistency( string const & fractureName,
                                       OrderedVariableToManyElementRelation const & elem2dToElems3d,
                                       FaceElementSubRegion::FaceMapType & elem2dToFaces )
{
  {
    localIndex const num2dElems = elem2dToFaces.size();
    // GEOS_ASSERT_EQ( elem2dToElems3d.size(), num2dElems );
    for( int e2d = 0; e2d < num2dElems; ++e2d )
    {
      std::set< localIndex > const sizes{
        elem2dToFaces[e2d].size(),
        elem2dToElems3d.m_toElementRegion[e2d].size(),
        elem2dToElems3d.m_toElementSubRegion[e2d].size(),
        elem2dToElems3d.m_toElementIndex[e2d].size()
      };

      if( sizes.size() != 1 || sizes.find( 2 ) == sizes.cend() )
      {
        continue;
      }

      localIndex const f0 = elem2dToFaces[e2d][0];
      localIndex const er0 = elem2dToElems3d.m_toElementRegion[e2d][0];
      localIndex const esr0 = elem2dToElems3d.m_toElementSubRegion[e2d][0];
      localIndex const ei0 = elem2dToElems3d.m_toElementIndex[e2d][0];
      auto const & faces0 = elem2dToElems3d.getElementRegionManager()->getRegion( er0 ).getSubRegion< CellElementSubRegion >( esr0 ).faceList()[ei0];

      localIndex const f1 = elem2dToFaces[e2d][1];
      localIndex const er1 = elem2dToElems3d.m_toElementRegion[e2d][1];
      localIndex const esr1 = elem2dToElems3d.m_toElementSubRegion[e2d][1];
      localIndex const ei1 = elem2dToElems3d.m_toElementIndex[e2d][1];
      auto const & faces1 = elem2dToElems3d.getElementRegionManager()->getRegion( er1 ).getSubRegion< CellElementSubRegion >( esr1 ).faceList()[ei1];

      bool const match00 = std::find( faces0.begin(), faces0.end(), f0 ) != faces0.end();
      bool const match11 = std::find( faces1.begin(), faces1.end(), f1 ) != faces1.end();
      bool const match01 = std::find( faces0.begin(), faces0.end(), f1 ) != faces0.end();
      bool const match10 = std::find( faces1.begin(), faces1.end(), f0 ) != faces1.end();

      bool const matchCrossed = !match00 && !match11 && match01 && match10;
      bool const matchStraight = match00 && match11 && !match01 && !match10;

      if( matchCrossed )
      {
        std::swap( elem2dToFaces[e2d][0], elem2dToFaces[e2d][1] );
      }
      else if( !matchStraight )
      {
        GEOS_ERROR( "Mapping neighbor inconsistency detected for fracture " << fractureName );
      }
    }
  }
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

  fixNeighborMappingsInconsistency( getName(), m_2dElemToElems, m_toFacesRelation );
}

/**
 * @brief Returns a mapping that links any collocated node to the collocated node with the lowest index.
 * @param elem2dToCollocatedNodesBuckets All the collocated nodes buckets.
 * @return The computed mapping.
 * @details For each collocated node, the returned mapping will provide
 * the lowest id among all the collocated nodes sharing the same position.
 * That way, it's possible to know if two nodes are collocated of each other by checking if they share the same lowest id.
 */
std::map< globalIndex, globalIndex > buildReferenceCollocatedNodes( ArrayOfArrays< array1d< globalIndex > > const & elem2dToCollocatedNodesBuckets )
{
  std::map< globalIndex, globalIndex > referenceCollocatedNodes;  // Will be returned.

  // Since some 2d elem may share some nodes, some of the collocated nodes buckets will be duplicated.
  // We want to remove this.
  std::set< std::set< globalIndex > > uniqueCollocatedNodesBuckets;
  for( localIndex e2d = 0; e2d < elem2dToCollocatedNodesBuckets.size(); ++e2d )
  {
    for( integer ni = 0; ni < elem2dToCollocatedNodesBuckets[e2d].size(); ++ni )
    {
      array1d< globalIndex > const & tmp = elem2dToCollocatedNodesBuckets( e2d, ni );
      std::set< globalIndex > const tmp2( tmp.begin(), tmp.end() );
      uniqueCollocatedNodesBuckets.insert( tmp2 );
    }
  }

  // Now we can perform simple loops to compute the mapping.
  for( std::set< globalIndex > const & bucket: uniqueCollocatedNodesBuckets )
  {
    globalIndex const & refNode = *std::min_element( bucket.cbegin(), bucket.cend() );
    for( globalIndex const & n: bucket )
    {
      referenceCollocatedNodes[n] = refNode;
    }
  }

  return referenceCollocatedNodes;
}


/**
 * @brief Returns a mapping that links any collocated edge to the collocated edge with the lowest index.
 * @param referenceCollocatedNodes The mapping that links any collocated node to its collocated node with the lowest index.
 * @param nl2g The local to glocal mapping for nodes.
 * @param edgeToNodes The edge to nodes mapping.
 * @param elem2dToEdges The 2d elem to edges mapping.
 * @param edgeGhostRanks The ghost rank of the edges.
 * @return The computed map.
 */
std::map< localIndex, localIndex > buildReferenceCollocatedEdges( std::map< globalIndex, globalIndex > const & referenceCollocatedNodes,
                                                                  arrayView1d< globalIndex const > const nl2g,
                                                                  EdgeManager::NodeMapType const & edgeToNodes,
                                                                  ArrayOfArraysView< localIndex const > const elem2dToEdges,
                                                                  arrayView1d< integer const > edgeGhostRanks )
{
  // `edgeIds` maps the nodes of the edges of the face element sub-region to its index.
  std::map< std::pair< globalIndex, globalIndex >, localIndex > edgesIds;
  for( int ei = 0; ei < elem2dToEdges.size(); ++ei )
  {
    for( localIndex const & edi: elem2dToEdges[ei] )
    {
      auto const nodes = edgeToNodes[edi];
      GEOS_ASSERT_EQ( nodes.size(), 2 );
      auto const p = std::minmax( { nl2g[nodes[0]], nl2g[nodes[1]] } );
      edgesIds[p] = edi;
    }
  }

  // The key of the `collocatedEdgeBuckets` map (i.e. `std::pair< globalIndex, globalIndex >`) represents the global indices of two nodes.
  // Those two nodes are the lowest index of collocated nodes. As such, those two nodes may not form an existing edge.
  // But this trick lets us define some kind of _hash_ that allows to compare the location of the edges:
  // edges sharing the same hash lie in the same position.
  std::map< std::pair< globalIndex, globalIndex >, std::set< localIndex > > collocatedEdgeBuckets;
  // The `collocatedEdgeIds` map gathers all the collocated edges together.
  for( auto const & p: edgesIds )
  {
    std::pair< globalIndex, globalIndex > const & nodes = p.first;
    localIndex const & edge = p.second;

    auto it0 = referenceCollocatedNodes.find( nodes.first );
    globalIndex const n0 = it0 != referenceCollocatedNodes.cend() ? it0->second : nodes.first;

    auto it1 = referenceCollocatedNodes.find( nodes.second );
    globalIndex const n1 = it1 != referenceCollocatedNodes.cend() ? it1->second : nodes.second;

    std::pair< globalIndex, globalIndex > const edgeHash = std::minmax( n0, n1 );
    collocatedEdgeBuckets[edgeHash].insert( edge );
  }

  std::map< localIndex, localIndex > referenceCollocatedEdges;

  // We want to consider in priority the edges that are owned by the rank.
  // So this comparator takes the ghost rank into account and favors lower values of the ghost rank.
  auto comp = [=]( int e,
                   int f ) -> bool
  {
    int const re = edgeGhostRanks[e] < 0 ? 0 : 1;
    int const rf = edgeGhostRanks[f] < 0 ? 0 : 1;
    return std::tie( re, e ) < std::tie( rf, f );
  };

  for( auto const & p: collocatedEdgeBuckets )
  {
    std::set< localIndex > const & collocatedEdges = p.second;
    localIndex const refEdge = *std::min_element( collocatedEdges.cbegin(), collocatedEdges.cend(), comp );
    for( localIndex const & collocatedEdge: collocatedEdges )
    {
      referenceCollocatedEdges[collocatedEdge] = refEdge;
    }
  }

  return referenceCollocatedEdges;
}


/**
 * @brief Returns a SortedArray with sequentially increasing values, starting with value and repetitively evaluating ++value
 * @param newSize The size of the returned array.
 * @param value The starting value.
 * @return The filled array.
 */
SortedArray< localIndex > makeSortedArrayIota( localIndex newSize, localIndex value = 0 )
{
  SortedArray< localIndex > result;
  result.reserve( newSize );
  for( localIndex i = value; i < value + newSize; ++i )
  {
    result.insert( i );
  }
  return result;
}


/**
 * @brief Builds the 2d face to 2d elements mapping.
 * @param num2dFaces The number of 2d faces.
 * @param num2dElems The number of 2d elements.
 * @param elem2dToEdges The 2d element (geometrical faces in 3d) to edges mapping.
 * @param edgesTo2dFaces The edges to 2d faces (geometrical edges in 3d) mapping.
 * @param referenceCollocatedEdges The mapping that, for a given edge index, returns the collocated edge with lowest index.
 * @return The computed mapping.
 */
ArrayOfArrays< geos::localIndex > build2dFaceTo2dElems( std::size_t num2dFaces,
                                                        localIndex num2dElems,
                                                        ArrayOfArraysView< localIndex const > const elem2dToEdges,
                                                        map< localIndex, localIndex > const & edgesTo2dFaces,
                                                        std::map< geos::localIndex, geos::localIndex > const & referenceCollocatedEdges )
{
  ArrayOfArrays< localIndex > m_2dFaceTo2dElems;
  // `tmp` contains the 2d face to 2d elements mappings as a `std` container.
  // Eventually, it's copied into an `LvArray` container.
  std::vector< std::vector< localIndex > > tmp( num2dFaces );
  for( auto i = 0; i < num2dElems; ++i )
  {
    for( auto const & e: elem2dToEdges[i] )
    {
      tmp[edgesTo2dFaces.at( referenceCollocatedEdges.at( e ) )].push_back( i );
    }
  }
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

  return m_2dFaceTo2dElems;
}


void FaceElementSubRegion::fixSecondaryMappings( NodeManager const & nodeManager,
                                                 EdgeManager const & edgeManager,
                                                 FaceManager const & faceManager,
                                                 ElementRegionManager const & elemManager )
{
  arrayView1d< globalIndex const > const nl2g = nodeManager.localToGlobalMap();

  // First let's create the reference mappings for both nodes and edges.
  std::map< globalIndex, globalIndex > const referenceCollocatedNodes = buildReferenceCollocatedNodes( m_2dElemToCollocatedNodesBuckets );
  std::map< localIndex, localIndex > const referenceCollocatedEdges =
    buildReferenceCollocatedEdges( referenceCollocatedNodes, nl2g, edgeManager.nodeList(), m_toEdgesRelation.toViewConst(), edgeManager.ghostRank().toViewConst() );

  localIndex const num2dElems = this->size();
  localIndex const num2dFaces = LvArray::integerConversion< localIndex >( referenceCollocatedEdges.size() );

  // For the `m_2dFaceToEdge`, we can select any values that we want.
  // But then we need to be consistent...
  m_2dFaceToEdge.clear();
  m_2dFaceToEdge.reserve( num2dFaces );
  for( auto const & p: referenceCollocatedEdges )
  {
    globalIndex const & refEdge = p.second;
    m_2dFaceToEdge.emplace_back( refEdge );
  }

  // `m_edgesTo2dFaces` is computed by the simple inversion of `m_2dFaceToEdge`
  m_edgesTo2dFaces.clear();
  for( localIndex i = 0; i < num2dFaces; ++i )
  {
    m_edgesTo2dFaces[m_2dFaceToEdge[i]] = i;
  }

  // Mark 2d elements and their connections as "new" so they get considered during the stencil computations.
  // This is mainly due to the dynamic process of the fracture propagation and the legacy unique way to import the fracture
  // by splitting the mesh during the first steps of the simulations.
  m_newFaceElements = makeSortedArrayIota( num2dElems );
  m_recalculateConnectionsFor2dFaces = makeSortedArrayIota( num2dFaces );

  m_2dFaceTo2dElems = build2dFaceTo2dElems( num2dFaces, num2dElems, m_toEdgesRelation.toViewConst(), m_edgesTo2dFaces, referenceCollocatedEdges );

  // When a fracture element has only one (or less!) neighbor, let's try to find the other one.
  // The `ElemPath` provides all the information of a given face: obviously its face index,
  // but also the element index that touch the face, and the node indices of the face as well.
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
  };

  // We are building the mapping that connects all the reference (collocated) nodes of any face to the elements those nodes are touching.
  // Using this nodal information will let use reconnect the fracture 2d element to its 3d neighbor.
  std::map< std::set< globalIndex >, std::set< ElemPath > > faceRefNodesToElems;
  elemManager.forElementSubRegionsComplete< CellElementSubRegion >( [&]( localIndex const er,
                                                                         localIndex const esr,
                                                                         ElementRegionBase const & GEOS_UNUSED_PARAM( region ),
                                                                         CellElementSubRegion const & subRegion )
  {
    auto const & elemToFaces = subRegion.faceList().base();
    auto const & faceToNodes = faceManager.nodeList();
    for( localIndex ei = 0; ei < elemToFaces.size( 0 ); ++ei )
    {
      for( auto const & face: elemToFaces[ei] )
      {
        // A set of the global indices of the nodes of the face is used as the "signature" of the face nodes.
        std::set< globalIndex > nodesOfFace;
        for( localIndex const & n: faceToNodes[face] )
        {
          auto const it = referenceCollocatedNodes.find( nl2g[n] );
          if( it != referenceCollocatedNodes.cend() )
          {
            nodesOfFace.insert( it->second );
          }
          else
          {
            // If all the nodes are not found, it makes no sense to keep on considering the current face as a candidate.
            // So we can exit the loop.
            break;
          }
        }
        // We still double-check that all the nodes of the face were found,
        // even if there's no chance for this entry to successfully match any request.
        if( nodesOfFace.size() == LvArray::integerConversion< std::size_t >( faceToNodes[face].size() ) )
        {
          std::vector< localIndex > const nodes( faceToNodes[face].begin(), faceToNodes[face].end() );
          faceRefNodesToElems[nodesOfFace].insert( ElemPath{ er, esr, ei, face, nodes } );
        }
      }
    }
  } );

  // The `misMatches` array will contain fracture element that did not find all of its neighbors.
  // This is used to display a more precise error message.
  std::vector< localIndex > misMatches;
  // Here we loop over all the elements of the fracture.
  // When there's neighbor missing, we search for a face that would lie on the collocated nodes of the fracture element.
  for( int e2d = 0; e2d < num2dElems; ++e2d )
  {
    if( m_2dElemToElems.m_toElementIndex.sizeOfArray( e2d ) >= 2 )  // All the neighbors are known.
    {
      continue;
    }

    std::set< globalIndex > refNodes;
    if( m_toNodesRelation[e2d].size() != 0 )
    {
      for( localIndex const & n: m_toNodesRelation[e2d] )
      {
        globalIndex const & gn = nl2g[n];
        auto const it = referenceCollocatedNodes.find( gn );
        if( it != referenceCollocatedNodes.cend() )
        {
          refNodes.insert( it->second );
        }
      }
    }
    else if( m_ghostRank[e2d] < 0 )
    {
      for( integer ni = 0; ni < m_2dElemToCollocatedNodesBuckets[e2d].size(); ++ni )
      {
        array1d< globalIndex > const & bucket = m_2dElemToCollocatedNodesBuckets( e2d, ni );
        for( globalIndex const & gni: bucket )
        {
          auto const it = referenceCollocatedNodes.find( gni );
          if( it != referenceCollocatedNodes.cend() )
          {
            refNodes.insert( it->second );
          }
        }
      }
    }

    auto const match = faceRefNodesToElems.find( refNodes );
    if( match != faceRefNodesToElems.cend() )
    {
      for( ElemPath const & path: match->second )
      {
        // This `if` prevents from storing the same data twice.
        if( m_2dElemToElems.m_toElementIndex.sizeOfArray( e2d ) == 0 || m_2dElemToElems.m_toElementIndex[e2d][0] != path.ei )
        {
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
    }
  }

  GEOS_ERROR_IF( !misMatches.empty(),
                 "Fracture " << this->getName() << " has elements {" << stringutilities::join( misMatches, ", " ) << "} without two neighbors." );

  // Checking that each face has two neighboring elements.
  // If not, we pop up an error.
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

  fixNeighborMappingsInconsistency( getName(), m_2dElemToElems, m_toFacesRelation );
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

std::set< std::set< globalIndex > > FaceElementSubRegion::getCollocatedNodes() const
{
  std::set< std::set< globalIndex > > result;

  for( localIndex e2d = 0; e2d < m_2dElemToCollocatedNodesBuckets.size(); ++e2d )
  {
    for( integer ni = 0; ni < m_2dElemToCollocatedNodesBuckets[e2d].size(); ++ni )
    {
      array1d< globalIndex > const & bucket = m_2dElemToCollocatedNodesBuckets( e2d, ni );
      std::set< globalIndex > const s( bucket.begin(), bucket.end() );
      result.insert( s );
    }
  }

  return result;
}

} /* namespace geos */
