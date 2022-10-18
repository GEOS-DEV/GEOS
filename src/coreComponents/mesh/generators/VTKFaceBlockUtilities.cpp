/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2020-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#include "VTKFaceBlockUtilities.hpp"

#include "dataRepository/Group.hpp"

#include <vtkCell.h>
#include <vtkDataArray.h>
#include <vtkFieldData.h>

#include <algorithm>

namespace geosx
{

/**
 * @brief Given a pair @p p of nodes, return the GEOSX edge index made of these 2 nodes.
 * @param p The pair of nodes in no specific order.
 * @param nodeToEdges The nodes to edges mapping.
 * @return The GEOSX edge index.
 */
int pairToEdge( std::pair< int, int > const & p,
                ArrayOfArraysView< localIndex const > nodeToEdges )
{
  // This functions takes all the edges that have one of the two nodes,
  // and find the unique edge with both nodes.
  // The current implementation can surely be improved,
  // but we want to challenge the feature first.
  localIndex const n0 = p.first, n1 = p.second;
  auto const & edges0 = nodeToEdges[n0];
  auto const & edges1 = nodeToEdges[n1];
  std::set< localIndex > const es0( edges0.begin(), edges0.end() );
  std::set< localIndex > const es1( edges1.begin(), edges1.end() );
  std::vector< localIndex > intersect;
  std::set_intersection( es0.cbegin(), es0.cend(), es1.cbegin(), es1.cend(), std::back_inserter( intersect ) );
  GEOSX_ERROR_IF_NE_MSG( intersect.size(), 1, "Internal error: we should only have one element" );
  return intersect.front();
}

/**
 * @brief Converts one fault/fracture raw information of the vtk field data input into a C++ class.
 */
class FaceData
{
public:
  /**
   * @brief Builds from the vtk input.
   * @param fieldData The vtk field data.
   * @param faceBlockName The face block name to search form in the @p fieldData.
   */
  FaceData( vtkFieldData * fieldData,
            string const & faceBlockName )
  {
    // TODO use a builder and not the constructor.
    // First, formatting the duplicated nodes' information.
    {
      vtkDataArray * duplicatedNodes = fieldData->GetArray( "duplicated_points_info" ); // TODO create constant.
      GEOSX_ERROR_IF( duplicatedNodes == nullptr, "Mandatory VTK field data array \"duplicated_points_info\" not found." );

      // vtk will provide the field as a (numComponents, numTuples) table.
      // numTuples is actually the number of duplicated nodes too.
      // But numComponents is just the maximum number of nodes multiplied nodes at one location.
      // For a T shaped fracture, the nodes at the T intersection will be tripled.
      // But all the other nodes will only get doubled.
      // Nevertheless, numComponents will equal 3 for all the nodes.
      // The dummy slot will be filled with negative dummy values.
      int const numComponents = duplicatedNodes->GetNumberOfComponents();
      int const numDuplicatedLocations = duplicatedNodes->GetNumberOfTuples();
      for( int i = 0; i < numDuplicatedLocations; ++i )
      {
        std::vector< double > rawTuple( numComponents );
        duplicatedNodes->GetTuple( i, rawTuple.data() );
        std::set< int > nodes;
        for( double const & v: rawTuple )
        {
          // This comparison filters our dummy values.
          if( v >= 0 )
          { nodes.insert( int( v ) ); }
        }
        // At that point, `nodes` contain all the duplicated nodes, sorted.
        for( const int & n: nodes )
        {
          auto const & it = m_duplicatedNodes.find( n );
          if( it == m_duplicatedNodes.end() )
          {
            m_duplicatedNodes[n] = *nodes.begin();
          }
          else
          {
            // For each duplicated node, we want to find the node with the lowest index.
            // This will help us find collocated edges, faces...
            it->second = std::min( it->second, *nodes.begin() );
          }
        }
      }
    }

    // Now formatting the faces' information.
    {
      vtkDataArray * array = fieldData->GetArray( faceBlockName.c_str() );
      GEOSX_ERROR_IF( array == 0, "Face block name \"" << faceBlockName << "\" was not found in the mesh." );

      m_num2dElements = array->GetNumberOfTuples();

      auto const numComponents = array->GetNumberOfComponents();
      GEOSX_ERROR_IF_NE_MSG( 4, numComponents, "FieldData \"" + faceBlockName + "\" should have 4 components. " << numComponents << " were found." );

      for( int i = 0; i < m_num2dElements; ++i )
      {
        std::vector< double > tuple( numComponents );
        array->GetTuple( i, tuple.data() );
        std::vector< globalIndex > res;
        res.reserve( numComponents );
        std::copy( tuple.cbegin(), tuple.cend(), std::back_inserter( res ) );
        m_data.push_back( res );
      }
    }
  }

  /**
   * @brief Number of 2d elements (geometrical faces in 3d).
   * @return The integer.
   */
  vtkIdType num2dElements() const
  {
    return m_num2dElements;
  }

  /**
   * @brief Given the vtk global index @p nodeId, returns the duplicated node with the lowest index, or @p nodeId itself if it's not duplicated.
   * @param nodeId The vtk global index.
   * @return The integer.
   *
   * This functions can be used to identify if two nodes are duplicated one of each other: they share the same lowest index.
   */
  int uniqueNode( int nodeId ) const
  {
    auto const it = m_duplicatedNodes.find( nodeId );
    return it == m_duplicatedNodes.cend() ? nodeId : it->second;
  }

  /**
   * @brief For @p elem2d, returns the global vtk index of @p neighbor element.
   * @param elem2d The local (to the "fracture") index of the 2d element.
   * @param neighbor The neighbor index, can be 0 or 1.
   * @return The global vtk index of the pointed cell.
   */
  int cell( int elem2d,
            int neighbor ) const
  {
    GEOSX_ERROR_IF( neighbor != 0 and neighbor != 1, "Wrong neighbor index." );
    return m_data[elem2d][2 * neighbor];
  }

  /**
   * @brief For @p elem2d, returns the local (to the vek cell) vtk index of @p neighbor element.
   * @param elem2d The local (to the "fracture") index of the 2d element.
   * @param neighbor The neighbor index, can be 0 or 1.
   * @return The local face vtk index of the pointed cell that touches the 2d element.
   */
  globalIndex face( int elem2d,
                    int neighbor ) const
  {
    GEOSX_ERROR_IF( neighbor != 0 and neighbor != 1, "Wrong neighbor index." );
    return m_data[elem2d][2 * neighbor + 1];
  }

private:
  /**
   * @brief Number of 2d elements (geometrical faces in 3d).
   */
  vtkIdType m_num2dElements;
  /**
   * @brief Raw vtk table that will be interpreted by @p FaceData.
   * Each row of the table contains the global index of an element,
   * then the element-local index of the face that makes belongs to the fracture.
   * And the same again for the second element.
   */
  std::vector< std::vector< globalIndex > > m_data;
  /**
   * @brief All for any duplicated node, among all duplicated nodes in the same "position", returns the node index with the lowest value.
   */
  std::map< int, int > m_duplicatedNodes;
};

namespace internal
{

/**
 * @brief This classes helps at building the global (vtk) element indices to the indices local to the cell blocks.
 * Then the global element to global faces can be be also computed.
 *
 * This implementation is mostly brittle and is work in progress.
 */
class ElementToFace
{
public:
  /**
   * @brief Builds for the cell blocks.
   * @param[in] cellBlocks The cell blocks group.
   */
  ElementToFace( dataRepository::Group const & cellBlocks )
    : m_localToGlobalMaps( cellBlocks.numSubGroups() ),
      m_elemToFacesMaps( cellBlocks.numSubGroups() )
  {
    // This implementation here is a bit brittle and dummy.
    // It's based on the fact that each cell block often spans a continuous range of indices.
    // This help us not to store a big map.
    // We check the ranges of the cell blocks and kill the process if it's not the case.
    // Something better should be done when more robustness is searched on the subject.
    localIndex const numCellBlocks = cellBlocks.numSubGroups();
    for( int i = 0; i < numCellBlocks; ++i )
    {
      CellBlock const * cb = cellBlocks.getGroupPointer< CellBlock >( i );
      m_localToGlobalMaps[i] = cb->localToGlobalMapConstView();
      m_elemToFacesMaps[i] = cb->getElemToFacesConstView();
    }

    for( int i = 0; i < numCellBlocks; ++i )
    {
      arrayView1d< globalIndex const > const l2g = m_localToGlobalMaps[i];
      auto const size = l2g.size();
      for( int j = 0; j < size - 1; ++j )
      {
        GEOSX_ERROR_IF_NE_MSG( l2g[i] + 1, l2g[i + 1], "Internal error: the code is assuming that the cell block indexing is continuously increasing." );
      }
      m_ranges.emplace_back( l2g.front(), l2g.back() );
    }
  }

  /**
   * @brief Given the global (vtk) id, returns the cell block index @p ei belongs to.
   * @param[in] ei global cell index.
   * @return The cell block index.
   */
  int getCellBlockIndex( int const & ei ) const
  {
    auto const it = std::find_if( m_ranges.cbegin(), m_ranges.cend(), [&]( auto const & p ) -> bool
    {
      return p.first <= ei and ei <= p.second;
    } );

    return std::distance( m_ranges.cbegin(), it );
  }

  /**
   * @brief From global index @p ei, returns the local element index in its cell block.
   * @param ei The global vtk element index.
   * @return The integer.
   */
  localIndex getElementIndexInCellBlock( int const & ei ) const
  {
    int const cbi = getCellBlockIndex( ei );
    int const offset = getElementCellBlockOffset( cbi );
    return ei - offset;
  }

  /**
   * @brief Given the global (vtk) id, returns the index local to the cell block @p ei belong to.
   * @param[in] ei global cell index.
   * @return The faces for element @p ei.
   */
  auto operator[]( int const & ei ) const
  {
    const localIndex iCellBlock = getCellBlockIndex( ei );
    arrayView2d< localIndex const > const e2f = m_elemToFacesMaps[iCellBlock];
    return e2f[ei - getElementCellBlockOffset( iCellBlock )];
  }

private:

  /**
   * @brief Returns the index offset for cell block of index @p iCellBlock
   * @param[in] iCellBlock The index of the cell block.
   * @return The offset.
   * Cell block indices span a continuous range of global indices.
   * The offset makes possible to compute the cell block local index using a simple subtraction.
   */
  int getElementCellBlockOffset( int const & iCellBlock ) const
  {
    return m_ranges[iCellBlock].first;
  }

  /**
   * @brief The cell indices ranges for each cell block.
   * he range is inclusive since min and max are in the cell block.
   */
  std::vector< std::pair< int, int > > m_ranges;
  /**
   * @brief The local to global map, for each cell block.
   * The indexing the the same as the cell block indexing.
   */
  std::vector< arrayView1d< globalIndex const > > m_localToGlobalMaps;
  /**
   * @brief The element to faces map, for each cell block.
   * The indexing the the same as the cell block indexing.
   */
  std::vector< arrayView2d< localIndex const > > m_elemToFacesMaps;
};

} // end of namespace internal

/**
 * @brief Computes the number of 2d faces (geometrical edges in 3d) for the fracture defined by @p faceData.
 * @param[in] vtkMesh The vtk mesh.
 * @param[in] faceData The parsed fault/fracture information.
 * @return The number of 2d faces.
 */
localIndex computeNum2dFaces( vtkSmartPointer< vtkDataSet > vtkMesh,
                              FaceData const & faceData )
{
  vtkIdType const num2dElements = faceData.num2dElements();
  // In order to compute the number of edges, we store all the edges nodes, and then we remove duplicates.
  // To deal with duplicated nodes, we do not store the given edge nodes, but the lowest duplicated node index instead.
  std::vector< std::pair< int, int > > allPairs;
  for( int i = 0; i < num2dElements; ++i )
  {
    for( int j = 0; j < 2; ++j )
    {
      vtkCell * cell = vtkMesh->GetCell( faceData.cell( i, j ) );
      vtkCell * face = cell->GetFace( faceData.face( i, j ) );
      for( int k = 0; k < face->GetNumberOfEdges(); ++k )
      {
        vtkCell * edge = face->GetEdge( k );
        vtkIdType const pt0 = edge->GetPointId( 0 ), pt1 = edge->GetPointId( 1 );
        allPairs.push_back( std::minmax( faceData.uniqueNode( pt0 ), faceData.uniqueNode( pt1 ) ) );
      }
    }
  }
  std::set< std::pair< int, int > > uniquePairs( allPairs.cbegin(), allPairs.cend() );
  return uniquePairs.size();
}

/**
 * @brief Compute the mapping from the 2d elements (geometrical faces in 3d) to their 3d element neighbors.
 * @param[in] faceData The parsed fault/fracture information.
 * @param[in] elemToFaces The element to face mapping, across all the cell blocks.
 * @return The mappings to both the 3d elements and the cell block they belong to.
 */
ToCellRelation< array2d< localIndex > > compute2dElemToElems( FaceData const & faceData,
                                                              internal::ElementToFace const & elemToFaces )
{
  vtkIdType const num2dElements = faceData.num2dElements();

  array2d< localIndex > elem2dToCellBlock( num2dElements, 2 );
  elem2dToCellBlock.setValues< parallelHostPolicy >( -1 );
  for( int i = 0; i < num2dElements; ++i )
  {
    for( int j = 0; j < 2; ++j )
    {
      elem2dToCellBlock[i][j] = elemToFaces.getCellBlockIndex( faceData.cell( i, j ) );
    }
  }

  array2d< localIndex > elem2dToElems( num2dElements, 2 );
  elem2dToElems.setValues< parallelHostPolicy >( -1 );
  for( int i = 0; i < num2dElements; ++i )
  {
    for( int j = 0; j < 2; ++j )
    {
      // The offset lets us convert the global index to an index local to the cell block.
      // It surely implements a pattern, which is indisputably bad and will need to be updated.
      elem2dToElems[i][j] = elemToFaces.getElementIndexInCellBlock( faceData.cell( i, j ) );
    }
  }

  return ToCellRelation< array2d< localIndex > >( elem2dToCellBlock, elem2dToElems );
}

/**
 * @brief Compute the 2d elements to nodes and faces
 * @param[in] vtkMesh The vtk mesh.
 * @param[in] faceData The parsed fault/fracture information.
 * @param[in] elemToFaces The element to face mapping, across all the cell blocks.
 * @param[in] faceToNodes The face to nodes mappings.
 * @return A tuple. First index being the node mapping, second being the face mapping.
 */
std::tuple< ArrayOfArrays< localIndex >, array2d< localIndex > >
compute2dElemToNodesAndFaces( vtkSmartPointer< vtkDataSet > vtkMesh,
                              FaceData const & faceData,
                              internal::ElementToFace const & elemToFaces,
                              ArrayOfArraysView< localIndex const > faceToNodes )
{
  // Small utility function to convert a vtkIdList into a std::set< int >.
  auto vtkIdListToSet = []( vtkIdList * in ) -> std::set< int >
  {
    std::set< int > result;
    for( int i = 0; i < in->GetNumberOfIds(); ++i )
    { result.insert( in->GetId( i ) ); }
    return result;
  };

  vtkIdType const num2dElements = faceData.num2dElements();

  ArrayOfArrays< localIndex > elem2dToNodes( num2dElements );
  array2d< localIndex > elem2dToFaces( num2dElements, 2 );

  // Standard element/face double loop to find/extract the nodes and faces indices required for the mappings.
  for( int i = 0; i < num2dElements; ++i )
  {
    std::vector< localIndex > nodes;
    std::vector< localIndex > faces;
    for( int j = 0; j < 2; ++j )
    {
      int const ei = faceData.cell( i, j );
      vtkCell * c = vtkMesh->GetCell( ei );
      vtkCell * f = c->GetFace( faceData.face( i, j ) );
      std::set< int > const vtkFaceNodes = vtkIdListToSet( f->GetPointIds() );
      // This loop somehow find which geosx face matches the given vtk local face index in `faceData`.
      for( int const & face: elemToFaces[ei] )
      {
        auto const tmp = faceToNodes[face];
        std::set< int > const faceNodes( tmp.begin(), tmp.end() );
        if( vtkFaceNodes == faceNodes )
        {
          std::copy( faceNodes.cbegin(), faceNodes.cend(), back_inserter( nodes ) );
          faces.push_back( face );
          break;
        }
      }
    }
    elem2dToNodes.resizeArray( i, nodes.size() );
    std::copy( nodes.cbegin(), nodes.cend(), elem2dToNodes[i].begin() );
    GEOSX_ERROR_IF_NE_MSG( faces.size(), 2, "Internal error: wrong number of faces in the fracture/fault mesh computation" );
    elem2dToFaces[i][0] = faces[0];
    elem2dToFaces[i][1] = faces[1];
  }

  return { elem2dToNodes, elem2dToFaces };
}

/**
 * @brief Compute the 2d face and 2d elem to edges mapping.
 * @param[in] vtkMesh The vtk mesh.
 * @param[in] faceData The parsed fault/fracture information.
 * @param[in] num2dFaces The number of 2d faces (geometrical edges in 3d).
 * @param[in] nodeToEdges The nodes to edges mapping.
 * @return A tuple. First index being the face mapping, second being the elem mapping.
 *
 * Actually, for both the 2d face to edges mapping,
 * it's important to note that there are actually 2 edges that match any 2d face.
 * But the returned mapping only has one edge value.
 * For the moment, this value is only used to define some kind of equivalent identification,
 * and therefore one value is enough.
 * From the two edges values, we only return the smallest value, to be sure it's unique.
 */
std::tuple< array1d< localIndex >, ArrayOfArrays< localIndex > >
compute2dFaceAnd2dElemToEdges( vtkSmartPointer< vtkDataSet > vtkMesh,
                               FaceData const & faceData,
                               localIndex const num2dFaces,
                               ArrayOfArraysView< localIndex const > nodeToEdges )
{
  // What we hereafter call a `pairSig` (for `pair signature`) the pair
  // with the lowest duplicated node indices for each of the two indices.
  // It's worth mentioning that such a pair may not be an existing edge.
  // It's a more kind of signature/hash to reconcile duplicated edges.
  // This function computes the pair id for any given pair of nodes (edge).
  auto const pairToPairSig = [&faceData]( std::pair< int, int > const & p ) -> std::pair< int, int >
  {
    return std::minmax( faceData.uniqueNode( p.first ), faceData.uniqueNode( p.second ) );
  };

  vtkIdType const num2dElements = faceData.num2dElements();

  // The pairs of `allEdges` do not deal with the duplication challenge: all the edges are there, even duplicates.
  // [2d elem index] -> [edges (stored as a pair of ints)]
  std::vector< std::vector< std::pair< int, int > > > allEdges;
  for( int i = 0; i < num2dElements; ++i )
  {
    allEdges.push_back( {} );
    for( int j = 0; j < 2; ++j )
    {
      vtkCell * cell = vtkMesh->GetCell( faceData.cell( i, j ) );
      vtkCell * face = cell->GetFace( faceData.face( i, j ) );
      for( int k = 0; k < face->GetNumberOfEdges(); ++k )
      {
        vtkCell * edge = face->GetEdge( k );
        vtkIdType const p0 = edge->GetPointId( 0 ), p1 = edge->GetPointId( 1 );
        allEdges[i].push_back( std::minmax( p0, p1 ) ); // FIXE is it useful to sort here?
      }
    }
  }

  // We now build the `pairSig` -> `edge` mapping.
  // `edgeId` is the minimal (existing) edge index for the given `pairSig` signature.
  // Two edges sharing the same `edgeId` are therefore duplicated.
  std::map< std::pair< int, int >, int > pairSigToEdge;
  for( auto const & edges: allEdges )
  {
    for( std::pair< int, int > const & p: edges )
    {
      std::pair< int, int > const pairSig = pairToPairSig( p );

      auto const & it = pairSigToEdge.find( pairSig );
      if( it == pairSigToEdge.end() )
      {
        pairSigToEdge[pairSig] = pairToEdge( p, nodeToEdges );
      }
      else
      {
        it->second = std::min( it->second, pairToEdge( p, nodeToEdges ) );
      }
    }
  }
  GEOSX_ERROR_IF_NE_MSG( LvArray::integerConversion< std::size_t >( num2dFaces ), pairSigToEdge.size(), "Internal error: inconsistency when computing the fault and fracture mappings." );

  // This loop build the `elem2dToEdges` mapping.
  // It uses the `pairSigToEdgeId` to remove duplicates.
  ArrayOfArrays< localIndex > elem2dToEdges( num2dElements );
  for( int i = 0; i < num2dElements; ++i )
  {
    std::set< int > edges;
    for( auto const & p: allEdges[i] )
    {
      edges.insert( pairSigToEdge.at( pairToPairSig( p ) ) );
    }
    elem2dToEdges.resizeArray( i, edges.size() );
    std::copy( edges.cbegin(), edges.cend(), elem2dToEdges[i].begin() );
  }

  // Here we compute the `face2dToEdges` mapping.
  // What matters is that each `face2d` matches the proper `edge`.
  // But we do not pay attention to the ordering of the `face2d`, it will be in the sorting order.
  auto const cmp = [&pairToPairSig]( std::pair< int, int > const & p0,
                                    std::pair< int, int > const & p1 ) -> bool
  {
    return pairToPairSig( p0 ) < pairToPairSig( p1 );
  };
  std::set< std::pair< int, int >, decltype( cmp ) > uniqueEdges( cmp );
  for( std::vector< std::pair< int, int > > const & edges: allEdges )
  {
    uniqueEdges.insert( edges.cbegin(), edges.cend() );
  }
  GEOSX_ERROR_IF_NE_MSG( LvArray::integerConversion< localIndex >( uniqueEdges.size() ), num2dFaces, "Internal error: wrong number of edges in the fracture/fault mesh computation." );

  array1d< localIndex > face2dToEdges;
  face2dToEdges.reserve( num2dFaces );
  for( std::pair< int, int > const & p: uniqueEdges )
  {
    face2dToEdges.emplace_back( pairToEdge( p, nodeToEdges ) );
  }

  return { face2dToEdges, elem2dToEdges };
}

/**
 * @brief Compute the 2d faces (geometrical edges in 2d) to 2d elems (geometrical faces in 3d) mapping.
 * @param[in] num2dFaces The number of 2d faces (geometrical edges in 3d).
 * @param[in] face2dToEdges The 2d faces (geometrical edges in 2d) to edges mapping.
 * @param[in] elem2dToEdges The 2d elem (geometrical faces in 3d) to edges mapping.
 * @return The mapping.
 */
ArrayOfArrays< localIndex > compute2dFaceToElems2d( localIndex const num2dFaces,
                                                    arrayView1d< localIndex const > face2dToEdges,
                                                    ArrayOfArraysView< localIndex const > elem2dToEdges )
{
  ArrayOfArrays< localIndex > face2dToElems2d( num2dFaces );

  // First, a mapping inversion
  std::map< int, std::vector< int > > edgesToElems2d;
  for( localIndex i = 0; i < elem2dToEdges.size(); ++i )
  {
    for( int const & e: elem2dToEdges[i] )
    {
      edgesToElems2d[e].push_back( i );
    }
  }
  // Then we fill the 2d face -> 2d elems mappings using the 2d face -> edges and edges -> 2d elems chain.
  for( int fi = 0; fi < num2dFaces; ++fi )
  {
    localIndex const & edge = face2dToEdges[fi];
    std::vector< int > const & elems2d = edgesToElems2d[edge];
    face2dToElems2d.resizeArray( fi, elems2d.size() );
    std::copy( elems2d.begin(), elems2d.end(), face2dToElems2d[fi].begin() );
  }

  return face2dToElems2d;
}

void importFracture( string const & faceBlockName,
                     vtkSmartPointer< vtkDataSet > vtkMesh,
                     CellBlockManager & cellBlockManager )
{
  vtkFieldData * fieldData = vtkMesh->GetFieldData();
  GEOSX_ERROR_IF( fieldData == nullptr, "Fault and fracture information is supposed to be defined in the VTK field data. None was found." );

  FaceData const faceData( fieldData, faceBlockName );

  vtkIdType const num2dElements = faceData.num2dElements();
  localIndex const num2dFaces = computeNum2dFaces( vtkMesh, faceData );

  internal::ElementToFace const elemToFaces( cellBlockManager.getCellBlocks() );

  ToCellRelation< array2d< localIndex > > elem2dToElems = compute2dElemToElems( faceData, elemToFaces );

  ArrayOfArrays< localIndex > const faceToNodes = cellBlockManager.getFaceToNodes();
  ArrayOfArrays< localIndex > const nodeToEdges = cellBlockManager.getNodeToEdges();

  auto tmp = compute2dElemToNodesAndFaces( vtkMesh, faceData, elemToFaces, faceToNodes.toViewConst() );
  ArrayOfArrays< localIndex > const elem2dToNodes( std::get< 0 >( tmp ) );
  array2d< localIndex > const elem2dToFaces( std::get< 1 >( tmp ) );

  auto tmp2 = compute2dFaceAnd2dElemToEdges( vtkMesh, faceData, num2dFaces, nodeToEdges.toViewConst() );
  array1d< localIndex > const face2dToEdges( std::get< 0 >( tmp2 ) );
  ArrayOfArrays< localIndex > const elem2dToEdges( std::get< 1 >( tmp2 ) );

  ArrayOfArrays< localIndex > const face2dToElems2d = compute2dFaceToElems2d( num2dFaces, face2dToEdges.toViewConst(), elem2dToEdges.toViewConst() );

  // Mappings are now computed. Just create the face block by value.
  FaceBlock & faceBlock = cellBlockManager.registerFaceBlock( faceBlockName );
  faceBlock.setNum2DElements( num2dElements );
  faceBlock.setNum2DFaces( num2dFaces );
  faceBlock.set2dElemToNodes( elem2dToNodes );
  faceBlock.set2dElemToEdges( elem2dToEdges );
  faceBlock.set2dElemToFaces( elem2dToFaces );
  faceBlock.set2dFaceTo2dElems( face2dToElems2d );
  faceBlock.set2dFaceToEdge( face2dToEdges );
  faceBlock.set2dElemToElems( elem2dToElems );
}

} // end of namespace
