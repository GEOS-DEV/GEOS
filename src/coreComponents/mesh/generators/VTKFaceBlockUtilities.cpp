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
#include <vtkSmartPointer.h>

#include <algorithm>

namespace geosx
{

/**
 * @brief Given a pair @p p of nodes, return the GEOSX edge index made of these 2 nodes.
 * @param p The pair of nodes in no specific order.
 * @param n2ed The nodes to edges mapping.
 * @return The GEOSX index
 */
int pairToEdge( std::pair< int, int > const & p,
                ArrayOfArrays< localIndex > const & n2ed )
{
  // This functions takes all the edges that have one of the two nodes,
  // and find the unique edge with both nodes.
  // The current implementation can surely be improved,
  // but for the moment it seems to be working.
  int const n0 = p.first, n1 = p.second;
  auto const & edges0 = n2ed[n0];
  auto const & edges1 = n2ed[n1];
  std::set< int > const es0( edges0.begin(), edges0.end() );
  std::set< int > const es1( edges1.begin(), edges1.end() );
  std::vector< int > intersect;
  std::set_intersection( es0.cbegin(), es0.cend(), es1.cbegin(), es1.cend(), std::back_inserter( intersect ) );
  GEOSX_ERROR_IF_NE_MSG( intersect.size(), 1, "We should only have one element" );
  return intersect.front();
}

class FaceData
{
public:
  FaceData( vtkFieldData * fieldData,
            string const & faceBlockName )
  {
    // TODO use a builder and not the constructor.
    // Formatting the duplicated nodes' information.
    {
      vtkDataArray * duplicatedNodes = fieldData->GetArray( "duplicated_points_info" ); // TODO create constant.
      GEOSX_ERROR_IF( duplicatedNodes == nullptr, "Mandatory VTK field data array \"duplicated_points_info\" not found." );
      auto const numComponents = duplicatedNodes->GetNumberOfComponents();

      int const numDuplicatedLocations = duplicatedNodes->GetNumberOfTuples();
      for( int i = 0; i < numDuplicatedLocations; ++i )
      {
        std::vector< double > tuple( numComponents );
        duplicatedNodes->GetTuple( i, tuple.data() );
        std::set< int > nodes;
        for( double const & v: tuple )
        {
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
            it->second = std::min( it->second, *nodes.begin() );
          }
        }
      }
    }

    // Formatting the faces' information.
    {
      vtkDataArray * array = fieldData->GetArray( faceBlockName.c_str() );
      GEOSX_ERROR_IF( array == 0, "Face block name \"" << faceBlockName << "\" was not found in the mesh." );

      m_num2dElements = array->GetNumberOfTuples();

      auto const numComponents = array->GetNumberOfComponents();
      GEOSX_ERROR_IF_NE_MSG( 4, numComponents, "FieldData \"" + faceBlockName + "\" should have 4 components." );

      for( int i = 0; i < m_num2dElements; ++i )
      {
        std::vector< double > tuple( numComponents );
        array->GetTuple( i, tuple.data() );
        std::vector< int > res;
        res.reserve( numComponents );
        std::copy( tuple.cbegin(), tuple.cend(), std::back_inserter( res ) );
        m_data.push_back( res );
      }
    }
  }

  vtkIdType num2dElements() const
  {
    return m_num2dElements;
  }

  /**
   * @brief Given the vtk global index @p nodeId, returns the duplicated node with the lowest index.
   * @param nodeId
   * @return
   *
   * This functions can be used to identify if two nodes are duplicated one of each other: they share the same lowest index.
   */
  int uniqueNode( int nodeId ) const
  {
    auto const it = m_duplicatedNodes.find( nodeId );
    return it == m_duplicatedNodes.cend() ? nodeId : it->second;
  }

  /**
   * @brief For @p elem2d, returns the global vtk index of @p neighbor.
   * @param elem2d The local (the the "fracture") index of the 2d element.
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
   * @brief For @p elem2d, returns the local (to the vek cell) vtk index of @p neighbor.
   * @param elem2d The local (the the "fracture") index of the 2d element.
   * @param neighbor The neighbor index, can be 0 or 1.
   * @return The local face vtk index of the pointed cell that touches the 2d element.
   */
  int face( int elem2d,
            int neighbor ) const
  {
    GEOSX_ERROR_IF( neighbor != 0 and neighbor != 1, "Wrong neighbor index." );
    return m_data[elem2d][2 * neighbor + 1];
  }

private:
  vtkIdType m_num2dElements;
  std::vector< std::vector< int > > m_data;
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
 */
class ElementToFace
{
public:
  /**
   * @brief Builds for the cell blocks.
   * @param cellBlocks The cell blocks group.
   */
  ElementToFace( dataRepository::Group const & cellBlocks )
    : m_localToGlobalMaps( cellBlocks.numSubGroups() ),
      m_elemToFacesMaps( cellBlocks.numSubGroups() )
  {
    // Checking the ranges of the cell blocks.
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
        GEOSX_ERROR_IF_NE_MSG( l2g[i] + 1, l2g[i + 1], "The code is assuming that the cell block indexing is continuously increasing." );
      }
      m_ranges.emplace_back( l2g.front(), l2g.back() );
    }
  }

  /**
   * @brief Given the global (vtk) id, returns the cell block index @p ei belongs to.
   * @param ei global cell index.
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
   * @brief Returns the index offset for cell block of index @p iCellBlock
   * @param iCellBlock The index of the cell block.
   * @return The offset.
   * Cell block indices span a continuous range of global indices.
   * The offset makes possible to compute the cell block local index using a simple subtraction.
   */
  int getElementCellBlockOffset( int const & iCellBlock ) const
  {
    return m_ranges[iCellBlock].first;
  }

  /**
   * @brief Given the global (vtk) id, returns the index local to the cell block @p ei belong to.
   * @param ei global cell index.
   * @return The faces for element @p ei.
   */
  auto operator[]( int const & ei ) const
  {
    const std::size_t iCellBlock = getCellBlockIndex( ei );
    arrayView2d< localIndex const > const e2f = m_elemToFacesMaps[iCellBlock];
    return e2f[ei - m_ranges[iCellBlock].first];
  }

private:

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
 * @param vtkMesh The vtk mesh.
 * @param faceData The parsed fault/fracture information.
 * @return The number of 2d faces.
 */
localIndex computeNum2dFaces( vtkSmartPointer< vtkDataSet > vtkMesh,
                              FaceData const & faceData )
{
  vtkIdType const num2dElements = faceData.num2dElements();
  // Storing all the edges nodes, then we remove duplicates to get the number of different edges.
  // To deal with duplicated nodes, we do not store the edge nodes, but the lowest duplicated node index instead.
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
 * @param faceData The parsed fault/fracture information.
 * @param e2f The element to face mapping, across all the cell blocks.
 * @return The mappings to both the 3d elements and the cell block they belong to.
 */
ToCellRelation< array2d< localIndex > > compute2dElemToElems( FaceData const & faceData,
                                                              internal::ElementToFace const & e2f )
{
  vtkIdType const num2dElements = faceData.num2dElements();

  array2d< localIndex > elem2dToCellBlock( num2dElements, 2 ); // API
  elem2dToCellBlock.setValues< parallelHostPolicy >( -1 );
  for( int i = 0; i < num2dElements; ++i )
  {
    for( int j = 0; j < 2; ++j )
    {
      elem2dToCellBlock[i][j] = e2f.getCellBlockIndex( faceData.cell( i, j ) );
    }
  }

  array2d< localIndex > elem2dToElems( num2dElements, 2 ); // API
  elem2dToElems.setValues< parallelHostPolicy >( -1 );
  for( int i = 0; i < num2dElements; ++i )
  {
    for( int j = 0; j < 2; ++j )
    {
      auto const offset = e2f.getElementCellBlockOffset( elem2dToCellBlock[i][j] );
      elem2dToElems[i][j] = faceData.cell( i, j ) - offset;
    }
  }

  return ToCellRelation< array2d< localIndex > >( elem2dToCellBlock, elem2dToElems );
}

/**
 * @brief Compute the 2d elements to nodes and faces
 * @param vtkMesh The vtk mesh.
 * @param faceData The parsed fault/fracture information.
 * @param e2f The element to face mapping, across all the cell blocks.
 * @param f2n The face to nodes mappings.
 * @return A tuple. First index being the node mapping, second being the face mapping.
 */
std::tuple< ArrayOfArrays< localIndex >, array2d< localIndex > >
compute2dElemToNodesAndFaces( vtkSmartPointer< vtkDataSet > vtkMesh,
                              FaceData const & faceData,
                              internal::ElementToFace const & e2f,
                              ArrayOfArrays< localIndex > const & f2n )
{
  // A utility function to convert a vtkIdList into a std::set< int >.
  auto vtkIdListToSet = []( vtkIdList * in ) -> std::set< int >
  {
    std::set< int > result;
    for( int i = 0; i < in->GetNumberOfIds(); ++i )
    { result.insert( in->GetId( i ) ); }
    return result;
  };

  vtkIdType const num2dElements = faceData.num2dElements();

  ArrayOfArrays< localIndex > elem2dToNodes( num2dElements ); // API
  array2d< localIndex > elem2dToFaces( num2dElements, 2 ); // API

  // Standard element/face double loop to find extract the nodes and faces indices required for the mappings.
  for( int i = 0; i < num2dElements; ++i )
  {
    std::vector< localIndex > nodeLine;
    std::vector< localIndex > faceLine;
    for( int j = 0; j < 2; ++j )
    {
      int const ei = faceData.cell( i, j );
      vtkCell * c = vtkMesh->GetCell( faceData.cell( i, j ) );
      vtkCell * f = c->GetFace( faceData.face( i, j ) );
      std::set< int > const vtkPoints = vtkIdListToSet( f->GetPointIds() );
      auto faces = e2f[ei];
      for( int const & face: faces )
      {
        auto const faceNodes = f2n[face];
        std::set< int > const nodes( faceNodes.begin(), faceNodes.end() );
        if( vtkPoints == nodes )
        {
          std::copy( nodes.cbegin(), nodes.cend(), back_inserter( nodeLine ) );
          faceLine.push_back( face );
          break;
        }
      }
    }
    elem2dToNodes.resizeArray( i, nodeLine.size() );
    std::copy( nodeLine.cbegin(), nodeLine.cend(), elem2dToNodes[i].begin() );
    GEOSX_ERROR_IF_NE_MSG( faceLine.size(), 2, "Internal error: wrong number of faces in the fracture/fault mesh computation" );
    elem2dToFaces[i][0] = faceLine[0];
    elem2dToFaces[i][1] = faceLine[1];
  }

  return { elem2dToNodes, elem2dToFaces };
}

/**
 * @brief Compute the 2d face and 2d elem to edges mapping.
 * @param vtkMesh The vtk mesh.
 * @param faceData The parsed fault/fracture information.
 * @param num2dFaces The number of 2d faces (geometrical edges in 3d).
 * @param n2ed The nodes to edges mapping.
 * @return A tuple. First index being the face mapping, second being the elem mapping.
 */
std::tuple< array1d< localIndex >, ArrayOfArrays< localIndex >>
compute2dFaceAnd2dElemToEdges( vtkSmartPointer< vtkDataSet > vtkMesh,
                               FaceData const & faceData,
                               localIndex const num2dFaces,
                               ArrayOfArrays< localIndex > const & n2ed )
{
  vtkIdType const num2dElements = faceData.num2dElements();

  array1d< localIndex > face2dToEdges; // API
  face2dToEdges.reserve( num2dFaces );
  ArrayOfArrays< localIndex > elem2dToEdges( num2dElements ); // API

  // [2d elem index] -> [left/right index] -> [edges (stored as a pair of ints)]
  // The pairs of `allEdges` do not deal with the duplication challenge: all the edges are there, even duplicates.
  std::vector< std::vector< std::vector< std::pair< int, int > > > > allEdges;
  for( int i = 0; i < num2dElements; ++i )
  {
    allEdges.push_back( {} );
    for( int j = 0; j < 2; ++j )
    {
      vtkCell * cell = vtkMesh->GetCell( faceData.cell( i, j ) );
      vtkCell * face = cell->GetFace( faceData.face( i, j ) );
      allEdges[i].push_back( {} );
      for( int k = 0; k < face->GetNumberOfEdges(); ++k )
      {
        vtkCell * edge = face->GetEdge( k );
        vtkIdType const p0 = edge->GetPointId( 0 ), p1 = edge->GetPointId( 1 );
        allEdges[i][j].push_back( { std::min( p0, p1 ), std::max( p0, p1 ) } ); // FIXE is it useful to sort here?
      }
    }
  }

  // `allEdgesFlat` is the flat list containing all the edges, stored as pairs of nodes.
  // Duplicated nodes challenge is not considered here. Thus, there will be "duplicated" edges.
  std::vector< std::pair< int, int > > allEdgesFlat;
  for( auto const & leftRight: allEdges )
  {
    for( auto const & edges: leftRight )
    {
      for( auto const & edge: edges )
      {
        allEdgesFlat.push_back( edge );
      }
    }
  }

  // What we hereafter call a `pairId` is the pair with the lowest duplicated node indices.
  // It's worth mentioning that such a pair may not be an existing edge.
  // It's a more kind of signature/hash to reconcile duplicated edges.
  // On the other hand, `edgeId` is the minimal (existing) edge index for the given `pairId` signature.
  // Two edges sharing the same `edgeId` are therefore duplicated.
  std::map< std::pair< int, int >, int > pairIdToEdgeId;
  for( std::pair< int, int > const & p: allEdgesFlat )
  {
    std::pair< int, int > const pairId( std::minmax( faceData.uniqueNode( p.first ), faceData.uniqueNode( p.second ) ) );

    auto const & it = pairIdToEdgeId.find( pairId );
    if( it == pairIdToEdgeId.end() )
    {
      pairIdToEdgeId[pairId] = pairToEdge( p, n2ed );
    }
    else
    {
      it->second = std::min( it->second, pairToEdge( p, n2ed ) );
    }
  }
  GEOSX_ERROR_IF_NE_MSG( std::size_t(num2dFaces), pairIdToEdgeId.size(), "Something wrong happened" );

  // This loop build the `elem2dToEdges` mapping.
  // It uses the `pairIdToEdgeId` to remove duplicates.
  for( int i = 0; i < num2dElements; ++i )
  {
    std::set< int > edges;
    for( std::size_t j = 0; j < allEdges[i].size(); ++j )
    {
      for( std::size_t k = 0; k < allEdges[i][j].size(); ++k )
      {
        std::pair< int, int > const & p = allEdges[i][j][k];
        std::pair< int, int > const pairId = std::minmax( faceData.uniqueNode( p.first ), faceData.uniqueNode( p.second ) );
        int const edge = pairIdToEdgeId.at( pairId );
        edges.insert( edge );
      }
    }
    elem2dToEdges.resizeArray( i, edges.size() );
    std::copy( edges.cbegin(), edges.cend(), elem2dToEdges[i].begin() );
  }

  // Here we compute the `face2dToEdges` mapping.
  // We do not pay attention to the ordering of the `face2d`, it will be in the sorting order.
  auto cmp = [&faceData]( std::pair< int, int > const & p0,
                          std::pair< int, int > const & p1 ) -> bool
  {
    std::pair< int, int > const tmp0( std::minmax( faceData.uniqueNode( p0.first ), faceData.uniqueNode( p0.second ) ) );
    std::pair< int, int > const tmp1( std::minmax( faceData.uniqueNode( p1.first ), faceData.uniqueNode( p1.second ) ) );
    return tmp0 < tmp1;
  };
  std::set< std::pair< int, int >, decltype( cmp ) > uniqueEdges( cmp );
  uniqueEdges.insert( allEdgesFlat.cbegin(), allEdgesFlat.cend() );
  GEOSX_ERROR_IF_NE_MSG( LvArray::integerConversion< localIndex >( uniqueEdges.size() ), num2dFaces, "Internal error: wrong number of edges in the fracture/fault mesh computation." );
  for( std::pair< int, int > const & p: uniqueEdges )
  {
    face2dToEdges.emplace_back( pairToEdge( p, n2ed ) );
  }

  return { face2dToEdges, elem2dToEdges };
}

/**
 * @brief Compute the 2d faces (geometrical edges in 2d) to 2d elems (geometrical faces in 3d) mapping.
 * @param num2dFaces The number of 2d faces (geometrical edges in 3d).
 * @param face2dToEdges The 2d faces (geometrical edges in 2d) to edges mapping.
 * @param elem2dToEdges The 2d elem (geometrical faces in 3d) to edges mapping.
 * @return
 */
ArrayOfArrays <localIndex> computeFace2dToElems2d( localIndex const num2dFaces,
                                                   array1d <localIndex> const & face2dToEdges,
                                                   ArrayOfArrays <localIndex> const & elem2dToEdges )
{
  ArrayOfArrays< localIndex > face2dToElems2d( num2dFaces ); // API

  // First, a mapping inversion TODO there are optimized implementations here and there.
  std::map< int, std::vector< int > > edgesToElems2d;
  for( localIndex i = 0; i < elem2dToEdges.size(); ++i )
  {
    for( int const & e: elem2dToEdges[i] )
    {
      edgesToElems2d[e].push_back( i );
    }
  }
  // Then we fill chain mappings...
  for( int fi = 0; fi < num2dFaces; ++fi )
  {
    localIndex const & edge = face2dToEdges[fi];
    std::vector< int > const & elems2d = edgesToElems2d[edge];
    face2dToElems2d.resizeArray( fi, elems2d.size() );
    std::copy( elems2d.begin(), elems2d.end(), face2dToElems2d[fi].begin() );
  }

  return face2dToElems2d;
}

/**
 * @brief Import face block @p faceBlockName from @p vtkMesh into the @p cellBlockManager.
 * @param faceBlockName
 * @param vtkMesh
 * @param cellBlockManager
 * @deprecated This is most likely not entirely the place to do this.
 */
void ImportFracture( string const & faceBlockName,
                     vtkSmartPointer< vtkDataSet > vtkMesh,
                     CellBlockManager & cellBlockManager )
{
  vtkFieldData * fieldData = vtkMesh->GetFieldData();
  GEOSX_ERROR_IF( fieldData == nullptr, "Fault and fracture information is supposed to be defined in the VTK field data. None was found." );

  FaceData const faceData( fieldData, faceBlockName );

  vtkIdType const num2dElements = faceData.num2dElements(); // API
  // Computing the number of 2d faces
  localIndex const num2dFaces = computeNum2dFaces( vtkMesh, faceData ); // API

  internal::ElementToFace const e2f( cellBlockManager.getCellBlocks() );

  ToCellRelation< array2d< localIndex > > elem2dToElems = compute2dElemToElems( faceData, e2f );

  auto tmp = compute2dElemToNodesAndFaces( vtkMesh, faceData, e2f, cellBlockManager.getFaceToNodes() );
  ArrayOfArrays< localIndex > const elem2dToNodes( std::get< 0 >( tmp ) );
  array2d< localIndex > const elem2dToFaces( std::get< 1 >( tmp ) );

  // The order is from the most little edge index to the largest.
  // Since two edges are in the same location, we always take the more little index of the two colocated edges.
  // Computing face 2d -> edges
  auto tmp2 = compute2dFaceAnd2dElemToEdges( vtkMesh, faceData, num2dFaces, cellBlockManager.getNodeToEdges() );
  array1d< localIndex > const face2dToEdges( std::get< 0 >( tmp2 ) ); // API
  ArrayOfArrays< localIndex > const elem2dToEdges( std::get< 1 >( tmp2 ) ); // API

  ArrayOfArrays< localIndex > const face2dToElems2d = computeFace2dToElems2d( num2dFaces, face2dToEdges, elem2dToEdges ); // API

  // Mappings are now computed. Just create the face block by value.
  FaceBlock & faceBlock = cellBlockManager.getFaceBlocks().registerGroup< FaceBlock >( faceBlockName );
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
