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

#include "mesh/generators/VTKUtilities.hpp"

#include "dataRepository/Group.hpp"

#include <vtkCell.h>
#include <vtkDataArray.h>
#include <vtkFieldData.h>

#include <vtkExtractEdges.h>
#include <vtkGeometryFilter.h>
#include <vtkPointData.h>
#include <vtkCellData.h>

#include <algorithm>

namespace geosx
{

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
      for( int j = 0, start = 0; j < size - 1; ++j )
      {
        if( l2g[j] + 1 != l2g[j + 1] )
        {
          m_ranges[{ l2g[start], l2g[j] + 1 }] = { i, -1 };
          start = j + 1;
        }
        else if( j + 1 == size - 1 )
        {
          m_ranges[{ l2g[start], l2g[j + 1] + 1 }] = { i, -1 };
        }
      }
    }

    std::vector< int > cbOffset( numCellBlocks, 0 );
    for( std::pair< Range const, CellBlockInfo > & p: m_ranges )
    {
      Range const & range = p.first;
      CellBlockInfo & cb = p.second;
      size_t const cbi = cb.index;
      int & offset = cb.offset;

      offset = cbOffset[cbi];
      for( size_t i = 0; i < cbOffset.size(); ++i )
      {
        if( i != cbi )
        {
          cbOffset[i] += range.max - range.min;
        }
      }
    }
  }

  /**
   * @brief Given the global (vtk) id, returns the cell block index @p ei belongs to.
   * @param[in] ei global cell index.
   * @return The cell block index.
   */
  int getCellBlockIndex( int const & ei ) const
  {
    return getCellBlockInfo( ei ).index;
  }

  /**
   * @brief From global index @p ei, returns the local element index in its cell block.
   * @param ei The global vtk element index.
   * @return The integer.
   */
  localIndex getElementIndexInCellBlock( int const & ei ) const
  {
    CellBlockInfo const & cbInfo = getCellBlockInfo( ei );
    return ei - cbInfo.offset;
  }

  /**
   * @brief Given the global (vtk) id, returns the index local to the cell block @p ei belong to. DAFUCK?
   * @param[in] ei global cell index.
   * @return The faces for element @p ei.
   */
  auto operator[]( int const & ei ) const
  {
    CellBlockInfo const & cbInfo = getCellBlockInfo( ei );
    arrayView2d< localIndex const > const & e2f = m_elemToFacesMaps[cbInfo.index];
    return e2f[ei - cbInfo.offset];
  }

private:

  struct Range
  {
    globalIndex min;
    globalIndex max;

    bool operator<( Range const & other ) const
    {
      return std::tie( min, max ) < std::tie( other.min, other.max );
    }
  };

  struct CellBlockInfo
  {
    int index;
    int offset;
  };

  /**
   * @brief Return the cell block information for element @p ei.
   * @param ei The global index of the element.
   * @return The cell block information.
   *
   * The cell block information contains which cell block (index) contains the global element that matches the @p range.
   * It also contains which offset should be subtracted in order to compute the local (to the cell block) indexing of the global element.
   */
  CellBlockInfo const & getCellBlockInfo( int const & ei ) const
  {
    auto const it = std::find_if( m_ranges.cbegin(), m_ranges.cend(), [&]( std::pair< Range, CellBlockInfo > const & p ) -> bool
    {
      return p.first.min <= ei && ei < p.first.max;
    } );

    return it->second;
  }

  /**
   * @brief Stores the cell block information for each given global element index range (lower included, upper excluded).
   */
  std::map< Range, CellBlockInfo > m_ranges;
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


class DuplicatedPoints
{
public:
  DuplicatedPoints( vtkSmartPointer< vtkDataSet > faceMesh )
  {
    vtkDataArray * duplicatedNodes = faceMesh->GetPointData()->GetArray( "duplicated_nodes" );  // TODO cast to int array?

    vtkIdType const numTuples = duplicatedNodes->GetNumberOfTuples();
    int const numComponents = duplicatedNodes->GetNumberOfComponents();
    m_duplicatedNodes.resize( numTuples );
    for( vtkIdType i = 0; i < numTuples; ++i )
    {
      m_duplicatedNodes[i].reserve( numComponents );
    }

    for( vtkIdType i = 0; i < numTuples; ++i )
    {
      double const * const tuple = duplicatedNodes->GetTuple( i );
      for( int j = 0; j < numComponents; ++j )
      {
        double const tmp = tuple[j];
        if( tmp > -1 )
        {
          m_duplicatedNodes[i].emplace_back( tmp );
        }
      }
    }
  }

  std::vector< vtkIdType > const & operator()( int i ) const // TODO operator[]
  {
    return m_duplicatedNodes[i];
  }

  std::size_t size() const
  {
    return m_duplicatedNodes.size();
  }

private:
  std::vector< std::vector< vtkIdType > > m_duplicatedNodes;
};

struct pairHash
{
  size_t operator()( const std::pair< vtkIdType, vtkIdType > & p ) const
  {
    auto hash1 = std::hash< vtkIdType >{}( p.first );
    auto hash2 = std::hash< vtkIdType >{}( p.second );

    // If hash1 == hash2, their XOR is zero.
    return hash1 != hash2 ? hash1 ^ hash2: hash1;
  }
};


/**
 * @brief Compute the 2d faces (geometrical edges in 2d) to 2d elems (geometrical faces in 3d) mapping.
 * @param edges The edges of the fracture mesh.
 * @param faceMesh The fracture mesh.
 * @return The mapping
 */
ArrayOfArrays< localIndex > computeFace2dToElems2d( vtkPolyData * edges,
                                                    vtkSmartPointer< vtkDataSet > faceMesh )
{
  std::unordered_map< std::pair< vtkIdType, vtkIdType >, int, pairHash > face2dIds;
  for( int i = 0; i < edges->GetNumberOfCells(); ++i )
  {
    vtkCell * c = edges->GetCell( i );
    vtkIdList * edgePointIds = c->GetPointIds();
    GEOSX_ASSERT( edgePointIds->GetNumberOfIds() == 2 );
    std::pair< vtkIdType, vtkIdType > const minMax = std::minmax( edgePointIds->GetId( 0 ), edgePointIds->GetId( 1 ) );
    face2dIds[minMax] = i;
  }

  ArrayOfArrays< localIndex > face2dToElems2d( edges->GetNumberOfCells() );
  for( int i = 0; i < faceMesh->GetNumberOfCells(); ++i )
  {
    vtkCell * c = faceMesh->GetCell( i );
    for( int j = 0; j < c->GetNumberOfEdges(); ++j )
    {
      vtkCell * e = c->GetEdge( j );
      vtkIdList * edgePointIds = e->GetPointIds();
      GEOSX_ASSERT( edgePointIds->GetNumberOfIds() == 2 );
      std::pair< vtkIdType, vtkIdType > const minMax = std::minmax( edgePointIds->GetId( 0 ), edgePointIds->GetId( 1 ) );
      face2dToElems2d.emplaceBack( face2dIds.at( minMax ), i ); // TODO cast
    }
  }

  return face2dToElems2d;
}

array1d< localIndex > computeFace2dToEdge( vtkPolyData * edges,
                                           DuplicatedPoints const & duplicatedPoints,
                                           ArrayOfArraysView< localIndex const > nodeToEdges )
{
  // Computing the face2dToEdges mapping. It does not require face2dToElems2d: split the function.
  auto const comp = []( std::pair< vtkIdType, int > const & l, std::pair< vtkIdType, int > const & r ) { return l.second < r.second; };
  array1d< localIndex > face2dToEdge( edges->GetNumberOfCells() );
  for( int i = 0; i < edges->GetNumberOfCells(); ++i )
  {
    std::vector< vtkIdType > allDuplicatedNodesOfEdge;
    vtkCell * edge = edges->GetCell( i );
    for( int j = 0; j < edge->GetNumberOfPoints(); ++j )
    {
      for( auto const & d: duplicatedPoints( edge->GetPointId( j ) ) )
      {
        allDuplicatedNodesOfEdge.emplace_back( d );
      }
    }
    std::map< vtkIdType, int > edgeCount;
    for( auto const & d: allDuplicatedNodesOfEdge )
    {
      auto const & es = nodeToEdges[d];
      for( int e = 0; e< es.size(); ++e )
      {
        edgeCount[es[e]]++;
      }
    }
    auto const res = std::max_element( edgeCount.cbegin(), edgeCount.cend(), comp );
    face2dToEdge[i] = res->first; // TODO cast
  }

  return face2dToEdge;
}

ArrayOfArrays< localIndex > computeElem2dToFace2d( vtkIdType num2dElements,
                                                   ArrayOfArraysView< localIndex const > face2dToElems2d )
{
  // Let's first invert face2dToElems2d into elem2dToFace2d
  ArrayOfArrays< localIndex > elem2dToFace2d( num2dElements );
  for( localIndex face2dIndex = 0; face2dIndex < face2dToElems2d.size(); ++face2dIndex )
  {
    for( localIndex const & elem2dIndex: face2dToElems2d[face2dIndex] )
    {
      elem2dToFace2d.emplaceBack(elem2dIndex, face2dIndex );  // TODO cast
    }
  }

  return elem2dToFace2d;
}

ArrayOfArrays< localIndex > computeElem2dToEdges( vtkIdType num2dElements,
                                                  arrayView1d< localIndex const > face2dToEdge,
                                                  ArrayOfArraysView< localIndex const > elem2dToFace2d )
{
  // Computing elem2dToEdges.
  ArrayOfArrays< localIndex > elem2dToEdges( num2dElements );
  // Let's first invert face2dToElems2d into elem2dToFace2d
  // Now that we have elem2dToFace2d, we can compose with face2dToEdge to create elem2dToEdges.
  for( localIndex elemIndex = 0; elemIndex < elem2dToFace2d.size(); ++elemIndex )
  {
    for( auto const & face2dIndex: elem2dToFace2d[elemIndex] )
    {
      elem2dToEdges.emplaceBack(elemIndex, face2dToEdge[face2dIndex] );
    }
  }

  return elem2dToEdges;
}

struct Elem2dTo3dInfo
{
  ToCellRelation< array2d< localIndex > > elem2dToElem3d;
  array2d< localIndex > elem2dToFaces;

  Elem2dTo3dInfo( ToCellRelation< array2d< localIndex > > && elem2dToElem3d_,
                  array2d< localIndex > && elem2dToFaces_ )
    : elem2dToElem3d( elem2dToElem3d_ ),
      elem2dToFaces( elem2dToFaces_ )
  { }
};

Elem2dTo3dInfo computeElem2dTo3dElemAndFaces( vtkSmartPointer< vtkDataSet > faceMesh,
                                              vtkSmartPointer< vtkDataSet > mesh,
                                              DuplicatedPoints const & duplicatedPoints,
                                              ArrayOfArraysView< localIndex const > faceToNodes,
                                              geosx::internal::ElementToFace const & elemToFaces )
{
  const char * const boundaryPointsName = "boundary points";
  const char * const boundaryCellsName = "boundary cells";

  auto boundaryExtractor = vtkSmartPointer< vtkGeometryFilter >::New();
  boundaryExtractor->PassThroughPointIdsOn();
  boundaryExtractor->PassThroughCellIdsOn();
  boundaryExtractor->FastModeOff();
  boundaryExtractor->SetOriginalPointIdsName( boundaryPointsName );
  boundaryExtractor->SetOriginalCellIdsName( boundaryCellsName );
  boundaryExtractor->SetInputData( mesh );
  boundaryExtractor->Update();
  vtkPolyData * boundary = boundaryExtractor->GetOutput();

  vtkIdTypeArray const * boundaryPoints = vtkIdTypeArray::FastDownCast( boundary->GetPointData()->GetArray( boundaryPointsName ) );
  vtkIdTypeArray const * boundaryCells = vtkIdTypeArray::FastDownCast( boundary->GetCellData()->GetArray( boundaryCellsName ) );

  // Building the elem2d to elem3d mapping. We need to find the 3d elements!
  std::map< vtkIdType, std::vector< vtkIdType > > nodesToCellsFull;
  for( vtkIdType i = 0; i < boundary->GetNumberOfCells(); ++i )
  {
    vtkIdType const cellId = boundaryCells->GetValue( i );
    vtkIdList * pointIds = boundary->GetCell( i )->GetPointIds();
    for( int j = 0; j < pointIds->GetNumberOfIds(); ++j )
    {
      vtkIdType const pointId = boundaryPoints->GetValue( pointIds->GetId( j ) );
      nodesToCellsFull[pointId].emplace_back( cellId );
    }
  }

  // Only keeping the nodes to cells for the duplicated nodes. Small optimisation.
  std::map< vtkIdType, std::set< vtkIdType > > nodesToCells;  // All the nodes in the duplicates and the cells they're "part" of.
  { // scope reduction
    std::set< vtkIdType > allDuplicatedNodes;
    for( std::size_t i = 0; i < duplicatedPoints.size(); ++i )
    {
      std::vector< vtkIdType > const & ns = duplicatedPoints( i );
      allDuplicatedNodes.insert( ns.cbegin(), ns.cend() );
    }

    for( vtkIdType const & n: allDuplicatedNodes )
    {
      std::vector< vtkIdType > const & tmp = nodesToCellsFull.at( n );
      std::set< vtkIdType > const cells{ tmp.cbegin(), tmp.cend() };
      nodesToCells[n] = cells;
    }
  }

  vtkIdType const num2dElements = faceMesh->GetNumberOfCells();
  array2d< localIndex > elem2dToElem3d( num2dElements, 2 );
  array2d< localIndex > elem2dToCellBlock( num2dElements, 2 );
  array2d< localIndex > elem2dToFaces( num2dElements, 2 );
  elem2dToElem3d.setValues< serialPolicy >( -1 );
  elem2dToCellBlock.setValues< serialPolicy >( -1 );
  elem2dToFaces.setValues< serialPolicy >( -1 );
  // Now we loop on all the 2d elements.
  for( int i = 0; i < num2dElements; ++i )
  {
    vtkIdList * pointIds = faceMesh->GetCell( i )->GetPointIds();
    std::size_t const elem2dNumPoints = pointIds->GetNumberOfIds();
    std::set< vtkIdType > duplicatedPointOfElem2d; // All the duplicated points of the 2d element. Note that we lose the collocation of the duplicated nodes.
    for( vtkIdType j = 0; j < pointIds->GetNumberOfIds(); ++j )
    {
      std::vector< vtkIdType > const & ns = duplicatedPoints( pointIds->GetId( j ) );
      duplicatedPointOfElem2d.insert( ns.cbegin(), ns.cend() );
    }

    // Now we search for the 3d elements (by finding its cell)
    // We first extract all the cells that are concerned by those duplicated nodes.
    std::map< vtkIdType, std::set< vtkIdType > > elem3dToDuplicatedNodes;
    for( vtkIdType const & n: duplicatedPointOfElem2d )
    {
      for( vtkIdType const & c: nodesToCells.at( n ) )
      {
        elem3dToDuplicatedNodes[c].insert( n );
      }
    }
    // Then, if we find a
    // Now cell2num will tell which cell is concerned by the duplicated nodes.
    for( auto const & e2n: elem3dToDuplicatedNodes )
    {
      // If the face of the element 3d has the same number of nodes than the elem 2d, it should be OK (the mesh is conformal). FIXME
      if( e2n.second.size() == elem2dNumPoints )
      {
        GEOSX_LOG( e2n.first << " is a bordering 3d element for 2d element " << i );
        // Now we know that the element 3d has a face that touches the element 2d. Let's find which one.
        localIndex const idx = elem2dToElem3d[i][0] == -1 ? 0 : 1;
        elem2dToElem3d[i][idx] = elemToFaces.getElementIndexInCellBlock( e2n.first );  // Warning about the element 3d being global id.
        // Computing the elem2dToFaces mapping.
        auto faces = elemToFaces[e2n.first];
        for( int faceIndex = 0; faceIndex < faces.size( 0 ); ++faceIndex )
        {
          localIndex const f = faces[faceIndex];
          auto nodes = faceToNodes[f];
          if( std::set< vtkIdType >( nodes.begin(), nodes.end() ) == e2n.second )
          {
            elem2dToFaces[i][idx] = f;
            elem2dToCellBlock[i][idx] = elemToFaces.getCellBlockIndex( e2n.first );
            break;
          }
        }
      }
    }
  }

  auto cellRelation = ToCellRelation< array2d< localIndex > >( std::move( elem2dToCellBlock ), std::move( elem2dToElem3d ) );
  return Elem2dTo3dInfo( std::move( cellRelation ), std::move( elem2dToFaces ) );
}


ArrayOfArrays< localIndex > computeElem2dToNodes( vtkIdType num2dElements,
                                                  ArrayOfArraysView< localIndex const > faceToNodes,
                                                  arrayView2d< localIndex const > elem2dToFaces )
{
  ArrayOfArrays< localIndex > elem2dToNodes( num2dElements );
  for( localIndex elem2dIndex = 0; elem2dIndex < elem2dToFaces.size( 0 ); ++elem2dIndex )
  {
    for( localIndex const & faceIndex: elem2dToFaces[elem2dIndex] )
    {
      std::set< localIndex > tmp;
      for( auto j = 0; j < faceToNodes[faceIndex].size(); ++j )
      {
        localIndex const & nodeIndex = faceToNodes[faceIndex][j];
        tmp.insert( nodeIndex );
      }
      for( localIndex const & nodeIndex: tmp )
      {
        elem2dToNodes.emplaceBack( elem2dIndex, nodeIndex );
      }
    }
  }

  return elem2dToNodes;
}


void importFractureNetwork( Path const & filePath,
                            string const & faceBlockName,
                            vtkSmartPointer< vtkDataSet > vtkMesh,
                            CellBlockManager & cellBlockManager )
{
  ArrayOfArrays< localIndex > const faceToNodes = cellBlockManager.getFaceToNodes();
  geosx::internal::ElementToFace const elemToFaces( cellBlockManager.getCellBlocks() );
  ArrayOfArrays< localIndex > const nodeToEdges = cellBlockManager.getNodeToEdges();

  vtkSmartPointer< vtkDataSet > faceMesh = vtk::loadMesh( filePath, faceBlockName );
  DuplicatedPoints const duplicatedPoints( faceMesh );
  // Add the appropriate validations (only 2d cells...)

  auto edgesExtractor = vtkSmartPointer< vtkExtractEdges >::New();
  edgesExtractor->SetInputData( faceMesh );
  edgesExtractor->UseAllPointsOn();  // Important: we want to prevent any node renumbering.
  edgesExtractor->Update();
  vtkPolyData * edges = edgesExtractor->GetOutput();

  vtkIdType const num2dFaces = edges->GetNumberOfCells();
  vtkIdType const num2dElements = faceMesh->GetNumberOfCells();
  // Now let's build the elem2dTo* mappings.
  Elem2dTo3dInfo elem2dTo3d = computeElem2dTo3dElemAndFaces( faceMesh, vtkMesh, duplicatedPoints, faceToNodes.toViewConst(), elemToFaces );
  ArrayOfArrays< localIndex > elem2dToNodes = computeElem2dToNodes( num2dElements, faceToNodes.toViewConst(), elem2dTo3d.elem2dToFaces.toViewConst() );

  ArrayOfArrays< localIndex > face2dToElems2d = computeFace2dToElems2d( edges, faceMesh );
  array1d< localIndex > face2dToEdge = computeFace2dToEdge( edges, duplicatedPoints, nodeToEdges.toViewConst() );
  ArrayOfArrays< localIndex > const elem2dToFace2d = computeElem2dToFace2d( num2dElements, face2dToElems2d.toViewConst() );
  ArrayOfArrays< localIndex > elem2DToEdges = computeElem2dToEdges( num2dElements, face2dToEdge.toViewConst(), elem2dToFace2d.toViewConst() );

  // Mappings are now computed. Just create the face block by value.
  FaceBlock & faceBlock = cellBlockManager.registerFaceBlock( faceBlockName );
  faceBlock.setNum2DElements( num2dElements );
  faceBlock.setNum2DFaces( num2dFaces );
  faceBlock.set2dElemToNodes( std::move( elem2dToNodes ) );
  faceBlock.set2dElemToEdges( std::move( elem2DToEdges ) );
  faceBlock.set2dFaceToEdge( std::move( face2dToEdge ) );
  faceBlock.set2dFaceTo2dElems( std::move( face2dToElems2d ) );
  faceBlock.set2dElemToFaces( std::move( elem2dTo3d.elem2dToFaces ) );
  faceBlock.set2dElemToElems( std::move( elem2dTo3d.elem2dToElem3d ) );
}

} // end of namespace
