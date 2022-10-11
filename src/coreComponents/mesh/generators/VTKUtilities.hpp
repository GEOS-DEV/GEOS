/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file VTKUtilities.hpp
 */

#ifndef GEOSX_MESH_GENERATORS_VTKUTILITIES_HPP
#define GEOSX_MESH_GENERATORS_VTKUTILITIES_HPP

#include "common/DataTypes.hpp"
#include "common/DataLayouts.hpp"
#include "common/MpiWrapper.hpp"
#include "common/TypeDispatch.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/MeshBody.hpp"
#include "mesh/generators/CellBlockManager.hpp"
#include "mesh/generators/VTKMeshGeneratorTools.hpp"
#include "mesh/generators/ParMETISInterface.hpp"
#ifdef GEOSX_USE_SCOTCH
#include "mesh/generators/PTScotchInterface.hpp"
#endif

#include <vtkSmartPointer.h>
#include <vtkDataSet.h>
#include <vtkCellArray.h>
#include <vtkMultiProcessController.h>
#include <vtkBoundingBox.h>
#include <vtkDataArray.h>
#include <vtkPointData.h>
#include <vtkPartitionedDataSet.h>
#include <vtkExtractCells.h>

#include <numeric>
#include <unordered_set>

namespace geosx
{

using namespace dataRepository;

/**
 * @brief Choice of advanced mesh partitioner
 */
enum class PartitionMethod : integer
{
  parmetis, ///< Use ParMETIS library
  ptscotch, ///< Use PTScotch library
};

/// Strings for VTKMeshGenerator::PartitionMethod enumeration
ENUM_STRINGS( PartitionMethod,
              "parmetis",
              "ptscotch" );

  /**
   * @brief Type of map used to store cell lists.
   *
   * This should be an unordered_map, but some outdated standard libraries on some systems
   * do not provide std::hash specialization for enums. This is not performance critical though.
   */
using CellMapType = std::map< ElementType, std::unordered_map< int, std::vector< vtkIdType > > >;

namespace vtk
{

/**
 * @brief Return a VTK controller for multiprocessing.
 */
vtkSmartPointer< vtkMultiProcessController > getController();

/**
 * @brief Load the VTK file into the VTK data structure
 * @param[in] filePath the Path of the file to load
 */
vtkSmartPointer< vtkDataSet > loadMesh( Path const & filePath );

/**
 * @brief Get the Cell Array object
 * @details Replaces GetCells() that exist only in vtkUnstructuredGrid
 * @param[in] mesh 
 * @return vtkSmartPointer< vtkCellArray > 
 */
vtkSmartPointer< vtkCellArray > GetCellArray( vtkDataSet & mesh );

vtkSmartPointer< vtkDataSet >
generateGlobalIDs( vtkDataSet & mesh );

vtkSmartPointer< vtkDataSet >
redistributeByCellGraph( vtkDataSet & mesh,
                         PartitionMethod const method,
                         MPI_Comm const comm,
                         int const numRefinements );

vtkSmartPointer< vtkDataSet >
redistributeByKdTree( vtkDataSet & mesh );

/**
 * @brief Compute the rank neighbor candidate list.
 * @param[in] boundingBoxes the bounding boxes used by the VTK partitioner for all ranks
 * @return the list of neighboring MPI ranks, will be updated
 */
std::vector< int >
findNeighborRanks( std::vector< vtkBoundingBox > boundingBoxes );

/**
 * @brief Generate global point/cell IDs and redistribute the mesh among MPI ranks.
 * @param[in] loadedMesh the mesh that was loaded on one or several MPI ranks
 * @param[in] comm the MPI communicator
 * @param[in] partitionRefinement number of graph partitioning refinement cycles
 * @param[in] useGlobalIds controls whether global id arrays from the vtk input should be used
 */
vtkSmartPointer< vtkDataSet >
redistributeMesh( vtkDataSet & loadedMesh,
                  MPI_Comm const comm,
                  PartitionMethod const method,
                  int const partitionRefinement,
                  int const useGlobalIds );

/**
 * @brief Identify the GEOSX type of the polyhedron
 *
 * @param cell The vtk cell VTK_POLYHEDRON
 * @return The geosx element type associated to VTK_POLYHEDRON
 */
geosx::ElementType buildGeosxPolyhedronType( vtkCell * const cell );

/**
 * @brief Get the GEOSX element type
 * @param[in] cellType The vtk cell type
 * @return The GEOSX element type
 */
ElementType convertVtkToGeosxElementType( VTKCellType const cellType );

std::map< ElementType, std::vector< vtkIdType > >
splitCellsByType( vtkDataSet & mesh );

CellMapType
splitCellsByTypeAndAttribute( std::map< ElementType, std::vector< vtkIdType > > & typeToCells,
                              vtkDataArray * const attributeDataArray );

void extendCellMapWithRemoteKeys( CellMapType & cellMap );

/**
 * @brief Get the tetrahedron node ordering from a VTK_POLYHEDRON
 *
 * @param cell The vtk cell, type VTK_POLYHEDRON
 * @return The node ordering
 */
std::vector< localIndex > getTetrahedronNodeOrderingFromPolyhedron( vtkCell * const cell );

/**
 * @brief Get the hexahedron node ordering from a VTK_POLYHEDRON
 *
 * @param cell The vtk cell, type VTK_POLYHEDRON
 * @return The node ordering
 *
 * It could be possible to use getPrismNodeOrderingFromPolyhedron< 4 > with additional
 * permutations. But at this point computationalGeometry::prismVolume< NUM_SIDES >
 * is not ready.
 */
std::vector< localIndex > getHexahedronNodeOrderingFromPolyhedron( vtkCell * const cell );

/**
 * @brief Get the wedge node ordering from a VTK_POLYHEDRON
 *
 * @param cell The vtk cell, type VTK_POLYHEDRON
 * @return The node ordering
 *
 * It could be possible to use getPrismNodeOrderingFromPolyhedron< 3 > with additional
 * permutations. But at this point computationalGeometry::prismVolume< NUM_SIDES >
 * is not ready.
 */
std::vector< localIndex > getWedgeNodeOrderingFromPolyhedron( vtkCell * const cell );

/**
 * @brief Get the pyramid node ordering from a VTK_POLYHEDRON
 *
 * @param cell The vtk cell, type VTK_POLYHEDRON
 * @return The node ordering
 */
std::vector< localIndex > getPyramidNodeOrderingFromPolyhedron( vtkCell * const cell );

/**
 * @brief Get the prism node ordering from a VTK_POLYHEDRON
 *
 * @tparam NUM_SIDES number of sides of the prism
 * @param cell The vtk cell, type VTK_POLYHEDRON
 * @return The node ordering
 */
template< integer NUM_SIDES >
std::vector< localIndex > getPrismNodeOrderingFromPolyhedron( vtkCell * const cell );

/**
 * @brief Collect lists of VTK cell indices organized by type and attribute value.
 * @param[in] mesh the vtkUnstructuredGrid or vtkStructuredGrid that is loaded
 * @param[in] attributeName name of the VTK data array containing the attribute, if any
 * @return A map from element type to a map of attribute to the associated cell ids for the current rank.
 *         The map contains entries for all types and attribute values across all MPI ranks,
 *         even if there are no cells on current rank (then the list will be empty).
 */
CellMapType buildCellMap( vtkDataSet & mesh,
                          string const & attributeName );

std::vector< int > getGeosxToVtkNodeOrdering( ElementType const elemType );

std::vector< int > getGeosxToVtkNodeOrdering( ElementType const elemType, vtkCell *cell );

/**
 * @brief Fill @p cellBlock with the appropriate nodes and local/global mappings.
 * @param[in] elemType the vtk cell type for cells of the CellBlock being written
 * @param[in] cellIds the cell indexes of cell type \p cellType within this region
 * @param[in] mesh the vtkUnstructuredGrid or vtkStructuredGrid that is loaded
 * @param[in,out] cellBlock The cell block to be written
 */
void fillCellBlock( vtkDataSet & mesh,
                    ElementType const elemType,
                    std::vector< vtkIdType > const & cellIds,
                    CellBlock & cellBlock );

/**
 * @brief Returns a string describing the element.
 * @param[in] type The element type.
 * @return The name.
 * @warning This information will be visible in the input file... Consider refactoring with great care.
 */
string getElementTypeName( ElementType const type );

/**
 * @brief Builds the cell block name.
 * @param[in] elementName The element name.
 * @param[in] regionId The region Id.
 * @return The name.
 * @warning This name will be visible in the input file... Consider refactoring with great care.
 */
string buildCellBlockName( ElementType const type, int const regionId );

/**
 * @brief Collect a set of material field names registered in a subregion.
 * @param subRegion the target subregion
 * @return a set of wrapper names
 */
std::unordered_set< string > getMaterialWrapperNames( ElementSubRegionBase const & subRegion );

/**
 * @brief Imports 2d and 3d arrays from @p vtkArray to @p wrapper, only for @p cellIds
 * @param cellIds The cells for which we should copy the data.
 * @param vtkArray The source.
 * @param wrapper The destination.
 */
void importMaterialField( std::vector< vtkIdType > const & cellIds,
                          vtkDataArray * vtkArray,
                          WrapperBase & wrapper );

/**
 * @brief Imports 1d and 2d arrays from @p vtkArray to @p wrapper, only for @p cellIds
 * @param cellIds The cells for which we should copy the data.
 * @param vtkArray The source.
 * @param wrapper The destination.
 */
void importRegularField( std::vector< vtkIdType > const & cellIds,
                         vtkDataArray * vtkArray,
                         WrapperBase & wrapper );

void printMeshStatistics( vtkDataSet & mesh,
                          CellMapType const & cellMap,
                          MPI_Comm const comm );

/**
 * @brief Collect the data to be imported.
 * @return A list of pointers to VTK data arrays.
 */
std::vector< vtkDataArray * >
findArraysForImport( vtkDataSet & mesh,
                     arrayView1d< string const > const & srcFieldNames );


/**
 * @brief Build the Elements to Nodes indexing
 * 
 * @tparam INDEX_TYPE 
 * @tparam POLICY 
 * @param mesh 
 * @param cells 
 * @return ArrayOfArrays< INDEX_TYPE, INDEX_TYPE > 
 */
template< typename INDEX_TYPE, typename POLICY >
ArrayOfArrays< INDEX_TYPE, INDEX_TYPE >
buildElemToNodesImpl( vtkDataSet & mesh,
                      vtkSmartPointer< vtkCellArray > const & cells )
{
  localIndex const numCells = LvArray::integerConversion< localIndex >( mesh.GetNumberOfCells() );
  array1d< INDEX_TYPE > nodeCounts( numCells );

  // GetCellSize() is always thread-safe, can run in parallel
  forAll< parallelHostPolicy >( numCells, [nodeCounts = nodeCounts.toView(), &cells] ( localIndex const cellIdx )
  {
    nodeCounts[cellIdx] = LvArray::integerConversion< INDEX_TYPE >( cells->GetCellSize( cellIdx ) );
  } );

  ArrayOfArrays< INDEX_TYPE, INDEX_TYPE > elemToNodes;
  elemToNodes.template resizeFromCapacities< parallelHostPolicy >( numCells, nodeCounts.data() );

  vtkIdTypeArray const & globalPointId = *vtkIdTypeArray::FastDownCast( mesh.GetPointData()->GetGlobalIds() );

  // GetCellAtId() is conditionally thread-safe, use POLICY argument
  forAll< POLICY >( numCells, [&cells, &globalPointId, elemToNodes = elemToNodes.toView()] ( localIndex const cellIdx )
  {
    vtkIdType numPts;
    vtkIdType const * points;
    cells->GetCellAtId( cellIdx, numPts, points );
    for( int a = 0; a < numPts; ++a )
    {
      vtkIdType const pointIdx = globalPointId.GetValue( points[a] );
      elemToNodes.emplaceBack( cellIdx, LvArray::integerConversion< INDEX_TYPE >( pointIdx ) );
    }
  } );

  return elemToNodes;
}

template< typename INDEX_TYPE >
ArrayOfArrays< INDEX_TYPE, INDEX_TYPE >
buildElemToNodes( vtkDataSet & mesh )
{
  vtkSmartPointer< vtkCellArray > const & cells = vtk::GetCellArray( mesh );
  // According to VTK docs, IsStorageShareable() indicates whether pointers extracted via
  // vtkCellArray::GetCellAtId() are pointers into internal storage rather than temp buffer
  // and thus results can be used in a thread-safe way.
  return cells->IsStorageShareable()
       ? buildElemToNodesImpl< INDEX_TYPE, parallelHostPolicy >( mesh, cells )
       : buildElemToNodesImpl< INDEX_TYPE, serialPolicy >( mesh, cells );
}

template< typename PART_INDEX >
vtkSmartPointer< vtkPartitionedDataSet >
splitMeshByPartition( vtkDataSet & mesh,
                      PART_INDEX const numParts,
                      arrayView1d< PART_INDEX const > const & part )
{
  array1d< localIndex > cellCounts( numParts );
  forAll< parallelHostPolicy >( part.size(), [part, cellCounts = cellCounts.toView()] ( localIndex const cellIdx )
  {
    RAJA::atomicInc< parallelHostAtomic >( &cellCounts[part[cellIdx]] );
  } );

  ArrayOfArrays< vtkIdType > cellsLists;
  cellsLists.resizeFromCapacities< serialPolicy >( numParts, cellCounts.data() );

  forAll< parallelHostPolicy >( part.size(), [part, cellsLists = cellsLists.toView()] ( localIndex const cellIdx )
  {
    cellsLists.emplaceBackAtomic< parallelHostAtomic >( LvArray::integerConversion< localIndex >( part[cellIdx] ),
                                                        LvArray::integerConversion< vtkIdType >( cellIdx ) );
  } );

  vtkNew< vtkPartitionedDataSet > result;
  result->SetNumberOfPartitions( LvArray::integerConversion< unsigned int >( numParts ) );

  vtkNew< vtkExtractCells > extractor;
  extractor->SetInputDataObject( &mesh );

  for( localIndex p = 0; p < numParts; ++p )
  {
    arraySlice1d< vtkIdType const > const cells = cellsLists[p];
    if( cells.size() > 0 )
    {
      extractor->SetCellIds( cells.dataIfContiguous(), LvArray::integerConversion< vtkIdType >( cells.size() ) );
      extractor->Update();

      vtkNew< vtkUnstructuredGrid > ug;
      ug->ShallowCopy( extractor->GetOutputDataObject( 0 ) );
      result->SetPartition( LvArray::integerConversion< unsigned int >( p ), ug );
    }
  }
  return result;
}

/**
 * @brief Gathers all the data from all ranks, merge them, sort them, and remove duplicates.
 * @tparam T Type of the exchanged data.
 * @param data The data to be exchanged.
 * @return The merged data.
 * @note This function makes MPI calls.
 */
template< typename T >
std::vector< T > collectUniqueValues( std::vector< T > const & data )
{
  // Exchange the sizes of the data across all ranks.
  array1d< int > dataSizes( MpiWrapper::commSize() );
  MpiWrapper::allGather( LvArray::integerConversion< int >( data.size() ), dataSizes, MPI_COMM_GEOSX );
  // `totalDataSize` contains the total data size across all the MPI ranks.
  int const totalDataSize = std::accumulate( dataSizes.begin(), dataSizes.end(), 0 );

  // Once the MPI exchange is done, `allData` will contain all the data of all the MPI ranks.
  // We want all ranks to get all the data. But each rank may have a different size of information.
  // Therefore, we use `allgatherv` that does not impose the same size across ranks like `allgather` does.
  std::vector< T > allData( totalDataSize );
  // `displacements` is the offset (relative to the receive buffer) to store the data for each rank.
  std::vector< int > displacements( MpiWrapper::commSize(), 0 );
  std::partial_sum( dataSizes.begin(), dataSizes.end() - 1, displacements.begin() + 1 );
  MpiWrapper::allgatherv( data.data(), data.size(), allData.data(), dataSizes.data(), displacements.data(), MPI_COMM_GEOSX );

  // Finalizing by sorting, removing duplicates and trimming the result vector at the proper size.
  std::sort( allData.begin(), allData.end() );
  auto newEnd = std::unique( allData.begin(), allData.end() );
  allData.erase( newEnd, allData.end() );

  return allData;
}


} // namespace vtk

/**
 * @brief Build all the vertex blocks.
 * @param[in] cellBlockManager The instance that stores the vertex blocks.
 */
real64 writeNodes( vtkDataSet & mesh,
                   string_array& nodesetNames,
                   CellBlockManager & cellBlockManager,
                   const geosx::R1Tensor& translate,
                   const geosx::R1Tensor& scale );


void importNodesets( vtkDataSet & mesh, 
                     string_array& nodesetNames,
                     CellBlockManager & cellBlockManager );

/**
 * @brief Build all the cell blocks.
 * @param[in] cellBlockManager The instance that stores the cell blocks.
 */
void writeCells( vtkDataSet & mesh,
                 const geosx::CellMapType &cellMap,
                 CellBlockManager & cellBlockManager );

/**
 * @brief Build the "surface" node sets from the surface information.
 * @param[in] mesh The vtkUnstructuredGrid or vtkStructuredGrid that is loaded
 * @param[in] surfacesIdsToCellsIds Map from the surfaces index to the list of cells in this surface in this rank.
 * @param[out] cellBlockManager The instance that stores the node sets.
 * @note @p surfacesIdsToCellsIds will contain all the surface ids across all the MPI ranks, but only its cell ids.
 * If the current MPI rank has no cell id for a given surface, then an empty set will be created.
 */
void writeSurfaces( vtkDataSet & mesh,
                    const geosx::CellMapType &cellMap,
                    CellBlockManager & cellBlockManager );

void importFieldOnCellElementSubRegion( int const regionId,
                                        ElementType const elemType,
                                        std::vector< vtkIdType > const & cellIds,
                                        ElementRegionManager & elemManager,
                                        arrayView1d< string const > const & fieldNames,
                                        std::vector< vtkDataArray * > const & srcArrays,
                                        FieldIdentifiers & fieldsToBeSync );


} // namespace geosx

#endif /* GEOSX_MESH_GENERATORS_VTKUTILITIES_HPP */

