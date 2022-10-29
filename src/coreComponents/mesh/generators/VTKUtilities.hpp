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
#include "mesh/DomainPartition.hpp"
#include "mesh/generators/CellBlockManager.hpp"

#include <vtkSmartPointer.h>
#include <vtkDataSet.h>
#include <vtkMultiProcessController.h>

#include <numeric>
#include <unordered_set>

namespace geosx
{

using namespace dataRepository;

namespace vtk
{

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

/**
 * @brief Return a VTK controller for multiprocessing.
 * @return Return a VTK controller for multiprocessing.
 */
vtkSmartPointer< vtkMultiProcessController > getController();

/**
 * @brief Load the VTK file into the VTK data structure
 * @param[in] filePath the Path of the file to load
 * @return a vtk mesh
 */
vtkSmartPointer< vtkDataSet > loadMesh( Path const & filePath );

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
 * @param[in] method the partitionning method
 * @param[in] partitionRefinement number of graph partitioning refinement cycles
 * @param[in] useGlobalIds controls whether global id arrays from the vtk input should be used
 * @return the vtk grid redistributed
 */
vtkSmartPointer< vtkDataSet >
redistributeMesh( vtkDataSet & loadedMesh,
                  MPI_Comm const comm,
                  PartitionMethod const method,
                  int const partitionRefinement,
                  int const useGlobalIds );

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

/**
 * @brief Print statistics for a vtk mesh
 *
 * @param[in] mesh an input mesh
 * @param cellMap a map of cell lists grouped by type
 * @param comm the MPI communicator
 */
void printMeshStatistics( vtkDataSet & mesh,
                          CellMapType const & cellMap,
                          MPI_Comm const comm );

/**
 * @brief Collect the data to be imported.
 * @param[in] mesh an input mesh
 * @param[in] srcFieldNames an array of field names
 * @return A list of pointers to VTK data arrays.
 */
std::vector< vtkDataArray * >
findArraysForImport( vtkDataSet & mesh,
                     arrayView1d< string const > const & srcFieldNames );

} // namespace vtk

/**
 * @brief Build all the vertex blocks.
 * @param[in] mesh The vtkUnstructuredGrid or vtkStructuredGrid that is loaded
 * @param[in] nodesetNames
 * @param[in] cellBlockManager The instance that stores the vertex blocks.
 * @param[in] translate translate the dataset
 * @param[in] scale scale the dataset
 * @return size of the dataset on x-axis
 */
real64 writeNodes( vtkDataSet & mesh,
                   string_array & nodesetNames,
                   CellBlockManager & cellBlockManager,
                   const geosx::R1Tensor & translate,
                   const geosx::R1Tensor & scale );

/**
 * @brief Build all the cell blocks.
 * @param[in] mesh The vtkUnstructuredGrid or vtkStructuredGrid that is loaded
 * @param[in] cellMap Map from the surfaces index to the list of cells in this surface in this rank.
 * @param[in] cellBlockManager The instance that stores the cell blocks.
 */
void writeCells( vtkDataSet & mesh,
                 const geosx::vtk::CellMapType & cellMap,
                 CellBlockManager & cellBlockManager );

/**
 * @brief Build the "surface" node sets from the surface information.
 * @param[in] mesh The vtkUnstructuredGrid or vtkStructuredGrid that is loaded
 * @param[in] cellMap Map from the surfaces index to the list of cells in this surface in this rank.
 * @param[out] cellBlockManager The instance that stores the node sets.
 * @note @p surfacesIdsToCellsIds will contain all the surface ids across all the MPI ranks, but only its cell ids.
 * If the current MPI rank has no cell id for a given surface, then an empty set will be created.
 */
void writeSurfaces( vtkDataSet & mesh,
                    const geosx::vtk::CellMapType & cellMap,
                    CellBlockManager & cellBlockManager );

/**
 * @brief Import data on 3d cells restricted to cells in a specific region
 *
 * @param regionId The id of the region
 * @param elemType The type of the element to store the data
 * @param cellIds The cell ids of the specific region
 * @param elemManager The instance that stores the elements.
 * @param fieldNames An array of the fields names
 * @param srcArrays an array of data to import
 * @param fieldsToBeSync Indentifies the fields to synchronize
 */
void importFieldOnCellElementSubRegion( int const regionId,
                                        ElementType const elemType,
                                        std::vector< vtkIdType > const & cellIds,
                                        ElementRegionManager & elemManager,
                                        arrayView1d< string const > const & fieldNames,
                                        std::vector< vtkDataArray * > const & srcArrays,
                                        FieldIdentifiers & fieldsToBeSync );


} // namespace geosx

#endif /* GEOSX_MESH_GENERATORS_VTKUTILITIES_HPP */
