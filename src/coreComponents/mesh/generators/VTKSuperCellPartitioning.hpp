/**
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#ifndef GEOS_MESH_GENERATORS_VTKSUPERCELLPARTITIONING_HPP
#define GEOS_MESH_GENERATORS_VTKSUPERCELLPARTITIONING_HPP


#include "common/MpiWrapper.hpp"
#include "common/TimingMacros.hpp"
#include "mesh/generators/VTKMeshGeneratorTools.hpp"  // for PartitionMethod


#include <vtkSmartPointer.h>

// GEOS types
#include "common/DataTypes.hpp"
#include "LvArray/src/ArrayOfArrays.hpp"

// VTK forward declarations
class vtkUnstructuredGrid;
class vtkDataSet;

// Get pmet_idx_t definition from ParMETIS interface
#include "mesh/generators/ParMETISInterface.hpp"  // Defines pmet_idx_t as int64_t



namespace geos
{

namespace vtk
{

/**
 * @brief Metadata about super-cells for partitioning
 * 
 * A super-cell is a group of 3D cells that must stay together during partitioning.
 * Typically these are cells connected by fractures.
 */
struct SuperCellInfo
{
  /// Map: SuperCellId -> vector of global cell IDs in that super-cell
  std::map< vtkIdType, std::vector< vtkIdType > > superCellToOriginalCells;
  
  /// Map: SuperCellId -> weight (number of cells in super-cell)
  std::map< vtkIdType, vtkIdType > vertexWeights;
  
  /// Set of SuperCellIds that contain multiple cells (atomic units)
  std::set< vtkIdType > atomicSuperCells;
};

/**
 * @brief Tag 3D cells with super-cell IDs based on fracture connectivity
 * 
 * Creates a "SuperCellId" cell data array where:
 * - Cells connected by fractures share the same super-cell ID
 * - Regular cells have their own unique super-cell ID (= their global ID)
 * 
 * This runs on rank 0 only, using the fracture neighbor information.
 * 
 * @param cells3D The 3D volumetric cells (modified in-place to add SuperCellId array)
 * @param fractureNeighbors Map of fracture name to neighbor mapping (fracture element → 3D cell neighbors)
 * @param comm MPI communicator
 * @return Super-cell metadata for partitioning
 */
SuperCellInfo tagCellsWithSuperCellIds(
  vtkSmartPointer< vtkUnstructuredGrid > cells3D,
  stdMap< string, ArrayOfArrays< vtkIdType, int64_t > > const & fractureNeighbors,
  MPI_Comm comm );

/**
 * @brief Reconstruct super-cell info from the SuperCellId array
 * 
 * After mesh redistribution, each rank needs to rebuild its local super-cell metadata
 * from the SuperCellId cell data array.
 * 
 * @param mesh The distributed mesh with SuperCellId array
 * @return Local super-cell metadata
 */
SuperCellInfo reconstructSuperCellInfo( vtkSmartPointer< vtkUnstructuredGrid > mesh );

/**
 * @brief Initial redistribution preserving super-cell integrity
 * @param cells3D Input mesh (only non-empty on rank 0)
 * @param comm MPI communicator
 * @return Redistributed mesh with SuperCellId array preserved
 * 
 * Uses simple round-robin assignment of super-cells to ranks.
 * This is faster than graph partitioning and suitable for initial distribution.
 */
vtkSmartPointer< vtkDataSet >
redistributeBySuperCellBlocks( vtkSmartPointer< vtkUnstructuredGrid > cells3D,
                               MPI_Comm comm );

/**
 * @brief Build a graph where nodes are super-cells (not individual cells)
 * 
 * Each super-cell becomes a single node in the graph. Edges connect super-cells
 * whose constituent cells are neighbors in the original mesh.
 * 
 * @param cells3D The tagged 3D mesh with SuperCellId array
 * @param baseGraph The original cell-to-cell adjacency graph
 * @param baseElemDist Element distribution for base graph (numRanks+1 array)
 * @param info Super-cell metadata
 * @param localStart Global index offset for this rank's cells (unused, for compatibility)
 * @param comm MPI communicator
 * @return Pair of (super-cell graph, super-cell vertex weights)
 */
std::pair< ArrayOfArrays< pmet_idx_t, pmet_idx_t >, array1d< pmet_idx_t > >
buildSuperCellGraph(
  vtkSmartPointer< vtkUnstructuredGrid > cells3D,
  ArrayOfArrays< pmet_idx_t, pmet_idx_t > const & baseGraph,
  arrayView1d< pmet_idx_t const > const & baseElemDist, 
  SuperCellInfo const & info,
  globalIndex const localStart,
  MPI_Comm comm );

/**
 * @brief Validate super-cell graph integrity before partitioning
 * 
 * Checks for:
 * - Self-loops
 * - Out-of-range neighbor indices
 * - Duplicate edges
 * - Isolated vertices
 * - Invalid vertex weights
 * 
 * @param superGraph The super-cell adjacency graph
 * @param superElemDist Element distribution array for super-cells
 * @param vertexWeights Vertex weights for load balancing
 * @param comm MPI communicator
 * @throws If graph has integrity errors
 */
void validateSuperCellGraph(
  ArrayOfArrays< pmet_idx_t, pmet_idx_t > const & superGraph,
  arrayView1d< pmet_idx_t const > const & superElemDist,
  arrayView1d< pmet_idx_t const > const & vertexWeights,
  MPI_Comm comm );

/**
 * @brief Unpack super-cell partitioning to individual cells
 * 
 * Maps the super-cell → rank assignments from ParMETIS back to individual cell assignments.
 * All cells in the same super-cell get assigned to the same rank.
 * 
 * @param cells3D The 3D mesh with SuperCellId array
 * @param superPartitioning Super-cell → rank assignments from ParMETIS
 * @param superCellIdToLocalIdx Mapping from SuperCellId to local super-cell index
 * @param comm MPI communicator
 * @return Cell-level partitioning array (cell → rank)
 */
array1d< int64_t >
unpackSuperCellPartitioning(
  vtkSmartPointer< vtkUnstructuredGrid > cells3D,
  array1d< int64_t > const & superPartitioning,
  stdMap< vtkIdType, localIndex > const & superCellIdToLocalIdx,
  MPI_Comm comm );

} // namespace vtk

} // namespace geos

#endif /* GEOS_MESH_GENERATORS_VTKSUPERCELLPARTITIONING_HPP */