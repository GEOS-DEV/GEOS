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

#include "common/DataTypes.hpp"
#include "common/MpiWrapper.hpp"
#include "common/TimingMacros.hpp"
#include "mesh/generators/VTKMeshGeneratorTools.hpp"
#include "mesh/generators/ParMETISInterface.hpp"
#include "LvArray/src/ArrayOfArrays.hpp"

#include <vtkSmartPointer.h>

class vtkUnstructuredGrid;
class vtkDataSet;

namespace geos
{
namespace vtk
{

/**
 * @brief Initial distribution strategy for super-cells
 */
enum class InitialDistributionStrategy
{
  MORTON,  ///< Morton Z-curve ordering for spatial locality
  BLOCK    ///< Contiguous block distribution
};

/**
 * @brief Super-cell metadata for constrained partitioning
 *
 * Groups cells that must remain co-located (e.g., fracture-connected cells).
 */
struct SuperCellInfo
{
  /// SuperCellId to global cell IDs
  std::map< vtkIdType, std::vector< vtkIdType > > superCellToOriginalCells;

  /// SuperCellId to vertex weight
  std::map< vtkIdType, vtkIdType > vertexWeights;

  /// SuperCellIds containing multiple cells
  std::set< vtkIdType > atomicSuperCells;
};

/**
 * @brief Assign super-cell IDs based on fracture connectivity
 *
 * Cells sharing a fracture get the same ID; isolated cells use their global ID.
 * Adds a "SuperCellId" array to the mesh. Runs on rank 0.
 *
 * @param cells3D Volumetric mesh (modified in-place)
 * @param fractureNeighbors Fracture element to 3D neighbor mappings
 * @param fractureWeight Load balancing weight boost for fracture super-cells
 * @return Super-cell metadata
 */
SuperCellInfo tagCellsWithSuperCellIds(
  vtkSmartPointer< vtkUnstructuredGrid > cells3D,
  stdMap< string, ArrayOfArrays< localIndex, int64_t > > const & fractureNeighbors,
  integer fractureWeight );

/**
 * @brief Rebuild super-cell metadata from SuperCellId array
 *
 * After redistribution, each rank reconstructs local super-cell info.
 * Applies fracture weights to account for contact computation cost.
 *
 * @param mesh Distributed mesh with SuperCellId array
 * @param fractureWeight Weight boost for multi-cell super-cells
 * @return Local super-cell metadata
 */
SuperCellInfo reconstructSuperCellInfo( vtkSmartPointer< vtkUnstructuredGrid > mesh,
                                        integer fractureWeight );

/**
 * @brief Initial redistribution preserving super-cell integrity
 *
 * Distributes super-cells across ranks without graph partitioning.
 * Faster than ParMETIS for initial scatter.
 *
 * @param cells3D Input mesh (non-empty on rank 0 only)
 * @param comm MPI communicator
 * @param strategy Distribution strategy (MORTON or BLOCK)
 * @return Redistributed mesh with preserved SuperCellId array
 */
vtkSmartPointer< vtkDataSet >
redistributeBySuperCellBlocks( vtkSmartPointer< vtkUnstructuredGrid > cells3D,
                               MPI_Comm comm,
                               InitialDistributionStrategy strategy = InitialDistributionStrategy::MORTON );

/**
 * @brief Build super-cell adjacency graph
 *
 * Collapses the cell graph into a super-cell graph where each vertex
 * represents a super-cell and edges connect neighboring super-cells.
 *
 * @param cells3D Mesh with SuperCellId array
 * @param baseGraph Original cell-to-cell adjacency
 * @param baseElemDist Element distribution for base graph
 * @param info Super-cell metadata
 * @param comm MPI communicator
 * @return (super-cell graph, vertex weights)
 */
std::pair< ArrayOfArrays< pmet_idx_t, pmet_idx_t >, array1d< pmet_idx_t > >
buildSuperCellGraph(
  vtkSmartPointer< vtkUnstructuredGrid > cells3D,
  ArrayOfArrays< pmet_idx_t, pmet_idx_t > const & baseGraph,
  arrayView1d< pmet_idx_t const > const & baseElemDist,
  SuperCellInfo const & info,
  MPI_Comm comm );

/**
 * @brief Validate super-cell graph before partitioning
 *
 * Checks for self-loops, invalid indices, duplicate edges,
 * isolated vertices, and invalid weights.
 *
 * @param superGraph Super-cell adjacency graph
 * @param superElemDist Element distribution array
 * @param vertexWeights Vertex weights
 * @param comm MPI communicator
 */
void validateSuperCellGraph(
  ArrayOfArrays< pmet_idx_t, pmet_idx_t > const & superGraph,
  arrayView1d< pmet_idx_t const > const & superElemDist,
  arrayView1d< pmet_idx_t const > const & vertexWeights,
  MPI_Comm comm );

/**
 * @brief Map super-cell partitions back to individual cells
 *
 * Converts super-cell -> rank assignments into cell -> rank assignments.
 * All cells in a super-cell go to the same rank.
 *
 * @param cells3D Mesh with SuperCellId array
 * @param superPartitioning Super-cell to rank assignments from ParMETIS
 * @param superCellIdToLocalIdx SuperCellId to local index mapping
 * @param comm MPI communicator
 * @return Cell-level partition assignments
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
