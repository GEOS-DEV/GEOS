/*
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

#ifndef GEOS_MESH_GENERATORS_VTKMESHSCATTERING_HPP
#define GEOS_MESH_GENERATORS_VTKMESHSCATTERING_HPP

#include "common/DataTypes.hpp"
#include "common/format/EnumStrings.hpp"
#include "common/MpiWrapper.hpp"

#include <vtkDataSet.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

namespace geos
{
namespace vtk
{

/**
 * @brief Method used to partition and scatter the mesh across MPI ranks.
 */
enum class ScatterMethod
{
  contiguous, ///< Assign contiguous blocks of cell indices to each rank (no geometry awareness).
  cartesian,  ///< Partition based on a user-specified Cartesian grid (nx × ny × nz).
  rcb,        ///< Recursive Coordinate Bisection along the longest bounding-box axis.
  kdtree      ///< VTK's built-in kd-tree redistribution (vtkRedistributeDataSetFilter). Retained for backward compatibility.
};

ENUM_STRINGS( ScatterMethod,
              "contiguous",
              "cartesian",
              "rcb",
              "kdtree" );

/**
 * @brief Compute the per-cell destination rank on rank 0 for a given scatter method.
 *
 * This is the assignment step of @ref scatterMesh, exposed so other code (e.g. the
 * super-cell aware scatter) can build its own assignment on top of the same algorithms.
 *
 * @param[in] method                the partitioning strategy. @p kdtree is not supported here
 *                                  (it does not produce a separable per-cell assignment).
 * @param[in] mesh                  the input mesh (must contain all cells on rank 0).
 * @param[in] cartesianPartitions   {nx, ny, nz}, only used for @p cartesian.
 * @param[in] numRanks              total number of MPI ranks.
 * @return rank-0 only: @c cellRanks[i] is the destination rank for cell @p i; empty on other ranks.
 */
stdVector< integer >
computeCellRanks( ScatterMethod method,
                  vtkDataSet & mesh,
                  arrayView1d< integer const > cartesianPartitions,
                  integer numRanks );

/**
 * @brief Ship cells from rank 0 to the ranks specified by @p assignment via a binary-tree exchange.
 *
 * This is the shipping step of @ref scatterMesh, exposed so other code can drive it with a
 * custom @p assignment (e.g. atom-aware partitioning that keeps fracture-connected cells together).
 *
 * @param[in] inputMesh   rank 0: source mesh with all cells. Other ranks: ignored (may be empty).
 * @param[in] assignment  rank 0 only: @c assignment[i] is the destination rank for cell @p i.
 * @param[in] comm        the MPI communicator.
 * @return the local subset of the mesh assigned to this rank.
 */
vtkSmartPointer< vtkUnstructuredGrid >
scatterByRankAssignment( vtkUnstructuredGrid * inputMesh,
                         stdVector< integer > assignment,
                         MPI_Comm comm );

/**
 * @brief Scatter a mesh held entirely on rank 0 to all MPI ranks.
 *
 * The mesh is first partitioned on rank 0 using the selected method, producing
 * a cell to rank assignment vector. The data is then distributed via a binary-tree
 * scatter pattern.
 *
 * @param[in] method     the partitioning strategy to use
 * @param[in] mesh       the input mesh (must contain all cells on rank 0; empty on other ranks)
 * @param[in] cartesianPartitions additional parameters for the cartesian partitioning method:
 *                        - For @p cartesian: must contain at least 3 values (nx, ny, nz)
 *                          with nx*ny*nz == MPI size.
 *                        - Ignored for other methods.
 * @param[in] comm       the MPI communicator
 * @return the local partition of the mesh on each rank
 */
vtkSmartPointer< vtkDataSet >
scatterMesh( ScatterMethod method,
             vtkDataSet & mesh,
             arrayView1d< integer const > cartesianPartitions,
             MPI_Comm comm );

} // namespace vtk
} // namespace geos

#endif // GEOS_MESH_GENERATORS_VTKMESHSCATTERING_HPP
