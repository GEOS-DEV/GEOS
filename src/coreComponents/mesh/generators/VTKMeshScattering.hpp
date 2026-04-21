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
