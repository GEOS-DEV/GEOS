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
 * @file VTKMeshGeneratorTools.hpp
 */

#ifndef GEOSX_VTK_MESH_GENERATORS_VTKMESHGENERATORTOOLS_HPP_
#define GEOSX_VTK_MESH_GENERATORS_VTKMESHGENERATORTOOLS_HPP_

#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPartitionedDataSet.h>

#include <mpi.h>

// NOTE: do NOT include anything from GEOSX here.
// In particular, nothing that directly or transitively includes "common/Format.hpp".
// The reason is "diy2" library includes an older version of {fmt} than the one used by GEOSX.
// Collision of includes leads to all kinds of impossible to fix compilation errors.
// Thankfully, no link errors, owing to namespace versioning employed by {fmt}.

namespace geosx
{
namespace vtk
{

/**
 * @brief Redistribute a dataset partitioned on each rank according to destination.
 * @param localParts the partitioned dataset (each part must be a vtkUnstructuredGrid)
 * @param mpiComm the MPI communicator to distribute over
 * @return the new, redistributed grid
 *
 * There must be exactly P partitions in @p localParts, where P is the size of @p comm.
 * Partition with index i represents a piece of mesh that must be shipped off to rank i.
 * Some partitions (usually most of them) can be empty, indicating nothing to send.
 * The return value on each rank is a combination of mesh pieces sent to current rank.
 */
vtkSmartPointer< vtkUnstructuredGrid >
redistribute( vtkPartitionedDataSet & localParts, MPI_Comm mpiComm );

} // namespace vtk
} // namespace geosx

#endif //GEOSX_VTK_MESH_GENERATORS_VTKMESHGENERATORTOOLS_HPP_
