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
 * @file ParMETISInterface.hpp
 */

#ifndef GEOSX_MESH_GENERATORS_PARMETISINTERFACE_HPP_
#define GEOSX_MESH_GENERATORS_PARMETISINTERFACE_HPP_

#include "common/DataTypes.hpp"
#include "common/MpiWrapper.hpp"

namespace geosx
{

#if defined(GEOSX_USE_HIP) // still need int32 hypre for the current hip-capable build
/// Typedef to allow us to specify required parmetis integer type.
using pmet_idx_t = int32_t;
#else
/// Typedef to allow us to specify required parmetis integer type.
using pmet_idx_t = int64_t;
#endif

namespace parmetis
{

/**
 * @brief Partition a mesh according to its dual graph.
 * @param elemToNodes the input mesh represented by its elem-node graph
 * @param comm the MPI communicator of processes to partition over
 * @param minCommonNodes minimum number of shared nodes to create an graph edge
 * @param numRefinements number of partition refinement iterations
 * @return an array of target partitions for each element in local mesh
 * @note the number of partitions will be equal to the size of @p comm.
 */
array1d< pmet_idx_t >
partMeshKway( ArrayOfArraysView< pmet_idx_t const, pmet_idx_t > const & elemToNodes,
              MPI_Comm comm,
              int minCommonNodes,
              int numRefinements );

} // namespace parmetis
} // namespace geosx

#endif //GEOSX_MESH_GENERATORS_PARMETISINTERFACE_HPP_
