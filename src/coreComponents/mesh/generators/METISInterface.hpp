/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file METISInterface.hpp
 */

#ifndef GEOS_MESH_GENERATORS_METISINTERFACE_HPP_
#define GEOS_MESH_GENERATORS_METISINTERFACE_HPP_

#include "common/DataTypes.hpp"
#include "mesh/generators/ParMETISInterface.hpp"

namespace geos
{
namespace metis
{

/**
 * @brief Partition a root-local weighted graph with serial METIS.
 * @param graph Symmetric CSR adjacency graph.
 * @param edgeWeights Positive edge weights in CSR order.
 * @param vertexWeights Vertex weights with shape [numVertices][numConstraints].
 * @param numParts Number of target partitions.
 * @param imbalance Relative imbalance tolerance for every constraint.
 * @param seed Deterministic METIS seed.
 * @return One target partition per graph vertex.
 */
array1d< pmet_idx_t >
partitionWeighted( ArrayOfArraysView< pmet_idx_t const, pmet_idx_t > const & graph,
                   arrayView1d< pmet_idx_t const > const & edgeWeights,
                   arrayView2d< pmet_idx_t const > const & vertexWeights,
                   pmet_idx_t numParts,
                   arrayView1d< real64 const > const & imbalance,
                   pmet_idx_t seed );

} // namespace metis
} // namespace geos

#endif // GEOS_MESH_GENERATORS_METISINTERFACE_HPP_
