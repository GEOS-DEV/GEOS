/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file ScotchInterface.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_MULTISCALE_SCOTCHINTERFACE_HPP_
#define GEOSX_LINEARALGEBRA_MULTISCALE_SCOTCHINTERFACE_HPP_

#include "common/DataTypes.hpp"
#include "linearAlgebra/utilities/LinearSolverParameters.hpp"

namespace geos
{
namespace scotch
{

/**
 * @brief Compute a graph partitioning using Scotch.
 * @param graph the input graph
 * @param params Scotch parameters
 * @param numPart number of partitions
 * @param partition writable view of pre-allocated output array
 */
void partition( CRSMatrixView< int64_t const, int64_t const, int64_t const > const & graph,
                LinearSolverParameters::Multiscale::Coarsening::Graph::Scotch const & params,
                int64_t const numPart,
                arrayView1d< int64_t > const & partition );

} // namespace scotch
} // namespace geos

#endif //GEOSX_LINEARALGEBRA_MULTISCALE_SCOTCHINTERFACE_HPP_
