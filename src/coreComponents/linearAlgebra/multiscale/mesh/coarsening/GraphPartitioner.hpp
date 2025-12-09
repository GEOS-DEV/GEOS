/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file GraphPartitioner.hpp
 */

#ifndef GEOS_LINEARALGEBRA_MULTISCALE_GRAPHPARTITIONER_HPP_
#define GEOS_LINEARALGEBRA_MULTISCALE_GRAPHPARTITIONER_HPP_

#include "PartitionerBase.hpp"

namespace geos
{
namespace multiscale
{

/**
 * @brief Graph-based mesh partitioner (uses METIS or SCOTCH).
 */
class GraphPartitioner final : public PartitionerBase
{
public:

  using PartitionerBase::PartitionerBase;

  virtual localIndex generate( multiscale::MeshLevel const & mesh,
                               arrayView1d< localIndex > const & partition ) override;
};

} // namespace multiscale
} // namespace geos

#endif //GEOS_LINEARALGEBRA_MULTISCALE_GRAPHPARTITIONER_HPP_
