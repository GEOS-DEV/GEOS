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
 * @file CartesianPartitioner.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_MULTISCALE_CARTESIANPARTITIONER_HPP_
#define GEOSX_LINEARALGEBRA_MULTISCALE_CARTESIANPARTITIONER_HPP_

#include "linearAlgebra/multiscale/mesh/coarsening/PartitionerBase.hpp"

namespace geos
{
namespace multiscale
{

/**
 * @brief Simple structured partitioner for Cartesian grids.
 */
class CartesianPartitioner final : public PartitionerBase
{
public:

  using PartitionerBase::PartitionerBase;

  virtual localIndex generate( multiscale::MeshLevel const & mesh,
                               arrayView1d< localIndex > const & partition ) override;

  virtual void setCoarseData( multiscale::MeshLevel & coarseMesh ) const override;

private:

  integer m_numPart[3]{};

};

} // namespace multiscale
} // namespace geos

#endif //GEOSX_LINEARALGEBRA_MULTISCALE_CARTESIANPARTITIONER_HPP_
