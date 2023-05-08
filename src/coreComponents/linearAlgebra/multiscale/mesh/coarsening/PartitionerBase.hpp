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
 * @file PartitionerBase.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_MULTISCALE_PARTITIONERBASE_HPP
#define GEOSX_LINEARALGEBRA_MULTISCALE_PARTITIONERBASE_HPP

#include "linearAlgebra/multiscale/mesh/MeshLevel.hpp"
#include "linearAlgebra/utilities/LinearSolverParameters.hpp"

#include <memory>

namespace geos
{
namespace multiscale
{

/**
 * @brief Base class for partitioner implementations.
 */
class PartitionerBase
{
public:

  /**
   * @brief Factory method for instantiating a partitioner based on parameters.
   * @param params coarsening parameters
   * @return owning pointer to the new partitioner object
   */
  static std::unique_ptr< PartitionerBase >
  create( LinearSolverParameters::Multiscale::Coarsening params );

  /**
   * @brief Constructor.
   * @param params coarsening parameters
   */
  explicit PartitionerBase( LinearSolverParameters::Multiscale::Coarsening params )
    : m_params( std::move( params ) )
  {}

  /**
   * @brief Destructor.
   */
  virtual ~PartitionerBase() = default;

  /**
   * @brief Generate a partitioning of fine-scale mesh cells.
   * @param mesh the fine-scale mesh
   * @param partition the partition index output array (that must be properly sized)
   * @return the number of partitions generated
   */
  virtual localIndex generate( multiscale::MeshLevel const & mesh,
                               arrayView1d< localIndex > const & partition ) = 0;

  /**
   * @brief Store auxiliary partitioning-related data on the coarse mesh.
   * @param coarseMesh the coarse mesh
   *
   * This function can be used to transfer auxiliary data used by the partitioner implementation
   * onto the coarse grid after it's been created based on the previously generated partition.
   * For example, a Cartesian partitioner may need to assign Cartesian indices to newly generated
   * coarse cells, so that it can be later applied to the coarse grid recursively.
   *
   * @note It is expected that the coarse mesh contains exactly as many cells as there are unique
   *       partition indices produced by generate().
   */
  virtual void setCoarseData( multiscale::MeshLevel & coarseMesh ) const
  {
    GEOS_UNUSED_VAR( coarseMesh );
  };

protected:

  /// Coarsening parameters
  LinearSolverParameters::Multiscale::Coarsening m_params;
};

} // namespace multiscale
} // namespace geos

#endif //GEOSX_LINEARALGEBRA_MULTISCALE_PARTITIONERBASE_HPP
