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

/**
 * @file GeometricPartitioner.hpp
 * @brief Base class for geometric (coordinate-based) domain partitioners
 */

#ifndef GEOS_MESH_MPICOMMUNICATIONS_GEOMETRICPARTITIONER_HPP_
#define GEOS_MESH_MPICOMMUNICATIONS_GEOMETRICPARTITIONER_HPP_

#include "DomainPartitioner.hpp"

namespace geos
{

/**
 * @class GeometricPartitioner
 * @brief Base class for partitioners that use geometric decomposition
 *
 * Geometric partitioners decompose the domain based on spatial coordinates rather than
 * mesh connectivity:
 * - Don't require loaded meshes
 * - Don't use graph partitioning algorithms
 * - Compute neighbors from geometric topology
 * - Can be used before mesh generation
 * - Provide spatial queries (isCoordInPartition)
 *
 * Concrete implementations:
 * - CartesianPartitioner: Regular X × Y × Z grid decomposition
 *
 * Workflow:
 * 1. Create partitioner from XML
 * 2. Call initializeDomain() with global bounding box
 * 3. Partitioner computes local subdomain and neighbors
 * 4. Use isCoordInPartition() to query point ownership
 *
 */
class GeometricPartitioner : public DomainPartitioner
{
public:

  /**
   * @brief Constructor
   * @param name The name of this partitioner instance
   * @param parent The parent group
   */
  explicit GeometricPartitioner( string const & name,
                                 dataRepository::Group * const parent );

  /**
   * @brief Destructor
   */
  virtual ~GeometricPartitioner() override;

  /**
   * @brief Initialize the partitioner with the global domain's bounding box
   *
   * This method:
   * 1. Computes local subdomain bounds for this rank
   * 2. Determines neighbor ranks from geometric topology
   * 3. Computes communication coloring
   *
   * @param globalMin The minimum coordinates (x, y, z) of the global domain
   * @param globalMax The maximum coordinates (x, y, z) of the global domain
   *
   * @post m_neighborsRank is set
   * @post m_color and m_numColors are set
   * @post Local subdomain bounds are computed
   *
   */
  virtual void initializeDomain( R1Tensor const & globalMin,
                                 R1Tensor const & globalMax ) = 0;

  /**
   * @brief Check if a coordinate falls within this rank's local partition
   *
   * @param coords The spatial coordinates (x, y, z) to test
   * @return true if the coordinate is within the local partition's boundaries
   *
   * @pre initializeDomain() must have been called
   */
  virtual bool isCoordInPartition( R1Tensor const & coords ) const = 0;

protected:

  /**
   * @brief Compute neighbors from geometric topology
   *
   * @post m_neighborsRank is populated
   *
   */
  virtual void computeNeighborsFromTopology() = 0;
};

} // namespace geos

#endif // GEOS_MESH_MPICOMMUNICATIONS_GEOMETRICPARTITIONER_HPP_
