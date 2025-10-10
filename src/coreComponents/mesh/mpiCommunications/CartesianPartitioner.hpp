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
 * @file CartesianPartitioner.hpp
 */

#ifndef GEOS_MESH_MPICOMMUNICATIONS_CARTESIANPARTITIONER_HPP_
#define GEOS_MESH_MPICOMMUNICATIONS_CARTESIANPARTITIONER_HPP_

#include "GeometricPartitioner.hpp"

namespace geos
{

/**
 * @class CartesianPartitioner
 * @brief Cartesian (regular grid) domain decomposition
 *
 * CartesianPartitioner divides a rectangular domain into a regular X × Y × Z
 * grid of subdomains. Each MPI rank is assigned one subdomain.
 *
 * Features:
 * - Checkerboard coloring (8 colors maximum)
 * - Support for periodic boundaries in any direction
 * - 26-connected neighbor topology (includes faces, edges, and corners)
 *
 */
class CartesianPartitioner : public GeometricPartitioner
{
public:

  /**
   * @brief Constructor
   * @param name The name of this partitioner instance
   * @param parent The parent group
   */
  explicit CartesianPartitioner( string const & name,
                                 dataRepository::Group * const parent );

  /**
   * @brief Destructor
   */
  virtual ~CartesianPartitioner() override = default;

  /**
   * @brief Catalog name for factory registration
   * @return The catalog key
   */
  static string catalogName() { return "CartesianPartitioner"; }

  /**
   * @brief Initialize the Cartesian grid partitioner
   *
   * Implementation of GeometricPartitioner::initializeDomain()
   *
   * Steps:
   * 1. Create MPI Cartesian communicator
   * 2. Compute local subdomain bounds
   * 3. Determine neighbor ranks from grid topology
   * 4. Compute checkerboard coloring
   *
   * @param globalMin Global domain minimum (x, y, z)
   * @param globalMax Global domain maximum (x, y, z)
   *
   * @post m_neighborsRank is populated
   * @post m_color and m_numColors are set
   * @post Local bounds (m_localMin, m_localMax) are set
   */
  void initializeDomain( R1Tensor const & globalMin,
                         R1Tensor const & globalMax ) override;

  /**
   * @brief Check if coordinate is in local Cartesian subdomain
   *
   * Handles periodic boundaries by mapping coordinates into
   * the periodic range before checking bounds.
   *
   * @param coords Spatial coordinates (x, y, z) to test
   * @return true if coordinate is in local partition
   *
   * @pre initializeDomain() must have been called
   */
  bool isCoordInPartition( R1Tensor const & coords ) const override;

  /**
   * @brief Check if coordinate is in partition along one axis
   *
   * Helper for isCoordInLocalPartition().
   * Handles periodic boundaries.
   *
   * @param coord Coordinate value to test
   * @param dir Direction (0=x, 1=y, 2=z)
   * @return true if coordinate is in range for this direction
   */
  bool isCoordInPartition( real64 const & coord, int const dir ) const;


  /**
   * @brief Process command-line partition overrides (-x, -y, -z flags)
   *
   * Priority:
   * 1. Command-line flags (if provided)
   * 2. XML partitionCounts (if specified)
   * 3. Default: 1 × 1 × N (1D decomposition along z-axis)
   *
   * @param xparCL X-direction partition count (0 = not specified)
   * @param yparCL Y-direction partition count (0 = not specified)
   * @param zparCL Z-direction partition count (0 = not specified)
   *
   * @post m_partitionCounts is set
   * @post m_numPartitions is set
   */
  void processCommandLineOverrides( unsigned int const xparCL,
                                    unsigned int const yparCL,
                                    unsigned int const zparCL ) override;

  /**
   * @brief Set partition counts in each direction
   *
   * @param xPartitions Number of partitions in X direction
   * @param yPartitions Number of partitions in Y direction
   * @param zPartitions Number of partitions in Z direction
   *
   * @post m_partitionCounts = {xPartitions, yPartitions, zPartitions}
   * @post m_numPartitions = xPartitions * yPartitions * zPartitions
   */
  void setPartitionCounts( unsigned int const xPartitions,
                           unsigned int const yPartitions,
                           unsigned int const zPartitions ) override;

  /**
   * @brief Get partition counts in all directions
   * @return Array of partition counts {nx, ny, nz}
   */
  array1d< int > const & getPartitionCounts() const { return m_partitionCounts; }

  /**
   * @brief Set periodicity for specific directions
   *
   * Used by mesh generators to enable
   * periodic boundaries in specific directions.
   *
   * @param dim Direction (0=x, 1=y, 2=z)
   * @param periodic Periodicity flag (1=periodic, 0=not periodic)
   */
  void setPeriodicity( int const dim, int const periodic )
  {
    GEOS_ERROR_IF( dim < 0 || dim >= m_ndim,
                   "Invalid dimension " << dim << ". Must be 0, 1, or 2." );
    m_periodic[dim] = periodic;
  }

  /**
   * @brief Get periodicity flags
   * @return Array of periodicity flags {px, py, pz} where 1 = periodic, 0 = not periodic
   */
  array1d< int > const & getPeriodicity() const { return m_periodic; }

  /**
   * @brief Compute checkerboard coloring
   *
   * @post m_color is set (0-7)
   * @post m_numColors is set (≤ 8)
   */
  void color() override;

  /**
   * @brief Get X-direction partition count
   * @return Number of partitions in X direction
   */
  unsigned int getXPartitions() const { return m_partitionCounts[0]; }

  /**
   * @brief Get Y-direction partition count
   * @return Number of partitions in Y direction
   */
  unsigned int getYPartitions() const { return m_partitionCounts[1]; }

  /**
   * @brief Get Z-direction partition count
   * @return Number of partitions in Z direction
   */
  unsigned int getZPartitions() const { return m_partitionCounts[2]; }

  real64 const * getLocalMin() const { return m_localMin;}

  real64 const * getLocalMax() const { return m_localMax; }

  real64 const * getGlobalMin() const { return m_globalGridMin; }

  real64 const * getGlobalMax() const { return m_globalGridMax; }


  /**
   * @brief Struct to serve as a container for variable strings and keys.
   * @struct viewKeyStruct
   */
  struct viewKeyStruct : GeometricPartitioner::viewKeyStruct
  {
    /// String key for partition counts
    static constexpr char const * partitionCountsString() { return "partitionCounts"; }
  };

protected:

  /**
   * @brief Compute neighbors from Cartesian grid topology
   *
   * Finds all neighbors in a 3×3×3 stencil (26-connected):
   * - 6 face neighbors
   * - 12 edge neighbors
   * - 8 corner neighbors
   *
   * Handles periodic boundaries.
   *
   * @post m_neighborsRank is populated
   */
  void computeNeighborsFromTopology() override;


private:

  /**
   * @brief Initialize MPI Cartesian communicator
   *
   * Creates an MPI_Cart communicator with the specified topology.
   * Sets m_cartRank and m_coords based on rank position in grid.
   *
   * @param[out] cartComm The created Cartesian communicator
   *
   * @pre m_partitionCounts must be set
   * @post m_cartRank is set
   * @post m_coords is set
   */
  void initializeCartesianCommunicator( MPI_Comm & cartComm );

  /**
   * @brief Recursively add all neighbors in 3×3×3 stencil
   *
   * @param idim Current dimension being explored (0, 1, or 2)
   * @param ncoords Neighbor coordinates being constructed
   * @param cartComm Cartesian communicator
   */
  void addNeighbors( unsigned int const idim, int * ncoords, MPI_Comm const cartComm );

  /**
   * @brief Set global domain values
   *
   * Stores global bounding box and computes global domain size.
   *
   * @param gridMin Global domain minimum
   * @param gridMax Global domain maximum
   *
   * @post m_globalGridMin, m_globalGridMax, m_globalGridSize are set
   */
  void setGlobalGridValues( R1Tensor const & gridMin, R1Tensor const & gridMax );

  /**
   * @brief Compute local subdomain bounds
   *
   * Divides global domain equally among partitions and computes
   * this rank's portion.
   *
   * @param gridMin Global domain minimum
   *
   * @post m_localMin, m_localMax, m_localSize are set
   */
  void setLocalPartitionValues( R1Tensor const & gridMin );

  /**
   * @brief Validate that partition count matches MPI size
   *
   * @param size MPI communicator size
   * @throw std::runtime_error if size != m_numPartitions
   */
  void validatePartitionSize( int const size ) const;

protected:
  /// Number of spatial dimensions (always 3)
  static constexpr int m_ndim = 3;

  /// Periodic boundary flags for each direction {x, y, z}
  array1d< int > m_periodic;

  /// Cartesian grid coordinates of this rank {i, j, k}
  array1d< int > m_coords;

  /// Number of partitions in each direction {nx, ny, nz}
  array1d< int > m_partitionCounts;

  /// Local subdomain minimum (x, y, z)
  real64 m_localMin[3];

  /// Local subdomain maximum (x, y, z)
  real64 m_localMax[3];

  /// Local subdomain size (x, y, z)
  real64 m_localSize[3];

  /// Global domain minimum (x, y, z)
  real64 m_globalGridMin[3];

  /// Global domain maximum (x, y, z)
  real64 m_globalGridMax[3];

  /// Global domain size (x, y, z)
  real64 m_globalGridSize[3];

  /// Rank in Cartesian communicator
  int m_cartRank;
};

} // namespace geos

#endif // GEOS_MESH_MPICOMMUNICATIONS_CARTESIANPARTITIONER_HPP_
