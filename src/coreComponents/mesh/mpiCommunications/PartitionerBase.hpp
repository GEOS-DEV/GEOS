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
 * @file PartitionerBase.hpp
 * @brief Defines the thin abstract base class for all partitioner objects.
 */

#ifndef GEOS_PARTITIONER_PARTITIONERBASE_HPP_
#define GEOS_PARTITIONER_PARTITIONERBASE_HPP_

#include "dataRepository/Group.hpp"


namespace geos
{

/**
 * @class PartitionerBase
 * @brief A thin abstract base class that defines the universal interface for all partitioners.
 *
 * This class provides a common foundation for different partitioning strategies.
 *
 * Specific partitioning algorithms (e.g., graph-based vs. geometric) are defined in
 * specialized derived abstract classes like GraphPartitioner and GeometricPartitioner,
 * which introduce the relevant partitioning methods.
 */
class PartitionerBase : public dataRepository::Group
{
public:

  explicit PartitionerBase( string const & name,
                            Group * const parent );

  virtual ~PartitionerBase() = default;


  /// using alias for templated Catalog PartitionerBase type
  using CatalogInterface = dataRepository::CatalogInterface< PartitionerBase, string const &, Group * const >;

  /**
   * @brief Accessor for the singleton Catalog object used by the factory.
   * @return A static reference to the Catalog object.
   */
  static CatalogInterface::CatalogType & getCatalog();


  /**
   * @brief Processes command-line overrides for partition counts.
   *
   * Each partitioner must implement this to handle user-specified layouts
   * (e.g., -x, -y, -z flags).
   */
  virtual void processCommandLineOverrides( unsigned int xparCL,
                                            unsigned int yparCL,
                                            unsigned int zparCL ) = 0;

  /**
   * @brief Sets the total number of partitions and their logical layout.
   *
   * The interpretation of x, y, and z partitions is implementation-dependent.
   */
  virtual void setPartitionCounts( unsigned int xPartitions,
                                   unsigned int yPartitions,
                                   unsigned int zPartitions ) = 0;

  /**
   * @brief Sets the ranks of neighboring partitions.
   *
   * This is typically used to build an adjacency graph for coloring.
   */
  virtual void setNeighborsRank( const std::vector< int > & neighborsRank ) = 0;

  /**
   * @brief Get the lightweight neighbor rank list (first-order neighbors only).
   * @return Vector of neighbor rank IDs.
   */
  std::vector< int > const & getNeighborsRank() const
  {
    return m_neighborsRank;
  }

  /**
   * @brief Get number of neighbors.
   */
  int getNumNeighbors() const
  {
    return static_cast< int >(m_neighborsRank.size());
  }

  /**
   * @brief Assigns a color to the current rank based on neighbor adjacencies.
   *
   * Ensures that no two adjacent partitions have the same color. This is essential
   * for scheduling non-conflicting communication and computation phases.
   */
  virtual void color() = 0;

  /**
   * @brief Gets the color assigned to this rank's partition.
   */
  int getColor() const
  {
    return m_color;
  }

  /**
   * @brief Gets the total number of colors used across all partitions.
   */
  int getNumColors() const
  {
    return m_numColors;
  }

  /**
   * @brief Gets the total number of partitions.
   */
  int getNumPartitions() const
  {
    return m_numPartitions;
  }

protected:
  /// List of neighbor rank IDs (first-order neighbors only).
  std::vector< int > m_neighborsRank;

  /// Total number of partitions in the simulation.
  int m_numPartitions;

  /// Total number of colors.
  int m_numColors;

  /// The color assigned to this specific rank/partition.
  int m_color;
};

} // namespace geos

#endif // GEOS_PARTITIONER_PARTITIONERBASE_HPP_
