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
 * @file DomainPartitioner.hpp
 */

#ifndef GEOS_MESH_MPICOMMUNICATIONS_DOMAINPARTITIONER_HPP_
#define GEOS_MESH_MPICOMMUNICATIONS_DOMAINPARTITIONER_HPP_

#include "dataRepository/Group.hpp"
#include "common/DataTypes.hpp"

namespace geos
{

/**
 * @class DomainPartitioner
 * @brief Abstract base class for all domain decomposition strategies
 *
 * ALL domain partitioners must:
 * - Compute MPI neighbor ranks
 * - Assign communication colors
 * - Support partition count configuration
 *
 * Concrete implementations:
 * - MeshPartitioner: Partitioners that operate on loaded meshes (uses GraphPartitionEngine)
 *   - CellGraphPartitioner: Standard cell-based dual graph partitioning
 *   - LayeredGraphPartitioner: Structured/layered mesh partitioning
 * - GeometricPartitioner: Partitioners based on spatial coordinates
 *   - CartesianPartitioner: Cartesian domain decomposition
 *   - ParticleCartesianPartitioner: Cartesian for particle meshes
 */
class DomainPartitioner : public dataRepository::Group
{
public:

  /**
   * @brief Constructor
   * @param name The name of this partitioner instance
   * @param parent The parent group
   */
  explicit DomainPartitioner( string const & name,
                              dataRepository::Group * const parent );

  /**
   * @brief Destructor
   */
  virtual ~DomainPartitioner() override;


  /**
   * @brief Get the catalog name for this partitioner type
   *
   */
  string getCatalogName() const
  {
    // Extract type name from the data context
    return getDataContext().toString();
  }

  /**
   * @brief Post-input initialization
   * @note Made public to allow initialization of dynamically created default partitioners
   */
  virtual void postInputInitialization() override;

  /**
   * @brief Process command-line partition count overrides (-x, -y, -z flags)
   *
   * This is called after XML parsing to allow command-line overrides.
   *
   * @param xparCL X-direction partition count from command line (0 = no override)
   * @param yparCL Y-direction partition count from command line (0 = no override)
   * @param zparCL Z-direction partition count from command line (0 = no override)
   */
  virtual void processCommandLineOverrides( unsigned int const xparCL,
                                            unsigned int const yparCL,
                                            unsigned int const zparCL ) = 0;

  /**
   * @brief Set target partition counts in each logical dimension
   *
   * @param xPartitions Number of partitions in X direction
   * @param yPartitions Number of partitions in Y direction
   * @param zPartitions Number of partitions in Z direction
   *
   * @post m_numPartitions = xPartitions * yPartitions * zPartitions
   */
  virtual void setPartitionCounts( unsigned int const xPartitions,
                                   unsigned int const yPartitions,
                                   unsigned int const zPartitions );

  /**
   * @brief Set the list of neighboring MPI ranks
   *
   * This is typically called internally by the partitioner after domain decomposition.
   *
   * @param neighborsRank Vector of neighbor rank IDs
   */
  void setNeighborsRank( std::vector< int > const & neighborsRank );

  /**
   * @brief Get the list of neighboring MPI ranks
   *
   * This is called by DomainPartition::setupBaseLevelMeshGlobalInfo() to configure
   * ghost communication.
   *
   * @return Vector of neighbor rank IDs
   *
   * @pre Partitioner must be fully initialized (neighbors computed)
   */
  std::vector< int > const & getNeighborsRank() const { return m_neighborsRank; }

  /**
   * @brief Get the number of neighboring ranks
   * @return Number of neighbors
   */
  int getNumNeighbors() const { return static_cast< int >( m_neighborsRank.size() ); }

  /**
   * @brief Compute graph coloring for communication scheduling
   *
   * Colors are assigned such that no two neighboring ranks have the same color.
   * This enables non-blocking communication patterns.
   *
   * @post m_color and m_numColors are set
   *
   * @note Must be called after neighbors are determined
   */
  virtual void color() = 0;

  /**
   * @brief Get the color assigned to this rank
   * @return Color ID [0, numColors)
   */
  int getColor() const { return m_color; }

  /**
   * @brief Get the total number of colors used
   * @return Number of colors
   */
  int getNumColors() const { return m_numColors; }

  /**
   * @brief Get the total number of partitions
   * @return Total partition count
   */
  int getNumPartitions() const { return m_numPartitions; }


  /**
   * @brief Print information about this partitioner to log
   */
  virtual void printInfo() const;

  /**
   * @brief Get a descriptive string about this partitioner
   */
  virtual string getInfoString() const;

  /**
   * @brief Accessor for the singleton Catalog object
   *
   * @return Reference to the Catalog object
   */
  using CatalogInterface = dataRepository::CatalogInterface< DomainPartitioner,
                                                             string const &,
                                                             dataRepository::Group * const >;

  static CatalogInterface::CatalogType & getCatalog();

protected:

  /// MPI neighbor ranks (set after domain decomposition)
  std::vector< int > m_neighborsRank;

  /// Total number of partitions
  int m_numPartitions;

  /// Total number of colors used for communication scheduling
  int m_numColors;

  /// Color assigned to this rank [0, numColors)
  int m_color;
};

} // namespace geos

#endif // GEOS_MESH_MPICOMMUNICATIONS_DOMAINPARTITIONER_HPP_
