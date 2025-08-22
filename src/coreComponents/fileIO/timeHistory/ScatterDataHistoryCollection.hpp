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
 * @file ScatterDataHistoryCollection.hpp
 */

#ifndef GEOS_FILEIO_TIMEHISTORY_SCATTERDATAHISTORYCOLLECTION_HPP_
#define GEOS_FILEIO_TIMEHISTORY_SCATTERDATAHISTORYCOLLECTION_HPP_

#include "HistoryCollectionBase.hpp"
#include "ScatterDataProvider.hpp"
#include "common/DataTypes.hpp"

namespace geos
{
class DomainPartition;


/**
 * @brief Collects time history data from scattered points
 *
 * This class collects time history data from solvers that implement the
 * ScatterDataProvider interface.
 */
class ScatterDataHistoryCollection : public HistoryCollectionBase
{
public:

  /// Constructor
  ScatterDataHistoryCollection( string const & name, Group * const parent );

  /// Destructor
  virtual ~ScatterDataHistoryCollection() = default;

  /// Non-copyable
  ScatterDataHistoryCollection( ScatterDataHistoryCollection const & ) = delete;
  ScatterDataHistoryCollection & operator=( ScatterDataHistoryCollection const & ) = delete;

  /// Movable: default move constructor, deleted move assignment operator
  ScatterDataHistoryCollection( ScatterDataHistoryCollection && ) = default;
  ScatterDataHistoryCollection & operator=( ScatterDataHistoryCollection && ) = delete;


  /**
   * @brief Get the catalog name for factory registration
   * @return The catalog name
   */
  static string catalogName() { return "ScatterDataHistoryCollection"; }

  /**
   * @brief Initialize the collection
   */
  virtual void initializePostInitialConditionsPreSubGroups() override;

  /**
   * @brief Get the metadata for a specific collection index
   * @param domain The domain partition
   * @param collectionIdx The collection index
   * @return HistoryMetadata describing the data to be collected
   */
  virtual HistoryMetadata getMetaData( DomainPartition const & domain, localIndex collectionIdx ) const override;

  /**
   * @brief Get the number of discrete collection operations
   * @return Number of collection operations
   */
  virtual localIndex numCollectors() const override;

  /**
   * @brief Get the name of the target being collected
   * @return The target name (solver name)
   */
  virtual const string & getTargetName() const override;

  /**
   * @brief Get the number of metadata collectors
   * @return Number of metadata collectors (always 0 for scatter data)
   */
  virtual localIndex numMetaDataCollectors() const override;

  /**
   * @brief Get a metadata collector (not applicable for scatter data)
   * @param metaIdx The metadata collector index
   * @return Reference to metadata collector
   */
  virtual HistoryCollection & getMetaDataCollector( localIndex metaIdx ) override;

protected:

  /**
   * @brief Update set indices (required by base class, no-op for scatter data)
   * @param domain The domain partition
   */
  virtual void updateSetsIndices( DomainPartition const & domain ) override;

  /**
   * @brief Collect data from the scatter data provider
   * @param domain The domain partition containing the data
   */
  virtual void collect( DomainPartition const & domain );

  /**
   * @brief Collect data from a specific collection index into a buffer
   * @param domain The domain partition containing the data
   * @param collectionIdx The index of the collection operation
   * @param buffer The buffer to write data into
   */
  virtual void collect( DomainPartition const & domain,
                        localIndex const collectionIdx,
                        buffer_unit_type * & buffer ) override;

private:

  /// Pointer to the scatter data provider (physics solver)
  ScatterDataProvider const * m_scatterDataProvider;

  /// Name of the physics solver that provides scatter data
  string m_solverName;

  /// Flag to include coordinates in output (default: true)
  integer m_includeCoordinates;

  /// Flag to include metadata in output (default: false)
  integer m_includeMetadata;

  /**
   * @brief Find and validate the scatter data provider
   */
  void findScatterDataProvider();

  struct viewKeyStruct
  {
    static constexpr char const * solverNameString() { return "solverName"; }
    static constexpr char const * includeCoordinatesString() { return "includeCoordinates"; }
    static constexpr char const * includeMetadataString() { return "includeMetadata"; }
  };
};

} // namespace geos

#endif /* GEOS_FILEIO_TIMEHISTORY_SCATTERDATAHISTORYCOLLECTION_HPP_ */
