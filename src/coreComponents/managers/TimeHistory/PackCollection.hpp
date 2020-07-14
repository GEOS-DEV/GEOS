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
 * @file PackCollection.hpp
 */

#ifndef GEOSX_PackCollection_HPP_
#define GEOSX_PackCollection_HPP_

#include "TimeHistoryCollection.hpp"

namespace geosx
{
/**
 * @class PackCollection
 *
 * A task class for serializing history information using the MPI communication packing routines.
 */
class PackCollection : public HistoryCollection
{
public:
  /**
   * @brief Constructor
   * @copydetails dataRepository::Group::Group( string const & name, Group * parent );
   */
  PackCollection ( string const & name, Group * parent );

  /**
   * @brief Catalog name interface
   * @return This type's catalog name
   */
  static string CatalogName() { return "PackCollection"; }

  /// @copydoc dataRepository::Group::InitializePostSubGroups
  void InitializePostSubGroups( Group * const group ) override;

  /// @copydoc geosx::HistoryCollection::GetMetadata
  virtual HistoryMetadata getMetadata( ProblemManager & problemManager ) override;

  /**
   * @brief Update the indices related to the sets being collected.
   * @param problemManager The problem manager cast to a group.
   * @note This is only required because we don't want to copy/move the
   *       indices each collection execution, becuase that causes data movement
   *       when collecting data from the device.
   * @note Refactoring the packing functions to allow direct usage of set indices
   *       from SortedArrayView instead of only ArrayViews will remove this
   *       duplication.
   */
  void UpdateSetsIndices( ProblemManager & problemManager );

  /**
   * @brief Count the total number of indices being collected by this process with this collector.
   * @param problemManager The ProblemManager.
   * @return The number of local indices being collected.
   */
  inline localIndex countLocalSetIndices( ProblemManager & problemManager );

  /**
   * @brief Count the number of indices being collected by this process for all sets up to and
   *        excluding the specified set index (see HistoryCollection::GetNumMetaCollectors).
   * @param problemManager The ProblemManager.
   * @param lastSetIdx The index of the Set to count all other local indices prior to.
   * @return The number of indices associate with all Sets locally prior to the specified Set.
   */
  localIndex countLocalSetIndicesExclusive( ProblemManager & problemManager, localIndex lastSetIdx = 0 );

  /// @copydoc geosx::HistoryCollection::GetNumMetaCollectors
  virtual localIndex getNumMetaCollectors( ) const override;

  /// @copydoc geosx::HistoryCollection::GetMetaCollector
  virtual std::unique_ptr< HistoryCollection > getMetaCollector( ProblemManager & problemManager, localIndex metaIdx, globalIndex metaRankOffset ) override;

  /// @copydoc geosx::HistoryCollection::Collect
  virtual void collect( Group * domain,
                        real64 const time_n,
                        real64 const dt,
                        buffer_unit_type * & buffer ) override;

  /// @cond DO_NOT_DOCUMENT
  struct viewKeysStruct
  {
    static constexpr auto objectPath = "objectPath";
    static constexpr auto fieldName = "fieldName";
    static constexpr auto setNames = "setNames";
  } keys;
  /// @endcond

protected:
  // todo : replace this with a vector of references to the actual set sortedarrays (after packing rework to allow sorted arrays to be used
  // for indexing)
  /// The indices for the specified sets to pack
  std::vector< array1d< localIndex > > m_setsIndices;
  /// The dataRepository name/path to get history data from
  string m_objectPath;
  /// The (packable) field associated with the specified object to get data from
  string m_fieldName;
  /// The names of the sets to collect history info from
  string_array m_setNames;
};

}
#endif
