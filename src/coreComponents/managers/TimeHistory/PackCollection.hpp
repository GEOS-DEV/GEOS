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
  virtual HistoryMetadata GetMetadata( Group * problem_group ) override;

  /**
   * @brief Update the indices related to the sets being collected.
   * @param problem_group The problem manager cast to a group.
   * @note This is only required because we don't want to copy/move the
   *       indices each collection execution, becuase that causes data movement
   *       when collecting data from the device.
   * @note Refactoring the packing functions to allow direct usage of set indices
   *       from SortedArrayView instead of only ArrayViews will remove this
   *       duplication.
   */
  void UpdateSetsIndices( Group * problem_group );

  /**
   * @brief Count the total number of indices being collected by this process with this collector.
   * @param problem_group The ProblemManager cast to a Group.
   * @return The number of local indices being collected.
   */
  inline localIndex CountLocalSetIndices( Group * problem_group );

  /**
   * @brief Count the number of indices being collected by this process for all sets up to and
   *        excluding the specified set index (see HistoryCollection::GetNumMetaCollectors).
   * @param problem_group The ProblemManager cast to a Group.
   * @param last_set_idx The index of the Set to count all other local indices prior to.
   * @return The number of indices associate with all Sets locally prior to the specified Set.
   */
  localIndex CountLocalSetIndicesExclusive( Group * problem_group, localIndex last_set_idx = 0 );

  /// @copydoc geosx::HistoryCollection::GetNumMetaCollectors
  virtual localIndex GetNumMetaCollectors( ) const override;

  /// @copydoc geosx::HistoryCollection::GetMetaCollector
  virtual std::unique_ptr< HistoryCollection > GetMetaCollector( Group * problem_group, localIndex meta_idx, globalIndex meta_rank_offset ) override;

  /// @copydoc geosx::HistoryCollection::Collect
  virtual void Collect( Group * domain,
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
  std::vector< array1d< localIndex > > m_sets_indices;
  /// The dataRepository name/path to get history data from
  string m_object_path;
  /// The (packable) field associated with the specified object to get data from
  string m_field_name;
  /// The names of the sets to collect history info from
  string_array m_set_names;
};

}
#endif
