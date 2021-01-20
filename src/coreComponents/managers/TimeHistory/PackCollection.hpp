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
   * @copydetails dataRepository::Group::Group(string const & name, Group * parent);
   */
  PackCollection (string const & name, Group * parent);

  /**
   * @brief Catalog name interface
   * @return This type's catalog name
   */
  static string CatalogName() {return "PackCollection";}

  /// @copydoc dataRepository::Group::InitializePostSubGroups
  void InitializePostSubGroups(Group * const group) override;

  /// @copydoc geosx::HistoryCollection::getMetadata
  virtual HistoryMetadata getMetadata(ProblemManager & problemManager, localIndex collectionIdx) override;

  /// @copydoc geosx::HistoryCollection::getTargetName
  virtual const string & getTargetName() const override
  {
    return m_fieldName;
  }

  /**
   * @brief Update the indices related to the sets being collected.
   * @param domain The domain partition.
   * @note This is only required because we don't want to copy/move the
   *       indices each collection execution, because that causes data movement
   *       when collecting data from the device.
   * @note Refactoring the packing functions to allow direct usage of set indices
   *       from SortedArrayView instead of only ArrayViews will remove this
   *       duplication.
   */
  virtual void updateSetsIndices(DomainPartition & domain) override final;

  /**
   * @brief Filters out ghost rank indices from setIndices to be collected.
   * @param setIndex which set (collection item) needs to be filtered
   * @param set the set of indices without the ghost ones
   * @param ghostRank the ghost rank of each index for the target object
   */
  void filterGhostIndices(localIndex const setIndex,
                           array1d<localIndex> & set,
                           arrayView1d<integer const> const & ghostRank);

  /// @cond DO_NOT_DOCUMENT
  struct viewKeysStruct
  {
    static constexpr auto objectPath = "objectPath";
    static constexpr auto fieldName = "fieldName";
    static constexpr auto setNames = "setNames";
    static constexpr auto minSetSize = "minSetSize";
  } keys;
  /// @endcond

protected:
  /**
   * @brief Get the target object from which to collect time history data.
   * @param domain The problem domain.
   * @return The target object as an ObjectManager.
   */
  ObjectManagerBase const * getTargetObject(DomainPartition & domain);

  /// @copydoc geosx::HistoryCollection::collect
  virtual void collect(DomainPartition & domain,
                        real64 const time_n,
                        real64 const dt,
                        localIndex const collectionIdx,
                        buffer_unit_type * & buffer) override;

private:
  // todo : replace this with a vector of references to the actual set sortedarrays (after packing rework to allow sorted arrays to be used
  // for indexing)
  /// The indices for the specified sets to pack
  std::vector<array1d<localIndex>> m_setsIndices;
  /// The dataRepository name/path to get history data from
  string m_objectPath;
  /// The (packable) field associated with the specified object to get data from
  string m_fieldName;
  /// The names of the sets to collect history info from
  string_array m_setNames;
  /// Set the minimum size of all the sets being collected on each process
  localIndex m_minimumSetSize;
};

}
#endif
