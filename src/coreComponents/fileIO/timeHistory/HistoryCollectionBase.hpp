/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2020-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#ifndef GEOSX_HISTORYCOLLECTIONBASE_HPP
#define GEOSX_HISTORYCOLLECTIONBASE_HPP

#include "HistoryCollection.hpp"

#if defined(GEOSX_USE_PYGEOSX)
#include "fileIO/python/PyHistoryCollectionType.hpp"
#endif

namespace geosx
{

/**
 * @brief Technical code factorisation after refactoring. Feel free to challenge this class.
 */
class HistoryCollectionBase : public HistoryCollection
{
public:
  /// @copydoc geosx::dataRepository::Group::Group(string const & name, Group * const parent)
  HistoryCollectionBase( string const & name, Group * parent ):
    HistoryCollection( name, parent ),
    m_targetIsMeshObject( false ),
    m_collectionCount( 1 ),
    m_timeBufferProvider(),
    m_bufferProviders()
  {   }

  void initializePostSubGroups() override;

  localIndex numCollectors() const override;

  virtual bool execute( real64 const time_n,
                        real64 const dt,
                        integer const cycleNumber,
                        integer const eventCounter,
                        real64 const eventProgress,
                        DomainPartition & domain ) override;

  void registerBufferProvider( localIndex collectionIdx, BufferProvider bufferProvider ) override;

  HistoryMetadata getTimeMetaData() const override;

  void registerTimeBufferProvider( TimeBufferProvider timeBufferProvider ) override;

  HistoryCollection & getMetaDataCollector( localIndex metaIdx ) override;

  /**
   * @brief Return PyHistoryCollection type.
   * @return Return PyHistoryCollection type.
   */
#if defined(GEOSX_USE_PYGEOSX)
  virtual PyTypeObject * getPythonType() const override
  { return python::getPyHistoryCollectionType(); }
#endif

protected:

  /**
   * @brief Update the indices from the sets being collected.
   * @param[in] domain The DomainPartition of the problem.
   */
  virtual void updateSetsIndices( DomainPartition const & domain ) = 0;

  /**
   * @brief Retrieve the target object from the data repository.
   * @param domain The DomainPartition of the problem.
   * @param objectPath The data repo path of the target object.
   * @return The target object as a Group
   * @note If the object path is absolute this returns the target object by calling getGroupByPath. If the
   *        object path is relative, it searches relative to the mesh in the same way the fieldSpecification does,
   *        so any objectPath that works in fieldSpecification will also work here.
   */
  dataRepository::Group const * getTargetObject( DomainPartition const & domain, string const & objectPath ) const;

  /**
   * @brief Collect history information into the provided buffer. Typically called from HistoryCollection::execute .
   * @param domain The DomainPartition to collect time history on.
   * @param collectionIdx The index of the collection operation to collect from the targeted collection event.
   * @param buffer A properly-sized buffer to serialize history data into.
   */
  virtual void collect( DomainPartition const & domain,
                        localIndex const collectionIdx,
                        buffer_unit_type * & buffer ) = 0;

  /// whether the target object is associated with mesh entities (fields, etc)
  bool m_targetIsMeshObject;

  /// The number of discrete collection operations described by metadata this collection collects.
  localIndex m_collectionCount;

  /// Callbacks to get the current time buffer head to write time data into
  TimeBufferProvider m_timeBufferProvider;

  /// Callbacks to get the current buffer head to write history data into
  std::vector< BufferProvider > m_bufferProviders;

  /**
   * @brief The set of metadata collectors for this collector
   * @note Currently only used to collect coordinates of mesh objects when collecting field data.
   */
  std::vector< std::unique_ptr< HistoryCollection > > m_metaDataCollectors;
};

}

#endif // include guard
