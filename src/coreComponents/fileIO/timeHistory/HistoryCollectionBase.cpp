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

#include "HistoryCollectionBase.hpp"

namespace geosx {

void HistoryCollectionBase::initializePostSubGroups()
{
  m_bufferProviders.resize( m_collectionCount );
}

localIndex HistoryCollectionBase::numCollectors() const
{
  return m_collectionCount;
}

void HistoryCollectionBase::registerBufferProvider( localIndex collectionIdx,
                                                    BufferProvider bufferProvider )
{
  GEOSX_ERROR_IF( collectionIdx < 0 || collectionIdx >= this->numCollectors(), "Invalid collection index specified." );
  m_bufferProviders[collectionIdx] = bufferProvider;
}

HistoryCollection & HistoryCollectionBase::getMetaDataCollector( localIndex metaIdx )
{
  GEOSX_ASSERT_MSG( metaIdx >= 0 && metaIdx < numMetaDataCollectors(), "Requesting nonexistent meta collector index." );
  return *m_metaDataCollectors[ metaIdx ].get( );
}

HistoryMetadata HistoryCollectionBase::getTimeMetaData() const
{
  return HistoryMetadata( this->getTargetName() + " Time", 1, std::type_index( typeid( real64 ) ) );
}

void HistoryCollectionBase::registerTimeBufferProvider( TimeBufferProvider timeBufferProvider )
{
  m_timeBufferProvider = timeBufferProvider;
}

bool HistoryCollectionBase::execute( real64 const time_n,
                                     real64 const dt,
                                     integer const cycleNumber,
                                     integer const eventCounter,
                                     real64 const eventProgress,
                                     DomainPartition & domain )
{
  for( localIndex collectionIdx = 0; collectionIdx < numCollectors(); ++collectionIdx )
  {
    // std::function defines the == and =! comparable against nullptr_t to check the
    // function pointer is actually assigned (an error would be thrown on the call attempt even so)
    GEOSX_ERROR_IF( m_bufferProviders[collectionIdx] == nullptr,
                    "History collection buffer retrieval function is unassigned, did you declare a related TimeHistoryOutput event?" );
    // using GEOSX_ERROR_IF_EQ caused type issues since the values are used in streams
    // we don't explicitly update the index sets here as they are updated in the buffer callback (see
    // TimeHistoryOutput.cpp::initCollectorParallel( ) )
    this->updateSetsIndices( domain );
    HistoryMetadata hmd = this->getMetaData( domain, collectionIdx );
    buffer_unit_type * buffer = m_bufferProviders[collectionIdx]( hmd.size( 0 ) );
    collect( domain, time_n, dt, collectionIdx, buffer );
  }
  for( auto & metaCollector: m_metaDataCollectors )
  {
    metaCollector->execute( time_n, dt, cycleNumber, eventCounter, eventProgress, domain );
  }
  int rank = MpiWrapper::commRank();
  if( rank == 0 && m_timeBufferProvider )
  {
    buffer_unit_type * timeBuffer = m_timeBufferProvider();
    memcpy( timeBuffer, &time_n, sizeof( time_n ) );
  }
  return false;
}

dataRepository::Group const * HistoryCollectionBase::getTargetObject( DomainPartition const & domain, string const & objectPath ) const
{
  // absolute objectPaths can be used to target anything in the data repo that is packable
  if( objectPath[0] == '/' )
  {
    return &Group::getGroupByPath( objectPath );
  }
  else // relative objectPaths use relative lookup identical to fieldSpecification to make xml input spec easier
  {
    string_array const targetTokens = stringutilities::tokenize( objectPath, "/" );
    localIndex const targetTokenLength = LvArray::integerConversion< localIndex >( targetTokens.size() );

    dataRepository::Group const * targetGroup = &domain.getMeshBody( 0 ).getMeshLevel( 0 );
    for( localIndex pathLevel = 0; pathLevel < targetTokenLength; ++pathLevel )
    {
      dataRepository::Group const * elemRegionSubGroup = targetGroup->getGroupPointer( ElementRegionManager::groupKeyStruct::elementRegionsGroup() );
      if( elemRegionSubGroup != nullptr )
      {
        targetGroup = elemRegionSubGroup;
      }
      dataRepository::Group const * elemSubRegionSubGroup = targetGroup->getGroupPointer( ElementRegionBase::viewKeyStruct::elementSubRegions() );
      if( elemSubRegionSubGroup != nullptr )
      {
        targetGroup = elemSubRegionSubGroup;
      }
      if( targetTokens[pathLevel] == ElementRegionManager::groupKeyStruct::elementRegionsGroup() ||
          targetTokens[pathLevel] == ElementRegionBase::viewKeyStruct::elementSubRegions() )
      {
        continue;
      }
      targetGroup = &targetGroup->getGroup( targetTokens[pathLevel] );
    }
    return targetGroup;
  }
}


}