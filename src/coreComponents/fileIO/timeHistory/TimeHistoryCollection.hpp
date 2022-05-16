/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file TimeHistoryCollection.hpp
 */

#ifndef GEOSX_FILEIO_TIMEHISTORY_TIMEHISTORYCOLLECTION_HPP_
#define GEOSX_FILEIO_TIMEHISTORY_TIMEHISTORYCOLLECTION_HPP_

#include "dataRepository/BufferOpsDevice.hpp"
#include "dataRepository/HistoryDataSpec.hpp"
#include "events/tasks/TaskBase.hpp"
#include "mesh/DomainPartition.hpp"

#if defined(GEOSX_USE_PYGEOSX)
#include "fileIO/python/PyHistoryCollectionType.hpp"
#endif

#include <functional>

namespace geosx
{

class DomainPartition;

using namespace dataRepository;

/**
 * @class HistoryCollection
 *
 * A task class for serializing time history data into a buffer for later I/O.
 */
class HistoryCollection : public TaskBase
{
public:
  /// @copydoc geosx::dataRepository::Group::Group(string const & name, Group * const parent)
  HistoryCollection( string const & name, Group * parent ):
    TaskBase( name, parent ),
    m_targetIsMeshObject( false ),
    m_collectionCount( 1 ),
    m_timeBufferCall(),
    m_bufferCalls()
  {   }


  void initializePostSubGroups() override
  {
    m_bufferCalls.resize( m_collectionCount );
  }

  /**
   * @brief Get the number of discrete collection operations this collector conducts.
   * @return The number of collection operations for this collector.
   */
  virtual localIndex getCollectionCount( ) const
  {
    return m_collectionCount;
  }

  /**
   * @brief Get the metadata for what this collector collects.
   * @param domain The DomainPartition.
   * @param collectionIdx Which collected item to get metadata for.
   * @return A HistoryMetadata object describing  the history data being collected by this collector.
   */
  virtual HistoryMetadata getMetadata( DomainPartition const & domain, localIndex collectionIdx )
  {
    GEOSX_UNUSED_VAR( domain );
    GEOSX_UNUSED_VAR( collectionIdx );
    return HistoryMetadata( );
  }

  /**
   * @brief Get the name of the object being targeted for collection.
   * @return The collection target's name
   */
  virtual const string & getTargetName( ) const = 0;

  /**
   * @brief Collects history data.
   * @copydoc EventBase::execute()
   */
  virtual bool execute( real64 const time_n,
                        real64 const dt,
                        integer const cycleNumber,
                        integer const eventCounter,
                        real64 const eventProgress,
                        DomainPartition & domain ) override
  {
    for( localIndex collectionIdx = 0; collectionIdx < getCollectionCount(); ++collectionIdx )
    {
      // std::function defines the == and =! comparable against nullptr_t to check the
      //  function pointer is actually assigned (an error would be thrown on the call attempt even so)
      GEOSX_ERROR_IF( m_bufferCalls[collectionIdx] == nullptr,
                      "History collection buffer retrieval function is unassigned, did you declare a related TimeHistoryOutput event?" );
      // using GEOSX_ERROR_IF_EQ caused type issues since the values are used in streams
      // we don't explicitly update the index sets here as they are updated in the buffer callback (see
      // TimeHistoryOutput.cpp::initCollectorParallel( ) )
      buffer_unit_type * buffer = m_bufferCalls[collectionIdx]( );
      collect( domain, time_n, dt, collectionIdx, buffer );
    }
    for( auto & metaCollector : m_metaCollectors )
    {
      metaCollector->execute( time_n, dt, cycleNumber, eventCounter, eventProgress, domain );
    }
    int rank = MpiWrapper::commRank();
    if( rank == 0 && m_timeBufferCall )
    {
      buffer_unit_type * timeBuffer = m_timeBufferCall();
      memcpy( timeBuffer, &time_n, sizeof(time_n) );
    }
    return false;
  }

  /**
   * @brief Register a callback that gives the current head of the time history data buffer.
   * @param collectionIdx Which collection item to register the buffer callback for.
   * @param bufferCall A functional that when invoked returns a pointer to the head of a buffer at least large enough to
   *                    serialize one timestep of history data into.
   * @note This is typically meant to callback to BufferedHistoryIO::GetBufferHead( )
   */
  void registerBufferCall( localIndex collectionIdx, std::function< buffer_unit_type *() > bufferCall )
  {
    GEOSX_ERROR_IF( collectionIdx < 0 || collectionIdx >= this->getCollectionCount( ), "Invalid collection index specified." );
    m_bufferCalls[collectionIdx] = bufferCall;
  }

  /**
   * @brief Get a metadata object relating the the Time variable itself.
   * @return A HistroyMetadata object describing the Time variable.
   */
  HistoryMetadata getTimeMetadata( ) const
  {
    return HistoryMetadata( this->getTargetName( ) + " Time", 1, std::type_index( typeid(real64) ) );
  }

  /**
   * @brief Register a callback that gives the current head of the time data buffer.
   * @param timeBufferCall A functional that when invoked returns a pointer to the head of a buffer at least large enough to
   *                    serialize one instance of the Time variable into.
   * @note This is typically meant to callback to BufferedHistoryIO::GetBufferHead( )
   */
  void registerTimeBufferCall( std::function< buffer_unit_type *() > timeBufferCall )
  {
    m_timeBufferCall = timeBufferCall;
  }

  /**
   * @brief Get the number of collectors of meta-information (set indices, etc) writing time-independent information during initialization.
   * @return The number of collectors of meta-information for this collector.
   */
  virtual localIndex getNumMetaCollectors( ) const
  {
    return 0;
  }

  /**
   * @brief Get a pointer to a collector of meta-information for this collector.
   * @param metaIdx Which of the meta-info collectors to return. (see HistoryCollection::getNumMetaCollectors( ) ).
   * @return A unique pointer to the HistoryCollection object used for meta-info collection. Intented to fall out of scope and desctruct
   * immediately
   *         after being used to perform output during simulation initialization.
   */
  virtual HistoryCollection & getMetaCollector( localIndex metaIdx )
  {
    GEOSX_ASSERT_MSG( metaIdx >= 0 && metaIdx < getNumMetaCollectors( ), "Requesting nonexistent meta collector index." );
    return *m_metaCollectors[ metaIdx ].get( );
  }
  /**
   * @brief Update the indices from the sets being collected.
   * @param domain The DomainPartition of the problem.
   */
  virtual void updateSetsIndices ( DomainPartition & domain ) = 0;

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
   * @brief Retrieve the target object from the data repository.
   * @param domain The DomainPartition of the problem.
   * @param objectPath The data repo path of the target object.
   * @return The target object as a Group
   * @note If the object path is absolute this returns the target object by calling getGroupByPath. If the
   *        object path is relative, it searches relative to the mesh in the same way the fieldSpecification does,
   *        so any objectPath that works in fieldSpecification will also work here.
   */
  inline dataRepository::Group const * getTargetObject( DomainPartition const & domain,
                                                        string const & objectPath )
  {
    // absolute objectPaths can be used to target anything in the data repo that is packable
    if( objectPath[0] == '/' )
    {
      return &Group::getGroupByPath( objectPath );
    }
    else // relative objectPaths use relative lookup identical to fieldSpecification to make xml input spec easier
    {
      string_array targetTokens = stringutilities::tokenize( objectPath, "/" );
      localIndex targetTokenLength = LvArray::integerConversion< localIndex >( targetTokens.size() );

      dataRepository::Group const * targetGroup = nullptr;
      int const numMeshBodies = domain.getMeshBodies().numSubGroups();


      if( numMeshBodies==1 )
      {
        string const singleMeshBodyName = domain.getMeshBody( 0 ).getName();
        if( targetTokens[0] != singleMeshBodyName )
        {
          ++targetTokenLength;
          targetTokens.insert( 0, &singleMeshBodyName, (&singleMeshBodyName)+1 );
        }
      }
      else
      {
        bool bodyFound = false;
        domain.forMeshBodies( [&]( MeshBody const & meshBody )
        {
          if( meshBody.getName()==targetTokens[0] )
          {
            bodyFound=true;
          }
        } );

        GEOSX_ERROR_IF( !bodyFound,
                        GEOSX_FMT( "MeshBody ({}) is specified, but not found.",
                                   targetTokens[0] ) );
      }



      string const meshBodyName = targetTokens[0];
      MeshBody const & meshBody = domain.getMeshBody( meshBodyName );

      // set mesh level in path
      localIndex const numMeshLevels = meshBody.getMeshLevels().numSubGroups();
      if( numMeshLevels==1 )
      {
        string const singleMeshLevelName = meshBody.getMeshLevels().getGroup< MeshLevel >( 0 ).getName();
        if( targetTokens[1] != singleMeshLevelName )
        {
          ++targetTokenLength;
          targetTokens.insert( 1, &singleMeshLevelName, (&singleMeshLevelName)+1 );
        }
        else
        {
          bool levelFound = false;
          meshBody.forMeshLevels( [&]( MeshLevel const & meshLevel )
          {
            if( meshLevel.getName()==targetTokens[1] )
            {
              levelFound=true;
            }
          } );

          GEOSX_ERROR_IF( !levelFound,
                          GEOSX_FMT( "MeshLevel ({}) is specified, but not found.",
                                     targetTokens[1] ) );
        }
      }
      else if( !meshBody.getMeshLevels().hasGroup< MeshLevel >( targetTokens[1] ) )
      {
        GEOSX_LOG_RANK_0( "In TimeHistoryCollection.hpp, Mesh Level Discretization not specified, "
                          "using baseDiscretizationString()." );

        string const baseMeshLevelName = MeshBody::groupStructKeys::baseDiscretizationString();
        ++targetTokenLength;
        targetTokens.insert( 1, &baseMeshLevelName, (&baseMeshLevelName)+1 );
      }

      string meshLevelName = targetTokens[1];
      MeshLevel const & meshLevel = meshBody.getMeshLevel( meshLevelName );
      targetGroup = &meshLevel;


      if( targetTokens[2]== MeshLevel::groupStructKeys::elemManagerString() )
      {
        ElementRegionManager const & elemRegionManager = meshLevel.getElemManager();

        string const elemRegionName = targetTokens[3];
        ElementRegionBase const & elemRegion = elemRegionManager.getRegion( elemRegionName );

        string const elemSubRegionName = targetTokens[4];
        ElementSubRegionBase const & elemSubRegion = elemRegion.getSubRegion( elemSubRegionName );
        targetGroup = &elemSubRegion;
      }
      else
      {
        for( localIndex pathLevel = 2; pathLevel < targetTokenLength; ++pathLevel )
        {
          targetGroup = targetGroup->getGroupPointer( targetTokens[pathLevel] );
        }
      }
      return targetGroup;
    }
  }

  /**
   * @brief Collect history information into the provided buffer. Typically called from HistoryCollection::execute .
   * @param domain The DomainPartition to collect time history on.
   * @param time_n The current simulation time.
   * @param dt The current simulation time delta.
   * @param collectionIdx The index of the collection operation to collect from the targeted collection event.
   * @param buffer A properly-sized buffer to serialize history data into.
   */
  virtual void collect( DomainPartition & domain, real64 const time_n, real64 const dt, localIndex const collectionIdx, buffer_unit_type * & buffer ) = 0;

protected:
  /// whether the target object is associated with mesh entities (fields, etc)
  bool m_targetIsMeshObject;
  /// The number of discrete collection operations described by metadata this collection collects.
  localIndex m_collectionCount;
  /// Callbacks to get the current time buffer head to write time data into
  std::function< buffer_unit_type *() > m_timeBufferCall;
  /// Callbacks to get the current buffer head to write history data into
  std::vector< std::function< buffer_unit_type *() > > m_bufferCalls;
  /// The set of metadata collectors for this collector ( currently only used to collect coordinates of mesh
  ///   objects when collecting field data )
  std::vector< std::unique_ptr< HistoryCollection > > m_metaCollectors;
};

}

#endif
