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
 * @brief Intermediate class for code factorisation.
 *        It mainly deals with collector and buffer management.
 *        It delegates the actual collection to derived classes.
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

#if defined(GEOSX_USE_PYGEOSX)
  /**
   * @brief Return PyHistoryCollection type.
   * @return Return PyHistoryCollection type.
   */
  virtual PyTypeObject * getPythonType() const override
  { return python::getPyHistoryCollectionType(); }
#endif

protected:

  ///@cond DO_NOT_DOCUMENT
  // Doxygen fails with error message `warning: documented empty return type of...`
  /**
   * @brief Update the indices from the sets being collected.
   * @param[in] domain The DomainPartition of the problem.
   */
  virtual void updateSetsIndices( DomainPartition const & domain ) = 0;
/// @endcond

  /**
   * @brief Retrieve the target object from the data repository.
   * @param[in] domain The DomainPartition of the problem.
   * @param[in] objectPath The data repo path of the target object.
   * @return The target object as a Group
   * @note If the object path is absolute this returns the target object by calling getGroupByPath. If the
   *        object path is relative, it searches relative to the mesh in the same way the fieldSpecification does,
   *        so any objectPath that works in fieldSpecification will also work here.
   */
  inline dataRepository::Group const * getTargetObject( DomainPartition const & domain,
                                                        string const & objectPath ) const
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
   * @param[in] domain The DomainPartition to collect time history on.
   * @param[in] collectionIdx The index of the collection operation to collect from the targeted collection event.
   * @param[in,out] buffer A properly-sized buffer to serialize history data into.
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
