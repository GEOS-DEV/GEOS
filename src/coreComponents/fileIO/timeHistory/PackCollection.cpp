#include "PackCollection.hpp"

namespace geosx
{
PackCollection::PackCollection ( string const & name, Group * parent )
  : HistoryCollectionBase( name, parent )
  , m_setsIndices( )
  , m_objectPath( )
  , m_fieldName( )
  , m_setNames( )
  , m_setChanged( true )
  , m_onlyOnSetChange( 0 )
  , m_disableCoordCollection( false )
  , m_initialized( false )
{
  registerWrapper( PackCollection::viewKeysStruct::objectPathString(), &m_objectPath ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "The name of the object from which to retrieve field values." );

  registerWrapper( PackCollection::viewKeysStruct::fieldNameString(), &m_fieldName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "The name of the (packable) field associated with the specified object to retrieve data from" );

  registerWrapper( PackCollection::viewKeysStruct::setNamesString(), &m_setNames ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "The set(s) for which to retrieve data." );

  registerWrapper( PackCollection::viewKeysStruct::onlyOnSetChangeString(), &m_onlyOnSetChange ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDefaultValue( 0 ).
    setDescription( "Whether or not to only collect when the collected sets of indices change in any way." );
}

void PackCollection::initializePostSubGroups( )
{
  if( !m_initialized )
  {
    DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
    localIndex numSets = m_setNames.size( );
    m_collectionCount = numSets == 0 ? 1 : numSets;
    // determine whether we're collecting from a mesh object manager
    Group const * const targetObject = this->getTargetObject( domain, m_objectPath );
    ObjectManagerBase const * const objectManagerTarget = dynamic_cast< ObjectManagerBase const * >( targetObject );
    m_targetIsMeshObject = objectManagerTarget != nullptr;
    // update sets after we know whether to filter ghost indices ( m_targetIsMeshObject )
    updateSetsIndices( domain );
    HistoryCollectionBase::initializePostSubGroups();
    // build any persistent meta collectors that are required
    buildMetaDataCollectors();
    for( auto & metaCollector : m_metaDataCollectors )
    {
      // coord meta collectors should have m_disableCoordCollection == true to avoid
      //  infinite recursive init calls here
      metaCollector->initializePostSubGroups();
    }
    m_initialized = true;
  }
}

HistoryMetadata PackCollection::getMetaData( DomainPartition const & domain, localIndex collectionIdx ) const
{
  Group const * targetObject = this->getTargetObject( domain, m_objectPath );
  WrapperBase const & targetField = targetObject->getWrapperBase( m_fieldName );
  if( m_setNames.empty() )
  {
    localIndex const packCount = m_setsIndices[0].size() == 0 ? 1 : m_setsIndices[0].size(); // TODO CHECK!
    return targetField.getHistoryMetadata( packCount );
//    return targetField.getHistoryMetadata(2);
  }
  else
  {
    GEOSX_ERROR_IF( collectionIdx < 0 || collectionIdx >= m_setNames.size(), "Invalid collection index specified." );
    localIndex collectionSize = m_setsIndices[ collectionIdx ].size( );
    if( (m_onlyOnSetChange != 0) && ( !m_setChanged ) ) // if we're only collecting when the set changes but the set hasn't changed
    {
      collectionSize = 0;
    }
    HistoryMetadata metadata = targetField.getHistoryMetadata( collectionSize );
    metadata.setName( metadata.getName() + " " + m_setNames[ collectionIdx ] );
    return metadata;
  }
}

// if we add a check for the setsizes changing we can use that data to determine whether the metadata collector
//  should only collect on size changes, otherwise not collect anything
void PackCollection::updateSetsIndices( DomainPartition const & domain )
{
  Group const * targetObject = this->getTargetObject( domain, m_objectPath );
  WrapperBase const & targetField = targetObject->getWrapperBase( m_fieldName );
  GEOSX_ERROR_IF( !targetField.isPackable( false ), "The object targeted for collection must be packable!" );
//  localIndex const numSets = m_setNames.size( );
//  std::vector< localIndex > oldSetSizes( numSets == 0 ? 1 : numSets );

  // If no set or "all" is specified we retrieve the entire field.
  // If sets are specified we retrieve the field only from those sets.
  bool const collectAll = m_setNames.empty() or std::find( m_setNames.begin(), m_setNames.end(), "all" ) != m_setNames.end();

//  std::set< string > setNames;
//  { // scope protection
//    Group const & setGroup = targetObject->getGroup( ObjectManagerBase::groupKeyStruct::setsString() );
//    auto predicate = [&setGroup]( string const & setName )
//    {
//      return setGroup.hasWrapper( setName );
//    };
//    std::copy_if( m_setNames.begin(), m_setNames.end(), std::inserter( setNames, setNames.end() ), predicate );
//  }

  localIndex const numSets = collectAll ? 1 : m_setNames.size();
//  localIndex const numSets = collectAll ? 1 : setNames.size();
  std::vector< localIndex > oldSetSizes( numSets );
  m_setsIndices.resize( numSets );

  if( collectAll )
  {
//    m_setsIndices.resize( numSets );
//    m_setsIndices[0].resize( targetField.size() );
//    for( localIndex ii = 0; ii < targetField.size(); ++ii )
    auto & si = m_setsIndices.front();
    oldSetSizes.front() = si.size( );
    si.resize( targetObject->size() );
    for( localIndex i = 0; i < targetObject->size(); ++i )
    {
      si[i] = i;
    }
    // this causes the ghost nodes to be filtered when collecting all
//    numSets = 1;
//    oldSetSizes[ 0 ] = si.size( );
//    oldSetSizes.front() = si.size( );
  }
  else
  {
    Group const & setGroup = targetObject->getGroup( ObjectManagerBase::groupKeyStruct::setsString() );
//    m_setsIndices.resize( numSets ); // TODO Warning refactoring because of this `numSets`, check twice.
    localIndex setIdx = 0;
    for( auto const & setName : m_setNames )
//    for( auto const & setName : setNames )
    {
      if( setGroup.hasWrapper( setName ) )
      {
        SortedArrayView< localIndex const > const & set = setGroup.getReference< SortedArray< localIndex > >( setName );
        oldSetSizes[ setIdx ] = m_setsIndices[ setIdx ].size( );
        localIndex setSize = set.size( );
        m_setsIndices[ setIdx ].resize( setSize );
        if( setSize > 0 )
        {
          for( localIndex idx = 0; idx < setSize; ++idx )
          {
            m_setsIndices[ setIdx ][ idx ] = set[ idx ];
          }
        }
      }
      setIdx++;
    }
  }

  // filter out the ghost indices immediately when we update the index sets
  if( m_targetIsMeshObject )
  {
    ObjectManagerBase const * objectManagerTarget = dynamic_cast< ObjectManagerBase const * >( targetObject );
    arrayView1d< integer const > const ghostRank = objectManagerTarget->ghostRank( );
//    for( localIndex setIdx = 0; setIdx < numSets; ++setIdx )
    for( std::size_t setIdx = 0; setIdx < m_setsIndices.size(); ++setIdx )
    {
      array1d< localIndex > ownedIndices( m_setsIndices[ setIdx ].size() );
      localIndex ownIdx = 0;
      for( localIndex idx = 0; idx < m_setsIndices[ setIdx ].size(); ++idx )
      {
        if( ghostRank[ m_setsIndices[ setIdx ][ idx ] ] < 0 )
        {
          ownedIndices[ ownIdx ] = m_setsIndices[ setIdx ][ idx ];
          ++ownIdx;
        }
      }
      localIndex newSetSize = ownIdx;
      if( oldSetSizes[ setIdx ] != newSetSize )
      {
        m_setChanged = true;
      }
      if( !m_setChanged )  // if the size hasn't changed check the individual values
      {
        for( localIndex idx = 0; idx < ownedIndices.size(); ++idx )
        {
          if( m_setsIndices[ setIdx ][ idx ] != ownedIndices[ idx ] )
          {
            m_setChanged = true;
            break;
          }
        }
      }
      m_setsIndices[ setIdx ].resize( newSetSize );
      for( localIndex idx = 0; idx < newSetSize; ++idx )
      {
        m_setsIndices[ setIdx ][ idx ] = ownedIndices[ idx ];
      }
    }
  }
}

localIndex PackCollection::numMetaDataCollectors( ) const
{
  return 1; // TODO This is unclear.
//  return m_targetIsMeshObject && !m_disableCoordCollection ? 1 : 0;
}

void PackCollection::buildMetaDataCollectors( )
{
  if( !m_disableCoordCollection )
  {
    char const * coordField = nullptr;
    if( m_objectPath.find( "nodeManager" ) != string::npos )
    {
      coordField = NodeManager::viewKeyStruct::referencePositionString();
    }
    else if( m_objectPath.find( "edgeManager" ) != string::npos )
    {
      GEOSX_ERROR( "Edge coordinate data collection is unimplemented." );
    }
    else if( m_objectPath.find( "faceManager" ) != string::npos )
    {
      coordField = FaceManager::viewKeyStruct::faceCenterString();
    }
    else if( m_objectPath.find( "ElementRegions" ) != string::npos )
    {
      coordField = ElementSubRegionBase::viewKeyStruct::elementCenterString();
    }
//    else
//    {
//      GEOSX_ERROR( "Data collection for \"" << m_objectPath << "\" is unimplemented." );
//    }

    string metaName( "coordinates" ); // FIXME Name?
    std::unique_ptr< PackCollection > collector = std::make_unique< PackCollection >( metaName, this );
    collector->m_objectPath = m_objectPath;
    collector->m_fieldName = coordField ? string( coordField ) : m_fieldName;
    collector->m_setNames = m_setNames;
    collector->m_onlyOnSetChange = true;
    // don't recursively keep creating metaDataCollectors
    collector->disableCoordCollection();
    m_metaDataCollectors.push_back( std::move( collector ) );
  }
}

void PackCollection::collect( DomainPartition const & domain,
                              localIndex const collectionIdx,
                              buffer_unit_type * & buffer )
{
  GEOSX_MARK_FUNCTION;
  GEOSX_ERROR_IF( collectionIdx < 0 || collectionIdx >= numCollectors(), "Attempting to collection from an invalid collection index!" );
  Group const * targetObject = this->getTargetObject( domain, m_objectPath );
  WrapperBase const & targetField = targetObject->getWrapperBase( m_fieldName );
  // if we have any indices to collect, and we're either collecting every time or we're only collecting when the set changes and the set has
  // changed
  parallelDeviceEvents events;
  if( m_setsIndices[ collectionIdx ].size() > 0 )
  {
    if( ( (m_onlyOnSetChange != 0) && m_setChanged ) || (m_onlyOnSetChange == 0) )
    {
      targetField.packByIndex< true >( buffer, m_setsIndices[collectionIdx], false, true, events );
    }
  }
  // if we're not collecting from a set of indices, we're collecting the entire object
  //  this will only happen when we're not targeting a mesh object since in that case while setnames size is 0,
  //  setsIndices[0] is the entire non-ghost index set, so all mesh object collection goes to packbyindex and any
  //  non-mesh objects (that don't somehow have index sets) are packed in their entirety
  else if( !m_targetIsMeshObject && m_setNames.size() == 0 )
  {
    targetField.pack< true >( buffer, false, true, events );
  }
  m_setChanged = false;
  GEOSX_ASYNC_WAIT( 6000000000, 10, testAllDeviceEvents( events ) );
}

REGISTER_CATALOG_ENTRY( TaskBase, PackCollection, string const &, Group * const )
}
