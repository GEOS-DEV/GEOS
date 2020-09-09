#include "PackCollection.hpp"

namespace geosx
{
PackCollection::PackCollection ( string const & name, Group * parent )
  : HistoryCollection( name, parent )
  , m_setsIndices( )
  , m_objectPath( )
  , m_fieldName( )
  , m_setNames( )
  , m_setSizeOverride( )
{
  registerWrapper( PackCollection::viewKeysStruct::objectPath, &m_objectPath )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "The name of the object from which to retrieve field values." );

  registerWrapper( PackCollection::viewKeysStruct::fieldName, &m_fieldName )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "The name of the (packable) field associated with the specified object to retrieve data from" );

  registerWrapper( PackCollection::viewKeysStruct::setNames, &m_setNames )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "The set(s) for which to retrieve data." );

  registerWrapper( PackCollection::viewKeysStruct::minSetsSize, &m_setSizeOverride)->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "The minimum size of the set(s) to be collected (use for sets that expand during the simulation).");
}

void PackCollection::InitializePostSubGroups( Group * const group )
{
  localIndex numSets = m_setNames.size( );
  m_collectionCount = numSets == 0 ? 1 : numSets;
  updateSetsIndices( dynamicCast< ProblemManager & >( *group ) );
  HistoryCollection::InitializePostSubGroups( group );
}

HistoryMetadata PackCollection::getMetadata( ProblemManager & pm, localIndex collectionIdx )
{
  DomainPartition & domain = *(pm.getDomainPartition( ));
  //MeshLevel const & meshLevel = *domain.getMeshBody( 0 )->getMeshLevel( 0 );
  Group const * target_object = this->getTargetObject( domain ); //meshLevel.GetGroupByPath( m_objectPath );
  WrapperBase const * target = target_object->getWrapperBase( m_fieldName );
  if( m_setNames.size() != 0 )
  {
    GEOSX_ERROR_IF( collectionIdx >= m_setNames.size(), "Invalid collection index specified." );
    localIndex num_indices = m_setSizeOverride > 0 ? m_setSizeOverride : m_setsIndices[ collectionIdx ].size( );
    HistoryMetadata metadata = target->getHistoryMetadata( num_indices );
    metadata.setName( metadata.getName() + " " + m_setNames[ collectionIdx ] );
    return metadata;
  }
  else
  {
    return target->getHistoryMetadata( m_setSizeOverride );
  }
}

void PackCollection::updateSetsIndices( ProblemManager & pm )
{
  DomainPartition & domain = *(pm.getDomainPartition( ));
  //MeshLevel const & meshLevel = *domain.getMeshBody( 0 )->getMeshLevel( 0 );
  //Group const * target_object = meshLevel.GetGroupByPath( m_objectPath );
  Group const * target_object = this->getTargetObject( domain );
  WrapperBase const * target = target_object->getWrapperBase( m_fieldName );
  GEOSX_ERROR_IF( !target->isPackable( false ), "The object targeted for collection must be packable!" );
  localIndex num_sets = m_setNames.size( );
  if( num_sets > 0 )
  {
    Group const * set_group = target_object->GetGroup( ObjectManagerBase::groupKeyStruct::setsString );
    m_setsIndices.resize( num_sets );
    localIndex set_idx = 0;
    for( auto & set_name : m_setNames )
    {
      dataRepository::Wrapper< SortedArray< localIndex > > const * const set_wrapper = set_group->getWrapper< SortedArray< localIndex > >( set_name );
      if( set_wrapper != nullptr )
      {
        SortedArrayView< localIndex const > const & set = set_wrapper->reference();
        m_setsIndices[ set_idx ].resize( 0 );
        if( set.size() > 0 )
        {
          m_setsIndices[ set_idx ].insert( 0, set.begin(), set.end() );
        }
      }
      set_idx++;
    }
  }
}

// TODO : once we add additional history collectors, this should likely be pulled into a super-class
Group const * PackCollection::getTargetObject( DomainPartition & domain )
{
  MeshLevel * const meshLevel = domain.getMeshBody( 0 )->getMeshLevel( 0 );
  dataRepository::Group * targetGroup = meshLevel;
  string_array const targetTokens = stringutilities::Tokenize( m_objectPath, "/" );
  localIndex const targetTokenLength = LvArray::integerConversion< localIndex >( targetTokens.size() );

  string processedPath;
  for( localIndex pathLevel = 0; pathLevel < targetTokenLength; ++pathLevel )
  {
    dataRepository::Group * const elemRegionSubGroup = targetGroup->GetGroup( ElementRegionManager::groupKeyStruct::elementRegionsGroup );
    if( elemRegionSubGroup != nullptr )
    {
      targetGroup = elemRegionSubGroup;
    }
    dataRepository::Group * const elemSubRegionSubGroup = targetGroup->GetGroup( ElementRegionBase::viewKeyStruct::elementSubRegions );
    if( elemSubRegionSubGroup != nullptr )
    {
      targetGroup = elemSubRegionSubGroup;
    }
    if( targetTokens[pathLevel] == ElementRegionManager::groupKeyStruct::elementRegionsGroup ||
        targetTokens[pathLevel] == ElementRegionBase::viewKeyStruct::elementSubRegions )
    {
      continue;
    }
    targetGroup = targetGroup->GetGroup( targetTokens[pathLevel] );
    processedPath += "/" + targetTokens[pathLevel];
    GEOSX_ERROR_IF( targetGroup == nullptr, "PackCollction::getTargetObject( ): Last entry in objectPath (" << processedPath << ") is not found" );
  }
  return targetGroup;
}

void PackCollection::collect( Group * domain_group,
                              real64 const GEOSX_UNUSED_PARAM( time_n ),
                              real64 const GEOSX_UNUSED_PARAM( dt ),
                              localIndex collectionIdx,
                              buffer_unit_type * & buffer )
{
  GEOSX_MARK_FUNCTION;
  GEOSX_ERROR_IF( collectionIdx >= getCollectionCount( ), "Attempting to collection from an invalid collection index!" );
  DomainPartition & domain = dynamicCast< DomainPartition & >( *domain_group );
  //MeshLevel const & meshLevel = *domain.getMeshBody( 0 )->getMeshLevel( 0 );
  //Group const * target_object = meshLevel.GetGroupByPath( m_objectPath );
  Group const * target_object = this->getTargetObject( domain );
  WrapperBase const * target = target_object->getWrapperBase( m_fieldName );
  if( m_setNames.size( ) > 0 )
  {
    localIndex sz = m_setsIndices.size( );
    if( sz > 0 )
    {
      // if we could directly transfer a sorted array to an array1d including on device this wouldn't require storing a copy of the indices
      target->PackByIndex( buffer, m_setsIndices[ collectionIdx ], false, true );
    }
  }
  else
  {
    target->Pack( buffer, false, true );
  }
}

REGISTER_CATALOG_ENTRY( TaskBase, PackCollection, std::string const &, Group * const )
}
