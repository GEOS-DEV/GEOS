#include "PackCollection.hpp"
#include "SetIndexCollection.hpp"

namespace geosx
{
PackCollection::PackCollection ( string const & name, Group * parent )
  : HistoryCollection( name, parent )
  , m_setsIndices( )
  , m_objectPath( )
  , m_fieldName( )
  , m_setNames( )
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
  DomainPartition const & domain = *(pm.getDomainPartition( ));
  MeshLevel const & meshLevel = *domain.getMeshBody( 0 )->getMeshLevel( 0 );
  Group const * target_object = meshLevel.GetGroupByPath( m_objectPath );
  WrapperBase const * target = target_object->getWrapperBase( m_fieldName );
  if( m_setNames.size() != 0 )
  {
    GEOSX_ERROR_IF( collectionIdx >= m_setNames.size(), "Invalid collection index specified." );
    localIndex num_indices = m_setsIndices[ collectionIdx ].size( );
    HistoryMetadata metadata = target->getHistoryMetadata( num_indices );
    metadata.setName( metadata.getName() + " " + m_setNames[ collectionIdx ] );
    return metadata;
  }
  else
  {
    return target->getHistoryMetadata( -1 );
  }
}

void PackCollection::updateSetsIndices( ProblemManager & pm )
{
  DomainPartition const & domain = *(pm.getDomainPartition( ));
  MeshLevel const & meshLevel = *domain.getMeshBody( 0 )->getMeshLevel( 0 );
  Group const * target_object = meshLevel.GetGroupByPath( m_objectPath );
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

void PackCollection::collect( Group * domain_group,
                              real64 const GEOSX_UNUSED_PARAM( time_n ),
                              real64 const GEOSX_UNUSED_PARAM( dt ),
			      localIndex collectionIdx,
                              buffer_unit_type * & buffer )
{
  GEOSX_MARK_FUNCTION;
  GEOSX_ERROR_IF( collectionIdx >= getCollectionCount( ), "Attempting to collection from an invalid collection index!" );
  DomainPartition const & domain = dynamicCast< DomainPartition const & >( *domain_group );
  MeshLevel const & meshLevel = *domain.getMeshBody( 0 )->getMeshLevel( 0 );
  Group const * target_object = meshLevel.GetGroupByPath( m_objectPath );
  WrapperBase const * target = target_object->getWrapperBase( m_fieldName );
  if( m_setNames.size( ) > 0 )
  {
    Group const * set_group = target_object->GetGroup( ObjectManagerBase::groupKeyStruct::setsString );
    dataRepository::Wrapper< SortedArray< localIndex > > const * const set_wrapper = set_group->getWrapper< SortedArray< localIndex > >( m_setNames[collectionIdx] );
    if( set_wrapper != nullptr )
    {
      SortedArrayView< localIndex const > const & set = set_wrapper->reference();
      localIndex sz = set.size( );
      if( sz > 0 )
      {
	// if we could directly transfer a sorted array to an array1d including on device this wouldn't require storing a copy of the
	// indices internally
	target->PackByIndex( buffer, m_setsIndices[ collectionIdx ], false, true );
      }
    }
  }
  else
  {
    target->Pack( buffer, false, true );
  }
}

REGISTER_CATALOG_ENTRY( TaskBase, PackCollection, std::string const &, Group * const )
}
