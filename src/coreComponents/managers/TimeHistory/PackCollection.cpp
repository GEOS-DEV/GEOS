#include "PackCollection.hpp"
#include "SetIndexCollection.hpp"

namespace geosx
{
PackCollection::PackCollection ( string const & name, Group * parent )
  : HistoryCollection( name, parent )
  , m_sets_indices( )
  , m_object_path( )
  , m_field_name( )
  , m_set_names( )
{
  registerWrapper( PackCollection::viewKeysStruct::objectPath, &m_object_path )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "The name of the object from which to retrieve field values." );

  registerWrapper( PackCollection::viewKeysStruct::fieldName, &m_field_name )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "The name of the (packable) field associated with the specified object to retrieve data from" );

  registerWrapper( PackCollection::viewKeysStruct::setNames, &m_set_names )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "The set(s) for which to retrieve data." );
}

void PackCollection::InitializePostSubGroups( Group * const group )
{
  UpdateSetsIndices( group );
}

/// @copydoc HistoryCollection::GetMetadata
HistoryMetadata PackCollection::GetMetadata( Group * problem_group )
{
  localIndex num_indices = CountLocalSetIndices( problem_group );
  ProblemManager const * pm = dynamicCast< ProblemManager const * >( problem_group );
  DomainPartition const & domain = *pm->getDomainPartition( );
  MeshLevel const & meshLevel = *domain.getMeshBody( 0 )->getMeshLevel( 0 );
  Group const * target_object = meshLevel.GetGroupByPath( m_object_path );
  WrapperBase const * target = target_object->getWrapperBase( m_field_name );
  return target->getBufferedIOMetadata( num_indices );
}

void PackCollection::UpdateSetsIndices( Group * problem_group )
{
  ProblemManager const * pm = dynamicCast< ProblemManager const * >( problem_group );
  DomainPartition const & domain = *pm->getDomainPartition( );
  MeshLevel const & meshLevel = *domain.getMeshBody( 0 )->getMeshLevel( 0 );
  Group const * target_object = meshLevel.GetGroupByPath( m_object_path );
  WrapperBase const * target = target_object->getWrapperBase( m_field_name );
  GEOSX_ERROR_IF( !target->isPackable( false ), "The object targeted for collection must be packable!" );
  localIndex num_sets = m_set_names.size( );
  if( num_sets > 0 )
  {
    Group const * set_group = target_object->GetGroup( ObjectManagerBase::groupKeyStruct::setsString );
    m_sets_indices.resize( num_sets );
    localIndex set_idx = 0;
    for( auto & set_name : m_set_names )
    {
      dataRepository::Wrapper< SortedArray< localIndex > > const * const set_wrapper = set_group->getWrapper< SortedArray< localIndex > >( set_name );
      if( set_wrapper != nullptr )
      {
        SortedArrayView< localIndex const > const & set = set_wrapper->reference();
        m_sets_indices[ set_idx ].insert( 0, set.begin(), set.end() );
      }
      set_idx++;
    }
  }
}

inline localIndex PackCollection::CountLocalSetIndices( Group * problem_group )
{
  return CountLocalSetIndicesExclusive( problem_group, m_set_names.size( ) );
}

localIndex PackCollection::CountLocalSetIndicesExclusive( Group * problem_group, localIndex last_set_idx )
{
  GEOSX_ERROR_IF( m_set_names.size() == 0 && last_set_idx > 0, "No set names specified in input, but trying to sum local set indices." );
  int num_sets = m_sets_indices.size( );
  if( num_sets == 0 )
  {
    UpdateSetsIndices( problem_group );
  }
  localIndex num_indices = 0;
  for( localIndex set_idx = 0; set_idx < last_set_idx; ++set_idx )
  {
    num_indices += m_sets_indices[ set_idx ].size( );
  }
  return num_indices;
}

localIndex PackCollection::GetNumMetaCollectors( ) const
{
  return m_set_names.size( );
}

std::unique_ptr< HistoryCollection > PackCollection::GetMetaCollector( Group * problem_group, localIndex meta_idx, globalIndex meta_rank_offset )
{
  globalIndex local_set_offset = LvArray::integerConversion< globalIndex >( CountLocalSetIndicesExclusive( problem_group, meta_idx ) );
  return std::unique_ptr< HistoryCollection >( new SetIndexCollection( m_object_path, m_set_names[meta_idx], meta_rank_offset + local_set_offset ));
}

void PackCollection::Collect( Group * domain_group,
                              real64 const GEOSX_UNUSED_PARAM( time_n ),
                              real64 const GEOSX_UNUSED_PARAM( dt ),
                              buffer_unit_type * & buffer )
{
  DomainPartition const & domain = dynamicCast< DomainPartition const & >( *domain_group );
  MeshLevel const & meshLevel = *domain.getMeshBody( 0 )->getMeshLevel( 0 );
  Group const * target_object = meshLevel.GetGroupByPath( m_object_path );
  WrapperBase const * target = target_object->getWrapperBase( m_field_name );
  if( m_set_names.size( ) > 0 )
  {
    Group const * set_group = target_object->GetGroup( ObjectManagerBase::groupKeyStruct::setsString );
    localIndex set_idx = 0;
    for( auto & set_name : m_set_names )
    {
      dataRepository::Wrapper< SortedArray< localIndex > > const * const set_wrapper = set_group->getWrapper< SortedArray< localIndex > >( set_name );
      if( set_wrapper != nullptr )
      {
        SortedArrayView< localIndex const > const & set = set_wrapper->reference();
        localIndex sz = set.size( );
        if( sz > 0 )
        {
          // if we could directly transfer a sorted array to an array1d including on device this wouldn't require storing a copy of the indices internally
          target->PackByIndex( buffer, m_sets_indices[ set_idx ], false, true );
        }
      }
      set_idx++;
    }
  }
  else
  {
    target->Pack( buffer, false, true );
  }
}

REGISTER_CATALOG_ENTRY( TaskBase, PackCollection, std::string const &, Group * const );
}
