#include "SetIndexCollection.hpp"

namespace geosx
{
SetIndexCollection::SetIndexCollection( string const & object_path, string const & set_name, globalIndex set_index_offset ):
  HistoryCollection( "SetIndexCollection", nullptr ),
  m_objectPath( object_path ),
  m_setName( set_name ),
  m_setIndexOffset( set_index_offset )
{ }

  HistoryMetadata SetIndexCollection::getMetadata( ProblemManager & pm, localIndex const collectionIdx )
{
  GEOSX_UNUSED_VAR( collectionIdx );
  DomainPartition const & domain = *(pm.getDomainPartition( ));
  MeshLevel const & meshLevel = *domain.getMeshBody( 0 )->getMeshLevel( 0 );
  Group const * target_object = meshLevel.GetGroupByPath( m_objectPath );
  Group const * set_group = target_object->GetGroup( ObjectManagerBase::groupKeyStruct::setsString );
  dataRepository::Wrapper< SortedArray< localIndex > > const * const set_wrapper = set_group->getWrapper< SortedArray< localIndex > >( m_setName );
  HistoryMetadata meta = set_wrapper->getHistoryMetadata( );
  string outname = meta.getName( );
  meta.setName( outname + " Indices" );
  meta.setType( std::type_index( typeid(globalIndex)) );
  return meta;
}

void SetIndexCollection::collect( Group * domain_group,
                                  real64 const GEOSX_UNUSED_PARAM( time_n ),
                                  real64 const GEOSX_UNUSED_PARAM( dt ),
				  localIndex const collectionIdx,
                                  buffer_unit_type * & buffer )
{
  GEOSX_UNUSED_VAR( collectionIdx );
  DomainPartition const & domain = dynamicCast< DomainPartition const & >( *domain_group );
  MeshLevel const & meshLevel = *domain.getMeshBody( 0 )->getMeshLevel( 0 );
  Group const * target_object = meshLevel.GetGroupByPath( m_objectPath );
  Group const * set_group = target_object->GetGroup( ObjectManagerBase::groupKeyStruct::setsString );
  dataRepository::Wrapper< SortedArray< localIndex > > const * const set_wrapper = set_group->getWrapper< SortedArray< localIndex > >( m_setName );
  if( set_wrapper != nullptr )
  {
    SortedArrayView< localIndex const > const & set = set_wrapper->reference();
    globalIndex sz = set.size( );
    if( sz != 0 )
    {
      std::vector< globalIndex > meta_idx( sz );
      for( localIndex ii = 0; ii < sz; ++ii )
      {
        meta_idx[ii] = m_setIndexOffset + LvArray::integerConversion< globalIndex >( ii );
      }
      memcpy( buffer, &meta_idx[0], sz * sizeof( globalIndex ));
    }
  }
}

}
