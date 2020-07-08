#include "SetIndexCollection.hpp"

namespace geosx
{
SetIndexCollection::SetIndexCollection( string const & object_path, string const & set_name, globalIndex set_index_offset ):
  HistoryCollection( "SetIndexCollection", nullptr ),
  m_object_path( object_path ),
  m_set_name( set_name ),
  m_set_index_offset( set_index_offset )
{ }

HistoryMetadata SetIndexCollection::GetMetadata( Group * problem_group )
{
  ProblemManager const * pm = dynamicCast< ProblemManager const * >( problem_group );
  DomainPartition const & domain = *pm->getDomainPartition( );
  MeshLevel const & meshLevel = *domain.getMeshBody( 0 )->getMeshLevel( 0 );
  Group const * target_object = meshLevel.GetGroupByPath( m_object_path );
  Group const * set_group = target_object->GetGroup( ObjectManagerBase::groupKeyStruct::setsString );
  dataRepository::Wrapper< SortedArray< localIndex > > const * const set_wrapper = set_group->getWrapper< SortedArray< localIndex > >( m_set_name );
  HistoryMetadata meta = set_wrapper->getBufferedIOMetadata( );
  string outname = meta.getName( );
  meta.setName( outname + " Indices" );
  meta.setType( std::type_index( typeid(globalIndex)) );
  return meta;
}

void SetIndexCollection::Collect( Group * domain_group,
                                  real64 const GEOSX_UNUSED_PARAM( time_n ),
                                  real64 const GEOSX_UNUSED_PARAM( dt ),
                                  buffer_unit_type * & buffer )
{
  DomainPartition const & domain = dynamicCast< DomainPartition const & >( *domain_group );
  MeshLevel const & meshLevel = *domain.getMeshBody( 0 )->getMeshLevel( 0 );
  Group const * target_object = meshLevel.GetGroupByPath( m_object_path );
  Group const * set_group = target_object->GetGroup( ObjectManagerBase::groupKeyStruct::setsString );
  dataRepository::Wrapper< SortedArray< localIndex > > const * const set_wrapper = set_group->getWrapper< SortedArray< localIndex > >( m_set_name );
  if( set_wrapper != nullptr )
  {
    SortedArrayView< localIndex const > const & set = set_wrapper->reference();
    globalIndex sz = set.size( );
    if( sz != 0 )
    {
      // this might wind up having to move the meta_idx array back and forth to device a lot
      array1d< globalIndex > meta_idx( sz );
      for( localIndex ii = 0; ii < sz; ++ii )
      {
        meta_idx[ii] = m_set_index_offset + LvArray::integerConversion< globalIndex >( ii );
      }
      bufferOps::PackDataDevice< true >( buffer, meta_idx.toViewConst() );
    }
  }
}

}
