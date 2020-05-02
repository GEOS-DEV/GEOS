#include "BufferAllocator.hpp"
#include "common/DataTypes.hpp"

#ifdef GEOSX_USE_CHAI
namespace geosx
{

bool prefer_pinned_buffer = true;

void setPreferPinned( bool p )
{
  prefer_pinned_buffer = p;
}

bool getPreferPinned( )
{
  return prefer_pinned_buffer;
}

template< typename T >
buffer_allocator< T >::buffer_allocator()
  : m_alloc( umpire::TypedAllocator< T >( umpire::ResourceManager::getInstance().getAllocator( umpire::resource::Host )))
  , m_prefer_pinned_l( prefer_pinned_buffer )
{
  auto & rm = umpire::ResourceManager::getInstance();
  if( rm.isAllocator( "PINNED" ) && m_prefer_pinned_l )
    m_alloc = umpire::TypedAllocator< T >( rm.getAllocator( umpire::resource::Pinned ));
}

template< typename T >
T * buffer_allocator< T >::allocate( size_t sz ) { return m_alloc.allocate( sz ); }

template< typename T >
void buffer_allocator< T >::deallocate( T * buffer, size_t sz )
{
  if( buffer != nullptr )
    m_alloc.deallocate( buffer, sz );
}

template< typename T >
bool buffer_allocator< T >::operator!=( const buffer_allocator< T > & )
{
  return false;
}

template< typename T >
bool buffer_allocator< T >::operator==( const buffer_allocator< T > & other )
{
  return !operator!=( other );
}

template class buffer_allocator< buffer_unit_type >;

}
#endif
