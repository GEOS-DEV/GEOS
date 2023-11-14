#include "BufferedHistoryIO.hpp"

namespace geos
{

namespace detail
{
  inline size_t getTypeSize( std::type_index const & type )
  {
    if( type == std::type_index( typeid(char)) )
    {
      return sizeof( char );
    }
    else if( type == std::type_index( typeid(signed char)) )
    {
      return sizeof( signed char );
    }
    else if( type == std::type_index( typeid(real32)) )
    {
      return sizeof( real32 );
    }
    else if( type == std::type_index( typeid(real64)) )
    {
      return sizeof( real64 );
    }
    else if( type == std::type_index( typeid(integer)) )
    {
      return sizeof( integer );
    }
    else if( type == std::type_index( typeid(localIndex)) )
    {
      return sizeof( localIndex );
    }
    else if( type == std::type_index( typeid(globalIndex)) )
    {
      return sizeof( globalIndex );
    }
    else
    {
      return sizeof( char );
    }
  }
}

BufferedHistoryIO::BufferedHistoryIO( std::type_index typeIdx,
                                      localIndex rank,
                                      std::vector< localIndex > const & dims ) :
  m_bufferedCount( 0 ),
  m_bufferHead( nullptr ),
  m_dataBuffer( 0 ),
  m_bufferedLocalIdxCounts( ),
  m_typeIdx( typeIdx ),
  m_typeSize( detail::getTypeSize( typeIdx ) ),
  m_typeCount( 1 ),
  m_rank( rank ),
  m_dims( rank )
{
  for( size_t dd = 0; dd < m_rank; ++dd )
  {
    m_dims[dd] = LvArray::integerConversion< size_t >( dims[dd] );
    m_typeCount *= m_dims[dd];
  }
  m_dataBuffer.resize( m_typeSize * m_typeCount );
  m_bufferHead = &m_dataBuffer[0];
}

buffer_unit_type * BufferedHistoryIO::getBufferHead()
{
  resizeBuffer();
  m_bufferedCount++;
  buffer_unit_type * const currentBufferHead = m_bufferHead;
  m_bufferHead += getRowBytes();
  return currentBufferHead;
}

size_t BufferedHistoryIO::getRowBytes()
{
  return m_typeCount * m_typeSize;
}

void BufferedHistoryIO::emptyBuffer()
{
  m_sizeChanged = false;
  m_bufferedCount = 0;
  m_bufferHead = &m_dataBuffer[0];
  m_bufferedLocalIdxCounts.clear();
}

void BufferedHistoryIO::resizeBuffer()
{
  // need to store the count every time we collect (which calls this)
  //  regardless of whether the size changes
  m_bufferedLocalIdxCounts.emplace_back( m_dims[0] );

  size_t const capacity = m_dataBuffer.size();
  size_t const inUse = m_bufferHead - &m_dataBuffer[0];
  size_t const nextRow = getRowBytes( );
  // if needed, resize the buffer
  if( inUse + nextRow > capacity )
  {
    // resize based on the ammount currently in use rather than on the capacity ( less aggressive w/ changing sizes )
    m_dataBuffer.resize( inUse + ( nextRow * m_overallocMultiple ) );
  }
  // reset the buffer head and advance based on count in case the underlying data buffer moves during a resize
  m_bufferHead = &m_dataBuffer[0] + inUse;
}

void BufferedHistoryIO::updateCollectingCount( localIndex count )
{
  if( LvArray::integerConversion< size_t >( count ) != m_dims[0] )
  {
    m_sizeChanged = true;
    m_dims[0] = count;
    m_typeCount = count;
    for( size_t dd = 1; dd < m_rank; ++dd )
    {
      m_typeCount *= m_dims[dd];
    }
  }
}

}