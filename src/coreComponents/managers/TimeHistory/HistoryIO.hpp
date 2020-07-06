#ifndef GEOSX_HISTORY_IO_HPP_
#define GEOSX_HISTORY_IO_HPP_

#include "dataRepository/Group.hpp"

namespace geosx
{
  class DataSpec;

  class BufferedHistoryIO
  {
  public:

    BufferedHistoryIO(): m_buffered_count(0), m_data_buffer(0) {}

    virtual ~BufferedHistoryIO() {}

    buffer_unit_type * GetBufferHead( )
    {
      size_t osize = m_data_buffer.size();
      resizeBuffer();
      m_buffered_count++;
      buffer_unit_type * head = nullptr;
      if ( m_data_buffer.size( ) )
      {
        head = &m_data_buffer[osize];
      }
      return head;
    }

    virtual void Init( bool exists_okay, bool once ) = 0;
    virtual void Write( ) = 0;
    virtual void CompressInFile( ) = 0;

    virtual globalIndex GetRankOffset( ) = 0;

    localIndex GetBufferedCount( ) { return m_buffered_count; }
  protected:
    virtual void resizeBuffer( ) = 0;

    void EmptyBuffer( )
    {
      m_buffered_count = 0;
      m_data_buffer = 0;
      m_data_buffer.resize( 0 );
    }

    localIndex m_buffered_count;
    array1d<buffer_unit_type> m_data_buffer;
  };

}
#endif
