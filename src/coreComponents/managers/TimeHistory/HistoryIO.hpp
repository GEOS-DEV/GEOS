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
      return &m_data_buffer[osize];
    }

    virtual void Init( bool exists_okay ) = 0;
    virtual void Write( ) = 0;

    localIndex GetBufferedCount( ) { return m_buffered_count; }
  protected:
    virtual void resizeBuffer( ) = 0;

    void EmptyBuffer( bool dealloc = false )
    {
      m_buffered_count = 0;
      if( dealloc )
      {
        m_data_buffer.resize(0);
      }
    }

    localIndex m_buffered_count;
    array1d<buffer_unit_type> m_data_buffer;
  };

}
#endif
