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
    // BufferedHistoryIO(size_t init_size) : m_buffered_count(0), m_data_buffer(init_size) {}

    buffer_unit_type * GetBufferHead( )
    {
      size_t osize = m_data_buffer.size();
      resizeBuffer();
      m_buffered_count++;
      return &m_data_buffer[osize];
    }

    virtual void Init( bool exists_okay ) = 0;
    virtual void Write( ) = 0;

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