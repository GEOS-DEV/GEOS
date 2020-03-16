#ifndef GEOSX_HISTORY_IO_HPP_
#define GEOSX_HISTORY_IO_HPP_

#include "dataRepository/Group.hpp"

namespace geosx
{
  class DataSpec;

  class BufferedHistoryIO
  {
  public:

    BufferedHistoryIO(): m_buffered_count(0), m_data_buffer() {}

    virtual buffer_unit_type * GetBufferHead( DataSpec const * spec )
    {
      size_t data_size = spec->getTotalDataSize(); //bytes
      m_data_buffer.resize(m_data_buffer.size() + data_size);
      buffer_unit_type * row_head = &m_data_buffer[m_buffered_count*data_size];
      m_buffered_count++;
      return row_head;
    }

    virtual void Init( const string & m_write_target_name, DataSpec * spec, bool exists_okay ) = 0;
    virtual void Write( const string & m_write_target_name, DataSpec const * spec ) = 0;
    virtual void ClearAfter( const string & m_write_target_name, DataSpec const * spec, localIndex last_good ) = 0;

  protected:
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