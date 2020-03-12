#ifndef GEOSX_HISTORY_IO_HPP_
#define GEOSX_HISTORY_IO_HPP_

#include "dataRepository/Group.hpp"

namespace geosx
{
  class DataSpec;

  class BufferedHistoryIO : public Group
  {
  public:
    BufferedHistoryIO( const string & name, Group * group )
    : Group(name,group)
    { }
    virtual void Init( const string & m_write_target_name, DataSpec * spec ) = 0;
    virtual buffer_unit_type * GetRowHead( DataSpec const * spec ) = 0;
    virtual void Write( const string & m_write_target_name, DataSpec const * spec ) = 0;
    virtual void ClearAfterHead( const string & m_write_target_name, DataSpec const * spec ) = 0;
  };

}
#endif