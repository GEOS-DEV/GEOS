#include "TimeHistoryCollector.hpp"

namespace geosx
{

  TimeHistoryCollector::TimeHistoryCollector( string const & name, Group * parent )
     : ExecutableGroup( name, parent )
     , m_data_title()
     , m_buffer_callback( )
  {
    registerWrapper(viewKeysStruct::dataTitleKey, &m_data_title, false)->
      setApplyDefaultValue("Time")->
      setInputFlag(InputFlags::REQUIRED)->
      setDescription("A title for the data collected by this collector.");
  }
}