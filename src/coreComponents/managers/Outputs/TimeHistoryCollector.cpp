#include "TimeHistoryCollector.hpp"

namespace geosx
{

  TimeHistoryCollector::TimeHistoryCollector( string const & name, Group * parent )
     : ExecutableGroup( name, parent )
     , m_data_title()
     , m_time_collector(new ScalarCollector<real64>("TimeCollector"))
  {
    registerWrapper(viewKeysStruct::dataTitleKey, &m_data_title, false)->
      setApplyDefaultValue("Time")->
      setInputFlag(InputFlags::REQUIRED)->
      setDescription("A title for the data collected by this collector.");

    //registerGroup()

   // add m_time_collector to the group
  }
  virtual void TimeHistoryCollector::AddToSpec( DataSpec & data_spec )
  {
    for( auto & collector : GetSubGroups<TimeHistoryCollector&>() )
    {
      collector.AddToSpec( data_spec );
    }
  }
  virtual void TimeHistoryCollector::Collect( real64 const time_n, real64 const dt, buffer_unit_type *& buffer )
  {
    for( auto & collector : GetSubGroups<TimeHistoryCollector&>() )
    {
      collector.Collect(time_n,dt,buffer);
    }
  }

}