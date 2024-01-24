#include "Events.hpp"

#include "Misc.hpp"

namespace geos::input::events
{

using namespace pugi;

void PeriodicEvent::fillEventsXmlNode( pugi::xml_node & eventsNode ) const
{
  pugi::xml_node periodicEvent = eventsNode.append_child( "PeriodicEvent" );
  periodicEvent.append_attribute( "target" ) = m_target.c_str();
  if( m_every.empty() )  // TODO Warning, this is a hack!
  {
    periodicEvent.append_attribute( "forceDt" ) = "1.0";
  }
  else
  {
    periodicEvent.append_attribute( "timeFrequency" ) = convertTime( m_every );
  }
}

void SoloEvent::fillEventsXmlNode( pugi::xml_node & eventsNode ) const
{
  pugi::xml_node soloEvent = eventsNode.append_child( "SoloEvent" );
  soloEvent.append_attribute( "targetTime" ) = convertTime( m_at );
  soloEvent.append_attribute( "target" ) = m_target.c_str();
}

} // end of namespace
