#ifndef GEOS_INPUT_EVENTS_HPP
#define GEOS_INPUT_EVENTS_HPP

#include "common/DataTypes.hpp"

#include <pugixml.hpp>
#include <yaml-cpp/yaml.h>

namespace geos::input::events
{

class Event
{
public:
  explicit Event( string const & target )
    : m_target( target )
  { }

  virtual ~Event() = default;

  virtual void fillEventsXmlNode( pugi::xml_node & eventsNode ) const = 0;

protected:
  string m_target;
};

class PeriodicEvent : public Event
{
public:
  /**
   * @brief
   * @param target
   * @param every If empty, then `forceDt = 1.0` is added to the xml.
   */
  PeriodicEvent( string const & target,
                 string const & every = "" )
    : Event( target ),
      m_every( every )
  { }

  void fillEventsXmlNode( pugi::xml_node & eventsNode ) const override;

private:
  string m_every;
};

class SoloEvent : public Event
{
public:
  SoloEvent( string const & target,
             string const & at )
    : Event( target ),
      m_at( at )
  { }

  void fillEventsXmlNode( pugi::xml_node & eventsNode ) const override;

private:
  string m_at;
};

} // end of namespace

#endif //GEOS_INPUT_EVENTS_HPP
