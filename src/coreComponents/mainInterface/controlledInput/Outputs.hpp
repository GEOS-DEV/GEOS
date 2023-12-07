#ifndef GEOS_INPUT_OUTPUTS_HPP
#define GEOS_INPUT_OUTPUTS_HPP

#include "Events.hpp"
#include "Misc.hpp"

#include "common/DataTypes.hpp"

#include <pugixml.hpp>
#include <yaml-cpp/yaml.h>

#include <memory>

namespace geos::input::outputs
{

class Output
{
public:
  Output( string name,
          string const & every,
          std::vector< string > const & at )
    : m_name( name ),
      m_every( every ),
      m_at( at )
  { }

  virtual ~Output() = default;

  virtual void fillOutputsXmlNode( pugi::xml_node & outputsNode ) const = 0;

  std::vector< std::shared_ptr< events::Event > > getEvents() const
  {
    std::vector< std::shared_ptr< events::Event > > result;

    if( !m_every.empty() )
    {
      result.emplace_back( std::make_shared< events::PeriodicEvent >( "/Outputs/" + m_name, m_every ) );
    }

    for( string const & at: m_at )
    {
      if( !at.empty() )
      {
        result.emplace_back( std::make_shared< events::SoloEvent >( "/Outputs/" + m_name, at ) );
      }
    }

    return result;
  };

protected:
  string m_name;

private:
  string m_every;
  std::vector< string > m_at;
};

void operator>>( const YAML::Node & node,
                 std::vector< std::shared_ptr< Output > > & outputs );


} // end of namespace

#endif //GEOS_INPUT_OUTPUTS_HPP
