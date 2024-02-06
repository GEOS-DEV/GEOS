#ifndef GEOS_FUNCTIONS_HPP
#define GEOS_FUNCTIONS_HPP

#include <yaml-cpp/yaml.h>

#include <pugixml.hpp>

#include <memory>


namespace geos::input::functions
{

class Function
{
public:
  virtual ~Function() = default;

  virtual void fillProblemXmlNode( pugi::xml_node & problemNode ) const = 0;
};

void operator>>( const YAML::Node & node,
                 std::vector< std::shared_ptr< Function > > & functions );

}

#endif //GEOS_FUNCTIONS_HPP
