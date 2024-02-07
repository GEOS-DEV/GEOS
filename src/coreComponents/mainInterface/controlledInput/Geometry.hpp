#ifndef GEOS_GEOMETRY_HPP
#define GEOS_GEOMETRY_HPP

#include <yaml-cpp/yaml.h>

#include <pugixml.hpp>

#include <memory>
#include <vector>


namespace geos::input::geometry
{

class Geometry
{
public:

  virtual ~Geometry() = default;

  virtual void fillProblemXmlNode( pugi::xml_node & problemNode ) const = 0;
};

void operator>>( const YAML::Node & node,
                 std::vector< std::shared_ptr< Geometry > > & geometries );

}

#endif //GEOS_GEOMETRY_HPP
