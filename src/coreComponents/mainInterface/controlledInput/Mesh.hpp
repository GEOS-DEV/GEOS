#ifndef GEOS_INPUT_MESH_HPP
#define GEOS_INPUT_MESH_HPP

#include "common/DataTypes.hpp"

#include <pugixml.hpp>
#include <yaml-cpp/yaml.h>

#include <memory>

namespace geos::input::meshes
{

class Mesh
{
public:
  virtual ~Mesh() = default;

  virtual std::vector< string > getDomains() const = 0;

  virtual void fillProblemXmlNode( pugi::xml_node & problemNode ) const = 0;
};

void operator>>( const YAML::Node & node,
                 std::shared_ptr< Mesh > & mesh );

} // end of namespaces

#endif //GEOS_INPUT_MESH_HPP
