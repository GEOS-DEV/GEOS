#ifndef GEOS_INPUT_MESH_HPP
#define GEOS_INPUT_MESH_HPP

#include <pugixml.hpp>
#include <yaml-cpp/yaml.h>

#include <memory>

namespace geos::input::meshes
{

class Mesh
{
public:
  virtual ~Mesh() = default;

  virtual void fillMeshXmlNode( pugi::xml_node & meshNode ) const = 0;
};

void operator>>( const YAML::Node & node,
                 std::shared_ptr< Mesh > & mesh );

} // end of namespaces

#endif //GEOS_INPUT_MESH_HPP
