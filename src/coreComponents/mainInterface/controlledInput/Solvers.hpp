#ifndef GEOS_INPUT_SOLVERS_HPP
#define GEOS_INPUT_SOLVERS_HPP

#include <yaml-cpp/yaml.h>

#include <pugixml.hpp>

#include <memory>
#include <vector>

namespace geos::input::solvers
{

class Solver
{
public:
  virtual ~Solver() = default;

  virtual void fillProblemXmlNode( pugi::xml_node & problemNode ) const = 0;
};

void operator>>( const YAML::Node & node,
                 std::vector< std::shared_ptr< Solver > > & solvers );

} // end of namespace

#endif //GEOS_INPUT_SOLVERS_HPP
