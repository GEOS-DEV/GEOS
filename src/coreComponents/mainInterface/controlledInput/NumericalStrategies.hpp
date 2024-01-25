#ifndef GEOS_NUMERICALSTRATEGIES_HPP
#define GEOS_NUMERICALSTRATEGIES_HPP

#include <yaml-cpp/yaml.h>

#include <pugixml.hpp>

#include <memory>

namespace geos::input::numericalStrategies
{

class NumericalStrategies
{
public:
  virtual ~NumericalStrategies() = default;

  virtual void fillProblemXmlNode( pugi::xml_node & problemNode ) const = 0;
};

void operator>>( const YAML::Node & node,
                 std::shared_ptr< NumericalStrategies > & numericalStrategy );

}

#endif //GEOS_NUMERICALSTRATEGIES_HPP
