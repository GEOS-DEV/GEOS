#ifndef GEOS_NUMERICALSTRATEGIES_HPP
#define GEOS_NUMERICALSTRATEGIES_HPP

#include <yaml-cpp/yaml.h>

#include <pugixml.hpp>

namespace geos::input::numericalStrategies
{

class NumericalStrategies
{
public:
  void fillProblemXmlNode( pugi::xml_node & problemNode ) const;
};

void operator>>( const YAML::Node & node,
                 NumericalStrategies & numericalStrategy );

}

#endif //GEOS_NUMERICALSTRATEGIES_HPP
