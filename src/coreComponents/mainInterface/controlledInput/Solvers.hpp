#ifndef GEOS_INPUT_SOLVERS_HPP
#define GEOS_INPUT_SOLVERS_HPP

#include "common/DataTypes.hpp"

#include <yaml-cpp/yaml.h>

#include <vector>

namespace geos::input::solvers
{

class Solver
{

};

class Variable
{
public:
  string functional_space;  // TOD Use an enum here?
  int fe_order;

};

class NamedVariable: public Variable
{
public:
  string name;
};

void operator>>( const YAML::Node & node,
                 std::vector< std::shared_ptr< Solver > > & solvers );

} // end of namespace

#endif //GEOS_INPUT_SOLVERS_HPP
