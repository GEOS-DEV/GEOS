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

void operator>>( const YAML::Node & node,
                 std::vector< std::shared_ptr< Solver > > & solvers );

} // end of namespace

#endif //GEOS_INPUT_SOLVERS_HPP
