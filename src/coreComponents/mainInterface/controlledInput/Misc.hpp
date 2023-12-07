#ifndef GEOS_INPUT_MISC_HPP
#define GEOS_INPUT_MISC_HPP

#include "codingUtilities/StringUtilities.hpp"

#include "common/DataTypes.hpp"

#include <yaml-cpp/yaml.h>

namespace geos::input
{

real64 convertTime( string const & time );

string convertYamlElementTypeToGeosElementType( string const yamlElementType );

template< class T >
string createGeosArray( T const & t )
{
  return "{ " + stringutilities::join( t, ", " ) + " }";
}

} // end of namespace

#endif //GEOS_INPUT_MISC_HPP
