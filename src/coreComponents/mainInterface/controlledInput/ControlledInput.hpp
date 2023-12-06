#ifndef GEOS_CONTROLLEDINPUT_HPP
#define GEOS_CONTROLLEDINPUT_HPP

#include "dataRepository/xmlWrapper.hpp"
#include "common/DataTypes.hpp"

namespace geos::input {

/**
 * @brief Converts the stable input file to its xml internal representation
 * @param stableInputFileName The yaml input file name.
 * @param doc The output xml document being fed.
 */
void convert( string const & stableInputFileName, xmlWrapper::xmlDocument & doc );

}




#endif //GEOS_CONTROLLEDINPUT_HPP
