/*
 * Logger.hpp
 *
 *  Created on: Jul 17, 2017
 *      Author: settgast1
 */

#ifndef SRC_COMPONENTS_CORE_SRC_COMMON_LOGGER_HPP_
#define SRC_COMPONENTS_CORE_SRC_COMMON_LOGGER_HPP_

#include <string>
#include "common/GeosxConfig.hpp"

#ifdef USE_ATK
#include "slic/slic.hpp"
#include "slic/GenericOutputStream.hpp"
#endif


namespace geosx
{
void geos_abort( std::string message );

//#ifdef USE_ATK
//#define GEOS_ERROR(msg) SLIC_ERROR(msg)
//#else
#define GEOS_ERROR(msg) \
  geos_abort(msg)
//#endif

}

#endif /* SRC_COMPONENTS_CORE_SRC_COMMON_LOGGER_HPP_ */
