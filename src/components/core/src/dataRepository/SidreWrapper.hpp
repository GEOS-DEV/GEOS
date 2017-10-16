/*
 * SidreWrapper.hpp
 *
 *  Created on: Jul 22, 2016
 *      Author: rrsettgast
 */

#ifndef COMPONENTS_CORE_SRC_DATAREPOSITORY_SIDREWRAPPER_HPP_
#define COMPONENTS_CORE_SRC_DATAREPOSITORY_SIDREWRAPPER_HPP_
#include "common/GeosxConfig.hpp"

#include "sidre/DataStore.hpp"

namespace geosx
{
namespace dataRepository
{

class SidreWrapper
{
public:
  SidreWrapper();
  ~SidreWrapper();

  static axom::sidre::DataStore& dataStore();

private:

};

}
} /* namespace geosx */

#endif /* COMPONENTS_CORE_SRC_DATAREPOSITORY_SIDREWRAPPER_HPP_ */
