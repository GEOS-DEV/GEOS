/*
 * SidreWrapper.hpp
 *
 *  Created on: Jul 22, 2016
 *      Author: rrsettgast
 */

#ifndef COMPONENTS_CORE_SRC_DATAREPOSITORY_SIDREWRAPPER_HPP_
#define COMPONENTS_CORE_SRC_DATAREPOSITORY_SIDREWRAPPER_HPP_

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

  static asctoolkit::sidre::DataStore& dataStore()
  {
    static asctoolkit::sidre::DataStore datastore;
    return datastore;
  }

private:

};

}
} /* namespace geosx */

#endif /* COMPONENTS_CORE_SRC_DATAREPOSITORY_SIDREWRAPPER_HPP_ */
