/*
 * SidreWrapper.cpp
 *
 *  Created on: Jul 22, 2016
 *      Author: rrsettgast
 */
#if ATK_FOUND
#include "SidreWrapper.hpp"

namespace geosx
{
namespace dataRepository
{

SidreWrapper::SidreWrapper()
{}

SidreWrapper::~SidreWrapper()
{}

axom::sidre::DataStore& SidreWrapper::dataStore()
{
  static axom::sidre::DataStore datastore;
  return datastore;
}

} /* end namespace dataRepository */
} /* end namespace geosx */
#endif /* ATK_FOUND */