/*
 * SidreWrapper.cpp
 *
 *  Created on: Jul 22, 2016
 *      Author: rrsettgast
 */

#include "SidreWrapper.hpp"

#ifdef USE_ATK
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

#endif
