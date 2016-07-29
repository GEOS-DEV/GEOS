/*
 * SidreWrapper.cpp
 *
 *  Created on: Jul 22, 2016
 *      Author: rrsettgast
 */

#include "SidreWrapper.hpp"

namespace geosx
{
namespace dataRepository
{

SidreWrapper::SidreWrapper()
{}

SidreWrapper::~SidreWrapper()
{}

asctoolkit::sidre::DataStore& SidreWrapper::dataStore()
{
  static asctoolkit::sidre::DataStore datastore;
  return datastore;
}

}
} /* namespace geosx */
