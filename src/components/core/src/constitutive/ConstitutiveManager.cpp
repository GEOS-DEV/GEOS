/*
 * ConstitutiveManager.cpp
 *
 *  Created on: Aug 4, 2016
 *      Author: rrsettgast
 */

#include "ConstitutiveManager.hpp"

namespace geosx
{

ConstitutiveManager::ConstitutiveManager( std::string const & name,
                                          WrapperCollection * const parent ):
  WrapperCollection(name,parent)
{
}

ConstitutiveManager::~ConstitutiveManager()
{
  // TODO Auto-generated destructor stub
}

} /* namespace geosx */
