/*
 * ConstitutiveBase.cpp
 *
 *  Created on: Jul 28, 2016
 *      Author: rrsettgast
 */

#include "ConstitutiveBase.hpp"

namespace geosx
{
using namespace dataRepository;
namespace constitutive
{

ConstitutiveBase::ConstitutiveBase( std::string const & name,
                                    ManagedGroup * const parent ):
    ManagedGroup(name,parent)
{
  this->RegisterGroup<ManagedGroup>("StateData");
  this->RegisterGroup<ManagedGroup>("ParameterData");
}

ConstitutiveBase::~ConstitutiveBase()
{}




ConstitutiveBase::CatalogInterface::CatalogType& ConstitutiveBase::GetCatalog()
{
  static ConstitutiveBase::CatalogInterface::CatalogType catalog;
  return catalog;
}


}
} /* namespace geosx */
