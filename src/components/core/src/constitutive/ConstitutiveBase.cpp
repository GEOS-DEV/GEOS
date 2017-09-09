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
  this->RegisterGroup<ManagedGroup>(keys::parameterData);
  this->RegisterGroup<ManagedGroup>(keys::stateData);
}

ConstitutiveBase::~ConstitutiveBase()
{}




ConstitutiveBase::CatalogInterface::CatalogType& ConstitutiveBase::GetCatalog()
{
  static ConstitutiveBase::CatalogInterface::CatalogType catalog;
  return catalog;
}

void ConstitutiveBase::resize( localIndex newsize )
{
  ManagedGroup::resize(newsize);
  this->GetGroup(keys::parameterData)->resize(newsize);
  this->GetGroup(keys::stateData)->resize(newsize);
}



}
} /* namespace geosx */
