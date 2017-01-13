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

ConstitutiveBase::ConstitutiveBase( std::string const & name )
{}

ConstitutiveBase::~ConstitutiveBase()
{}




ConstitutiveBase::CatalogInterface::CatalogType& ConstitutiveBase::GetCatalog()
{
  static ConstitutiveBase::CatalogInterface::CatalogType catalog;
  return catalog;
}


}
} /* namespace geosx */
