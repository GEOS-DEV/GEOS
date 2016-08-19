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
                                    SynchronizedGroup * const parent ):
  SynchronizedGroup(name,parent)
{
}

ConstitutiveBase::~ConstitutiveBase()
{
}


void ConstitutiveBase::Registration( dataRepository::SynchronizedGroup * const )
{

}


ConstitutiveBase::CatalogInterface::CatalogType& ConstitutiveBase::GetCatalog()
{
  static ConstitutiveBase::CatalogInterface::CatalogType catalog;
  return catalog;
}


}
} /* namespace geosx */
