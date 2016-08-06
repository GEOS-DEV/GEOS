/*
 * ConstitutiveManager.cpp
 *
 *  Created on: Aug 4, 2016
 *      Author: rrsettgast
 */

#include "ConstitutiveManager.hpp"

namespace geosx
{

using namespace dataRepository;

namespace constitutive
{


ConstitutiveManager::ConstitutiveManager( std::string const & name,
                                          WrapperCollection * const parent ):
  WrapperCollection(name,parent)
{
}

ConstitutiveManager::~ConstitutiveManager()
{}

void ConstitutiveManager::ReadXMLInput()
{
//  this->RegisterChildWrapperCollection<ConstitutiveBase>(keys::ConstitutiveBase);
  std::string newName = "matmodel";

  std::unique_ptr<ConstitutiveBase> temp = ConstitutiveBase::CatalogInterface::Factory("HypoElasticLinear", newName, this );
  auto& matmodel = this->RegisterChildWrapperCollection( newName, std::move(temp) );

}

}

} /* namespace geosx */
