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
                                          SynchronizedGroup * const parent ):
  SynchronizedGroup(name,parent)
{
}

ConstitutiveManager::~ConstitutiveManager()
{}

void ConstitutiveManager::ReadXMLInput()
{
//  this->RegisterChildWrapperCollection<ConstitutiveBase>(keys::ConstitutiveBase);
  std::string newName = "matmodel";

//  RegisterGroup( newName,
  //               ConstitutiveBase::CatalogInterface::Factory("HypoElasticLinear", newName, this ) );

  RegisterGroup<ConstitutiveBase>(newName,"HypoElasticLinear");
}

}

} /* namespace geosx */
