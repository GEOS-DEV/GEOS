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
                                          ManagedGroup * const parent ) :
  ManagedGroup(name,parent)
{}

ConstitutiveManager::~ConstitutiveManager()
{}

void ConstitutiveManager::ReadXMLsub( pugi::xml_node const & targetNode )
{
  for (pugi::xml_node childNode=targetNode.first_child(); childNode; childNode=childNode.next_sibling())
  {
      std::string materialName = childNode.attribute("name").value();
      std::string materialKey = childNode.name();
      std::cout<<materialName<<std::endl;
      std::unique_ptr<ConstitutiveBase> material = ConstitutiveBase::CatalogInterface::Factory( materialKey, materialName, this );
      ConstitutiveBase & newMaterial = this->RegisterGroup<ConstitutiveBase>( materialName, std::move(material) );
      newMaterial.SetDocumentationNodes( nullptr );
      newMaterial.ReadXML( childNode );
  }
}

}

} /* namespace geosx */
