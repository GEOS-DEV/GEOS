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

void ConstitutiveManager::FillDocumentationNode( dataRepository::ManagedGroup * const group )
{
  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();
  docNode->setName("Constitutive");
  docNode->setSchemaType("Node");
}

void ConstitutiveManager::ReadXMLsub( xmlWrapper::xmlNode const & targetNode )
{
  for (xmlWrapper::xmlNode childNode=targetNode.first_child(); childNode; childNode=childNode.next_sibling())
  {
      std::string materialName = childNode.attribute("name").value();
      std::string materialKey = childNode.name();
//      std::cout<<materialName<<std::endl;
      std::unique_ptr<ConstitutiveBase> material = ConstitutiveBase::CatalogInterface::Factory( materialKey, materialName, this );
      ConstitutiveBase & newMaterial = this->RegisterGroup<ConstitutiveBase>( materialName, std::move(material) );
      newMaterial.SetDocumentationNodes( nullptr );
//      newMaterial.RegisterDocumentationNodes();
      newMaterial.ReadXML( childNode );
  }
//  std::cout<<this->GetSubGroups().size()<<std::endl;
//  for( auto const & material : this->GetSubGroups() )
//  {
//    std::cout<<material.first<<std::endl;
//  }
}

ConstitutiveManager::constitutiveMaps & ConstitutiveManager::GetMaps( int32 const reinit ) const
{
//  static array<ManagedGroup *> map0;
//  static map<string,int32> map1;
  static constitutiveMaps rval;
  auto & map0 = rval.first;
  auto & map1 = rval.second;

  if( reinit==1 )
  {
    map0.clear();
    map1.clear();
    for( auto & material : this->GetSubGroups() )
    {
      map0.push_back( material.second );
      map1.insert({material.first,map0.size()-1});
    }
  }

  return rval;
}


}

} /* namespace geosx */
