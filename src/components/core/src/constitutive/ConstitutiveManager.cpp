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
                                          ManagedGroup * const parent ):
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

void ConstitutiveManager::CreateChild( string const & childKey, string const & childName )
{
  std::unique_ptr<ConstitutiveBase> material = ConstitutiveBase::CatalogInterface::Factory( childKey, childName, this );
  ConstitutiveBase * newMaterial = this->RegisterGroup<ConstitutiveBase>( childName, std::move(material) );
}

//ConstitutiveManager::constitutiveMaps & ConstitutiveManager::GetMaps( integer
// const reinit ) const
//{
//  static constitutiveMaps rval;
//  auto & map0 = rval.first;
//  auto & map1 = rval.second;
//
//  if( reinit==1 )
//  {
//    map0.clear();
//    map1.clear();
//    for( auto const & material : this->GetSubGroups() )
//    {
//      map0.push_back( &(*material.second) );
//      map1.insert({material.first,map0.size()-1});
//    }
//  }
//
//  return rval;
//}


}

} /* namespace geosx */
