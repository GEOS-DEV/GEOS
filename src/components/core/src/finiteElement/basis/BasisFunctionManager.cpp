/*
 * BasisFunctionManager.cpp
 *
 *  Created on: Dec 5, 2017
 *      Author: sherman
 */

#include "BasisFunctionManager.hpp"
#include "BasisBase.hpp"

namespace geosx
{
using namespace dataRepository;

BasisFunctionManager::BasisFunctionManager( string const & name, ManagedGroup * const parent ):
  ManagedGroup(name,parent)
{}

BasisFunctionManager::~BasisFunctionManager()
{
  // TODO Auto-generated destructor stub
}

void BasisFunctionManager::FillDocumentationNode()
{
  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();
  docNode->setName("BasisFunctions");
  docNode->setSchemaType("UniqueNode");
}

void BasisFunctionManager::CreateChild( string const & childKey, string const & childName )
{
  std::unique_ptr<BasisBase> basis = BasisBase::CatalogInterface::Factory( childKey );
  this->RegisterViewWrapper( childName, std::move(basis) );
}


} /* namespace geosx */
