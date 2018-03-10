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
  std::cout << "Basis Function: " << childKey << ", " << childName << std::endl;
  std::unique_ptr<BasisBase> basis = BasisBase::CatalogInterface::Factory( childKey );
  this->RegisterViewWrapper( childName, std::move(basis) )->setWriteToRestart(false);
}

// Basis Base is not derived from ManagedGroup, so we need to do this manually:
void BasisFunctionManager::ReadXMLsub( xmlWrapper::xmlNode const & targetNode )
{
  for (xmlWrapper::xmlNode childNode=targetNode.first_child() ; childNode ; childNode=childNode.next_sibling())
  {
    std::string childName = childNode.attribute("name").value();
    BasisBase * basis = this->getData<BasisBase>(childName);

    if (basis != nullptr)
    {
      basis->ReadXML(childNode);
    }
  }
}


} /* namespace geosx */
