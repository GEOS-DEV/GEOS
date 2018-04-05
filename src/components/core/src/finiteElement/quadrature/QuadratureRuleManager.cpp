/*
 * QuadratureRuleManager.cpp
 *
 *  Created on: Dec 5, 2017
 *      Author: sherman
 */

#include "QuadratureRuleManager.hpp"
#include "QuadratureBase.hpp"

namespace geosx
{
using namespace dataRepository;

QuadratureRuleManager::QuadratureRuleManager( string const & name, ManagedGroup * const parent ):
  ManagedGroup(name,parent)
{}

QuadratureRuleManager::~QuadratureRuleManager()
{
  // TODO Auto-generated destructor stub
}

void QuadratureRuleManager::FillDocumentationNode()
{
  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();
  docNode->setName("QuadratureRules");
  docNode->setSchemaType("UniqueNode");
}

void QuadratureRuleManager::CreateChild( string const & childKey, string const & childName )
{
  std::cout << "Quadrature Rule: " << childKey << ", " << childName << std::endl;
  std::unique_ptr<QuadratureBase> quadrature = QuadratureBase::CatalogInterface::Factory( childKey );
  this->RegisterViewWrapper( childName, std::move(quadrature) )->setRestartFlags(RestartFlags::NO_WRITE);
}

// Basis Base is not derived from ManagedGroup, so we need to do this manually:
void QuadratureRuleManager::ReadXMLsub( xmlWrapper::xmlNode const & targetNode )
{
  for (xmlWrapper::xmlNode childNode=targetNode.first_child() ; childNode ; childNode=childNode.next_sibling())
  {
    std::string childName = childNode.attribute("name").value();
    QuadratureBase * quadrature = this->getData<QuadratureBase>(childName);

    if (quadrature != nullptr)
    {
      quadrature->ReadXML(childNode);
    }
  }
}


} /* namespace geosx */
