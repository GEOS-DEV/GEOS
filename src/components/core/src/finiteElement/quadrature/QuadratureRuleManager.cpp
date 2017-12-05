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
  std::unique_ptr<QuadratureBase> quadrature = QuadratureBase::CatalogInterface::Factory( childKey );
  this->RegisterViewWrapper( childName, std::move(quadrature) );
}


} /* namespace geosx */
