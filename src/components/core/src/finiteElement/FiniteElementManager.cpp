/*
 * FiniteElementManager.cpp
 *
 *  Created on: Apr 18, 2017
 *      Author: rrsettgast
 */

#include "FiniteElementManager.hpp"
#include "basis/BasisBase.hpp"
#include "quadrature/QuadratureBase.hpp"

namespace geosx
{
using namespace dataRepository;

FiniteElementManager::FiniteElementManager( string const & name, ManagedGroup * const parent ):
  ManagedGroup(name,parent)
{
  this->RegisterGroup<ManagedGroup>(keys::basisFunctions);
  this->RegisterGroup<ManagedGroup>(keys::quadratureRules);
  this->RegisterGroup<ManagedGroup>(keys::finiteElementSpace);

}

FiniteElementManager::~FiniteElementManager()
{
  // TODO Auto-generated destructor stub
}

void FiniteElementManager::FillDocumentationNode()
{
  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();
  docNode->setName("NumericalMethods");
  docNode->setSchemaType("Node");
}

void FiniteElementManager::CreateChild( string const & childKey, string const & childName )
{
  if (strcmp(childKey.c_str(), keys::basisFunctions.c_str()) == 0)
  {
    // This will not place the pointer in a place where ReadXML will find:
    // Maybe place it on this?
    ManagedGroup * basisFunctions = this->GetGroup(keys::basisFunctions);
    std::unique_ptr<BasisBase> basis = BasisBase::CatalogInterface::Factory( childKey );
    basisFunctions->RegisterViewWrapper( childName, std::move(basis) );
  }
  else if (strcmp(childKey.c_str(), keys::quadratureRules.c_str()) == 0)
  {
    // This will not place the pointer in a place where ReadXML will find:
    // Maybe place it on this?
    ManagedGroup * quadratureRules = this->GetGroup(keys::quadratureRules);
    std::unique_ptr<QuadratureBase> quadrature = QuadratureBase::CatalogInterface::Factory( childKey );
    quadratureRules->RegisterViewWrapper( childName, std::move(quadrature) );
  }
  else
  {
    std::unique_ptr<ManagedGroup> fem = ManagedGroup::CatalogInterface::Factory( childKey, childName, this );
    ManagedGroup * feSpace = this->RegisterGroup( childName, std::move(fem) );
  }
}



} /* namespace geosx */
