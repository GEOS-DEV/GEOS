/*
 * FiniteElementManager.cpp
 *
 *  Created on: Apr 18, 2017
 *      Author: rrsettgast
 */

#include "FiniteElementManager.hpp"
#include "basis/BasisFunctionManager.hpp"
#include "quadrature/QuadratureRuleManager.hpp"
#include "FiniteElementSpaceManager.hpp"


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
}



} /* namespace geosx */
