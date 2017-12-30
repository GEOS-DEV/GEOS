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
  this->RegisterGroup<BasisFunctionManager>(keys::basisFunctions);
  this->RegisterGroup<QuadratureRuleManager>(keys::quadratureRules);
  this->RegisterGroup<FiniteElementSpaceManager>(keys::finiteElementSpaces);
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
