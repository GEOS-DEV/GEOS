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

void FiniteElementManager::FillDocumentationNode( dataRepository::ManagedGroup * const group )
{
  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();
  docNode->setName("NumericalMethods");
  docNode->setSchemaType("Node");
}

void FiniteElementManager::ReadXMLsub( xmlWrapper::xmlNode const & node )
{
  std::cout << "Reading Components for Numerical Methods:" << std::endl;
  if ( node == nullptr )
  {
    throw std::invalid_argument("Numerical Methods block not present in input xml file!");
  }
  else
  {
    xmlWrapper::xmlNode basisNode = node.child(keys::basisFunctions.c_str());
    if( basisNode != nullptr )
    {
      ManagedGroup * basisFunctions = this->GetGroup(keys::basisFunctions);

      for (xmlWrapper::xmlNode childNode=basisNode.first_child(); childNode; childNode=childNode.next_sibling())
      {
        string catalogName = childNode.name();
        string name = childNode.attribute("name").value();
        std::cout <<childNode.name()<<", "<<childNode.attribute("name").value()<< std::endl;

        std::unique_ptr<BasisBase> basis = BasisBase::CatalogInterface::Factory( catalogName );
        basis->ReadXML( childNode );
        basisFunctions->RegisterViewWrapper( name, std::move(basis) );
      }
    }

    xmlWrapper::xmlNode quadratureNode = node.child(keys::quadratureRules.c_str());
    if( quadratureNode != nullptr )
    {
      ManagedGroup * quadratureRules = this->GetGroup(keys::quadratureRules);

      for (xmlWrapper::xmlNode childNode=quadratureNode.first_child(); childNode; childNode=childNode.next_sibling())
      {
        string catalogName = childNode.name();
        string name = childNode.attribute("name").value();
        std::cout <<childNode.name()<<", "<<childNode.attribute("name").value()<< std::endl;

        std::unique_ptr<QuadratureBase> quadrature = QuadratureBase::CatalogInterface::Factory( catalogName );
        quadrature->ReadXML(childNode);
        quadratureRules->RegisterViewWrapper( name, std::move(quadrature) );
      }

    }

    xmlWrapper::xmlNode finiteElementNode = node.child(keys::finiteElements.c_str());
    if( finiteElementNode != nullptr )
    {
//      ManagedGroup & feSpaces = RegisterGroup(keys::FE_Space);
      for (xmlWrapper::xmlNode childNode=finiteElementNode.first_child(); childNode; childNode=childNode.next_sibling())
      {
        string catalogName = childNode.name();
        string name = childNode.attribute("name").value();
        std::cout <<childNode.name()<<", "<<childNode.attribute("name").value()<< std::endl;

        std::unique_ptr<ManagedGroup> fem = ManagedGroup::CatalogInterface::Factory( catalogName, name, this );
        fem->SetDocumentationNodes(nullptr);
        fem->RegisterDocumentationNodes();
        ManagedGroup & feSpace = this->RegisterGroup( name, std::move(fem) );
        feSpace.ReadXML(childNode);

      }

    }

  }
}



} /* namespace geosx */
