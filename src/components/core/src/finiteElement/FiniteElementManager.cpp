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

void FiniteElementManager::ReadXMLsub( pugi::xml_node const & node )
{
  std::cout << "Reading Components for Numerical Methods:" << std::endl;
  if ( node == nullptr )
  {
    throw std::invalid_argument("Numerical Methods block not present in input xml file!");
  }
  else
  {
    pugi::xml_node basisNode = node.child(keys::basisFunctions.c_str());
    if( basisNode != nullptr )
    {
      ManagedGroup & basisFunctions = this->GetGroup(keys::basisFunctions);

      for (pugi::xml_node childNode=basisNode.first_child(); childNode; childNode=childNode.next_sibling())
      {
        string catalogName = childNode.name();
        string name = childNode.attribute("name").value();
        std::cout <<childNode.name()<<", "<<childNode.attribute("name").value()<< std::endl;

        std::unique_ptr<BasisBase> basis = BasisBase::CatalogInterface::Factory( catalogName );
        basis->ReadXML( childNode );
        basisFunctions.RegisterViewWrapper( name, std::move(basis) );
      }
    }

    pugi::xml_node quadratureNode = node.child(keys::quadratureRules.c_str());
    if( quadratureNode != nullptr )
    {
      ManagedGroup & quadratureRules = this->GetGroup(keys::quadratureRules);

      for (pugi::xml_node childNode=quadratureNode.first_child(); childNode; childNode=childNode.next_sibling())
      {
        string catalogName = childNode.name();
        string name = childNode.attribute("name").value();
        std::cout <<childNode.name()<<", "<<childNode.attribute("name").value()<< std::endl;

        std::unique_ptr<QuadratureBase> quadrature = QuadratureBase::CatalogInterface::Factory( catalogName );
        quadrature->ReadXML(childNode);
        quadratureRules.RegisterViewWrapper( name, std::move(quadrature) );
      }

    }

  }
}


} /* namespace geosx */
