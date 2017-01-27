/*
 * FiniteElementSpace.cpp
 *
 *  Created on: Aug 4, 2016
 *      Author: rrsettgast
 */

#include "FiniteElementSpace.hpp"
#include "managers/DomainPartition.hpp"
#include "managers/ObjectManagerBase.hpp"
#include "managers/NodeManager.hpp"

namespace geosx
{
using namespace dataRepository;



FiniteElementSpace::FiniteElementSpace( std::string const & name, ManagedGroup * const parent ) :
  ManagedGroup(name,parent)
{}

FiniteElementSpace::~FiniteElementSpace()
{}


void FiniteElementSpace::BuildDataStructure( dataRepository::ManagedGroup * const parent )
{
  m_nodeManager    = dynamic_cast<NodeManager*>(&(parent->RegisterGroup<ObjectManagerBase>( keys::FEM_Nodes, NodeManager::CatalogName() ) ) );
  m_elementManager = dynamic_cast<ElementManager*>(&(parent->RegisterGroup<ObjectManagerBase>(keys::FEM_Elements, ElementManager::CatalogName())));



}



ManagedGroup & FiniteElementSpace::getNodeManager()
{
  return GetGroup<DomainPartition>(keys::FEM_Nodes);
}

ManagedGroup & FiniteElementSpace::getEdgeManager()
{
  return GetGroup<ManagedGroup>(keys::FEM_Nodes);
}

ManagedGroup & FiniteElementSpace::getFaceManager()
{
  return GetGroup<ManagedGroup>(keys::FEM_Faces);
}

ManagedGroup & FiniteElementSpace::getElementManager()
{
  return GetGroup<ManagedGroup>(keys::FEM_Elements);
}

REGISTER_CATALOG_ENTRY( ManagedGroup, FiniteElementSpace, std::string const &, ManagedGroup * const )

} /* namespace geosx */
