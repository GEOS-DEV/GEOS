/*
 * FiniteElementSpace.cpp
 *
 *  Created on: Aug 4, 2016
 *      Author: rrsettgast
 */

#include "FiniteElementSpace.hpp"
#include "managers/DomainPartition.hpp"

namespace geosx
{
using namespace dataRepository;



FiniteElementSpace::FiniteElementSpace( std::string const & name, ManagedGroup * const parent ) :
  ManagedGroup(name,parent)
{}

FiniteElementSpace::~FiniteElementSpace()
{}


void FiniteElementSpace::Registration( dataRepository::ManagedGroup * const parent )
{
  parent->RegisterGroup<ManagedGroup>(keys::FEM_Nodes);
  parent->RegisterGroup<ManagedGroup>(keys::FEM_Elements);
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


} /* namespace geosx */
