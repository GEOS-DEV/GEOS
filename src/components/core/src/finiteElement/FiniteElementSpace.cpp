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



FiniteElementSpace::FiniteElementSpace( std::string const & name, SynchronizedGroup * const parent ):
    SynchronizedGroup(name,parent)
{}

FiniteElementSpace::~FiniteElementSpace()
{}


void FiniteElementSpace::Registration( dataRepository::SynchronizedGroup * const parent )
{
  parent->RegisterGroup<SynchronizedGroup>(keys::FEM_Nodes);
  parent->RegisterGroup<SynchronizedGroup>(keys::FEM_Elements);
}



SynchronizedGroup & FiniteElementSpace::getNodeManager()
{
  return GetGroup<DomainPartition>(keys::FEM_Nodes);
}

SynchronizedGroup & FiniteElementSpace::getEdgeManager()
{
  return GetGroup<SynchronizedGroup>(keys::FEM_Nodes);
}

SynchronizedGroup & FiniteElementSpace::getFaceManager()
{
  return GetGroup<SynchronizedGroup>(keys::FEM_Faces);
}

SynchronizedGroup & FiniteElementSpace::getElementManager()
{
  return GetGroup<SynchronizedGroup>(keys::FEM_Elements);
}


} /* namespace geosx */
