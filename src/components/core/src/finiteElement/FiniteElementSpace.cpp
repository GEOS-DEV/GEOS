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



FiniteElementSpace::FiniteElementSpace( std::string const & name, WrapperCollection * const parent ):
    WrapperCollection(name,parent)
{}

FiniteElementSpace::~FiniteElementSpace()
{}


void FiniteElementSpace::Registration( dataRepository::WrapperCollection * const parent )
{
  parent->RegisterChildWrapperCollection<WrapperCollection>(keys::FEM_Nodes);
  parent->RegisterChildWrapperCollection<WrapperCollection>(keys::FEM_Elements);
}



WrapperCollection & FiniteElementSpace::getNodeManager()
{
  return GetChildWrapperCollection<DomainPartition>(keys::FEM_Nodes);
}

WrapperCollection & FiniteElementSpace::getEdgeManager()
{
  return GetChildWrapperCollection<WrapperCollection>(keys::FEM_Nodes);
}

WrapperCollection & FiniteElementSpace::getFaceManager()
{
  return GetChildWrapperCollection<WrapperCollection>(keys::FEM_Faces);
}

WrapperCollection & FiniteElementSpace::getElementManager()
{
  return GetChildWrapperCollection<WrapperCollection>(keys::FEM_Elements);
}


} /* namespace geosx */
