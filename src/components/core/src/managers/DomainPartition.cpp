/*
 * DomainPartition.cpp
 *
 *  Created on: Jul 28, 2016
 *      Author: rrsettgast
 */

#include "DomainPartition.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "NodeManager.hpp"

namespace geosx
{
using namespace dataRepository;

DomainPartition::DomainPartition( std::string const & name,
                                  ManagedGroup * const parent ) :
  ManagedGroup( name, parent )
{}

DomainPartition::~DomainPartition()
{}


void DomainPartition::BuildDataStructure( ManagedGroup * const )
{
  this->RegisterGroup<constitutive::ConstitutiveManager>(keys::ConstitutiveManager);

  this->RegisterGroup<NodeManager>(keys::FEM_Nodes);
  this->RegisterGroup<ElementManager>(keys::FEM_Elements);
//  this->RegisterGroup<FaceManager,ObjectManagerBase>(keys::FEM_Faces);
}

} /* namespace geosx */
