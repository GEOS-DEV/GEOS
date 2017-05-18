/*
 * DomainPartition.cpp
 *
 *  Created on: Jul 28, 2016
 *      Author: rrsettgast
 */

#include "DomainPartition.hpp"

#include "../MPI_Communications/SpatialPartition.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "NodeManager.hpp"
#include "managers/ElementRegionManager.hpp"

namespace geosx
{
using namespace dataRepository;

DomainPartition::DomainPartition( std::string const & name,
                                  ManagedGroup * const parent ) :
  ManagedGroup( name, parent )
{
  m_partition = new SpatialPartition();
}

DomainPartition::~DomainPartition()
{}


void DomainPartition::BuildDataStructure( ManagedGroup * const )
{
  this->RegisterGroup<constitutive::ConstitutiveManager>(keys::ConstitutiveManager);

  this->RegisterGroup<NodeManager>(keys::FEM_Nodes);
  this->RegisterGroup<ElementRegionManager>(keys::FEM_Elements);
  this->RegisterGroup<CellBlockManager>(keys::cellManager);
//  this->RegisterGroup<FaceManager,ObjectManagerBase>(keys::FEM_Faces);
}

} /* namespace geosx */
