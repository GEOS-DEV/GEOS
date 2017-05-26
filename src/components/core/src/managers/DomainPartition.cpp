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

void DomainPartition::InitializationOrder( string_array & order )
{
  set<string> usedNames;
  {
    order.push_back(keys::ConstitutiveManager);
    usedNames.insert(keys::ConstitutiveManager);
  }

  {
    order.push_back(keys::FEM_Elements);
    usedNames.insert(keys::FEM_Elements);
  }


  for( auto const & subGroup : this->GetSubGroups() )
  {
    if( usedNames.count(subGroup.first) == 0 )
    {
      order.push_back(subGroup.first);
    }
  }
}


} /* namespace geosx */
