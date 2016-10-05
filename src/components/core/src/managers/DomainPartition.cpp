/*
 * DomainPartition.cpp
 *
 *  Created on: Jul 28, 2016
 *      Author: rrsettgast
 */

#include "DomainPartition.hpp"
#include "constitutive/ConstitutiveManager.hpp"

namespace geosx
{
using namespace dataRepository;

DomainPartition::DomainPartition(  std::string const & name,
                                   ManagedGroup * const parent ) :
  ManagedGroup( name, parent )
{}

DomainPartition::~DomainPartition()
{}


void DomainPartition::BuildDataStructure( ManagedGroup * const )
{
  this->RegisterGroup<constitutive::ConstitutiveManager>(keys::ConstitutiveManager);

}

} /* namespace geosx */
