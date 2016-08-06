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
                                   WrapperCollection * const parent ):
    WrapperCollection( name, parent )
{
}

DomainPartition::~DomainPartition()
{
}


void DomainPartition::Registration( dataRepository::WrapperCollection * const )
{
  this->RegisterChildWrapperCollection<constitutive::ConstitutiveManager>(keys::ConstitutiveManager);
}

} /* namespace geosx */
