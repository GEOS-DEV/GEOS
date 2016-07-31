/*
 * DomainPartition.cpp
 *
 *  Created on: Jul 28, 2016
 *      Author: rrsettgast
 */

#include "DomainPartition.hpp"

namespace geosx
{

DomainPartition::DomainPartition(  std::string const & name,
                                   WrapperCollection * const parent ):
    WrapperCollection( name, parent )
{
  // TODO Auto-generated constructor stub

}

DomainPartition::~DomainPartition()
{
  // TODO Auto-generated destructor stub
}


void DomainPartition::Registration( dataRepository::WrapperCollection * const )
{
  RegisterChildWrapperCollection<DomainPartition>("domain");
  RegisterChildWrapperCollection<WrapperCollection>("solvers");


  RegisterWrapper< std::unordered_map<string,string> >("simulationParameterMap");

}

} /* namespace geosx */
