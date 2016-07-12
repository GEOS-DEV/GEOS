/*
 * NewComponent.cpp
 *
 *  Created on: Jun 8, 2016
 *      Author: settgast
 */

#include "NewComponent.hpp"

namespace geosx
{

NewComponent::NewComponent( std::string const & name ):
    SolverBase(name)
{

}

NewComponent::~NewComponent()
{
  // TODO Auto-generated destructor stub
}

void RegisterDataObjects( dataRepository::WrapperCollection& /*domain*/ )
{}


void TimeStep( real64 const& /*time_n*/,
               real64 const& /*dt*/,
               int32 const /*cycleNumber*/,
               dataRepository::WrapperCollection& /*domain*/ )
{}

REGISTER_FACTORY( NewComponent, SolverBase, std::string )

} /* namespace geosx */
