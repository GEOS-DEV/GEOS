/*
 * NewComponent.cpp
 *
 *  Created on: Jun 8, 2016
 *      Author: settgast
 */

#include "NewComponent.hpp"

namespace geosx
{

NewComponent::NewComponent( std::string const & name,
                            WrapperCollection * const parent ):
    SolverBase(name,parent)
{

}

NewComponent::~NewComponent()
{}

void NewComponent::Registration( dataRepository::WrapperCollection& /*domain*/ )
{}


void NewComponent::TimeStep( real64 const& /*time_n*/,
               real64 const& /*dt*/,
               int32 const /*cycleNumber*/,
               dataRepository::WrapperCollection& /*domain*/ )
{}

REGISTER_CATALOG_ENTRY( SolverBase, NewComponent, std::string const &, dataRepository::WrapperCollection * const )

} /* namespace geosx */
