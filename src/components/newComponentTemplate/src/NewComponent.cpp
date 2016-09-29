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
                            ManagedGroup * const parent ):
    SolverBase(name,parent)
{

}

NewComponent::~NewComponent()
{}


void NewComponent::Initialize( dataRepository::ManagedGroup& /*domain*/ )
{}


void NewComponent::ReadXML( pugi::xml_node const & /*solverNode*/  )
{}


void NewComponent::Registration( dataRepository::ManagedGroup * const /*domain*/ )
{}


void NewComponent::TimeStep( real64 const& /*time_n*/,
               real64 const& /*dt*/,
               int32 const /*cycleNumber*/,
               dataRepository::ManagedGroup& /*domain*/ )
{}

REGISTER_CATALOG_ENTRY( SolverBase, NewComponent, std::string const &, dataRepository::ManagedGroup * const )

} /* namespace geosx */
