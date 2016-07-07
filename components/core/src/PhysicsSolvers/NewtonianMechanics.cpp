/*
 * NewtonianMechanics.cpp
 *
 *  Created on: Dec 4, 2014
 *      Author: rrsettgast
 */

#include "NewtonianMechanics.hpp"
#include <vector>

#include "../dataRepository/intrinsic/DataObjectManager.hpp"
namespace geosx
{

NewtonianMechanics::NewtonianMechanics( const std::string& name ):
    SolverBase(name)
{

}




NewtonianMechanics::~NewtonianMechanics()
{
  // TODO Auto-generated destructor stub
}

void NewtonianMechanics::RegisterDataObjects( DataObjectManager& domain )
{

  DataObjectManager& nodes = domain.RegisterChildDataObjectManager<DataObjectManager >("FEM_Nodes");
  nodes.RegisterDataObject<real64_array>("ReferencePosition");
  nodes.RegisterDataObject<real64_array>("TotalDisplacement");
  nodes.RegisterDataObject<real64_array>("IncrementalDisplacement");
  nodes.RegisterDataObject<real64_array>("Velocity");
  nodes.RegisterDataObject<real64_array>("Acceleration");
}

void NewtonianMechanics::TimeStep( real64 const& time_n,
                                   real64 const& dt,
                                   const int cycleNumber,
                                   DataObjectManager& domain )
{
  TimeStepExplicit( time_n, dt, cycleNumber, domain );
}

void NewtonianMechanics::TimeStepExplicit( real64 const& time_n,
                                           real64 const& dt,
                                           const int cycleNumber,
                                           DataObjectManager& domain )
{
  std::cout<<"breakpoint 2: "<<LOCATION<<std::endl;
  DataObjectManager& nodes = domain.GetDataObjectManager<DataObjectManager>("FEM_Nodes");

  std::size_t numNodes = nodes.size();
  DataObject<real64_array>::rtype    X = nodes.GetDataObjectData<real64_array>("ReferencePosition");
  DataObject<real64_array>::rtype    u = nodes.GetDataObjectData<real64_array>("TotalDisplacement");
  DataObject<real64_array>::rtype uhat = nodes.GetDataObjectData<real64_array>("IncrementalDisplacement");
  DataObject<real64_array>::rtype vel  = nodes.GetDataObjectData<real64_array>("Velocity");
  DataObject<real64_array>::rtype acc  = nodes.GetDataObjectData<real64_array>("Acceleration");


  Integration::OnePoint( acc, vel, dt, numNodes );
  Integration::OnePoint( vel, uhat, u, dt, numNodes );


  (void) time_n;
  (void) cycleNumber;
}


//REGISTER_FACTORY( NewtonianMechanics, SolverBase, std::string )
//REGISTER_SOLVER(NewtonianMechanics)

REGISTER_CATALOGUE_ENTRY( SolverBase, NewtonianMechanics )
} /* namespace ANST */
