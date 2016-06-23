/*
 * NewtonianMechanics.cpp
 *
 *  Created on: Dec 4, 2014
 *      Author: rrsettgast
 */

#include "NewtonianMechanics.hpp"
#include "dataRepository/DataObjectManager.hpp"
//#include "DomainPartition.hpp"
//#include <StringUtilities.hpp>

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
  /*
  DataObjectManager& nodes = domain.RegisterChildDataObjectManager<DataObjectManager >("FEM_Nodes");

  nodes.RegisterDataObject<rArray1d>("ReferencePosition");
  nodes.RegisterDataObject<rArray1d>("TotalDisplacement");
  nodes.RegisterDataObject<rArray1d>("IncrementalDisplacement");
  nodes.RegisterDataObject<rArray1d>("Velocity");
  nodes.RegisterDataObject<rArray1d>("Acceleration");
*/
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
/*
  DataObjectManager& nodes = domain.GetDataObjectManager<DataObjectManager>("FEM_Nodes");

  std::size_t numNodes = nodes.size();
  rArray1d const &    X = nodes.GetDataObjectData<rArray1d>("ReferencePosition");
  rArray1d &    u = nodes.GetDataObjectData<rArray1d>("TotalDisplacement");
  rArray1d & uhat = nodes.GetDataObjectData<rArray1d>("IncrementalDisplacement");
  rArray1d & vel  = nodes.GetDataObjectData<rArray1d>("Velocity");
  rArray1d const & acc  = nodes.GetDataObjectData<rArray1d>("Acceleration");

  Integration::OnePoint( acc.data(), vel.data(), dt, numNodes );
  Integration::OnePoint( vel.data(), uhat.data(), u.data(), dt, numNodes );


  (void) time_n;
  (void) cycleNumber;
  */
}


REGISTER_FACTORY( NewtonianMechanics, SolverBase, std::string )
} /* namespace ANST */
