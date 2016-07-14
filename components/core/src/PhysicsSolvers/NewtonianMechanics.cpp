/*
 * NewtonianMechanics.cpp
 *
 *  Created on: Dec 4, 2014
 *      Author: rrsettgast
 */

#include "NewtonianMechanics.hpp"
#include <vector>

#include "../dataRepository/intrinsic/WrapperCollection.hpp"
#include "RAJA/RAJA.hxx"
#include "../dataRepository/DataTypes.hpp"

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

void NewtonianMechanics::RegisterDataObjects( dataRepository::WrapperCollection& domain )
{

  dataRepository::WrapperCollection& nodes = domain.RegisterChildDataObjectManager<dataRepository::WrapperCollection >("FEM_Nodes");
  dataRepository::WrapperCollection& elems = domain.RegisterChildDataObjectManager<dataRepository::WrapperCollection >("FEM_Elements");

  nodes.RegisterDataObject<real64_array>("TotalDisplacement");
  nodes.RegisterDataObject<real64_array>("IncrementalDisplacement");
  nodes.RegisterDataObject<real64_array>("Velocity");
  nodes.RegisterDataObject<real64_array>("Acceleration");

  dataRepository::Wrapper<real64_array>::rtype    X = nodes.RegisterDataObject<real64_array>("ReferencePosition")->getObjectData<real64_array>();
  dataRepository::Wrapper<real64_array>::rtype mass = nodes.RegisterDataObject<real64_array>("Mass")->getObjectData<real64_array>();;
  std::cout<<"breakpoint 2"<<std::endl;

  elems.RegisterDataObject<real64_array>("Strain");
  elems.RegisterDataObject("Force",  rtTypes::TypeIDs::real64_array_id );
  std::cout<<"breakpoint 3"<<std::endl;

  // HACK
  nodes.resize(101);
  elems.resize(100);
  std::cout<<"breakpoint 4"<<std::endl;

  for( uint64 a=0 ; a<nodes.size() ; ++a )
  {
    mass[a] = 1;
    X[a] = a * ( 1.0/nodes.size() );
  }
  std::cout<<"breakpoint 5"<<std::endl;


}

void NewtonianMechanics::TimeStep( real64 const& time_n,
                                   real64 const& dt,
                                   const int cycleNumber,
                                   dataRepository::WrapperCollection& domain )
{
  TimeStepExplicit( time_n, dt, cycleNumber, domain );
}

void NewtonianMechanics::TimeStepExplicit( real64 const& time_n,
                                           real64 const& dt,
                                           const int cycleNumber,
                                           dataRepository::WrapperCollection& domain )
{
  dataRepository::WrapperCollection& nodes = domain.GetDataObjectManager<dataRepository::WrapperCollection>("FEM_Nodes");
  dataRepository::WrapperCollection& elems = domain.GetDataObjectManager<dataRepository::WrapperCollection>("FEM_Elements");

  std::size_t const numNodes = nodes.size();
  std::size_t const numElems = elems.size();

  dataRepository::Wrapper<real64_array>::rtype    X = nodes.GetDataObjectData<real64_array>("ReferencePosition");
  dataRepository::Wrapper<real64_array>::rtype    u = nodes.GetDataObjectData<real64_array>("TotalDisplacement");
  dataRepository::Wrapper<real64_array>::rtype uhat = nodes.GetDataObjectData<real64_array>("IncrementalDisplacement");
  dataRepository::Wrapper<real64_array>::rtype vel  = nodes.GetDataObjectData<real64_array>("Velocity");
  dataRepository::Wrapper<real64_array>::rtype acc  = nodes.GetDataObjectData<real64_array>("Acceleration");
  dataRepository::Wrapper<real64_array>::rtype mass = nodes.GetDataObjectData<real64_array>("Mass");

  dataRepository::Wrapper<real64_array>::rtype    Felem = elems.GetDataObjectData<real64_array>("Force");
  dataRepository::Wrapper<real64_array>::rtype   Strain = elems.GetDataObjectData<real64_array>("Strain");

  Integration::OnePoint( acc, vel, dt/2, numNodes );
  Integration::OnePoint( vel, uhat, u, dt, numNodes );


  for( uint64 a=0 ; a<numElems ; ++a )
  {
    acc[a] = 0.0;
  }
  real64 Ey = 1e9;
  for( uint64 k=0 ; k<numElems ; ++k )
  {
    Strain[k] = ( u[k+1] - u[k] ) / ( X[k+1] - X[k] ) ;
    Felem[k] = Ey * Strain[k];
    acc[k]   -= Felem[k];
    acc[k+1] += Felem[k];
  }

  for( uint64 a=0 ; a<numNodes ; ++a )
  {
    acc[a] /= mass[a];
  }
  Integration::OnePoint( acc, vel, dt/2, numNodes );


  for( uint64 a=0 ; a<numNodes ; ++a )
  {
    std::cout<<vel[a]<<std::endl;
  }

  (void) time_n;
  (void) cycleNumber;
}


//REGISTER_FACTORY( NewtonianMechanics, SolverBase, std::string )
//REGISTER_SOLVER(NewtonianMechanics)

REGISTER_CATALOGUE_ENTRY( SolverBase, NewtonianMechanics )
} /* namespace ANST */
