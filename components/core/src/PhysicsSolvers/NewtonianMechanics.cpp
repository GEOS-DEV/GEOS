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
using namespace dataRepository;

NewtonianMechanics::NewtonianMechanics( const std::string& name ):
    SolverBase(name)
{

}




NewtonianMechanics::~NewtonianMechanics()
{
  // TODO Auto-generated destructor stub
}

void NewtonianMechanics::Registration( WrapperCollection& domain )
{

  WrapperCollection& nodes = domain.RegisterChildWrapperCollection<WrapperCollection >("FEM_Nodes");
  WrapperCollection& elems = domain.RegisterChildWrapperCollection<WrapperCollection >("FEM_Elements");

  nodes.RegisterWrapper<real64_array>("TotalDisplacement");
  nodes.RegisterWrapper<real64_array>("IncrementalDisplacement");
  nodes.RegisterWrapper<real64_array>("Velocity");
  nodes.RegisterWrapper<real64_array>("Acceleration");

  Wrapper<real64_array>::rtype    X = nodes.RegisterWrapper<real64_array>("ReferencePosition")->data<real64_array>();
  Wrapper<real64_array>::rtype mass = nodes.RegisterWrapper<real64_array>("Mass")->data<real64_array>();;

  elems.RegisterWrapper<real64_array>("Strain");
  elems.RegisterWrapper("Force",  rtTypes::TypeIDs::real64_array_id );

  // HACK
  nodes.resize(101);
  elems.resize(100);
  std::cout<<"breakpoint 4"<<std::endl;
  std::cout<<nodes.size()<<std::endl;
//  std::cout<<mass.size()<<std::cout;
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
                                   WrapperCollection& domain )
{
  TimeStepExplicit( time_n, dt, cycleNumber, domain );
}

void NewtonianMechanics::TimeStepExplicit( real64 const& time_n,
                                           real64 const& dt,
                                           const int cycleNumber,
                                           WrapperCollection& domain )
{
  WrapperCollection& nodes = domain.GetDataObjectManager<WrapperCollection>("FEM_Nodes");
  WrapperCollection& elems = domain.GetDataObjectManager<WrapperCollection>("FEM_Elements");

  std::size_t const numNodes = nodes.size();
  std::size_t const numElems = elems.size();

  Wrapper<real64_array>::rtype    X = nodes.getWrappedObjectData<real64_array>("ReferencePosition");
  Wrapper<real64_array>::rtype    u = nodes.getWrappedObjectData<real64_array>("TotalDisplacement");
  Wrapper<real64_array>::rtype uhat = nodes.getWrappedObjectData<real64_array>("IncrementalDisplacement");
  Wrapper<real64_array>::rtype vel  = nodes.getWrappedObjectData<real64_array>("Velocity");
//  Wrapper<real64_array>::rtype acc  = nodes.getWrappedObjectData<real64_array>("Acceleration");
//  Wrapper<real64_array>::rtype mass = nodes.getWrappedObjectData<real64_array>("Mass");


  Wrapper<real64_array>::rtype acc  = nodes.getWrappedObjectData<real64_array>("Acceleration");


  Wrapper<real64_array>::rtype mass = nodes.getWrapper<real64_array>("Mass").data();

  Wrapper<real64_array>::rtype    Felem = elems.getWrappedObjectData<real64_array>("Force");
  Wrapper<real64_array>::rtype   Strain = elems.getWrappedObjectData<real64_array>("Strain");

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
