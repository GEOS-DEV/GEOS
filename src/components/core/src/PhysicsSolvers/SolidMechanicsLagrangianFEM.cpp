/*
 * NewtonianMechanics.cpp
 *
 *  Created on: Dec 4, 2014
 *      Author: rrsettgast
 */

#include "SolidMechanicsLagrangianFEM.hpp"

#include <vector>
#include <math.h>

#include "../dataRepository/WrapperCollection.hpp"
#include "RAJA/RAJA.hxx"
#include "common/DataTypes.hpp"


namespace geosx
{

namespace dataRepository
{
namespace keys
{
std::string const K = "K";
std::string const Ey = "Ey";
}
}

using namespace dataRepository;

SolidMechanics_LagrangianFEM::SolidMechanics_LagrangianFEM( const std::string& name,
                                                            WrapperCollection * const parent ) :
  SolverBase( name, parent )
{}



SolidMechanics_LagrangianFEM::~SolidMechanics_LagrangianFEM()
{
  // TODO Auto-generated destructor stub
}

void SolidMechanics_LagrangianFEM::Registration( WrapperCollection * const domain )
{

  WrapperCollection& nodes = domain->RegisterChildWrapperCollection<WrapperCollection >(keys::FEM_Nodes);
  WrapperCollection& elems = domain->RegisterChildWrapperCollection<WrapperCollection >(keys::FEM_Elements);

  nodes.RegisterWrapper<real64_array>(keys::TotalDisplacement);
  nodes.RegisterWrapper<real64_array>(keys::IncrementalDisplacement);
  nodes.RegisterWrapper<real64_array>(keys::Velocity);
  nodes.RegisterWrapper<real64_array>(keys::Acceleration);

  elems.RegisterWrapper<real64_array>(keys::Strain);
  elems.RegisterWrapper(keys::Force,  rtTypes::TypeIDs::real64_array_id );

  // HACK
  nodes.resize(11);
  elems.resize(10);

  Wrapper<real64_array>::rtype    X = nodes.RegisterWrapper<real64_array>(keys::ReferencePosition).data();
  Wrapper<real64_array>::rtype mass = nodes.RegisterWrapper<real64_array>(keys::Mass).data();
  real64& Ey = *(elems.RegisterWrapper<real64>(keys::Ey).data());
  Wrapper<real64_array>::rtype K = elems.RegisterWrapper<real64_array>(keys::K).data();

  Ey = 10e9;

  double rho = 2700;
  double A = 1;
  double L = 1.0;


  std::cout<<"sound speed = "<<sqrt(Ey/rho)<<1.0/0.0<<std::endl;
  for( uint64 a=0 ; a<nodes.size() ; ++a )
  {
    X[a] = a * ( L/nodes.size() );
  }

  for( uint64 k=0 ; k<elems.size() ; ++k )
  {
    double dx = L / elems.size();
    mass[k] += rho * A * dx / 2;
    mass[k+1] += rho * A * dx / 2;
    K[k] = Ey*A*dx;
  }

}

void SolidMechanics_LagrangianFEM::TimeStep( real64 const& time_n,
                                             real64 const& dt,
                                             const int cycleNumber,
                                             WrapperCollection& domain )
{
  TimeStepExplicit( time_n, dt, cycleNumber, domain );
}

void SolidMechanics_LagrangianFEM::TimeStepExplicit( real64 const& time_n,
                                                     real64 const& dt,
                                                     const int cycleNumber,
                                                     WrapperCollection& domain )
{
  WrapperCollection& nodes = domain.GetChildWrapperCollection<WrapperCollection>(keys::FEM_Nodes);
  WrapperCollection& elems = domain.GetChildWrapperCollection<WrapperCollection>(keys::FEM_Elements);

  std::size_t const numNodes = nodes.size();
  std::size_t const numElems = elems.size();

  Wrapper<real64_array>::rtype          X = nodes.getData<real64_array>(keys::ReferencePosition);
  Wrapper<real64_array>::rtype          u = nodes.getData<real64_array>(keys::TotalDisplacement);
  Wrapper<real64_array>::rtype       uhat = nodes.getData<real64_array>(keys::IncrementalDisplacement);
  Wrapper<real64_array>::rtype       vel  = nodes.getData<real64_array>(keys::Velocity);
  Wrapper<real64_array>::rtype       acc  = nodes.getData<real64_array>(keys::Acceleration);
  Wrapper<real64_array>::rtype_const mass = nodes.getWrapper<real64_array>(keys::Mass).data();

  Wrapper<real64_array>::rtype    Felem = elems.getData<real64_array>(keys::Force);
  Wrapper<real64_array>::rtype   Strain = elems.getData<real64_array>(keys::Strain);
  Wrapper<real64_array>::rtype_const  K = elems.getData<real64_array>(keys::K);

  Integration::OnePoint( acc, vel, dt/2, numNodes );
  vel[0] = 1.0;
  Integration::OnePoint( vel, uhat, u, dt, numNodes );

  for( uint64 a=0 ; a<numElems ; ++a )
  {
    acc[a] = 0.0;
  }

  for( uint64 k=0 ; k<numElems ; ++k )
  {
    Strain[k] = ( u[k+1] - u[k] ) / ( X[k+1] - X[k] );
    Felem[k] = K[k] * Strain[k];
    acc[k]   += Felem[k];
    acc[k+1] -= Felem[k];
  }

  for( uint64 a=0 ; a<numNodes ; ++a )
  {
    acc[a] /= mass[a];
  }
  Integration::OnePoint( acc, vel, dt/2, numNodes );
  vel[0] = 1.0;


  for( uint64 a=0 ; a<numNodes ; ++a )
  {
    std::cout<<vel[a]<<std::endl;
  }

  (void) time_n;
  (void) cycleNumber;
}


REGISTER_CATALOG_ENTRY( SolverBase, SolidMechanics_LagrangianFEM, std::string const &, WrapperCollection * const )
} /* namespace ANST */

