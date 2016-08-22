/*
 * NewtonianMechanics.cpp
 *
 *  Created on: Dec 4, 2014
 *      Author: rrsettgast
 */

#include "SolidMechanicsLagrangianFEM.hpp"

#include <vector>
#include <math.h>

#include "dataRepository/SynchronizedGroup.hpp"
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
                                                            SynchronizedGroup * const parent ) :
  SolverBase( name, parent )
{}



SolidMechanics_LagrangianFEM::~SolidMechanics_LagrangianFEM()
{
  // TODO Auto-generated destructor stub
}

void SolidMechanics_LagrangianFEM::Registration( SynchronizedGroup * const domain )
{

  SynchronizedGroup& nodes = domain->RegisterGroup<SynchronizedGroup >(keys::FEM_Nodes);
  SynchronizedGroup& elems = domain->RegisterGroup<SynchronizedGroup >(keys::FEM_Elements);

  nodes.RegisterViewWrapper<real64_array>(keys::TotalDisplacement);
  nodes.RegisterViewWrapper<real64_array>(keys::IncrementalDisplacement);
  nodes.RegisterViewWrapper<real64_array>(keys::Velocity);
  nodes.RegisterViewWrapper<real64_array>(keys::Acceleration);

  elems.RegisterViewWrapper<real64_array>(keys::Strain);
  elems.RegisterViewWrapper(keys::Force,  rtTypes::TypeIDs::real64_array_id );

  // HACK
  nodes.resize(11);
  elems.resize(10);

  ViewWrapper<real64_array>::rtype    X = nodes.RegisterViewWrapper<real64_array>(keys::ReferencePosition).data();
  ViewWrapper<real64_array>::rtype mass = nodes.RegisterViewWrapper<real64_array>(keys::Mass).data();
  real64& Ey = *(elems.RegisterViewWrapper<real64>(keys::Ey).data());
  ViewWrapper<real64_array>::rtype K = elems.RegisterViewWrapper<real64_array>(keys::K).data();

  Ey = 10e9;

  double rho = 2700;
  double A = 1;
  double L = 1.0;


  std::cout<<"sound speed = "<<sqrt(Ey/rho)<<std::endl;
//  std::cout<<1.0/0.0<<std::endl;
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
                                             SynchronizedGroup& domain )
{
  TimeStepExplicit( time_n, dt, cycleNumber, domain );
}

void SolidMechanics_LagrangianFEM::TimeStepExplicit( real64 const& time_n,
                                                     real64 const& dt,
                                                     const int cycleNumber,
                                                     SynchronizedGroup& domain )
{
  SynchronizedGroup& nodes = domain.GetGroup<SynchronizedGroup>(keys::FEM_Nodes);
  SynchronizedGroup& elems = domain.GetGroup<SynchronizedGroup>(keys::FEM_Elements);

  std::size_t const numNodes = nodes.size();
  std::size_t const numElems = elems.size();

  ViewWrapper<real64_array>::rtype          X = nodes.getData<real64_array>(keys::ReferencePosition);
  ViewWrapper<real64_array>::rtype          u = nodes.getData<real64_array>(keys::TotalDisplacement);
  ViewWrapper<real64_array>::rtype       uhat = nodes.getData<real64_array>(keys::IncrementalDisplacement);
  ViewWrapper<real64_array>::rtype       vel  = nodes.getData<real64_array>(keys::Velocity);
  ViewWrapper<real64_array>::rtype       acc  = nodes.getData<real64_array>(keys::Acceleration);
  ViewWrapper<real64_array>::rtype_const mass = nodes.getWrapper<real64_array>(keys::Mass).data();

  ViewWrapper<real64_array>::rtype    Felem = elems.getData<real64_array>(keys::Force);
  ViewWrapper<real64_array>::rtype   Strain = elems.getData<real64_array>(keys::Strain);
  ViewWrapper<real64_array>::rtype_const  K = elems.getData<real64_array>(keys::K);

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


REGISTER_CATALOG_ENTRY( SolverBase, SolidMechanics_LagrangianFEM, std::string const &, SynchronizedGroup * const )
} /* namespace ANST */

