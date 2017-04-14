/*
 * NewtonianMechanics.cpp
 *
 *  Created on: Dec 4, 2014
 *      Author: rrsettgast
 */

#include "SolidMechanicsLagrangianFEM.hpp"

#include <vector>
#include <math.h>

// #include "RAJA/RAJA.hxx"
#include "dataRepository/ManagedGroup.hpp"
#include "common/DataTypes.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "managers/NodeManager.hpp"

namespace geosx
{

namespace dataRepository
{
namespace keys
{
std::string const K = "K";
std::string const Ey = "Ey";
std::string const rho = "rho";
std::string const area = "area";
std::string const barLength = "barLength";
std::string const nElements = "nElements";
}
}

using namespace dataRepository;
using namespace constitutive;

SolidMechanics_LagrangianFEM::SolidMechanics_LagrangianFEM( const std::string& name,
                                                            ManagedGroup * const parent ) :
  SolverBase( name, parent )
{}



SolidMechanics_LagrangianFEM::~SolidMechanics_LagrangianFEM()
{
  // TODO Auto-generated destructor stub
}


void SolidMechanics_LagrangianFEM::FillDocumentationNode( dataRepository::ManagedGroup * const domain )
{
  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();
  SolverBase::FillDocumentationNode( domain );

  NodeManager& nodes    = domain->GetGroup<NodeManager>(keys::FEM_Nodes);
  ElementManager& elems = domain->GetGroup<ElementManager>(keys::FEM_Elements);


  docNode->setName(this->CatalogName());
  docNode->setSchemaType("Node");
  docNode->setShortDescription("An example solid mechanics solver");
  
  docNode->AllocateChildNode( "area",
                              keys::area,
                              -1,
                              "real64",
                              "real64",
                              "cross section area",
                              "cross section area",
                              "1.0",
                              "",
                              0,
                              1,
                              0 );

  nodes.getDocumentationNode()->AllocateChildNode( keys::TotalDisplacement,
                                                   keys::TotalDisplacement,
                                                   -1,
                                                   "real64_array",
                                                   "real64_array",
                                                   "Total Displacement",
                                                   "Total Displacement",
                                                   "0.0",
                                                   keys::nodeManager,
                                                   1,
                                                   0,
                                                   0 );

  nodes.getDocumentationNode()->AllocateChildNode( keys::IncrementalDisplacement,
                                                   keys::IncrementalDisplacement,
                                                   -1,
                                                   "real64_array",
                                                   "real64_array",
                                                   "Incremental Displacement",
                                                   "Incremental Displacement",
                                                   "0.0",
                                                   keys::nodeManager,
                                                   1,
                                                   0,
                                                   0 );

  nodes.getDocumentationNode()->AllocateChildNode( keys::Velocity,
                                                   keys::Velocity,
                                                   -1,
                                                   "real64_array",
                                                   "real64_array",
                                                   "Velocity",
                                                   "Velocity",
                                                   "0.0",
                                                   keys::nodeManager,
                                                   1,
                                                   0,
                                                   0 );

  nodes.getDocumentationNode()->AllocateChildNode( keys::Acceleration,
                                                   keys::Acceleration,
                                                   -1,
                                                   "real64_array",
                                                   "real64_array",
                                                   "Acceleration",
                                                   "Acceleration",
                                                   "0.0",
                                                   keys::nodeManager,
                                                   1,
                                                   0,
                                                   0 );

  nodes.getDocumentationNode()->AllocateChildNode( keys::Mass,
                                                   keys::Mass,
                                                   -1,
                                                   "real64_array",
                                                   "real64_array",
                                                   "Acceleration",
                                                   "Acceleration",
                                                   "0.0",
                                                   keys::nodeManager,
                                                   1,
                                                   0,
                                                   0 );

  elems.getDocumentationNode()->AllocateChildNode( keys::Strain,
                                                   keys::Strain,
                                                   -1,
                                                   "real64_array",
                                                   "real64_array",
                                                   "Acceleration",
                                                   "Acceleration",
                                                   "0.0",
                                                   keys::nodeManager,
                                                   1,
                                                   0,
                                                   0 );


  elems.getDocumentationNode()->AllocateChildNode( keys::Force,
                                                   keys::Force,
                                                   -1,
                                                   "real64_array",
                                                   "real64_array",
                                                   "Acceleration",
                                                   "Acceleration",
                                                   "0.0",
                                                   keys::nodeManager,
                                                   1,
                                                   0,
                                                   0 );



}


void SolidMechanics_LagrangianFEM::BuildDataStructure( ManagedGroup * const domain )
{
  SolverBase::BuildDataStructure( domain );

  // Test auto-registration:
  RegisterDocumentationNodes();

}


void SolidMechanics_LagrangianFEM::Initialize( dataRepository::ManagedGroup& domain )
{
  ManagedGroup& nodes = domain.GetGroup<ManagedGroup >(keys::FEM_Nodes);
  ElementManager& elems = domain.GetGroup<ElementManager >(keys::FEM_Elements);
  ConstitutiveManager & constitutiveManager = domain.GetGroup<ConstitutiveManager >(keys::ConstitutiveManager);

  ViewWrapper<r1_array>::rtype    X = nodes.getData<r1_array>(keys::ReferencePosition);
  ViewWrapper<real64_array>::rtype mass = nodes.getData<real64_array>(keys::Mass);
//  ViewWrapper<real64_array>::rtype K = elems.getData<real64_array>(keys::K);


  elems.forRegions([ this, &X, &mass ]( ElementRegion& elemRegion ) -> void
  {
    int32_array const & constitutiveMap = elemRegion.getData<int32_array>(keys::constitutiveMap);
    real64 rho = *(elemRegion.getData<real64>(keys::rho));
    lArray2d const & elemsToNodes = elemRegion.getData<lArray2d>(keys::nodeList);
    real64 area = 1;
    for( localIndex k=0 ; k<elemRegion.size() ; ++k )
    {
      localIndex const * nodeList = elemsToNodes[k];
      real64 dx = X[nodeList[1]][0] - X[nodeList[0]][0];
      mass[k]   += rho * area * dx / 2;
      mass[k+1] += rho * area * dx / 2;    }
  });
}

void SolidMechanics_LagrangianFEM::TimeStep( real64 const& time_n,
                                             real64 const& dt,
                                             const int cycleNumber,
                                             ManagedGroup& domain )
{
  TimeStepExplicit( time_n, dt, cycleNumber, domain );
}

void SolidMechanics_LagrangianFEM::TimeStepExplicit( real64 const& time_n,
                                                     real64 const& dt,
                                                     const int cycleNumber,
                                                     ManagedGroup& domain )
{
  ManagedGroup& nodes = domain.GetGroup<ManagedGroup>(keys::FEM_Nodes);
  ManagedGroup& elems = domain.GetGroup<ManagedGroup>(keys::FEM_Elements);

  localIndex const numNodes = nodes.size();
  localIndex const numElems = elems.size();

  ViewWrapper<real64_array>::rtype          X = nodes.getData<real64_array>(keys::ReferencePosition);
  ViewWrapper<real64_array>::rtype          u = nodes.getData<real64_array>(keys::TotalDisplacement);
  ViewWrapper<real64_array>::rtype       uhat = nodes.getData<real64_array>(keys::IncrementalDisplacement);
  ViewWrapper<real64_array>::rtype       vel  = nodes.getData<real64_array>(keys::Velocity);
  ViewWrapper<real64_array>::rtype       acc  = nodes.getData<real64_array>(keys::Acceleration);
  ViewWrapper<real64_array>::rtype_const mass = nodes.getWrapper<real64_array>(keys::Mass).data();

  ViewWrapper<real64_array>::rtype    Felem = elems.getData<real64_array>(keys::Force);
  ViewWrapper<real64_array>::rtype   Strain = elems.getData<real64_array>(keys::Strain);
  ViewWrapper<real64_array>::rtype_const  K = elems.getData<real64_array>(keys::K);


//  ViewWrapper<real64_array>::rtype          X2 = nodes.GetData(keys::ReferencePosition);


  Integration::OnePoint( acc, vel, dt/2, numNodes );
  vel[0] = 1.0;
  Integration::OnePoint( vel, uhat, u, dt, numNodes );

  for( localIndex a=0 ; a<numElems ; ++a )
  {
    acc[a] = 0.0;
  }

  for( localIndex k=0 ; k<numElems ; ++k )
  {
    Strain[k] = ( u[k+1] - u[k] ) / ( X[k+1] - X[k] );
    Felem[k] = K[k] * Strain[k];
    acc[k]   += Felem[k];
    acc[k+1] -= Felem[k];
  }

  for( localIndex a=0 ; a<numNodes ; ++a )
  {
    acc[a] /= mass[a];
  }
  Integration::OnePoint( acc, vel, dt/2, numNodes );
  vel[0] = 1.0;


  printf(" %6.5f : ", time_n + dt );
  for( localIndex a=0 ; a<numNodes ; ++a )
  {
    printf(" %4.2f ",vel[a] );
//    std::cout<<vel[a]<<std::endl;
  }
  printf("\n" );

  (void) time_n;
  (void) cycleNumber;
}


REGISTER_CATALOG_ENTRY( SolverBase, SolidMechanics_LagrangianFEM, std::string const &, ManagedGroup * const )
} /* namespace ANST */
