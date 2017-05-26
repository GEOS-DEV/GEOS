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
#include "constitutive/LinearElasticIsotropic.hpp"
#include "managers/NodeManager.hpp"
#include "managers/ElementRegionManager.hpp"
#include "finiteElement/FiniteElementManager.hpp"
#include "finiteElement/ElementLibrary/FiniteElement.h"
#include "finiteElement/Kinematics.h"
//#include "finiteElement/ElementLibrary/FiniteElementUtilities.h"

namespace geosx
{

namespace dataRepository
{
namespace keys
{
std::string const area = "area";
}
}

using namespace dataRepository;
using namespace constitutive;

void Integrate( const R2SymTensor& fieldvar,
                const R1Tensor* const dNdX,
                const realT& detJ,
                const realT& detF,
                const R2Tensor& Finv,
                const int numPoints,
                R1Tensor* const result)
{
  const realT integrationFactor = detJ * detF;

  R2Tensor P;
  P.AijBkj( fieldvar,Finv);
  P *= integrationFactor;

  for( int a=0 ; a<numPoints ; ++a )  // loop through all shape functions in element
  {
    result[a].minusAijBj( P , dNdX[a] );
  }


}


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

  docNode->setName(this->CatalogName());
  docNode->setSchemaType("Node");
  docNode->setShortDescription("An example solid mechanics solver");


  nodes.getDocumentationNode()->AllocateChildNode( keys::TotalDisplacement,
                                                   keys::TotalDisplacement,
                                                   -1,
                                                   "r1_array",
                                                   "r1_array",
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
                                                   "r1_array",
                                                   "r1_array",
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
                                                   "r1_array",
                                                   "r1_array",
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
                                                   "r1_array",
                                                   "r1_array",
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

}


void SolidMechanics_LagrangianFEM::BuildDataStructure( ManagedGroup * const domain )
{
  SolverBase::BuildDataStructure( domain );

  // Test auto-registration:
  RegisterDocumentationNodes();

}


void SolidMechanics_LagrangianFEM::InitializePreSubGroups( ManagedGroup& problemManager )
{
  ManagedGroup & domain = problemManager.GetGroup(keys::domain);
  ManagedGroup& nodes = domain.GetGroup<ManagedGroup >(keys::FEM_Nodes);
  ElementRegionManager& elementRegions = domain.GetGroup<ElementRegionManager >(keys::FEM_Elements);
  ConstitutiveManager & constitutiveManager = domain.GetGroup<ConstitutiveManager >(keys::ConstitutiveManager);
//  ConstitutiveManager::constitutiveMaps const & constitutiveMaps = constitutiveManager.GetMaps(0);

  ViewWrapper<real64_array>::rtype mass = nodes.getData<real64_array>(keys::Mass);
//  ViewWrapper<real64_array>::rtype K = elems.getData<real64_array>(keys::K);


  array< ConstitutiveWrapper< view_rtype<real64> > >  rho = constitutiveManager.GetParameterData<real64>(keys::density);

  elementRegions.forCellBlocks([&]( ManagedGroup& cellBlock ) -> void
  {
    auto const & detJ            = cellBlock.RegisterViewWrapper< Array2dT<real64> >(keys::detJ).data();
    auto const & constitutiveMap = cellBlock.getData< Array2dT<mapPair> >(keys::constitutiveMap);
    lArray2d const & elemsToNodes = cellBlock.getWrapper<lArray2d>(keys::nodeList).reference();// getData<lArray2d>(keys::nodeList);

    for( localIndex k=0 ; k<cellBlock.size() ; ++k )
    {
      localIndex const * nodeList = elemsToNodes[k];
      real64 const * detJq = detJ[k];
      for( localIndex q=0 ; q<constitutiveMap.Dimension(1) ; ++q )
      {
        mass[nodeList[q]] += rho[ constitutiveMap(k,q).first ].m_object[0] * detJq[q];
      }
    }
  });
  
  real64 totalMass = 0;
  for( localIndex a=0 ; a<mass.size() ; ++a )
  {
    std::cout<<"mass["<<a<<"] = "<<mass[a]<<std::endl;
    totalMass += mass[a];
  }
  std::cout<<"totalMass = "<<totalMass<<std::endl;
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
  ElementRegionManager & elemManager = domain.GetGroup<ElementRegionManager>( keys::FEM_Elements );
  FiniteElementManager const & numericalMethodManager = domain.getParent()->GetGroup<FiniteElementManager>(keys::finiteElementManager);
  ConstitutiveManager & constitutiveManager = domain.GetGroup<ConstitutiveManager >(keys::ConstitutiveManager);
  ConstitutiveManager::constitutiveMaps const & constitutiveMaps = constitutiveManager.GetMaps(0);

  localIndex const numNodes = nodes.size();

  view_rtype_const<r1_array> X = nodes.getData<r1_array>(keys::ReferencePosition);
  view_rtype<r1_array>       u = nodes.getData<r1_array>(keys::TotalDisplacement);
  view_rtype<r1_array>       uhat = nodes.getData<r1_array>(keys::IncrementalDisplacement);
  view_rtype<r1_array>       vel  = nodes.getData<r1_array>(keys::Velocity);
  view_rtype<r1_array>       acc  = nodes.getData<r1_array>(keys::Acceleration);
  view_rtype_const<real64_array> mass = nodes.getWrapper<real64_array>(keys::Mass).data();

  array< ConstitutiveWrapper< view_rtype<real64> > > const & E = constitutiveManager.GetParameterData<real64>(keys::youngsModulus);

  Integration::OnePoint( acc, vel, dt/2, numNodes );
  vel[0][0] = 1.0;
  vel[1][0] = 1.0;
  vel[2][0] = 1.0;
  vel[3][0] = 1.0;

  for( localIndex a=0 ; a<numNodes ; ++a )
  {
    vel[a][1] = 0.0;
    vel[a][2] = 0.0;
  }
  Integration::OnePoint( vel, uhat, u, dt, numNodes );

  for( localIndex a=0 ; a<numNodes ; ++a )
  {
    acc[a] = 0.0;
  }

  
  
  elemManager.forElementRegions([&]( ElementRegion & elementRegion )
  {
    auto const & numMethodName = elementRegion.getData<string>(keys::numericalMethod);
    FiniteElementSpace const & feSpace = numericalMethodManager.GetGroup<FiniteElementSpace>(numMethodName);

    elementRegion.forCellBlocks( [&] ( CellBlockSubRegion& cellBlock )
    {
      auto const & dNdX            = cellBlock.RegisterViewWrapper< Array1dT< Array2dT<R1Tensor> > >(keys::dNdX).data();
      auto const & detJ            = cellBlock.RegisterViewWrapper< Array2dT<real64> >(keys::detJ).data();

      view_rtype_const< Array2dT<mapPair> > constitutiveMap = cellBlock.getData< Array2dT<mapPair> >(keys::constitutiveMap);
      auto const & constitutiveGrouping = cellBlock.getReference< map< string, int32_array > >(dataRepository::keys::constitutiveGrouping);
      lArray2d const & elemsToNodes = cellBlock.getWrapper<lArray2d>(keys::nodeList).reference();// getData<lArray2d>(keys::nodeList);

      r1_array uhat_local( elemsToNodes.Dimension(1) );
      r1_array u_local( elemsToNodes.Dimension(1) );
      r1_array f_local( elemsToNodes.Dimension(1) );

      for( auto const & constitutiveGroup : constitutiveGrouping )
      {
        string const constitutiveName = constitutiveGroup.first;
        int32_array const & elementList = constitutiveGroup.second;


//        std::cout<<constitutiveManager.GetSubGroups().size()<<std::endl;
//        for( auto const & material : constitutiveManager.GetSubGroups() )
//        {
//          std::cout<<material.first<<std::endl;
//        }

        ConstitutiveBase & constitutiveModel = constitutiveManager.GetGroup<ConstitutiveBase>( constitutiveName );

        for( auto const k : elementList )
        {
          std::cout<<"Region "<<elementRegion.getName()<<", CellBlockSubRegion "<<cellBlock.getName()<<", Material "<<constitutiveName<<", index "<<k<<std::endl;

          f_local = 0;

          localIndex const * nodelist = elemsToNodes[k];

          CopyGlobalToLocal( nodelist,
                             u, uhat,
                             u_local, uhat_local );

          std::pair<int32,int32> const * constitutiveMapQuadrature = constitutiveMap[k];
          for( auto q=0 ; q<feSpace.m_finiteElement->n_quadrature_points() ; ++q )
          {
//            std::cout<<"quadrature point "<<q<<std::endl;
            R2Tensor dUhatdX, dUdX;
            CalculateGradient( dUhatdX ,uhat_local, dNdX[k][q] );
            CalculateGradient( dUdX ,u_local, dNdX[k][q] );
//            std::cout<<"dUhatdX"<<std::endl;
//            std::cout<<dUhatdX<<std::endl;

            R2Tensor F,L, Finv;
            {
            // calculate dv/dX
            R2Tensor dvdX = dUhatdX;
            dvdX *= 1.0 / dt;

            // calculate du/dX
            F = dUhatdX;
            F *= 0.5;
            F += dUdX;
            F.PlusIdentity(1.0);

            // calculate dX/du
            Finv.Inverse(F);

            // chain rule: calculate dv/du = dv/dX * dX/du
            L.AijBjk(dvdX, Finv);
            }

            // calculate gradient (end of step)
            F = dUhatdX;
            F += dUdX;
            F.PlusIdentity(1.0);

            realT detF = F.Det();

          // calculate element volume
  //        detJ_np1(k,q) = detJ(k,q) * detF;
  //        volume[k] += detJ_np1(k,q);
  //        initVolume += detJ(k,q);

            Finv.Inverse(F);

            // Calculate Rate of deformation tensor and rotationAxis tensor
            R2Tensor A, Rot;
            R2SymTensor Dadt;
            A.AijBjk(dUhatdX,Finv);
            IncrementalKinematics(A,Dadt,Rot);

//            std::cout<<"Dadt"<<std::endl;
//            std::cout<<Dadt<<std::endl;
//
//            std::cout<<"ConstitutiveIndex = "<<constitutiveMapQuadrature[q].second<<std::endl;
            R2SymTensor const totalStress = constitutiveModel.StateUpdatePoint(Dadt,Rot,constitutiveMapQuadrature[q].second,0);

//            std::cout<<"stress"<<q<<std::endl;
//            std::cout<<totalStress<<std::endl;
            Integrate( totalStress,
                                               dNdX[k][q],
                                               detJ(k,q),
                                               detF,
                                               Finv,
                                               f_local.size(),
                                               f_local.data() );

          }

          AddLocalToGlobal( nodelist, f_local, acc);
          
          std::cout<<"nodal forces"<<std::endl;
          for( auto force : f_local )
          {
            printf(" %4.2f ",force[0] );
//            std::cout<<force[0]<<std::endl;
          }
          printf("\n" );
          
        }
      }


    });


  });



  for( localIndex a=0 ; a<numNodes ; ++a )
  {
    acc[a] /= mass[a];
  }
  Integration::OnePoint( acc, vel, dt/2, numNodes );
  vel[0][0] = 1.0;
  vel[1][0] = 1.0;
  vel[2][0] = 1.0;
  vel[3][0] = 1.0;
  for( localIndex a=0 ; a<numNodes ; ++a )
  {
    vel[a][1] = 0.0;
    vel[a][2] = 0.0;
  }

  printf("\n            X at time %6.5f : ", time_n + dt );
  for( localIndex a=0 ; a<numNodes ; ++a )
  {
    printf(" %7.2f ",X[a][0] );
    //    std::cout<<vel[a]<<std::endl;
  }
  printf("\nAccelerations at time %6.5f : ", time_n + dt );
  for( localIndex a=0 ; a<numNodes ; ++a )
  {
    printf(" %7.2f ",acc[a][0] );
    //    std::cout<<vel[a]<<std::endl;
  }
  printf("\nVelocities    at time %6.5f : ", time_n + dt );
  for( localIndex a=0 ; a<numNodes ; ++a )
  {
    printf(" %7.2f ",vel[a][0] );
    //    std::cout<<vel[a]<<std::endl;
  }
  printf("\n" );

  (void) time_n;
  (void) cycleNumber;
}


REGISTER_CATALOG_ENTRY( SolverBase, SolidMechanics_LagrangianFEM, std::string const &, ManagedGroup * const )
} /* namespace ANST */
