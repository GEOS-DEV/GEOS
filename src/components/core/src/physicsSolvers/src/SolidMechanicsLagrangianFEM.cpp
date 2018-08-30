/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

#include "SolidMechanicsLagrangianFEM.hpp"
#include "SolidMechanicsLagrangianFEMKernels_impl.hpp"
#include "../miniApps/SolidMechanicsLagrangianFEM-MiniApp/Layout.hpp"
#include "../miniApps/SolidMechanicsLagrangianFEM-MiniApp/ConstitutiveUpdate_impl.hpp"

#include <vector>
#include <math.h>


#include "common/TimingMacros.hpp"

#include "dataRepository/ManagedGroup.hpp"
#include <common/DataTypes.hpp>
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/LinearElasticIsotropic.hpp"
#include "managers/NumericalMethodsManager.hpp"
#include "finiteElement/FiniteElementSpaceManager.hpp"
#include "finiteElement/ElementLibrary/FiniteElement.h"
#include "finiteElement/Kinematics.h"
//#include "finiteElement/ElementLibrary/FiniteElementUtilities.h"
#include "managers/BoundaryConditions/BoundaryConditionManager.hpp"

#include "codingUtilities/Utilities.hpp"

#include "managers/DomainPartition.hpp"
#include "meshUtilities/ComputationalGeometry.hpp"
#include "MPI_Communications/CommunicationTools.hpp"
#include "../../rajaInterface/GEOS_RAJA_Interface.hpp"

//#define verbose 0 //Need to move this somewhere else

namespace geosx
{

namespace dataRepository
{
namespace keys
{}
}

using namespace dataRepository;
using namespace constitutive;
using namespace systemSolverInterface;

void Integrate( const R2SymTensor& fieldvar,
                const R1Tensor* const dNdX,
                real64 const& detJ,
                real64 const& detF,
                const R2Tensor& Finv,
                const localIndex numPoints,
                R1Tensor* const result)
{
  real64 const integrationFactor = detJ * detF;

  R2Tensor P;
  P.AijBkj( fieldvar,Finv);
  P *= integrationFactor;

  for( int a=0 ; a<numPoints ; ++a )  // loop through all shape functions in
                                      // element
  {
    result[a].minusAijBj( P, dNdX[a] );
  }

}


inline void HughesWinget(R2Tensor &Rot, R2SymTensor &Dadt, R2Tensor L, real64 dt){

  R2Tensor Omega;

  real64 * Dadt_data = Dadt.Data();
  real64 * L_data = L.Data();

  real64 *Omega_data = Omega.Data();

  //Omega = 0.5*(L - LT); 
  Omega_data[0] = 0.5*(L_data[0] - L_data[0]);
  Omega_data[1] = 0.5*(L_data[1] - L_data[3]);
  Omega_data[2] = 0.5*(L_data[2] - L_data[6]);

  Omega_data[3] = 0.5*(L_data[3] - L_data[1]);
  Omega_data[4] = 0.5*(L_data[4] - L_data[4]);
  Omega_data[5] = 0.5*(L_data[5] - L_data[7]);

  Omega_data[6] = 0.5*(L_data[6] - L_data[2]);
  Omega_data[7] = 0.5*(L_data[7] - L_data[5]);
  Omega_data[8] = 0.5*(L_data[8] - L_data[8]);


  //Dadt = 0.5*(L + LT)*dt;
  Dadt_data[0] = L_data[0]*dt;
  
  Dadt_data[1] = 0.5*(L_data[1] + L_data[3])*dt;
  Dadt_data[2] = L_data[4]*dt;
  
  Dadt_data[3] = 0.5*(L_data[6] + L_data[2])*dt;
  Dadt_data[4] = 0.5*(L_data[7] + L_data[5])*dt;
  Dadt_data[5] = L_data[8] *dt;
  
  
  R2Tensor IpR, ImR, ImRinv;
  IpR  = Omega;
  IpR *= 0.5;
  IpR.PlusIdentity(1.0);
  
  ImR  = Omega;
  ImR *= -0.5;
  ImR.PlusIdentity(1.0);                           
  
  ImRinv.Inverse(ImR);
  
  Rot.AijBjk(ImRinv,ImR);  
}

inline void LinearElasticIsotropic_Kernel(R2SymTensor &Dadt, R2SymTensor &TotalStress, R2Tensor &Rot,
                                          localIndex i, real64 bulkModulus, real64 shearModulus,
                                          array_view<real64,1> meanStress,
                                          array_view<R2SymTensor,1> devStress)
{
  real64 volumeStrain = Dadt.Trace();
  TotalStress = Dadt; 

  meanStress[i] += volumeStrain * bulkModulus;
  TotalStress *= 2.0 * shearModulus;
  
  devStress[i] += TotalStress;
  
  TotalStress.QijAjkQlk(devStress[i],Rot);
  devStress[i] = TotalStress;
  
  TotalStress.PlusIdentity(meanStress[i]);
  
}




SolidMechanics_LagrangianFEM::SolidMechanics_LagrangianFEM( const std::string& name,
                                                            ManagedGroup * const parent ):
  SolverBase( name, parent )
{
  getLinearSystemRepository()->
    SetBlockID( BlockIDs::displacementBlock, this->getName() );

}



SolidMechanics_LagrangianFEM::~SolidMechanics_LagrangianFEM()
{
  // TODO Auto-generated destructor stub
}


void SolidMechanics_LagrangianFEM::FillDocumentationNode()
{
  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();
  SolverBase::FillDocumentationNode();

  docNode->setName(this->CatalogName());
  docNode->setSchemaType("Node");
  docNode->setShortDescription("An example solid mechanics solver");


  docNode->AllocateChildNode( viewKeys.newmarkGamma.Key(),
                              viewKeys.newmarkGamma.Key(),
                              -1,
                              "real64",
                              "real64",
                              "newmark method Gamma",
                              "newmark method Gamma",
                              "0.5",
                              "",
                              0,
                              1,
                              1 );

  // correct default for this value is pow(newmarkGamma+0.5,2.0)/4.0
  docNode->AllocateChildNode( viewKeys.newmarkBeta.Key(),
                              viewKeys.newmarkBeta.Key(),
                              -1,
                              "real64",
                              "real64",
                              "newmark method Beta",
                              "newmark method Beta",
                              "0.25",
                              "",
                              0,
                              1,
                              1 );

  docNode->AllocateChildNode( viewKeys.massDamping.Key(),
                              viewKeys.massDamping.Key(),
                              -1,
                              "real64",
                              "real64",
                              "newmark method Beta",
                              "newmark method Beta",
                              "0",
                              "",
                              0,
                              1,
                              1 );

  docNode->AllocateChildNode( viewKeys.stiffnessDamping.Key(),
                              viewKeys.stiffnessDamping.Key(),
                              -1,
                              "real64",
                              "real64",
                              "newmark method Beta",
                              "newmark method Beta",
                              "0",
                              "",
                              0,
                              1,
                              1 );


  docNode->AllocateChildNode( viewKeys.timeIntegrationOption.Key(),
                              viewKeys.timeIntegrationOption.Key(),
                              -1,
                              "string",
                              "string",
                              "option for default time integration method",
                              "option for default time integration method",
                              "ExplicitDynamic",
                              "",
                              0,
                              1,
                              1 );

  docNode->AllocateChildNode( viewKeys.useVelocityEstimateForQS.Key(),
                              viewKeys.useVelocityEstimateForQS.Key(),
                              -1,
                              "integer",
                              "integer",
                              "option to use quasi-static deformation rate as a velocity in estimate of next solution",
                              "option to use quasi-static deformation rate as a velocity in estimate of next solution",
                              "ExplicitDynamic",
                              "",
                              0,
                              1,
                              1 );
}



void SolidMechanics_LagrangianFEM::FillOtherDocumentationNodes( dataRepository::ManagedGroup * const rootGroup )
{
  DomainPartition * domain  = rootGroup->GetGroup<DomainPartition>(keys::domain);

  for( auto & mesh : domain->getMeshBodies()->GetSubGroups() )
  {
    NodeManager * const nodes = mesh.second->group_cast<MeshBody*>()->getMeshLevel(0)->getNodeManager();
    cxx_utilities::DocumentationNode * const docNode = nodes->getDocumentationNode();

    docNode->AllocateChildNode( viewKeys.vTilde.Key(),
                                viewKeys.vTilde.Key(),
                                -1,
                                "r1_array",
                                "r1_array",
                                "intermediate velocity",
                                "intermediate velocity",
                                "0.0",
                                NodeManager::CatalogName(),
                                1,
                                0,
                                1 );

    docNode->AllocateChildNode( viewKeys.uhatTilde.Key(),
                                viewKeys.uhatTilde.Key(),
                                -1,
                                "r1_array",
                                "r1_array",
                                "intermediate incremental displacement",
                                "intermediate incremental displacement",
                                "0.0",
                                NodeManager::CatalogName(),
                                1,
                                0,
                                1 );

    docNode->AllocateChildNode( keys::TotalDisplacement,
                                keys::TotalDisplacement,
                                -1,
                                "r1_array",
                                "r1_array",
                                "Total Displacement",
                                "Total Displacement",
                                "0.0",
                                NodeManager::CatalogName(),
                                1,
                                0,
                                0 );

    docNode->AllocateChildNode( keys::IncrementalDisplacement,
                                keys::IncrementalDisplacement,
                                -1,
                                "r1_array",
                                "r1_array",
                                "Incremental Displacement",
                                "Incremental Displacement",
                                "0.0",
                                NodeManager::CatalogName(),
                                1,
                                0,
                                2 );

    docNode->AllocateChildNode( keys::Velocity,
                                keys::Velocity,
                                -1,
                                "r1_array",
                                "r1_array",
                                "Velocity",
                                "Velocity",
                                "0.0",
                                NodeManager::CatalogName(),
                                1,
                                0,
                                0 );

    docNode->AllocateChildNode( keys::Acceleration,
                                keys::Acceleration,
                                -1,
                                "r1_array",
                                "r1_array",
                                "Acceleration",
                                "Acceleration",
                                "0.0",
                                NodeManager::CatalogName(),
                                1,
                                0,
                                2 );

    docNode->AllocateChildNode( keys::Mass,
                                keys::Mass,
                                -1,
                                "real64_array",
                                "real64_array",
                                "Acceleration",
                                "Acceleration",
                                "0.0",
                                NodeManager::CatalogName(),
                                1,
                                0,
                                1 );

    docNode->AllocateChildNode( viewKeys.trilinosIndex.Key(),
                                viewKeys.trilinosIndex.Key(),
                                -1,
                                "globalIndex_array",
                                "globalIndex_array",
                                "Acceleration",
                                "Acceleration",
                                "-1",
                                NodeManager::CatalogName(),
                                1,
                                0,
                                1 );
  }

}

void SolidMechanics_LagrangianFEM::ReadXML_PostProcess()
{
  string tiOption = this->getReference<string>(viewKeys.timeIntegrationOption);

  if( tiOption == "ExplicitDynamic" )
  {
    this->m_timeIntegrationOption = timeIntegrationOption::ExplicitDynamic;
  }
  else if( tiOption == "ImplicitDynamic" )
  {
    this->m_timeIntegrationOption = timeIntegrationOption::ImplicitDynamic;
  }
  else if ( tiOption == "QuasiStatic" )
  {
    this->m_timeIntegrationOption = timeIntegrationOption::QuasiStatic;
  }
  else
  {
    GEOS_ERROR("invalid time integration option");
  }
}

void SolidMechanics_LagrangianFEM::FinalInitialization( ManagedGroup * const problemManager )
{
  DomainPartition * domain = problemManager->GetGroup<DomainPartition>(keys::domain);
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);

  NodeManager * const nodes = mesh->getNodeManager();
  FaceManager * const faceManager = mesh->getFaceManager();


  ElementRegionManager * elementRegionManager = mesh->getElemManager();
  ConstitutiveManager * constitutiveManager = domain->GetGroup<ConstitutiveManager >(keys::ConstitutiveManager);
//  ConstitutiveManager::constitutiveMaps const & constitutiveMaps =
// constitutiveManager->GetMaps(0);

  ViewWrapper<real64_array>::rtype mass = nodes->getData<real64_array>(keys::Mass);
//  ViewWrapper<real64_array>::rtype K = elems.getData<real64_array>(keys::K);

  ElementRegionManager::MaterialViewAccessor< array2d<real64> >
  rho = elementRegionManager->ConstructMaterialViewAccessor< array2d<real64> >("density",
                                                                               constitutiveManager);

  for( localIndex er=0 ; er<elementRegionManager->numRegions() ; ++er )
  {
    ElementRegion const * const elemRegion = elementRegionManager->GetRegion(er);
    for( localIndex esr=0 ; esr<elemRegion->numSubRegions() ; ++esr )
    {
      CellBlockSubRegion const * const cellBlock = elemRegion->GetSubRegion(esr);

      auto const & detJ            = cellBlock->getReference< array2d<real64> >(keys::detJ);
      auto const & constitutiveMap = cellBlock->getReference< std::pair< array2d<localIndex>,array2d<localIndex> > >(cellBlock->viewKeys().constitutiveMap);
      FixedOneToManyRelation const & elemsToNodes = cellBlock->getWrapper<FixedOneToManyRelation>(cellBlock->viewKeys().nodeList)->reference();// getData<array2d<localIndex>>(keys::nodeList);
//      array2d<real64> & rho = cellBlock->getReference< array2d<real64> >( string("density"));

      for( localIndex k=0 ; k < cellBlock->size() ; ++k )
      {
        arrayView1d<localIndex const> const nodeList = elemsToNodes[k];
        arrayView1d<real64 const> detJq = detJ[k];
        for( localIndex q=0 ; q<constitutiveMap.second.size(1) ; ++q )
        {
          mass[nodeList[q]] += rho[er][esr][0][k][q] * detJq[q];
        }
      }
    }
  }

  real64 totalMass = 0;
  for( localIndex a=0 ; a<nodes->size() ; ++a )
  {
    // std::cout<<"mass["<<a<<"] = "<<mass[a]<<std::endl;
    totalMass += mass[a];
  }
  std::cout<<"totalMass = "<<totalMass<<std::endl;


  real64_array & faceArea  = faceManager->getReference<real64_array>(FaceManager::
                                                                     viewKeyStruct::
                                                                     faceAreaString);
  array1d<array1d<localIndex>> const & facesToNodes = faceManager->nodeList();
  r1_array const & X = nodes->referencePosition();


  R1Tensor faceNormal, faceCenter;

  // loop over faces and calculate faceArea, faceNormal and faceCenter
  for (localIndex kf = 0; kf < faceManager->size(); ++kf)
  {
    faceArea[kf] = computationalGeometry::Centroid_3DPolygon(facesToNodes[kf],
                                                             X,
                                                             faceCenter,
                                                             faceNormal);
  }
}

real64 SolidMechanics_LagrangianFEM::SolverStep( real64 const& time_n,
                                             real64 const& dt,
                                             const int cycleNumber,
                                             DomainPartition * domain )
{
  real64 dtReturn = dt;
  if( m_timeIntegrationOption == timeIntegrationOption::ExplicitDynamic )
  {
    dtReturn = ExplicitStep( time_n, dt, cycleNumber, ManagedGroup::group_cast<DomainPartition*>(domain) );
  }
  else if( m_timeIntegrationOption == timeIntegrationOption::ImplicitDynamic ||
           m_timeIntegrationOption == timeIntegrationOption::QuasiStatic )
  {

    ImplicitStepSetup( time_n, dt, domain, getLinearSystemRepository() );

    dtReturn = NonlinearImplicitStep( time_n,
                                      dt,
                                      cycleNumber,
                                      domain->group_cast<DomainPartition*>(),
                                      getLinearSystemRepository() );

    ImplicitStepComplete( time_n, dt,  domain );

  }
return dtReturn;
}

real64 SolidMechanics_LagrangianFEM::ExplicitStep( real64 const& time_n,
                                                   real64 const& dt,
                                                   const int cycleNumber,
                                                   DomainPartition * const domain )
{

  GEOS_MARK_BEGIN(initialization);

  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  NodeManager * const nodes = mesh->getNodeManager();
  ElementRegionManager * elemManager = mesh->getElemManager();
  NumericalMethodsManager const * numericalMethodManager = domain->getParent()->GetGroup<NumericalMethodsManager>(keys::numericalMethodsManager);
  FiniteElementSpaceManager const * feSpaceManager = numericalMethodManager->GetGroup<FiniteElementSpaceManager>(keys::finiteElementSpaces);
  ConstitutiveManager * constitutiveManager = domain->GetGroup<ConstitutiveManager >(keys::ConstitutiveManager);

  BoundaryConditionManager * bcManager = BoundaryConditionManager::get();
  localIndex const numNodes = nodes->size();

  r1_array const &        X = nodes->getReference<r1_array>(nodes->viewKeys.referencePosition);
  real64_array const & mass = nodes->getReference<real64_array>(keys::Mass);
  r1_array &           vel  = nodes->getReference<r1_array>(keys::Velocity);

#if !defined(OBJECT_OF_ARRAYS_LAYOUT)

  r1_array &    u = nodes->getReference<r1_array>(keys::TotalDisplacement);
  r1_array & uhat = nodes->getReference<r1_array>(keys::IncrementalDisplacement);
  r1_array & acc  = nodes->getReference<r1_array>(keys::Acceleration);

#elif defined(OBJECT_OF_ARRAYS_LAYOUT)

  static geosxData acc_x = new double[numNodes];
  static geosxData acc_y = new double[numNodes];
  static geosxData acc_z = new double[numNodes];

  static geosxData uhat_x = new double[numNodes]; 
  static geosxData uhat_y = new double[numNodes]; 
  static geosxData uhat_z = new double[numNodes]; 

  static geosxData u_x = new double[numNodes]; 
  static geosxData u_y = new double[numNodes]; 
  static geosxData u_z = new double[numNodes];

  static bool setIc = true;
  if(setIc){    
    std::memset(acc_x, 0, numNodes*sizeof(double));
    std::memset(acc_y, 0, numNodes*sizeof(double));
    std::memset(acc_z, 0, numNodes*sizeof(double));

    std::memset(uhat_x, 0, numNodes*sizeof(double));
    std::memset(uhat_y, 0, numNodes*sizeof(double));
    std::memset(uhat_z, 0, numNodes*sizeof(double));

    std::memset(u_x, 0, numNodes*sizeof(double));
    std::memset(u_y, 0, numNodes*sizeof(double));
    std::memset(u_z, 0, numNodes*sizeof(double));
    setIc = false;
  }

#else
  GEOS_ERROR("Invalid data layout");
#endif

  GEOS_MARK_END(initialization);  

  GEOS_MARK_BEGIN(BC1);
#if !defined(OBJECT_OF_ARRAYS_LAYOUT)  
  bcManager->ApplyBoundaryConditionToField( time_n,
                                            domain,
                                            "nodeManager",
                                            keys::Acceleration );
#endif    
  GEOS_MARK_END(BC1);

  //3: v^{n+1/2} = v^{n} + a^{n} dt/2
  GEOS_CXX_MARK_LOOP_BEGIN(onepointloop,onepointloop1);
#if !defined(OBJECT_OF_ARRAYS_LAYOUT)  
  SolidMechanicsLagrangianFEMKernels::OnePoint( acc, vel, dt/2, numNodes );
#else
  SolidMechanicsLagrangianFEMKernels::OnePoint( acc_x, acc_y, acc_z,
                                                vel, dt/2, numNodes );
#endif  
  GEOS_CXX_MARK_LOOP_END(onepointloop);


  GEOS_MARK_BEGIN(BC2);
#if !defined(OBJECT_OF_ARRAYS_LAYOUT)  
  //  bcManager->ApplyBoundaryCondition( nodes, keys::Velocity, time_n + dt/2);

  bcManager->ApplyBoundaryConditionToField( time_n,
                                            domain,
                                            "nodeManager",
                                            keys::Velocity );

#endif  
  GEOS_MARK_END(BC2);


  //                     dydx, dy,   y, dx, length
  //4. x^{n+1} = x^{n} + v^{n+{1}/{2}} dt (x is displacement)
  GEOS_CXX_MARK_LOOP_BEGIN(onepointloop2,onepointloop2);
#if !defined(OBJECT_OF_ARRAYS_LAYOUT)  
  SolidMechanicsLagrangianFEMKernels::OnePoint( vel, uhat, u, dt, numNodes );
#else
  SolidMechanicsLagrangianFEMKernels::OnePoint(vel,uhat_x,uhat_y,uhat_z,
                                               u_x, u_y, u_z, dt, numNodes );
#endif  
  GEOS_CXX_MARK_LOOP_END(onepointloop2);


  GEOS_MARK_BEGIN(BC3);
#if !defined(OBJECT_OF_ARRAYS_LAYOUT)  
  //  bcManager->ApplyBoundaryCondition( this, &SolidMechanics_LagrangianFEM::ApplyDisplacementBC_explicit,
  //                                     nodes, keys::TotalDisplacement, time_n + dt, dt, u, uhat, vel );

  bcManager->ApplyBoundaryConditionToField( time_n+dt,
                                            domain,
                                            "nodeManager",
                                            keys::TotalDisplacement,
                                            [&]( BoundaryConditionBase const * const bc,
                                                set<localIndex> const & targetSet )->void
                                                {
    integer const component = bc->GetComponent();
    for( auto a : targetSet )
    {
      uhat[a][component] = u[a][component] - u[a][component];
      vel[a][component]  = uhat[a][component] / dt;
    }
                                                });


#endif  
  GEOS_MARK_END(BC3);

  //Set memory to zero
  GEOS_CXX_MARK_LOOP_BEGIN(memset,memset);

  FORALL_NODES( a, 0, numNodes )
  {
#if !defined(OBJECT_OF_ARRAYS_LAYOUT)
    acc[a] = 0;
#else
    acc_x[a] = 0;
    acc_y[a] = 0;
    acc_z[a] = 0; 
#endif    
  } END_FOR
  GEOS_CXX_MARK_LOOP_END(memset);

  ElementRegionManager::MaterialViewAccessor< array2d<real64> >
  meanStress = elemManager->ConstructMaterialViewAccessor< array2d<real64> >("MeanStress",
                                                                             constitutiveManager);

  ElementRegionManager::MaterialViewAccessor< array2d<R2SymTensor> >
  devStress = elemManager->ConstructMaterialViewAccessor< array2d<R2SymTensor> >("DeviatorStress",
                                                                                 constitutiveManager);

  ElementRegionManager::ConstitutiveRelationAccessor<ConstitutiveBase>
  constitutiveRelations = elemManager->ConstructConstitutiveAccessor<ConstitutiveBase>(constitutiveManager);


  //Step 5. Calculate deformation input to constitutive model and update state to
  // Q^{n+1}
  for( localIndex er=0 ; er<elemManager->numRegions() ; ++er )
  {
    ElementRegion * const elementRegion = elemManager->GetRegion(er);

    auto const & numMethodName = elementRegion->getData<string>(keys::numericalMethod);
    FiniteElementSpace const * feSpace = feSpaceManager->GetGroup<FiniteElementSpace>(numMethodName);

    for( localIndex esr=0 ; esr<elementRegion->numSubRegions() ; ++esr )
    {
      CellBlockSubRegion * const cellBlock = elementRegion->GetSubRegion(esr);
      //always needed
      array3d< R1Tensor > const & dNdX = cellBlock->getReference< array3d< R1Tensor > >(keys::dNdX);

      array2d<real64> const & detJ            = cellBlock->getReference< array2d<real64> >(keys::detJ);

      array_view<localIndex,2> const elemsToNodes = cellBlock->getWrapper<FixedOneToManyRelation>(cellBlock->viewKeys().nodeList)->reference().View();// getData<array2d<localIndex>>(keys::nodeList);

      localIndex const numNodesPerElement = elemsToNodes.size(1);

      localIndex const numQuadraturePoints = feSpace->m_finiteElement->n_quadrature_points();

      //Storage for holding intermediate results
#if defined(EXTERNAL_KERNELS) && defined(THREE_KERNEL_UPDATE) 
      static geosxData Dadt = new double[localMatSz*inumQuadraturePoints*elementList.size()];
      static geosxData Rot  = new double[localMatSz*inumQuadraturePoints*elementList.size()];
      static geosxData detF = new double[inumQuadraturePoints*elementList.size()];
      static geosxData inverseF = new double[localMatSz*inumQuadraturePoints*elementList.size()];
#endif

      //
      //Internal GEOSX Kernel
      //
      GEOS_CXX_MARK_LOOP_BEGIN(elemLoop,elemLoop);
#if !defined(EXTERNAL_KERNELS)         


      //          geosx::forall_in_set<elemPolicy>(elementList.data(), elementList.size(), GEOSX_LAMBDA ( globalIndex k) {
      for( localIndex k=0 ; k<cellBlock->size() ; ++k )
      {
        r1_array uhat_local( numNodesPerElement );
        r1_array u_local( numNodesPerElement );
        r1_array f_local( numNodesPerElement );

        f_local = R1Tensor(0.0);
        arrayView1d<localIndex const> const nodelist = elemsToNodes[k];

        CopyGlobalToLocal( nodelist,
                           u, uhat,
                           u_local.data(), uhat_local.data(), numNodesPerElement );


        //Compute Quadrature
        for(auto q = 0 ; q<numQuadraturePoints ; ++q)
        {

          R2Tensor dUhatdX, dUdX;
          CalculateGradient( dUhatdX,uhat_local, dNdX[k][q] );
          CalculateGradient( dUdX,u_local, dNdX[k][q] );

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
          real64 detF = F.Det();


          // calculate element volume
          //        detJ_np1(k,q) = detJ(k,q) * detF;
          //        volume[k] += detJ_np1(k,q);
          //        initVolume += detJ(k,q);


          Finv.Inverse(F);


          //-------------------------[Incremental Kinematics]----------------------------------
          R2Tensor Rot;
          R2SymTensor Dadt;
          HughesWinget(Rot, Dadt, L, dt);
          //-----------------------[Compute Total Stress - Linear Elastic Isotropic]-----------

          constitutiveRelations[er][esr][0]->StateUpdatePoint( Dadt, Rot, k, q, 0);

          R2SymTensor TotalStress;
          TotalStress = devStress[er][esr][0][k][q];
          TotalStress.PlusIdentity( meanStress[er][esr][0][k][q] );

          //----------------------

          Integrate( TotalStress, dNdX[k][q], detJ(k,q), detF, Finv, f_local.size(), f_local.data() );

        }//quadrature loop


        AddLocalToGlobal(nodelist, f_local.data(), acc, numNodesPerElement);

      } //Element loop
#else// defined(EXTERNAL_KERNELS) 

      //
      // Setup for external kernels
      //

      //Setup pointers
      const real64 *Xptr = static_cast<const real64 *>(X[0].Data());
      geosxData imeanStress   = static_cast<real64*>(&(meanStress[er][esr][0][0][0])); //Symmetric tensors
      geosxData idevStress    = devStress[er][esr][0][0][0].Data();
      localIndex const *  iconstitutiveMap = reinterpret_cast<localIndex const*>(constitutiveMap.second.data());
      const real64 * idetJ  = detJ.data();

      //Setup pointer for the external constitutive update
      void (*externConstitutiveUpdate)(real64 D[local_dim][local_dim], real64 Rot[local_dim][local_dim],
          localIndex m, localIndex q, globalIndex k, geosxData devStressData2,
          geosxData meanStress2, real64 shearModulus2, real64 bulkModulus2, localIndex NoElem);
      externConstitutiveUpdate = UpdateStatePoint;

#if defined(COMPUTE_SHAPE_FUN)
      static bool computeP = true;
      static P_Wrapper P;
      if(computeP)
      {
        generateP(P, 8, 8);
        computeP = false;
      }
#endif


#if defined(ARRAY_OF_OBJECTS_LAYOUT)

      //Setup pointers
      geosxData  ivel = static_cast<real64*>(vel[0].Data());
      geosxData iuhat = static_cast<real64*>(uhat[0].Data());
      geosxData iu = static_cast<real64*>(u[0].Data());
      geosxData idNdX = const_cast<real64*>(dNdX[0][0][0].Data());

      //Calculation is done in a monolithic kernel
#if !defined(THREE_KERNEL_UPDATE) && !defined(COMPUTE_SHAPE_FUN)


      SolidMechanicsLagrangianFEMKernels::ArrayOfObjectsKernel<elemPolicy>(elementList.size(),elementList.data(), dt,
                                                                           elemsToNodes.data(), iu, iuhat, idNdX,
                                                                           iconstitutiveMap, idevStress, imeanStress,
                                                                           shearModulus, bulkModulus, detJ.data(), iacc, externConstitutiveUpdate);
#elif !defined(THREE_KERNEL_UPDATE) && defined(COMPUTE_SHAPE_FUN)

      SolidMechanicsLagrangianFEMKernels::ArrayOfObjectsKernel_Shape<elemPolicy>(elementList.size(),elementList, dt,
                                                                                 elemsToNodes.data(), iu, iuhat, Xptr, P,
                                                                                 iconstitutiveMap, idevStress, imeanStress,
                                                                                 shearModulus, bulkModulus, detJ.data(), iacc, externConstitutiveUpdate);
      //Calculation is split across three kernels
#elif defined(THREE_KERNEL_UPDATE)

      //Kinematic step
      SolidMechanicsLagrangianFEMKernels::ArrayOfObjects_KinematicKernel<elemPolicy>(elementList.size(),elementList, dt, elemsToNodes.data(), iu, iuhat, idNdX,
                                                                                     iconstitutiveMap, idevStress, imeanStress,
                                                                                     shearModulus, bulkModulus, detJ.data(), iacc, Dadt, Rot, detF, inverseF);
      //Constitutive step
      SolidMechanicsLagrangianFEMKernels::ConstitutiveUpdateKernel<elemPolicy>(elementList.size(), elementList, Dadt, Rot, iconstitutiveMap, idevStress,
                                                                               imeanStress, shearModulus, bulkModulus);

      //Integration step
      SolidMechanicsLagrangianFEMKernels::ArrayOfObjects_IntegrationKernel<elemPolicy>(elementList.size(),elementList, dt, elemsToNodes.data(), iu, iuhat, idNdX,
                                                                                       iconstitutiveMap, idevStress, imeanStress,
                                                                                       shearModulus, bulkModulus, detJ.data(), iacc, Dadt, Rot, detF, inverseF);
#else

      //Throw error if invalid kernel is asked for
      GEOS_ERROR("Invalid External Kernel");

#endif // THREE_KERNEL UPDATE


#elif defined(OBJECT_OF_ARRAYS_LAYOUT)

      //Split the shape function derivatives into three arrays
      static geosxData dNdX_x = new double[inumNodesPerElement*inumQuadraturePoints*elementList.size()];
      static geosxData dNdX_y = new double[inumNodesPerElement*inumQuadraturePoints*elementList.size()];
      static geosxData dNdX_z = new double[inumNodesPerElement*inumQuadraturePoints*elementList.size()];

      static bool copy = true;
      if(copy){

        //Generate shape function derivatives
        //const real64 *Xptr = static_cast<const real64 *>(X[0].Data());
        make_dNdX(dNdX_x, dNdX_y, dNdX_z,
                  Xptr, elemsToNodes.data(),elementList.size(), inumNodesPerElement, inumQuadraturePoints);

        //Verify correctness
        for(localIndex k = 0; k < elementList.size(); ++k){
          for(localIndex a = 0; a < inumNodesPerElement; ++a){
            for(localIndex q = 0; q < inumQuadraturePoints; ++q){
              localIndex id = q + inumQuadraturePoints*(a + inumNodesPerElement*k);
              assert( std::abs( dNdX_x[id] - dNdX[k][a][q][0]) < 1e-12);
              assert( std::abs( dNdX_y[id] - dNdX[k][a][q][1]) < 1e-12);
              assert( std::abs( dNdX_z[id] - dNdX[k][a][q][2]) < 1e-12);
            }
          }
        }
        std::cout<<"Successful copy !"<<std::endl;
        copy = false;
      }

#if !defined(THREE_KERNEL_UPDATE) && !defined(COMPUTE_SHAPE_FUN)

      //Carry out computation in a monolithic kernel
      SolidMechanicsLagrangianFEMKernels::ObjectOfArraysKernel<elemPolicy>(elementList.size(), elementList, dt,elemsToNodes.data(),
                                                                           u_x, u_y, u_z, uhat_x, uhat_y, uhat_z, dNdX_x, dNdX_y, dNdX_z,
                                                                           iconstitutiveMap, idevStress, imeanStress, shearModulus,
                                                                           bulkModulus, detJ.data(), acc_x, acc_y, acc_z, externConstitutiveUpdate);
#elif !defined(THREE_KERNEL_UPDATE) && defined(COMPUTE_SHAPE_FUN)

      SolidMechanicsLagrangianFEMKernels::ObjectOfArraysKernel_Shape<elemPolicy>(elementList.size(), elementList, dt,elemsToNodes.data(),
                                                                                 u_x, u_y, u_z, uhat_x, uhat_y, uhat_z, Xptr, P,
                                                                                 iconstitutiveMap, idevStress, imeanStress, shearModulus,
                                                                                 bulkModulus, detJ.data(), acc_x, acc_y, acc_z, externConstitutiveUpdate);

#elif defined(THREE_KERNEL_UPDATE)
      //Kinematic step
      SolidMechanicsLagrangianFEMKernels::ObjectOfArrays_KinematicKernel<elemPolicy>(elementList.size(),elementList, dt, elemsToNodes.data(), u_x,u_y,u_z,
                                                                                     uhat_x, uhat_y, uhat_z, dNdX_x, dNdX_y, dNdX_z,
                                                                                     iconstitutiveMap, idevStress, imeanStress,
                                                                                     shearModulus, bulkModulus, detJ.data(),
                                                                                     acc_x, acc_y, acc_z,
                                                                                     Dadt, Rot, detF, inverseF);
      //Constitutive step
      SolidMechanicsLagrangianFEMKernels::ConstitutiveUpdateKernel<elemPolicy>(elementList.size(), elementList, Dadt, Rot, iconstitutiveMap, idevStress,
                                                                               imeanStress, shearModulus, bulkModulus);

      //Integration step
      SolidMechanicsLagrangianFEMKernels::ObjectOfArrays_IntegrationKernel<elemPolicy>(elementList.size(),elementList, dt, elemsToNodes.data(), u_x,u_y,u_z,
                                                                                       uhat_x, uhat_y, uhat_z, dNdX_x, dNdX_y, dNdX_z,
                                                                                       iconstitutiveMap, idevStress, imeanStress,
                                                                                       shearModulus, bulkModulus, detJ.data(),
                                                                                       acc_x, acc_y, acc_z,
                                                                                       Dadt, Rot, detF, inverseF);
#endif //defined THREE_KERNEL_UPDATE


#else

      GEOS_ERROR("Invalid External Kernel");

#endif //if defined(OBJECT_OF_ARRAYS_LAYOUT)


#endif// If !defined(EXTERNAL_KERNELS)
      GEOS_CXX_MARK_LOOP_END(elemLoop);

    } //Element Region

  } //Element Manager


//Compute Force : Point-wise computations
GEOS_CXX_MARK_LOOP_BEGIN(computeForce,computeForce);
FORALL_NODES( a, 0, numNodes )
{
#if !defined(OBJECT_OF_ARRAYS_LAYOUT)    
  acc[a] /=mass[a];
#else    
  acc_x[a] /=mass[a];
  acc_y[a] /=mass[a];
  acc_z[a] /=mass[a];
#endif
} END_FOR
GEOS_CXX_MARK_LOOP_END(computeForce);


//Integration::OnePoint( acc, vel, dt/2, numNodes );
GEOS_CXX_MARK_LOOP_BEGIN(onepointloop3,onepointloop3);
#if !defined(OBJECT_OF_ARRAYS_LAYOUT)      
SolidMechanicsLagrangianFEMKernels::OnePoint(acc, vel, (dt/2), numNodes);
#else  
SolidMechanicsLagrangianFEMKernels::OnePoint(acc_x, acc_y, acc_z, vel, (dt/2), numNodes);
#endif  
GEOS_CXX_MARK_LOOP_END(onepointloop3);


GEOS_MARK_BEGIN(BC4);
#if !defined(OBJECT_OF_ARRAYS_LAYOUT)
//bcManager->ApplyBoundaryCondition( nodes, keys::Velocity, time_n + dt);
bcManager->ApplyBoundaryConditionToField( time_n, domain, "nodeManager", keys::Velocity );

#endif
GEOS_MARK_END(BC4);


(void) cycleNumber;

return dt;
}


void SolidMechanics_LagrangianFEM::ApplyDisplacementBC_implicit( real64 const time,
                                                                 DomainPartition & domain,
                                                                 EpetraBlockSystem & blockSystem )
{

  BoundaryConditionManager const * const bcManager = BoundaryConditionManager::get();

  bcManager->ApplyBoundaryCondition( time,
                                     &domain,
                                     "nodeManager",
                                     keys::TotalDisplacement,
                                     [&]( BoundaryConditionBase const * const bc,
                                          string const &,
                                          set<localIndex> const & targetSet,
                                          ManagedGroup * const targetGroup,
                                          string const fieldName )->void
    {
    bc->ApplyBoundaryConditionToSystem<BcEqual>( targetSet,
                                                 time,
                                                 targetGroup,
                                                 fieldName,
                                                 viewKeys.trilinosIndex.Key(),
                                                 3,
                                                 &blockSystem,
                                                 BlockIDs::displacementBlock );
  });
}


void SolidMechanics_LagrangianFEM::ApplyTractionBC( DomainPartition * const domain,
                                                    real64 const time,
                                                    systemSolverInterface::EpetraBlockSystem & blockSystem )
{
  BoundaryConditionManager * const bcManager = BoundaryConditionManager::get();
  NewFunctionManager * const functionManager = NewFunctionManager::Instance();

  FaceManager * const faceManager = domain->getMeshBody(0)->getMeshLevel(0)->getFaceManager();
  NodeManager * const nodeManager = domain->getMeshBody(0)->getMeshLevel(0)->getNodeManager();

  real64_array const & faceArea  = faceManager->getReference<real64_array>("faceArea");
  array1d<localIndex_array> const & facesToNodes = faceManager->nodeList();

  globalIndex_array const &
  blockLocalDofNumber = nodeManager->getReference<globalIndex_array>(viewKeys.trilinosIndex);

  Epetra_FEVector * const rhs = blockSystem.GetResidualVector( BlockIDs::displacementBlock );



  bcManager->ApplyBoundaryCondition( time,
                                     domain,
                                     "faceManager",
                                     string("Traction"),
                                     [&]( BoundaryConditionBase const * const bc,
                                         string const &,
                                         set<localIndex> const & targetSet,
                                         ManagedGroup * const targetGroup,
                                         string const fieldName ) -> void
  {
    string const functionName = bc->getData<string>( BoundaryConditionBase::viewKeyStruct::functionNameString);

    globalIndex_array nodeDOF;
    real64_array nodeRHS;
    integer const component = bc->GetComponent();

    if( functionName.empty() )
    {
      integer counter=0;
      for( auto kf : targetSet )
      {
        localIndex const numNodes = facesToNodes[kf].size();
        nodeDOF.resize( numNodes );
        nodeRHS.resize( numNodes );
        for( localIndex a=0 ; a<numNodes ; ++a )
        {
          nodeDOF[a] = 3*blockLocalDofNumber[facesToNodes[kf][a]]+component;
          nodeRHS[a] = bc->GetScale() * faceArea[kf] / numNodes;
        }
        rhs->SumIntoGlobalValues( integer_conversion<int>(numNodes), nodeDOF.data(), nodeRHS.data() );
      }
    }
    else
    {
      FunctionBase const * const function  = functionManager->GetGroup<FunctionBase>(functionName);
      GEOS_ASSERT( function!=nullptr, "SolidMechanicsLagrangianFEM::ApplyTractionBC() application function not found");

        if( function->isFunctionOfTime()==2 )
        {
          real64 value = bc->GetScale() * function->Evaluate( &time );
          for( auto kf : targetSet )
          {
            localIndex const numNodes = facesToNodes[kf].size();
            nodeDOF.resize( numNodes );
            nodeRHS.resize( numNodes );
            for( localIndex a=0 ; a<numNodes ; ++a )
            {
              nodeDOF[a] = blockLocalDofNumber[facesToNodes[kf][a]]+component;
              nodeRHS[a] = value;
            }
            rhs->SumIntoGlobalValues( integer_conversion<int>(nodeDOF.size()), nodeDOF.data(), nodeRHS.data() );
          }
        }
        else
        {
          real64_array result;
          result.resize( targetSet.size() );
          function->Evaluate( faceManager, time, targetSet, result );

          integer counter=0;
          for( auto kf : targetSet )
          {
            localIndex const numNodes = facesToNodes[kf].size();
            nodeDOF.resize( numNodes );
            nodeRHS.resize( numNodes );
            for( localIndex a=0 ; a<numNodes ; ++a )
            {
              nodeDOF[a] = blockLocalDofNumber[facesToNodes[kf][a]]+component;
              nodeRHS[a] = result[a];
            }
            rhs->SumIntoGlobalValues( integer_conversion<int>(nodeDOF.size()), nodeDOF.data(), nodeRHS.data() );
          }
      }
    }
  });

}

void
SolidMechanics_LagrangianFEM::
ImplicitStepSetup( real64 const& time_n,
                   real64 const& dt,
                   DomainPartition * const domain,
                   systemSolverInterface::EpetraBlockSystem * const blockSystem )
{

  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ManagedGroup * const nodeManager = mesh->getNodeManager();

  if( this->m_timeIntegrationOption == timeIntegrationOption::ImplicitDynamic )
  {
    view_rtype_const<r1_array> v_n = nodeManager->getData<r1_array>(keys::Velocity);
    view_rtype_const<r1_array> a_n = nodeManager->getData<r1_array>(keys::Acceleration);
    view_rtype<r1_array> vtilde   = nodeManager->getData<r1_array>(viewKeys.vTilde);
    view_rtype<r1_array> uhatTilde   = nodeManager->getData<r1_array>(viewKeys.uhatTilde);

    view_rtype<r1_array> uhat  = nodeManager->getData<r1_array>(keys::IncrementalDisplacement);
    view_rtype<r1_array> disp = nodeManager->getData<r1_array>(keys::TotalDisplacement);

    localIndex const numNodes = nodeManager->size();
    real64 const newmarkGamma = this->getReference<real64>(viewKeys.newmarkGamma);
    real64 const newmarkBeta = this->getReference<real64>(viewKeys.newmarkBeta);

    for( auto a = 0 ; a < numNodes ; ++a )
    {
      for( int i=0 ; i<3 ; ++i )
      {
        vtilde[a][i] = v_n[a][i] + (1.0-newmarkGamma) * a_n[a][i] * dt;
        uhatTilde[a][i] = ( v_n[a][i] + 0.5 * ( 1.0 - 2.0*newmarkBeta ) * a_n[a][i] * dt ) *dt;
        uhat[a][i] = uhatTilde[a][i];
        disp[a][i] += uhatTilde[a][i];
      }
    }
  }
  else if( this->m_timeIntegrationOption == timeIntegrationOption::QuasiStatic  )
  {

    view_rtype<r1_array> uhat  = nodeManager->getData<r1_array>(keys::IncrementalDisplacement);
    integer const useVelocityEstimateForQS = this->getReference<integer>(viewKeys.useVelocityEstimateForQS);
    localIndex const numNodes = nodeManager->size();

    if( useVelocityEstimateForQS==1 )
    {
      view_rtype_const<r1_array> v_n = nodeManager->getData<r1_array>(keys::Velocity);
      view_rtype<r1_array> disp = nodeManager->getData<r1_array>(keys::TotalDisplacement);

      for( auto a = 0 ; a < numNodes ; ++a )
      {
        for( int i=0 ; i<3 ; ++i )
        {
          uhat[a][i] = v_n[a][i] * dt;
          disp[a][i] += uhat[a][i];
        }
      }
    }
    else
    {
      for( auto a = 0 ; a < numNodes ; ++a )
      {
        uhat[a] = 0.0;
      }
    }
  }

  SetupSystem( domain, blockSystem );
}

void SolidMechanics_LagrangianFEM::ImplicitStepComplete( real64 const & time_n,
                                                             real64 const & dt,
                                                             DomainPartition * const domain)
{
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ManagedGroup * const nodeManager = mesh->getNodeManager();
  localIndex const numNodes = nodeManager->size();

  view_rtype<r1_array> v_n = nodeManager->getData<r1_array>(keys::Velocity);
  view_rtype<r1_array> uhat  = nodeManager->getData<r1_array>(keys::IncrementalDisplacement);

  if( this->m_timeIntegrationOption == timeIntegrationOption::ImplicitDynamic )
  {
    view_rtype<r1_array> a_n = nodeManager->getData<r1_array>(keys::Acceleration);
    view_rtype<r1_array> vtilde    = nodeManager->getData<r1_array>(viewKeys.vTilde);
    view_rtype<r1_array> uhatTilde = nodeManager->getData<r1_array>(viewKeys.uhatTilde);
    real64 const newmarkGamma = this->getReference<real64>(viewKeys.newmarkGamma);
    real64 const newmarkBeta = this->getReference<real64>(viewKeys.newmarkBeta);

    for( auto a = 0 ; a < numNodes ; ++a )
    {
      for( int i=0 ; i<3 ; ++i )
      {
        //        real64 a_np1 = 4.0 / (dt*dt) * ( uhat[a][i] - uhatTilde[a][i]
        // );
        a_n[a][i] = 1.0 / ( newmarkBeta * dt*dt) * ( uhat[a][i] - uhatTilde[a][i] );
        v_n[a][i] = vtilde[a][i] + newmarkGamma * a_n[a][i] * dt;
      }
    }
  }
  else if( this->m_timeIntegrationOption == timeIntegrationOption::QuasiStatic && dt > 0.0)
  {
    for( auto a = 0 ; a < numNodes ; ++a )
    {
      for( int i=0 ; i<3 ; ++i )
      {
        v_n[a][i] = uhat[a][i] / dt;
      }
    }
  }
}

void SolidMechanics_LagrangianFEM::SetNumRowsAndTrilinosIndices( ManagedGroup * const nodeManager,
                                                                 localIndex & numLocalRows,
                                                                 globalIndex & numGlobalRows,
                                                                 localIndex_array& localIndices,
                                                                 localIndex offset )
{
//  dim =
// domain.m_feElementManager.m_ElementRegions.begin()->second.m_ElementDimension;
  int dim = 3;

  int n_mpi_processes;
  MPI_Comm_size( MPI_COMM_WORLD, &n_mpi_processes );

  int this_mpi_process = 0;
  MPI_Comm_rank( MPI_COMM_WORLD, &this_mpi_process );

  std::vector<int> gather(n_mpi_processes);

  int intNumLocalRows = static_cast<int>(numLocalRows);
  m_linearSolverWrapper.m_epetraComm.GatherAll( &intNumLocalRows,
                                                &gather.front(),
                                                1 );
  numLocalRows = intNumLocalRows;

  localIndex first_local_row = 0;
  numGlobalRows = 0;

  for( integer p=0 ; p<n_mpi_processes ; ++p)
  {
    numGlobalRows += gather[p];
    if(p<this_mpi_process)
      first_local_row += gather[p];
  }

  // create trilinos dof indexing

  globalIndex_array& trilinos_index = nodeManager->getReference<globalIndex_array>(viewKeys.trilinosIndex);
  integer_array& is_ghost       = nodeManager->getReference<integer_array>(viewKeys.ghostRank);


  trilinos_index = -1;

  integer local_count = 0;
  for(integer r=0 ; r<trilinos_index.size() ; ++r )
  {
    if(is_ghost[r] < 0)
    {
      trilinos_index[r] = first_local_row+local_count+offset;
      localIndices.push_back(trilinos_index[r]);
      local_count++;
    }
    else
    {
      trilinos_index[r] = -INT_MAX;
    }
  }

  assert(local_count == numLocalRows );

//  partition.SynchronizeFields(m_syncedFields, CommRegistry::lagrangeSolver02);

}


void SolidMechanics_LagrangianFEM :: SetupSystem ( DomainPartition * const domain,
                                                   EpetraBlockSystem * const blockSystem )
{
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  NodeManager * const nodeManager = mesh->getNodeManager();

//  using namespace BoundaryConditionFunctions;

  // determine the global/local degree of freedom distribution.

  //  const auto& kinematicSibling =
  // domain.m_feNodeManager.GetOneToOneMap("NodeToKinematicSibling");
  //  const auto& isNodeDead =
  // domain.m_feNodeManager.GetFieldData<int>("isDead");
  //  array1d<integer> kinematicSibling(domain.m_feNodeManager.m_numNodes);

  localIndex dim = 3;//domain.m_feElementManager.m_ElementRegions.begin()->second.m_ElementDimension;
  localIndex n_ghost_rows  = nodeManager->GetNumberOfGhosts();
  localIndex n_local_rows  = nodeManager->size()-n_ghost_rows;
  globalIndex n_global_rows = 0;

  localIndex_array displacementIndices;
  SetNumRowsAndTrilinosIndices( nodeManager,
                                n_local_rows,
                                n_global_rows,
                                displacementIndices,
                                0 );


  std::map<string, string_array > fieldNames;
  fieldNames["node"].push_back(viewKeys.trilinosIndex.Key());

  CommunicationTools::SynchronizeFields(fieldNames,
                              mesh,
                              domain->getReference< array1d<NeighborCommunicator> >( domain->viewKeys.neighbors ) );

  // create epetra map


  Epetra_Map * const rowMap = blockSystem->SetRowMap( BlockIDs::displacementBlock,
                                                      std::make_unique<Epetra_Map>( dim*n_global_rows,
                                                                                    dim*n_local_rows,
                                                                                    0,
                                                                                    m_linearSolverWrapper.m_epetraComm ) );

  Epetra_FECrsGraph * const sparsity = blockSystem->SetSparsity( BlockIDs::displacementBlock,
                                                                 BlockIDs::displacementBlock,
                                                                 std::make_unique<Epetra_FECrsGraph>(Copy,*rowMap,0) );



//  integer_array dummyDof;


//  const array1d<int>* isDetachedFromSolidMesh =
// domain.m_feNodeManager.GetFieldDataPointer<int> ("isDetachedFromSolidMesh");

//  if(domain.m_externalFaces.m_contactActive &&
// domain.m_contactManager.m_use_contact_search)
//  {
//    UpdateContactDataStructures(domain, false);
//  }

//  if(domain.m_externalFaces.m_contactActive)
//  {
//    InsertGlobalIndices( domain);
//    const bool planeStress =
//  SolidMechanics_LagrangianFEM::m_2dOption==SolidMechanics_LagrangianFEM::PlaneStress
// ;
//    if(domain.m_contactManager.m_nitsche_active)
//      domain.m_externalFaces.GetProjectionTensorAndWeightingAndStabilizationParameters(dim,
// planeStress, domain);
//  }

#ifdef SRC_INTERNAL2
  if(domain.m_xfemManager != nullptr)
  {
    SetSparsityPatternXFEM(domain);
  }
  else
#endif
  {
    SetSparsityPattern( domain, sparsity );
  }

  sparsity->GlobalAssemble();
  sparsity->OptimizeStorage();

  blockSystem->SetMatrix( BlockIDs::displacementBlock,
                          BlockIDs::displacementBlock,
                          std::make_unique<Epetra_FECrsMatrix>(Copy,*sparsity) );

  blockSystem->SetSolutionVector( BlockIDs::displacementBlock,
                                  std::make_unique<Epetra_FEVector>(*rowMap) );

  blockSystem->SetResidualVector( BlockIDs::displacementBlock,
                                  std::make_unique<Epetra_FEVector>(*rowMap) );

//  std::cout<<"m_sparsity->NumGlobalRows()     =
// "<<m_sparsity->NumGlobalRows()<<std::endl;
//  std::cout<<"m_sparsity->NumGlobalNonzeros() =
// "<<m_sparsity->NumGlobalNonzeros()<<std::endl;
}

void SolidMechanics_LagrangianFEM::SetSparsityPattern( DomainPartition const * const domain,
                                                       Epetra_FECrsGraph * const sparsity )
{
  int dim=3;
  MeshLevel const * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ManagedGroup const * const nodeManager = mesh->getNodeManager();

  globalIndex_array const & trilinos_index = nodeManager->getReference<globalIndex_array>(viewKeys.trilinosIndex);
  ElementRegionManager const * const elemManager = mesh->getElemManager();


  elemManager->forElementRegions([&](ElementRegion const * const elementRegion)
    {
      auto const & numMethodName = elementRegion->getData<string>(keys::numericalMethod);

      elementRegion->forCellBlocks([&](CellBlockSubRegion const * const cellBlock)
      {
        localIndex const numElems = cellBlock->size();
        array2d<localIndex> const & elemsToNodes = cellBlock->getWrapper<FixedOneToManyRelation>(cellBlock->viewKeys().nodeList)->reference();// getData<array2d<localIndex>>(keys::nodeList);
        localIndex const numNodesPerElement = elemsToNodes.size(1);

        globalIndex_array elementLocalDofIndex (dim*numNodesPerElement);

        for( localIndex k=0 ; k<numElems ; ++k )
        {
          arrayView1d<localIndex const> const localNodeIndices = elemsToNodes[k];

          for( localIndex a=0 ; a<numNodesPerElement ; ++a )
          {
            for(localIndex i=0 ; i<numNodesPerElement ; ++i)
            {
              for( int d=0 ; d<dim ; ++d )
              {
                elementLocalDofIndex[i*dim+d] = dim*static_cast<int>(trilinos_index[localNodeIndices[i]])+d;
              }
            }

            sparsity->InsertGlobalIndices(static_cast<int>(elementLocalDofIndex.size()),
                                          elementLocalDofIndex.data(),
                                          static_cast<int>(elementLocalDofIndex.size()),
                                          elementLocalDofIndex.data());
          }

        }
      });
    });
}



void SolidMechanics_LagrangianFEM::AssembleSystem ( DomainPartition * const  domain,
                                                EpetraBlockSystem * const blockSystem,
                                                real64 const time_n,
                                                real64 const dt )
{
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ManagedGroup * const nodeManager = mesh->getNodeManager();
  ConstitutiveManager  * const constitutiveManager = domain->GetGroup<ConstitutiveManager >(keys::ConstitutiveManager);
  ElementRegionManager * const elemManager = mesh->getElemManager();
  NumericalMethodsManager const * numericalMethodManager = domain->getParent()->GetGroup<NumericalMethodsManager>(keys::numericalMethodsManager);
  FiniteElementSpaceManager const * feSpaceManager = numericalMethodManager->GetGroup<FiniteElementSpaceManager>(keys::finiteElementSpaces);

  auto fluidPressure = elemManager->ConstructViewAccessor<real64_array>("fluidPressure");


  Epetra_FECrsMatrix * const matrix = blockSystem->GetMatrix( BlockIDs::displacementBlock,
                                                              BlockIDs::displacementBlock );
  Epetra_FEVector * const rhs = blockSystem->GetResidualVector( BlockIDs::displacementBlock );
  Epetra_FEVector * const solution = blockSystem->GetSolutionVector( BlockIDs::displacementBlock );

  matrix->Scale(0.0);
  rhs->Scale(0.0);

  real64 maxForce = 0.0;

  view_rtype_const<r1_array> disp = nodeManager->getData<r1_array>(keys::TotalDisplacement);
  view_rtype_const<r1_array> uhat = nodeManager->getData<r1_array>(keys::IncrementalDisplacement);
  view_rtype_const<r1_array> vel  = nodeManager->getData<r1_array>(keys::Velocity);


  view_rtype_const<r1_array> uhattilde = nullptr;
  view_rtype_const<r1_array> vtilde = nullptr;

  globalIndex_array const & trilinos_index = nodeManager->getReference<globalIndex_array>(viewKeys.trilinosIndex);

  static array1d< R1Tensor > u_local(8);
  static array1d< R1Tensor > uhat_local(8);
  static array1d< R1Tensor > vtilde_local(8);
  static array1d< R1Tensor > uhattilde_local(8);

  ElementRegionManager::ConstitutiveRelationAccessor<ConstitutiveBase>
  constitutiveRelations = elemManager->ConstructConstitutiveAccessor<ConstitutiveBase>(constitutiveManager);

  ElementRegionManager::MaterialViewAccessor< real64 > const
  density = elemManager->ConstructMaterialViewAccessor< real64 >( "density0",
                                                                  constitutiveManager );


  // begin region loop
  for( localIndex er=0 ; er<elemManager->numRegions() ; ++er )
  {
    ElementRegion * const elementRegion = elemManager->GetRegion(er);
    auto const & numMethodName = elementRegion->getData<string>(keys::numericalMethod);
    FiniteElementSpace const * feSpace = feSpaceManager->GetGroup<FiniteElementSpace>(numMethodName);

    for( localIndex esr=0 ; esr<elementRegion->numSubRegions() ; ++esr )
    {
      CellBlockSubRegion * const cellBlock = elementRegion->GetSubRegion(esr);

      multidimensionalArray::ManagedArray<R1Tensor, 3> const &
      dNdX = cellBlock->getReference< multidimensionalArray::ManagedArray<R1Tensor, 3> >(keys::dNdX);

      array2d<real64> const & detJ = cellBlock->getReference< array2d<real64> >(keys::detJ);

      array2d< localIndex > const & elemsToNodes = cellBlock->nodeList();
      localIndex const numNodesPerElement = elemsToNodes.size(1);

      u_local.resize(numNodesPerElement);
      uhat_local.resize(numNodesPerElement);
      vtilde_local.resize(numNodesPerElement);
      uhattilde_local.resize(numNodesPerElement);


      // space for element matrix and rhs
      int dim = 3;
      Epetra_LongLongSerialDenseVector  elementLocalDofIndex   (dim*static_cast<int>(numNodesPerElement));
      Epetra_SerialDenseVector     element_rhs     (dim*static_cast<int>(numNodesPerElement));
      Epetra_SerialDenseMatrix     element_matrix  (dim*static_cast<int>(numNodesPerElement),
                                                    dim*static_cast<int>(numNodesPerElement));
      Epetra_SerialDenseVector     element_dof_np1 (dim*static_cast<int>(numNodesPerElement));

      array1d<integer> const & elemGhostRank = cellBlock->m_ghostRank;

      GEOS_CXX_MARK_LOOP_BEGIN(elemLoop,elemLoop);

      for( localIndex k=0 ; k<cellBlock->size() ; ++k )
      {

        real64 stiffness[6][6];
        constitutiveRelations[er][esr][0]->GetStiffness( stiffness );

        if(elemGhostRank[k] < 0)
        {
          arrayView1d<localIndex const> const localNodeIndices = elemsToNodes[k];

          for( localIndex a=0 ; a<numNodesPerElement ; ++a)
          {

            localIndex localNodeIndex = localNodeIndices[a];

            for( int i=0 ; i<dim ; ++i )
            {
              elementLocalDofIndex[static_cast<int>(a)*dim+i] = dim*static_cast<int>(trilinos_index[localNodeIndex])+i;

              // TODO must add last solution estimate for this to be valid
              element_dof_np1(static_cast<int>(a)*dim+i) = disp[localNodeIndex][i];
            }
          }

          if( this->m_timeIntegrationOption == timeIntegrationOption::ImplicitDynamic )
          {
            CopyGlobalToLocal( localNodeIndices,
                               disp, uhat, vtilde, uhattilde,
                               u_local.data(), uhat_local.data(), vtilde_local.data(), uhattilde_local.data(),
                               numNodesPerElement );
          }
          else
          {
            CopyGlobalToLocal( localNodeIndices,
                               disp, uhat,
                               u_local.data(), uhat_local.data(),
                               numNodesPerElement );
          }

          R2SymTensor referenceStress;
          if( fluidPressure[er][esr].isValid() )
          {
            referenceStress.PlusIdentity( fluidPressure[er][esr][k] );
          }
          real64 maxElemForce = CalculateElementResidualAndDerivative( density[er][esr][0],
                                                                       feSpace->m_finiteElement,
                                                                       dNdX[k],
                                                                       detJ[k],
                                                                       &referenceStress,
                                                                       u_local,
                                                                       uhat_local,
                                                                       uhattilde_local,
                                                                       vtilde_local,
                                                                       dt,
                                                                       element_matrix,
                                                                       element_rhs,
                                                                       stiffness);


          //            if( maxElemForce > maxForce )
          //              maxForce = maxElemForce;

          matrix->SumIntoGlobalValues( elementLocalDofIndex,
                                       element_matrix);


          rhs->SumIntoGlobalValues( elementLocalDofIndex,
                                    element_rhs);
        }
      }
    }
  }


  // Global assemble
  matrix->GlobalAssemble(true);
  rhs->GlobalAssemble();

  if( verboseLevel() >= 2 )
  {
    matrix->Print(std::cout);
    rhs->Print(std::cout);
  }

 // return maxForce;
}

void
SolidMechanics_LagrangianFEM::
ApplyBoundaryConditions( DomainPartition * const domain,
                         systemSolverInterface::EpetraBlockSystem * const blockSystem,
                         real64 const time_n,
                         real64 const dt )
{

  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);

  FaceManager * const faceManager = mesh->getFaceManager();
  NodeManager * const nodeManager = mesh->getNodeManager();

  BoundaryConditionManager * bcManager = BoundaryConditionManager::get();
//  bcManager->ApplyBoundaryCondition( this, &SolidMechanics_LagrangianFEM::ForceBC,
//                                     nodeManager, keys::Force, time_n + dt, *blockSystem );

  bcManager->ApplyBoundaryCondition( time_n+dt,
                                     domain,
                                     "nodeManager",
                                     keys::Force,
                                     [&]( BoundaryConditionBase const * const bc,
                                          string const &,
                                          set<localIndex> const & targetSet,
                                          ManagedGroup * const targetGroup,
                                          string const fieldName )->void
  {
    bc->ApplyBoundaryConditionToSystem<BcAdd>( targetSet,
                                               time_n+dt,
                                               targetGroup,
                                               keys::TotalDisplacement, // TODO fix use of dummy name for
                                               viewKeys.trilinosIndex.Key(),
                                               3,
                                               blockSystem,
                                               BlockIDs::displacementBlock );
  });

  ApplyDisplacementBC_implicit( time_n + dt, *domain, *blockSystem );
//  bcManager->ApplyBoundaryCondition( this, &,
//                                     nodeManager, keys::TotalDisplacement, time_n + dt, *blockSystem );

  ApplyTractionBC( domain,
                   time_n+dt,
                   *blockSystem );


  Epetra_FECrsMatrix const * const matrix = blockSystem->GetMatrix( BlockIDs::displacementBlock,
                                                                    BlockIDs::displacementBlock );
  Epetra_FEVector const * const rhs = blockSystem->GetResidualVector( BlockIDs::displacementBlock );

  if( verboseLevel() >= 2 )
  {
    matrix->Print(std::cout);
    rhs->Print(std::cout);
  }

}

real64
SolidMechanics_LagrangianFEM::
CalculateResidualNorm(systemSolverInterface::EpetraBlockSystem const *const blockSystem, DomainPartition *const domain)
{

  Epetra_FEVector const * const
  residual = blockSystem->GetResidualVector( BlockIDs::displacementBlock );

  real64 localResidual = 0.0;
//  residual->Norm2(&scalarResidual);

  real64 * residualData = nullptr;
  int length;
  residual->ExtractView(&residualData,&length);
  for( localIndex i=0 ; i<length ; ++i )
  {
    localResidual += residualData[i]*residualData[i];
  }
  realT globalResidualNorm;
  MPI_Allreduce (&localResidual,&globalResidualNorm,1,MPI_DOUBLE,MPI_SUM ,MPI_COMM_WORLD);


  return sqrt(globalResidualNorm);

}

realT SolidMechanics_LagrangianFEM::CalculateElementResidualAndDerivative( real64 const density,
                                                                           FiniteElementBase const * const fe,
                                                                           const array_view<R1Tensor,2>& dNdX,
                                                                           const realT* const detJ,
                                                                           R2SymTensor const * const refStress,
                                                                           array1d<R1Tensor> const & u,
                                                                           array1d<R1Tensor> const & uhat,
                                                                           array1d<R1Tensor> const & uhattilde,
                                                                           array1d<R1Tensor> const & vtilde,
                                                                           realT const dt,
                                                                           Epetra_SerialDenseMatrix& dRdU,
                                                                           Epetra_SerialDenseVector& R,
                                                                           real64 c[6][6] )
{
  const integer dim = 3;
  realT maxForce = 0;
  realT amass = *this->getData<real64>(viewKeys.massDamping);
  realT astiff = *this->getData<real64>(viewKeys.stiffnessDamping);
  real64 const newmarkBeta = *(getData<real64>(viewKeys.newmarkBeta));
  real64 const newmarkGamma = *(getData<real64>(viewKeys.newmarkGamma));


//  if( LagrangeSolverBase::m_2dOption==LagrangeSolverBase::PlaneStress )
//  {
//    lambda = 2*lambda*G / ( lambda + 2*G );
//  }

  dRdU.Scale(0);
  R.Scale(0);

  Epetra_SerialDenseVector R_InertiaMassDamping(R);
  Epetra_SerialDenseMatrix dRdU_InertiaMassDamping(dRdU);

  Epetra_SerialDenseVector R_StiffnessDamping(R);
  Epetra_SerialDenseMatrix dRdU_StiffnessDamping(dRdU);

  dRdU_InertiaMassDamping.Scale(0);
  R_InertiaMassDamping.Scale(0);

  dRdU_StiffnessDamping.Scale(0);
  R_StiffnessDamping.Scale(0);

  R1Tensor dNdXa;
  R1Tensor dNdXb;


  for( integer q=0 ; q<fe->n_quadrature_points() ; ++q )
  {
    const realT detJq = detJ[q];
    std::vector<double> const & N = fe->values(q);

    for( integer a=0 ; a<fe->dofs_per_element() ; ++a )
    {
//      realT const * const dNdXa = dNdX(q,a).Data();
      dNdXa = dNdX(q,a);

      for( integer b=0 ; b<fe->dofs_per_element() ; ++b )
      {
//        realT const * const dNdXb = dNdX(q,b).Data();
        dNdXb = dNdX(q,b);

        if( dim==3 )
        {
          dRdU(a*dim+0,b*dim+0) -= ( c[0][0]*dNdXa[0]*dNdXb[0] + c[5][5]*dNdXa[1]*dNdXb[1] + c[4][4]*dNdXa[2]*dNdXb[2] ) * detJq;
          dRdU(a*dim+0,b*dim+1) -= ( c[5][5]*dNdXa[1]*dNdXb[0] + c[0][1]*dNdXa[0]*dNdXb[1] ) * detJq;
          dRdU(a*dim+0,b*dim+2) -= ( c[4][4]*dNdXa[2]*dNdXb[0] + c[0][2]*dNdXa[0]*dNdXb[2] ) * detJq;

          dRdU(a*dim+1,b*dim+0) -= ( c[0][1]*dNdXa[1]*dNdXb[0] + c[5][5]*dNdXa[0]*dNdXb[1] ) * detJq;
          dRdU(a*dim+1,b*dim+1) -= ( c[5][5]*dNdXa[0]*dNdXb[0] + c[1][1]*dNdXa[1]*dNdXb[1] + c[3][3]*dNdXa[2]*dNdXb[2] ) * detJq;
          dRdU(a*dim+1,b*dim+2) -= ( c[3][3]*dNdXa[2]*dNdXb[1] + c[1][2]*dNdXa[1]*dNdXb[2] ) * detJq;

          dRdU(a*dim+2,b*dim+0) -= ( c[0][2]*dNdXa[2]*dNdXb[0] + c[4][4]*dNdXa[0]*dNdXb[2] ) * detJq;
          dRdU(a*dim+2,b*dim+1) -= ( c[1][2]*dNdXa[2]*dNdXb[1] + c[3][3]*dNdXa[1]*dNdXb[2] ) * detJq;
          dRdU(a*dim+2,b*dim+2) -= ( c[4][4]*dNdXa[0]*dNdXb[0] + c[3][3]*dNdXa[1]*dNdXb[1] + c[2][2]*dNdXa[2]*dNdXb[2] ) * detJq;


          if( this->m_timeIntegrationOption == timeIntegrationOption::ImplicitDynamic )
          {

            double integrationFactor = density * N[a] * N[b] * detJq;
            double temp1 = ( amass * newmarkGamma/( newmarkBeta * dt ) + 1.0 / ( newmarkBeta * dt * dt ) )* integrationFactor;

            for( int i=0 ; i<dim ; ++i )
            {
              realT const acc = 1.0 / ( newmarkBeta * dt * dt ) * ( uhat[b][i] - uhattilde[b][i] );
              realT const velb = vtilde[b][i] + newmarkGamma/( newmarkBeta * dt ) *( uhat[b][i] - uhattilde[b][i] );

              dRdU_InertiaMassDamping(a*dim+i,b*dim+i) -= temp1;
              R_InertiaMassDamping(a*dim+i) -= ( amass * velb + acc ) * integrationFactor;
            }
          }
        }
        else if( dim==2 )
        {
//          dRdU(a*dim+0,b*dim+0) -= ( dNdXa[1]*dNdXb[1]*G +
// dNdXa[0]*dNdXb[0]*(2*G + lambda) ) * detJq;
//          dRdU(a*dim+0,b*dim+1) -= ( dNdXa[1]*dNdXb[0]*G +
// dNdXa[0]*dNdXb[1]*lambda ) * detJq;
//
//          dRdU(a*dim+1,b*dim+0) -= ( dNdXa[0]*dNdXb[1]*G +
// dNdXa[1]*dNdXb[0]*lambda ) * detJq;
//          dRdU(a*dim+1,b*dim+1) -= ( dNdXa[0]*dNdXb[0]*G +
// dNdXa[1]*dNdXb[1]*(2*G + lambda) ) * detJq;
        }
      }
    }
  }



  if( refStress!=nullptr )
  {
    R1Tensor temp;
    for( integer q=0 ; q<fe->n_quadrature_points() ; ++q )
    {
      const realT detJq = detJ[q];
      R2SymTensor stress0 = *refStress;
      stress0 *= detJq;
      for( integer a=0 ; a<fe->dofs_per_element() ; ++a )
      {
        dNdXa = dNdX(q,a);

        temp.AijBj(stress0,dNdXa);
        realT maxf = temp.MaxVal();
        if( maxf > maxForce )
        {
          maxForce = maxf;
        }

        R(a*dim+0) -= temp[0];
        R(a*dim+1) -= temp[1];
        R(a*dim+2) -= temp[2];
      }
    }
  }


// TODO It is simpler to do this...try it.
//  dRdU.Multiply(dof_np1,R);
  for( integer a=0 ; a<fe->dofs_per_element() ; ++a )
  {
    realT nodeForce = 0;
    for( integer b=0 ; b<fe->dofs_per_element() ; ++b )
    {
      for( int i=0 ; i<dim ; ++i )
      {
        for( int j=0 ; j<dim ; ++j )
        {
          R(a*dim+i) += dRdU(a*dim+i,b*dim+j) * u[b][j];
        }
      }

      if( this->m_timeIntegrationOption == timeIntegrationOption::ImplicitDynamic )
      {
        for( int i=0 ; i<dim ; ++i )
        {
          for( int j=0 ; j<dim ; ++j )
          {
            R_StiffnessDamping(a*dim+i) += astiff * dRdU(a*dim+i,b*dim+j) * ( vtilde[b][j] + newmarkGamma/(newmarkBeta * dt)*(uhat[b][j]-uhattilde[b][j]) );
          }
        }
      }

    }

    if (dim ==3)
    {
      nodeForce = std::max( std::max( R(a*dim+0), R(a*dim+1) ),  R(a*dim+2) );
    }
    else
    {
      nodeForce = std::max( R(a*dim+0), R(a*dim+1));
    }
//    std::cout<<"nodeForce["<<a<<"] = "<<nodeForce<<std::endl;
    if( fabs(nodeForce) > maxForce )
    {
      maxForce = fabs(nodeForce);
    }
  }


  if( this->m_timeIntegrationOption == timeIntegrationOption::ImplicitDynamic )
  {
    dRdU_StiffnessDamping = dRdU;
    dRdU_StiffnessDamping.Scale( astiff * newmarkGamma / ( newmarkBeta * dt ) );

    dRdU += dRdU_InertiaMassDamping;
    dRdU += dRdU_StiffnessDamping;
    R    += R_InertiaMassDamping;
    R    += R_StiffnessDamping;
  }


  return maxForce;
}



void SolidMechanics_LagrangianFEM::ApplySystemSolution( EpetraBlockSystem const * const blockSystem,
                                                        real64 const scalingFactor,
                                                        DomainPartition * const domain )
{
  NodeManager * const nodeManager = domain->getMeshBody(0)->getMeshLevel(0)->getNodeManager();

  Epetra_Map const * const rowMap        = blockSystem->GetRowMap( BlockIDs::displacementBlock );
  Epetra_FEVector const * const solution = blockSystem->GetSolutionVector( BlockIDs::displacementBlock );

  int solutionLength;
  double* local_solution = nullptr;
  solution->ExtractView(&local_solution,&solutionLength);


  view_rtype_const<globalIndex_array> trilinos_index = nodeManager->getData< globalIndex_array >(viewKeys.trilinosIndex);

  view_rtype<r1_array> X        = nodeManager->getData<r1_array>(nodeManager->viewKeys.referencePosition);
  view_rtype<r1_array> disp     = nodeManager->getData<r1_array>(keys::TotalDisplacement);
  view_rtype<r1_array> incdisp  = nodeManager->getData<r1_array>(keys::IncrementalDisplacement);

  localIndex const numNodes = nodeManager->size();
  realT maxpos = 0.0;
  realT maxdisp = 0.0;

  integer const dim = 3;

  for(integer r=0 ; r<numNodes ; ++r)
  {
    {
      for( int d=0 ; d<dim ; ++d )
      {
        int lid = rowMap->LID( static_cast<int>(dim*trilinos_index[r]) + d );

        if( lid >=0 )
        {
          incdisp[r][d] -= scalingFactor*local_solution[lid];
          disp[r][d] -= scalingFactor*local_solution[lid];
          maxpos = std::max( maxpos, fabs(X[r][d]+disp[r][d]) );
          maxdisp = std::max( maxdisp, fabs(disp[r][d]) );
        }
      }
    }
  }
//  m_maxDofVal = maxpos;
//  std::cout<<"Maximum DeltaDisplacement, Position = "<<maxinc<<",
// "<<maxpos<<", "<<maxinc/maxpos<<std::endl;

}

void SolidMechanics_LagrangianFEM::SolveSystem( EpetraBlockSystem * const blockSystem,
                                        SystemSolverParameters const * const params )
{
  Epetra_FEVector * const
  solution = blockSystem->GetSolutionVector( BlockIDs::displacementBlock );

  Epetra_FEVector * const
  residual = blockSystem->GetResidualVector( BlockIDs::displacementBlock );
//  residual->Scale(-1.0);

  solution->Scale(0.0);

  m_linearSolverWrapper.SolveSingleBlockSystem( blockSystem,
                                                 params,
                                                 BlockIDs::displacementBlock );

  if( verboseLevel() >= 2 )
  {
    solution->Print(std::cout);
  }

}

void SolidMechanics_LagrangianFEM::ResetStateToBeginningOfStep( DomainPartition * const domain )
{
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  NodeManager * const nodeManager = mesh->getNodeManager();

  view_rtype<r1_array> incdisp  = nodeManager->getData<r1_array>(keys::IncrementalDisplacement);

  // TODO need to finish this rewind
  FORALL_NODES( a, 0, nodeManager->size() )
  {
    incdisp[a] = 0.0;
  } END_FOR
}

REGISTER_CATALOG_ENTRY( SolverBase, SolidMechanics_LagrangianFEM, std::string const &, ManagedGroup * const )
} /* namespace ANST */
