/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
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

#include <sys/time.h>

#include "common/TimingMacros.hpp"

#include "dataRepository/ManagedGroup.hpp"
#include <common/DataTypes.hpp>
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/LinearElasticIsotropic.hpp"
#include "managers/NumericalMethodsManager.hpp"
#include "finiteElement/ElementLibrary/FiniteElement.h"
#include "finiteElement/FiniteElementDiscretizationManager.hpp"
#include "finiteElement/Kinematics.h"
#include "managers/BoundaryConditions/BoundaryConditionManager.hpp"

#include "codingUtilities/Utilities.hpp"

#include "managers/DomainPartition.hpp"
#include "meshUtilities/ComputationalGeometry.hpp"
#include "MPI_Communications/CommunicationTools.hpp"
#include "../../rajaInterface/GEOS_RAJA_Interface.hpp"
#include "MPI_Communications/NeighborCommunicator.hpp"


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
                arraySlice1d<R1Tensor> const & dNdX,
                real64 const& detJ,
                real64 const& detF,
                const R2Tensor& Finv,
                const localIndex numPoints,
                arraySlice1d<R1Tensor> & result)
{
  real64 const integrationFactor = detJ * detF;

  R2Tensor P;
  P.AijBkj( fieldvar,Finv);
  P *= integrationFactor;

  for( int a=0 ; a<numPoints ; ++a )  // loop through all shape functions in  element
  {
    result[a].minusAijBj( P, dNdX[a] );
  }
}

template< int N >
void Integrate( const R2SymTensor& fieldvar,
                arraySlice1d<R1Tensor const> const & dNdX,
                real64 const& detJ,
                real64 const& detF,
                const R2Tensor& Finv,
                arraySlice1d<R1Tensor> & result)
{
  real64 const integrationFactor = detJ * detF;

  R2Tensor P;
  P.AijBkj( fieldvar,Finv);
  P *= integrationFactor;

  for( int a=0 ; a<N ; ++a )  // loop through all shape functions in element
  {
    result[a].minusAijBj( P, dNdX[a] );
  }
}


inline void HughesWinget( R2Tensor &Rot, R2SymTensor & Dadt, R2Tensor const & G)
{

  real64 * restrict const Dadt_data = Dadt.Data();
  real64 * restrict const Rot_data = Rot.Data();
  real64 const * restrict const G_data = G.Data();


  //Dadt = 0.5*(G + GT);
  Dadt_data[0] = G_data[0];
  
  Dadt_data[1] = 0.5*(G_data[1] + G_data[3]);
  Dadt_data[2] = G_data[4];
  
  Dadt_data[3] = 0.5*(G_data[6] + G_data[2]);
  Dadt_data[4] = 0.5*(G_data[7] + G_data[5]);
  Dadt_data[5] = G_data[8];


  //Omega = 0.5*(G - GT);
  real64 const w12 = 0.5*(G_data[1] - G_data[3]);
  real64 const w13 = 0.5*(G_data[2] - G_data[6]);
  real64 const w23 = 0.5*(G_data[5] - G_data[7]);

  real64 const w12w12div4 = 0.25*w12*w12;
  real64 const w13w13div4 = 0.25*w13*w13;
  real64 const w23w23div4 = 0.25*w23*w23;
  real64 const w12w13div2 = 0.5*(w12*w13);
  real64 const w12w23div2 = 0.5*(w12*w23);
  real64 const w13w23div2 = 0.5*(w13*w23);
  real64 const invDetIplusOmega = 1.0 / ( 1 + ( w12w12div4 + w13w13div4 + w23w23div4 ) );

  Rot_data[0] = ( 1.0 + (-w12w12div4 - w13w13div4 + w23w23div4) ) * invDetIplusOmega;
  Rot_data[1] = ( w12 - w13w23div2 ) * invDetIplusOmega;
  Rot_data[2] = ( w13 + w12w23div2 ) * invDetIplusOmega;

  Rot_data[3] = (-w12 - w13w23div2 ) * invDetIplusOmega;
  Rot_data[4] = ( 1.0 + (-w12w12div4 + w13w13div4 - w23w23div4) ) * invDetIplusOmega;
  Rot_data[5] = ( w23 - w12w13div2 ) * invDetIplusOmega;

  Rot_data[6] = (-w13 + w12w23div2 ) * invDetIplusOmega;
  Rot_data[7] = (-w23 - w12w13div2 ) * invDetIplusOmega;
  Rot_data[8] = ( 1.0 + ( w12w12div4 - w13w13div4 - w23w23div4) ) * invDetIplusOmega;
  
  
}

inline void LinearElasticIsotropic_Kernel(R2SymTensor & Dadt, R2SymTensor & TotalStress, R2Tensor & Rot,
                                          localIndex i, real64 bulkModulus, real64 shearModulus,
                                          arrayView1d<real64> & meanStress,
                                          arrayView1d<R2SymTensor> & devStress)
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
  SolverBase( name, parent ),
  m_newmarkGamma(0.5),
  m_newmarkBeta(0.25),
  m_massDamping(0.0),
  m_stiffnessDamping(0.0),
  m_timeIntegrationOptionString(),
  m_timeIntegrationOption(timeIntegrationOption::ExplicitDynamic),
  m_useVelocityEstimateForQS(0),
  m_maxForce(0.0),
  m_elemsAttachedToSendOrReceiveNodes(),
  m_elemsNotAttachedToSendOrReceiveNodes(),
  m_sendOrRecieveNodes(),
  m_nonSendOrRecieveNodes(),
  m_icomm()
{
  // To generate the schema, multiple solvers of that use this command are constructed
  // Doing this can cause an error in the block setup, so move it to InitializePreSubGroups
  // getLinearSystemRepository()->SetBlockID( BlockIDs::displacementBlock, this->getName() );


  RegisterViewWrapper(viewKeyStruct::newmarkGammaString, &m_newmarkGamma, false )->
    setApplyDefaultValue(0.5)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Value of \\gamma in the Newmark Method for time integration");

  RegisterViewWrapper(viewKeyStruct::newmarkBetaString, &m_newmarkBeta, false )->
    setApplyDefaultValue(0.25)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Value of \\beta in the Newmark Method for time integration. "
          "This should be pow(newmarkGamma+0.5,2.0)/4.0 unless you know what you are doing.");

  RegisterViewWrapper(viewKeyStruct::massDampingString, &m_massDamping, false )->
    setApplyDefaultValue(0.0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Value of mass based damping coefficient in equations of motion. ");

  RegisterViewWrapper(viewKeyStruct::stiffnessDampingString, &m_stiffnessDamping, false )->
    setApplyDefaultValue(0.0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Value of stiffness based damping coefficient in equations of motion. ");

  RegisterViewWrapper(viewKeyStruct::timeIntegrationOptionString, &m_timeIntegrationOption, false )->
    setInputFlag(InputFlags::FALSE)->
    setDescription("Time integration method. Options are QuasiStatic, ImplicitDynamic, ExplicitDynamic");

  RegisterViewWrapper(viewKeyStruct::timeIntegrationOptionStringString, &m_timeIntegrationOptionString, false )->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Time integration method. Options are QuasiStatic, ImplicitDynamic, ExplicitDynamic");

  RegisterViewWrapper(viewKeyStruct::useVelocityEstimateForQSString, &m_useVelocityEstimateForQS, false )->
    setApplyDefaultValue(0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Flag to indicate the use of the incremental displacement from the previous step as an "
          "initial estimate for the incremental displacement of the current step.");
}

void SolidMechanics_LagrangianFEM::PostProcessInput()
{
  if( !m_timeIntegrationOptionString.empty() )
  {
    SetTimeIntegrationOption( m_timeIntegrationOptionString );
  }
}

SolidMechanics_LagrangianFEM::~SolidMechanics_LagrangianFEM()
{
  // TODO Auto-generated destructor stub
}


void SolidMechanics_LagrangianFEM::RegisterDataOnMesh( ManagedGroup * const MeshBodies )
{
  for( auto & mesh : MeshBodies->GetSubGroups() )
  {
    NodeManager * const nodes = mesh.second->group_cast<MeshBody*>()->getMeshLevel(0)->getNodeManager();
    nodes->RegisterViewWrapper<array1d<R1Tensor> >( viewKeyStruct::vTildeString );
    nodes->RegisterViewWrapper<array1d<R1Tensor> >( viewKeyStruct::uhatTildeString );
    nodes->RegisterViewWrapper<array1d<R1Tensor> >( keys::TotalDisplacement )->setPlotLevel(PlotLevel::LEVEL_0);
    nodes->RegisterViewWrapper<array1d<R1Tensor> >( keys::IncrementalDisplacement )->setPlotLevel(PlotLevel::LEVEL_2);
    nodes->RegisterViewWrapper<array1d<R1Tensor> >( keys::Velocity )->setPlotLevel(PlotLevel::LEVEL_0);
    nodes->RegisterViewWrapper<array1d<R1Tensor> >( keys::Acceleration )->setPlotLevel(PlotLevel::LEVEL_1);
    nodes->RegisterViewWrapper<array1d<real64> >( keys::Mass )->setPlotLevel(PlotLevel::LEVEL_0);
    nodes->RegisterViewWrapper<array1d<globalIndex> >( viewKeyStruct::trilinosIndexString )->setPlotLevel(PlotLevel::LEVEL_1);

  }
}


void SolidMechanics_LagrangianFEM::InitializePreSubGroups(ManagedGroup * const rootGroup)
{
  getLinearSystemRepository()->SetBlockID( BlockIDs::displacementBlock, this->getName() );
}


void SolidMechanics_LagrangianFEM::InitializePostInitialConditions_PreSubGroups( ManagedGroup * const problemManager )
{
  DomainPartition * domain = problemManager->GetGroup<DomainPartition>(keys::domain);
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);

  NodeManager * const nodes = mesh->getNodeManager();
  FaceManager * const faceManager = mesh->getFaceManager();


  ElementRegionManager * elementRegionManager = mesh->getElemManager();
  ConstitutiveManager * constitutiveManager = domain->GetGroup<ConstitutiveManager >(keys::ConstitutiveManager);

  arrayView1d<real64> & mass = nodes->getReference<array1d<real64>>(keys::Mass);

  ElementRegionManager::MaterialViewAccessor< arrayView2d<real64> > rho =
    elementRegionManager->ConstructMaterialViewAccessor< array2d<real64>, arrayView2d<real64> >("density", constitutiveManager);

  m_elemsAttachedToSendOrReceiveNodes.resize( elementRegionManager->numRegions() );
  m_elemsNotAttachedToSendOrReceiveNodes.resize( elementRegionManager->numRegions() );

  for( localIndex er=0 ; er<elementRegionManager->numRegions() ; ++er )
  {
    ElementRegion const * const elemRegion = elementRegionManager->GetRegion(er);
    m_elemsAttachedToSendOrReceiveNodes[er].resize( elemRegion->numSubRegions() );
    m_elemsNotAttachedToSendOrReceiveNodes[er].resize( elemRegion->numSubRegions() );

    for( localIndex esr=0 ; esr<elemRegion->numSubRegions() ; ++esr )
    {
      CellBlockSubRegion const * const cellBlock = elemRegion->GetSubRegion(esr);

      arrayView2d<real64> const & detJ = cellBlock->getReference< array2d<real64> >(keys::detJ);

      arrayView2d<localIndex> const & elemsToNodes =
        cellBlock->getWrapper<FixedOneToManyRelation>(cellBlock->viewKeys().nodeList)->reference();// getReference<array2d<localIndex>>(keys::nodeList);

      for( localIndex k=0 ; k < elemsToNodes.size(0) ; ++k )
      {
        for( localIndex q=0 ; q< elemsToNodes.size(1) ; ++q )
        {
          mass[elemsToNodes[k][q]] += rho[er][esr][0][k][q] * detJ[k][q];
        }

        bool isAttachedToGhostNode = false;
        for( localIndex a=0 ; a<cellBlock->numNodesPerElement() ; ++a )
        {
          if( nodes->GhostRank()[elemsToNodes[k][a]] >= -1 )
          {
            isAttachedToGhostNode = true;
            m_sendOrRecieveNodes.insert( elemsToNodes[k][a] );
          }
          else
          {
            m_nonSendOrRecieveNodes.insert( elemsToNodes[k][a] );
          }
        }

        if( isAttachedToGhostNode )
        {
          m_elemsAttachedToSendOrReceiveNodes[er][esr].insert(k);
        }
        else
        {
          m_elemsNotAttachedToSendOrReceiveNodes[er][esr].insert(k);
        }
      }
    }
  }
}

real64 SolidMechanics_LagrangianFEM::SolverStep( real64 const& time_n,
                                             real64 const& dt,
                                             const int cycleNumber,
                                             DomainPartition * domain )
{
  real64 dtReturn = dt;

  SolverBase * const surfaceGenerator =  this->getParent()->GetGroup<SolverBase>("SurfaceGen");

  if( m_timeIntegrationOption == timeIntegrationOption::ExplicitDynamic )
  {
    dtReturn = ExplicitStep( time_n, dt, cycleNumber, ManagedGroup::group_cast<DomainPartition*>(domain) );
  }
  else if( m_timeIntegrationOption == timeIntegrationOption::ImplicitDynamic ||
           m_timeIntegrationOption == timeIntegrationOption::QuasiStatic )
  {
    ImplicitStepSetup( time_n, dt, domain, getLinearSystemRepository() );

    dtReturn = NonlinearImplicitStep( time_n, dt, cycleNumber, domain->group_cast<DomainPartition*>(),
                                      getLinearSystemRepository() );

    ImplicitStepComplete( time_n, dt,  domain );

    if( surfaceGenerator!=nullptr )
    {
      surfaceGenerator->SolverStep( time_n, dt, cycleNumber, domain );
    }

  }

  return dtReturn;
}

real64 SolidMechanics_LagrangianFEM::ExplicitStep( real64 const& time_n,
                                                   real64 const& dt,
                                                   const int cycleNumber,
                                                   DomainPartition * const domain )
{
  GEOSX_MARK_FUNCTION;

  static real64 minTimes[10] = {1.0e9,1.0e9,1.0e9,1.0e9,1.0e9,1.0e9,1.0e9,1.0e9,1.0e9,1.0e9};
  static real64 maxTimes[10] = {0.0};

  GEOSX_GET_TIME( t0 );

  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  NodeManager * const nodes = mesh->getNodeManager();
  ElementRegionManager * elemManager = mesh->getElemManager();
  NumericalMethodsManager const * numericalMethodManager = domain->getParent()->GetGroup<NumericalMethodsManager>(keys::numericalMethodsManager);
  FiniteElementDiscretizationManager const * feDiscretizationManager = numericalMethodManager->GetGroup<FiniteElementDiscretizationManager>(keys::finiteElementDiscretizations);
  ConstitutiveManager * constitutiveManager = domain->GetGroup<ConstitutiveManager >(keys::ConstitutiveManager);

  BoundaryConditionManager * bcManager = BoundaryConditionManager::get();
  localIndex const numNodes = nodes->size();

  arrayView1d<R1Tensor> const & X = nodes->getReference<array1d<R1Tensor>>(nodes->viewKeys.referencePosition);
  arrayView1d<real64> const & mass = nodes->getReference<array1d<real64>>(keys::Mass);
  arrayView1d<R1Tensor> & vel = nodes->getReference<array1d<R1Tensor>>(keys::Velocity);

  arrayView1d<R1Tensor> & u = nodes->getReference<array1d<R1Tensor>>(keys::TotalDisplacement);
  arrayView1d<R1Tensor> & uhat = nodes->getReference<array1d<R1Tensor>>(keys::IncrementalDisplacement);
  arrayView1d<R1Tensor> & acc = nodes->getReference<array1d<R1Tensor>>(keys::Acceleration);

  array1d<NeighborCommunicator> & neighbors = domain->getReference< array1d<NeighborCommunicator> >( domain->viewKeys.neighbors );
  std::map<string, string_array > fieldNames;
  fieldNames["node"].push_back("Velocity");

  CommunicationTools::SynchronizePackSendRecvSizes( fieldNames, mesh, neighbors, m_icomm );

  bcManager->ApplyBoundaryConditionToField( time_n, domain, "nodeManager", keys::Acceleration );

  GEOSX_MARK_BEGIN(firstVelocityUpdate);

  //3: v^{n+1/2} = v^{n} + a^{n} dt/2
  SolidMechanicsLagrangianFEMKernels::OnePoint( acc, vel, dt/2, numNodes );
  GEOSX_MARK_END(firstVelocityUpdate);

  bcManager->ApplyBoundaryConditionToField( time_n, domain, "nodeManager", keys::Velocity );

  //4. x^{n+1} = x^{n} + v^{n+{1}/{2}} dt (x is displacement)
  SolidMechanicsLagrangianFEMKernels::OnePoint( vel, uhat, u, dt, numNodes );


  bcManager->ApplyBoundaryConditionToField( time_n + dt, domain, "nodeManager", keys::TotalDisplacement,
    [&]( BoundaryConditionBase const * const bc, set<localIndex> const & targetSet )->void
    {
      integer const component = bc->GetComponent();
      for( auto const a : targetSet )
      {
        uhat[a][component] = u[a][component] - u[a][component];
        vel[a][component]  = uhat[a][component] / dt;
      }
    }
  );

  forall_in_range(0, numNodes, GEOSX_LAMBDA (localIndex a) mutable
  {
    acc[a] = 0;
  });

  ElementRegionManager::MaterialViewAccessor< arrayView2d<real64> > meanStress =
    elemManager->ConstructMaterialViewAccessor< array2d<real64>, arrayView2d<real64> >("MeanStress", constitutiveManager);

  ElementRegionManager::MaterialViewAccessor< arrayView2d<R2SymTensor> > const devStress =
    elemManager->ConstructMaterialViewAccessor< array2d<R2SymTensor>, arrayView2d<R2SymTensor> >("DeviatorStress", constitutiveManager);

  ElementRegionManager::ConstitutiveRelationAccessor<ConstitutiveBase> constitutiveRelations =
    elemManager->ConstructConstitutiveAccessor<ConstitutiveBase>(constitutiveManager);

  GEOSX_GET_TIME( t1 );

  //Step 5. Calculate deformation input to constitutive model and update state to
  // Q^{n+1}
  for( localIndex er=0 ; er<elemManager->numRegions() ; ++er )
  {
    ElementRegion * const elementRegion = elemManager->GetRegion(er);
    FiniteElementDiscretization const * feDiscretization = feDiscretizationManager->GetGroup<FiniteElementDiscretization>(m_discretizationName);

    for( localIndex esr=0 ; esr<elementRegion->numSubRegions() ; ++esr )
    {
      CellBlockSubRegion * const cellBlock = elementRegion->GetSubRegion(esr);

      arrayView3d< R1Tensor > const & dNdX = cellBlock->getReference< array3d< R1Tensor > >(keys::dNdX);

      arrayView2d<real64> const & detJ = cellBlock->getReference< array2d<real64> >(keys::detJ);

      arrayView2d<localIndex> const & elemsToNodes = cellBlock->getWrapper<FixedOneToManyRelation>(cellBlock->viewKeys().nodeList)->reference();

      localIndex const numNodesPerElement = elemsToNodes.size(1);

      localIndex const numQuadraturePoints = feDiscretization->m_finiteElement->n_quadrature_points();

      GEOSX_MARK_BEGIN(externalElemsLoop);

      ElementKernelSelector( er, 
                             esr,
                             this->m_elemsAttachedToSendOrReceiveNodes[er][esr],
                             elemsToNodes,
                             dNdX,
                             detJ,
                             u,
                             uhat,
                             acc,
                             constitutiveRelations,
                             meanStress,
                             devStress,
                             dt,
                             numNodesPerElement,
                             numQuadraturePoints );

      GEOSX_MARK_END(externalElemsLoop);

    } //Element Region

  } //Element Manager

  //Compute Force : Point-wise computations
  forall_in_set(m_sendOrRecieveNodes.data(), m_sendOrRecieveNodes.size(), GEOSX_LAMBDA (localIndex a) mutable
  {
    acc[a] /=mass[a];
  });

  // apply this over a set
  SolidMechanicsLagrangianFEMKernels::OnePoint( acc, vel, dt / 2, m_sendOrRecieveNodes.data(), m_sendOrRecieveNodes.size() );

  bcManager->ApplyBoundaryConditionToField( time_n, domain, "nodeManager", keys::Velocity );

  GEOSX_GET_TIME( t2 );
  CommunicationTools::SynchronizePackSendRecv( fieldNames, mesh, neighbors, m_icomm );
  GEOSX_GET_TIME( t3 );

  for( localIndex er=0 ; er<elemManager->numRegions() ; ++er )
  {
    ElementRegion * const elementRegion = elemManager->GetRegion(er);

    FiniteElementDiscretization const * feDiscretization = feDiscretizationManager->GetGroup<FiniteElementDiscretization>(m_discretizationName);

    for( localIndex esr=0 ; esr<elementRegion->numSubRegions() ; ++esr )
    {
      CellBlockSubRegion * const cellBlock = elementRegion->GetSubRegion(esr);

      arrayView3d< R1Tensor > const & dNdX = cellBlock->getReference< array3d< R1Tensor > >(keys::dNdX);

      arrayView2d<real64> const & detJ = cellBlock->getReference< array2d<real64> >(keys::detJ);

      arrayView2d<localIndex> const & elemsToNodes = cellBlock->getWrapper<FixedOneToManyRelation>(cellBlock->viewKeys().nodeList)->reference();

      localIndex const numNodesPerElement = elemsToNodes.size(1);

      localIndex const numQuadraturePoints = feDiscretization->m_finiteElement->n_quadrature_points();

      GEOSX_MARK_BEGIN(internalElemsLoop);

      ElementKernelSelector( er,
                             esr,
                             this->m_elemsNotAttachedToSendOrReceiveNodes[er][esr],
                             elemsToNodes,
                             dNdX,
                             detJ,
                             u,
                             uhat,
                             acc,
                             constitutiveRelations,
                             meanStress,
                             devStress,
                             dt,
                             numNodesPerElement,
                             numQuadraturePoints );

      GEOSX_MARK_END(internalElemsLoop);

    } //Element Region

  } //Element Manager

  GEOSX_GET_TIME( t4 );

  //Compute Force : Point-wise computations
  forall_in_set(m_nonSendOrRecieveNodes.data(), m_nonSendOrRecieveNodes.size(), GEOSX_LAMBDA (localIndex a)
  {
    acc[a] /=mass[a];
  });

  // apply this over a set
  SolidMechanicsLagrangianFEMKernels::OnePoint( acc, vel, dt / 2, m_nonSendOrRecieveNodes.data(), m_nonSendOrRecieveNodes.size());

  bcManager->ApplyBoundaryConditionToField( time_n, domain, "nodeManager", keys::Velocity );

  CommunicationTools::SynchronizeUnpack( mesh, neighbors, m_icomm );


  GEOSX_GET_TIME( tf );

#ifdef GEOSX_USE_TIMERS
  GEOSX_MARK_BEGIN("MPI_Barrier");
  MPI_Barrier(MPI_COMM_GEOSX);
  GEOSX_MARK_END("MPI_Barrier");

  minTimes[0] = std::min( minTimes[0], t2-t1 );
  minTimes[1] = std::min( minTimes[1], t4-t3 );
  minTimes[2] = std::min( minTimes[2], tf-t0 );


  maxTimes[0] = std::max( maxTimes[0], t2-t1 );
  maxTimes[1] = std::max( maxTimes[1], t4-t3 );
  maxTimes[2] = std::max( maxTimes[2], tf-t0 );


  GEOS_LOG_RANK( "     outer loop: "<< (t2-t1)<<", "<<minTimes[0]<<", "<<maxTimes[0] );
  GEOS_LOG_RANK( "     inner loop: "<< (t4-t3)<<", "<<minTimes[1]<<", "<<maxTimes[1] );
  GEOS_LOG_RANK( "     total time: "<< (tf-t0)<<", "<<minTimes[2]<<", "<<maxTimes[2] );
#endif

  return dt;
}


real64 SolidMechanics_LagrangianFEM::ElementKernelSelector( localIndex const er,
                                                            localIndex const esr,
                                                            set<localIndex> const & elementList,
                                                            arrayView2d<localIndex> const & elemsToNodes,
                                                            arrayView3d< R1Tensor > const & dNdX,
                                                            arrayView2d<real64> const & detJ,
                                                            arrayView1d<R1Tensor> const & u,
                                                            arrayView1d<R1Tensor> const & uhat,
                                                            arrayView1d<R1Tensor> & acc,
                                                            ElementRegionManager::ConstitutiveRelationAccessor<constitutive::ConstitutiveBase>& constitutiveRelations,
                                                            ElementRegionManager::MaterialViewAccessor< arrayView2d<real64> > const & meanStress,
                                                            ElementRegionManager::MaterialViewAccessor< arrayView2d<R2SymTensor> > const & devStress,
                                                            real64 const dt,
                                                            localIndex NUM_NODES_PER_ELEM,
                                                            localIndex NUM_QUADRATURE_POINTS )
{

  real64 rval = 0;

  if( NUM_NODES_PER_ELEM==8 && NUM_QUADRATURE_POINTS==8 )
  {
    rval =
    ExplicitElementKernel<8,8>( er,
                                esr,
                                elementList,
                                elemsToNodes,
                                dNdX,
                                detJ,
                                u,
                                uhat,
                                acc,
                                constitutiveRelations,
                                meanStress,
                                devStress,
                                dt );
  }

  return rval;
}

template< localIndex NUM_NODES_PER_ELEM, localIndex NUM_QUADRATURE_POINTS >
real64 SolidMechanics_LagrangianFEM::ExplicitElementKernel( localIndex const er,
                                                            localIndex const esr,
                                                            set<localIndex> const & elementList,
                                                            arrayView2d<localIndex> const & elemsToNodes,
                                                            arrayView3d< R1Tensor > const & dNdX,
                                                            arrayView2d<real64> const & detJ,
                                                            arrayView1d<R1Tensor> const & u,
                                                            arrayView1d<R1Tensor> const & uhat,
                                                            arrayView1d<R1Tensor> & acc,
                                                            ElementRegionManager::ConstitutiveRelationAccessor<ConstitutiveBase> constitutiveRelations,
                                                            ElementRegionManager::MaterialViewAccessor< arrayView2d<real64> > const & meanStress,
                                                            ElementRegionManager::MaterialViewAccessor< arrayView2d<R2SymTensor> > const & devStress,
                                                            real64 const dt )
{

  ConstitutiveBase::UpdateFunctionPointer update = constitutiveRelations[er][esr][0]->GetStateUpdateFunctionPointer();
  void * data = nullptr;
  constitutiveRelations[er][esr][0]->SetParamStatePointers( data );
  forall_in_set<elemPolicy>( elementList.data(),
                                   elementList.size(),
                                   GEOSX_LAMBDA ( localIndex k) mutable
  {
    r1_array uhat_local( NUM_NODES_PER_ELEM );
    r1_array u_local( NUM_NODES_PER_ELEM );
    r1_array f_local( NUM_NODES_PER_ELEM );

    CopyGlobalToLocal<R1Tensor,NUM_NODES_PER_ELEM>( elemsToNodes[k],
                                                    u, uhat,
                                                    u_local, uhat_local );

    //Compute Quadrature
    for( localIndex q = 0 ; q<NUM_QUADRATURE_POINTS ; ++q)
    {

      R2Tensor dUhatdX, dUdX;
      CalculateGradients<NUM_NODES_PER_ELEM>( dUhatdX, dUdX, uhat_local, u_local, dNdX[k][q]);

      R2Tensor F,Ldt, Finv;

      // calculate du/dX
      F = dUhatdX;
      F *= 0.5;
      F += dUdX;
      F.PlusIdentity(1.0);
      Finv.Inverse(F);

      // chain rule: calculate dv/du = dv/dX * dX/du
      Ldt.AijBjk(dUhatdX, Finv);

      // calculate gradient (end of step)
      F = dUhatdX;
      F += dUdX;
      F.PlusIdentity(1.0);
      real64 detF = F.Det();
      Finv.Inverse(F);


      R2Tensor Rot;
      R2SymTensor Dadt;
      HughesWinget(Rot, Dadt, Ldt);

      constitutiveRelations[er][esr][0]->StateUpdatePoint( Dadt, Rot, k, q, 0);

      R2SymTensor TotalStress;
      TotalStress = devStress[er][esr][0][k][q];
      TotalStress.PlusIdentity( meanStress[er][esr][0][k][q] );

      Integrate<NUM_NODES_PER_ELEM>( TotalStress, dNdX[k][q], detJ[k][q], detF, Finv, f_local );
    }//quadrature loop


    AddLocalToGlobal( elemsToNodes[k], f_local, acc, NUM_NODES_PER_ELEM );
  });

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
                                                 viewKeyStruct::trilinosIndexString,
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

  arrayView1d<globalIndex> const & blockLocalDofNumber =
    nodeManager->getReference<globalIndex_array>(solidMechanicsViewKeys.trilinosIndex);

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
    string const & functionName = bc->getReference<string>( BoundaryConditionBase::viewKeyStruct::functionNameString);

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
      assert( function!=nullptr);

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
              nodeDOF[a] = 3*blockLocalDofNumber[facesToNodes[kf][a]]+component;
              nodeRHS[a] = value * faceArea[kf] / numNodes;
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
              nodeDOF[a] = 3*blockLocalDofNumber[facesToNodes[kf][a]]+component;
              nodeRHS[a] = result[kf] * faceArea[kf] / numNodes;
            }
            rhs->SumIntoGlobalValues( integer_conversion<int>(nodeDOF.size()), nodeDOF.data(), nodeRHS.data() );
          }
      }
    }
  });
}

void
SolidMechanics_LagrangianFEM::
ImplicitStepSetup( real64 const& time_n, real64 const& dt, DomainPartition * const domain,
                   systemSolverInterface::EpetraBlockSystem * const blockSystem )
{
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ManagedGroup * const nodeManager = mesh->getNodeManager();

  if( this->m_timeIntegrationOption == timeIntegrationOption::ImplicitDynamic )
  {
    arrayView1d<R1Tensor> const& v_n = nodeManager->getReference<array1d<R1Tensor>>(keys::Velocity);
    arrayView1d<R1Tensor> const& a_n = nodeManager->getReference<array1d<R1Tensor>>(keys::Acceleration);
    arrayView1d<R1Tensor>& vtilde   = nodeManager->getReference<array1d<R1Tensor>>(solidMechanicsViewKeys.vTilde);
    arrayView1d<R1Tensor>& uhatTilde   = nodeManager->getReference<array1d<R1Tensor>>(solidMechanicsViewKeys.uhatTilde);

    arrayView1d<R1Tensor>& uhat  = nodeManager->getReference<array1d<R1Tensor>>(keys::IncrementalDisplacement);
    arrayView1d<R1Tensor>& disp = nodeManager->getReference<array1d<R1Tensor>>(keys::TotalDisplacement);

    localIndex const numNodes = nodeManager->size();
    real64 const newmarkGamma = this->getReference<real64>(solidMechanicsViewKeys.newmarkGamma);
    real64 const newmarkBeta = this->getReference<real64>(solidMechanicsViewKeys.newmarkBeta);

    for( localIndex a = 0 ; a < numNodes ; ++a )
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

    arrayView1d<R1Tensor>& uhat = nodeManager->getReference<array1d<R1Tensor>>(keys::IncrementalDisplacement);
    integer const useVelocityEstimateForQS = this->getReference<integer>(solidMechanicsViewKeys.useVelocityEstimateForQS);
    localIndex const numNodes = nodeManager->size();

    if( useVelocityEstimateForQS==1 )
    {
      arrayView1d<R1Tensor> const& v_n = nodeManager->getReference<array1d<R1Tensor>>(keys::Velocity);
      arrayView1d<R1Tensor>& disp = nodeManager->getReference<array1d<R1Tensor>>(keys::TotalDisplacement);

      for( localIndex a = 0 ; a < numNodes ; ++a )
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
      for( localIndex a = 0 ; a < numNodes ; ++a )
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

  r1_array& v_n = nodeManager->getReference<r1_array>(keys::Velocity);
  r1_array& uhat  = nodeManager->getReference<r1_array>(keys::IncrementalDisplacement);

  if( this->m_timeIntegrationOption == timeIntegrationOption::ImplicitDynamic )
  {
    r1_array& a_n = nodeManager->getReference<r1_array>(keys::Acceleration);
    r1_array& vtilde    = nodeManager->getReference<r1_array>(solidMechanicsViewKeys.vTilde);
    r1_array& uhatTilde = nodeManager->getReference<r1_array>(solidMechanicsViewKeys.uhatTilde);
    real64 const newmarkGamma = this->getReference<real64>(solidMechanicsViewKeys.newmarkGamma);
    real64 const newmarkBeta = this->getReference<real64>(solidMechanicsViewKeys.newmarkBeta);

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
  MPI_Comm_size( MPI_COMM_GEOSX, &n_mpi_processes );

  int this_mpi_process = 0;
  MPI_Comm_rank( MPI_COMM_GEOSX, &this_mpi_process );

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

  globalIndex_array& trilinos_index = nodeManager->getReference<globalIndex_array>(solidMechanicsViewKeys.trilinosIndex);
  integer_array& is_ghost       = nodeManager->getReference<integer_array>( ObjectManagerBase::viewKeyStruct::ghostRankString);

  trilinos_index = -1;

  integer local_count = 0;
  for(integer r=0 ; r<trilinos_index.size() ; ++r )
  {
    if(is_ghost[r] < 0)
    {
      trilinos_index[r] = first_local_row+local_count+offset;
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

  localIndex dim = 3;
  localIndex n_ghost_rows  = nodeManager->GetNumberOfGhosts();
  localIndex n_local_rows  = nodeManager->size()-n_ghost_rows;
  globalIndex n_global_rows = 0;

  localIndex_array displacementIndices;
  SetNumRowsAndTrilinosIndices( nodeManager, n_local_rows, n_global_rows, displacementIndices, 0 );

  std::map<string, string_array > fieldNames;
  fieldNames["node"].push_back(viewKeyStruct::trilinosIndexString);

  CommunicationTools::SynchronizeFields( fieldNames, mesh,
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

  SetSparsityPattern( domain, sparsity );

  sparsity->GlobalAssemble();
  sparsity->OptimizeStorage();

  blockSystem->SetMatrix( BlockIDs::displacementBlock,
                          BlockIDs::displacementBlock,
                          std::make_unique<Epetra_FECrsMatrix>(Copy,*sparsity) );

  blockSystem->SetSolutionVector( BlockIDs::displacementBlock,
                                  std::make_unique<Epetra_FEVector>(*rowMap) );

  blockSystem->SetResidualVector( BlockIDs::displacementBlock,
                                  std::make_unique<Epetra_FEVector>(*rowMap) );
}

void SolidMechanics_LagrangianFEM::SetSparsityPattern( DomainPartition const * const domain,
                                                       Epetra_FECrsGraph * const sparsity )
{
  int dim=3;
  MeshLevel const * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ManagedGroup const * const nodeManager = mesh->getNodeManager();

  arrayView1d<globalIndex> const & trilinos_index = nodeManager->getReference<array1d<globalIndex>>(solidMechanicsViewKeys.trilinosIndex);

  ElementRegionManager const * const elemManager = mesh->getElemManager();

  for( localIndex er=0 ; er<elemManager->numRegions() ; ++er )
  {
    ElementRegion const * const elementRegion = elemManager->GetRegion(er);

    for( localIndex esr=0 ; esr<elementRegion->numSubRegions() ; ++esr )
    {

      CellBlockSubRegion const * const cellBlock = elementRegion->GetSubRegion(esr);
      localIndex const numElems = cellBlock->size();
      arrayView2d<localIndex> const & elemsToNodes = cellBlock->getWrapper<FixedOneToManyRelation>(cellBlock->viewKeys().nodeList)->reference();// getReference<array2d<localIndex>>(keys::nodeList);
      localIndex const numNodesPerElement = elemsToNodes.size(1);

      globalIndex_array elementLocalDofIndex(dim * numNodesPerElement);

      for( localIndex k=0 ; k<numElems ; ++k )
      {
        for( localIndex a=0 ; a<numNodesPerElement ; ++a )
        {
          for(localIndex i=0 ; i<numNodesPerElement ; ++i)
          {
            for( int d=0 ; d<dim ; ++d )
            {
              elementLocalDofIndex[i * dim + d] = dim * trilinos_index[elemsToNodes[k][i]] + d;
            }
          }

          sparsity->InsertGlobalIndices(static_cast<int>(elementLocalDofIndex.size()),
                                        elementLocalDofIndex.data(),
                                        static_cast<int>(elementLocalDofIndex.size()),
                                        elementLocalDofIndex.data());
        }
      }
    }
  }
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
  FiniteElementDiscretizationManager const * feDiscretizationManager = numericalMethodManager->GetGroup<FiniteElementDiscretizationManager>(keys::finiteElementDiscretizations);

  ElementRegionManager::MaterialViewAccessor<real64> const biotCoefficient =
    elemManager->ConstructMaterialViewAccessor<real64>( "BiotCoefficient", constitutiveManager);

  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> const fluidPres =
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>("pressure");

  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> const dPres =
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>("deltaPressure");

  Epetra_FECrsMatrix * const matrix = blockSystem->GetMatrix( BlockIDs::displacementBlock,
                                                              BlockIDs::displacementBlock );
  Epetra_FEVector * const rhs = blockSystem->GetResidualVector( BlockIDs::displacementBlock );
  Epetra_FEVector * const solution = blockSystem->GetSolutionVector( BlockIDs::displacementBlock );

  matrix->Scale(0.0);
  rhs->Scale(0.0);

  r1_array const& disp = nodeManager->getReference<r1_array>(keys::TotalDisplacement);
  r1_array const& uhat = nodeManager->getReference<r1_array>(keys::IncrementalDisplacement);
  r1_array const& vel  = nodeManager->getReference<r1_array>(keys::Velocity);

  r1_array const uhattilde;
  r1_array const vtilde;

  globalIndex_array const & trilinos_index = nodeManager->getReference<globalIndex_array>(solidMechanicsViewKeys.trilinosIndex);

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

    FiniteElementDiscretization const * feDiscretization = feDiscretizationManager->GetGroup<FiniteElementDiscretization>(m_discretizationName);

    for( localIndex esr=0 ; esr<elementRegion->numSubRegions() ; ++esr )
    {
      CellBlockSubRegion * const cellBlock = elementRegion->GetSubRegion(esr);

      array3d<R1Tensor> const &
      dNdX = cellBlock->getReference< array3d<R1Tensor> >(keys::dNdX);

      arrayView2d<real64> const & detJ = cellBlock->getReference< array2d<real64> >(keys::detJ);

      arrayView2d< localIndex > const & elemsToNodes = cellBlock->nodeList();
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


      GEOSX_MARK_LOOP_BEGIN(elemLoop,elemLoop);

      for( localIndex k=0 ; k<cellBlock->size() ; ++k )
      {

        real64 stiffness[6][6];
        constitutiveRelations[er][esr][0]->GetStiffness( stiffness );

        if(elemGhostRank[k] < 0)
        {
          for( localIndex a=0 ; a<numNodesPerElement ; ++a)
          {

            localIndex localNodeIndex = elemsToNodes[k][a];

            for( int i=0 ; i<dim ; ++i )
            {
              elementLocalDofIndex[static_cast<int>(a)*dim+i] = dim*trilinos_index[localNodeIndex]+i;

              // TODO must add last solution estimate for this to be valid
              element_dof_np1(static_cast<int>(a)*dim+i) = disp[localNodeIndex][i];
            }
          }

          if( this->m_timeIntegrationOption == timeIntegrationOption::ImplicitDynamic )
          {
            GEOS_ERROR("Option not supported");
            CopyGlobalToLocal<R1Tensor>( elemsToNodes[k],
                               disp, uhat, vtilde, uhattilde, u_local, uhat_local, vtilde_local, uhattilde_local,
                               numNodesPerElement );
          }
          else
          {
            CopyGlobalToLocal<R1Tensor>( elemsToNodes[k], disp, uhat, u_local, uhat_local, numNodesPerElement );
          }

          R2SymTensor referenceStress;
          if( !fluidPres[er][esr].empty() )
          {
            referenceStress.PlusIdentity( - biotCoefficient[er][esr][0] * (fluidPres[er][esr][k] + dPres[er][esr][k]));
          }
          real64 maxElemForce = CalculateElementResidualAndDerivative( density[er][esr][0],
                                                                       feDiscretization->m_finiteElement,
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


          if( maxElemForce > m_maxForce )
          {
            m_maxForce = maxElemForce;
          }

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
                                               viewKeyStruct::trilinosIndexString,
                                               3,
                                               blockSystem,
                                               BlockIDs::displacementBlock );
  });

  ApplyTractionBC( domain,
                   time_n+dt,
                   *blockSystem );

  ApplyDisplacementBC_implicit( time_n + dt, *domain, *blockSystem );
//  bcManager->ApplyBoundaryCondition( this, &,
//                                     nodeManager, keys::TotalDisplacement, time_n + dt, *blockSystem );

  Epetra_FECrsMatrix * const matrix = blockSystem->GetMatrix( BlockIDs::displacementBlock,
                                                                    BlockIDs::displacementBlock );
  Epetra_FEVector * const rhs = blockSystem->GetResidualVector( BlockIDs::displacementBlock );

  if( verboseLevel() >= 2 )
  {
    matrix->Print(std::cout);
    rhs->Print(std::cout);
  }

  matrix->GlobalAssemble(true);
  rhs->GlobalAssemble();


}

real64
SolidMechanics_LagrangianFEM::
CalculateResidualNorm(systemSolverInterface::EpetraBlockSystem const *const blockSystem, DomainPartition *const domain)
{

  Epetra_FEVector const * const
  residual = blockSystem->GetResidualVector( BlockIDs::displacementBlock );

  real64 localResidual[2] = {0.0, this->m_maxForce};
//  residual->Norm2(&scalarResidual);

  real64 * residualData = nullptr;
  int length;
  residual->ExtractView(&residualData,&length);
  for( localIndex i=0 ; i<length ; ++i )
  {
    localResidual[0] += residualData[i]*residualData[i];
  }


  real64 globalResidualNorm[2] = {0,0};
//  MPI_Allreduce (&localResidual,&globalResidualNorm,1,MPI_DOUBLE,MPI_SUM ,MPI_COMM_GEOSX);


  int rank, size;
  MPI_Comm_rank(MPI_COMM_GEOSX, &rank);
  MPI_Comm_size(MPI_COMM_GEOSX, &size);
  array1d<real64> globalValues( size * 2 );
  globalValues = 0;
  MPI_Gather( localResidual,
              2,
              MPI_DOUBLE,
              globalValues.data(),
              2,
              MPI_DOUBLE,
              0,
              MPI_COMM_GEOSX );

  if( rank==0 )
  {
    for( int r=0 ; r<size ; ++r )
    {
      globalResidualNorm[0] += globalValues[r*2];

      if( globalResidualNorm[1] < globalValues[r*2+1] )
      {
        globalResidualNorm[1] = globalValues[r*2+1];
      }
    }
  }

  MPI_Bcast( globalResidualNorm, 2, MPI_DOUBLE, 0, MPI_COMM_GEOSX );



  return sqrt(globalResidualNorm[0])/(globalResidualNorm[1]+1);

}

realT SolidMechanics_LagrangianFEM::CalculateElementResidualAndDerivative( real64 const density,
                                                                           FiniteElementBase const * const fe,
                                                                           arraySlice2d<R1Tensor const> const& dNdX,
                                                                           arraySlice1d<realT const> const& detJ,
                                                                           R2SymTensor const * const refStress,
                                                                           r1_array const& u,
                                                                           r1_array const& uhat,
                                                                           r1_array const& uhattilde,
                                                                           r1_array const& vtilde,
                                                                           realT const dt,
                                                                           Epetra_SerialDenseMatrix& dRdU,
                                                                           Epetra_SerialDenseVector& R,
                                                                           real64 c[6][6] )
{
  const integer dim = 3;
  realT maxForce = 0;
  realT amass = getReference<real64>(solidMechanicsViewKeys.massDamping);
  realT astiff = getReference<real64>(solidMechanicsViewKeys.stiffnessDamping);
  real64 const newmarkBeta = getReference<real64>(solidMechanicsViewKeys.newmarkBeta);
  real64 const newmarkGamma = getReference<real64>(solidMechanicsViewKeys.newmarkGamma);


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
      dNdXa = dNdX[q][a];

      for( integer b=0 ; b<fe->dofs_per_element() ; ++b )
      {
//        realT const * const dNdXb = dNdX(q,b).Data();
        dNdXb = dNdX[q][b];

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
        dNdXa = dNdX[q][a];

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


  globalIndex_array const& trilinos_index = nodeManager->getReference< globalIndex_array >(solidMechanicsViewKeys.trilinosIndex);

  r1_array& X        = nodeManager->getReference<r1_array>(nodeManager->viewKeys.referencePosition);
  r1_array& disp     = nodeManager->getReference<r1_array>(keys::TotalDisplacement);
  r1_array& incdisp  = nodeManager->getReference<r1_array>(keys::IncrementalDisplacement);

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

  std::map<string, string_array > fieldNames;
  fieldNames["node"].push_back(keys::IncrementalDisplacement);
  fieldNames["node"].push_back(keys::TotalDisplacement);

  CommunicationTools::SynchronizeFields( fieldNames,
                                         domain->getMeshBody(0)->getMeshLevel(0),
                                         domain->getReference< array1d<NeighborCommunicator> >( domain->viewKeys.neighbors ) );



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

  r1_array& incdisp  = nodeManager->getReference<r1_array>(keys::IncrementalDisplacement);

  // TODO need to finish this rewind
  forall_in_range(0, nodeManager->size(), GEOSX_LAMBDA (localIndex a) mutable
  {
    incdisp[a] = 0.0;
  });
}

REGISTER_CATALOG_ENTRY( SolverBase, SolidMechanics_LagrangianFEM, std::string const &, ManagedGroup * const )
} /* namespace ANST */
