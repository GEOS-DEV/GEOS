/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file SolidMechanicsLagrangianFEM.cpp
 */

#include "SolidMechanicsLagrangianFEM.hpp"

#include <vector>
#include <math.h>
#include <sys/time.h>

#include "common/TimingMacros.hpp"
#include "managers/FieldSpecification/FieldSpecificationManager.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "managers/NumericalMethodsManager.hpp"
#include "finiteElement/ElementLibrary/FiniteElement.h"
#include "finiteElement/FiniteElementDiscretizationManager.hpp"
#include "finiteElement/Kinematics.h"
#include "codingUtilities/Utilities.hpp"
#include "mesh/FaceElementSubRegion.hpp"

#include "managers/DomainPartition.hpp"
#include "meshUtilities/ComputationalGeometry.hpp"
#include "MPI_Communications/CommunicationTools.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"
#include "MPI_Communications/NeighborCommunicator.hpp"


//#define verbose 0 //Need to move this somewhere else

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;

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

SolidMechanicsLagrangianFEM::SolidMechanicsLagrangianFEM( const std::string& name,
                                                          Group * const parent ):
  SolverBase( name, parent ),
  m_newmarkGamma(0.5),
  m_newmarkBeta(0.25),
  m_massDamping(0.0),
  m_stiffnessDamping(0.0),
  m_timeIntegrationOptionString(),
  m_timeIntegrationOption(timeIntegrationOption::ExplicitDynamic),
  m_useVelocityEstimateForQS(0),
  m_maxForce(0.0),
  m_maxNumResolves(10),
  m_strainTheory(0),
  m_solidMaterialName(""),
  m_solidMaterialFullIndex(0),
  m_elemsAttachedToSendOrReceiveNodes(),
  m_elemsNotAttachedToSendOrReceiveNodes(),
  m_sendOrReceiveNodes(),
  m_nonSendOrReceiveNodes(),
  m_iComm()
{
  m_sendOrReceiveNodes.setUserCallBack("SolidMechanicsLagrangianFEM::m_sendOrReceiveNodes");
  m_nonSendOrReceiveNodes.setUserCallBack("SolidMechanicsLagrangianFEM::m_nonSendOrReceiveNodes");

  registerWrapper(viewKeyStruct::newmarkGammaString, &m_newmarkGamma, false )->
    setApplyDefaultValue(0.5)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Value of :math:`\\gamma` in the Newmark Method for Implicit Dynamic time integration option");

  registerWrapper(viewKeyStruct::newmarkBetaString, &m_newmarkBeta, false )->
    setApplyDefaultValue(0.25)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Value of :math:`\\beta` in the Newmark Method for Implicit Dynamic time integration option. "
                   "This should be pow(newmarkGamma+0.5,2.0)/4.0 unless you know what you are doing.");

  registerWrapper(viewKeyStruct::massDampingString, &m_massDamping, false )->
    setApplyDefaultValue(0.0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Value of mass based damping coefficient. ");

  registerWrapper(viewKeyStruct::stiffnessDampingString, &m_stiffnessDamping, false )->
    setApplyDefaultValue(0.0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Value of stiffness based damping coefficient. ");

  registerWrapper(viewKeyStruct::timeIntegrationOptionString, &m_timeIntegrationOption, false )->
    setInputFlag(InputFlags::FALSE)->
    setDescription("Time integration enum class value.");

  registerWrapper(viewKeyStruct::timeIntegrationOptionStringString, &m_timeIntegrationOptionString, false )->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Time integration method. Options are: \n QuasiStatic \n ImplicitDynamic \n ExplicitDynamic");

  registerWrapper(viewKeyStruct::useVelocityEstimateForQSString, &m_useVelocityEstimateForQS, false )->
    setApplyDefaultValue(0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription( "Flag to indicate the use of the incremental displacement from the previous step as an "
                    "initial estimate for the incremental displacement of the current step.");

  registerWrapper(viewKeyStruct::maxNumResolvesString, &m_maxNumResolves, false )->
    setApplyDefaultValue(10)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription( "Value to indicate how many resolves may be executed after some other event is executed. "
                    "For example, if a SurfaceGenerator is specified, it will be executed after the mechanics solve. "
                    "However if a new surface is generated, then the mechanics solve must be executed again due to the "
                    "change in topology.");

  registerWrapper(viewKeyStruct::strainTheoryString, &m_strainTheory, false )->
    setApplyDefaultValue(0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription( "Indicates whether or not to use "
                    "`Infinitesimal Strain Theory <https://en.wikipedia.org/wiki/Infinitesimal_strain_theory>`_, or "
                    "`Finite Strain Theory <https://en.wikipedia.org/wiki/Finite_strain_theory>`_. Valid Inputs are:\n"
                    " 0 - Infinitesimal Strain \n"
                    " 1 - Finite Strain");

  registerWrapper(viewKeyStruct::solidMaterialNameString, &m_solidMaterialName, false )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription( "The name of the material that should be used in the constitutive updates");


}

void SolidMechanicsLagrangianFEM::PostProcessInput()
{
  SolverBase::PostProcessInput();

  m_linearSolverParameters.amg.isSymmetric = true;
  m_linearSolverParameters.dofsPerNode = 3;

  if( !m_timeIntegrationOptionString.empty() )
  {
    SetTimeIntegrationOption( m_timeIntegrationOptionString );
  }
}

SolidMechanicsLagrangianFEM::~SolidMechanicsLagrangianFEM()
{
  // TODO Auto-generated destructor stub
}


void SolidMechanicsLagrangianFEM::RegisterDataOnMesh( Group * const MeshBodies )
{
  for( auto & mesh : MeshBodies->GetSubGroups() )
  {
    NodeManager * const nodes = mesh.second->group_cast<MeshBody*>()->getMeshLevel(0)->getNodeManager();


    nodes->registerWrapper<array1d<R1Tensor> >( keys::TotalDisplacement )->
      setPlotLevel(PlotLevel::LEVEL_0)->
      setRegisteringObjects(this->getName())->
      setDescription( "An array that holds the total displacements on the nodes.");

    nodes->registerWrapper<array1d<R1Tensor> >( keys::IncrementalDisplacement )->
      setPlotLevel(PlotLevel::LEVEL_3)->
      setRegisteringObjects(this->getName())->
      setDescription( "An array that holds the incremental displacements for the current time step on the nodes.");

    nodes->registerWrapper<array1d<R1Tensor> >( keys::Velocity )->
      setPlotLevel(PlotLevel::LEVEL_0)->
      setRegisteringObjects(this->getName())->
      setDescription( "An array that holds the current velocity on the nodes.");

    nodes->registerWrapper<array1d<R1Tensor> >( keys::Acceleration )->
      setPlotLevel(PlotLevel::LEVEL_1)->
      setRegisteringObjects(this->getName())->
      setDescription( "An array that holds the current acceleration on the nodes. This array also is used "
                      "to hold the summation of nodal forces resulting from the governing equations.");

    nodes->registerWrapper<array1d<R1Tensor> >( viewKeyStruct::forceExternal )->
      setPlotLevel(PlotLevel::LEVEL_0)->
      setRegisteringObjects(this->getName())->
      setDescription( "An array that holds the external forces on the nodes. This includes any boundary"
                      " conditions as well as coupling forces such as hydraulic forces.");

    nodes->registerWrapper<array1d<real64> >( keys::Mass )->setPlotLevel(PlotLevel::LEVEL_0)->
        setPlotLevel(PlotLevel::LEVEL_0)->
        setRegisteringObjects(this->getName())->
        setDescription( "An array that holds the mass on the nodes.");

    nodes->registerWrapper<array1d<R1Tensor> >( viewKeyStruct::vTildeString )->
      setPlotLevel(PlotLevel::NOPLOT)->
      setRegisteringObjects(this->getName())->
      setDescription( "An array that holds the velocity predictors on the nodes.");

    nodes->registerWrapper<array1d<R1Tensor> >( viewKeyStruct::uhatTildeString )->
      setPlotLevel(PlotLevel::NOPLOT)->
      setRegisteringObjects(this->getName())->
      setDescription( "An array that holds the incremental displacement predictors on the nodes.");

  }
}


void SolidMechanicsLagrangianFEM::InitializePreSubGroups(Group * const rootGroup)
{
  SolverBase::InitializePreSubGroups(rootGroup);

  DomainPartition * domain = rootGroup->GetGroup<DomainPartition>(keys::domain);
  ConstitutiveManager const * const cm = domain->getConstitutiveManager();

  ConstitutiveBase const * const solid  = cm->GetConstitutiveRelation<ConstitutiveBase>( m_solidMaterialName );
  GEOS_ERROR_IF( solid == nullptr, "constitutive model " + m_solidMaterialName + " not found" );
  m_solidMaterialFullIndex = solid->getIndexInParent();

}

void SolidMechanicsLagrangianFEM::updateIntrinsicNodalData( DomainPartition * const domain )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  NodeManager * const nodes = mesh->getNodeManager();
  //FaceManager * const faceManager = mesh->getFaceManager();

  ElementRegionManager const * const elementRegionManager = mesh->getElemManager();
  ConstitutiveManager const * const constitutiveManager = domain->GetGroup<ConstitutiveManager >(keys::ConstitutiveManager);

  arrayView1d<real64> & mass = nodes->getReference<array1d<real64>>(keys::Mass);
  mass = 0.0;

  NumericalMethodsManager const *
  numericalMethodManager = domain->getParent()->GetGroup<NumericalMethodsManager>(keys::numericalMethodsManager);

  FiniteElementDiscretizationManager const *
  feDiscretizationManager = numericalMethodManager->GetGroup<FiniteElementDiscretizationManager>(keys::finiteElementDiscretizations);

  FiniteElementDiscretization const *
  feDiscretization = feDiscretizationManager->GetGroup<FiniteElementDiscretization>(m_discretizationName);

  ElementRegionManager::MaterialViewAccessor< arrayView2d<real64> > rho =
    elementRegionManager->ConstructFullMaterialViewAccessor< array2d<real64>, arrayView2d<real64> >("density", constitutiveManager);

  for( localIndex er=0 ; er<elementRegionManager->numRegions() ; ++er )
  {
    ElementRegionBase const * const elemRegion = elementRegionManager->GetRegion(er);

    elemRegion->forElementSubRegionsIndex<CellElementSubRegion>([&]( localIndex const esr, CellElementSubRegion const * const elementSubRegion )
    {
      arrayView2d<real64> const & detJ = elementSubRegion->getReference< array2d<real64> >(keys::detJ);
      arrayView2d<localIndex> const & elemsToNodes = elementSubRegion->nodeList();

      std::unique_ptr<FiniteElementBase>
      fe = feDiscretization->getFiniteElement( elementSubRegion->GetElementTypeString() );

      for( localIndex k=0 ; k < elemsToNodes.size(0) ; ++k )
      {

        // TODO this integration needs to be be carried out properly.
        real64 elemMass = 0;
        for( localIndex q=0 ; q<fe->n_quadrature_points() ; ++q )
        {
          elemMass += rho[er][esr][m_solidMaterialFullIndex][k][q] * detJ[k][q];
        }
        for( localIndex a=0 ; a< elemsToNodes.size(1) ; ++a )
        {
          mass[elemsToNodes[k][a]] += elemMass/elemsToNodes.size(1);
        }

        for( localIndex a=0 ; a<elementSubRegion->numNodesPerElement() ; ++a )
        {
          if( nodes->GhostRank()[elemsToNodes[k][a]] >= -1 )
          {
            m_sendOrReceiveNodes.insert( elemsToNodes[k][a] );
          }
          else
          {
            m_nonSendOrReceiveNodes.insert( elemsToNodes[k][a] );
          }
        }
      }
    });
  }
}

void SolidMechanicsLagrangianFEM::InitializePostInitialConditions_PreSubGroups( Group * const problemManager )
{
  DomainPartition * domain = problemManager->GetGroup<DomainPartition>(keys::domain);
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);

  NodeManager * const nodes = mesh->getNodeManager();
  //FaceManager * const faceManager = mesh->getFaceManager();


  ElementRegionManager * elementRegionManager = mesh->getElemManager();
  ConstitutiveManager * constitutiveManager = domain->GetGroup<ConstitutiveManager >(keys::ConstitutiveManager);

  arrayView1d<real64> & mass = nodes->getReference<array1d<real64>>(keys::Mass);

  ElementRegionManager::MaterialViewAccessor< arrayView2d<real64> > rho =
    elementRegionManager->ConstructFullMaterialViewAccessor< array2d<real64>, arrayView2d<real64> >("density", constitutiveManager);


  NumericalMethodsManager const *
  numericalMethodManager = domain->getParent()->GetGroup<NumericalMethodsManager>(keys::numericalMethodsManager);

  FiniteElementDiscretizationManager const *
  feDiscretizationManager = numericalMethodManager->GetGroup<FiniteElementDiscretizationManager>(keys::finiteElementDiscretizations);

  FiniteElementDiscretization const *
  feDiscretization = feDiscretizationManager->GetGroup<FiniteElementDiscretization>(m_discretizationName);

  m_elemsAttachedToSendOrReceiveNodes.resize( elementRegionManager->numRegions() );
  m_elemsNotAttachedToSendOrReceiveNodes.resize( elementRegionManager->numRegions() );

  for( localIndex er=0 ; er<elementRegionManager->numRegions() ; ++er )
  {
    ElementRegionBase const * const elemRegion = elementRegionManager->GetRegion(er);
    m_elemsAttachedToSendOrReceiveNodes[er].resize( elemRegion->numSubRegions() );
    m_elemsNotAttachedToSendOrReceiveNodes[er].resize( elemRegion->numSubRegions() );

    elemRegion->forElementSubRegionsIndex<CellElementSubRegion>([&]( localIndex const esr, CellElementSubRegion const * const elementSubRegion )
    {
      m_elemsAttachedToSendOrReceiveNodes[er][esr].setUserCallBack(
        "SolidMechanicsLagrangianFEM::m_elemsAttachedToSendOrReceiveNodes["
        + std::to_string(er) + "][" + std::to_string(esr) + "]" );

      m_elemsNotAttachedToSendOrReceiveNodes[er][esr].setUserCallBack(
        "SolidMechanicsLagrangianFEM::m_elemsNotAttachedToSendOrReceiveNodes["
        + std::to_string(er) + "][" + std::to_string(esr) + "]" );

      arrayView2d<real64> const & detJ = elementSubRegion->getReference< array2d<real64> >(keys::detJ);
      arrayView2d<localIndex> const & elemsToNodes = elementSubRegion->nodeList();

      std::unique_ptr<FiniteElementBase>
      fe = feDiscretization->getFiniteElement( elementSubRegion->GetElementTypeString() );

      for( localIndex k=0 ; k < elemsToNodes.size(0) ; ++k )
      {

        // TODO this integration needs to be be carried out properly.
        real64 elemMass = 0;
        for( localIndex q=0 ; q<fe->n_quadrature_points() ; ++q )
        {
          elemMass += rho[er][esr][m_solidMaterialFullIndex][k][q] * detJ[k][q];
        }
        for( localIndex a=0 ; a< elemsToNodes.size(1) ; ++a )
        {
          mass[elemsToNodes[k][a]] += elemMass/elemsToNodes.size(1);
        }

        bool isAttachedToGhostNode = false;
        for( localIndex a=0 ; a<elementSubRegion->numNodesPerElement() ; ++a )
        {
          if( nodes->GhostRank()[elemsToNodes[k][a]] >= -1 )
          {
            isAttachedToGhostNode = true;
            m_sendOrReceiveNodes.insert( elemsToNodes[k][a] );
          }
          else
          {
            m_nonSendOrReceiveNodes.insert( elemsToNodes[k][a] );
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
    });
  }
}

real64 SolidMechanicsLagrangianFEM::SolverStep( real64 const& time_n,
                                                real64 const& dt,
                                                const int cycleNumber,
                                                DomainPartition * domain )
{
  real64 dtReturn = dt;

  SolverBase * const surfaceGenerator =  this->getParent()->GetGroup<SolverBase>("SurfaceGen");

  if( m_timeIntegrationOption == timeIntegrationOption::ExplicitDynamic )
  {
    dtReturn = ExplicitStep( time_n, dt, cycleNumber, Group::group_cast<DomainPartition*>(domain) );

    if( surfaceGenerator!=nullptr )
    {
      surfaceGenerator->SolverStep( time_n, dt, cycleNumber, domain );
    }

  }
  else if( m_timeIntegrationOption == timeIntegrationOption::ImplicitDynamic ||
           m_timeIntegrationOption == timeIntegrationOption::QuasiStatic )
  {
    int const maxNumResolves = m_maxNumResolves;
    int locallyFractured = 0;
    int globallyFractured = 0;
    for( int solveIter=0 ; solveIter<maxNumResolves ; ++solveIter )
    {
      ImplicitStepSetup( time_n, dt, domain, m_dofManager, m_matrix, m_rhs, m_solution );

      dtReturn = NonlinearImplicitStep( time_n, dt, cycleNumber, domain->group_cast<DomainPartition *>(), m_dofManager,
                                        m_matrix, m_rhs, m_solution );
      if( surfaceGenerator!=nullptr )
      {
        if( surfaceGenerator->SolverStep( time_n, dt, cycleNumber, domain ) > 0 )
        {
          locallyFractured = 1;
        }
        MPI_Allreduce( &locallyFractured,
                       &globallyFractured,
                       1,
                       MPI_INT,
                       MPI_MAX,
                       MPI_COMM_GEOSX);
      }
      if( globallyFractured == 0 )
      {
        break;
      }
    }
    ImplicitStepComplete( time_n, dt,  domain );
  }

  return dtReturn;
}

real64 SolidMechanicsLagrangianFEM::ExplicitStep( real64 const& time_n,
                                                  real64 const& dt,
                                                  const int GEOSX_UNUSED_ARG( cycleNumber ),
                                                  DomainPartition * const domain )
{
  GEOSX_MARK_FUNCTION;

  // updateIntrinsicNodalData(domain);

  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  NodeManager * const nodes = mesh->getNodeManager();
  ElementRegionManager * const elemManager = mesh->getElemManager();
  NumericalMethodsManager const * const numericalMethodManager = domain->getParent()->GetGroup<NumericalMethodsManager>(keys::numericalMethodsManager);
  FiniteElementDiscretizationManager const * const feDiscretizationManager = numericalMethodManager->GetGroup<FiniteElementDiscretizationManager>(keys::finiteElementDiscretizations);
  ConstitutiveManager * const constitutiveManager = domain->GetGroup<ConstitutiveManager >(keys::ConstitutiveManager);

  FieldSpecificationManager * const fsManager = FieldSpecificationManager::get();

  arrayView1d<real64 const> const & mass = nodes->getReference<array1d<real64>>(keys::Mass);
  array1d<R1Tensor> & velocityArray = nodes->getReference<array1d<R1Tensor>>(keys::Velocity);
  arrayView1d<R1Tensor> const & vel = velocityArray;

  arrayView1d<R1Tensor> const & u = nodes->getReference<array1d<R1Tensor>>(keys::TotalDisplacement);
  arrayView1d<R1Tensor> const & uhat = nodes->getReference<array1d<R1Tensor>>(keys::IncrementalDisplacement);
  arrayView1d<R1Tensor> const & acc = nodes->getReference<array1d<R1Tensor>>(keys::Acceleration);

  array1d<NeighborCommunicator> & neighbors = domain->getReference< array1d<NeighborCommunicator> >( domain->viewKeys.neighbors );
  std::map<string, string_array > fieldNames;
  fieldNames["node"].push_back("Velocity");

  CommunicationTools::SynchronizePackSendRecvSizes( fieldNames, mesh, neighbors, m_iComm );

  fsManager->ApplyFieldValue< parallelDevicePolicy< 1024 > >( time_n, domain, "nodeManager", keys::Acceleration );

  //3: v^{n+1/2} = v^{n} + a^{n} dt/2
  SolidMechanicsLagrangianFEMKernels::velocityUpdate( acc, vel, dt/2 );

  fsManager->ApplyFieldValue< parallelDevicePolicy< 1024 > >( time_n, domain, "nodeManager", keys::Velocity );

  //4. x^{n+1} = x^{n} + v^{n+{1}/{2}} dt (x is displacement)
  SolidMechanicsLagrangianFEMKernels::displacementUpdate( vel, uhat, u, dt );

  fsManager->ApplyFieldValue( time_n + dt, domain, "nodeManager", keys::TotalDisplacement,
    [&]( FieldSpecificationBase const * const bc, SortedArrayView<localIndex const> const & targetSet )->void
    {
      integer const component = bc->GetComponent();
      forall_in_range< parallelDevicePolicy< 1024 > >(0, targetSet.size(),
        GEOSX_DEVICE_LAMBDA( localIndex const i )
        {
          localIndex const a = targetSet[ i ];
          vel[a][component] = u[a][component];
        }
      );
    },
    [&]( FieldSpecificationBase const * const bc, SortedArrayView<localIndex const> const & targetSet )->void
    {
      integer const component = bc->GetComponent();
      forall_in_range< parallelDevicePolicy< 1024 > >(0, targetSet.size(),
        GEOSX_DEVICE_LAMBDA( localIndex const i )
        {
          localIndex const a = targetSet[ i ];
          uhat[a][component] = u[a][component] - vel[a][component];
          vel[a][component]  = uhat[a][component] / dt;
        }
      );
    }
  );

  ElementRegionManager::MaterialViewAccessor< arrayView2d<real64> >
  meanStress = elemManager->ConstructFullMaterialViewAccessor< array2d<real64>,
                                                               arrayView2d<real64> >("MeanStress",
                                                                                     constitutiveManager);

  ElementRegionManager::MaterialViewAccessor< arrayView2d<R2SymTensor> > const
  devStress = elemManager->ConstructFullMaterialViewAccessor< array2d<R2SymTensor>,
                                                              arrayView2d<R2SymTensor> >("DeviatorStress",
                                                                                         constitutiveManager);

  ElementRegionManager::ConstitutiveRelationAccessor<ConstitutiveBase> constitutiveRelations =
    elemManager->ConstructFullConstitutiveAccessor<ConstitutiveBase>(constitutiveManager);

  //Step 5. Calculate deformation input to constitutive model and update state to
  // Q^{n+1}
  for( localIndex er=0 ; er<elemManager->numRegions() ; ++er )
  {
    ElementRegionBase * const elementRegion = elemManager->GetRegion(er);
    FiniteElementDiscretization const * feDiscretization = feDiscretizationManager->GetGroup<FiniteElementDiscretization>(m_discretizationName);

    elementRegion->forElementSubRegionsIndex<CellElementSubRegion>([&]( localIndex const esr, CellElementSubRegion const * const elementSubRegion )
    {
      arrayView3d< R1Tensor > const & dNdX = elementSubRegion->getReference< array3d< R1Tensor > >(keys::dNdX);

      arrayView2d<real64> const & detJ = elementSubRegion->getReference< array2d<real64> >(keys::detJ);

      arrayView2d<localIndex> const & elemsToNodes = elementSubRegion->nodeList();

      localIndex const numNodesPerElement = elemsToNodes.size(1);

      localIndex const numQuadraturePoints = feDiscretization->m_finiteElement->n_quadrature_points();

      ExplicitElementKernelLaunch( numNodesPerElement,
                                   numQuadraturePoints,
                                   constitutiveRelations[er][esr][m_solidMaterialFullIndex],
                                   this->m_elemsAttachedToSendOrReceiveNodes[er][esr],
                                   elemsToNodes,
                                   dNdX,
                                   detJ,
                                   u,
                                   vel,
                                   acc,
                                   meanStress[er][esr][m_solidMaterialFullIndex],
                                   devStress[er][esr][m_solidMaterialFullIndex],
                                   dt );

    }); //Element Region

  } //Element Manager

  // apply this over a set
  SolidMechanicsLagrangianFEMKernels::velocityUpdate( acc, mass, vel, dt / 2, m_sendOrReceiveNodes );

  fsManager->ApplyFieldValue< parallelDevicePolicy< 1024 > >( time_n, domain, "nodeManager", keys::Velocity );

  // HACK: Move velocity back to the CPU to be packed. It is not modified so we don't touch it.
  if (neighbors.size() > 0) velocityArray.move(chai::CPU, false);
  CommunicationTools::SynchronizePackSendRecv( fieldNames, mesh, neighbors, m_iComm );

  for( localIndex er=0 ; er<elemManager->numRegions() ; ++er )
  {
    ElementRegionBase * const elementRegion = elemManager->GetRegion(er);

    FiniteElementDiscretization const * feDiscretization = feDiscretizationManager->GetGroup<FiniteElementDiscretization>(m_discretizationName);

    elementRegion->forElementSubRegionsIndex<CellElementSubRegion>([&]( localIndex const esr, CellElementSubRegion const * const elementSubRegion )
    {
      arrayView3d< R1Tensor > const & dNdX = elementSubRegion->getReference< array3d< R1Tensor > >(keys::dNdX);

      arrayView2d<real64> const & detJ = elementSubRegion->getReference< array2d<real64> >(keys::detJ);

      arrayView2d<localIndex> const & elemsToNodes = elementSubRegion->nodeList();

      localIndex const numNodesPerElement = elemsToNodes.size(1);

      localIndex const numQuadraturePoints = feDiscretization->m_finiteElement->n_quadrature_points();

      ExplicitElementKernelLaunch( numNodesPerElement,
                                   numQuadraturePoints,
                                   constitutiveRelations[er][esr][m_solidMaterialFullIndex],
                                   this->m_elemsNotAttachedToSendOrReceiveNodes[er][esr],
                                   elemsToNodes,
                                   dNdX,
                                   detJ,
                                   u,
                                   vel,
                                   acc,
                                   meanStress[er][esr][m_solidMaterialFullIndex],
                                   devStress[er][esr][m_solidMaterialFullIndex],
                                   dt );
    }); //Element Region

  } //Element Manager

  // apply this over a set
  SolidMechanicsLagrangianFEMKernels::velocityUpdate( acc, mass, vel, dt / 2, m_nonSendOrReceiveNodes );

  fsManager->ApplyFieldValue< parallelDevicePolicy< 1024 > >( time_n, domain, "nodeManager", keys::Velocity );

  // HACK: Move velocity back to the CPU to be unpacked. It is modified so we touch it.
  if (neighbors.size() > 0) velocityArray.move(chai::CPU, true);
  CommunicationTools::SynchronizeUnpack( mesh, neighbors, m_iComm );

  return dt;
}



void SolidMechanicsLagrangianFEM::ApplyDisplacementBC_implicit( real64 const time,
                                                                DofManager const & dofManager,
                                                                DomainPartition & domain,
                                                                ParallelMatrix & matrix,
                                                                ParallelVector & rhs )
{
  string const dofKey = dofManager.getKey( keys::TotalDisplacement );

  FieldSpecificationManager const * const fsManager = FieldSpecificationManager::get();

  fsManager->Apply( time,
                     &domain,
                     "nodeManager",
                     keys::TotalDisplacement,
                     [&]( FieldSpecificationBase const * const bc,
                          string const &,
                          set<localIndex> const & targetSet,
                          Group * const targetGroup,
                          string const fieldName )
  {
    bc->ApplyBoundaryConditionToSystem<FieldSpecificationEqual, LAInterface>( targetSet,
                                                                              false,
                                                                              time,
                                                                              targetGroup,
                                                                              fieldName,
                                                                              dofKey,
                                                                              3,
                                                                              matrix,
                                                                              rhs );
  } );
}


void SolidMechanicsLagrangianFEM::ApplyTractionBC( real64 const time,
                                                   DofManager const & dofManager,
                                                   DomainPartition * const domain,
                                                   ParallelVector & rhs )
{
  FieldSpecificationManager * const fsManager = FieldSpecificationManager::get();
  NewFunctionManager * const functionManager = NewFunctionManager::Instance();

  FaceManager * const faceManager = domain->getMeshBody(0)->getMeshLevel(0)->getFaceManager();
  NodeManager * const nodeManager = domain->getMeshBody(0)->getMeshLevel(0)->getNodeManager();

  real64_array const & faceArea  = faceManager->getReference<real64_array>("faceArea");
  ArrayOfArraysView< localIndex const > const & faceToNodeMap = faceManager->nodeList();

  string const dofKey = dofManager.getKey( keys::TotalDisplacement );

  arrayView1d<globalIndex> const &
  blockLocalDofNumber = nodeManager->getReference<globalIndex_array>( dofKey );

  arrayView1d<integer const> const & faceGhostRank = faceManager->GhostRank();
  fsManager->Apply( time,
                    domain,
                    "faceManager",
                    string("Traction"),
                    [&]( FieldSpecificationBase const * const bc,
                    string const &,
                    set<localIndex> const & targetSet,
                    Group * const GEOSX_UNUSED_ARG( targetGroup ),
                    string const GEOSX_UNUSED_ARG( fieldName ) ) -> void
  {
    string const & functionName = bc->getReference<string>( FieldSpecificationBase::viewKeyStruct::functionNameString);

    globalIndex_array nodeDOF;
    real64_array nodeRHS;
    integer const component = bc->GetComponent();

    if( functionName.empty() )
    {
      for( auto kf : targetSet )
      {
        if( faceGhostRank[kf] < 0 )
        {
          localIndex const numNodes = faceToNodeMap.sizeOfArray( kf );
          nodeDOF.resize( numNodes );
          nodeRHS.resize( numNodes );
          for( localIndex a=0 ; a<numNodes ; ++a )
          {
            nodeDOF[a] = blockLocalDofNumber[ faceToNodeMap( kf, a ) ] + component;
            nodeRHS[a] = bc->GetScale() * faceArea[kf] / numNodes;
          }
          rhs.add( nodeDOF, nodeRHS );
        }
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
            if( faceGhostRank[kf] < 0 )
            {
              localIndex const numNodes = faceToNodeMap.sizeOfArray( kf );
              nodeDOF.resize( numNodes );
              nodeRHS.resize( numNodes );
              for( localIndex a=0 ; a<numNodes ; ++a )
              {
                nodeDOF[a] = blockLocalDofNumber[ faceToNodeMap( kf, a ) ] + component;
                nodeRHS[a] = value * faceArea[kf] / numNodes;
              }
              rhs.add( nodeDOF, nodeRHS );
            }
          }
        }
        else
        {
          real64_array result;
          result.resize( targetSet.size() );
          function->Evaluate( faceManager, time, targetSet, result );

          for( auto kf : targetSet )
          {
            if( faceGhostRank[kf] < 0 )
            {
              localIndex const numNodes = faceToNodeMap.sizeOfArray( kf );
              nodeDOF.resize( numNodes );
              nodeRHS.resize( numNodes );
              for( localIndex a=0 ; a<numNodes ; ++a )
              {
                nodeDOF[a] = blockLocalDofNumber[ faceToNodeMap( kf, a ) ] + component;
                nodeRHS[a] = result[kf] * faceArea[kf] / numNodes;
              }
              rhs.add( nodeDOF, nodeRHS );
            }
          }
      }
    }
  });
}

void SolidMechanicsLagrangianFEM::ApplyChomboPressure( DofManager const & dofManager,
                                                       DomainPartition * const domain,
                                                       ParallelVector & rhs )
{
  FaceManager * const faceManager = domain->getMeshBody(0)->getMeshLevel(0)->getFaceManager();
  NodeManager * const nodeManager = domain->getMeshBody(0)->getMeshLevel(0)->getNodeManager();

  arrayView1d<real64 const> const & faceArea  = faceManager->faceArea();
  arrayView1d<R1Tensor const> const & faceNormal  = faceManager->faceNormal();
  ArrayOfArraysView< localIndex const > const & faceToNodeMap = faceManager->nodeList();

  string const dofKey = dofManager.getKey( keys::TotalDisplacement );

  arrayView1d<globalIndex> const &
  blockLocalDofNumber =  nodeManager->getReference<globalIndex_array>( dofKey );

  arrayView1d<real64 const> const & facePressure = faceManager->getReference< array1d<real64> >("ChomboPressure");

  for( localIndex kf=0 ; kf<faceManager->size() ; ++kf )
  {
    globalIndex nodeDOF[20];
    real64 nodeRHS[20];

    int const numNodes = integer_conversion<int>(faceToNodeMap.sizeOfArray(kf));
    for( int a=0 ; a<numNodes ; ++a )
    {
      for( int component=0 ; component<3 ; ++component )
      {
        nodeDOF[3*a+component] = blockLocalDofNumber[faceToNodeMap(kf,a)] + component;
        nodeRHS[3*a+component] = - facePressure[kf] * faceNormal[kf][component] * faceArea[kf] / numNodes;
      }
    }
    rhs.add( nodeDOF, nodeRHS, numNodes*3 );
  }

}



void
SolidMechanicsLagrangianFEM::
ImplicitStepSetup( real64 const & GEOSX_UNUSED_ARG( time_n ),
                   real64 const & dt,
                   DomainPartition * const domain,
                   DofManager & dofManager,
                   ParallelMatrix & matrix,
                   ParallelVector & rhs,
                   ParallelVector & solution )
{
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  Group * const nodeManager = mesh->getNodeManager();

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

  SetupSystem( domain, dofManager, matrix, rhs, solution );
}

void SolidMechanicsLagrangianFEM::ImplicitStepComplete( real64 const & GEOSX_UNUSED_ARG( time_n ),
                                                        real64 const & dt,
                                                        DomainPartition * const domain)
{
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  Group * const nodeManager = mesh->getNodeManager();
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

void SolidMechanicsLagrangianFEM::SetupDofs( DomainPartition const * const GEOSX_UNUSED_ARG( domain ),
                                             DofManager & dofManager ) const
{
  dofManager.addField( keys::TotalDisplacement,
                       DofManager::Location::Node,
                       DofManager::Connectivity::Elem,
                       3 );
}

void SolidMechanicsLagrangianFEM::AssembleSystem( real64 const GEOSX_UNUSED_ARG( time_n ),
                                                  real64 const dt,
                                                  DomainPartition * const domain,
                                                  DofManager const & dofManager,
                                                  ParallelMatrix & matrix,
                                                  ParallelVector & rhs )
{
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  Group * const nodeManager = mesh->getNodeManager();
  ConstitutiveManager  * const constitutiveManager = domain->GetGroup<ConstitutiveManager >(keys::ConstitutiveManager);
  ElementRegionManager * const elemManager = mesh->getElemManager();
  NumericalMethodsManager const * numericalMethodManager = domain->getParent()->GetGroup<NumericalMethodsManager>(keys::numericalMethodsManager);
  FiniteElementDiscretizationManager const * feDiscretizationManager = numericalMethodManager->GetGroup<FiniteElementDiscretizationManager>(keys::finiteElementDiscretizations);

  ElementRegionManager::MaterialViewAccessor<real64> const biotCoefficient =
    elemManager->ConstructFullMaterialViewAccessor<real64>( "BiotCoefficient", constitutiveManager);

  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> const fluidPres =
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>("pressure");

  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> const dPres =
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>("deltaPressure");

  matrix.zero();
  rhs.zero();

  matrix.open();
  rhs.open();

  r1_array const& disp = nodeManager->getReference<r1_array>(keys::TotalDisplacement);
  r1_array const& uhat = nodeManager->getReference<r1_array>(keys::IncrementalDisplacement);
  //r1_array const& vel  = nodeManager->getReference<r1_array>(keys::Velocity);

  r1_array const uhattilde;
  r1_array const vtilde;

  string const dofKey = dofManager.getKey( keys::TotalDisplacement );

  globalIndex_array const & dofNumber = nodeManager->getReference<globalIndex_array>( dofKey );


  ElementRegionManager::ConstitutiveRelationAccessor<ConstitutiveBase>
  constitutiveRelations = elemManager->ConstructFullConstitutiveAccessor<ConstitutiveBase>(constitutiveManager);

  ElementRegionManager::MaterialViewAccessor< real64 > const
  density = elemManager->ConstructFullMaterialViewAccessor< real64 >( "density0",
                                                                  constitutiveManager );

  // begin region loop
  for( localIndex er=0 ; er<elemManager->numRegions() ; ++er )
  {
    ElementRegionBase * const elementRegion = elemManager->GetRegion(er);

    FiniteElementDiscretization const *
    feDiscretization = feDiscretizationManager->GetGroup<FiniteElementDiscretization>(m_discretizationName);

    elementRegion->forElementSubRegionsIndex<CellElementSubRegion>([&]( localIndex const esr,
                                                                        CellElementSubRegion const * const elementSubRegion )
    {
      array3d<R1Tensor> const &
      dNdX = elementSubRegion->getReference< array3d<R1Tensor> >(keys::dNdX);

      arrayView2d<real64> const & detJ = elementSubRegion->getReference< array2d<real64> >(keys::detJ);

      arrayView2d< localIndex > const & elemsToNodes = elementSubRegion->nodeList();
      localIndex const numNodesPerElement = elemsToNodes.size(1);

      std::unique_ptr<FiniteElementBase>
      fe = feDiscretization->getFiniteElement( elementSubRegion->GetElementTypeString() );

      // space for element matrix and rhs

      m_maxForce = ImplicitElementKernelLaunch( numNodesPerElement,
                                                fe->n_quadrature_points(),
                                                constitutiveRelations[er][esr][m_solidMaterialFullIndex],
                                                elementSubRegion->size(),
                                                dt,
                                                dNdX,
                                                detJ,
                                                fe.get(),
                                                elementSubRegion->m_ghostRank,
                                                elemsToNodes,
                                                dofNumber,
                                                disp,
                                                uhat,
                                                vtilde,
                                                uhattilde,
                                                density[er][esr],
                                                fluidPres[er][esr],
                                                dPres[er][esr],
                                                biotCoefficient[er][esr],
                                                m_timeIntegrationOption,
                                                this->m_stiffnessDamping,
                                                this->m_massDamping,
                                                this->m_newmarkBeta,
                                                this->m_newmarkGamma,
                                                &dofManager,
                                                &matrix,
                                                &rhs );

    });
  }


  matrix.close();
  rhs.close();

  if( verboseLevel() >= 2 )
  {
    GEOS_LOG_RANK_0( "After SolidMechanicsLagrangianFEM::AssembleSystem" );
    GEOS_LOG_RANK_0( "\nMatrix:\n" );
    matrix.print(std::cout);
    GEOS_LOG_RANK_0( "\nRhs:\n" );
    rhs.print(std::cout);
  }
}

void
SolidMechanicsLagrangianFEM::
ApplyBoundaryConditions( real64 const time_n,
                         real64 const dt,
                         DomainPartition * const domain,
                         DofManager const & dofManager,
                         ParallelMatrix & matrix,
                         ParallelVector & rhs )
{

  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);

  FaceManager * const faceManager = mesh->getFaceManager();
  FieldSpecificationManager * fsManager = FieldSpecificationManager::get();
//  fsManager->ApplyBoundaryCondition( this, &SolidMechanics_LagrangianFEM::ForceBC,
//                                     nodeManager, keys::Force, time_n + dt, *blockSystem );

  string const dofKey = dofManager.getKey( keys::TotalDisplacement );

  matrix.open();
  rhs.open();

  fsManager->Apply( time_n+dt,
                    domain,
                    "nodeManager",
                    keys::Force,
                    [&]( FieldSpecificationBase const * const bc,
                         string const &,
                         set< localIndex > const & targetSet,
                         Group * const targetGroup,
                         string const GEOSX_UNUSED_ARG( fieldName ) )
  {
    bc->ApplyBoundaryConditionToSystem<FieldSpecificationAdd, LAInterface>( targetSet,
                                                                            false,
                                                                            time_n + dt,
                                                                            targetGroup,
                                                                            keys::TotalDisplacement, // TODO fix use of dummy name for
                                                                            dofKey,
                                                                            3,
                                                                            matrix,
                                                                            rhs );
  });

  ApplyTractionBC(
    time_n + dt,
    dofManager, domain,
    rhs );

  ApplyDisplacementBC_implicit( time_n + dt, dofManager, *domain, matrix, rhs );

  if( faceManager->hasView("ChomboPressure") )
  {
    fsManager->ApplyFieldValue( time_n, domain, "faceManager", "ChomboPressure" );
    ApplyChomboPressure( dofManager, domain, rhs );
  }

  matrix.close();
  rhs.close();

  if( verboseLevel() >= 2 )
  {
    GEOS_LOG_RANK_0( "After SolidMechanicsLagrangianFEM::AssembleSystem" );
    GEOS_LOG_RANK_0("\nJacobian:\n");
    std::cout << matrix;
    GEOS_LOG_RANK_0("\nResidual:\n");
    std::cout << rhs;
  }


}

real64
SolidMechanicsLagrangianFEM::
CalculateResidualNorm( DomainPartition const * const GEOSX_UNUSED_ARG( domain ),
                       DofManager const & GEOSX_UNUSED_ARG( dofManager ),
                       ParallelVector const & rhs )
{
  real64 const * localResidual = rhs.extractLocalVector();

  real64 localResidualNorm[2] = { 0.0, this->m_maxForce };

  for( localIndex i=0 ; i<rhs.localSize() ; ++i )
  {
    localResidualNorm[0] += localResidual[i] * localResidual[i];
  }


  real64 globalResidualNorm[2] = {0,0};
//  MPI_Allreduce (&localResidual,&globalResidualNorm,1,MPI_DOUBLE,MPI_SUM ,MPI_COMM_GEOSX);


  int rank, size;
  MPI_Comm_rank(MPI_COMM_GEOSX, &rank);
  MPI_Comm_size(MPI_COMM_GEOSX, &size);
  array1d<real64> globalValues( size * 2 );
  globalValues = 0;
  MPI_Gather( localResidualNorm,
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



void
SolidMechanicsLagrangianFEM::ApplySystemSolution( DofManager const & dofManager,
                                                  ParallelVector const & solution,
                                                  real64 const scalingFactor,
                                                  DomainPartition * const domain )
{
  MeshLevel * const mesh = domain->getMeshBody( 0 )->getMeshLevel( 0 );
  NodeManager * const nodeManager = mesh->getNodeManager();

  string const fieldName = keys::TotalDisplacement;
  dofManager.addVectorToField( solution, fieldName, -scalingFactor, nodeManager, keys::IncrementalDisplacement );
  dofManager.addVectorToField( solution, fieldName, -scalingFactor, nodeManager, keys::TotalDisplacement );

  std::map<string, string_array > fieldNames;
  fieldNames["node"].push_back( keys::IncrementalDisplacement );
  fieldNames["node"].push_back( keys::TotalDisplacement );

  CommunicationTools::SynchronizeFields( fieldNames, mesh,
                                         domain->getReference< array1d<NeighborCommunicator> >( domain->viewKeys.neighbors ) );
}

void SolidMechanicsLagrangianFEM::SolveSystem( DofManager const & dofManager,
                                               ParallelMatrix & matrix,
                                               ParallelVector & rhs,
                                               ParallelVector & solution )
{
  solution.zero();

  SolverBase::SolveSystem( dofManager, matrix, rhs, solution );

  if( verboseLevel() >= 2 )
  {
    GEOS_LOG_RANK_0( "After SolidMechanicsLagrangianFEM::SolveSystem" );
    GEOS_LOG_RANK_0( "\nSolution:\n");
    std::cout << solution;
  }

}

void SolidMechanicsLagrangianFEM::ResetStateToBeginningOfStep( DomainPartition * const domain )
{
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  NodeManager * const nodeManager = mesh->getNodeManager();

  r1_array& incdisp  = nodeManager->getReference<r1_array>(keys::IncrementalDisplacement);
  r1_array& disp     = nodeManager->getReference<r1_array>(keys::TotalDisplacement);

  // TODO need to finish this rewind
  forall_in_range(0, nodeManager->size(), GEOSX_LAMBDA (localIndex a) mutable
  {
    disp[a] -= incdisp[a];
    incdisp[a] = 0.0;
  });
}




REGISTER_CATALOG_ENTRY( SolverBase, SolidMechanicsLagrangianFEM, string const &, dataRepository::Group * const )
}
