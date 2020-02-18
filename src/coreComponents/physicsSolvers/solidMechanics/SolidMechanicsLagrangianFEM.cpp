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

#include "codingUtilities/Utilities.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/contact/ContactRelationBase.hpp"
#include "finiteElement/ElementLibrary/FiniteElement.h"
#include "finiteElement/FiniteElementDiscretizationManager.hpp"
#include "finiteElement/Kinematics.h"
#include "managers/DomainPartition.hpp"
#include "managers/FieldSpecification/FieldSpecificationManager.hpp"
#include "managers/NumericalMethodsManager.hpp"
#include "mesh/FaceElementSubRegion.hpp"
#include "meshUtilities/ComputationalGeometry.hpp"
#include "mpiCommunications/CommunicationTools.hpp"
#include "mpiCommunications/NeighborCommunicator.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"


namespace geosx
{

using namespace dataRepository;
using namespace constitutive;

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
  m_sendOrReceiveNodes.setName("SolidMechanicsLagrangianFEM::m_sendOrReceiveNodes");
  m_nonSendOrReceiveNodes.setName("SolidMechanicsLagrangianFEM::m_nonSendOrReceiveNodes");

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

  registerWrapper(viewKeyStruct::contactRelationNameString, &m_contactRelationName, 0)->
    setApplyDefaultValue(viewKeyStruct::noContactRelationNameString)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Name of contact relation to enforce constraints on fracture boundary.");

  registerWrapper(viewKeyStruct::maxForce, &m_maxForce, false )->
    setInputFlag(InputFlags::FALSE)->
    setDescription( "The maximum force contribution in the problem domain.");
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


    nodes->registerWrapper< array2d< real64, nodes::TOTAL_DISPLACEMENT_PERM > >( keys::TotalDisplacement )->
      setPlotLevel(PlotLevel::LEVEL_0)->
      setRegisteringObjects(this->getName())->
      setDescription( "An array that holds the total displacements on the nodes.");
    nodes->totalDisplacement().resizeDimension< 1 >( 3 );

    nodes->registerWrapper< array2d< real64, nodes::INCR_DISPLACEMENT_PERM > >( keys::IncrementalDisplacement )->
      setPlotLevel(PlotLevel::LEVEL_3)->
      setRegisteringObjects(this->getName())->
      setDescription( "An array that holds the incremental displacements for the current time step on the nodes.");
    nodes->incrementalDisplacement().resizeDimension< 1 >( 3 );

    nodes->registerWrapper< array2d< real64, nodes::VELOCITY_PERM > >( keys::Velocity )->
      setPlotLevel(PlotLevel::LEVEL_0)->
      setRegisteringObjects(this->getName())->
      setDescription( "An array that holds the current velocity on the nodes.");
    nodes->velocity().resizeDimension< 1 >( 3 );

    nodes->registerWrapper< array2d< real64, nodes::ACCELERATION_PERM > >( keys::Acceleration )->
      setPlotLevel(PlotLevel::LEVEL_1)->
      setRegisteringObjects(this->getName())->
      setDescription( "An array that holds the current acceleration on the nodes. This array also is used "
                      "to hold the summation of nodal forces resulting from the governing equations.");
    nodes->acceleration().resizeDimension< 1 >( 3 );

    nodes->registerWrapper<array1d<R1Tensor> >( viewKeyStruct::forceExternal )->
      setPlotLevel(PlotLevel::LEVEL_0)->
      setRegisteringObjects(this->getName())->
      setDescription( "An array that holds the external forces on the nodes. This includes any boundary"
                      " conditions as well as coupling forces such as hydraulic forces.");

    nodes->registerWrapper<array1d<real64> >( keys::Mass )->
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

    nodes->registerWrapper<array1d<R1Tensor> >( viewKeyStruct::contactForceString )->
      setPlotLevel(PlotLevel::LEVEL_0)->
      setRegisteringObjects(this->getName())->
      setDescription( "An array that holds the contact force.");

    ElementRegionManager * const
    elementRegionManager = mesh.second->group_cast<MeshBody*>()->getMeshLevel(0)->getElemManager();
    elementRegionManager->forElementRegions<CellElementRegion>( m_targetRegions,
                                                                [&]( CellElementRegion * const elemRegion )
    {
      elemRegion->forElementSubRegions<CellElementSubRegion>([&]( CellElementSubRegion * const subRegion )
      {
        subRegion->registerWrapper<array2d<R2SymTensor> >( viewKeyStruct::stress_n )->
          setPlotLevel(PlotLevel::NOPLOT)->
          setRestartFlags(RestartFlags::NO_WRITE)->
          setRegisteringObjects(this->getName())->
          setDescription("Array to hold the beginning of step stress for implicit problem rewinds");
      });
    });


  }
}


void SolidMechanicsLagrangianFEM::InitializePreSubGroups(Group * const rootGroup)
{
  SolverBase::InitializePreSubGroups(rootGroup);

  DomainPartition * domain = rootGroup->GetGroup<DomainPartition>(keys::domain);
  ConstitutiveManager const * const cm = domain->getConstitutiveManager();

  ConstitutiveBase const * const solid  = cm->GetConstitutiveRelation<ConstitutiveBase>( m_solidMaterialName );
  GEOSX_ERROR_IF( solid == nullptr, "constitutive model " + m_solidMaterialName + " not found" );
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
      arrayView2d<localIndex const, cells::NODE_MAP_USD> const & elemsToNodes = elementSubRegion->nodeList();

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
      m_elemsAttachedToSendOrReceiveNodes[er][esr].setName(
        "SolidMechanicsLagrangianFEM::m_elemsAttachedToSendOrReceiveNodes["
        + std::to_string(er) + "][" + std::to_string(esr) + "]" );

      m_elemsNotAttachedToSendOrReceiveNodes[er][esr].setName(
        "SolidMechanicsLagrangianFEM::m_elemsNotAttachedToSendOrReceiveNodes["
        + std::to_string(er) + "][" + std::to_string(esr) + "]" );

      arrayView2d<real64> const & detJ = elementSubRegion->getReference< array2d<real64> >(keys::detJ);
      arrayView2d<localIndex const, cells::NODE_MAP_USD> const & elemsToNodes = elementSubRegion->nodeList();

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
  GEOSX_MARK_FUNCTION;
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
    ImplicitStepSetup( time_n, dt, domain, m_dofManager, m_matrix, m_rhs, m_solution );
    for( int solveIter=0 ; solveIter<maxNumResolves ; ++solveIter )
    {
      SetupSystem( domain, m_dofManager, m_matrix, m_rhs, m_solution );

      if( solveIter>0 )
      {
        ResetStressToBeginningOfStep(domain);
      }

      dtReturn = NonlinearImplicitStep( time_n, dt, cycleNumber, domain->group_cast<DomainPartition *>(), m_dofManager,
                                        m_matrix, m_rhs, m_solution );

      updateStress( domain );
      if( surfaceGenerator!=nullptr )
      {
        if( surfaceGenerator->SolverStep( time_n, dt, cycleNumber, domain ) > 0 )
        {
          locallyFractured = 1;
        }
        MpiWrapper::allReduce( &locallyFractured,
                               &globallyFractured,
                               1,
                               MPI_MAX,
                               MPI_COMM_GEOSX);
      }
      if( globallyFractured == 0 )
      {
        break;
      }
      else
      {
        GEOSX_LOG_RANK_0("Fracture Occurred. Resolve");
      }
    }
    ImplicitStepComplete( time_n, dt,  domain );
  }

  return dtReturn;
}

real64 SolidMechanicsLagrangianFEM::ExplicitStep( real64 const& time_n,
                                                  real64 const& dt,
                                                  const int GEOSX_UNUSED_PARAM( cycleNumber ),
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

  FieldSpecificationManager & fsManager = FieldSpecificationManager::get();

  arrayView1d<real64 const> const & mass = nodes->getReference<array1d<real64>>(keys::Mass);
  arrayView2d<real64, nodes::VELOCITY_USD> const & vel = nodes->velocity();

  arrayView2d<real64, nodes::TOTAL_DISPLACEMENT_USD> const & u = nodes->totalDisplacement();
  arrayView2d<real64, nodes::INCR_DISPLACEMENT_USD> const & uhat = nodes->incrementalDisplacement();
  arrayView2d<real64, nodes::ACCELERATION_USD> const & acc = nodes->acceleration();

  array1d<NeighborCommunicator> & neighbors = domain->getReference< array1d<NeighborCommunicator> >( domain->viewKeys.neighbors );
  std::map<string, string_array > fieldNames;
  fieldNames["node"].push_back( keys::Velocity);
  fieldNames["node"].push_back( keys::Acceleration);

  CommunicationTools::SynchronizePackSendRecvSizes( fieldNames, mesh, neighbors, m_iComm, true );

  fsManager.ApplyFieldValue< parallelDevicePolicy< 1024 > >( time_n, domain, "nodeManager", keys::Acceleration );

  //3: v^{n+1/2} = v^{n} + a^{n} dt/2
  SolidMechanicsLagrangianFEMKernels::velocityUpdate( acc, vel, dt/2 );

  fsManager.ApplyFieldValue< parallelDevicePolicy< 1024 > >( time_n, domain, "nodeManager", keys::Velocity );

  //4. x^{n+1} = x^{n} + v^{n+{1}/{2}} dt (x is displacement)
  SolidMechanicsLagrangianFEMKernels::displacementUpdate( vel, uhat, u, dt );

  fsManager.ApplyFieldValue( time_n + dt, domain, "nodeManager", keys::TotalDisplacement,
    [&]( FieldSpecificationBase const * const bc, SortedArrayView<localIndex const> const & targetSet )
    {
      integer const component = bc->GetComponent();
      forall_in_range< parallelDevicePolicy< 1024 > >(0, targetSet.size(),
        GEOSX_DEVICE_LAMBDA( localIndex const i )
        {
          localIndex const a = targetSet[ i ];
          vel( a, component ) = u( a, component );
        }
      );
    },
    [&]( FieldSpecificationBase const * const bc, SortedArrayView<localIndex const> const & targetSet )
    {
      integer const component = bc->GetComponent();
      forall_in_range< parallelDevicePolicy< 1024 > >(0, targetSet.size(),
        GEOSX_DEVICE_LAMBDA( localIndex const i )
        {
          localIndex const a = targetSet[ i ];
          uhat( a, component ) = u( a, component ) - vel( a, component );
          vel( a, component )  = uhat( a, component ) / dt;
        }
      );
    }
  );

  ElementRegionManager::MaterialViewAccessor< arrayView3d<real64, solid::STRESS_USD> > const
  stress = elemManager->ConstructFullMaterialViewAccessor< array3d<real64, solid::STRESS_PERMUTATION>,
                                                           arrayView3d<real64, solid::STRESS_USD> >( SolidBase::viewKeyStruct::stressString,
                                                                                       constitutiveManager);

  ElementRegionManager::ConstitutiveRelationAccessor<ConstitutiveBase>
  constitutiveRelations = elemManager->ConstructFullConstitutiveAccessor<ConstitutiveBase>(constitutiveManager);

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

      arrayView2d<localIndex const, cells::NODE_MAP_USD> const & elemsToNodes = elementSubRegion->nodeList();

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
                                   stress[er][esr][m_solidMaterialFullIndex],
                                   dt );

    }); //Element Region

  } //Element Manager

  // apply this over a set
  SolidMechanicsLagrangianFEMKernels::velocityUpdate( acc, mass, vel, dt / 2, m_sendOrReceiveNodes );

  fsManager.ApplyFieldValue< parallelDevicePolicy< 1024 > >( time_n, domain, "nodeManager", keys::Velocity );

  CommunicationTools::SynchronizePackSendRecv( fieldNames, mesh, neighbors, m_iComm, true );

  for( localIndex er=0 ; er<elemManager->numRegions() ; ++er )
  {
    ElementRegionBase * const elementRegion = elemManager->GetRegion(er);

    FiniteElementDiscretization const * feDiscretization = feDiscretizationManager->GetGroup<FiniteElementDiscretization>(m_discretizationName);

    elementRegion->forElementSubRegionsIndex<CellElementSubRegion>([&]( localIndex const esr, CellElementSubRegion const * const elementSubRegion )
    {
      arrayView3d< R1Tensor > const & dNdX = elementSubRegion->getReference< array3d< R1Tensor > >(keys::dNdX);

      arrayView2d<real64> const & detJ = elementSubRegion->getReference< array2d<real64> >(keys::detJ);

      arrayView2d<localIndex const, cells::NODE_MAP_USD> const & elemsToNodes = elementSubRegion->nodeList();

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
                                   stress[er][esr][m_solidMaterialFullIndex],
                                   dt );
    }); //Element Region

  } //Element Manager

  // apply this over a set
  SolidMechanicsLagrangianFEMKernels::velocityUpdate( acc, mass, vel, dt / 2, m_nonSendOrReceiveNodes );

  fsManager.ApplyFieldValue< parallelDevicePolicy< 1024 > >( time_n, domain, "nodeManager", keys::Velocity );

  CommunicationTools::SynchronizeUnpack( mesh, neighbors, m_iComm, true );

  return dt;
}



void SolidMechanicsLagrangianFEM::ApplyDisplacementBC_implicit( real64 const time,
                                                                DofManager const & dofManager,
                                                                DomainPartition & domain,
                                                                ParallelMatrix & matrix,
                                                                ParallelVector & rhs )
{
  string const dofKey = dofManager.getKey( keys::TotalDisplacement );

  FieldSpecificationManager const & fsManager = FieldSpecificationManager::get();

  fsManager.Apply( time,
                   &domain,
                   "nodeManager",
                   keys::TotalDisplacement,
                   [&]( FieldSpecificationBase const * const bc,
                        string const &,
                        SortedArray<localIndex> const & targetSet,
                        Group * const targetGroup,
                        string const fieldName )
  {
    bc->ApplyBoundaryConditionToSystem<FieldSpecificationEqual, LAInterface>( targetSet,
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
  FieldSpecificationManager & fsManager = FieldSpecificationManager::get();
  FunctionManager & functionManager = FunctionManager::Instance();

  FaceManager * const faceManager = domain->getMeshBody(0)->getMeshLevel(0)->getFaceManager();
  NodeManager * const nodeManager = domain->getMeshBody(0)->getMeshLevel(0)->getNodeManager();

  real64_array const & faceArea  = faceManager->getReference<real64_array>("faceArea");
  ArrayOfArraysView< localIndex const > const & faceToNodeMap = faceManager->nodeList();

  string const dofKey = dofManager.getKey( keys::TotalDisplacement );

  arrayView1d<globalIndex> const &
  blockLocalDofNumber = nodeManager->getReference<globalIndex_array>( dofKey );

  arrayView1d<integer const> const & faceGhostRank = faceManager->GhostRank();
  fsManager.Apply( time,
                   domain,
                   "faceManager",
                   string("Traction"),
                   [&]( FieldSpecificationBase const * const bc,
                        string const &,
                        SortedArray<localIndex> const & targetSet,
                        Group * const GEOSX_UNUSED_PARAM( targetGroup ),
                        string const GEOSX_UNUSED_PARAM( fieldName ) )
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
      FunctionBase const * const function = functionManager.GetGroup<FunctionBase>(functionName);
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
ImplicitStepSetup( real64 const & GEOSX_UNUSED_PARAM( time_n ),
                   real64 const & dt,
                   DomainPartition * const domain,
                   DofManager & GEOSX_UNUSED_PARAM( dofManager ),
                   ParallelMatrix & GEOSX_UNUSED_PARAM( matrix ),
                   ParallelVector & GEOSX_UNUSED_PARAM( rhs ),
                   ParallelVector & GEOSX_UNUSED_PARAM( solution ) )
{
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  NodeManager * const nodeManager = mesh->getNodeManager();

  arrayView2d<real64 const, nodes::VELOCITY_USD> const & v_n = nodeManager->velocity();
  arrayView2d<real64, nodes::INCR_DISPLACEMENT_USD> const & uhat = nodeManager->incrementalDisplacement();
  arrayView2d<real64, nodes::TOTAL_DISPLACEMENT_USD> const & disp = nodeManager->totalDisplacement();

  if( this->m_timeIntegrationOption == timeIntegrationOption::ImplicitDynamic )
  {
    arrayView2d<real64 const, nodes::ACCELERATION_USD> const & a_n = nodeManager->acceleration();
    arrayView1d<R1Tensor> const & vtilde   = nodeManager->getReference<array1d<R1Tensor>>(solidMechanicsViewKeys.vTilde);
    arrayView1d<R1Tensor> const & uhatTilde   = nodeManager->getReference<array1d<R1Tensor>>(solidMechanicsViewKeys.uhatTilde);

    localIndex const numNodes = nodeManager->size();
    real64 const newmarkGamma = this->getReference<real64>(solidMechanicsViewKeys.newmarkGamma);
    real64 const newmarkBeta = this->getReference<real64>(solidMechanicsViewKeys.newmarkBeta);

    RAJA::forall< parallelHostPolicy >( RAJA::TypedRangeSegment< localIndex >( 0, numNodes ),
                                        GEOSX_LAMBDA ( localIndex const a )
    {
      for( int i=0 ; i<3 ; ++i )
      {
        vtilde[a][i] = v_n( a, i ) + (1.0-newmarkGamma) * a_n( a, i ) * dt;
        uhatTilde[a][i] = ( v_n( a, i ) + 0.5 * ( 1.0 - 2.0*newmarkBeta ) * a_n( a, i ) * dt ) *dt;
        uhat( a, i ) = uhatTilde[a][i];
        disp( a, i ) += uhatTilde[a][i];
      }
    } );
  }
  else if( this->m_timeIntegrationOption == timeIntegrationOption::QuasiStatic  )
  {
    localIndex const numNodes = nodeManager->size();

    if( m_useVelocityEstimateForQS==1 )
    {
      RAJA::forall< parallelHostPolicy >( RAJA::TypedRangeSegment< localIndex >( 0, numNodes ),
                                          GEOSX_LAMBDA ( localIndex const a )
      {
        for( int i=0 ; i<3 ; ++i )
        {
          uhat( a, i ) = v_n( a, i ) * dt;
          disp( a, i ) += uhat( a, i );
        }
      } );
    }
    else
    {
      RAJA::forall< parallelHostPolicy >( RAJA::TypedRangeSegment< localIndex >( 0, numNodes ),
                                          GEOSX_LAMBDA ( localIndex const a )
      {
        for( int i=0 ; i<3 ; ++i )
        {
          uhat( a, i ) = 0.0;
        }
      } );
    }
  }

  ElementRegionManager * const elementRegionManager = mesh->getElemManager();
  ConstitutiveManager  * const
  constitutiveManager = domain->GetGroup<ConstitutiveManager >(dataRepository::keys::ConstitutiveManager);
  ElementRegionManager::ConstitutiveRelationAccessor<ConstitutiveBase>
  constitutiveRelations = elementRegionManager->ConstructFullConstitutiveAccessor<ConstitutiveBase>(constitutiveManager);

  elementRegionManager->
  forElementSubRegionsComplete<CellElementSubRegion>( m_targetRegions,
                                                      [&]( localIndex const er,
                                                           localIndex const esr,
                                                           ElementRegionBase * const,
                                                           CellElementSubRegion * const subRegion )
  {
    SolidBase * const
    constitutiveRelation = constitutiveRelations[er][esr][m_solidMaterialFullIndex]->group_cast<SolidBase*>();

    arrayView3d<real64 const, solid::STRESS_USD> const & stress = constitutiveRelation->getStress();

    array2d<R2SymTensor> &
    stress_n = subRegion->getReference<array2d<R2SymTensor>>(viewKeyStruct::stress_n);
    stress_n.resize( stress.size(0), stress.size(1) );

    for( localIndex k=0 ; k<stress.size(0) ; ++k )
    {
      for( localIndex a=0 ; a<stress.size(1) ; ++a )
      {
        stress_n(k,a) = stress[k][a];
      }
    }
  });



}

void SolidMechanicsLagrangianFEM::ImplicitStepComplete( real64 const & GEOSX_UNUSED_PARAM( time_n ),
                                                        real64 const & dt,
                                                        DomainPartition * const domain)
{
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  NodeManager * const nodeManager = mesh->getNodeManager();
  localIndex const numNodes = nodeManager->size();

  arrayView2d< real64, nodes::VELOCITY_USD > const & v_n = nodeManager->velocity();
  arrayView2d< real64 const, nodes::INCR_DISPLACEMENT_USD > const & uhat  = nodeManager->incrementalDisplacement();

  if( this->m_timeIntegrationOption == timeIntegrationOption::ImplicitDynamic )
  {
    arrayView2d< real64, nodes::ACCELERATION_USD > const & a_n = nodeManager->acceleration();
    arrayView1d< R1Tensor const > const & vtilde    = nodeManager->getReference<r1_array>(solidMechanicsViewKeys.vTilde);
    arrayView1d< R1Tensor const > const & uhatTilde = nodeManager->getReference<r1_array>(solidMechanicsViewKeys.uhatTilde);
    real64 const newmarkGamma = this->getReference<real64>(solidMechanicsViewKeys.newmarkGamma);
    real64 const newmarkBeta = this->getReference<real64>(solidMechanicsViewKeys.newmarkBeta);

    RAJA::forall< parallelHostPolicy >( RAJA::TypedRangeSegment< localIndex >( 0, numNodes ),
                                        GEOSX_LAMBDA ( localIndex const a )
    {
      for( int i=0 ; i<3 ; ++i )
      {
        a_n( a, i ) = 1.0 / ( newmarkBeta * dt*dt) * ( uhat( a, i ) - uhatTilde[a][i] );
        v_n[a][i] = vtilde[a][i] + newmarkGamma * a_n( a, i ) * dt;
      }
    } );
  }
  else if( this->m_timeIntegrationOption == timeIntegrationOption::QuasiStatic && dt > 0.0)
  {
    RAJA::forall< parallelHostPolicy >( RAJA::TypedRangeSegment< localIndex >( 0, numNodes ),
                                        GEOSX_LAMBDA ( localIndex const a )
    {
      for( int i=0 ; i<3 ; ++i )
      {
        v_n[a][i] = uhat( a, i ) / dt;
      }
    } );
  }
}

void SolidMechanicsLagrangianFEM::SetupDofs( DomainPartition const * const GEOSX_UNUSED_PARAM( domain ),
                                             DofManager & dofManager ) const
{
  dofManager.addField( keys::TotalDisplacement,
                       DofManager::Location::Node,
                       3 );

  dofManager.addCoupling( keys::TotalDisplacement,
                          keys::TotalDisplacement,
                          DofManager::Connectivity::Elem );
}

void SolidMechanicsLagrangianFEM::AssembleSystem( real64 const GEOSX_UNUSED_PARAM( time_n ),
                                                  real64 const dt,
                                                  DomainPartition * const domain,
                                                  DofManager const & dofManager,
                                                  ParallelMatrix & matrix,
                                                  ParallelVector & rhs )
{
  GEOSX_MARK_FUNCTION;
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  NodeManager const * const nodeManager = mesh->getNodeManager();
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

  arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const & disp = nodeManager->totalDisplacement();
  arrayView2d< real64 const, nodes::INCR_DISPLACEMENT_USD > const & uhat = nodeManager->incrementalDisplacement();

  r1_array const uhattilde;
  r1_array const vtilde;

  string const dofKey = dofManager.getKey( keys::TotalDisplacement );

  globalIndex_array const & dofNumber = nodeManager->getReference<globalIndex_array>( dofKey );


  ElementRegionManager::ConstitutiveRelationAccessor<ConstitutiveBase>
  constitutiveRelations = elemManager->ConstructFullConstitutiveAccessor<ConstitutiveBase>(constitutiveManager);

  ElementRegionManager::MaterialViewAccessor< arrayView2d<real64> > const
  density = elemManager->ConstructFullMaterialViewAccessor< array2d<real64>,
                                                            arrayView2d<real64> >( SolidBase::viewKeyStruct::densityString,
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

      arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes = elementSubRegion->nodeList();
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
                                                density[er][esr][m_solidMaterialFullIndex],
                                                fluidPres[er][esr],
                                                dPres[er][esr],
                                                biotCoefficient[er][esr],
                                                m_timeIntegrationOption,
                                                this->m_stiffnessDamping,
                                                this->m_massDamping,
                                                this->m_newmarkBeta,
                                                this->m_newmarkGamma,
                                                gravityVector(),
                                                &dofManager,
                                                &matrix,
                                                &rhs );

    });
  }


  ApplyContactConstraint( dofManager,
                          *domain,
                          &matrix,
                          &rhs );

  matrix.close();
  rhs.close();

  // Debug for logLevel >= 2
  GEOSX_LOG_LEVEL_RANK_0( 2, "After SolidMechanicsLagrangianFEM::AssembleSystem" );
  GEOSX_LOG_LEVEL_RANK_0( 2, "\nJacobian:\n" << matrix );
  GEOSX_LOG_LEVEL_RANK_0( 2, "\nResidual:\n" << rhs );
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
  GEOSX_MARK_FUNCTION;
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);

  FaceManager * const faceManager = mesh->getFaceManager();
  FieldSpecificationManager & fsManager = FieldSpecificationManager::get();

  string const dofKey = dofManager.getKey( keys::TotalDisplacement );

  matrix.open();
  rhs.open();

  fsManager.Apply( time_n + dt,
                   domain,
                   "nodeManager",
                   keys::Force,
                   [&]( FieldSpecificationBase const * const bc,
                        string const &,
                        SortedArray< localIndex > const & targetSet,
                        Group * const targetGroup,
                        string const GEOSX_UNUSED_PARAM( fieldName ) )
  {
    bc->ApplyBoundaryConditionToSystem<FieldSpecificationAdd, LAInterface>( targetSet,
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

  if( faceManager->hasWrapper( "ChomboPressure" ) )
  {
    fsManager.ApplyFieldValue( time_n, domain, "faceManager", "ChomboPressure" );
    ApplyChomboPressure( dofManager, domain, rhs );
  }

  matrix.close();
  rhs.close();

  // Debug for logLevel >= 2
  GEOSX_LOG_LEVEL_RANK_0( 2, "After SolidMechanicsLagrangianFEM::AssembleSystem" );
  GEOSX_LOG_LEVEL_RANK_0( 2, "\nJacobian:\n" << matrix );
  GEOSX_LOG_LEVEL_RANK_0( 2, "\nResidual:\n" << rhs );


}

real64
SolidMechanicsLagrangianFEM::
CalculateResidualNorm( DomainPartition const * const GEOSX_UNUSED_PARAM( domain ),
                       DofManager const & GEOSX_UNUSED_PARAM( dofManager ),
                       ParallelVector const & rhs )
{
  GEOSX_MARK_FUNCTION;
  real64 const * localResidual = rhs.extractLocalVector();

  real64 localResidualNorm[2] = { 0.0, this->m_maxForce };
  //real64 localResInfNorm[2] = {}

  for( localIndex i=0 ; i<rhs.localSize() ; ++i )
  {
  // sum(rhs^2) on each rank.
    localResidualNorm[0] += localResidual[i] * localResidual[i];
  }


  // globalResidualNorm[0]: the sum of all the local sum(rhs^2).
  // globalResidualNorm[1]: max of max force of each rank. Basically max force globally
  real64 globalResidualNorm[2] = {0,0};

  int const rank = MpiWrapper::Comm_rank(MPI_COMM_GEOSX);
  int const size = MpiWrapper::Comm_size(MPI_COMM_GEOSX);
  array1d<real64> globalValues( size * 2 );
  globalValues = 0;

  // Everything is done on rank 0
  MpiWrapper::gather( localResidualNorm,
                      2,
                      globalValues.data(),
                      2,
                      0,
                      MPI_COMM_GEOSX );

  if( rank==0 )
  {
    for( int r=0 ; r<size ; ++r )
    {
      // sum across all ranks
      globalResidualNorm[0] += globalValues[r*2];

      // check if it is greater than the other ranks.
      // If yes, change the entry of globalResidualNorm[1] (new max)
      if( globalResidualNorm[1] < globalValues[r*2+1] )
      {
        globalResidualNorm[1] = globalValues[r*2+1];
      }
    }
  }

  MpiWrapper::bcast( globalResidualNorm, 2, 0, MPI_COMM_GEOSX );


  real64 const residual = sqrt(globalResidualNorm[0])/(globalResidualNorm[1]+1);  // the + 1 is for the first time-step when maxForce = 0;

  if( getLogLevel() >= 1 && logger::internal::rank==0 )
  {
    char output[200] = {0};
    sprintf( output,
             "( RSolid ) = (%4.2e) ; ",
             residual);
    std::cout<<output;
  }


  return residual;
}



void
SolidMechanicsLagrangianFEM::ApplySystemSolution( DofManager const & dofManager,
                                                  ParallelVector const & solution,
                                                  real64 const scalingFactor,
                                                  DomainPartition * const domain )
{
  dofManager.addVectorToField( solution, keys::TotalDisplacement, keys::IncrementalDisplacement, -scalingFactor );
  dofManager.addVectorToField( solution, keys::TotalDisplacement, keys::TotalDisplacement, -scalingFactor );

  std::map<string, string_array > fieldNames;
  fieldNames["node"].push_back( keys::IncrementalDisplacement );
  fieldNames["node"].push_back( keys::TotalDisplacement );

  CommunicationTools::SynchronizeFields( fieldNames,
                                         domain->getMeshBody( 0 )->getMeshLevel( 0 ),
                                         domain->getReference< array1d<NeighborCommunicator> >( domain->viewKeys.neighbors ) );
}

void SolidMechanicsLagrangianFEM::SolveSystem( DofManager const & dofManager,
                                               ParallelMatrix & matrix,
                                               ParallelVector & rhs,
                                               ParallelVector & solution )
{
  solution.zero();

  SolverBase::SolveSystem( dofManager, matrix, rhs, solution );

  // Debug for logLevel >= 2
  GEOSX_LOG_LEVEL_RANK_0( 2, "After SolidMechanicsLagrangianFEM::SolveSystem" );
  GEOSX_LOG_LEVEL_RANK_0( 2, "\nSolution:\n" << solution );
}

void SolidMechanicsLagrangianFEM::ResetStateToBeginningOfStep( DomainPartition * const domain )
{
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  NodeManager * const nodeManager = mesh->getNodeManager();

  arrayView2d<real64, nodes::INCR_DISPLACEMENT_USD> const & incdisp  = nodeManager->incrementalDisplacement();
  arrayView2d<real64, nodes::TOTAL_DISPLACEMENT_USD> const & disp = nodeManager->totalDisplacement();

  // TODO need to finish this rewind
  forAll< serialPolicy >( nodeManager->size(), GEOSX_LAMBDA (localIndex const a)
  {
    for ( localIndex i = 0; i < 3; ++i )
    {
      disp( a, i ) -= incdisp( a, i );
      incdisp( a, i ) = 0.0;
    }
  });

  ResetStressToBeginningOfStep(domain);
}

void SolidMechanicsLagrangianFEM::ResetStressToBeginningOfStep( DomainPartition * const domain )
{
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);

  ElementRegionManager * const elementRegionManager = mesh->getElemManager();
  ConstitutiveManager  * const
  constitutiveManager = domain->GetGroup<ConstitutiveManager >(dataRepository::keys::ConstitutiveManager);
  ElementRegionManager::ConstitutiveRelationAccessor<ConstitutiveBase>
  constitutiveRelations = elementRegionManager->ConstructFullConstitutiveAccessor<ConstitutiveBase>(constitutiveManager);

  elementRegionManager->
  forElementSubRegionsComplete<CellElementSubRegion>( m_targetRegions,
                                                      [&]( localIndex const er,
                                                           localIndex const esr,
                                                           ElementRegionBase * const,
                                                           CellElementSubRegion * const subRegion )
  {
    SolidBase * const
    constitutiveRelation = constitutiveRelations[er][esr][m_solidMaterialFullIndex]->group_cast<SolidBase*>();

    arrayView3d<real64, solid::STRESS_USD> const & stress = constitutiveRelation->getStress();

    array2d<R2SymTensor> &
    stress_n = subRegion->getReference<array2d<R2SymTensor>>(viewKeyStruct::stress_n);

    for( localIndex k=0 ; k<stress.size(0) ; ++k )
    {
      for( localIndex a=0 ; a<stress.size(1) ; ++a )
      {
        for ( localIndex i = 0; i < 6; ++i )
        {
          stress( k, a, i ) = stress_n( k, a ).Data()[ i ];
        }
      }
    }
  });
}


void SolidMechanicsLagrangianFEM::ApplyContactConstraint( DofManager const & dofManager,
                                                          DomainPartition & domain,
                                                          ParallelMatrix * const matrix,
                                                          ParallelVector * const rhs )
{
  GEOSX_MARK_FUNCTION;

  if( m_contactRelationName != viewKeyStruct::noContactRelationNameString )
  {
    MeshLevel * const mesh = domain.getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
    FaceManager const * const faceManager = mesh->getFaceManager();
    NodeManager * const nodeManager = mesh->getNodeManager();
    ElementRegionManager * const elemManager = mesh->getElemManager();


    ConstitutiveManager const * const
    constitutiveManager = domain.GetGroup<ConstitutiveManager>(keys::ConstitutiveManager);

    ContactRelationBase const * const
    contactRelation = constitutiveManager->GetGroup<ContactRelationBase>(m_contactRelationName);

    real64 const contactStiffness = contactRelation->stiffness();

    arrayView2d<real64 const, nodes::TOTAL_DISPLACEMENT_USD> const & u = nodeManager->totalDisplacement();
    arrayView1d<R1Tensor> const & fc = nodeManager->getReference< array1d<R1Tensor> >( viewKeyStruct::contactForceString );
    fc = {0,0,0};

    arrayView1d<R1Tensor const> const & faceNormal = faceManager->faceNormal();
    ArrayOfArraysView<localIndex const> const & facesToNodes = faceManager->nodeList();

    string const dofKey = dofManager.getKey( keys::TotalDisplacement );
    arrayView1d<globalIndex> const & nodeDofNumber = nodeManager->getReference<globalIndex_array>( dofKey );


    elemManager->forElementSubRegions<FaceElementSubRegion>([&]( FaceElementSubRegion * const subRegion )->void
    {
        arrayView1d<integer const> const & ghostRank = subRegion->GhostRank();
        arrayView1d<real64> const & area = subRegion->getElementArea();
        arrayView2d< localIndex const > const & elemsToFaces = subRegion->faceList();

        forall_in_range<serialPolicy>( 0,
                                       subRegion->size(),
                                       GEOSX_LAMBDA ( localIndex const kfe )
        {

          if( ghostRank[kfe] < 0 )
          {
            R1Tensor Nbar = faceNormal[elemsToFaces[kfe][0]];
            Nbar -= faceNormal[elemsToFaces[kfe][1]];
            Nbar.Normalize();

            localIndex const kf0 = elemsToFaces[kfe][0];
            localIndex const kf1 = elemsToFaces[kfe][1];
            localIndex const numNodesPerFace=facesToNodes.sizeOfArray(kf0);
            real64 const Ja = area[kfe] / numNodesPerFace;


            globalIndex rowDOF[24] = {0};
            real64 nodeRHS[24] = {0};
            stackArray2d<real64, (4*3*2)*(4*3*2)> dRdP(numNodesPerFace*3*2, numNodesPerFace*3*2);

            for( localIndex a=0 ; a<numNodesPerFace ; ++a )
            {
              R1Tensor penaltyForce = Nbar;
              localIndex const node0 = facesToNodes[kf0][a];
              localIndex const node1 = facesToNodes[kf1][ a==0 ? a : numNodesPerFace-a ];
              R1Tensor gap = u[node1];
              gap -= u[node0];
              real64 const gapNormal = Dot(gap,Nbar);

              if( gapNormal < 0 )
              {
                penaltyForce *= -contactStiffness * gapNormal * Ja;
                for( int i=0 ; i<3 ; ++i )
                {
                  rowDOF[3*a+i]                     = nodeDofNumber[node0]+i;
                  rowDOF[3*(numNodesPerFace + a)+i] = nodeDofNumber[node1]+i;


                  fc[node0] -= penaltyForce;
                  fc[node1] += penaltyForce;
                  nodeRHS[3*a+i]                     -= penaltyForce[i];
                  nodeRHS[3*(numNodesPerFace + a)+i] += penaltyForce[i];

                  dRdP(3*a+i,3*a+i)                                         -= contactStiffness * Ja * Nbar[i] * Nbar[i] ;
                  dRdP(3*a+i,3*(numNodesPerFace + a)+i)                     += contactStiffness * Ja * Nbar[i] * Nbar[i] ;
                  dRdP(3*(numNodesPerFace + a)+i,3*a+i)                     += contactStiffness * Ja * Nbar[i] * Nbar[i] ;
                  dRdP(3*(numNodesPerFace + a)+i,3*(numNodesPerFace + a)+i) -= contactStiffness * Ja * Nbar[i] * Nbar[i] ;
                }
              }
            }

            rhs->add( rowDOF,
                      nodeRHS,
                      numNodesPerFace*3*2 );

            matrix->add( rowDOF,
                         rowDOF,
                         dRdP.data(),
                         numNodesPerFace * 3 *2,
                         numNodesPerFace * 3 *2 );
          }
        });
    });
  }
}

real64
SolidMechanicsLagrangianFEM::ScalingForSystemSolution( DomainPartition const * const GEOSX_UNUSED_PARAM( domain ),
                                                       DofManager const & GEOSX_UNUSED_PARAM( dofManager ),
                                                       ParallelVector const & GEOSX_UNUSED_PARAM( solution ) )
{
  GEOSX_MARK_FUNCTION;
  real64 scalingFactor = 1.0;
//  MeshLevel const * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
//  FaceManager const * const faceManager = mesh->getFaceManager();
//  NodeManager const * const nodeManager = mesh->getNodeManager();
//  ElementRegionManager const * const elemManager = mesh->getElemManager();
//
//
//  arrayView1d<R1Tensor const> const & u = nodeManager->getReference< array1d<R1Tensor> >( keys::TotalDisplacement );
//
//  arrayView1d<R1Tensor const> const & faceNormal = faceManager->faceNormal();
//  array1d<localIndex_array> const & facesToNodes = faceManager->nodeList();
//
//  string const dofKey = dofManager.getKey( keys::TotalDisplacement );
//  arrayView1d<globalIndex> const & nodeDofNumber = nodeManager->getReference<globalIndex_array>( dofKey );
//
//  real64 const * soln = nullptr;
//  solution.extractLocalVector( &( const_cast<real64*&>(soln) ) );
//  globalIndex const rankOffset = dofManager.rankOffset();
//
//
//
//  elemManager->forElementSubRegions<FaceElementSubRegion>([&]( FaceElementSubRegion const * const subRegion )->void
//  {
//      arrayView1d<integer const> const & ghostRank = subRegion->GhostRank();
//      arrayView1d<real64 const > const & area = subRegion->getElementArea();
//      arrayView2d< localIndex const > const & elemsToFaces = subRegion->faceList();
//
//      RAJA::ReduceMin<RAJA::seq_reduce, real64> newScaleFactor(1.0);
//      forall_in_range<serialPolicy>( 0,
//                                     subRegion->size(),
//                                     GEOSX_LAMBDA ( localIndex const kfe )
//      {
//
//        if( ghostRank[kfe] < 0 )
//        {
//          R1Tensor Nbar = faceNormal[elemsToFaces[kfe][0]];
//          Nbar -= faceNormal[elemsToFaces[kfe][1]];
//          Nbar.Normalize();
//
//          localIndex const kf0 = elemsToFaces[kfe][0];
//          localIndex const kf1 = elemsToFaces[kfe][1];
//          localIndex const numNodesPerFace=facesToNodes[kf0].size();
//          localIndex const * const nodelist0 = facesToNodes[kf0];
//          localIndex const * const nodelist1 = facesToNodes[kf1];
//
//          for( localIndex a=0 ; a<numNodesPerFace ; ++a )
//          {
//            R1Tensor penaltyForce = Nbar;
//            localIndex const node0 = facesToNodes[kf0][a];
//            localIndex const node1 = facesToNodes[kf1][ a==0 ? a : numNodesPerFace-a ];
//
//            localIndex const lid0 = nodeDofNumber[node0] - rankOffset;
//            localIndex const lid1 = nodeDofNumber[node1] - rankOffset;
//            GEOSX_ASSERT( lid0 >= 0 && lid1 >=0 ); // since vectors are partitioned same as the mesh
//
//            R1Tensor gap = u[node1];
//            gap -= u[node0];
//            real64 const gapNormal = Dot(gap,Nbar);
//
//            R1Tensor deltaGap ;
//            for( int i=0 ; i<3 ; ++i )
//            {
//              deltaGap[i] = soln[lid1] - soln[lid0];
//            }
//            real64 const deltaGapNormal = Dot( deltaGap, Nbar );
//
//            if( ( gapNormal + deltaGapNormal ) * gapNormal < 0 )
//            {
//              newScaleFactor.min( - gapNormal / deltaGapNormal * 1.05 );
//              newScaleFactor.min( 1.0 );
//              std::cout<< "gapNormal, deltaGapNormal, scaleFactor, newGap = "<<
//                  gapNormal<<", "<<deltaGapNormal<<", "<<newScaleFactor.get()<<", "<<gapNormal+newScaleFactor.get()*deltaGapNormal<<std::endl;
//            }
//          }
//        }
//      });
//      scalingFactor = newScaleFactor.get();
//  });

  return scalingFactor;
}

void
SolidMechanicsLagrangianFEM::updateStress( DomainPartition * const GEOSX_UNUSED_PARAM(domain) )
{

}

REGISTER_CATALOG_ENTRY( SolverBase, SolidMechanicsLagrangianFEM, string const &, dataRepository::Group * const )
}
