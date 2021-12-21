/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file SolidMechanicsMPM.cpp
 */

#include "SolidMechanicsMPM.hpp"
#include "SolidMechanicsSmallStrainQuasiStaticKernel.hpp"
#include "SolidMechanicsSmallStrainImplicitNewmarkKernel.hpp"
#include "SolidMechanicsSmallStrainExplicitNewmarkKernel.hpp"
#include "SolidMechanicsFiniteStrainExplicitNewmarkKernel.hpp"

#include "codingUtilities/Utilities.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/contact/ContactBase.hpp"
#include "finiteElement/FiniteElementDiscretizationManager.hpp"
#include "finiteElement/Kinematics.h"
#include "LvArray/src/output.hpp"
#include "mesh/DomainPartition.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "discretizationMethods/NumericalMethodsManager.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "fieldSpecification/TractionBoundaryCondition.hpp"
#include "mesh/FaceElementSubRegion.hpp"
#include "mesh/utilities/ComputationalGeometry.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "mesh/mpiCommunications/NeighborCommunicator.hpp"
#include "common/GEOS_RAJA_Interface.hpp"


namespace geosx
{

using namespace dataRepository;
using namespace constitutive;

SolidMechanicsMPM::SolidMechanicsMPM( const string & name,
                                                          Group * const parent ):
  SolverBase( name, parent ),
  m_newmarkGamma( 0.5 ),
  m_newmarkBeta( 0.25 ),
  m_massDamping( 0.0 ),
  m_stiffnessDamping( 0.0 ),
  m_timeIntegrationOption( TimeIntegrationOption::ExplicitDynamic ),
  m_useVelocityEstimateForQS( 0 ),
  m_maxForce( 0.0 ),
  m_maxNumResolves( 10 ),
  m_strainTheory( 0 ),
//  m_elemsAttachedToSendOrReceiveNodes(),
//  m_elemsNotAttachedToSendOrReceiveNodes(),
//  m_sendOrReceiveNodes(),
//  m_nonSendOrReceiveNodes(),
//  m_targetNodes(),
  m_iComm( CommunicationTools::getInstance().getCommID() )
{
//  m_sendOrReceiveNodes.setName( "SolidMechanicsMPM::m_sendOrReceiveNodes" );
//  m_nonSendOrReceiveNodes.setName( "SolidMechanicsMPM::m_nonSendOrReceiveNodes" );
//  m_targetNodes.setName( "SolidMechanicsMPM::m_targetNodes" );

  registerWrapper( viewKeyStruct::newmarkGammaString(), &m_newmarkGamma ).
    setApplyDefaultValue( 0.5 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Value of :math:`\\gamma` in the Newmark Method for Implicit Dynamic time integration option" );

  registerWrapper( viewKeyStruct::newmarkBetaString(), &m_newmarkBeta ).
    setApplyDefaultValue( 0.25 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Value of :math:`\\beta` in the Newmark Method for Implicit Dynamic time integration option. "
                    "This should be pow(newmarkGamma+0.5,2.0)/4.0 unless you know what you are doing." );

  registerWrapper( viewKeyStruct::massDampingString(), &m_massDamping ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Value of mass based damping coefficient. " );

  registerWrapper( viewKeyStruct::stiffnessDampingString(), &m_stiffnessDamping ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Value of stiffness based damping coefficient. " );

  registerWrapper( viewKeyStruct::timeIntegrationOptionString(), &m_timeIntegrationOption ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_timeIntegrationOption ).
    setDescription( "Time integration method. Options are:\n* " + EnumStrings< TimeIntegrationOption >::concat( "\n* " ) );

  registerWrapper( viewKeyStruct::useVelocityEstimateForQSString(), &m_useVelocityEstimateForQS ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flag to indicate the use of the incremental displacement from the previous step as an "
                    "initial estimate for the incremental displacement of the current step." );

  registerWrapper( viewKeyStruct::maxNumResolvesString(), &m_maxNumResolves ).
    setApplyDefaultValue( 10 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Value to indicate how many resolves may be executed after some other event is executed. "
                    "For example, if a SurfaceGenerator is specified, it will be executed after the mechanics solve. "
                    "However if a new surface is generated, then the mechanics solve must be executed again due to the "
                    "change in topology." );

  registerWrapper( viewKeyStruct::strainTheoryString(), &m_strainTheory ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Indicates whether or not to use "
                    "`Infinitesimal Strain Theory <https://en.wikipedia.org/wiki/Infinitesimal_strain_theory>`_, or "
                    "`Finite Strain Theory <https://en.wikipedia.org/wiki/Finite_strain_theory>`_. Valid Inputs are:\n"
                    " 0 - Infinitesimal Strain \n"
                    " 1 - Finite Strain" );

  registerWrapper( viewKeyStruct::solidMaterialNamesString(), &m_solidMaterialNames ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "The name of the material that should be used in the constitutive updates" );

  registerWrapper( viewKeyStruct::contactRelationNameString(), &m_contactRelationName ).
    setApplyDefaultValue( viewKeyStruct::noContactRelationNameString() ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Name of contact relation to enforce constraints on fracture boundary." );

  registerWrapper( viewKeyStruct::maxForceString(), &m_maxForce ).
    setInputFlag( InputFlags::FALSE ).
    setDescription( "The maximum force contribution in the problem domain." );

}

void SolidMechanicsMPM::postProcessInput()
{
  SolverBase::postProcessInput();

  checkModelNames( m_solidMaterialNames, viewKeyStruct::solidMaterialNamesString() );

  LinearSolverParameters & linParams = m_linearSolverParameters.get();
  linParams.isSymmetric = true;
  linParams.dofsPerNode = 3;
  linParams.amg.separateComponents = true;
}

SolidMechanicsMPM::~SolidMechanicsMPM()
{
  // TODO Auto-generated destructor stub
}


void SolidMechanicsMPM::registerDataOnMesh( Group & meshBodies )
{
  forMeshTargets( meshBodies, [&] ( MeshLevel & meshLevel,
                                    arrayView1d<string const> const & regionNames )
  {
    NodeManager & nodes = meshLevel.getNodeManager();

    nodes.registerWrapper< array2d< real64, nodes::TOTAL_DISPLACEMENT_PERM > >( keys::TotalDisplacement ).
      setPlotLevel( PlotLevel::LEVEL_0 ).
      setRegisteringObjects( this->getName()).
      setDescription( "An array that holds the total displacements on the nodes." ).
      reference().resizeDimension< 1 >( 3 );

    nodes.registerWrapper< array2d< real64, nodes::INCR_DISPLACEMENT_PERM > >( keys::IncrementalDisplacement ).
      setPlotLevel( PlotLevel::LEVEL_3 ).
      setRegisteringObjects( this->getName()).
      setDescription( "An array that holds the incremental displacements for the current time step on the nodes." ).
      reference().resizeDimension< 1 >( 3 );

    nodes.registerWrapper< array2d< real64, nodes::VELOCITY_PERM > >( keys::Velocity ).
      setPlotLevel( PlotLevel::LEVEL_0 ).
      setRegisteringObjects( this->getName()).
      setDescription( "An array that holds the current velocity on the nodes." ).
      reference().resizeDimension< 1 >( 3 );

    nodes.registerWrapper< array2d< real64, nodes::ACCELERATION_PERM > >( keys::Acceleration ).
      setPlotLevel( PlotLevel::LEVEL_1 ).
      setRegisteringObjects( this->getName()).
      setDescription( "An array that holds the current acceleration on the nodes. This array also is used "
                      "to hold the summation of nodal forces resulting from the governing equations." ).
      reference().resizeDimension< 1 >( 3 );

    nodes.registerWrapper< array2d< real64 > >( viewKeyStruct::forceExternalString() ).
      setPlotLevel( PlotLevel::LEVEL_0 ).
      setRegisteringObjects( this->getName()).
      setDescription( "An array that holds the external forces on the nodes. This includes any boundary"
                      " conditions as well as coupling forces such as hydraulic forces." ).
      reference().resizeDimension< 1 >( 3 );

    nodes.registerWrapper< array1d< real64 > >( keys::Mass ).
      setPlotLevel( PlotLevel::LEVEL_0 ).
      setRegisteringObjects( this->getName()).
      setDescription( "An array that holds the mass on the nodes." );

    nodes.registerWrapper< array2d< real64 > >( viewKeyStruct::vTildeString() ).
      setPlotLevel( PlotLevel::NOPLOT ).
      setRegisteringObjects( this->getName()).
      setDescription( "An array that holds the velocity predictors on the nodes." ).
      reference().resizeDimension< 1 >( 3 );

    nodes.registerWrapper< array2d< real64 > >( viewKeyStruct::uhatTildeString() ).
      setPlotLevel( PlotLevel::NOPLOT ).
      setRegisteringObjects( this->getName()).
      setDescription( "An array that holds the incremental displacement predictors on the nodes." ).
      reference().resizeDimension< 1 >( 3 );

    nodes.registerWrapper< array2d< real64 > >( viewKeyStruct::contactForceString() ).
      setPlotLevel( PlotLevel::LEVEL_0 ).
      setRegisteringObjects( this->getName()).
      setDescription( "An array that holds the contact force." ).
      reference().resizeDimension< 1 >( 3 );

    Group & nodeSets = nodes.sets();
    nodeSets.registerWrapper<SortedArray<localIndex>>( viewKeyStruct::sendOrRecieveNodesString() ).
      setPlotLevel( PlotLevel::NOPLOT ).
      setRestartFlags( RestartFlags::NO_WRITE );

    nodeSets.registerWrapper<SortedArray<localIndex>>( viewKeyStruct::nonSendOrReceiveNodesString() ).
      setPlotLevel( PlotLevel::NOPLOT ).
      setRestartFlags( RestartFlags::NO_WRITE );

    nodeSets.registerWrapper<SortedArray<localIndex>>( viewKeyStruct::targetNodesString() ).
      setPlotLevel( PlotLevel::NOPLOT ).
      setRestartFlags( RestartFlags::NO_WRITE );

    ElementRegionManager & elementRegionManager = meshLevel.getElemManager();
    elementRegionManager.forElementSubRegions< CellElementSubRegion >( regionNames,
                                                                       [&]( localIndex const,
                                                                            CellElementSubRegion & subRegion )
    {
      subRegion.registerWrapper< SortedArray< localIndex > >( viewKeyStruct::elemsAttachedToSendOrReceiveNodesString() ).
        setPlotLevel( PlotLevel::NOPLOT ).
        setRestartFlags( RestartFlags::NO_WRITE );

      subRegion.registerWrapper< SortedArray< localIndex > >( viewKeyStruct::elemsNotAttachedToSendOrReceiveNodesString() ).
        setPlotLevel( PlotLevel::NOPLOT ).
        setRestartFlags( RestartFlags::NO_WRITE );
    } );

  } );
}



void SolidMechanicsMPM::initializePreSubGroups()
{
  SolverBase::initializePreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );

  GEOSX_ERROR_IF(domain.getMeshBodies().numSubGroups() != 2, "The GEOSX MPM solver requires two mesh bodies corresponding to the background grid and particles!");

  // Validate solid models in target regions
  domain.forMeshBodies( [&]( MeshBody const & meshBody )
  {
    MeshLevel const & meshLevel = meshBody.getMeshLevel( 0 );
    validateModelMapping< SolidBase >( meshLevel.getElemManager(), m_solidMaterialNames );
  } );

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();

  FiniteElementDiscretizationManager const &
  feDiscretizationManager = numericalMethodManager.getFiniteElementDiscretizationManager();

  FiniteElementDiscretization const &
  feDiscretization = feDiscretizationManager.getGroup< FiniteElementDiscretization >( m_discretizationName );
  GEOSX_UNUSED_VAR( feDiscretization );
}



template< typename ... PARAMS >
real64 SolidMechanicsMPM::explicitKernelDispatch( MeshLevel & mesh,
                                                            arrayView1d< string const > const & targetRegions,
                                                            string const & finiteElementName,
                                                            arrayView1d< string const > const & constitutiveNames,
                                                            real64 const dt,
                                                            std::string const & elementListName )
{
  GEOSX_MARK_FUNCTION;
  real64 rval = 0;
  auto kernelFactory = SolidMechanicsLagrangianFEMKernels::ExplicitFiniteStrainFactory( dt, elementListName );
  rval = finiteElement::
           regionBasedKernelApplication< parallelDevicePolicy< 32 >,
                                         constitutive::SolidBase,
                                         CellElementSubRegion >( mesh,
                                                                 targetRegions,
                                                                 finiteElementName,
                                                                 constitutiveNames,
                                                                 kernelFactory );
  return rval;
}



real64 SolidMechanicsMPM::solverStep( real64 const & time_n,
                                                real64 const & dt,
                                                const int cycleNumber,
                                                DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;
  real64 dtReturn = dt;

  SolverBase * const surfaceGenerator = this->getParent().getGroupPointer< SolverBase >( "SurfaceGen" );

  if( m_timeIntegrationOption == TimeIntegrationOption::ExplicitDynamic )
  {
    dtReturn = explicitStep( time_n, dt, cycleNumber, domain );

    if( surfaceGenerator!=nullptr )
    {
      surfaceGenerator->solverStep( time_n, dt, cycleNumber, domain );
    }
  }
  else
  {
    GEOSX_LOG_RANK_0( "MPM solver only currently supports explicit time stepping!" );
  }

  return dtReturn;
}

real64 SolidMechanicsMPM::explicitStep( real64 const & time_n,
                                                  real64 const & dt,
                                                  const int GEOSX_UNUSED_PARAM( cycleNumber ),
                                                  DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  #define USE_PHYSICS_LOOP

  Group & meshBodies = domain.getMeshBodies();





//  This is the explicit step function for the Lagrangian FEM solver which I'm keeping here for reference.
//  for( localIndex a = 0; a < meshBodies.numSubGroups(); ++a )
//  {
//    MeshLevel & mesh = domain.getMeshBody( a ).getMeshLevel( 0 );
//    NodeManager & nodes = mesh.getNodeManager();
//    Group const & nodeSets = nodes.sets();
//
//    SortedArrayView< localIndex const > const &
//    m_sendOrReceiveNodes = nodeSets.getReference<SortedArray<localIndex>>( viewKeyStruct::sendOrRecieveNodesString() ).toViewConst();
//
//    SortedArrayView< localIndex const > const &
//    m_nonSendOrReceiveNodes = nodeSets.getReference<SortedArray<localIndex>>( viewKeyStruct::nonSendOrReceiveNodesString() ).toViewConst();
//
//    SortedArrayView< localIndex const> const &
//    m_targetNodes = nodeSets.getReference<SortedArray<localIndex>>( viewKeyStruct::targetNodesString() ).toViewConst();
//
//
//    // save previous constitutive state data in preparation for next timestep
//    forTargetSubRegions< CellElementSubRegion >( mesh, [&]( localIndex const targetIndex,
//                                                            CellElementSubRegion & subRegion )
//    {
//      SolidBase & constitutiveRelation = getConstitutiveModel< SolidBase >( subRegion, m_solidMaterialNames[targetIndex] );
//      constitutiveRelation.saveConvergedState();
//    } );
//
//
//
//    // Print out nodal reference positions?
//    arrayView2d< real64, nodes::REFERENCE_POSITION_USD > const & x0 = nodes.referencePosition();
//    if(time_n==0.0)
//    {
//      std::cout << "Printing... List size is: " << x0.size()/3 << std::endl; // The list type comes from a class that probably does fancy things when you print it naively like I'm doing here... it's stored as a 1D array in memory, hence having to divide by 3.
//      for(localIndex i=0; i<x0.size()/3; i++)
//      {
//        std::cout << x0[i] << std::endl;
//      }
//    }
//
//
//
//
//    FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();
//
//    arrayView1d< real64 const > const & mass = nodes.getReference< array1d< real64 > >( keys::Mass );
//    arrayView2d< real64, nodes::VELOCITY_USD > const & vel = nodes.velocity();
//
//    arrayView2d< real64, nodes::TOTAL_DISPLACEMENT_USD > const & u = nodes.totalDisplacement();
//    arrayView2d< real64, nodes::INCR_DISPLACEMENT_USD > const & uhat = nodes.incrementalDisplacement();
//    arrayView2d< real64, nodes::ACCELERATION_USD > const & acc = nodes.acceleration();
//
//    std::map< string, string_array > fieldNames;
//    fieldNames["node"].emplace_back( keys::Velocity );
//    fieldNames["node"].emplace_back( keys::Acceleration );
//
//    m_iComm.resize( domain.getNeighbors().size() );
//    CommunicationTools::getInstance().synchronizePackSendRecvSizes( fieldNames, mesh, domain.getNeighbors(), m_iComm, true );
//
//    fsManager.applyFieldValue< parallelDevicePolicy< 1024 > >( time_n, domain, "nodeManager", keys::Acceleration );
//
//    //3: v^{n+1/2} = v^{n} + a^{n} dt/2
//    SolidMechanicsLagrangianFEMKernels::velocityUpdate( acc, vel, dt/2 );
//
//    fsManager.applyFieldValue< parallelDevicePolicy< 1024 > >( time_n, domain, "nodeManager", keys::Velocity );
//
//    //4. x^{n+1} = x^{n} + v^{n+{1}/{2}} dt (x is displacement)
//    SolidMechanicsLagrangianFEMKernels::displacementUpdate( vel, uhat, u, dt );
//
//    fsManager.applyFieldValue( time_n + dt,
//                               domain, "nodeManager",
//                               NodeManager::viewKeyStruct::totalDisplacementString(),
//                               [&]( FieldSpecificationBase const & bc,
//                                    SortedArrayView< localIndex const > const & targetSet )
//    {
//      integer const component = bc.getComponent();
//      GEOSX_ERROR_IF_LT_MSG( component, 0, "Component index required for displacement BC " << bc.getName() );
//
//      forAll< parallelDevicePolicy< 1024 > >( targetSet.size(),
//                                              [=] GEOSX_DEVICE ( localIndex const i )
//      {
//        localIndex const a = targetSet[ i ];
//        vel( a, component ) = u( a, component );
//      } );
//    },
//                               [&]( FieldSpecificationBase const & bc,
//                                    SortedArrayView< localIndex const > const & targetSet )
//    {
//      integer const component = bc.getComponent();
//      GEOSX_ERROR_IF_LT_MSG( component, 0, "Component index required for displacement BC " << bc.getName() );
//
//      forAll< parallelDevicePolicy< 1024 > >( targetSet.size(),
//                                              [=] GEOSX_DEVICE ( localIndex const i )
//      {
//        localIndex const a = targetSet[ i ];
//        uhat( a, component ) = u( a, component ) - vel( a, component );
//        vel( a, component )  = uhat( a, component ) / dt;
//      } );
//    } );
//
//    //Step 5. Calculate deformation input to constitutive model and update state to
//    // Q^{n+1}
//    explicitKernelDispatch( mesh,
//                            targetRegionNames(),
//                            this->getDiscretizationName(),
//                            m_solidMaterialNames,
//                            dt,
//                            string( viewKeyStruct::elemsAttachedToSendOrReceiveNodesString() ) );
//
//    // apply this over a set
//    SolidMechanicsLagrangianFEMKernels::velocityUpdate( acc, mass, vel, dt / 2, m_sendOrReceiveNodes.toViewConst() );
//
//    fsManager.applyFieldValue< parallelDevicePolicy< 1024 > >( time_n, domain, "nodeManager", keys::Velocity );
//
//    parallelDeviceEvents packEvents;
//    CommunicationTools::getInstance().asyncPack( fieldNames, mesh, domain.getNeighbors(), m_iComm, true, packEvents );
//
//    waitAllDeviceEvents( packEvents );
//
//    CommunicationTools::getInstance().asyncSendRecv( domain.getNeighbors(), m_iComm, true, packEvents );
//
//    explicitKernelDispatch( mesh,
//                            targetRegionNames(),
//                            this->getDiscretizationName(),
//                            m_solidMaterialNames,
//                            dt,
//                            string( viewKeyStruct::elemsNotAttachedToSendOrReceiveNodesString() ) );
//
//    // apply this over a set
//    SolidMechanicsLagrangianFEMKernels::velocityUpdate( acc, mass, vel, dt / 2, m_nonSendOrReceiveNodes.toViewConst() );
//    fsManager.applyFieldValue< parallelDevicePolicy< 1024 > >( time_n, domain, "nodeManager", keys::Velocity );
//
//    // this includes  a device sync after launching all the unpacking kernels
//    parallelDeviceEvents unpackEvents;
//    CommunicationTools::getInstance().finalizeUnpack( mesh, domain.getNeighbors(), m_iComm, true, unpackEvents );
//  }



  return dt;
}



REGISTER_CATALOG_ENTRY( SolverBase, SolidMechanicsMPM, string const &, dataRepository::Group * const )
}
