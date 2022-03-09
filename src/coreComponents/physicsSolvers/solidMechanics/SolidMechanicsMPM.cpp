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

//  checkModelNames( m_solidMaterialNames, viewKeyStruct::solidMaterialNamesString() );

  LinearSolverParameters & linParams = m_linearSolverParameters.get();
  linParams.isSymmetric = true;
  linParams.dofsPerNode = 3;
  linParams.amg.separateComponents = true;
}

SolidMechanicsMPM::~SolidMechanicsMPM()
{
  // TODO Auto-generated destructor stub
}


void SolidMechanicsMPM::registerDataOnMesh( Group & meshBodies ) // Apparently I wasted time on this and it's not actually called by anything... SJP
{
  ExecutableGroup::registerDataOnMesh( meshBodies );

  forMeshTargets( meshBodies, [&] ( string const &,
                                    MeshLevel & meshLevel,
                                    arrayView1d< string const > const & regionNames )
  {
    ParticleManager & particleManager = meshLevel.getParticleManager();

    MeshBody const & meshBody = dynamicCast< MeshBody const & >( meshLevel.getParent().getParent() );

    // Set constitutive names on particles
    if(meshBody.m_hasParticles)
    {
      particleManager.forParticleSubRegions< ParticleSubRegionBase >( regionNames,
                                                                      [&]( localIndex const,
                                                                          ParticleSubRegionBase & subRegion )
      {
        setConstitutiveNamesCallSuper( subRegion );
        setConstitutiveNames( subRegion );
      } );
    }

  } );

  forMeshTargets( meshBodies, [&] ( string const &,
                                    MeshLevel & meshLevel,
                                    arrayView1d<string const> const & regionNames )
  {
    MeshBody const & meshBody = dynamicCast< MeshBody const & >( meshLevel.getParent().getParent() );

    if(!meshBody.m_hasParticles) // Background grid field registration
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
      nodeSets.registerWrapper<SortedArray<localIndex>>( viewKeyStruct::sendOrReceiveNodesString() ).
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
    }
    else // Particle field registration
    {
      std::cout << "Registering particle fields" << std::endl;
    }

  } );
}



void SolidMechanicsMPM::initializePreSubGroups()
{
  SolverBase::initializePreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );

  forMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                MeshLevel & meshLevel,
                                                arrayView1d<string const> const & regionNames )
  {
    MeshBody const & meshBody = dynamicCast< MeshBody const & >( meshLevel.getParent().getParent() );

    if(meshBody.m_hasParticles) // Only particle regions will hold actual materials. Background grid currently holds a null material so that the input file parser doesn't complain, but we don't need to actually do anything with it.
    {
      ParticleManager & particleManager = meshLevel.getParticleManager();
      particleManager.forParticleSubRegions< ParticleSubRegion >( regionNames, [&]( localIndex const,
                                                                                    ParticleSubRegion & subRegion )
      {
        string & solidMaterialName = subRegion.getReference<string>( viewKeyStruct::solidMaterialNamesString() );
        solidMaterialName = SolverBase::getConstitutiveName<SolidBase>( subRegion );
      });
    }
  });

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
                                                  real64 const dt,
                                                  std::string const & elementListName )
{
  GEOSX_MARK_FUNCTION;
  real64 rval = 0;
  if( m_strainTheory==0 )
  {
    auto kernelFactory = SolidMechanicsLagrangianFEMKernels::ExplicitSmallStrainFactory( dt, elementListName );
    rval = finiteElement::
             regionBasedKernelApplication< parallelDevicePolicy< 32 >,
                                           constitutive::SolidBase,
                                           CellElementSubRegion >( mesh,
                                                                   targetRegions,
                                                                   finiteElementName,
                                                                   viewKeyStruct::solidMaterialNamesString(),
                                                                   kernelFactory );
  }
  else if( m_strainTheory==1 )
  {
    auto kernelFactory = SolidMechanicsLagrangianFEMKernels::ExplicitFiniteStrainFactory( dt, elementListName );
    rval = finiteElement::
             regionBasedKernelApplication< parallelDevicePolicy< 32 >,
                                           constitutive::SolidBase,
                                           CellElementSubRegion >( mesh,
                                                                   targetRegions,
                                                                   finiteElementName,
                                                                   viewKeyStruct::solidMaterialNamesString(),
                                                                   kernelFactory );
  }
  else
  {
    GEOSX_ERROR( "Invalid option for strain theory (0 = infinitesimal strain, 1 = finite strain" );
  }

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
    GEOSX_ERROR( "MPM solver only currently supports explicit time stepping!" );
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

  // Get node and particle managers. ***** We implicitly assume that there are exactly two mesh bodies, and that one has particles and one does not. *****
  Group & meshBodies = domain.getMeshBodies();
  MeshBody & meshBody1 = meshBodies.getGroup< MeshBody >(0);
  MeshBody & meshBody2 = meshBodies.getGroup< MeshBody >(1);
  ParticleManager & particleManager = meshBody1.m_hasParticles ? meshBody1.getMeshLevel(0).getParticleManager() : meshBody2.getMeshLevel(0).getParticleManager();
  NodeManager & nodeManager = !meshBody1.m_hasParticles ? meshBody1.getMeshLevel(0).getNodeManager() : meshBody2.getMeshLevel(0).getNodeManager();

  // Loop over nodes
  arrayView2d< real64, nodes::REFERENCE_POSITION_USD > & X = nodeManager.referencePosition();
  std::vector<double> xMin{DBL_MAX,DBL_MAX,DBL_MAX}, xMax{DBL_MIN,DBL_MIN,DBL_MIN}, domainL{0.0,0.0,0.0};

  for(int i=0; i<X.size()/3; i++)
  {
    for(int j=0; j<3; j++)
    {
      xMin[j] = std::fmin(xMin[j],X[i][j]);
      xMax[j] = std::fmax(xMin[j],X[i][j]);
      domainL[j] = xMax[j] - xMin[j];
    }
  }

  // Loop over particles
  particleManager.forParticleRegions( [&]( auto & region )
  {
    region.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
    {
      arrayView2d< real64 > const particleCenter = subRegion.getParticleCenter();
      for(int i=0; i<subRegion.size(); i++)
      {
        particleCenter[i][0] += 0.15708*cos(0.314159*(time_n + 0.5));
//        intArray nodes = getNodes(particle);
//        realArray weights = getWeights(particle);
//        for(nodes : node)
//        {
//          momentum(node) += particle contribution;
//        }
      }
      // forall with policy - use atomics to avoid race conditions
      // use array2d's 2nd template argument like node.referencePosition()
    } );
  } );

  return dt;
}

void SolidMechanicsMPM::setConstitutiveNamesCallSuper( ParticleSubRegionBase & subRegion ) const
{
  SolverBase::setConstitutiveNamesCallSuper( subRegion );

  subRegion.registerWrapper< string >( viewKeyStruct::solidMaterialNamesString() ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setSizedFromParent( 0 );

  string & solidMaterialName = subRegion.getReference< string >( viewKeyStruct::solidMaterialNamesString() );
  solidMaterialName = SolverBase::getConstitutiveName< SolidBase >( subRegion );
  GEOSX_ERROR_IF( solidMaterialName.empty(), GEOSX_FMT( "SolidBase model not found on subregion {}", subRegion.getName() ) );
}

void SolidMechanicsMPM::setConstitutiveNames( ParticleSubRegionBase & subRegion ) const
{
  GEOSX_UNUSED_VAR( subRegion );
}


REGISTER_CATALOG_ENTRY( SolverBase, SolidMechanicsMPM, string const &, dataRepository::Group * const )
}
