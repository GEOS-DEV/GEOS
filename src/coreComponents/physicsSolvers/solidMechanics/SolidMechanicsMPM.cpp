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
#include "mesh/mpiCommunications/SpatialPartition.hpp"
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
    else // Particle field registration? TODO: What goes here?
    {
      GEOSX_LOG_RANK_0("Registering particle fields, I guess");
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



//template< typename ... PARAMS >
//real64 SolidMechanicsMPM::explicitKernelDispatch( MeshLevel & mesh,
//                                                  arrayView1d< string const > const & targetRegions,
//                                                  string const & finiteElementName,
//                                                  real64 const dt,
//                                                  std::string const & elementListName )
//{
//  GEOSX_MARK_FUNCTION;
//  real64 rval = 0;
//  if( m_strainTheory==0 )
//  {
//    auto kernelFactory = SolidMechanicsLagrangianFEMKernels::ExplicitSmallStrainFactory( dt, elementListName );
//    rval = finiteElement::
//             regionBasedKernelApplication< parallelDevicePolicy< 32 >,
//                                           constitutive::SolidBase,
//                                           CellElementSubRegion >( mesh,
//                                                                   targetRegions,
//                                                                   finiteElementName,
//                                                                   viewKeyStruct::solidMaterialNamesString(),
//                                                                   kernelFactory );
//  }
//  else if( m_strainTheory==1 )
//  {
//    auto kernelFactory = SolidMechanicsLagrangianFEMKernels::ExplicitFiniteStrainFactory( dt, elementListName );
//    rval = finiteElement::
//             regionBasedKernelApplication< parallelDevicePolicy< 32 >,
//                                           constitutive::SolidBase,
//                                           CellElementSubRegion >( mesh,
//                                                                   targetRegions,
//                                                                   finiteElementName,
//                                                                   viewKeyStruct::solidMaterialNamesString(),
//                                                                   kernelFactory );
//  }
//  else
//  {
//    GEOSX_ERROR( "Invalid option for strain theory (0 = infinitesimal strain, 1 = finite strain" );
//  }
//
//  return rval;
//}

bool SolidMechanicsMPM::execute( real64 const time_n,
                                 real64 const dt,
                                 integer const cycleNumber,
                                 integer const GEOSX_UNUSED_PARAM( eventCounter ),
                                 real64 const GEOSX_UNUSED_PARAM( eventProgress ),
                                 DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  m_nextDt = solverStep( time_n,
                         dt,
                         cycleNumber,
                         domain );

  return false;
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

void SolidMechanicsMPM::initialize(arrayView2d< real64, nodes::REFERENCE_POSITION_USD > const & g_X, ParticleManager & particleManager, SpatialPartition & partition)
{
  // Get global domain extent
  for(int i=0; i<g_X.size()/3; i++)
  {
    for(int j=0; j<3; j++)
    {
      m_xLocalMin[j] = std::min(m_xLocalMin[j],g_X[i][j]);
      m_xLocalMax[j] = std::max(m_xLocalMax[j],g_X[i][j]);
    }
  }
  for(int i=0; i<3; i++)
  {
    m_xLocalMinNoGhost[i] = partition.getLocalMin()[i];
    m_xLocalMaxNoGhost[i] = partition.getLocalMax()[i];
    m_xGlobalMin[i] = partition.getGlobalMin()[i];
    m_xGlobalMax[i] = partition.getGlobalMax()[i];
    m_domainL[i] = m_xLocalMax[i] - m_xLocalMin[i];
  }

  // Get element size
  for(int i=0; i<g_X.size()/3; i++)
  {
    for(int j=0; j<3; j++)
    {
      real64 test = g_X[i][j] - m_xLocalMin[j]; // By definition, this should always be positive. Furthermore, the g_X should only be those on the local partition
      if(test > 0.0) // We're looking for the smallest nonzero distance from the "min" node. TODO: Could be vulnerable to a finite precision bug.
      {
        m_hx[j] = std::fmin(test, m_hx[j]);
      }
    }
  }

  // Get number of elements in each direction
  for(int i=0; i<3; i++)
  {
    m_nEl[i] = std::round(m_domainL[i]/m_hx[i]);
  }

  // Create element map
  m_ijkMap.resize( m_nEl[0] + 1 );
  for( int i = 0 ; i <= m_nEl[0] ; i++ )
  {
    m_ijkMap[i].resize( m_nEl[1] + 1 );
    for( int j = 0 ; j <= m_nEl[1] ; j++ )
    {
      m_ijkMap[i][j].resize( m_nEl[2] + 1 );
    }
  }
  for( int ii = 0 ; ii < g_X.size()/3 ; ii++ )
  {
    int i = std::round( ( g_X[ii][0] - m_xLocalMin[0] ) / m_hx[0] ) ;
    int j = std::round( ( g_X[ii][1] - m_xLocalMin[1] ) / m_hx[1] ) ;
    int k = std::round( ( g_X[ii][2] - m_xLocalMin[2] ) / m_hx[2] ) ;
    m_ijkMap[i][j][k] = ii;
  }

  // Set particle masses based on their volume and density. Set stress to zero. Set deformation gradient to identity;
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    string const & solidMaterialName = subRegion.template getReference< string >( viewKeyStruct::solidMaterialNamesString() );
    SolidBase & constitutiveRelation = getConstitutiveModel< SolidBase >( subRegion, solidMaterialName ); // For the time being we restrict our attention to elastic isotropic solids. TODO: Have all constitutive models automatically calculate a wave speed.
    arrayView2d< real64 > const particleDensity = constitutiveRelation.getDensity(); // 2d array because there's a density for each quadrature point, we just access with [particle][0]
    //arrayView3d< real64, solid::STRESS_USD > const particleStress = constitutiveRelation.getStress();
    arrayView1d< real64 > const particleVolume = subRegion.getParticleVolume();
    arrayView1d< real64 > const particleMass = subRegion.getParticleMass();
    arrayView3d< real64 > const particleDeformationGradient = subRegion.getParticleDeformationGradient();

    // mass
    for(int i=0; i<subRegion.size(); i++)
    {
      particleMass[i] = particleDensity[i][0]*particleVolume[i]; // TODO: This should probably be done in ParticleMeshGenerator...
    }

    // stress
    //particleStress.zero(); // This is actually handled by SolidBase.cpp - stress gets initialized to zero by default in the constructor

    // deformation gradient - TODO: there's probably a LvArray function that makes this a one-liner - I don't think the ParticleSubRegionBase constructor can easily initialize this to identity
    for(int p=0; p<subRegion.size(); p++)
    {
      for(int i=0; i<3; i++)
      {
        for(int j=0; j<3; j++)
        {
          if(i==j)
          {
            particleDeformationGradient[p][i][j] = 1.0;
          }
          else
          {
            particleDeformationGradient[p][i][j] = 0.0;
          }
        }
      }
    }

  } );

}

real64 SolidMechanicsMPM::explicitStep( real64 const & time_n,
                                        real64 const & dt,
                                        const int cycleNumber,
                                        DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  #define USE_PHYSICS_LOOP


  // Constitutive manager
  //ConstitutiveManager & constitutiveManager = domain.getConstitutiveManager();

  // Spatial Partition
  SpatialPartition & partition = dynamic_cast< SpatialPartition & >(domain.getReference< PartitionBase >( keys::partitionManager ) );

  // Get node and particle managers. ***** We implicitly assume that there are exactly two mesh bodies, and that one has particles and one does not. *****
  Group & meshBodies = domain.getMeshBodies();

  MeshBody & meshBody1 = meshBodies.getGroup< MeshBody >(0);
  MeshBody & meshBody2 = meshBodies.getGroup< MeshBody >(1);
  MeshBody & particles = meshBody1.m_hasParticles ? meshBody1 : meshBody2;
  MeshBody & grid = !meshBody1.m_hasParticles ? meshBody1 : meshBody2;

  ParticleManager & particleManager = particles.getMeshLevel(0).getParticleManager();
  MeshLevel & mesh = grid.getMeshLevel(0);
  NodeManager & nodeManager = mesh.getNodeManager();


  // Get nodal fields
  arrayView2d< real64, nodes::REFERENCE_POSITION_USD > & g_X = nodeManager.referencePosition();
  arrayView2d< real64, nodes::VELOCITY_USD > & g_V = nodeManager.velocity(); // Velocity, initially overloaded to be momentum
  arrayView2d< real64, nodes::ACCELERATION_USD > & g_A = nodeManager.acceleration(); // Acceleration, initially overloaded to be internal force
  arrayView1d< real64 > & g_M = nodeManager.getReference< array1d< real64 > >( keys::Mass );


  // Zero out nodal fields
  g_V.zero();
  g_A.zero();
  g_M.zero();


  // At time step zero, perform initialization calculations
  if(cycleNumber == 0)
  {
    initialize(g_X, particleManager, partition);
  }


  // Particle to grid interpolation
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    arrayView2d< real64 > const particleCenter = subRegion.getParticleCenter();
    arrayView2d< real64 > const particleVelocity = subRegion.getParticleVelocity();
    arrayView1d< real64 > const particleMass = subRegion.getParticleMass();
    arrayView1d< real64 > const particleVolume = subRegion.getParticleVolume();

    string const & solidMaterialName = subRegion.template getReference< string >( viewKeyStruct::solidMaterialNamesString() );
    SolidBase & constitutiveRelation = getConstitutiveModel< SolidBase >( subRegion, solidMaterialName );
    arrayView3d< real64, solid::STRESS_USD > const particleStress = constitutiveRelation.getStress();
    if(cycleNumber ==0)
    {
      particleStress.zero(); // Zero out particle stress on first time step
    }

    for(int p=0; p<subRegion.size(); p++)
    {
      auto const & p_x = particleCenter[p]; // auto = LvArray::ArraySlice<double, 1, 0, long>
      auto const & p_v = particleVelocity[p]; // auto = LvArray::ArraySlice<double, 1, 0, long>
      real64 const & p_m = particleMass[p];
      real64 const & p_Vol = particleVolume[p];
      auto const & p_stress = particleStress[p][0];

      // Get interpolation kernel
      std::vector<int> nodeIDs; // nodes that the particle maps to
      std::vector<real64> weights; // shape function value for each node
      std::vector< std::vector<real64> > gradWeights; // shape function gradient value for each node; 1st index = direction, 2nd index = node
      gradWeights.resize(3);
      subRegion.getAllWeights(p,
                              p_x,
                              m_xLocalMin,
                              m_hx,
                              m_ijkMap,
                              g_X,
                              nodeIDs,      // output
                              weights,      // output
                              gradWeights); // output

      // Update grid values
      for(size_t i=0; i<nodeIDs.size(); i++)
      {
        int g = nodeIDs[i];
        g_M[g] += p_m*weights[i];
        for(int j=0; j<3; j++)
        {
          g_V[g][j] += p_m*p_v[j]*weights[i];
          for(int k=0; k<3; k++)
          {
            int voigt = m_voigtMap[k][j];
            g_A[g][j] -= p_stress[voigt]*gradWeights[k][i]*p_Vol;
          }
        }
      }

    } // particle loop
  } ); // subregion loop


  // Grid MPI operations

  // (1) Initialize
  std::map< string, string_array > fieldNames;
  fieldNames["node"].emplace_back( keys::Mass );
  fieldNames["node"].emplace_back( keys::Velocity );
  fieldNames["node"].emplace_back( keys::Acceleration );
  std::vector< NeighborCommunicator > & neighbors = domain.getNeighbors();
  m_iComm.resize( neighbors.size() );

  // (2) Additive sync
  CommunicationTools::getInstance().synchronizePackSendRecvSizes( fieldNames, mesh, neighbors, m_iComm, true );
  parallelDeviceEvents packEvents;
  CommunicationTools::getInstance().asyncPack( fieldNames, mesh, neighbors, m_iComm, true, packEvents );
  waitAllDeviceEvents( packEvents );
  CommunicationTools::getInstance().asyncSendRecv( neighbors, m_iComm, true, packEvents );
  parallelDeviceEvents unpackEvents;
  CommunicationTools::getInstance().finalizeUnpack( mesh, neighbors, m_iComm, true, unpackEvents, MPI_SUM ); // needs an extra argument to indicate additive unpack

  // (3) Swap send and receive indices
  for(size_t n=0; n<neighbors.size(); n++)
  {
    int const neighborRank = neighbors[n].neighborRank();
    array1d< localIndex > & nodeGhostsToReceive = nodeManager.getNeighborData( neighborRank ).ghostsToReceive();
    array1d< localIndex > & nodeGhostsToSend = nodeManager.getNeighborData( neighborRank ).ghostsToSend();
    array1d< localIndex > temp = nodeGhostsToSend;
    nodeGhostsToSend = nodeGhostsToReceive;
    nodeGhostsToReceive = temp;
  }

  // (4) Perform sync
  CommunicationTools::getInstance().synchronizePackSendRecvSizes( fieldNames, mesh, neighbors, m_iComm, true );
  parallelDeviceEvents packEvents2;
  CommunicationTools::getInstance().asyncPack( fieldNames, mesh, neighbors, m_iComm, true, packEvents2 );
  waitAllDeviceEvents( packEvents2 );
  CommunicationTools::getInstance().asyncSendRecv( neighbors, m_iComm, true, packEvents2 );
  parallelDeviceEvents unpackEvents2;
  CommunicationTools::getInstance().finalizeUnpack( mesh, neighbors, m_iComm, true, unpackEvents2 );


  // Grid update
  for(int i=0; i<nodeManager.size(); i++)
  {
    if(g_M[i] > 1.0e-12) // small mass threshold
    {
      for(int j=0; j<3; j++)
      {
        g_A[i][j] /= g_M[i]; // g_A holds the nodal forces before this update, divide by mass to obtain acceleration
        g_V[i][j] /= g_M[i]; // g_V holds the nodal momenta before this update, divide by mass to obtain velocity
        g_V[i][j] += g_A[i][j]*dt;
      }
//      // hard-coded rotation about z-axis passing thru mesh center
//      double xRel = g_X[i][0] - 0.5*(m_xGlobalMax[0] - m_xGlobalMin[0]);
//      double yRel = g_X[i][1] - 0.5*(m_xGlobalMax[1] - m_xGlobalMin[1]);
//      double theta = atan2(yRel,xRel);
//      double r = hypot(xRel,yRel);
//      g_V[i][0] = -r*sin(theta)*50.0;
//      g_V[i][1] = r*cos(theta)*50.0;
//      g_V[i][2] = 0.0;
//      // hard-coded simple shear
//      g_V[i][0] = 50.0*(g_X[i][1] - 0.5*(m_xGlobalMax[1] - m_xGlobalMin[1]));
//      g_V[i][1] = 25.0*(g_X[i][0] - 0.5*(m_xGlobalMax[0] - m_xGlobalMin[0]));
//      g_V[i][2] = 0.0;
    }
    else
    {
      for(int j=0; j<3; j++)
      {
        g_A[i][j] = 0.0;
        g_V[i][j] = 0.0;
      }
    }
  }


  // Grid to particle interpolation
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    arrayView2d< real64 > const particleCenter = subRegion.getParticleCenter();
    arrayView2d< real64 > const particleVelocity = subRegion.getParticleVelocity();
    arrayView1d< real64 > const particleMass = subRegion.getParticleMass();
    arrayView1d< real64 > const particleVolume = subRegion.getParticleVolume();
    arrayView1d< real64 > const particleVolume0 = subRegion.getParticleVolume0();
    arrayView3d< real64 > const particleDeformationGradient = subRegion.getParticleDeformationGradient();

    string const & solidMaterialName = subRegion.template getReference< string >( viewKeyStruct::solidMaterialNamesString() );
    ElasticIsotropic & constitutiveRelation = getConstitutiveModel< ElasticIsotropic >( subRegion, solidMaterialName ); // again, limiting to elastic isotropic for now
    arrayView2d< real64 > const particleDensity = constitutiveRelation.getDensity();
    arrayView3d< real64, solid::STRESS_USD > const particleStress = constitutiveRelation.getStress();
    arrayView1d< real64 > const shearModulus = constitutiveRelation.shearModulus();
    arrayView1d< real64 > const bulkModulus = constitutiveRelation.bulkModulus();

    // Particle loop - we might be able to get rid of this someday and have everything happen via MPMParticleSubRegion methods
    for(int p=0; p<subRegion.size(); p++)
    {
      auto const & p_x = particleCenter[p];
      auto const & p_v = particleVelocity[p];
      real64 & p_m = particleMass[p];
      real64 & p_Vol = particleVolume[p];
      real64 const & p_Vol0 = particleVolume0[p];
      real64 & p_rho = particleDensity[p][0];
      auto const & p_F = particleDeformationGradient[p]; // auto = LvArray::ArraySlice<double, 2, 1, long>
      auto const & p_stress = particleStress[p][0];
      real64 p_L[3][3] = { {0} }; // Velocity gradient
      real64 p_FOld[3][3] = { {0} }; // Old particle F
      real64 EG[3][3] = { {0} }; // Green-Lagrange Strain
      real64 PK2[3][3] = { {0} }; // 2nd Piola-Kirchhoff Stress
      real64 sigTemp[3][3] = { {0} }; // Temporary stress-like object
      real64 detF = 0.0; // Material Jacobian

      // Store the old particle F
      for(int i=0; i<3; i++)
      {
        for(int j=0; j<3; j++)
        {
          p_FOld[i][j] = p_F[i][j];
        }
      }

      // Get interpolation kernel
      std::vector<int> nodeIDs; // nodes that the particle maps to
      std::vector<real64> weights; // shape function value for each node
      std::vector< std::vector<real64> > gradWeights; // shape function gradient value for each node; 1st index = direction, 2nd index = node
      gradWeights.resize(3);
      subRegion.getAllWeights(p,
                              p_x,
                              m_xLocalMin,
                              m_hx,
                              m_ijkMap,
                              g_X,
                              nodeIDs,      // output
                              weights,      // output
                              gradWeights); // output

      // Particle-to-grid map
      for(size_t i=0; i<nodeIDs.size(); i++)
      {
        int g = nodeIDs[i];
        for(int j=0; j<3; j++)
        {
          p_x[j] += g_V[g][j]*dt*weights[i];
          p_v[j] += g_A[g][j]*dt*weights[i]; // FLIP
          for(int k=0; k<3; k++)
          {
            p_L[j][k] += g_V[g][j]*gradWeights[k][i];
          }
        }
      }

      // Particle kinematic update - TODO: surely there's a nicer way to do this with LvArray
      // Add identity tensor to velocity gradient and multiply by dt
      for(int i=0; i<3; i++)
      {
        for(int j=0; j<3; j++)
        {
          if(i==j)
            p_L[i][j] = p_L[i][j]*dt + 1.0;
          else
            p_L[i][j] *= dt;
        }
      }

      // Get new F
      for(int i=0; i<3; i++)
      {
        for(int j=0; j<3; j++)
        {
          p_F[i][j] = p_L[i][0]*p_FOld[0][j] + p_L[i][1]*p_FOld[1][j] + p_L[i][2]*p_FOld[2][j]; // matrix multiply
        }
      }

      // Get det(F), update volume and r-vectors
      detF = -p_F[0][2]*p_F[1][1]*p_F[2][0] + p_F[0][1]*p_F[1][2]*p_F[2][0] + p_F[0][2]*p_F[1][0]*p_F[2][1] - p_F[0][0]*p_F[1][2]*p_F[2][1] - p_F[0][1]*p_F[1][0]*p_F[2][2] + p_F[0][0]*p_F[1][1]*p_F[2][2];
      p_Vol = p_Vol0*detF;
      p_rho = p_m/p_Vol;
      subRegion.updateRVectors(p, p_F);

      // Particle constitutive update - Elastic Isotropic model doesn't have a hyperelastic update yet (waiting on strain and stress measure confirmation?) so we implement our own - St. Venant-Kirchhoff
      // Get Green-Lagrange strain
      for(int i=0; i<3; i++)
      {
       for(int j=0; j<3; j++)
       {
         if(i == j)
         {
           EG[i][j] = (p_F[0][i]*p_F[0][j] + p_F[1][i]*p_F[1][j] + p_F[2][i]*p_F[2][j]) - 1.0;
         }
         else
         {
           EG[i][j] = p_F[0][i]*p_F[0][j] + p_F[1][i]*p_F[1][j] + p_F[2][i]*p_F[2][j];
         }
         EG[i][j] *= 0.5;
       }
      }

      // Get PK2 stress
      for(int i=0; i<3; i++)
      {
        for(int j=0; j<3; j++)
        {
          if(i == j)
          {
            real64 lambda = bulkModulus[p] - (2.0/3.0)*shearModulus[p];
            PK2[i][j] = lambda*(EG[0][0] + EG[1][1] + EG[2][2]) + 2.0*shearModulus[p]*EG[i][j];
          }
          else
          {
            PK2[i][j] = 2.0*shearModulus[p]*EG[i][j];
          }
        }
      }

      // Partially convert to Cauchy stress
      for(int i=0; i<3; i++)
      {
        for(int j=0; j<3; j++)
        {
          sigTemp[i][j] = (p_F[i][0]*PK2[0][j] + p_F[i][1]*PK2[1][j] + p_F[i][2]*PK2[2][j])/detF;
        }
      }

      // Finish conversion to Cauchy stress
      for(int i=0; i<3; i++)
      {
        for(int j=i; j<3; j++) // symmetric update, yes this works, I checked it
        {
          int voigt = m_voigtMap[i][j];
          p_stress[voigt] = (sigTemp[i][0]*p_F[j][0] + sigTemp[i][1]*p_F[j][1] + sigTemp[i][2]*p_F[j][2]);
        }
      }

    } // particle loop
  } ); // subregion loop


  // Particle repartitioning
  partition.repartitionMasterParticlesToNeighbors( domain, m_iComm );

  // Calculate stable time step
  real64 wavespeed = 0.0;
  real64 length = std::fmin(m_hx[0],std::fmin(m_hx[1],m_hx[2]));

  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    string const & solidMaterialName = subRegion.template getReference< string >( viewKeyStruct::solidMaterialNamesString() );
    ElasticIsotropic & constitutiveRelation = getConstitutiveModel< ElasticIsotropic >( subRegion, solidMaterialName ); // For the time being we restrict our attention to elastic isotropic solids. TODO: Have all constitutive models automatically calculate a wave speed.
    arrayView2d< real64 > const rho = constitutiveRelation.getDensity();
    arrayView1d< real64 > const g = constitutiveRelation.shearModulus();
    arrayView1d< real64 > const k = constitutiveRelation.bulkModulus();
    for(int p=0; p<subRegion.size(); p++)
    {
      wavespeed = std::max(wavespeed,sqrt((k[p]+(4.0/3.0)*g[p])/rho[p][0]));
    }
  } );

  real64 dtReturn = wavespeed > 1.0e-12 ? m_cflFactor*length/wavespeed : DBL_MAX; // This partitions's dt, make it huge if wavespeed=0.0 (this happens when there are no particles on this partition)
  return dtReturn;
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
