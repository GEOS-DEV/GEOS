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

#include "chrono"
#include "thread"

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
  m_timeIntegrationOption( TimeIntegrationOption::ExplicitDynamic ),
  m_iComm( CommunicationTools::getInstance().getCommID() ),
  m_numContactGroups(),
  m_numContactFlags(),
  m_numVelocityFields(),
  m_damageFieldPartitioning( false ),
  m_hEl{DBL_MAX, DBL_MAX, DBL_MAX},
  m_xLocalMin{DBL_MAX, DBL_MAX, DBL_MAX},
  m_xLocalMax{DBL_MIN, DBL_MIN, DBL_MIN},
  m_xLocalMinNoGhost{0.0, 0.0, 0.0},
  m_xLocalMaxNoGhost{0.0, 0.0, 0.0},
  m_xGlobalMin{0.0, 0.0, 0.0},
  m_xGlobalMax{0.0, 0.0, 0.0},
  m_domainLengths{0.0, 0.0, 0.0},
  m_nEl{0, 0, 0},
  m_ijkMap(),
  m_voigtMap{ {0, 5, 4}, {5, 1, 3}, {4, 3, 2} }
{
  registerWrapper( viewKeyStruct::timeIntegrationOptionString(), &m_timeIntegrationOption ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_timeIntegrationOption ).
    setDescription( "Time integration method. Options are:\n* " + EnumStrings< TimeIntegrationOption >::concat( "\n* " ) );

  registerWrapper( "hEl", &m_hEl ).
    setInputFlag( InputFlags::FALSE ).
    setDescription( "Element dimensions" );

  registerWrapper( "xLocalMin", &m_xLocalMin ).
    setInputFlag( InputFlags::FALSE ).
    setDescription( "Local minimum grid extent including ghosts" );

  registerWrapper( "xLocalMax", &m_xLocalMax ).
    setInputFlag( InputFlags::FALSE ).
    setDescription( "Local maximum grid extent including ghosts" );

  registerWrapper( "xLocalMinNoGhost", &m_xLocalMinNoGhost ).
    setInputFlag( InputFlags::FALSE ).
    setDescription( "Local minimum grid extent excluding ghosts" );

  registerWrapper( "xLocalMaxNoGhost", &m_xLocalMaxNoGhost ).
    setInputFlag( InputFlags::FALSE ).
    setDescription( "Local maximum grid extent excluding ghosts" );

  registerWrapper( "xGlobalMin", &m_xGlobalMin ).
    setInputFlag( InputFlags::FALSE ).
    setDescription( "Global minimum grid extent including ghosts" );

  registerWrapper( "xGlobalMax", &m_xGlobalMax ).
    setInputFlag( InputFlags::FALSE ).
    setDescription( "Global maximum grid extent including ghosts" );

  registerWrapper( "domainLengths", &m_domainLengths ).
    setInputFlag( InputFlags::FALSE ).
    setDescription( "Lengths of each side of the computational domain" );

  registerWrapper( "nEl", &m_nEl ).
    setInputFlag( InputFlags::FALSE ).
    setDescription( "Number of elements in each spatial dimension" );

  registerWrapper( "ijkMap", &m_ijkMap ).
    setInputFlag( InputFlags::FALSE ).
    setDescription( "Map from cell-spaced coordinates to cell ID" );
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


void SolidMechanicsMPM::registerDataOnMesh( Group & meshBodies )
{
  ExecutableGroup::registerDataOnMesh( meshBodies );

  forMeshTargets( meshBodies, [&] ( string const & meshBodyName,
                                    MeshLevel & meshLevel,
                                    arrayView1d< string const > const & regionNames )
  {
    ParticleManager & particleManager = meshLevel.getParticleManager();

    MeshBody const & meshBody = meshBodies.getGroup< MeshBody >( meshBodyName );

    // Set constitutive names on particles
    if( meshBody.hasParticles() )
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

  forMeshTargets( meshBodies, [&] ( string const & meshBodyName,
                                    MeshLevel & meshLevel,
                                    arrayView1d<string const> const & GEOSX_UNUSED_PARAM( regionNames ) )
  {
    MeshBody const & meshBody = meshBodies.getGroup< MeshBody >( meshBodyName );
    if( meshBody.hasParticles() ) // Particle field registration? TODO: What goes here?
    {
      GEOSX_LOG_RANK_0("Registering particle fields.");
    }
    else // Background grid field registration
    {
      NodeManager & nodes = meshLevel.getNodeManager();

      nodes.registerWrapper< array2d< real64 > >( keys::Mass ).
        setPlotLevel( PlotLevel::LEVEL_0 ).
        setRegisteringObjects( this->getName()).
        setDescription( "An array that holds the mass on the nodes." );
      
      nodes.registerWrapper< array3d< real64 > >( keys::Velocity ).
        setPlotLevel( PlotLevel::LEVEL_0 ).
        setRegisteringObjects( this->getName()).
        setDescription( "An array that holds the current velocity on the nodes." );

      nodes.registerWrapper< array3d< real64 > >( viewKeyStruct::momentumString() ).
        setPlotLevel( PlotLevel::LEVEL_0 ).
        setRegisteringObjects( this->getName()).
        setDescription( "An array that holds the current momentum on the nodes." );

      nodes.registerWrapper< array3d< real64 > >( keys::Acceleration ).
        setPlotLevel( PlotLevel::LEVEL_1 ).
        setRegisteringObjects( this->getName()).
        setDescription( "An array that holds the current acceleration on the nodes." );

      nodes.registerWrapper< array3d< real64 > >( viewKeyStruct::forceExternalString() ).
        setPlotLevel( PlotLevel::LEVEL_0 ).
        setRegisteringObjects( this->getName()).
        setDescription( "An array that holds the external forces on the nodes. This includes any boundary"
                        " conditions as well as coupling forces such as hydraulic forces." );

      nodes.registerWrapper< array3d< real64 > >( viewKeyStruct::forceInternalString() ).
        setPlotLevel( PlotLevel::LEVEL_0 ).
        setRegisteringObjects( this->getName()).
        setDescription( "An array that holds the internal forces on the nodes." );

      nodes.registerWrapper< array3d< real64 > >( viewKeyStruct::forceContactString() ).
        setPlotLevel( PlotLevel::LEVEL_0 ).
        setRegisteringObjects( this->getName()).
        setDescription( "An array that holds the contact force." );

      // Group & nodeSets = nodes.sets(); // Do we need to register BC nodes on this?
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

    if( meshBody.hasParticles() ) // Only particle regions will hold actual materials. Background grid currently holds a null material so that the input file parser doesn't complain, but we don't need to actually do anything with it.
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

void SolidMechanicsMPM::initialize( NodeManager & nodeManager,
                                    ParticleManager & particleManager,
                                    SpatialPartition & partition )
{
  // Get grid quantites
  int numNodes = nodeManager.size();
  arrayView2d< real64, nodes::REFERENCE_POSITION_USD > const & gridReferencePosition = nodeManager.referencePosition();
  
  // Get global domain extent
  for(int g=0; g<numNodes; g++)
  {
    for(int i=0; i<3; i++)
    {
      m_xLocalMin[i] = std::min(m_xLocalMin[i],gridReferencePosition[g][i]);
      m_xLocalMax[i] = std::max(m_xLocalMax[i],gridReferencePosition[g][i]);
    }
  }
  for(int i=0; i<3; i++)
  {
    m_xLocalMinNoGhost[i] = partition.getLocalMin()[i];
    m_xLocalMaxNoGhost[i] = partition.getLocalMax()[i];
    m_xGlobalMin[i] = partition.getGlobalMin()[i];
    m_xGlobalMax[i] = partition.getGlobalMax()[i];
    m_domainLengths[i] = m_xLocalMax[i] - m_xLocalMin[i];
  }

  // Get element size
  for(int g=0; g<numNodes; g++)
  {
    for(int i=0; i<3; i++)
    {
      real64 test = gridReferencePosition[g][i] - m_xLocalMin[i]; // By definition, this should always be positive. Furthermore, the gridReferencePosition should only be those on the local partition
      if(test > 0.0) // We're looking for the smallest nonzero distance from the "min" node. TODO: Could be vulnerable to a finite precision bug.
      {
        m_hEl[i] = std::fmin(test, m_hEl[i]);
      }
    }
  }

  // Get number of elements in each direction
  for(int i=0; i<3; i++)
  {
    m_nEl[i] = std::round(m_domainLengths[i]/m_hEl[i]);
  }

  // Create element map
  m_ijkMap.resize( m_nEl[0] + 1, m_nEl[1] + 1, m_nEl[2] + 1 );
  for( int g = 0 ; g < numNodes ; g++ )
  {
    int i = std::round( ( gridReferencePosition[g][0] - m_xLocalMin[0] ) / m_hEl[0] ) ;
    int j = std::round( ( gridReferencePosition[g][1] - m_xLocalMin[1] ) / m_hEl[1] ) ;
    int k = std::round( ( gridReferencePosition[g][2] - m_xLocalMin[2] ) / m_hEl[2] ) ;
    m_ijkMap[i][j][k] = g;
  }

  // Set particle masses based on their volume and density. Set deformation gradient to identity;
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    string const & solidMaterialName = subRegion.template getReference< string >( viewKeyStruct::solidMaterialNamesString() );
    SolidBase & constitutiveRelation = getConstitutiveModel< SolidBase >( subRegion, solidMaterialName ); // For the time being we restrict our attention to elastic isotropic solids. TODO: Have all constitutive models automatically calculate a wave speed.
    arrayView2d< real64 > const particleDensity = constitutiveRelation.getDensity(); // 2d array because there's a density for each quadrature point, we just access with [particle][0]
    arrayView1d< real64 > const particleVolume = subRegion.getParticleVolume();
    arrayView1d< real64 > const particleMass = subRegion.getParticleMass();
    arrayView3d< real64 > const particleDeformationGradient = subRegion.getParticleDeformationGradient();

    // mass
    for(int p=0; p<subRegion.size(); p++)
    {
      particleMass[p] = particleDensity[p][0] * particleVolume[p]; // TODO: This should probably be done in ParticleMeshGenerator...
    }

    // deformation gradient - TODO: there's probably a LvArray function that makes this a one-liner - I don't think the ParticleSubRegionBase constructor can easily initialize this to identity
    for(int p=0; p<subRegion.size(); p++)
    {
      for(int i=0; i<3; i++)
      {
        for(int j=0; j<3; j++)
        {
          particleDeformationGradient[p][i][j] = i==j ? 1.0 : 0.0;
        }
      }
    }

  } );

}

real64 SolidMechanicsMPM::explicitStep( real64 const & GEOSX_UNUSED_PARAM( time_n ),
                                        real64 const & dt,
                                        const int cycleNumber,
                                        DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  #define USE_PHYSICS_LOOP

  // Spatial Partition
  SpatialPartition & partition = dynamic_cast< SpatialPartition & >(domain.getReference< PartitionBase >( keys::partitionManager ) );

  // Get node and particle managers. ***** We implicitly assume that there are exactly two mesh bodies, and that one has particles and one does not. *****
  Group & meshBodies = domain.getMeshBodies();

  MeshBody & meshBody1 = meshBodies.getGroup< MeshBody >(0);
  MeshBody & meshBody2 = meshBodies.getGroup< MeshBody >(1);
  MeshBody & particles = meshBody1.hasParticles() ? meshBody1 : meshBody2;
  MeshBody & grid = !meshBody1.hasParticles() ? meshBody1 : meshBody2;

  ParticleManager & particleManager = particles.getMeshLevel(0).getParticleManager();
  MeshLevel & mesh = grid.getMeshLevel(0);
  NodeManager & nodeManager = mesh.getNodeManager();

  // At time step zero, perform initialization calculations
  if(cycleNumber == 0)
  {
    initialize(nodeManager, particleManager, partition);
  }

  // Get nodal fields
  int numNodes = nodeManager.size();
  arrayView2d< real64, nodes::REFERENCE_POSITION_USD > & gridReferencePosition = nodeManager.referencePosition();
  array2d< real64 > & gridMass = nodeManager.getReference< array2d< real64 > >( keys::Mass );
  array3d< real64 > & gridVelocity = nodeManager.getReference< array3d< real64 > >( keys::Velocity );
  array3d< real64 > & gridMomentum = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::momentumString() );
  array3d< real64 > & gridAcceleration = nodeManager.getReference< array3d< real64 > >( keys::Acceleration );
  array3d< real64 > & gridInternalForce = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::forceInternalString() );
  array3d< real64 > & gridExternalForce = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::forceExternalString() );
  array3d< real64 > & gridContactForce = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::forceContactString() );

  // Get number of contact groups
  int maxLocalGroupNumber = 0; // Maximum contact group number on this partition.
  int maxGlobalGroupNumber; // Maximum contact group number on global domain.
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    arrayView1d< int > const particleGroup = subRegion.getParticleGroup();
    int subregionMaxGroupNumber = 0;
    for(int p=0; p<subRegion.size(); p++)
    {
      if(particleGroup[p] > subregionMaxGroupNumber)
      {
        subregionMaxGroupNumber = particleGroup[p];
      }
    }
    if( subregionMaxGroupNumber > maxLocalGroupNumber)
    {
      maxLocalGroupNumber = subregionMaxGroupNumber;
    }
  } );  
  MPI_Allreduce( &maxLocalGroupNumber,
                 &maxGlobalGroupNumber,
                 1,
                 MPI_INT,
                 MPI_MAX,
                 MPI_COMM_WORLD );

  // Number of contact groups
  m_numContactGroups = maxGlobalGroupNumber + 1;

  // Specified number of damage flags.
  m_numContactFlags = m_damageFieldPartitioning ? 2 : 1;

  // Total number of velocity fields:
  m_numVelocityFields = m_numContactGroups * m_numContactFlags;

  // Resize grid field arrays
  gridMass.resize( numNodes, m_numVelocityFields );
  gridVelocity.resize( numNodes, m_numVelocityFields, 3 );
  gridMomentum.resize( numNodes, m_numVelocityFields, 3 );
  gridAcceleration.resize( numNodes, m_numVelocityFields, 3 );
  gridInternalForce.resize( numNodes, m_numVelocityFields, 3 );
  gridExternalForce.resize( numNodes, m_numVelocityFields, 3 );
  gridContactForce.resize( numNodes, m_numVelocityFields, 3 );

  // Set grid multi-field labels to avoid a VTK output bug
  // TODO: Only doing this on the 1st cycle breaks restarts, can we do better? Why aren't labels part of the restart data?
  std::vector< std::string > keys = { keys::Velocity, viewKeyStruct::momentumString(),keys::Acceleration, viewKeyStruct::forceInternalString(), viewKeyStruct::forceExternalString(), viewKeyStruct::forceContactString() };
  std::vector< std::string > labels1(m_numVelocityFields);
  std::generate( labels1.begin(), labels1.end(), [i=0]() mutable { return std::to_string( i++ ); } );
  string const labels2[] = { "0", "1", "2" };
  for( size_t gridField=0; gridField<keys.size(); gridField++)
  {
    WrapperBase & wrapper = nodeManager.getWrapper< array3d< real64 > >( keys[gridField] );
    wrapper.setDimLabels( 1, labels1 );
    wrapper.setDimLabels( 2, labels2 );
  }

  // Zero out grid fields
  gridMass.zero();
  gridVelocity.zero();
  gridMomentum.zero();
  gridAcceleration.zero();
  gridInternalForce.zero();
  gridExternalForce.zero();
  gridContactForce.zero();

  // Particle to grid interpolation
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    arrayView2d< real64 > const particleCenter = subRegion.getParticleCenter();
    arrayView2d< real64 > const particleVelocity = subRegion.getParticleVelocity();
    arrayView1d< real64 > const particleMass = subRegion.getParticleMass();
    arrayView1d< real64 > const particleVolume = subRegion.getParticleVolume();
    arrayView1d< int > const particleGroup = subRegion.getParticleGroup();

    string const & solidMaterialName = subRegion.template getReference< string >( viewKeyStruct::solidMaterialNamesString() );
    SolidBase & constitutiveRelation = getConstitutiveModel< SolidBase >( subRegion, solidMaterialName );
    arrayView3d< real64, solid::STRESS_USD > const particleStress = constitutiveRelation.getStress();

    for(int p=0; p<subRegion.size(); p++)
    {
      auto const & p_x = particleCenter[p]; // auto = LvArray::ArraySlice<double, 1, 0, int>
      auto const & p_v = particleVelocity[p]; // auto = LvArray::ArraySlice<double, 1, 0, int>
      real64 const & p_m = particleMass[p];
      real64 const & p_vol = particleVolume[p];
      auto const & p_stress = particleStress[p][0];
      int const & p_group = particleGroup[p];

      // Get interpolation kernel
      std::vector<int> nodeIDs; // nodes that the particle maps to
      std::vector<real64> weights; // shape function value for each node
      std::vector< std::vector<real64> > gradWeights; // shape function gradient value for each node; 1st index = direction, 2nd index = node
      gradWeights.resize(3);
      subRegion.getAllWeights( p,
                               p_x,
                               m_xLocalMin,
                               m_hEl,
                               m_ijkMap,
                               gridReferencePosition,
                               nodeIDs,       // output
                               weights,       // output
                               gradWeights ); // output

      // Map to grid
      for(size_t g=0; g<nodeIDs.size(); g++)
      {
        int mappedNode = nodeIDs[g];
        int fieldIndex = p_group;
        gridMass[mappedNode][fieldIndex] += p_m*weights[g];
        for(int i=0; i<3; i++)
        {
          gridMomentum[mappedNode][fieldIndex][i] += p_m*p_v[i]*weights[g];
          for(int k=0; k<3; k++)
          {
            int voigt = m_voigtMap[k][i];
            gridInternalForce[mappedNode][fieldIndex][i] -= p_stress[voigt]*gradWeights[k][g]*p_vol;
          }
        }
      }

    } // particle loop
  } ); // subregion loop


  // Grid MPI operations

  // (1) Initialize
  FieldIdentifiers fieldsToBeSynced;
  fieldsToBeSynced.addFields( FieldLocation::Node, { keys::Mass, viewKeyStruct::momentumString(), viewKeyStruct::forceInternalString(), viewKeyStruct::forceExternalString() } );
  std::vector< NeighborCommunicator > & neighbors = domain.getNeighbors();
  m_iComm.resize( neighbors.size() );

  // (2) Swap send and receive indices so we can sum from ghost to master
  for(size_t n=0; n<neighbors.size(); n++)
  {
    int const neighborRank = neighbors[n].neighborRank();
    array1d< localIndex > & nodeGhostsToReceive = nodeManager.getNeighborData( neighborRank ).ghostsToReceive();
    array1d< localIndex > & nodeGhostsToSend = nodeManager.getNeighborData( neighborRank ).ghostsToSend();
    array1d< localIndex > temp = nodeGhostsToSend;
    nodeGhostsToSend = nodeGhostsToReceive;
    nodeGhostsToReceive = temp;
  }

  // (3) Additive sync
  CommunicationTools::getInstance().synchronizePackSendRecvSizes( fieldsToBeSynced, mesh, neighbors, m_iComm, true );
  parallelDeviceEvents packEvents;
  CommunicationTools::getInstance().asyncPack( fieldsToBeSynced, mesh, neighbors, m_iComm, true, packEvents );
  waitAllDeviceEvents( packEvents );
  CommunicationTools::getInstance().asyncSendRecv( neighbors, m_iComm, true, packEvents );
  parallelDeviceEvents unpackEvents;
  CommunicationTools::getInstance().finalizeUnpack( mesh, neighbors, m_iComm, true, unpackEvents, MPI_SUM ); // needs an extra argument to indicate additive unpack

  // (4) Swap send and receive indices back so we can sync from master to ghost
  for(size_t n=0; n<neighbors.size(); n++)
  {
    int const neighborRank = neighbors[n].neighborRank();
    array1d< localIndex > & nodeGhostsToReceive = nodeManager.getNeighborData( neighborRank ).ghostsToReceive();
    array1d< localIndex > & nodeGhostsToSend = nodeManager.getNeighborData( neighborRank ).ghostsToSend();
    array1d< localIndex > temp = nodeGhostsToSend;
    nodeGhostsToSend = nodeGhostsToReceive;
    nodeGhostsToReceive = temp;
  }

  // (5) Perform sync
  CommunicationTools::getInstance().synchronizePackSendRecvSizes( fieldsToBeSynced, mesh, neighbors, m_iComm, true );
  parallelDeviceEvents packEvents2;
  CommunicationTools::getInstance().asyncPack( fieldsToBeSynced, mesh, neighbors, m_iComm, true, packEvents2 );
  waitAllDeviceEvents( packEvents2 );
  CommunicationTools::getInstance().asyncSendRecv( neighbors, m_iComm, true, packEvents2 );
  parallelDeviceEvents unpackEvents2;
  CommunicationTools::getInstance().finalizeUnpack( mesh, neighbors, m_iComm, true, unpackEvents2 );


  // Determine updated trial velocities before contact enforcement
  for(int fieldIndex=0; fieldIndex<m_numVelocityFields; fieldIndex++)
  {
    for(int g=0; g<numNodes; g++)
    {
      if(gridMass[g][fieldIndex] > 1.0e-12) // small mass threshold
      {
        for(int i=0; i<3; i++)
        {
          real64 totalForce = gridInternalForce[g][fieldIndex][i] + gridExternalForce[g][fieldIndex][i];
          gridAcceleration[g][fieldIndex][i] = totalForce/gridMass[g][fieldIndex];
          gridVelocity[g][fieldIndex][i] = gridMomentum[g][fieldIndex][i]/gridMass[g][fieldIndex];
          gridVelocity[g][fieldIndex][i] += gridAcceleration[g][fieldIndex][i]*dt;
        }
      }
      else
      {
        for(int i=0; i<3; i++)
        {
          gridAcceleration[g][fieldIndex][i] = 0.0;
          gridVelocity[g][fieldIndex][i] = 0.0;
        }
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
    arrayView1d< real64 > const particleInitialVolume = subRegion.getParticleInitialVolume();
    arrayView3d< real64 > const particleDeformationGradient = subRegion.getParticleDeformationGradient();
    arrayView1d< int > const particleGroup = subRegion.getParticleGroup();

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
      real64 & p_vol = particleVolume[p];
      real64 const & p_vol0 = particleInitialVolume[p];
      real64 & p_rho = particleDensity[p][0];
      auto const & p_F = particleDeformationGradient[p]; // auto = LvArray::ArraySlice<double, 2, 1, int>
      auto const & p_stress = particleStress[p][0];
      int const & p_group = particleGroup[p];
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

      // Get interpolation kernel  - TODO: Seems dumb to have to do this twice
      std::vector<int> nodeIDs; // nodes that the particle maps to
      std::vector<real64> weights; // shape function value for each node
      std::vector< std::vector<real64> > gradWeights; // shape function gradient value for each node; 1st index = direction, 2nd index = node
      gradWeights.resize(3);
      subRegion.getAllWeights( p,
                               p_x,
                               m_xLocalMin,
                               m_hEl,
                               m_ijkMap,
                               gridReferencePosition,
                               nodeIDs,       // output
                               weights,       // output
                               gradWeights ); // output

      // Particle-to-grid map
      for(size_t g=0; g<nodeIDs.size(); g++)
      {
        int mappedNode = nodeIDs[g];
        int fieldIndex = p_group;
        for(int i=0; i<3; i++)
        {
          p_x[i] += gridVelocity[mappedNode][fieldIndex][i]*dt*weights[g];
          p_v[i] += gridAcceleration[mappedNode][fieldIndex][i]*dt*weights[g]; // FLIP
          for(int j=0; j<3; j++)
          {
            p_L[i][j] += gridVelocity[mappedNode][fieldIndex][i]*gradWeights[j][g];
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
      p_vol = p_vol0*detF;
      p_rho = p_m/p_vol;
      subRegion.updateRVectors(p, p_F);

      // Particle constitutive update - Elastic Isotropic model doesn't have a hyperelastic update yet (waiting on strain and stress measure confirmation?) so we implement our own - St. Venant-Kirchhoff
      // Get Green-Lagrange strain
      for(int i=0; i<3; i++)
      {
       for(int j=0; j<3; j++)
       {
         EG[i][j] = 0.5*( (p_F[0][i]*p_F[0][j] + p_F[1][i]*p_F[1][j] + p_F[2][i]*p_F[2][j]) - ( i==j ? 1.0 : 0.0 ) );
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
  real64 length = std::fmin(m_hEl[0],std::fmin(m_hEl[1],m_hEl[2]));

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
