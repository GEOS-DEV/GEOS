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
namespace tOps = LvArray::tensorOps; // We call this A LOT, make an alias



SolidMechanicsMPM::SolidMechanicsMPM( const string & name,
                                      Group * const parent ):
  SolverBase( name, parent ),
  m_solverProfiling( 0 ),
  m_timeIntegrationOption( TimeIntegrationOption::ExplicitDynamic ),
  m_iComm( CommunicationTools::getInstance().getCommID() ),
  m_boundaryConditionTypes(),
  m_smallMass(),
  m_numContactGroups(),
  m_numContactFlags(),
  m_numVelocityFields(),
  m_damageFieldPartitioning( false ),
  m_contactGapCorrection( 0 ),
  m_frictionCoefficient( 0.0 ),
  m_planeStrain( false ),
  m_numDims( 3 ),
  m_hEl{ DBL_MAX, DBL_MAX, DBL_MAX },
  m_xLocalMin{ DBL_MAX, DBL_MAX, DBL_MAX },
  m_xLocalMax{ DBL_MIN, DBL_MIN, DBL_MIN },
  m_xLocalMinNoGhost{ 0.0, 0.0, 0.0 },
  m_xLocalMaxNoGhost{ 0.0, 0.0, 0.0 },
  m_xGlobalMin{ 0.0, 0.0, 0.0 },
  m_xGlobalMax{ 0.0, 0.0, 0.0 },
  m_domainLengths{ 0.0, 0.0, 0.0 },
  m_nEl{ 0, 0, 0 },
  m_ijkMap(),
  m_voigtMap{ {0, 5, 4}, {5, 1, 3}, {4, 3, 2} }
{
  registerWrapper( "solverProfiling", &m_solverProfiling ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flag for timing subroutines in the solver" );
  
  registerWrapper( viewKeyStruct::timeIntegrationOptionString(), &m_timeIntegrationOption ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_timeIntegrationOption ).
    setDescription( "Time integration method. Options are:\n* " + EnumStrings< TimeIntegrationOption >::concat( "\n* " ) );

  registerWrapper( "boundaryConditionTypes", &m_boundaryConditionTypes ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Boundary conditions on x-, x+, y-, y+, z- and z+ faces. Options are:\n* " + EnumStrings< BoundaryConditionOption >::concat( "\n* " ) );

  registerWrapper( "numContactGroups", &m_numContactGroups ).
    setInputFlag( InputFlags::FALSE ).
    setDescription( "Number of prescribed contact groups" );

  registerWrapper( "numContactFlags", &m_numContactFlags ).
    setInputFlag( InputFlags::FALSE ).
    setDescription( "Number of contact flags that can appear due to damage" );

  registerWrapper( "numVelocityFields", &m_numVelocityFields ).
    setInputFlag( InputFlags::FALSE ).
    setDescription( "Number of velocity fields" );

  registerWrapper( "damageFieldPartitioning", &m_damageFieldPartitioning ).
    setApplyDefaultValue( false ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flag for using the gradient of the particle damage field to partition material into separate velocity fields" );

  registerWrapper( "contactGapCorrection", &m_contactGapCorrection ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flag for mitigating contact gaps" );

  registerWrapper( "frictionCoefficient", &m_frictionCoefficient ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Coefficient of friction, currently assumed to be the same everywhere" );

  registerWrapper( "planeStrain", &m_planeStrain ).
    setApplyDefaultValue( false ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flag for performing plane strain calculations" );

  registerWrapper( "numDims", &m_numDims ).
    setApplyDefaultValue( 3 ).
    setInputFlag( InputFlags::FALSE ).
    setDescription( "The number of active spatial dimensions, 2 for plane strain, 3 otherwise" );

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
    setDescription( "Map from indices in each spatial dimension to local node ID" );
}

void SolidMechanicsMPM::postProcessInput()
{
  SolverBase::postProcessInput();

  // Set number of active dimensions based on m_planeStrain
  m_numDims = m_planeStrain ? 2 : 3;

  // Throw error if boundary conditions are incorrectly specified
  GEOSX_ERROR_IF( m_boundaryConditionTypes.size() != 6 && m_boundaryConditionTypes.size() > 0,
                  "boundaryConditionTypes must be of length 6. "
                  "The 6 entries correspond to BCs on the x-, x+, y-, y+, z- and z+ faces." );

  // Initialize boundary condition types if they're not specified by the user
  if( m_boundaryConditionTypes.size() == 0 )
  {
    m_boundaryConditionTypes.resize(6);
    tOps::fill< 6 >( m_boundaryConditionTypes, 0 );
  }
}

SolidMechanicsMPM::~SolidMechanicsMPM()
{
  // TODO Auto-generated destructor stub
}


void SolidMechanicsMPM::registerDataOnMesh( Group & meshBodies )
{
  ExecutableGroup::registerDataOnMesh( meshBodies );

  forDiscretizationOnMeshTargets( meshBodies, [&] ( string const & meshBodyName,
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

  forDiscretizationOnMeshTargets( meshBodies, [&] ( string const & meshBodyName,
                                                    MeshLevel & meshLevel,
                                                    arrayView1d< string const > const & GEOSX_UNUSED_PARAM( regionNames ) )
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
        setRegisteringObjects( this->getName() ).
        setDescription( "An array that holds the mass on the nodes." );
      
      nodes.registerWrapper< array3d< real64 > >( keys::Velocity ).
        setPlotLevel( PlotLevel::LEVEL_0 ).
        setRegisteringObjects( this->getName() ).
        setDescription( "An array that holds the current velocity on the nodes." );

      nodes.registerWrapper< array3d< real64 > >( viewKeyStruct::momentumString() ).
        setPlotLevel( PlotLevel::LEVEL_1 ).
        setRegisteringObjects( this->getName() ).
        setDescription( "An array that holds the current momentum on the nodes." );

      nodes.registerWrapper< array3d< real64 > >( keys::Acceleration ).
        setPlotLevel( PlotLevel::LEVEL_1 ).
        setRegisteringObjects( this->getName() ).
        setDescription( "An array that holds the current acceleration on the nodes." );

      nodes.registerWrapper< array3d< real64 > >( viewKeyStruct::forceExternalString() ).
        setPlotLevel( PlotLevel::LEVEL_1 ).
        setRegisteringObjects( this->getName() ).
        setDescription( "An array that holds the external forces on the nodes. This includes any boundary"
                        " conditions as well as coupling forces such as hydraulic forces." );

      nodes.registerWrapper< array3d< real64 > >( viewKeyStruct::forceInternalString() ).
        setPlotLevel( PlotLevel::LEVEL_1 ).
        setRegisteringObjects( this->getName() ).
        setDescription( "An array that holds the internal forces on the nodes." );

      nodes.registerWrapper< array3d< real64 > >( viewKeyStruct::forceContactString() ).
        setPlotLevel( PlotLevel::LEVEL_1 ).
        setRegisteringObjects( this->getName() ).
        setDescription( "An array that holds the contact force on the nodes." );

      nodes.registerWrapper< array2d< real64 > >( viewKeyStruct::damageString() ).
        setPlotLevel( PlotLevel::LEVEL_1 ).
        setRegisteringObjects( this->getName() ).
        setDescription( "An array that holds the result of mapping particle damage to the nodes." );

      nodes.registerWrapper< array2d< real64 > >( viewKeyStruct::maxDamageString() ).
        setPlotLevel( PlotLevel::LEVEL_1 ).
        setRegisteringObjects( this->getName() ).
        setDescription( "An array that holds the maximum damage of any particle mapping to a given node." );

      nodes.registerWrapper< array3d< real64 > >( viewKeyStruct::surfaceNormalString() ).
        setPlotLevel( PlotLevel::LEVEL_1 ).
        setRegisteringObjects( this->getName() ).
        setDescription( "An array that holds the contact surface normals on the nodes." );

      nodes.registerWrapper< array3d< real64 > >( viewKeyStruct::materialPositionString() ).
        setPlotLevel( PlotLevel::LEVEL_1 ).
        setRegisteringObjects( this->getName() ).
        setDescription( "An array that holds the result of mapping particle positions to the nodes." );



      Group & nodeSets = nodes.sets();

      nodeSets.registerWrapper< array1d< SortedArray< localIndex > > >( viewKeyStruct::boundaryNodesString() ).
        setPlotLevel( PlotLevel::NOPLOT ).
        setRestartFlags( RestartFlags::NO_WRITE );

      nodeSets.registerWrapper< array1d< SortedArray< localIndex > > >( viewKeyStruct::bufferNodesString() ).
        setPlotLevel( PlotLevel::NOPLOT ).
        setRestartFlags( RestartFlags::NO_WRITE );
    }
  } );
}

void SolidMechanicsMPM::initializePreSubGroups()
{
  SolverBase::initializePreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );

  Group & meshBodies = domain.getMeshBodies();

  forDiscretizationOnMeshTargets( meshBodies, [&] ( string const & meshBodyName,
                                                                MeshLevel & meshLevel,
                                                                arrayView1d<string const> const & regionNames )
  {
    MeshBody const & meshBody = meshBodies.getGroup< MeshBody >( meshBodyName );

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
  array2d< real64 > & gridMass = nodeManager.getReference< array2d< real64 > >( keys::Mass );
  array2d< real64 > & gridDamage = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::damageString() );
  array2d< real64 > & gridMaxDamage = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::maxDamageString() );
  array3d< real64 > & gridVelocity = nodeManager.getReference< array3d< real64 > >( keys::Velocity );
  array3d< real64 > & gridMomentum = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::momentumString() );
  array3d< real64 > & gridAcceleration = nodeManager.getReference< array3d< real64 > >( keys::Acceleration );
  array3d< real64 > & gridInternalForce = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::forceInternalString() );
  array3d< real64 > & gridExternalForce = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::forceExternalString() );
  array3d< real64 > & gridContactForce = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::forceContactString() );
  array3d< real64 > & gridSurfaceNormal = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::surfaceNormalString() );
  array3d< real64 > & gridMaterialPosition = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::materialPositionString() );
  
  // Get local domain extent
  for(int g=0; g<numNodes; g++)
  {
    for(int i=0; i<3; i++)
    {
      m_xLocalMin[i] = std::fmin( m_xLocalMin[i], gridReferencePosition[g][i] );
      m_xLocalMax[i] = std::fmax( m_xLocalMax[i], gridReferencePosition[g][i] );
    }
  }
  for(int i=0; i<3; i++)
  {
    m_xLocalMinNoGhost[i] = partition.getLocalMin()[i];
    m_xLocalMaxNoGhost[i] = partition.getLocalMax()[i];
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
        m_hEl[i] = std::fmin( test, m_hEl[i] );
      }
    }
  }

  // Get global domain extent excluding buffer nodes
  for(int i=0; i<3; i++)
  {
    m_xGlobalMin[i] = partition.getGlobalMin()[i] + m_hEl[i];
    m_xGlobalMax[i] = partition.getGlobalMax()[i] - m_hEl[i];
  }

  // Get number of elements in each direction
  for(int i=0; i<3; i++)
  {
    m_nEl[i] = std::round( m_domainLengths[i] / m_hEl[i] );
  }

  // Create element map
  m_ijkMap.resize( m_nEl[0] + 1, m_nEl[1] + 1, m_nEl[2] + 1 );
  for( int g=0 ; g<numNodes ; g++ )
  {
    int i = std::round( ( gridReferencePosition[g][0] - m_xLocalMin[0] ) / m_hEl[0] ) ;
    int j = std::round( ( gridReferencePosition[g][1] - m_xLocalMin[1] ) / m_hEl[1] ) ;
    int k = std::round( ( gridReferencePosition[g][2] - m_xLocalMin[2] ) / m_hEl[2] ) ;
    m_ijkMap[i][j][k] = g;
  }

  // Identify node sets for applying boundary conditions. We need boundary nodes and buffer nodes.
  Group & nodeSets = nodeManager.sets();
  array1d< SortedArray< localIndex > > & m_boundaryNodes = nodeSets.getReference< array1d< SortedArray< localIndex > > >( viewKeyStruct::boundaryNodesString() );
  array1d< SortedArray< localIndex > > & m_bufferNodes = nodeSets.getReference< array1d< SortedArray< localIndex > > >( viewKeyStruct::bufferNodesString() );
  m_boundaryNodes.resize(6);
  m_bufferNodes.resize(6);

  for( int face=0; face<6; face++ )
  {
    std::set< localIndex > tmpBoundaryNodes;
    std::set< localIndex > tmpBufferNodes;

    int dir0 = face / 2;           // 0, 0, 1, 1, 2, 2 (x-, x+, y-, y+, z-, z+)

    int positiveNormal = face % 2; // even => (-) => 0, odd => (+) => 1
    int sign = 2 * positiveNormal - 1;
    real64 minOrMax = positiveNormal == 0 ? m_xGlobalMin[dir0] : m_xGlobalMax[dir0];

    for( localIndex g=0; g<numNodes; g++ )
    {
      real64 positionRelativeToBoundary = gridReferencePosition[g][dir0] - minOrMax;
      real64 tolerance = 1.0e-12 * m_hEl[dir0]; // small multiple of element dimension

      if( std::fabs( positionRelativeToBoundary ) < tolerance )
      {
        tmpBoundaryNodes.insert(g);
      }

      if( sign * positionRelativeToBoundary > 0 ) // basically a dot product with the face normal
      {
        tmpBufferNodes.insert(g);
      }
    }

    m_boundaryNodes[face].insert( tmpBoundaryNodes.begin(), tmpBoundaryNodes.end() );
    m_bufferNodes[face].insert( tmpBufferNodes.begin(), tmpBufferNodes.end() );
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

    // Set particle masses and small mass threshold
    real64 minMass = 0.0;
    for(int p=0; p<subRegion.size(); p++)
    {
      particleMass[p] = particleDensity[p][0] * particleVolume[p]; // TODO: This should probably be done in ParticleMeshGenerator...
      minMass = particleMass[p] < minMass ? particleMass[p] : minMass;
    }
    m_smallMass = minMass * 1.0e-12;

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


  // Resize grid arrays according to number of velocity fields

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
  gridDamage.resize( numNodes, m_numVelocityFields );
  gridMaxDamage.resize( numNodes, m_numVelocityFields );
  gridVelocity.resize( numNodes, m_numVelocityFields, 3 );
  gridMomentum.resize( numNodes, m_numVelocityFields, 3 );
  gridAcceleration.resize( numNodes, m_numVelocityFields, 3 );
  gridInternalForce.resize( numNodes, m_numVelocityFields, 3 );
  gridExternalForce.resize( numNodes, m_numVelocityFields, 3 );
  gridContactForce.resize( numNodes, m_numVelocityFields, 3 );
  gridSurfaceNormal.resize( numNodes, m_numVelocityFields, 3 );
  gridMaterialPosition.resize( numNodes, m_numVelocityFields, 3 );
}

real64 SolidMechanicsMPM::explicitStep( real64 const & GEOSX_UNUSED_PARAM( time_n ),
                                        real64 const & dt,
                                        const int cycleNumber,
                                        DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  #define USE_PHYSICS_LOOP

  //#######################################################################################
  solverProfiling( "Get spatial partition, get node and particle managers" );
  //#######################################################################################
  SpatialPartition & partition = dynamic_cast< SpatialPartition & >(domain.getReference< PartitionBase >( keys::partitionManager ) );
  
  // ***** We implicitly assume that there are exactly two mesh bodies, and that one has particles and one does not. *****
  Group & meshBodies = domain.getMeshBodies();

  MeshBody & meshBody1 = meshBodies.getGroup< MeshBody >(0);
  MeshBody & meshBody2 = meshBodies.getGroup< MeshBody >(1);
  MeshBody & particles = meshBody1.hasParticles() ? meshBody1 : meshBody2;
  MeshBody & grid = !meshBody1.hasParticles() ? meshBody1 : meshBody2;

  ParticleManager & particleManager = particles.getBaseDiscretization().getParticleManager();
  MeshLevel & mesh = grid.getBaseDiscretization();
  NodeManager & nodeManager = mesh.getNodeManager();

  //#######################################################################################
  solverProfiling( "Get nodal fields" );
  //#######################################################################################
  // TODO: If we decide to have all subroutines only take particleManager and nodeManager as inputs,
  //       we won't need these at all as they'll all be defined locally in each subroutine.
  int numNodes = nodeManager.size();
  arrayView2d< real64, nodes::REFERENCE_POSITION_USD > & gridReferencePosition = nodeManager.referencePosition();
  array2d< real64 > & gridMass = nodeManager.getReference< array2d< real64 > >( keys::Mass );
  array2d< real64 > & gridDamage = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::damageString() );
  array2d< real64 > & gridMaxDamage = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::maxDamageString() );
  array3d< real64 > & gridVelocity = nodeManager.getReference< array3d< real64 > >( keys::Velocity );
  array3d< real64 > & gridMomentum = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::momentumString() );
  array3d< real64 > & gridAcceleration = nodeManager.getReference< array3d< real64 > >( keys::Acceleration );
  array3d< real64 > & gridInternalForce = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::forceInternalString() );
  array3d< real64 > & gridExternalForce = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::forceExternalString() );
  array3d< real64 > & gridContactForce = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::forceContactString() );
  array3d< real64 > & gridSurfaceNormal = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::surfaceNormalString() );
  array3d< real64 > & gridMaterialPosition = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::materialPositionString() );

  //#######################################################################################
  solverProfiling( "At time step zero, perform initialization calculations" );
  //#######################################################################################
  if(cycleNumber == 0)
  {
    initialize( nodeManager, particleManager, partition );
  }
  

  //#######################################################################################
  solverProfiling( "Set grid multi-field labels to avoid a VTK output bug" );
  //#######################################################################################
  // Must be done every time step despite grid fields being registered
  // TODO: Only doing this on the 1st cycle breaks restarts, can we do better? Why aren't labels part of the restart data?
  setGridFieldLabels( nodeManager );


  //#######################################################################################
  solverProfiling( "Set grid fields to zero" );
  //#######################################################################################
  gridMass.zero();
  gridDamage.zero();
  gridMaxDamage.zero();
  gridVelocity.zero();
  gridMomentum.zero();
  gridAcceleration.zero();
  gridInternalForce.zero();
  gridExternalForce.zero();
  gridContactForce.zero();
  gridSurfaceNormal.zero();
  gridMaterialPosition.zero();


  //#######################################################################################
  solverProfiling( "Particle-to-grid interpolation" );
  //#######################################################################################
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
        gridMass[mappedNode][fieldIndex] += p_m * weights[g];
        for(int i=0; i<3; i++)
        {
          gridMomentum[mappedNode][fieldIndex][i] += p_m * p_v[i] * weights[g];
          gridMaterialPosition[mappedNode][fieldIndex][i] += p_m * ( p_x[i] - gridReferencePosition[mappedNode][i] ) * weights[g];
          for(int k=0; k<3; k++)
          {
            int voigt = m_voigtMap[k][i];
            gridInternalForce[mappedNode][fieldIndex][i] -= p_stress[voigt] * gradWeights[k][g] * p_vol;
          }
        }
      }

    } // particle loop
  } ); // subregion loop


  //#######################################################################################
  solverProfiling( "Grid MPI operations" );
  //#######################################################################################
  std::vector< std::string > fieldNames = { keys::Mass,
                                            viewKeyStruct::momentumString(),
                                            viewKeyStruct::forceInternalString(),
                                            viewKeyStruct::forceExternalString() };
  syncGridFields( fieldNames, domain, nodeManager, mesh );


  //#######################################################################################
  solverProfiling( "Determine trial momenta and velocities based on acceleration due to internal and external forces, but before contact enforcement" );
  //#######################################################################################
  for( int fieldIndex=0; fieldIndex<m_numVelocityFields; fieldIndex++ )
  {
    for( int g=0; g<numNodes; g++ )
    {
      if( gridMass[g][fieldIndex] > m_smallMass ) // small mass threshold
      {
        for( int i=0; i<3; i++ )
        {
          real64 totalForce = gridInternalForce[g][fieldIndex][i] + gridExternalForce[g][fieldIndex][i];
          gridAcceleration[g][fieldIndex][i] = totalForce / gridMass[g][fieldIndex];
          gridMomentum[g][fieldIndex][i] += totalForce * dt;
          gridVelocity[g][fieldIndex][i] = gridMomentum[g][fieldIndex][i] / gridMass[g][fieldIndex];
          gridMaterialPosition[g][fieldIndex][i] /= gridMass[g][fieldIndex];
        }
      }
      else
      {
        for( int i=0; i<3; i++ )
        {
          gridAcceleration[g][fieldIndex][i] = 0.0;
          gridVelocity[g][fieldIndex][i] = 0.0;
          gridMomentum[g][fieldIndex][i] = 0.0;
          gridMaterialPosition[g][fieldIndex][i] = gridReferencePosition[g][i];
        }
      }
    }
  }


  //#######################################################################################
  solverProfiling( "Contact enforcement" );
  //#######################################################################################
  if( m_numVelocityFields > 1 )
  {
    // Compute grid surface normals
    computeGridSurfaceNormals( particleManager, gridReferencePosition, gridSurfaceNormal );

    // Sync surface normals
    syncGridFields( { viewKeyStruct::surfaceNormalString() }, domain, nodeManager, mesh );

    // Apply symmetry boundary conditions to surface normals
    enforceGridVectorFieldSymmetryBC( gridSurfaceNormal, gridReferencePosition, nodeManager.sets() );

    // Normalize grid surface normals
    normalizeGridSurfaceNormals( gridMass, gridSurfaceNormal );

    // Compute contact forces
    computeContactForces( dt,
                          gridMass,
                          gridDamage,
                          gridMaxDamage,
                          gridVelocity,
                          gridMomentum,
                          gridSurfaceNormal,
                          gridMaterialPosition,
                          gridContactForce ); // output

    // Update grid momenta and velocities based on contact forces
    for( int fieldIndex=0; fieldIndex<m_numVelocityFields; fieldIndex++ )
    {
      for( int g=0; g<numNodes; g++ )
      {
        if( gridMass[g][fieldIndex] > m_smallMass ) // small mass threshold
        {
          for( int i=0; i<3; i++ )
          {
            real64 contactForce = gridContactForce[g][fieldIndex][i];
            gridAcceleration[g][fieldIndex][i] += contactForce / gridMass[g][fieldIndex];
            gridMomentum[g][fieldIndex][i] += contactForce * dt;
            gridVelocity[g][fieldIndex][i] = gridMomentum[g][fieldIndex][i] / gridMass[g][fieldIndex];
          }
        }
      }
    }
  }

  // Read F-Table?


  //#######################################################################################
  solverProfiling( "Apply essential boundary conditions" );
  //#######################################################################################
  enforceGridVectorFieldSymmetryBC( gridVelocity, gridReferencePosition, nodeManager.sets() );
  enforceGridVectorFieldSymmetryBC( gridAcceleration, gridReferencePosition, nodeManager.sets() );


  //#######################################################################################
  solverProfiling( "Grid-to-particle interpolation" );
  //#######################################################################################
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
      real64 p_D[3][3] = { {0} }; // Rate of deformation
      real64 p_Diso[3][3] = { {0} }; // Rate of deformation
      real64 p_Ddev[3][3] = { {0} }; // Rate of deformation
      real64 p_FOld[3][3] = { {0} }; // Old particle F
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
          p_x[i] += gridVelocity[mappedNode][fieldIndex][i] * dt * weights[g];
          p_v[i] += gridAcceleration[mappedNode][fieldIndex][i] * dt * weights[g]; // FLIP
          for(int j=0; j<3; j++)
          {
            p_L[i][j] += gridVelocity[mappedNode][fieldIndex][i] * gradWeights[j][g];
          }
        }
      }

      // Get D
      for(int i=0; i<3; i++)
      {
        for(int j=0; j<3; j++)
        {
          p_D[i][j] = 0.5 * ( p_L[i][j] + p_L[j][i] );
        }
      }

      // Get Diso
      for(int i=0; i<3; i++)
      {
        for(int j=0; j<3; j++)
        {
          p_Diso[i][j] = i == j ? ( p_D[0][0] + p_D[1][1] + p_D[2][2] ) / 3.0 : 0.0;
        }
      }

      // Get Ddev
      for(int i=0; i<3; i++)
      {
        for(int j=0; j<3; j++)
        {
          p_Ddev[i][j] = p_D[i][j] - p_Diso[i][j];
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

      // Particle constitutive update - hyperelastic model from vortex MMS paper
      real64 lambda = bulkModulus[p] - 2.0 * shearModulus[p] / 3.0;
      real64 stressFull[3][3] = { { 0.0 } };

      // Populate full stress tensor
      for( int i = 0 ; i < 3 ; i++ )
      {
        for( int j = 0 ; j < 3 ; j++ )
        {
          stressFull[i][j] = i == j ? ( lambda * log( detF ) / detF - ( shearModulus[p] / detF ) ) : 0.0;

          for( int k = 0 ; k < 3 ; k++ )
          {
            stressFull[i][j] += ( shearModulus[p] / detF ) * p_F[i][k] * p_F[j][k];
          }
        }
      }

      // Assign full stress tensor to the symmetric version
      for(int i=0; i<3; i++)
      {
        for(int j=i; j<3; j++) // symmetric update, yes this works, I checked it
        {
          int voigt = m_voigtMap[i][j];
          p_stress[voigt] = stressFull[i][j];
        }
      }

    } // particle loop
  } ); // subregion loop


  // Resize grid based on F-table


  //#######################################################################################
  solverProfiling( "Particle repartitioning" );
  //#######################################################################################
  partition.repartitionMasterParticlesToNeighbors( domain, m_iComm );


  //#######################################################################################
  solverProfiling( "Delete particles that map outside the domain (including buffer cells)" );
  //#######################################################################################
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    subRegion.deleteOutOfRangeParticles( m_xGlobalMin, m_xGlobalMax, m_hEl );
  } );


  //#######################################################################################
  solverProfiling( "Calculate stable time step" );
  //#######################################################################################
  real64 wavespeed = 0.0;
  real64 length = std::fmin( m_hEl[0], std::fmin( m_hEl[1], m_hEl[2] ) );

  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    string const & solidMaterialName = subRegion.template getReference< string >( viewKeyStruct::solidMaterialNamesString() );
    ElasticIsotropic & constitutiveRelation = getConstitutiveModel< ElasticIsotropic >( subRegion, solidMaterialName ); // For the time being we restrict our attention to elastic isotropic solids. TODO: Have all constitutive models automatically calculate a wave speed.
    arrayView2d< real64 > const rho = constitutiveRelation.getDensity();
    arrayView1d< real64 > const g = constitutiveRelation.shearModulus();
    arrayView1d< real64 > const k = constitutiveRelation.bulkModulus();
    for(int p=0; p<subRegion.size(); p++)
    {
      wavespeed = std::max( wavespeed, sqrt( ( k[p] + (4.0/3.0) * g[p] ) / rho[p][0] ) );
    }
  } );

  real64 dtReturn = wavespeed > 1.0e-12 ? m_cflFactor*length/wavespeed : DBL_MAX; // This partitions's dt, make it huge if wavespeed=0.0 (this happens when there are no particles on this partition)


  //#######################################################################################
  solverProfiling( "End of explicitStep" );
  //#######################################################################################
  if( m_solverProfiling )
  {
    // Use MPI reduction to get the average elapsed time for each step on all partitions
    int rank = MpiWrapper::commRank( MPI_COMM_GEOSX );
    unsigned int numIntervals = m_profilingTimes.size() - 1;
    std::vector< real64 > timeIntervalsAllRanks( numIntervals );    

    // Get total CPU time for the entire time step
    real64 totalStepTimeThisRank = m_profilingTimes[numIntervals] - m_profilingTimes[0];
    real64 totalStepTimeAllRanks;
    MpiWrapper::allReduce< real64 >( &totalStepTimeThisRank,
                                     &totalStepTimeAllRanks,
                                     1,
                                     MPI_SUM,
                                     MPI_COMM_GEOSX );

    // Get total CPU times for each queried time interval
    for( unsigned int i = 0 ; i < numIntervals; i++ )
    {
      real64 timeIntervalThisRank = ( m_profilingTimes[i+1] - m_profilingTimes[i] );
      real64 timeIntervalAllRanks;
      MpiWrapper::allReduce< real64 >( &timeIntervalThisRank,
                                       &timeIntervalAllRanks,
                                       1,
                                       MPI_SUM,
                                       MPI_COMM_GEOSX );
      if( rank == 0 )
      {
        timeIntervalsAllRanks[i] = timeIntervalAllRanks;
      }
    }

    // Print out solver profiling
    MPI_Barrier( MPI_COMM_GEOSX );
    if( rank == 0 )
    {
      std::cout << "---------------------------------------------" << std::endl;
      std::cout << "Fraction of total time for one step: " << std::endl;
      for( unsigned int i = 0 ; i < numIntervals ; i++ )
      {
        std::cout << " (" << i << ") ";
        std::cout << std::fixed;
        std::cout << std::showpoint;
        std::cout << std::setprecision(6) << timeIntervalsAllRanks[i] / totalStepTimeAllRanks;
        std::cout << ": " << m_profilingLabels[i] << std::endl;
      }
      std::cout << " ** Total step CPU time:  " << totalStepTimeAllRanks << " s **" << std::endl;
      std::cout << "---------------------------------------------" << std::endl;
    }

    // Reset profiling arrays
    m_profilingTimes.resize(0);
    m_profilingLabels.resize(0);
  }

  return dtReturn;
}

void SolidMechanicsMPM::syncGridFields( std::vector< std::string > const & fieldNames,
                                        DomainPartition & domain,
                                        NodeManager & nodeManager,
                                        MeshLevel & mesh )
{
  // (1) Initialize
  FieldIdentifiers fieldsToBeSynced;
  fieldsToBeSynced.addFields( FieldLocation::Node, fieldNames );
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
}

void SolidMechanicsMPM::enforceGridVectorFieldSymmetryBC( array3d< real64 > & vectorMultiField,
                                                          arrayView2d< real64, nodes::REFERENCE_POSITION_USD > const & gridReferencePosition,
                                                          Group & nodeSets )
{
  array1d< SortedArray< localIndex > > & m_boundaryNodes = nodeSets.getReference< array1d< SortedArray< localIndex > > >( viewKeyStruct::boundaryNodesString() );
  array1d< SortedArray< localIndex > > & m_bufferNodes = nodeSets.getReference< array1d< SortedArray< localIndex > > >( viewKeyStruct::bufferNodesString() );

  for( int fieldIndex=0; fieldIndex<m_numVelocityFields; fieldIndex++ )
  {
    for( int face=0; face<6; face++ )
    {
      if( m_boundaryConditionTypes[face] == 1 )
      {
        // Face-assocaited quantities
        int dir0 = face / 2;           // 0, 0, 1, 1, 2, 2 (x-, x+, y-, y+, z-, z+)
        int dir1 = ( dir0 + 1 ) % 3;   // 1, 1, 2, 2, 0, 0
        int dir2 = ( dir0 + 2 ) % 3;   // 2, 2, 0, 0, 1, 1
        int positiveNormal = face % 2; // even => (-) => 0, odd => (+) => 1
        
        // Enforce BCs on boundary nodes
        for( const auto & g : m_boundaryNodes[face] )
        {
          vectorMultiField[g][fieldIndex][dir0] = 0.0;
        }
      
        // Perform field reflection on buffer nodes
        for( const auto & g : m_bufferNodes[face] )
        {
          int ijk[3];
          ijk[dir0] = positiveNormal * ( m_nEl[dir0] - 2 ) + ( 1 - positiveNormal ) * ( 2 );
          ijk[dir1] = std::round( ( gridReferencePosition[g][dir1] - m_xLocalMin[dir1] ) / m_hEl[dir1] );
          ijk[dir2] = std::round( ( gridReferencePosition[g][dir2] - m_xLocalMin[dir2] ) / m_hEl[dir2] );

          localIndex gFrom = m_ijkMap[ijk[0]][ijk[1]][ijk[2]];

          vectorMultiField[g][fieldIndex][dir0] = -vectorMultiField[gFrom][fieldIndex][dir0]; // Negate component aligned with surface normal
          vectorMultiField[g][fieldIndex][dir1] =  vectorMultiField[gFrom][fieldIndex][dir1];
          vectorMultiField[g][fieldIndex][dir2] =  vectorMultiField[gFrom][fieldIndex][dir2];
        }
      }
    }
  }
}

void SolidMechanicsMPM::computeGridSurfaceNormals( ParticleManager & particleManager,
                                                   arrayView2d< real64, nodes::REFERENCE_POSITION_USD > const & gridReferencePosition,
                                                   array3d< real64 > & gridSurfaceNormal )
{
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    arrayView1d< real64 > const particleVolume = subRegion.getParticleVolume();
    arrayView1d< int > const particleGroup = subRegion.getParticleGroup();

    for(int p=0; p<subRegion.size(); p++)
    {
      // Get interpolation kernel
      std::vector<int> nodeIDs; // nodes that the particle maps to
      std::vector<real64> weights; // shape function value for each node
      std::vector< std::vector<real64> > gradWeights; // shape function gradient value for each node; 1st index = direction, 2nd index = node
      gradWeights.resize(3);
      subRegion.getAllWeights( p,
                               m_xLocalMin,
                               m_hEl,
                               m_ijkMap,
                               gridReferencePosition,
                               nodeIDs,       // output
                               weights,       // output
                               gradWeights ); // output

      // Map to grid
      for( size_t g=0; g<nodeIDs.size(); g++ )
      {
        int mappedNode = nodeIDs[g];
        int fieldIndex = particleGroup[p];
        for(int i=0; i<3; i++)
        {
          gridSurfaceNormal[mappedNode][fieldIndex][i] += gradWeights[i][g] * particleVolume[p];
        }
      }

    } // particle loop
  } ); // subregion loop
}

void SolidMechanicsMPM::normalizeGridSurfaceNormals( array2d< real64 > & gridMass,
                                                     array3d< real64 > & gridSurfaceNormal )
{
  for( localIndex g=0 ; g<gridSurfaceNormal.size(0) ; g++ )
  {
    for( localIndex fieldIndex = 0 ; fieldIndex < m_numVelocityFields; fieldIndex++ )
    {
      if( gridMass[g][fieldIndex] > m_smallMass ) // small mass threshold
      {
        real64 norm = tOps::l2Norm< 3 >( gridSurfaceNormal[g][fieldIndex] );
        if( norm > 1.0e-12 )
        {
          tOps::scale< 3 >( gridSurfaceNormal[g][fieldIndex], 1.0/norm );
        }
        else
        {
          tOps::fill< 3 >( gridSurfaceNormal[g][fieldIndex], 0.0 );
        }
      }
      else
      {
        tOps::fill< 3 >( gridSurfaceNormal[g][fieldIndex], 0.0 );
      }
    }
  }
}

void SolidMechanicsMPM::computeContactForces( const real64 dt,
                                              const array2d< real64 > & gridMass,
                                              const array2d< real64 > & GEOSX_UNUSED_PARAM( gridDamage ),
                                              const array2d< real64 > & GEOSX_UNUSED_PARAM( gridMaxDamage ),
                                              const array3d< real64 > & gridVelocity,
                                              const array3d< real64 > & gridMomentum,
                                              const array3d< real64 > & gridSurfaceNormal,
                                              const array3d< real64 > & gridMaterialPosition,
                                              array3d< real64 > & gridContactForce )
{
  // Get number of nodes
  int numNodes = gridMass.size(0);

  // Initialize contact force to zero
  gridContactForce.zero();

  for( localIndex g=0 ; g<numNodes ; g++ )
  {
    // Loop over all possible field pairings and enforce contact on each pair.
    // gridContactForce will be gradually updated due to a '+=' in computePairwiseNodalContactForce
    for( localIndex A = 0 ; A < m_numVelocityFields - 1 ; A++ )
    {
      for( localIndex B = A + 1 ; B < m_numVelocityFields ; B++ )
      {
        // Make sure both fields in the pair are active
        bool active = ( gridMass[g][A] > m_smallMass ) && ( tOps::l2NormSquared< 3 >( gridSurfaceNormal[g][A] ) > 1.0e-12 )
                      and
                      ( gridMass[g][B] > m_smallMass ) && ( tOps::l2NormSquared< 3 >( gridSurfaceNormal[g][B] ) > 1.0e-12 );

        if( active )
        {
          // Evaluate the separability criterion for the contact pair.
  //      int separable = evaluateSeparabilityCriterion( A,
  //                                                     B,
  //                                                     gridDamage[g][A],
  //                                                     gridDamage[g][B],
  //                                                     gridMaxDamage[A],
  //                                                     gridMaxDamage[B] );
          int separable = 1;

          computePairwiseNodalContactForce( separable,
                                            dt,
                                            gridMass[g][A],
                                            gridMass[g][B],
                                            gridVelocity[g][A],
                                            gridVelocity[g][B],
                                            gridMomentum[g][A],
                                            gridMomentum[g][B],
                                            gridSurfaceNormal[g][A],
                                            gridSurfaceNormal[g][B],
                                            gridMaterialPosition[g][A],
                                            gridMaterialPosition[g][B],
                                            gridContactForce[g][A],
                                            gridContactForce[g][B] );
        }
      }
    }
  }
}

void SolidMechanicsMPM::computePairwiseNodalContactForce( const int & separable,
                                                          const real64 & dt,
                                                          const real64 & mA,
                                                          const real64 & mB,
                                                          const arraySlice1d< real64 > vA,
                                                          const arraySlice1d< real64 > GEOSX_UNUSED_PARAM( vB ),
                                                          const arraySlice1d< real64 > qA,
                                                          const arraySlice1d< real64 > qB,
                                                          const arraySlice1d< real64 > nA,
                                                          const arraySlice1d< real64 > nB,
                                                          const arraySlice1d< real64 > xA, // Position of field A
                                                          const arraySlice1d< real64 > xB, // Position of field B
                                                          arraySlice1d< real64 > fA,
                                                          arraySlice1d< real64 > fB )
{
  // Total mass for the contact pair.
  real64 mAB = mA + mB;

  // Outward normal of field A with respect to field B.

  // Use the surface normal for whichever field has more mass.
  // This should be good for corners against flat surfaces.

  // nAB = mA > mB ? nA : -nB;
  array1d< real64 > nAB(3);
  if( mA > mB )
  {
    tOps::copy< 3 >( nAB, nA );
  }
  else
  {
    tOps::scaledCopy< 3 >( nAB, nB, -1.0 );
  }
  tOps::normalize< 3 >( nAB );

  // gap = ( xB - xA ) * nAB - minCellLength;
  real64 minCellLength = m_planeStrain ? std::min( m_hEl[0], m_hEl[1] ) : std::min( m_hEl[0], std::min( m_hEl[1], m_hEl[2] ) );
  real64 gap = subtractDot( xB, xA, nAB ) - minCellLength;

  // Total momentum for the contact pair.
  array1d< real64 > qAB(3);
  tOps::copy< 3 >( qAB, qA );
  tOps::add< 3 >( qAB, qB );

  // Center-of-mass velocity for the contact pair.
  array1d< real64 > vAB(3);
  tOps::scaledCopy< 3 >( vAB, qAB, 1.0/mAB );

  // Compute s1AB and s2AB, to form an orthonormal basis.  This uses the method by E. Herbold, to ensure
  // consistency between surfaces.
  array1d< real64 > s1AB(3), s2AB(3); // Tangential vectors for the contact pair
  computeOrthonormalBasis( nAB, s1AB, s2AB );

  // Check for separability, and enforce either slip, or no-slip contact, accordingly
  if( separable == 0 )
  {
    // Surfaces are bonded, treat as single velocity field by applying normal force to prevent
    // interpenetration, and tangential force to prevent slip.
    real64 fnor =  ( mA / dt ) * subtractDot( vAB, vA, nAB ),
           ftan1 = ( mA / dt ) * subtractDot( vAB, vA, s1AB ),
           ftan2 = ( mA / dt ) * subtractDot( vAB, vA, s2AB );

    tOps::scaledCopy< 3 >( fA, nAB, fnor );
    tOps::scaledAdd< 3 >( fA, s1AB, ftan1 );
    tOps::scaledAdd< 3 >( fA, s2AB, ftan2 );
    tOps::scaledCopy< 3 >( fB, fA, -1.0 );
  }
  else
  {
    // Surfaces are separable, for frictional contact, apply a normal force to
    // prevent interpenetration, and tangential force to prevent slip unless f_tan > mu*f_nor

    real64 contact = 0.0;
    if( m_contactGapCorrection == 1 )
    {
      contact = ( subtractDot( vA, vAB, nAB ) > 0 ) && ( gap <= 0 ) ? 1. : 0.;
    }
    else if( m_contactGapCorrection == 2 )
    {
      if( subtractDot( vA, vAB, nAB ) > 0 )
      {
        if( gap <= 0.0 )
        {
          contact = 1.0;
        }
        else if( gap < minCellLength )
        {
          contact = 1.0 - gap/minCellLength;
        }
      }
    }
    else
    {
      contact = subtractDot( vA, vAB, nAB ) > 0 ? 1. : 0.;
    }

    // Determine forces for no normal interpenetration and tangential sticking
    real64 fnor = contact * ( mA / dt ) * subtractDot( vAB, vA, nAB ),
           ftan1 = ( mA / dt ) * subtractDot( vAB, vA, s1AB ),
           ftan2 = ( mA / dt ) * subtractDot( vAB, vA, s2AB );

    // Get direction of tangential contact force
    array1d< real64 > sAB(3);
    tOps::scaledCopy< 3 >( sAB, s1AB, ftan1 );
    tOps::scaledAdd< 3 >( sAB, s2AB, ftan2 );
    tOps::normalize< 3 >( sAB );

    // Update fA and fB - tangential force is bounded by friction
    tOps::scaledAdd< 3 >( fA, nAB, fnor );
    tOps::scaledAdd< 3 >( fA, sAB, std::min( m_frictionCoefficient * std::abs( fnor ), sqrt( ftan1 * ftan1 + ftan2 * ftan2 ) ) );
    tOps::scaledAdd< 3 >( fB, fA, -1.0 );
  }
}

real64 SolidMechanicsMPM::subtractDot( const arraySlice1d< real64 > & u,
                                       const arraySlice1d< real64 > & v,
                                       const arraySlice1d< real64 > & n )
{
  // Calculates (u-v).n
  // This is done so much in the contact subroutines that it's nice to have a separate function for it
  array1d< real64 > u_minus_v(3);
  tOps::copy< 3 >( u_minus_v, u );
  tOps::subtract< 3 >( u_minus_v, v );
  return tOps::AiBi< 3 >( u_minus_v, n );
}

void SolidMechanicsMPM::computeOrthonormalBasis( const array1d< real64 > & e1, // input "normal" unit vector.
                                                 array1d< real64 > & e2,       // output "tangential" unit vector.
                                                 array1d< real64 > & e3 )      // output "tangential" unit vector.
{
  // This routine takes in a normalized vector and gives two orthogonal vectors.
  // It is rather arbitrary, in general, how this is done; however, this routine
  // is written in a consistent way. This is important, for example, if you
  // have a face on one processor and a ghost face on another processor and
  // you try to define a basis that is EXACTLY the same for debugging purposes.
  // This came up when trying to debug shear traction components and this
  // routine helped a lot.

  real64 e1x = e1[0], e1y = e1[1], e1z = e1[2],
         e2x, e2y, e2z;

  // find the maximum pair of values from the first vector.
  real64 maxp1 = std::abs( e1x ) + std::abs( e1y ),
         maxp2 = std::abs( e1y ) + std::abs( e1z ),
         maxp3 = std::abs( e1x ) + std::abs( e1z );

  // This part amounts to setting the smallest component of the input vector to 0.
  // Then, the orthogonal vector is a simple operation in the 2D plane.
  if( maxp1 > maxp2 && maxp1 > maxp3 )
  {
    e2x = e1y;
    e2y = -e1x;
    e2z = 0.0;
  }
  else if( maxp2 > maxp1 && maxp2 > maxp3 )
  {
    e2y = e1z;
    e2z = -e1y;
    e2x = 0.0;
  }
  else
  {
    e2z = e1x;
    e2x = -e1z;
    e2y = 0.0;
  }

  // Normalize the first base vector
  real64 e2norm = sqrt( e2x * e2x + e2y * e2y + e2z * e2z );
  e2[0] = e2x / e2norm;
  e2[1] = e2y / e2norm;
  e2[2] = e2z / e2norm;

  // The second base vector is simply a cross product of the first two.  For
  // safety sake, this may need to be normalized to combat issues of numerical
  // precision.
  e3[0] =  e1y * e2z - e1z * e2y;
  e3[1] = -e1x * e2z + e1z * e2x;
  e3[2] =  e1x * e2y - e1y * e2x;
  tOps::normalize< 3 >( e3 );
}

void SolidMechanicsMPM::setGridFieldLabels( NodeManager & nodeManager )
{
  // Generate labels
  std::vector< std::string > labels1(m_numVelocityFields);
  std::generate( labels1.begin(), labels1.end(), [i=0]() mutable { return "velocityField" + std::to_string( i++ ); } );
  string const labels2[] = { "X", "Y", "Z" };

  // Apply labels to scalar multi-fields
  std::vector< std::string > keys2d = { keys::Mass,
                                        viewKeyStruct::damageString(),
                                        viewKeyStruct::maxDamageString()};
  for( size_t gridField=0; gridField<keys2d.size(); gridField++)
  {
    WrapperBase & wrapper = nodeManager.getWrapper< array2d< real64 > >( keys2d[gridField] );
    wrapper.setDimLabels( 1, labels1 );
  }

  // Apply labels to vector multi-fields
  std::vector< std::string > keys3d = { keys::Velocity,
                                        viewKeyStruct::momentumString(),
                                        keys::Acceleration,
                                        viewKeyStruct::forceInternalString(),
                                        viewKeyStruct::forceExternalString(),
                                        viewKeyStruct::forceContactString(),
                                        viewKeyStruct::surfaceNormalString(),
                                        viewKeyStruct::materialPositionString() };
  for( size_t gridField=0; gridField<keys3d.size(); gridField++)
  {
    WrapperBase & wrapper = nodeManager.getWrapper< array3d< real64 > >( keys3d[gridField] );
    wrapper.setDimLabels( 1, labels1 );
    wrapper.setDimLabels( 2, labels2 );
  }
}

void SolidMechanicsMPM::solverProfiling( std::string label )
{
  if( m_solverProfiling )
  {
    MPI_Barrier( MPI_COMM_GEOSX );
    m_profilingTimes.push_back( MPI_Wtime() );
    m_profilingLabels.push_back( label );
  }
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

