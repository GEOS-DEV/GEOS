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
#include "kernels/ImplicitSmallStrainNewmark.hpp"
#include "kernels/ImplicitSmallStrainQuasiStatic.hpp"
#include "kernels/ExplicitSmallStrain.hpp"
#include "kernels/ExplicitFiniteStrain.hpp"

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
#include "MPMSolverBaseFields.hpp"


namespace geosx
{

using namespace dataRepository;
using namespace constitutive;
namespace tOps = LvArray::tensorOps;

SolidMechanicsMPM::SolidMechanicsMPM( const string & name,
                                      Group * const parent ):
  SolverBase( name, parent ),
  m_solverProfiling( 0 ),
  m_timeIntegrationOption( TimeIntegrationOption::ExplicitDynamic ),
  m_iComm( CommunicationTools::getInstance().getCommID() ),
  m_prescribedBcTable( 0 ),
  m_boundaryConditionTypes(),
  m_bcTable(),
  m_prescribedBoundaryFTable( 0 ),
  m_fTableInterpType( 0 ),
  m_fTable(),
  m_needsNeighborList( 0 ),
  m_neighborRadius( -1.0 ),
  m_binSizeMultiplier( 1 ),
  m_useDamageAsSurfaceFlag( 0 ),
  m_cpdiDomainScaling( 0 ),
  m_smallMass( DBL_MAX ),
  m_numContactGroups(),
  m_numContactFlags(),
  m_numVelocityFields(),
  m_damageFieldPartitioning( 0 ),
  m_contactGapCorrection( 0 ),
  m_frictionCoefficient( 0.0 ),
  m_planeStrain( 0 ),
  m_numDims( 3 ),
  m_hEl{ DBL_MAX, DBL_MAX, DBL_MAX },
  m_xLocalMin{ DBL_MAX, DBL_MAX, DBL_MAX },
  m_xLocalMax{ DBL_MIN, DBL_MIN, DBL_MIN },
  m_xLocalMinNoGhost{ 0.0, 0.0, 0.0 },
  m_xLocalMaxNoGhost{ 0.0, 0.0, 0.0 },
  m_xGlobalMin{ 0.0, 0.0, 0.0 },
  m_xGlobalMax{ 0.0, 0.0, 0.0 },
  m_partitionExtent{ 0.0, 0.0, 0.0 },
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

  registerWrapper( "prescribedBcTable", &m_prescribedBcTable ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flag for whether to have time-dependent boundary condition types" );

  registerWrapper( "boundaryConditionTypes", &m_boundaryConditionTypes ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Boundary conditions on x-, x+, y-, y+, z- and z+ faces. Options are:\n* " + EnumStrings< BoundaryConditionOption >::concat( "\n* " ) );

  registerWrapper( "bcTable", &m_bcTable ).
    setInputFlag( InputFlags::FALSE ).
    setDescription( "Array that stores time-dependent bc types on x-, x+, y-, y+, z- and z+ faces." );

  registerWrapper( "prescribedBoundaryFTable", &m_prescribedBoundaryFTable ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flag for whether to have time-dependent boundary conditions described by a global background grid F" );

  registerWrapper( "fTableInterpType", &m_fTableInterpType ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "The type of F table interpolation. Options are 0 (linear), 1 (cosine), 2 (quintic polynomial)." );

  registerWrapper( "fTable", &m_fTable ).
    setInputFlag( InputFlags::FALSE ).
    setDescription( "Array that stores time-dependent grid-aligned stretches interpreted as a global background grid F." );

  registerWrapper( "needsNeighborList", &m_needsNeighborList ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0 ).
    setDescription( "Flag for whether to construct neighbor list" );

  registerWrapper( "neighborRadius", &m_neighborRadius ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( -1.0 ).
    setDescription( "Neighbor radius for SPH-type calculations" );

  registerWrapper( "binSizeMultiplier", &m_binSizeMultiplier ).
    setInputFlag( InputFlags::FALSE ).
    setApplyDefaultValue( 1 ).
    setDescription( "Multiplier for setting bin size, used to speed up particle neighbor sorting" );

  registerWrapper( "useDamageAsSurfaceFlag", &m_useDamageAsSurfaceFlag ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0 ).
    setDescription( "Indicates whether particle damage at the beginning of the simulation should be interpreted as a surface flag" );

  registerWrapper( "cpdiDomainScaling", &m_cpdiDomainScaling ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDefaultValue( 0 ).
    setDescription( "Option for CPDI domain scaling" );

  registerWrapper( "smallMass", &m_smallMass ).
    setInputFlag( InputFlags::FALSE ).
    setApplyDefaultValue( DBL_MAX ).
    setDescription( "The small mass threshold for ignoring extremely low-mass nodes." );

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
    setApplyDefaultValue( 0 ).
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
    setDescription( "Global minimum grid extent excluding buffer cells" );

  registerWrapper( "xGlobalMax", &m_xGlobalMax ).
    setInputFlag( InputFlags::FALSE ).
    setDescription( "Global maximum grid extent excluding buffer cells" );

  registerWrapper( "domainLengths", &m_partitionExtent ).
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

  // Activate neighbor list if necessary
  if( m_damageFieldPartitioning == 1 )
  {
    m_needsNeighborList = 1;
  }

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
  using namespace fields::mpm;
  
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
                                                    arrayView1d< string const > const & regionNames )
  {
    MeshBody const & meshBody = meshBodies.getGroup< MeshBody >( meshBodyName );
    if( meshBody.hasParticles() ) // Particle field registration? TODO: What goes here?
    {
      ParticleManager & particleManager = meshLevel.getParticleManager();
      
      particleManager.forParticleSubRegions< ParticleSubRegion >( regionNames,
                                                                  [&]( localIndex const,
                                                                       ParticleSubRegion & subRegion )
      {
        // Registration automatically allocates the 1st dimension (i.e. particle index) of the associated arrays.
        // Vector/tensor fields need to have their other dimensions specified.
        
        string const voightLabels[6] = { "XX", "YY", "ZZ", "YZ", "XZ", "XY" };

        // Single-indexed fields (scalars)
        subRegion.registerField< isBad >( getName() );
        subRegion.registerField< particleMass >( getName() );
        subRegion.registerField< particleInitialVolume >( getName() );
        subRegion.registerField< particleDensity >( getName() );
        //subRegion.registerField< particleDamage >( getName() ); // This is an intrinsic field for now since it comes from ParticleMeshGenerator
        subRegion.registerField< particleSurfaceFlag >( getName() );

        // Double-indexed fields (vectors and symmetric tensors stored in Voigt notation)
        subRegion.registerField< particleStress >( getName() ).setDimLabels( 1, voightLabels ).reference().resizeDimension< 1 >( 6 );
        subRegion.registerField< particleDamageGradient >( getName() ).reference().resizeDimension< 1 >( 3 );

        // Triple-indexed fields (vectors of vectors, non-symmetric tensors)
        subRegion.registerField< particleInitialRVectors >( getName() ).reference().resizeDimension< 1, 2 >( 3, 3 );
        subRegion.registerField< particleDeformationGradient >( getName() ).reference().resizeDimension< 1, 2 >( 3, 3 );
      } );
    }
    else // Background grid field registration
    {
      NodeManager & nodes = meshLevel.getNodeManager();

      nodes.registerWrapper< array2d< real64 > >( viewKeyStruct::massString() ).
        setPlotLevel( PlotLevel::LEVEL_0 ).
        setRegisteringObjects( this->getName() ).
        setDescription( "An array that holds the mass on the nodes." );
      
      nodes.registerWrapper< array3d< real64 > >( viewKeyStruct::velocityString() ).
        setPlotLevel( PlotLevel::LEVEL_0 ).
        setRegisteringObjects( this->getName() ).
        setDescription( "An array that holds the current velocity on the nodes." );

      nodes.registerWrapper< array3d< real64 > >( viewKeyStruct::momentumString() ).
        setPlotLevel( PlotLevel::LEVEL_1 ).
        setRegisteringObjects( this->getName() ).
        setDescription( "An array that holds the current momentum on the nodes." );

      nodes.registerWrapper< array3d< real64 > >( viewKeyStruct::accelerationString() ).
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
        setRestartFlags( RestartFlags::WRITE_AND_READ );

      nodeSets.registerWrapper< array1d< SortedArray< localIndex > > >( viewKeyStruct::bufferNodesString() ).
        setPlotLevel( PlotLevel::NOPLOT ).
        setRestartFlags( RestartFlags::WRITE_AND_READ );
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
      } );
    }
  } );

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
  // Initialize neighbor lists
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    OrderedVariableToManyParticleRelation & neighborList = subRegion.neighborList();
    neighborList.setParticleManager( particleManager );
  } );
  
  // Read and distribute BC table
  if( m_prescribedBcTable == 1 )
  {
    int rank = MpiWrapper::commRank( MPI_COMM_GEOSX );
    int BCTableSize;
    std::vector<double> BCTable1D; // Need 1D version of BC table for MPI broadcast

    if( rank == 0 ) // Rank 0 process parses the BC table file
    {
      std::ifstream fp( "BCTable.dat" );
      double BCTableEntry = 0;
      while( fp >> BCTableEntry )
      {
        BCTable1D.push_back( BCTableEntry ); // Push values into BCTable1D
      }
      BCTableSize = BCTable1D.size();
    }

    MPI_Bcast( &BCTableSize, 1, MPI_INT, 0, MPI_COMM_GEOSX ); // Broadcast the size of BCTable1D to other processes
    if( rank != 0 ) // All processes except for root resize their versions of BCTable1D
    {
      BCTable1D.resize( BCTableSize );
    }
    MPI_Bcast( BCTable1D.data(), BCTableSize, MPI_DOUBLE, 0, MPI_COMM_GEOSX ); // Broadcast BCTable1D to other processes

    // Technically don't need to reshape BCTable1D into a 2D array, but it makes things more readable and should have little runtime penalty
    m_bcTable.resize( BCTableSize / 7, 7 ); // Initialize size of m_BCTable
    for( int i = 0 ; i < BCTableSize ; i++ ) // Populate m_BCTable
    {
      m_bcTable[i / 7][i % 7] = BCTable1D[i]; // Taking advantage of integer division rounding towards zero
    }
  }

  // Initialize domain F and L, then read and distribute F table 
  m_domainF.resize(3);
  m_domainL.resize(3);
  for( int i=0; i<3; i++ )
  {
    m_domainF[i] = 1.0;
    m_domainL[i] = 0.0;
  }
  if( m_prescribedBoundaryFTable == 1 )
  {
    int rank = MpiWrapper::commRank( MPI_COMM_GEOSX );
    int FTableSize;
    std::vector<double> FTable1D; // Need 1D version of F table for MPI broadcast

    if( rank == 0 ) // Rank 0 process parses the F table file
    {
      std::ifstream fp( "FTable.dat" );
      double FTableEntry = 0;
      while( fp >> FTableEntry )
      {
        FTable1D.push_back( FTableEntry ); // Push values into FTable1D
      }
      FTableSize = FTable1D.size();
    }

    MPI_Bcast( &FTableSize, 1, MPI_INT, 0, MPI_COMM_GEOSX ); // Broadcast the size of FTable1D to other processes
    if( rank != 0 ) // All processes except for root resize their versions of FTable1D
    {
      FTable1D.resize( FTableSize );
    }
    MPI_Bcast( FTable1D.data(), FTableSize, MPI_DOUBLE, 0, MPI_COMM_GEOSX ); // Broadcast FTable1D to other processes

    // Techinically don't need to reshape FTable1D into a 2D array, but it makes things more readable and should have little runtime penalty
    m_fTable.resize( FTableSize / 4, 4 ); // Initialize size of m_fTable
    for( int i = 0 ; i < FTableSize ; i++ ) // Populate m_fTable
    {
      m_fTable[i / 4][i % 4] = FTable1D[i]; // Taking advantage of integer division rounding towards zero
    }
  }
  
  // Get grid quantites
  int numNodes = nodeManager.size();
  arrayView2d< real64, nodes::REFERENCE_POSITION_USD > const gridPosition = nodeManager.referencePosition();
  array2d< real64 > & gridMass = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::massString() );
  array2d< real64 > & gridDamage = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::damageString() );
  array2d< real64 > & gridMaxDamage = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::maxDamageString() );
  array3d< real64 > & gridVelocity = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::velocityString() );
  array3d< real64 > & gridMomentum = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::momentumString() );
  array3d< real64 > & gridAcceleration = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::accelerationString() );
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
      m_xLocalMin[i] = std::fmin( m_xLocalMin[i], gridPosition[g][i] );
      m_xLocalMax[i] = std::fmax( m_xLocalMax[i], gridPosition[g][i] );
    }
  }
  for(int i=0; i<3; i++)
  {
    m_xLocalMinNoGhost[i] = partition.getLocalMin()[i];
    m_xLocalMaxNoGhost[i] = partition.getLocalMax()[i];
    m_partitionExtent[i] = m_xLocalMax[i] - m_xLocalMin[i];
  }

  // Get element size
  for(int g=0; g<numNodes; g++)
  {
    for(int i=0; i<3; i++)
    {
      real64 test = gridPosition[g][i] - m_xLocalMin[i]; // By definition, this should always be positive. Furthermore, the gridPosition should only be those on the local partition
      if(test > 0.0) // We're looking for the smallest nonzero distance from the "min" node. TODO: Could be vulnerable to a finite precision bug.
      {
        m_hEl[i] = std::fmin( test, m_hEl[i] );
      }
    }
  }

  // Set SPH neighbor radius if necessary
  if( m_neighborRadius <= 0.0 )
  {
    if( m_planeStrain )
    {
      m_neighborRadius *= -1.0 * fmin( m_hEl[0], m_hEl[1] );
    }
    else
    {
      m_neighborRadius *= -1.0 * fmin( m_hEl[0], fmin( m_hEl[1], m_hEl[2] ) );
    }
  }

  // Get global domain extent excluding buffer nodes
  for(int i=0; i<3; i++)
  {
    m_xGlobalMin[i] = partition.getGlobalMin()[i] + m_hEl[i];
    m_xGlobalMax[i] = partition.getGlobalMax()[i] - m_hEl[i];
    m_domainExtent[i] = m_xGlobalMax[i] - m_xGlobalMin[i];
  }

  // Get number of elements in each direction
  for(int i=0; i<3; i++)
  {
    m_nEl[i] = std::round( m_partitionExtent[i] / m_hEl[i] );
  }

  // Create element map
  m_ijkMap.resize( m_nEl[0] + 1, m_nEl[1] + 1, m_nEl[2] + 1 );
  for( int g=0 ; g<numNodes ; g++ )
  {
    int i = std::round( ( gridPosition[g][0] - m_xLocalMin[0] ) / m_hEl[0] ) ;
    int j = std::round( ( gridPosition[g][1] - m_xLocalMin[1] ) / m_hEl[1] ) ;
    int k = std::round( ( gridPosition[g][2] - m_xLocalMin[2] ) / m_hEl[2] ) ;
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

    int dir0 = face / 2; // 0, 0, 1, 1, 2, 2 (x-, x+, y-, y+, z-, z+)

    int positiveNormal = face % 2; // even => (-) => 0, odd => (+) => 1
    int sign = 2 * positiveNormal - 1;
    real64 minOrMax = positiveNormal == 0 ? m_xGlobalMin[dir0] : m_xGlobalMax[dir0];

    for( localIndex g=0; g<numNodes; g++ )
    {
      real64 positionRelativeToBoundary = gridPosition[g][dir0] - minOrMax;
      real64 tolerance = 1.0e-12 * m_hEl[dir0]; // small multiple of element dimension

      if( std::fabs( positionRelativeToBoundary ) < tolerance )
      {
        tmpBoundaryNodes.insert(g);
      }

      if( sign * positionRelativeToBoundary > 0 && std::fabs( positionRelativeToBoundary ) > tolerance ) // basically a dot product with the face normal
      {
        tmpBufferNodes.insert(g);
      }
    }

    m_boundaryNodes[face].insert( tmpBoundaryNodes.begin(), tmpBoundaryNodes.end() );
    m_bufferNodes[face].insert( tmpBufferNodes.begin(), tmpBufferNodes.end() );
  }

  // Initialize reaction force history file and write its header
  if( MpiWrapper::commRank( MPI_COMM_GEOSX ) == 0 )
  {
    std::ofstream file;
    file.open( "reactionHistory.csv", std::ios::out | std::ios::app );
    if( file.fail() )
      throw std::ios_base::failure( std::strerror( errno ) );
    //make sure write fails with exception if something is wrong
    file.exceptions( file.exceptions() | std::ios::failbit | std::ifstream::badbit );
    file << "time, F00, F11, F22, length_x, length_y, length_z, Rx-, Rx+, Ry-, Ry+, Rz-, Rz+, L00, L11, L22" << std::endl;
    file << std::setprecision( std::numeric_limits<long double>::digits10 )
         << 0.0 << ","
         << 1.0 << "," << 1.0 << "," << 1.0 << ","
         << m_domainExtent[0] << "," << m_domainExtent[1] << "," << m_domainExtent[2] << ","
         << 0.0 << "," << 0.0 << "," << 0.0 << "," << 0.0 << "," << 0.0 << "," << 0.0 << ","
         << 0.0 << "," << 0.0 << "," << 0.0
         << std::endl;
  }

  // Set particle masses based on their volume and density. Set deformation gradient to identity;
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    // Getters
    string const & solidMaterialName = subRegion.template getReference< string >( viewKeyStruct::solidMaterialNamesString() );
    SolidBase & constitutiveRelation = getConstitutiveModel< SolidBase >( subRegion, solidMaterialName ); // For the time being we restrict our attention to elastic isotropic solids. TODO: Have all constitutive models automatically calculate a wave speed.
    arrayView2d< real64 > const constitutiveDensity = constitutiveRelation.getDensity();
    arrayView1d< real64 > const particleDensity = subRegion.getField< fields::mpm::particleDensity >();
    arrayView1d< real64 > const particleVolume = subRegion.getParticleVolume();
    arrayView1d< real64 > const particleDamage = subRegion.getParticleDamage();
    arrayView3d< real64 > const particleRVectors = subRegion.getParticleRVectors();
    arrayView1d< real64 > const particleMass = subRegion.getField< fields::mpm::particleMass >();
    arrayView3d< real64 > const particleDeformationGradient = subRegion.getField< fields::mpm::particleDeformationGradient >();
    arrayView1d< real64 > const particleInitialVolume = subRegion.getField< fields::mpm::particleInitialVolume >();
    arrayView3d< real64 > const particleInitialRVectors = subRegion.getField< fields::mpm::particleInitialRVectors >();
    arrayView1d< int > const particleSurfaceFlag = subRegion.getField< fields::mpm::particleSurfaceFlag >();

    // Set initial volume and R-vectors
    for(int p=0; p<subRegion.size(); p++)
    {
      particleInitialVolume[p] = particleVolume[p];
      for(int i=0; i<3; i++)
      {
        for(int j=0; j<3; j++)
        {
          particleInitialRVectors[p][i][j] = particleRVectors[p][i][j];
        }
      }
    }

    // Pull density from constitutive model, set particle masses and small mass threshold
    real64 localMinMass = 0.0;
    real64 globalMinMass;
    for(int p=0; p<subRegion.size(); p++)
    {
      particleDensity[p] = constitutiveDensity[p][0];
      particleMass[p] = particleDensity[p] * particleVolume[p];
      localMinMass = particleMass[p] < localMinMass ? particleMass[p] : localMinMass;
    }
    if( subRegion.size() == 0 ) // Handle empty partitions
    {
      localMinMass = DBL_MAX;
    }
    MPI_Allreduce( &localMinMass,
                   &globalMinMass,
                   1,
                   MPI_DOUBLE,
                   MPI_MIN,
                   MPI_COMM_GEOSX );
    m_smallMass = fmin( globalMinMass * 1.0e-12, m_smallMass );

    // Initialize deformation gradient
    for(int p=0; p<subRegion.size(); p++)
    {
      tOps::addIdentity< 3 >( particleDeformationGradient[p], 1.0 );
    }

    // Set surface flags
    for(int p=0; p<subRegion.size(); p++)
    {
      particleSurfaceFlag[p] = 0;
      if( particleDamage[p] > 0.0 && m_useDamageAsSurfaceFlag == 1 )
      {
        particleDamage[p] = 0.0;
        particleSurfaceFlag[p] = 2;
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
                 MPI_COMM_GEOSX );

  // Number of contact groups
  m_numContactGroups = maxGlobalGroupNumber + 1;

  // Specified number of damage flags.
  m_numContactFlags = m_damageFieldPartitioning == 1 ? 2 : 1;

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

real64 SolidMechanicsMPM::explicitStep( real64 const & time_n,
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
  arrayView2d< real64, nodes::REFERENCE_POSITION_USD > const gridPosition = nodeManager.referencePosition();
  arrayView1d< int > const gridGhostRank = nodeManager.ghostRank();
  array2d< real64 > & gridMass = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::massString() );
  array2d< real64 > & gridDamage = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::damageString() );
  array2d< real64 > & gridMaxDamage = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::maxDamageString() );
  array3d< real64 > & gridVelocity = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::velocityString() );
  array3d< real64 > & gridMomentum = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::momentumString() );
  array3d< real64 > & gridAcceleration = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::accelerationString() );
  array3d< real64 > & gridInternalForce = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::forceInternalString() );
  array3d< real64 > & gridExternalForce = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::forceExternalString() );
  array3d< real64 > & gridContactForce = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::forceContactString() );
  array3d< real64 > & gridSurfaceNormal = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::surfaceNormalString() );
  array3d< real64 > & gridMaterialPosition = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::materialPositionString() );


  //#######################################################################################
  solverProfilingIf( "At time step zero, perform initialization calculations", cycleNumber == 0 );
  //#######################################################################################
  if( cycleNumber == 0 )
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
  solverProfiling( "Update global-to-local map" );
  //#######################################################################################
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    subRegion.updateMaps();
  } );


  //#######################################################################################
  solverProfilingIf( "Perform particle ghosting", MpiWrapper::commSize( MPI_COMM_GEOSX ) > 1 && m_needsNeighborList == 1 );
  //#######################################################################################
  if( MpiWrapper::commSize( MPI_COMM_GEOSX ) > 1 && m_needsNeighborList == 1 )
  {
    partition.getGhostParticlesFromNeighboringPartitions( domain, m_iComm, m_neighborRadius );
  }


  //#######################################################################################
  solverProfiling( "Get indices of non-ghost particles" );
  //#######################################################################################
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    subRegion.setNonGhostIndices();
  } );


  //#######################################################################################
  solverProfilingIf( "Construct neighbor list", m_needsNeighborList == 1 );
  //#######################################################################################
  if( m_needsNeighborList == 1 )
  {
    if( cycleNumber == 0)
    {
      optimizeBinSort( particleManager );
    }
    else
    {
      (void) computeNeighborList( particleManager );
    }
  }


  //#######################################################################################
  solverProfilingIf( "Compute damage field gradient", m_damageFieldPartitioning == 1 );
  //#######################################################################################
  if( m_damageFieldPartitioning == 1 )
  {
    // Get volume, position, damage, surface flag
    ParticleManager::ParticleViewAccessor< arrayView1d< real64 const > > particleVolumeAccessor = particleManager.constructArrayViewAccessor< real64, 1 >( "particleVolume" );
    ParticleManager::ParticleViewAccessor< arrayView2d< real64 const > > particlePositionAccessor = particleManager.constructArrayViewAccessor< real64, 2 >( "particleCenter" );
    ParticleManager::ParticleViewAccessor< arrayView1d< real64 const > > particleDamageAccessor = particleManager.constructArrayViewAccessor< real64, 1 >( "particleDamage" );
    ParticleManager::ParticleViewAccessor< arrayView1d< int const > > particleSurfaceFlagAccessor = particleManager.constructFieldAccessor< fields::mpm::particleSurfaceFlag >();

    // Perform neighbor operations
    particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
    {
      // Get neighbor list
      OrderedVariableToManyParticleRelation & neighborList = subRegion.neighborList();
      arrayView1d< localIndex const > numNeighborsAll = neighborList.m_numParticles.toViewConst();
      ArrayOfArraysView< localIndex const > neighborRegions = neighborList.m_toParticleRegion.toViewConst();
      ArrayOfArraysView< localIndex const > neighborSubRegions = neighborList.m_toParticleSubRegion.toViewConst();
      ArrayOfArraysView< localIndex const > neighborIndices = neighborList.m_toParticleIndex.toViewConst();

      // Get particle position and damage field gradient
      arrayView2d< real64 > const particlePosition = subRegion.getParticleCenter();
      arrayView2d< real64 > const particleDamageGradient = subRegion.getField< fields::mpm::particleDamageGradient >();

      // Declare arrays for holding neighbor data
      array1d< real64 > neighborVolumes;
      array2d< real64 > neighborPositions;
      array1d< real64 > neighborDamages;

      // Loop over neighbors
      for( localIndex& p: subRegion.nonGhostIndices() )
      {
        // Get number of neighbors and accessor indices
        localIndex numNeighbors = numNeighborsAll[p];
        arraySlice1d< localIndex const > regionIndices = neighborRegions[p];
        arraySlice1d< localIndex const > subRegionIndices = neighborSubRegions[p];
        arraySlice1d< localIndex const > particleIndices = neighborIndices[p];

        // Size neighbor data arrays
        neighborVolumes.resize(numNeighbors);
        neighborPositions.resize(numNeighbors,3);
        neighborDamages.resize(numNeighbors);

        // Populate neighbor data arrays
        for( localIndex neighborIndex = 0; neighborIndex < numNeighbors; neighborIndex++ )
        {
          localIndex regionIndex = regionIndices[neighborIndex];
          localIndex subRegionIndex = subRegionIndices[neighborIndex];
          localIndex particleIndex = particleIndices[neighborIndex];
          neighborVolumes[neighborIndex] = particleVolumeAccessor[regionIndex][subRegionIndex][particleIndex];
          neighborPositions[neighborIndex][0] = particlePositionAccessor[regionIndex][subRegionIndex][particleIndex][0];
          neighborPositions[neighborIndex][1] = particlePositionAccessor[regionIndex][subRegionIndex][particleIndex][1];
          neighborPositions[neighborIndex][2] = particlePositionAccessor[regionIndex][subRegionIndex][particleIndex][2];
          if( particleSurfaceFlagAccessor[regionIndex][subRegionIndex][particleIndex] == 2 )
          {
            neighborDamages[neighborIndex] = 1.0;
          }
          else
          {
            neighborDamages[neighborIndex] = particleDamageAccessor[regionIndex][subRegionIndex][particleIndex];
          }
        }

        // Call kernel field gradient function
        computeKernelFieldGradient( particlePosition[p],          // input
                                    neighborPositions,            // input
                                    neighborVolumes,              // input
                                    neighborDamages,              // input
                                    particleDamageGradient[p] );  // OUTPUT
      }
    } );
  }


  //#######################################################################################
  solverProfilingIf( "Update BCs based on bcTable", m_prescribedBcTable == 1 );
  //#######################################################################################
  if( m_prescribedBcTable == 1)
  {
    int bcInterval;

    for( localIndex i = 0 ; i < m_bcTable.size(0) ; i++ ) // Naive method for determining what part of BC table we're currently in, can optimize later (TODO)
    {
      if( time_n + 0.5 * dt > m_bcTable[i][0] )
      {
        bcInterval = i;
      }
    }

    for( int i=0; i<6; i++ )
    {
      m_boundaryConditionTypes[i] = m_bcTable[bcInterval][i+1];
    }
  }


  //#######################################################################################
  solverProfilingIf( "Perform r-vector scaling (CPDI domain scaling)", m_cpdiDomainScaling == 1 );
  //#######################################################################################
  if( m_cpdiDomainScaling == 1 )
  {
    particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
    {
      if( subRegion.getParticleType() == ParticleType::CPDI )
      {
        real64 lCrit = m_planeStrain ? 0.49999 * fmin( m_hEl[0], m_hEl[1] ) : 0.49999 * fmin( m_hEl[0], fmin( m_hEl[1], m_hEl[2] ) );
        subRegion.cpdiDomainScaling( lCrit, m_planeStrain );
      }
    } );
  }


  //#######################################################################################
  solverProfiling( "Particle-to-grid interpolation" );
  //#######################################################################################
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    // subregion fields
    arrayView2d< real64 > const particleCenter = subRegion.getParticleCenter();
    arrayView2d< real64 > const particleVelocity = subRegion.getParticleVelocity();
    arrayView1d< real64 > const particleMass = subRegion.getField< fields::mpm::particleMass >();
    arrayView1d< real64 > const particleVolume = subRegion.getParticleVolume();
    arrayView1d< int > const particleGroup = subRegion.getParticleGroup();

    // solver fields
    arrayView2d< real64 > const particleStress = subRegion.getField< fields::mpm::particleStress >();

    // initialize mapping helpers
    int numNodesMappedTo = subRegion.numNodesMappedTo();
    array1d< int > nodeIDs(numNodesMappedTo); // nodes that the particle maps to
    array1d< real64 > weights(numNodesMappedTo); // shape function value for each node
    array2d< real64 > gradWeights(numNodesMappedTo, 3); // shape function gradient value for each node; 1st index = direction, 2nd index = node

    for( localIndex& p: subRegion.nonGhostIndices() )
    {
      arraySlice1d< real64 > const p_x = particleCenter[p]; // auto = LvArray::ArraySlice<double, 1, 0, int>
      arraySlice1d< real64 > const p_v = particleVelocity[p]; // auto = LvArray::ArraySlice<double, 1, 0, int>
      real64 const & p_m = particleMass[p];
      real64 const & p_vol = particleVolume[p];
      arraySlice1d< real64 > const p_stress = particleStress[p];
      int const & p_group = particleGroup[p];

      // Get interpolation kernel
      subRegion.getAllWeights( p,
                               m_xLocalMin,
                               m_hEl,
                               m_ijkMap,
                               gridPosition,
                               nodeIDs,       // output
                               weights,       // output
                               gradWeights ); // output

      // Map to grid
      for(int g=0; g<numNodesMappedTo; g++)
      {
        int mappedNode = nodeIDs[g];
        int fieldIndex = p_group;
        gridMass[mappedNode][fieldIndex] += p_m * weights[g];
        for(int i=0; i<m_numDims; i++)
        {
          gridMomentum[mappedNode][fieldIndex][i] += p_m * p_v[i] * weights[g];
          gridMaterialPosition[mappedNode][fieldIndex][i] += p_m * (p_x[i] - gridPosition[mappedNode][i]) * weights[g]; // TODO: Switch to volume weighting?
          for(int k=0; k<3; k++)
          {
            int voigt = m_voigtMap[k][i];
            gridInternalForce[mappedNode][fieldIndex][i] -= p_stress[voigt] * gradWeights[g][k] * p_vol;
          }
        }
      }

    } // particle loop
  } ); // subregion loop


  //#######################################################################################
  solverProfiling( "Grid MPI operations" );
  //#######################################################################################
  std::vector< std::string > fieldNames = { viewKeyStruct::massString(),
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
        for( int i=0; i<m_numDims; i++ )
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
        for( int i=0; i<m_numDims; i++ )
        {
          gridAcceleration[g][fieldIndex][i] = 0.0;
          gridVelocity[g][fieldIndex][i] = 0.0;
          gridMomentum[g][fieldIndex][i] = 0.0;
          gridMaterialPosition[g][fieldIndex][i] = gridPosition[g][i];
        }
      }
    }
  }


  //#######################################################################################
  solverProfilingIf( "Contact enforcement", m_numVelocityFields > 1 );
  //#######################################################################################
  if( m_numVelocityFields > 1 )
  {
    // Compute grid surface normals
    computeGridSurfaceNormals( particleManager, gridPosition, gridSurfaceNormal );

    // Sync surface normals
    syncGridFields( { viewKeyStruct::surfaceNormalString() }, domain, nodeManager, mesh );

    // Apply symmetry boundary conditions to surface normals
    enforceGridVectorFieldSymmetryBC( gridSurfaceNormal, gridPosition, nodeManager.sets() );

    // Normalize grid surface normals
    normalizeGridSurfaceNormals( gridMass, gridSurfaceNormal );

    // TODO: Enforce symmetry on grid mass and material position?

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
          for( int i=0; i<m_numDims; i++ )
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


  //#######################################################################################
  solverProfilingIf( "Interpolate F table", m_prescribedBoundaryFTable == 1 );
  //#######################################################################################
  if( m_prescribedBoundaryFTable == 1 )
  {
    double Fii_new;
    double Fii_dot;
    int fInterval;
    double timeInterval;
    double timePast;

    if( time_n + dt < m_fTable[m_fTable.size(0) - 1][0] ) // If within F table bounds...
    {
      for( localIndex i = 0 ; i < m_fTable.size(0) ; i++ ) // Naive method for determining what part of F table we're currently in, can optimize later (TODO)
      {
        if( time_n + dt / 2 > m_fTable[i][0] )
        {
          fInterval = i;
        }
      }

      timeInterval = m_fTable[fInterval + 1][0] - m_fTable[fInterval][0]; // Time fInterval for current part of F table we're in
      timePast = time_n + dt - m_fTable[fInterval][0]; // Time elapsed since switching intervals in F table

      for( int i = 0 ; i < m_numDims ; i++ ) // Update L and F
      {
        // smooth-step interpolation with cosine, zero endpoint velocity
        if( m_fTableInterpType == 1 )
        {
          Fii_new = m_fTable[fInterval][i + 1] - 0.5 * ( m_fTable[fInterval + 1][i + 1] - m_fTable[fInterval][i + 1] ) * ( cos(
              3.141592653589793 * timePast / timeInterval ) - 1.0 );
        }
        // smooth-step interpolation with 5th order polynomial, zero endpoint velocity and acceleration
        else if( m_fTableInterpType == 2 )
        {
          Fii_new = m_fTable[fInterval][i+1] + (m_fTable[fInterval+1][i+1] - m_fTable[fInterval][i+1])*(10.0*pow(timePast/timeInterval,3) - 15.0*pow(timePast/timeInterval,4) + 6.0*pow(timePast/timeInterval,5));
        }
        // default linear interpolation
        else
        {
          Fii_new = m_fTable[fInterval][i + 1] * ( timeInterval - timePast ) / timeInterval + m_fTable[fInterval + 1][i + 1] * ( timePast ) / timeInterval;
        }

        // // drive F22 and F33 using F11 to get pure shear
        // if(m_fTable_pure_shear && i > 0)
        // {
        //   Fii_new = m_planeStrain ? 1.0/m_domainF[0] : 1.0/sqrt(m_domainF[0]);
        // }

        Fii_dot = ( Fii_new - m_domainF[i] ) / dt;
        m_domainL[i] = Fii_dot / Fii_new; // L = Fdot.Finv
        m_domainF[i] = Fii_new;
      }
    }
    else if( time_n + dt >= m_fTable[m_fTable.size(0) - 1][0] ) // Else (i.e. if we exceed F table upper bound)...
    {
      for( int i = 0 ; i < m_numDims ; i++ )
      {
        Fii_new = m_fTable[m_fTable.size(0) - 1][i + 1]; // Set Fii_new to the last prescribed F table value
        Fii_dot = ( Fii_new - m_domainF[i] ) / dt;
        m_domainL[i] = Fii_dot / Fii_new; // L = Fdot.Finv
        m_domainF[i] = Fii_new;
      }
    }
  }

  //#######################################################################################
  solverProfiling( "Apply essential boundary conditions" );
  //#######################################################################################
  applyEssentialBCs( dt, time_n, gridVelocity, gridAcceleration, gridMass, gridGhostRank, gridPosition, nodeManager.sets() );


  //#######################################################################################
  solverProfiling( "Grid-to-particle interpolation" );
  //#######################################################################################
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    // Registered by subregion
    arrayView2d< real64 > const particleCenter = subRegion.getParticleCenter();
    arrayView2d< real64 > const particleVelocity = subRegion.getParticleVelocity();
    arrayView1d< real64 > const particleVolume = subRegion.getParticleVolume();
    arrayView1d< int > const particleGroup = subRegion.getParticleGroup();
    arrayView1d< globalIndex > particleID = subRegion.getParticleID();

    // Registered by constitutive model
    string const & solidMaterialName = subRegion.template getReference< string >( viewKeyStruct::solidMaterialNamesString() );
    ElasticIsotropic & constitutiveRelation = getConstitutiveModel< ElasticIsotropic >( subRegion, solidMaterialName ); // again, limiting to elastic isotropic for now
    arrayView1d< real64 > const shearModulus = constitutiveRelation.shearModulus();
    arrayView1d< real64 > const bulkModulus = constitutiveRelation.bulkModulus();

    // Registered by MPM solver
    arrayView1d< int > const isBad = subRegion.getField< fields::mpm::isBad >();
    arrayView1d< real64 > const particleMass = subRegion.getField< fields::mpm::particleMass >();
    arrayView1d< real64 > const particleInitialVolume = subRegion.getField< fields::mpm::particleInitialVolume >();
    arrayView3d< real64 > const particleDeformationGradient = subRegion.getField< fields::mpm::particleDeformationGradient >();
    arrayView3d< real64 > const particleInitialRVectors = subRegion.getField< fields::mpm::particleInitialRVectors >();
    arrayView1d< real64 > const particleDensity = subRegion.getField< fields::mpm::particleDensity >();
    arrayView2d< real64 > const particleStress = subRegion.getField< fields::mpm::particleStress >();

    // Initialize mapping helpers
    int numNodesMappedTo = subRegion.numNodesMappedTo();
    array1d< int > nodeIDs(numNodesMappedTo); // nodes that the particle maps to
    array1d< real64 > weights(numNodesMappedTo); // shape function value for each node
    array2d< real64 > gradWeights(numNodesMappedTo, 3); // shape function gradient value for each node; 1st index = direction, 2nd index = node

    // Particle loop - we might be able to get rid of this someday and have everything happen via MPMParticleSubRegion methods
    for( localIndex& p: subRegion.nonGhostIndices() )
    {
      arraySlice1d< real64 > const p_x = particleCenter[p];
      arraySlice1d< real64 > const p_v = particleVelocity[p];
      real64 & p_m = particleMass[p];
      real64 & p_vol = particleVolume[p];
      real64 const & p_vol0 = particleInitialVolume[p];
      real64 & p_rho = particleDensity[p];
      arraySlice2d< real64 > const p_F = particleDeformationGradient[p]; // auto = LvArray::ArraySlice<double, 2, 1, int>
      arraySlice1d< real64 > const p_stress = particleStress[p];
      int const & p_group = particleGroup[p];
      arraySlice2d< real64 > const p_rvec0 = particleInitialRVectors[p];
      real64 p_L[3][3] = { {0} }; // Velocity gradient
      real64 p_D[3][3] = { {0} }; // Rate of deformation
      real64 p_Diso[3][3] = { {0} }; // Rate of deformation
      real64 p_Ddev[3][3] = { {0} }; // Rate of deformation
      real64 p_FOld[3][3] = { {0} }; // Old particle F
      real64 detF; // Material Jacobian

      // Store the old particle F
      for(int i=0; i<3; i++)
      {
        for(int j=0; j<3; j++)
        {
          p_FOld[i][j] = p_F[i][j];
        }
      }

      // Get interpolation kernel  - TODO: Seems dumb to have to do this twice
      subRegion.getAllWeights( p,
                               m_xLocalMin,
                               m_hEl,
                               m_ijkMap,
                               gridPosition,
                               nodeIDs,       // output
                               weights,       // output
                               gradWeights ); // output

      // Particle-to-grid map
      for(int g=0; g<numNodesMappedTo; g++)
      {
        int mappedNode = nodeIDs[g];
        int fieldIndex = p_group;
        for(int i=0; i<m_numDims; i++)
        {
          p_x[i] += gridVelocity[mappedNode][fieldIndex][i] * dt * weights[g];
          p_v[i] += gridAcceleration[mappedNode][fieldIndex][i] * dt * weights[g]; // FLIP
          for(int j=0; j<m_numDims; j++)
          {
            p_L[i][j] += gridVelocity[mappedNode][fieldIndex][i] * gradWeights[g][j]; // Technically wrong, the best kind of wrong (end-of-step velocity with beginning-of-step gradient)
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
      if( detF <= 0.1 || detF >= 10.0 )
      {
        isBad[p] = 1;
        GEOSX_LOG_RANK( "Flagging particle with unreasonable Jacobian (<0.1 or >10) for deletion! Global particle ID: " << particleID[p] );
        continue; // move on to the next particle since log(detF) will throw an error
      }
      p_vol = p_vol0 * detF;
      p_rho = p_m / p_vol;
      subRegion.computeRVectors( p, p_F, p_rvec0 );

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
        for(int j=i; j<3; j++)
        {
          int voigt = m_voigtMap[i][j];
          p_stress[voigt] = stressFull[i][j];
        }
      }

    } // particle loop
  } ); // subregion loop


  //#######################################################################################
  solverProfiling( "Calculate stable time step" );
  //#######################################################################################
  real64 wavespeed = 0.0;
  real64 length = m_planeStrain == 1 ? std::fmin( m_hEl[0], m_hEl[1] ) : std::fmin( m_hEl[0], std::fmin( m_hEl[1], m_hEl[2] ) );

  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    arrayView2d< real64 > const particleVelocity = subRegion.getParticleVelocity();
    string const & solidMaterialName = subRegion.template getReference< string >( viewKeyStruct::solidMaterialNamesString() );
    ElasticIsotropic & constitutiveRelation = getConstitutiveModel< ElasticIsotropic >( subRegion, solidMaterialName ); // For the time being we restrict our attention to elastic isotropic solids. TODO: Have all constitutive models automatically calculate a wave speed.
    arrayView1d< real64 > const rho = subRegion.getField< fields::mpm::particleDensity >();
    arrayView1d< real64 > const g = constitutiveRelation.shearModulus();
    arrayView1d< real64 > const k = constitutiveRelation.bulkModulus();
    for( localIndex& p: subRegion.nonGhostIndices() )
    {
      wavespeed = std::max( wavespeed, sqrt( ( k[p] + (4.0/3.0) * g[p] ) / rho[p] ) + tOps::l2Norm< 3 >( particleVelocity[p] ) );
    }
  } );

  real64 dtReturn = wavespeed > 1.0e-12 ? m_cflFactor * length / wavespeed : DBL_MAX; // This partitions's dt, make it huge if wavespeed=0.0 (this happens when there are no particles on this partition)


  //#######################################################################################
  solverProfiling( "Delete bad particles" );
  //#######################################################################################
  // Cases covered:
  // 1.) Particles that map outside the domain (including buffer cells)
  // 2.) Particles with unacceptable Jacobian (<0.1 or >10)
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    // Get relevant particle arrays
    arrayView1d< int > const isBad = subRegion.getField< fields::mpm::isBad >();
    arrayView1d< int > const particleRank = subRegion.getParticleRank();

    // Initialize the set of particles to delete
    std::set< localIndex > indicesToErase;

    subRegion.flagOutOfRangeParticles( m_xGlobalMin, m_xGlobalMax, m_hEl, isBad ); // This skips ghost particles

    for( int p=subRegion.size()-1; p>=0; p-- ) // TODO: Looping over a set containing indices ordered largest->smallest would be more elegant and probably faster.
    {                                          //       We could also do away with the rank check since only master particles are ever evaluated for deletion.
      if( isBad[p] == 1 && particleRank[p] == MpiWrapper::commRank( MPI_COMM_GEOSX ) )
      {
        indicesToErase.insert(p);
      }
      subRegion.erase(indicesToErase);
    }
  } );


  //#######################################################################################
  solverProfilingIf( "Particle repartitioning", MpiWrapper::commSize( MPI_COMM_GEOSX ) > 1 );
  //#######################################################################################
  if( MpiWrapper::commSize( MPI_COMM_GEOSX ) > 1 )
  {
    particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
    {
      partition.repartitionMasterParticles( subRegion, m_iComm );
    } );
  }


  //#######################################################################################
  solverProfilingIf( "Resize grid based on F-table", m_prescribedBoundaryFTable == 1 );
  //#######################################################################################
  if( m_prescribedBoundaryFTable == 1 )
  {
    resizeGrid( partition, gridPosition, dt );
  }


  //#######################################################################################
  solverProfiling( "End of explicitStep" );
  //#######################################################################################
  if( m_solverProfiling >= 1 )
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

  // Return stable time step
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
  for( size_t n=0; n<neighbors.size(); n++ )
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

void SolidMechanicsMPM::singleFaceVectorFieldSymmetryBC( const int face,
                                                         array3d< real64 > & vectorMultiField,
                                                         arrayView2d< real64, nodes::REFERENCE_POSITION_USD > const gridPosition,
                                                         Group & nodeSets )
{
  // This is a helper function for enforcing symmetry BCs on a single face and is meant to be called by other functions, not directly by the solver:
  //   * enforceGridVectorFieldSymmetryBC calls this on all grid faces
  //   * applyEssentialBCs calls this on faces that aren't moving (moving faces due to F-table need special treatment)
  
  array1d< SortedArray< localIndex > > & m_boundaryNodes = nodeSets.getReference< array1d< SortedArray< localIndex > > >( viewKeyStruct::boundaryNodesString() );
  array1d< SortedArray< localIndex > > & m_bufferNodes = nodeSets.getReference< array1d< SortedArray< localIndex > > >( viewKeyStruct::bufferNodesString() );

  for( int fieldIndex=0; fieldIndex<m_numVelocityFields; fieldIndex++ )
  {
    // Face-associated quantities
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
      ijk[dir1] = std::round( ( gridPosition[g][dir1] - m_xLocalMin[dir1] ) / m_hEl[dir1] );
      ijk[dir2] = std::round( ( gridPosition[g][dir2] - m_xLocalMin[dir2] ) / m_hEl[dir2] );

      localIndex gFrom = m_ijkMap[ijk[0]][ijk[1]][ijk[2]];

      vectorMultiField[g][fieldIndex][dir0] = -vectorMultiField[gFrom][fieldIndex][dir0]; // Negate component aligned with surface normal
      vectorMultiField[g][fieldIndex][dir1] =  vectorMultiField[gFrom][fieldIndex][dir1];
      vectorMultiField[g][fieldIndex][dir2] =  vectorMultiField[gFrom][fieldIndex][dir2];
    }
  }
}                                                          

void SolidMechanicsMPM::enforceGridVectorFieldSymmetryBC( array3d< real64 > & vectorMultiField,
                                                          arrayView2d< real64, nodes::REFERENCE_POSITION_USD > const gridPosition,
                                                          Group & nodeSets )
{
  for( int face=0; face<6; face++ )
  {
    if( m_boundaryConditionTypes[face] == 1 || m_boundaryConditionTypes[face] == 2 )
    {
      singleFaceVectorFieldSymmetryBC( face, vectorMultiField, gridPosition, nodeSets );
    }
  }
}

void SolidMechanicsMPM::applyEssentialBCs( const real64 dt,
                                           const real64 time_n,
                                           array3d<real64> & velocity,
                                           array3d<real64> & acceleration,
                                           array2d< real64 > & mass,
                                           arrayView1d< int > const ghostRank,
                                           arrayView2d<real64, nodes::REFERENCE_POSITION_USD> const gridPosition,
                                           Group & nodeSets )
{
  // Get node sets
  array1d< SortedArray< localIndex > > & m_boundaryNodes = nodeSets.getReference< array1d< SortedArray< localIndex > > >( viewKeyStruct::boundaryNodesString() );
  array1d< SortedArray< localIndex > > & m_bufferNodes = nodeSets.getReference< array1d< SortedArray< localIndex > > >( viewKeyStruct::bufferNodesString() );

  // Impose BCs on each face while gathering reaction forces
  real64 localFaceReactions[6] = {0.0};
  for( int face = 0; face < 6; face++ )
  {
    if( m_boundaryConditionTypes[face] == 1 )
    {
      singleFaceVectorFieldSymmetryBC( face, velocity, gridPosition, nodeSets );
      singleFaceVectorFieldSymmetryBC( face, acceleration, gridPosition, nodeSets );
    }
    else if( m_boundaryConditionTypes[face] == 2 && m_prescribedBoundaryFTable == 1 )
    {
      for( int fieldIndex = 0; fieldIndex < m_numVelocityFields; fieldIndex++ )
      {
        // Face-associated quantities
        int dir0 = face / 2;           // 0, 0, 1, 1, 2, 2 (x-, x+, y-, y+, z-, z+)
        int dir1 = (dir0 + 1) % 3;     // 1, 1, 2, 2, 0, 0
        int dir2 = (dir0 + 2) % 3;     // 2, 2, 0, 0, 1, 1
        int positiveNormal = face % 2; // even => (-) => 0, odd => (+) => 1

        // Enforce BCs on boundary nodes using F-table
        for( const auto &g : m_boundaryNodes[face] )
        {
          real64 prescribedVelocity = m_domainL[dir0] * gridPosition[g][dir0];
          real64 accelerationForBC = (prescribedVelocity - velocity[g][fieldIndex][dir0]) / dt; // acceleration needed to satisfy BC
          velocity[g][fieldIndex][dir0] = prescribedVelocity;
          acceleration[g][fieldIndex][dir0] += accelerationForBC;
          if( ghostRank[g] <= -1 ) // so we don't double count reactions at partition boundaries
          {
            localFaceReactions[face] += accelerationForBC * mass[g][fieldIndex];
          }
        }

        // Perform field reflection on buffer nodes - accounts for moving boundary effects
        for( const auto &g : m_bufferNodes[face] )
        {
          // Initialize grid ijk indices
          int ijk[3];
          ijk[dir1] = std::round((gridPosition[g][dir1] - m_xLocalMin[dir1]) / m_hEl[dir1]);
          ijk[dir2] = std::round((gridPosition[g][dir2] - m_xLocalMin[dir2]) / m_hEl[dir2]);

          // Grab the node index that we're copying from
          ijk[dir0] = positiveNormal * (m_nEl[dir0] - 2) + (1 - positiveNormal) * (2);
          localIndex gFrom = m_ijkMap[ijk[0]][ijk[1]][ijk[2]];

          // Grab the associated boundary node index for moving boundary correction
          ijk[dir0] = positiveNormal * (m_nEl[dir0] - 1) + (1 - positiveNormal) * (1);
          localIndex gBoundary = m_ijkMap[ijk[0]][ijk[1]][ijk[2]];

          // Calculate velocity
          velocity[g][fieldIndex][dir0] = -velocity[gFrom][fieldIndex][dir0] + 2.0 * velocity[gBoundary][fieldIndex][dir0]; // Negate component aligned with surface normal and correct for moving boundary
          velocity[g][fieldIndex][dir1] = velocity[gFrom][fieldIndex][dir1];
          velocity[g][fieldIndex][dir2] = velocity[gFrom][fieldIndex][dir2];

          // Calculate acceleration
          acceleration[g][fieldIndex][dir0] = -acceleration[gFrom][fieldIndex][dir0] + 2.0 * acceleration[gBoundary][fieldIndex][dir0]; // Negate component aligned with surface normal and correct for moving boundary
          acceleration[g][fieldIndex][dir1] = acceleration[gFrom][fieldIndex][dir1];
          acceleration[g][fieldIndex][dir2] = acceleration[gFrom][fieldIndex][dir2];
        }
      }
    }
  }

  // Reduce reaction forces from all partitions
  real64 globalFaceReactions[6];
  for( int face = 0 ; face < 6 ; face++ )
  {
    MPI_Allreduce( &localFaceReactions[face],
                   &globalFaceReactions[face],
                   1,
                   MPI_DOUBLE,
                   MPI_SUM,
                   MPI_COMM_GEOSX );
  }

  // Get end-of-step domain dimensions - note that m_domainExtent is updated later
  real64 length, width, height;
  length = m_domainExtent[0] * (1.0 + m_domainL[0] * dt);
  width  = m_domainExtent[1] * (1.0 + m_domainL[1] * dt);
  height = m_domainExtent[2] * (1.0 + m_domainL[2] * dt);

  // Write global reactions to file
  if( MpiWrapper::commRank( MPI_COMM_GEOSX ) == 0 )
  {
    std::ofstream file;
    // can't enable exception now because of gcc bug that raises ios_base::failure with useless message
    // file.exceptions(file.exceptions() | std::ios::failbit);
    file.open( "reactionHistory.csv", std::ios::out | std::ios::app );
    if( file.fail() )
    {
      throw std::ios_base::failure( std::strerror( errno ) );
    }
    // make sure write fails with exception if something is wrong
    file.exceptions( file.exceptions() | std::ios::failbit | std::ifstream::badbit );
    file << std::setprecision( std::numeric_limits<long double>::digits10 )
         << time_n + dt << ","
         << m_domainF[0] << ","
         << m_domainF[1] << ","
         << m_domainF[2] << ","
         << length << ","
         << width << ","
         << height << ","
         << globalFaceReactions[0] << ","
         << globalFaceReactions[1] << ","
         << globalFaceReactions[2] << ","
         << globalFaceReactions[3] << ","
         << globalFaceReactions[4] << ","
         << globalFaceReactions[5] << ","
         << m_domainL[0] << ","
         << m_domainL[1] << ","
         << m_domainL[2] << std::endl;
    file.close();
  }
}

void SolidMechanicsMPM::computeGridSurfaceNormals( ParticleManager & particleManager,
                                                   arrayView2d< real64, nodes::REFERENCE_POSITION_USD > const gridPosition,
                                                   array3d< real64 > & gridSurfaceNormal )
{
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    // Registered by subregion
    arrayView1d< real64 > const particleVolume = subRegion.getParticleVolume();
    arrayView1d< int > const particleGroup = subRegion.getParticleGroup();

    // Initialize mapping helpers
    int numNodesMappedTo = subRegion.numNodesMappedTo();
    array1d< int > nodeIDs(numNodesMappedTo); // nodes that the particle maps to
    array1d< real64 > weights(numNodesMappedTo); // shape function value for each node
    array2d< real64 > gradWeights(numNodesMappedTo, 3); // shape function gradient value for each node; 1st index = direction, 2nd index = node

    for( localIndex& p: subRegion.nonGhostIndices() )
    {
      // Get interpolation kernel
      subRegion.getAllWeights( p,
                               m_xLocalMin,
                               m_hEl,
                               m_ijkMap,
                               gridPosition,
                               nodeIDs,       // output
                               weights,       // output
                               gradWeights ); // output

      // Map to grid
      for( int g=0; g<numNodesMappedTo; g++ )
      {
        int mappedNode = nodeIDs[g];
        int fieldIndex = particleGroup[p];
        for(int i=0; i<m_numDims; i++)
        {
          gridSurfaceNormal[mappedNode][fieldIndex][i] += gradWeights[g][i] * particleVolume[p];
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
        arraySlice1d< real64 > const surfaceNormal = gridSurfaceNormal[g][fieldIndex];
        real64 norm = m_planeStrain == 1 ? sqrt( surfaceNormal[0] * surfaceNormal[0] + surfaceNormal[1] * surfaceNormal[1] ) : tOps::l2Norm< 3 >( surfaceNormal );
        if( norm > 0.0 ) // TODO: Set a finite threshold?
        {
          tOps::scale< 3 >( surfaceNormal, 1.0/norm );
        }
        else
        {
          surfaceNormal[0] = 1.0;
          surfaceNormal[1] = 0.0;
          surfaceNormal[2] = 0.0;
        }
      }
      else
      {
        gridSurfaceNormal[g][fieldIndex][0] = 1.0;
        gridSurfaceNormal[g][fieldIndex][1] = 0.0;
        gridSurfaceNormal[g][fieldIndex][2] = 0.0;
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
  real64 nAB[3];
  if( mA > mB )
  {
    nAB[0] = nA[0];
    nAB[1] = nA[1];
    nAB[2] = nA[2];
  }
  else
  {
    nAB[0] = -nB[0];
    nAB[1] = -nB[1];
    nAB[2] = -nB[2];
  }
  real64 norm = sqrt(nAB[0] * nAB[0] + nAB[1] * nAB[1] + nAB[2] * nAB[2]);
  nAB[0] /= norm;
  nAB[1] /= norm;
  nAB[2] /= norm;

  real64 weightedSpacingProjection = m_planeStrain ?
                                     0.65 * ( fabs(nAB[0]) * m_hEl[0] + fabs(nAB[1]) * m_hEl[1] ) :
                                     0.65 * ( fabs(nAB[0]) * m_hEl[0] + fabs(nAB[1]) * m_hEl[1] + fabs(nAB[2]) * m_hEl[2] ); // TODO: What should these weights be?
  real64 gap = (xB[0] - xA[0]) * nAB[0] + (xB[1] - xA[1]) * nAB[1] + (xB[2] - xA[2]) * nAB[2] - weightedSpacingProjection;

  // Total momentum for the contact pair.
  real64 qAB[3];
  qAB[0] = qA[0] + qB[0];
  qAB[1] = qA[1] + qB[1];
  qAB[2] = qA[2] + qB[2];

  // Center-of-mass velocity for the contact pair.
  real64 vAB[3];
  vAB[0] = qAB[0] / mAB;
  vAB[1] = qAB[1] / mAB;
  vAB[2] = qAB[2] / mAB;

  // Compute s1AB and s2AB, to form an orthonormal basis. This uses the method by E. Herbold
  // to ensure consistency between surfaces.
  real64 s1AB[3], s2AB[3]; // Tangential vectors for the contact pair
  computeOrthonormalBasis( nAB, s1AB, s2AB );

  // Compute force decomposition, declare increment in fA from "this" contact pair
  real64 fnor =  ( mA / dt ) * ( (vAB[0] - vA[0]) * nAB[0]  + (vAB[1] - vA[1]) * nAB[1]  + (vAB[2] - vA[2]) * nAB[2] ),
         ftan1 = ( mA / dt ) * ( (vAB[0] - vA[0]) * s1AB[0] + (vAB[1] - vA[1]) * s1AB[1] + (vAB[2] - vA[2]) * s1AB[2] ),
         ftan2 = ( mA / dt ) * ( (vAB[0] - vA[0]) * s2AB[0] + (vAB[1] - vA[1]) * s2AB[1] + (vAB[2] - vA[2]) * s2AB[2] );
  real64 dfA[3];

  // Check for separability, and enforce either slip, or no-slip contact, accordingly
  if( separable == 0 )
  {
    // Surfaces are bonded, treat as single velocity field by applying normal force to prevent
    // interpenetration, and tangential force to prevent slip.
    dfA[0] = fnor * nAB[0] + ftan1 * s1AB[0] + ftan2 * s2AB[0];
    dfA[1] = fnor * nAB[1] + ftan1 * s1AB[1] + ftan2 * s2AB[1];
    dfA[2] = fnor * nAB[2] + ftan1 * s1AB[2] + ftan2 * s2AB[2];
    fA[0] += dfA[0];
    fA[1] += dfA[1];
    fA[2] += dfA[2];
    fB[0] -= dfA[0];
    fB[1] -= dfA[1];
    fB[2] -= dfA[2];
  }
  else
  {
    // Surfaces are separable. For frictional contact, apply a normal force to
    // prevent interpenetration, and tangential force to prevent slip unless f_tan > mu*f_nor
    real64 contact;
    real64 test = (vA[0] - vAB[0]) * nAB[0] + (vA[1] - vAB[1]) * nAB[1] + (vA[2] - vAB[2]) * nAB[2];
    if( m_contactGapCorrection == 1 ) // TODO: Implement soft/ramped contact option?
    {
      contact = test > 0.0 && gap < 0.0 ? 1.0 : 0.0;
    }
    else
    {
      contact = test > 0.0 ? 1.0 : 0.0;
    }

    // Modify normal contact force
    fnor *= contact;

    // Determine force for tangential sticking
    real64 ftanMag = sqrt( ftan1 * ftan1 + ftan2 * ftan2 );

    // Get direction of tangential contact force
    real64 sAB[3];
    if( ftanMag > 0.0 )
    {
      sAB[0] = (s1AB[0] * ftan1 + s2AB[0] * ftan2) / ftanMag;
      sAB[0] = (s1AB[1] * ftan1 + s2AB[1] * ftan2) / ftanMag;
      sAB[0] = (s1AB[2] * ftan1 + s2AB[2] * ftan2) / ftanMag;
    }
    else
    {
      sAB[0] = 0.0;
      sAB[1] = 0.0;
      sAB[2] = 0.0;
    }

    // Update fA and fB - tangential force is friction bounded by sticking force
    real64 ftan = std::min( m_frictionCoefficient * std::abs( fnor ), ftanMag ); // This goes to zero when contact=0 due to the std::min
    dfA[0] = fnor * nAB[0] + ftan * sAB[0];
    dfA[1] = fnor * nAB[1] + ftan * sAB[1];
    dfA[2] = fnor * nAB[2] + ftan * sAB[2];
    fA[0] += dfA[0];
    fA[1] += dfA[1];
    fA[2] += dfA[2];
    fB[0] -= dfA[0];
    fB[1] -= dfA[1];
    fB[2] -= dfA[2];
  }
}

void SolidMechanicsMPM::computeOrthonormalBasis( const real64* e1, // input "normal" unit vector.
                                                 real64* e2,       // output "tangential" unit vector.
                                                 real64* e3 )      // output "tangential" unit vector.
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
  if( maxp1 >= maxp2 && maxp1 >= maxp3 )
  {
    e2x = e1y;
    e2y = -e1x;
    e2z = 0.0;
  }
  else if( maxp2 >= maxp1 && maxp2 >= maxp3 )
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
  real64 e3norm = sqrt( e3[0] * e3[0] + e3[1] * e3[1] + e3[2] * e3[2] );
  e3[0] /= e3norm;
  e3[1] /= e3norm;
  e3[2] /= e3norm;
}

void SolidMechanicsMPM::setGridFieldLabels( NodeManager & nodeManager )
{
  // Generate labels
  std::vector< std::string > labels1(m_numVelocityFields);
  std::generate( labels1.begin(), labels1.end(), [i=0]() mutable { return "velocityField" + std::to_string( i++ ); } );
  string const labels2[] = { "X", "Y", "Z" };

  // Apply labels to scalar multi-fields
  std::vector< std::string > keys2d = { viewKeyStruct::massString(),
                                        viewKeyStruct::damageString(),
                                        viewKeyStruct::maxDamageString()};
  for( size_t gridField=0; gridField<keys2d.size(); gridField++)
  {
    WrapperBase & wrapper = nodeManager.getWrapper< array2d< real64 > >( keys2d[gridField] );
    wrapper.setDimLabels( 1, labels1 );
  }

  // Apply labels to vector multi-fields
  std::vector< std::string > keys3d = { viewKeyStruct::velocityString(),
                                        viewKeyStruct::momentumString(),
                                        viewKeyStruct::accelerationString(),
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

void SolidMechanicsMPM::resizeGrid( SpatialPartition & partition,
                                    arrayView2d< real64, nodes::REFERENCE_POSITION_USD > const gridPosition,
                                    real64 const dt )
{
  // Modify SpatialPartition class members
  partition.updateSizes( m_domainL, dt );

  // Modify SolidMechanicsMPM class members
  for(int i=0; i<3; i++)
  {
    // Incremental stretch
    real64 ratio = 1.0 + m_domainL[i] * dt;

    // Modify SolidMechanicsMPM class members
    m_hEl[i] *= ratio;
    m_xLocalMin[i] *= ratio;
    m_xLocalMax[i] *= ratio;
    m_xLocalMinNoGhost[i] *= ratio;
    m_xLocalMaxNoGhost[i] *= ratio;
    m_xGlobalMin[i] *= ratio;
    m_xGlobalMax[i] *= ratio;
    m_partitionExtent[i] *= ratio;
    m_domainExtent[i] *= ratio;

    // Update nodal positions
    for( int g=0; g<gridPosition.size(0); g++ )
    {
      gridPosition[g][i] *= ratio;
    }
  }
}

void SolidMechanicsMPM::solverProfiling( std::string label )
{
  if( m_solverProfiling >= 1 )
  {
    MPI_Barrier( MPI_COMM_GEOSX );
    GEOSX_LOG_RANK_IF( m_solverProfiling == 2, label );    
    m_profilingTimes.push_back( MPI_Wtime() );
    m_profilingLabels.push_back( label );
  }
}

void SolidMechanicsMPM::solverProfilingIf( std::string label, bool condition )
{
  if( condition )
  {
    solverProfiling( label );
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

real64 SolidMechanicsMPM::computeNeighborList( ParticleManager & particleManager )
{
  // Time this function
  real64 tStart = MPI_Wtime();
  
  // Expand bin limits by neighbor radius to account for the buffer zone of ghost particles outside the patch limits
  real64 neighborRadiusSquared = m_neighborRadius * m_neighborRadius;
  real64 xmin = m_xLocalMinNoGhost[0] - m_neighborRadius,
         xmax = m_xLocalMaxNoGhost[0] + m_neighborRadius,
         ymin = m_xLocalMinNoGhost[1] - m_neighborRadius,
         ymax = m_xLocalMaxNoGhost[1] + m_neighborRadius,
         zmin = m_xLocalMinNoGhost[2] - m_neighborRadius,
         zmax = m_xLocalMaxNoGhost[2] + m_neighborRadius;

  // Initialize bin sort
  real64 binWidth = m_binSizeMultiplier * m_neighborRadius;
  int nxbins = std::ceil( ( xmax - xmin ) / binWidth ),
      nybins = std::ceil( ( ymax - ymin ) / binWidth ),
      nzbins = m_planeStrain ? 1 : std::ceil( ( zmax - zmin ) / binWidth );
  int nbins = nxbins * nybins * nzbins;
  real64 dx = ( xmax - xmin ) / nxbins,
         dy = ( ymax - ymin ) / nybins,
         dz = ( zmax - zmin ) / nzbins;

  // Declare and add particles to bins
  std::vector<std::vector< ArrayOfArrays<int> >> bins; // This is horrible
  bins.resize( particleManager.numRegions() );
  particleManager.forParticleRegions< ParticleRegion >( [&]( ParticleRegion & region ) // idk why this requires a template argument and the subregion loops don't
  {
    localIndex regionIndex = region.getIndexInParent();
    bins[regionIndex].resize( region.numSubRegions() );
    region.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
    {
      localIndex subRegionIndex = subRegion.getIndexInParent();
      bins[regionIndex][subRegionIndex].resize( nbins );
      arrayView2d< real64 > const p_x = subRegion.getParticleCenter();
      for( localIndex p=0; p<subRegion.size(); p++ )
      {
        // Particle bin ijk indices
        int i, j, k;
        i = std::floor( ( p_x[p][0] - xmin ) / dx ),
        j = std::floor( ( p_x[p][1] - ymin ) / dy ),
        k = std::floor( ( p_x[p][2] - zmin ) / dz );

        // Bin number
        int b = i + j * nxbins + k * nxbins * nybins;
        bins[regionIndex][subRegionIndex].emplaceBack( b, p );
      }
    } );
  } );

  // Perform neighbor search over appropriate bins
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegionA )
  {
    // Declare relative position    
    real64 xBA[3];

    // Get and initialize neighbor list
    OrderedVariableToManyParticleRelation & neighborList = subRegionA.neighborList();
    neighborList.resize(0); // Clear the existing neighbor list. It would be better to resize each particle's array to zero to avoid the next line
    neighborList.resize(subRegionA.size());

    // Get 'this' particle's location
    arrayView2d< real64 > const xA = subRegionA.getParticleCenter();

    // Find neighbors of 'this' particle
    for( localIndex& a: subRegionA.nonGhostIndices() )
    {
      // Bin ijk indices bounding a sphere of radius m_neighborRadius centered at 'this' particle
      int imin, imax, jmin, jmax, kmin, kmax;
      imin = std::floor( ( xA[a][0] - m_neighborRadius - xmin ) / dx ),
      jmin = std::floor( ( xA[a][1] - m_neighborRadius - ymin ) / dy ),
      kmin = std::floor( ( xA[a][2] - m_neighborRadius - zmin ) / dz );
      imax = std::floor( ( xA[a][0] + m_neighborRadius - xmin ) / dx ),
      jmax = std::floor( ( xA[a][1] + m_neighborRadius - ymin ) / dy ),
      kmax = std::floor( ( xA[a][2] + m_neighborRadius - zmin ) / dz );

      // Adjust bin ijk indices if necessary
      imin = std::max( imin, 0 );
      imax = std::min( imax, nxbins-1 );
      jmin = std::max( jmin, 0 );
      jmax = std::min( jmax, nybins-1 );
      kmin = std::max( kmin, 0 );
      kmax = std::min( kmax, nzbins-1 );

      // Inner subregion loop
      particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegionB )
      {
        // Get region and subregion indices
        ParticleRegion & region = dynamicCast< ParticleRegion & >( subRegionB.getParent().getParent() );
        localIndex regionIndex = region.getIndexInParent();
        localIndex subRegionIndex = subRegionB.getIndexInParent();

        // Get 'other' particle location
        arrayView2d< real64 > const xB = subRegionB.getParticleCenter();

        // Loop over bins
        for( int iBin=imin; iBin<=imax; iBin++ )
        {
          for( int jBin=jmin; jBin<=jmax; jBin++ )
          {
            for( int kBin=kmin; kBin<=kmax; kBin++ )
            {
              int bin = iBin + jBin * nxbins + kBin * nxbins * nybins;
              for( localIndex & b: bins[regionIndex][subRegionIndex][bin] )
              {
                xBA[0] = xB[b][0] - xA[a][0];
                xBA[1] = xB[b][1] - xA[a][1];
                xBA[2] = xB[b][2] - xA[a][2];
                real64 rSquared = xBA[0] * xBA[0] + xBA[1] * xBA[1] + xBA[2] * xBA[2];
                if( rSquared <= neighborRadiusSquared ) // Would you be my neighbor?
                {
                  insert( neighborList,
                          a,
                          regionIndex,
                          subRegionIndex,
                          b );
                }
              }
            }
          }
        }
      } );
    }
  } );

  return( MPI_Wtime() - tStart );
}

void SolidMechanicsMPM::optimizeBinSort( ParticleManager & particleManager )
{
  // Each partition determines its optimal multiplier which results in the minimum time for neighbor list construction
  // The global multiplier is set by a weighted average of each partition's multiplier, with the number of particles
  // on the partition being the weight factor.

  // Start
  int optimalMultiplier;

  // Identify the largest possible multiplier - we can limit this if it's prohibitive to check larger multipliers in large 3D sims
  real64 xL = m_xLocalMaxNoGhost[0] - m_xLocalMinNoGhost[0] + 2 * m_neighborRadius;
  real64 yL = m_xLocalMaxNoGhost[1] - m_xLocalMinNoGhost[1] + 2 * m_neighborRadius;
  real64 zL = m_xLocalMaxNoGhost[2] - m_xLocalMinNoGhost[2] + 2 * m_neighborRadius;
  int maxMultiplier = std::max( std::ceil( xL / m_neighborRadius), std::max( std::ceil( yL / m_neighborRadius), std::ceil( zL / m_neighborRadius) ));
  maxMultiplier = std::max( maxMultiplier, 1 );

  // Identify this partition's optimal multiplier
  real64 minTime = DBL_MAX;
  for( int multiplier=1; multiplier<=maxMultiplier; multiplier++)
  {
    m_binSizeMultiplier = multiplier;
    real64 sortingTime = computeNeighborList( particleManager );
    if( sortingTime < minTime )
    {
      minTime = sortingTime;
      optimalMultiplier = multiplier;
    }
  }

  // MPI comms
  int globalNumberOfParticles;
  real64 globalWeightedMultiplier;
  int localNumberOfParticles = particleManager.getNumberOfParticles();
  real64 localWeightedMultiplier = optimalMultiplier * localNumberOfParticles;
  MPI_Allreduce( &localWeightedMultiplier,
                 &globalWeightedMultiplier,
                 1,
                 MPI_DOUBLE,
                 MPI_SUM,
                 MPI_COMM_GEOSX );
  MPI_Allreduce( &localNumberOfParticles,
                 &globalNumberOfParticles,
                 1,
                 MPI_INT,
                 MPI_SUM,
                 MPI_COMM_GEOSX );

  // Set bin size multiplier
  m_binSizeMultiplier = std::max( (int) std::round( globalWeightedMultiplier / globalNumberOfParticles ), 1 );
}

real64 SolidMechanicsMPM::kernel( const real64 & r ) // distance from particle to query point.
{
  // Compute the value of a particle kernel function at some point.
  
  const real64 R = m_neighborRadius; // kernel Radius
  const real64 norm = m_planeStrain ? 1.0610329539459689051/(R*R) : 1.1936620731892150183/(R*R*R); // kernel should integrate to unity

  if( r < R )
  {
    return norm * ( 1 - 3.0 * r * r / ( R * R ) + 2.0 * r * r * r / ( R * R * R ) );
  }
  else
  {
    return ( 0.0 );
  }
}

void SolidMechanicsMPM::kernelGradient( arraySlice1d< real64 > const x,  // query point
                                        arraySlice1d< real64 > const xp, // particle location
                                        const real64 & r,                // distance from particle to query point.
                                        real64* result )
{
  // Compute the value of a particle kernel function gradient at some point.
  
  real64 s;
  const real64 R = m_neighborRadius; // kernel Radius
  const real64 norm = m_planeStrain ? 1.0610329539459689051/(R*R) : 1.1936620731892150183/(R*R*R); // kernel should integrate to unity

  if( r < R )
  {
    s = norm * 6.0 * ( r - R ) / ( R * R * R );
    result[0] = s * (x[0] - xp[0]);
    result[1] = s * (x[1] - xp[1]);
    result[2] = s * (x[2] - xp[2]);
  }
  else
  {
    result[0] = 0.0;
    result[1] = 0.0;
    result[2] = 0.0;
  }
}

real64 SolidMechanicsMPM::computeKernelField( arraySlice1d< real64 > const x,  // query point
                                              arrayView2d< real64 > const xp,  // List of neighbor particle locations.
                                              arrayView1d< real64 > const Vp,  // List of neighbor particle volumes.
                                              arrayView1d< real64 > const fp ) // scalar field values (e.g. damage) at neighbor particles
{
  // Compute the kernel scalar field at a point, for a given list of neighbor particles.
  // The lists xp, fp, and the length np could refer to all the particles in the patch,
  // but generally this function will be evaluated with x equal to some particle center,
  // and xp, fp, will be lists for just the neighbors of the particle.

  // Initialize
  real64 relativePosition[3];
  real64 kernelVal,
         f = 0.0,
         k = 0.0,
         r;

  // Sum
  for( localIndex p = 0 ; p < Vp.size() ; ++p )
  {
    relativePosition[0] = x[0] - xp[p][0];
    r = sqrt( relativePosition[0] * relativePosition[0] + relativePosition[1] * relativePosition[1] + relativePosition[2] * relativePosition[2] );
    kernelVal = kernel( r );
    k += Vp[p] * kernelVal;
    f += Vp[p] * fp[p] * kernelVal;
  }

  // Return the normalized kernel field (which eliminates edge effects)
  if( k > 0.0 )
  {
    return ( f / k );
  }
  else
  {
    return ( 0.0 );
  }
}

void SolidMechanicsMPM::computeKernelFieldGradient( arraySlice1d< real64 > const x,  // query point
                                                    arrayView2d< real64 > const xp,  // List of neighbor particle locations.
                                                    arrayView1d< real64 > const Vp,  // List of neighbor particle volumes.
                                                    arrayView1d< real64 > const fp,  // scalar field values (e.g. damage) at neighbor particles
                                                    arraySlice1d< real64 > const result )
{
  // Compute the kernel scalar field at a point, for a given list of neighbor particles.
  // The lists xp, fp, and the length np could refer to all the particles in the patch,
  // but generally this function will be evaluated with x equal to some particle center,
  // and xp, fp, will be lists for just the neighbors of the particle.
  // TODO: Modify to also "return" the kernel field value

  // Scalar kernel field values
  real64 kernelVal,
         f = 0.0,
         k = 0.0,
         r;

  // Gradient of the scalar field
  real64 relativePosition[3],
         fGrad[3] = {0.0, 0.0, 0.0},
         kGrad[3] = {0.0, 0.0, 0.0},
         kernelGradVal[3];

  for( localIndex p = 0 ; p < Vp.size() ; ++p )
  {
    relativePosition[0] = x[0] - xp[p][0];
    relativePosition[1] = x[1] - xp[p][1];
    relativePosition[2] = x[2] - xp[p][2];
    r = sqrt( relativePosition[0] * relativePosition[0] + relativePosition[1] * relativePosition[1] + relativePosition[2] * relativePosition[2] );

    kernelVal = kernel( r );
    k += Vp[p] * kernelVal;
    f += Vp[p] * fp[p] * kernelVal;

    kernelGradient( x, xp[p], r, kernelGradVal );
    kGrad[0] += kernelGradVal[0] * Vp[p];
    kGrad[1] += kernelGradVal[1] * Vp[p];
    kGrad[2] += kernelGradVal[2] * Vp[p];
    fGrad[0] += kernelGradVal[0] * Vp[p] * fp[p];
    fGrad[1] += kernelGradVal[1] * Vp[p] * fp[p];
    fGrad[2] += kernelGradVal[2] * Vp[p] * fp[p];
  }

  // Return the normalized kernel field gradient (which eliminates edge effects)
  if( k > 0.0 )
  {
    //kernelField = f/k;
    result[0] = fGrad[0] / k - f * kGrad[0] / (k * k);
    result[1] = fGrad[1] / k - f * kGrad[1] / (k * k);
    result[2] = fGrad[2] / k - f * kGrad[2] / (k * k);
  }
  else
  {
    //kernelField = 0.0;
    result[0] = 0.0;
    result[1] = 0.0;
    result[2] = 0.0;
  }
}

REGISTER_CATALOG_ENTRY( SolverBase, SolidMechanicsMPM, string const &, dataRepository::Group * const )
}

