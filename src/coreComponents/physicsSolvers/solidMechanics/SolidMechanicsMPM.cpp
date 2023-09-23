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
#include "constitutive/ConstitutivePassThruHandler.hpp"

#include "events/mpmEvents/MPMEventBase.hpp"
#include "events/mpmEvents/MaterialSwapMPMEvent.hpp"
#include "events/mpmEvents/AnnealMPMEvent.hpp"
#include "events/mpmEvents/HealMPMEvent.hpp"
#include "events/mpmEvents/CrystalHealMPMEvent.hpp"
#include "events/mpmEvents/InsertPeriodicContactSurfacesMPMEvent.hpp"
#include "events/mpmEvents/MachineSampleMPMEvent.hpp"
#include "events/mpmEvents/FrictionCoefficientSwapMPMEvent.hpp"

namespace geos
{

using namespace dataRepository;
using namespace constitutive;

SolidMechanicsMPM::SolidMechanicsMPM( const string & name,
                                      Group * const parent ):
  SolverBase( name, parent ),
  m_solverProfiling( 0 ),
  m_timeIntegrationOption( TimeIntegrationOption::ExplicitDynamic ),
  m_updateMethod( UpdateMethodOption::FLIP ),
  m_iComm( CommunicationTools::getInstance().getCommID() ),
  m_prescribedBcTable( 0 ),
  m_boundaryConditionTypes(),
  m_bcTable(),
  m_prescribedFTable( 0 ),
  m_prescribedBoundaryFTable( 0 ),
  m_fTableInterpType( 0 ),
  m_fTable(),
  m_domainF(),
  m_domainL(),
  m_bodyForce(),
  m_stressControl(),
  m_stressTableInterpType( 0 ),
  m_stressControlKp( 0.1 ),
  m_stressControlKi( 0 ),
  m_stressControlKd( 0 ),
  m_boxAverageHistory( 0 ),
  m_reactionHistory( 0 ),
  m_needsNeighborList( 0 ),
  m_neighborRadius( -1.0 ),
  m_binSizeMultiplier( 1 ),
  m_cpdiDomainScaling( 0 ),
  m_smallMass( DBL_MAX ),
  m_numContactGroups(),
  m_numContactFlags(),
  m_numVelocityFields(),
  m_separabilityMinDamage( 0.5 ),
  m_treatFullyDamagedAsSingleField( 1 ),
  m_surfaceDetection( 0 ),
  m_damageFieldPartitioning( 0 ),
  m_contactGapCorrection( 0 ),
  // m_directionalOverlapCorrection( 0 ),
  m_frictionCoefficient(),
  m_frictionCoefficientTable(),
  m_planeStrain( 0 ),
  m_numDims( 3 ),
  m_hEl{ DBL_MAX, DBL_MAX, DBL_MAX },
  m_xLocalMin{ DBL_MAX, DBL_MAX, DBL_MAX },
  m_xLocalMax{ -DBL_MAX, -DBL_MAX, -DBL_MAX },
  m_xLocalMinNoGhost{ 0.0, 0.0, 0.0 },
  m_xLocalMaxNoGhost{ 0.0, 0.0, 0.0 },
  m_xGlobalMin{ 0.0, 0.0, 0.0 },
  m_xGlobalMax{ 0.0, 0.0, 0.0 },
  m_partitionExtent{ 0.0, 0.0, 0.0 },
  m_nEl{ 0, 0, 0 },
  m_ijkMap(),
  m_mpmEventManager( nullptr ),
  m_surfaceHealing( false )
{
  // setInputFlags( InputFlags::OPTIONAL );

  registerWrapper( "solverProfiling", &m_solverProfiling ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flag for timing subroutines in the solver" );

  registerWrapper( viewKeyStruct::timeIntegrationOptionString(), &m_timeIntegrationOption ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_timeIntegrationOption ).
    setDescription( "Time integration method. Options are:\n* " + EnumStrings< TimeIntegrationOption >::concat( "\n* " ) );

  registerWrapper( "updateMethod", &m_updateMethod ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_updateMethod ).
    setDescription( "Update method. Options are:\n* " + EnumStrings< UpdateMethodOption >::concat( "\n* ") );

  registerWrapper( "updateOrder", &m_updateOrder ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 2 ).
    setDescription("Order for update method, only applies to XPIC and FMPM");

  registerWrapper( "prescribedBcTable", &m_prescribedBcTable ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0 ).
    setDescription( "Flag for whether to have time-dependent boundary condition types" );

  registerWrapper( "boxAverageHistory", &m_boxAverageHistory ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0 ).
    setDescription( "Flag for whether to output box average history" );

  registerWrapper( "boxAverageWriteInterval", &m_boxAverageWriteInterval ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0 ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Interval between writing box averages to files" );

  registerWrapper( "boxAverageMin", &m_boxAverageMin).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Minimum corner position of box average" );

  registerWrapper( "boxAverageMax", &m_boxAverageMax).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Maximum corner position of box average" );

  registerWrapper( "nextBoxAverageWriteTime", &m_nextBoxAverageWriteTime ).
    setInputFlag( InputFlags::FALSE ).
    setApplyDefaultValue( 0 ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Next time to write box averages" );

  registerWrapper( "reactionHistory", &m_reactionHistory ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0 ).
    setDescription( "Flag for whether to output face reaction history" );

  registerWrapper( "reactionWriteInterval", &m_reactionWriteInterval ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0 ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Interval between writing reactions to files" );

  registerWrapper( "nextReactionWriteTime", &m_nextReactionWriteTime ).
    setInputFlag( InputFlags::FALSE ).
    setApplyDefaultValue( 0 ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Next time to write reactions" );

  registerWrapper( "boundaryConditionTypes", &m_boundaryConditionTypes ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Boundary conditions on x-, x+, y-, y+, z- and z+ faces. Options are:\n* " + EnumStrings< BoundaryConditionOption >::concat( "\n* " ) );

  registerWrapper( "bodyForce", &m_bodyForce ).
    setInputFlag(InputFlags::OPTIONAL).
    setDescription( "Array that stores uniform body force" );

  registerWrapper( "bcTable", &m_bcTable ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Array that stores time-dependent bc types on x-, x+, y-, y+, z- and z+ faces." );

  registerWrapper( "prescribedFTable", &m_prescribedFTable ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flag for whether to have time-dependent superimposed velocity gradient for triply periodic simulations" );

  registerWrapper( "prescribedBoundaryFTable", &m_prescribedBoundaryFTable ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flag for whether to have time-dependent boundary conditions described by a global background grid F" );

  registerWrapper( "fTableInterpType", &m_fTableInterpType ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "The type of F table interpolation. Options are 0 (linear), 1 (cosine), 2 (quintic polynomial)." );

  registerWrapper( "fTable", &m_fTable ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Array that stores time-dependent grid-aligned stretches interpreted as a gloabl background grid F read from the XML file." );

  registerWrapper( "stressControl" , &m_stressControl).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flag for whether stress control using box averages is enabled" );

  registerWrapper( "stressTableInterpType", &m_stressTableInterpType ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "The type of stress table interpolation. Options are 0 (linear), 1 (cosine), 2 (quintic polynomial)." );

  registerWrapper( "stressTable", &m_stressTable ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Array that stores the time-depended grid aligned stresses" );

  registerWrapper( "stressControlKp", &m_stressControlKp ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Proportional gain of stress PID controller" );

  registerWrapper( "stressControlKi", &m_stressControlKi ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Integral gain of stress PID controller" );

   registerWrapper( "stressControlKd", &m_stressControlKd ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Derivative gain of stress PID controller" );

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

  registerWrapper( "generalizedVortexMMS", &m_generalizedVortexMMS ).
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

  registerWrapper( "separabilityMinDamage", &m_separabilityMinDamage ).
    setApplyDefaultValue( 0.5 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Damage threshold for field separability" );

  registerWrapper( "treatFullyDamagedAsSingleField", &m_treatFullyDamagedAsSingleField ).
    setApplyDefaultValue( 1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Whether to consolidate fully damaged fields into a single field. Nice for modeling damaged mush." );

  registerWrapper( "surfaceDetection", &m_surfaceDetection ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flag for automatic surface detection on the 1st cycle" );

  registerWrapper( "damageFieldPartitioning", &m_damageFieldPartitioning ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flag for using the gradient of the particle damage field to partition material into separate velocity fields" );

  registerWrapper( "contactGapCorrection", &m_contactGapCorrection ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flag for mitigating contact gaps" );

  // registerWrapper( "directionalOverlapCorrection", &m_directionalOverlapCorrection ).
  //   setApplyDefaultValue( 0 ).
  //   setInputFlag( InputFlags::OPTIONAL ).
  //   setDescription( "Flag for mitigating pile-up of particles at contact interfaces" );

  registerWrapper( "frictionCoefficient", &m_frictionCoefficient ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Coefficient of friction, currently assumed to be the same everywhere" );

  registerWrapper( "frictionCoefficientTable", &m_frictionCoefficientTable ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Friction coefficient table for different groups" );

  registerWrapper( "planeStrain", &m_planeStrain ).
    setApplyDefaultValue( false ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flag for performing plane strain calculations" );

  registerWrapper( "numDims", &m_numDims ).
    setApplyDefaultValue( 3 ).
    setInputFlag( InputFlags::FALSE ).
    setDescription( "The number of active spatial dimensions, 2 for plane strain, 3 otherwise" );

  registerWrapper( "m_ijkMap", &m_ijkMap ).
    setInputFlag( InputFlags::FALSE ).
    setDescription( "Map from indices in each spatial dimension to local node ID" );

  registerWrapper("useEvents", &m_useEvents).
    setInputFlag( InputFlags::OPTIONAL ).
    setDefaultValue( 0 ).
    setDescription( "Enable events" );

  registerWrapper( "debugFlag", &m_debugFlag ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDefaultValue( 0 ).
    setDescription( "Enables debugging of MPM explicit timestep" );

  m_mpmEventManager = &registerGroup< MPMEventManager >( groupKeys.mpmEventManager );
}


void SolidMechanicsMPM::postProcessInput()
{
  SolverBase::postProcessInput();

  // Activate neighbor list if necessary
  if( m_damageFieldPartitioning == 1 || m_surfaceDetection == 1 /*|| m_directionalOverlapCorrection == 1*/ )
  {
    m_needsNeighborList = 1;
  }

  // Set number of active dimensions based on m_planeStrain
  m_numDims = m_planeStrain ? 2 : 3;

  // Throw error if boundary conditions are incorrectly specified
  GEOS_ERROR_IF( m_boundaryConditionTypes.size() != 6 && m_boundaryConditionTypes.size() > 0,
                 "boundaryConditionTypes must be of length 6. "
                 "The 6 entries correspond to BCs on the x-, x+, y-, y+, z- and z+ faces." );

  // Initialize boundary condition types if they're not specified by the user
  if( m_boundaryConditionTypes.size() == 0 )
  {
    m_boundaryConditionTypes.resize( 6 );
    LvArray::tensorOps::fill< 6 >( m_boundaryConditionTypes, 0 );
  }

  // Throw error if boundary conditions are incorrectly specified
  GEOS_ERROR_IF( m_bodyForce.size() != 3 && m_bodyForce.size() > 0,
                 "bodyForce must be of length 3. ");

  //Initialize body force if they're not specified by the user
  if(m_bodyForce.size() == 0)
  {
    m_bodyForce.resize(3);
    LvArray::tensorOps::fill< 3 >(m_bodyForce, 0.0 );
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
      GEOS_LOG_RANK("SolidMechanicsMPM::registerDataOnMesh, " << regionNames);
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
        subRegion.registerField< particleCrystalHealFlag >( getName() );
        subRegion.registerField< particleMass >( getName() );
        subRegion.registerField< particleInitialVolume >( getName() );
        subRegion.registerField< particleDensity >( getName() );
        subRegion.registerField< particleOverlap >( getName() );

        // Double-indexed fields (vectors and symmetric tensors stored in Voigt notation)
        subRegion.registerField< particleBodyForce >( getName() ).reference().resizeDimension< 1 >( 3 );
        subRegion.registerField< particleStress >( getName() ).setDimLabels( 1, voightLabels ).reference().resizeDimension< 1 >( 6 );
        subRegion.registerField< particlePlasticStrain >( getName() ).setDimLabels( 1, voightLabels ).reference().resizeDimension< 1 >( 6 );
        subRegion.registerField< particleDamageGradient >( getName() ).reference().resizeDimension< 1 >( 3 );
        subRegion.registerField< particleReferencePosition >( getName() ).reference().resizeDimension< 1 >( 3 );

        // Triple-indexed fields (vectors of vectors, non-symmetric tensors)
        subRegion.registerField< particleInitialRVectors >( getName() ).reference().resizeDimension< 1, 2 >( 3, 3 );
        subRegion.registerField< particleDeformationGradient >( getName() ).reference().resizeDimension< 1, 2 >( 3, 3 );
        subRegion.registerField< particleFDot >( getName() ).reference().resizeDimension< 1, 2 >( 3, 3 );
        subRegion.registerField< particleVelocityGradient >( getName() ).reference().resizeDimension< 1, 2 >( 3, 3 );
        subRegion.registerField< particleSphF >( getName() ).reference().resizeDimension< 1, 2 >( 3, 3 );
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

      nodes.registerWrapper< array2d< real64 > >( viewKeyStruct::damageGradientString() ).
        setPlotLevel( PlotLevel::LEVEL_1 ).
        setRegisteringObjects( this->getName() ).
        setDescription( "An array that holds the result of mapping particle damage gradients to the nodes." );

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
                                                    arrayView1d< string const > const & regionNames )
  {
    MeshBody const & meshBody = meshBodies.getGroup< MeshBody >( meshBodyName );

    if( meshBody.hasParticles() ) // Only particle regions will hold actual materials. Background grid currently holds a null material so
                                  // that the input file parser doesn't complain, but we don't need to actually do anything with it.
    {
      ParticleManager & particleManager = meshLevel.getParticleManager();
      particleManager.forParticleSubRegions< ParticleSubRegion >( regionNames, [&]( localIndex const,
                                                                                    ParticleSubRegion & subRegion )
      {
        string & solidMaterialName = subRegion.getReference< string >( viewKeyStruct::solidMaterialNamesString() );
        solidMaterialName = SolverBase::getConstitutiveName< SolidBase >( subRegion );
      } );
    }
  } );

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();

  FiniteElementDiscretizationManager const &
  feDiscretizationManager = numericalMethodManager.getFiniteElementDiscretizationManager();

  FiniteElementDiscretization const &
  feDiscretization = feDiscretizationManager.getGroup< FiniteElementDiscretization >( m_discretizationName );
  GEOS_UNUSED_VAR( feDiscretization );
}


bool SolidMechanicsMPM::execute( real64 const time_n,
                                 real64 const dt,
                                 integer const cycleNumber,
                                 integer const GEOS_UNUSED_PARAM( eventCounter ),
                                 real64 const GEOS_UNUSED_PARAM( eventProgress ),
                                 DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

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
  GEOS_MARK_FUNCTION;
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
    GEOS_ERROR( "MPM solver only currently supports explicit time stepping!" );
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
    // Reads the FTable directly from the xml
    int numRows = m_bcTable.size( 0 );

    GEOS_ERROR_IF(numRows == 0, "Prescribed boundary conditions is enabled but no bcTable was specified.");
    for(int i = 0; i < numRows; ++i){
      GEOS_ERROR_IF(m_bcTable[i].size() != 7, "BCtable row " << i+1 << " must have 7 elements.");
      
      GEOS_ERROR_IF(m_bcTable[i][0] < 0, "BCTable times must be positive.");

      GEOS_ERROR_IF(roundf(m_bcTable[i][1]) != m_bcTable[i][1] || 
                    roundf(m_bcTable[i][2]) != m_bcTable[i][1] || 
                    roundf(m_bcTable[i][3]) != m_bcTable[i][1] || 
                    roundf(m_bcTable[i][4]) != m_bcTable[i][1] || 
                    roundf(m_bcTable[i][5]) != m_bcTable[i][1] || 
                    roundf(m_bcTable[i][6]) != m_bcTable[i][1], "Only integer boundary condition types are permitted.");

      GEOS_ERROR_IF( ( m_bcTable[i][1] < 0 || 
                       m_bcTable[i][2] < 0 || 
                       m_bcTable[i][3] < 0 || 
                       m_bcTable[i][4] < 0 || 
                       m_bcTable[i][5] < 0 || 
                       m_bcTable[i][6] < 0 ) && 
                     ( m_bcTable[i][1] > 3 || 
                       m_bcTable[i][2] > 3 || 
                       m_bcTable[i][3] > 3 || 
                       m_bcTable[i][4] > 3 || 
                       m_bcTable[i][5] > 3 || 
                       m_bcTable[i][6] > 3 ), "Boundary types must be 0, 1, 2 or 3");
    }
  }

  // Initialize domain F and L, then read and distribute F table
  m_domainF.resize( 3 );
  m_domainL.resize( 3 );
  for( int i=0; i < 3; i++ )
  {
    m_domainF[i] = 1.0;
    m_domainL[i] = 0.0;
  }
  if( m_prescribedBoundaryFTable == 1 && m_prescribedFTable == 1 )
  {
    // Reads the FTable directly from the xml
    int numRows = m_fTable.size( 0 );
    GEOS_ERROR_IF(numRows == 0, "Prescribed boundary deformation is enabled but no fTable was specified.");
    for(int i = 0; i < numRows; ++i){
      GEOS_ERROR_IF(m_fTable[i].size() != 4, "Ftable row " << i+1 << " must have 4 elements.");
      GEOS_ERROR_IF(m_fTable[i][0] < 0, "FTable times must be positive.");

      if(i == 0)
      {
        GEOS_ERROR_IF( m_fTable[i][1] != 1.0 || 
                       m_fTable[i][2] != 1.0 || 
                       m_fTable[i][3] != 1.0 , "Deformation of first row of FTable must be 1." );
      }
      else
      {
        GEOS_ERROR_IF( m_fTable[i][1] < 0 || 
                       m_fTable[i][2] < 0 || 
                       m_fTable[i][3] < 0, "Deformations of FTable must be positive." );
      }
    }
  }

  // Check stress control
  if( m_stressControl.size() == 0 ){
    m_stressControl.resize( 3 );
    LvArray::tensorOps::fill<3>( m_stressControl, 0 );
  }

  if( m_stressControl[0] == 1 || m_stressControl[1] == 1 || m_stressControl[2] == 1 )
  {
    m_domainStress.resize( 3 );
    m_stressControlLastError.resize(3);
    LvArray::tensorOps::fill< 3 >( m_stressControlLastError, 0.0 );
    m_stressControlITerm.resize(3);
    LvArray::tensorOps::fill< 3 >( m_stressControlITerm, 0.0 );
    

    GEOS_ERROR_IF(m_stressTable.size(0) == 0, "Stress table cannot be empty is stress control is enabled");
    GEOS_ERROR_IF(m_stressTable.size(1) == 0, "Stress table must have 4 columns");

    for(int i = 0; i < m_stressTable.size(0) ; i++)
    {
      GEOS_ERROR_IF(m_stressTable[i][0] < 0.0, "Stress table times must be positive");
    }
  }

  // Get nodal position
  arrayView1d< int const > const periodic = partition.getPeriodic();
  int numNodes = nodeManager.size();
  arrayView2d< real64, nodes::REFERENCE_POSITION_USD > const & gridPosition = nodeManager.referencePosition();

  for(int i =0; i < 3; i++)
  {
    if( periodic[i] && (partition.m_coords[i] == 0 || partition.m_coords[i] == partition.m_Partitions[i]-1) )
    {
      real64 xExtent = partition.getGlobalMax()[i] - partition.getGlobalMin()[i];
      for(int g=0; g<nodeManager.size(); g++)
      {
        // if (gridPosition[g][i] < partition.getLocalMin()[i] && gridPosition[g][i] > partition.getLocalMax()[i] ){
          //Partition is on positive face
          if( partition.m_coords[i] == partition.m_Partitions[i]-1) // CC: Does this need to be toleranced?
          {
            if(gridPosition[g][i] < partition.getLocalMin()[i] - xExtent/2)
            {
              gridPosition[g][i] += xExtent; //Do I nee to subtract two cells that are ghost? Shouldn't have those if periodic boundaries are on
            }
          }

          //Partition is on negative face
          if( partition.m_coords[i] == 0){
            if(gridPosition[g][i] > partition.getLocalMax()[i] + xExtent/2) // CC: Does this need to be toleranced?
            {
              gridPosition[g][i] -= xExtent; //Do I nee to subtract two cells that are ghost? Shouldn't have those if periodic boundaries are on
            }
          }
        // }

      }
    }
  }

  // Get local domain extent
  for( int g=0; g<numNodes; g++ )
  {
    for( int i=0; i<3; i++ )
    {
      m_xLocalMin[i] = std::fmin( m_xLocalMin[i], gridPosition[g][i] );
      m_xLocalMax[i] = std::fmax( m_xLocalMax[i], gridPosition[g][i] );
    }
  }
  for( int i=0; i<3; i++ )
  {
    m_xLocalMinNoGhost[i] = partition.getLocalMin()[i];
    m_xLocalMaxNoGhost[i] = partition.getLocalMax()[i];
    m_partitionExtent[i] = m_xLocalMax[i] - m_xLocalMin[i];
  }

  // CC: why not compute element size directly from domain extent and number of cpps across direction?
  // Get element size
  for( int g=0; g<numNodes; g++ )
  {
    for( int i=0; i<3; i++ )
    {
      real64 test = gridPosition[g][i] - m_xLocalMin[i]; // By definition, this should always be positive. Furthermore, the gridPosition
                                                         // should only be those on the local partition
      if( test > 0.0 ) // We're looking for the smallest nonzero distance from the "min" node. TODO: Could be vulnerable to a finite
                       // precision bug.
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
      m_neighborRadius *= -1.0 * sqrt( m_hEl[0] * m_hEl[0] + m_hEl[1] * m_hEl[1] );
    }
    else
    {
      m_neighborRadius *= -1.0 * sqrt( m_hEl[0] * m_hEl[0] + m_hEl[1] * m_hEl[1] + m_hEl[2] * m_hEl[2] );
    }
  }

  // Get global domain extent excluding buffer nodes
  for( int i=0; i<3; i++ )
  {
    m_xGlobalMin[i] = partition.getGlobalMin()[i];
    m_xGlobalMax[i] = partition.getGlobalMax()[i];
    if( !periodic[i] )
    {
      m_xGlobalMin[i] += m_hEl[i];
      m_xGlobalMax[i] -= m_hEl[i];  
    }

    m_domainExtent[i] = m_xGlobalMax[i] - m_xGlobalMin[i];
  }

  // Get number of elements in each direction
  for( int i=0; i<3; i++ )
  {
    m_nEl[i] = std::round( m_partitionExtent[i] / m_hEl[i] );
  }

  // Create element map
  m_ijkMap.resize( m_nEl[0] + 1, m_nEl[1] + 1, m_nEl[2] + 1 );
  for( int g=0; g<numNodes; g++ )
  {
    int i = std::round( ( gridPosition[g][0] - m_xLocalMin[0] ) / m_hEl[0] );
    int j = std::round( ( gridPosition[g][1] - m_xLocalMin[1] ) / m_hEl[1] );
    int k = std::round( ( gridPosition[g][2] - m_xLocalMin[2] ) / m_hEl[2] );
    m_ijkMap[i][j][k] = g;
    
  }

  // Identify node sets for applying boundary conditions. We need boundary nodes and buffer nodes.
  Group & nodeSets = nodeManager.sets();
  array1d< SortedArray< localIndex > > & m_boundaryNodes = nodeSets.getReference< array1d< SortedArray< localIndex > > >( viewKeyStruct::boundaryNodesString() );
  array1d< SortedArray< localIndex > > & m_bufferNodes = nodeSets.getReference< array1d< SortedArray< localIndex > > >( viewKeyStruct::bufferNodesString() );
  m_boundaryNodes.resize( 6 );
  m_bufferNodes.resize( 6 );

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
        tmpBoundaryNodes.insert( g );
      }

      if( sign * positionRelativeToBoundary > 0 && std::fabs( positionRelativeToBoundary ) > tolerance ) // basically a dot product with the
                                                                                                         // face normal
      {
        tmpBufferNodes.insert( g );
      }
    }

    m_boundaryNodes[face].insert( tmpBoundaryNodes.begin(), tmpBoundaryNodes.end() );
    m_bufferNodes[face].insert( tmpBufferNodes.begin(), tmpBufferNodes.end() );
  }

  // Initialize reaction force history file and write its header
  if( MpiWrapper::commRank( MPI_COMM_GEOSX ) == 0 && m_reactionHistory == 1 )
  {
    std::ofstream file;
    file.open( "reactionHistory.csv", std::ios::out); // | std::ios::app );
    if( file.fail() )
      throw std::ios_base::failure( std::strerror( errno ) );
    //make sure write fails with exception if something is wrong
    file.exceptions( file.exceptions() | std::ios::failbit | std::ifstream::badbit );
    file << "time, F00, F11, F22, length_x, length_y, length_z, Rx-, Rx+, Ry-, Ry+, Rz-, Rz+, L00, L11, L22" << std::endl;
    file << std::setprecision( std::numeric_limits< long double >::digits10 )
         << 0.0 << ","
         << 1.0 << "," << 1.0 << "," << 1.0 << ","
         << m_domainExtent[0] << "," << m_domainExtent[1] << "," << m_domainExtent[2] << ","
         << 0.0 << "," << 0.0 << "," << 0.0 << "," << 0.0 << "," << 0.0 << "," << 0.0 << ","
         << 0.0 << "," << 0.0 << "," << 0.0
         << std::endl;
  }

  // Initialize particle fields that weren't intialized by reading the particle input file
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    // Getters
    string const & solidMaterialName = subRegion.template getReference< string >( viewKeyStruct::solidMaterialNamesString() );
    SolidBase & constitutiveRelation = getConstitutiveModel< SolidBase >( subRegion, solidMaterialName );
    arrayView2d< real64 const > const constitutiveDensity = constitutiveRelation.getDensity();

    arrayView2d< real64 const > const particlePosition = subRegion.getParticleCenter();
    arrayView1d< real64 const > const particleVolume = subRegion.getParticleVolume();
    arrayView3d< real64 const > const particleRVectors = subRegion.getParticleRVectors();

    arrayView2d< real64 > const particleDisplacement = subRegion.getParticleDisplacement();
    arrayView2d< real64 > const particleBodyForce = subRegion.getField< fields::mpm::particleBodyForce >();
    arrayView2d< real64 > const particlePlasticStrain = subRegion.getField< fields::mpm::particlePlasticStrain >(); 
    arrayView1d< real64 > const particleDensity = subRegion.getField< fields::mpm::particleDensity >();
    arrayView1d< real64 > const particleMass = subRegion.getField< fields::mpm::particleMass >();
    arrayView3d< real64 > const particleDeformationGradient = subRegion.getField< fields::mpm::particleDeformationGradient >();
    arrayView3d< real64 > const particleFDot = subRegion.getField< fields::mpm::particleFDot >();
    arrayView3d< real64 > const particleVelocityGradient = subRegion.getField< fields::mpm::particleVelocityGradient >();
    arrayView1d< real64 > const particleInitialVolume = subRegion.getField< fields::mpm::particleInitialVolume >();
    arrayView3d< real64 > const particleInitialRVectors = subRegion.getField< fields::mpm::particleInitialRVectors >();
    arrayView3d< real64 > const particleSphF = subRegion.getField< fields::mpm::particleSphF >();
    arrayView2d< real64 > const particleReferencePosition = subRegion.getField< fields::mpm::particleReferencePosition >();
    arrayView1d< int > const particleCrystalHealFlag = subRegion.getField< fields::mpm::particleCrystalHealFlag >();

    // Set reference position, volume and R-vectors
    for( int p=0; p<subRegion.size(); p++ )
    {
      particleInitialVolume[p] = particleVolume[p];
      particleCrystalHealFlag[p] = 0;
      for( int i=0; i<3; i++ )
      {
        particleBodyForce[p][i] = 0;
        particleDisplacement[p][i] = 0;
        particleReferencePosition[p][i] = particlePosition[p][i];
        for( int j=0; j<3; j++ )
        {
          particleInitialRVectors[p][i][j] = particleRVectors[p][i][j];
        }
      }
    }

    // Pull initial density from constitutive model, set particle masses and small mass threshold
    real64 localMinMass = 0.0;
    real64 globalMinMass;
    for( int p=0; p<subRegion.size(); p++ )
    {
      particleDensity[p] = constitutiveDensity[p][0];
      particleMass[p] = particleDensity[p] * particleVolume[p];
      localMinMass = particleMass[p] < localMinMass ? particleMass[p] : localMinMass;

      for (int j=0; j < 6; ++j){
        particlePlasticStrain[p][j] = 0;
      }
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

    // Initialize deformation gradient and velocity gradient
    for( int p=0; p<subRegion.size(); p++ )
    {
      for( int i=0; i<3; i++ )
      {
        for( int j=0; j<3; j++ )
        {
          particleDeformationGradient[p][i][j] = i == j ? 1.0 : 0.0;
          particleSphF[p][i][j] = i == j ? 1.0 : 0.0;
          particleFDot[p][i][j] = 0.0;
          particleVelocityGradient[p][i][j] = 0.0;
        }
      }
    }

  } );

  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    subRegion.setActiveParticleIndices(); // Needed for computeAndWriteBoxAverage().
  } );

  // Initialize box average history file and write its header
  if( m_boxAverageHistory == 1 )
  {

    if( MpiWrapper::commRank( MPI_COMM_GEOSX ) == 0 )
    {
      std::ofstream file;
      file.open( "boxAverageHistory.csv", std::ios::out );
      if( file.fail() )
      {
        throw std::ios_base::failure( std::strerror( errno ) );
      }
      file.exceptions( file.exceptions() | std::ios::failbit | std::ifstream::badbit );
      file << "Time, Sxx, Syy, Szz, Syz, Sxz, Sxy, Density, Damage, epxx, epyy, epzz, epyz, epxz, epxy" << std::endl;
    }
    MpiWrapper::barrier( MPI_COMM_GEOSX ); // wait for the header to be written

    // Initialize box average extent
    if(  m_boxAverageMin.size() == 0  )
    {
      m_boxAverageMin.resize( 3 );
      LvArray::tensorOps::copy< 3 >( m_boxAverageMin, m_xGlobalMin );
    }
    else
    {
      GEOS_ERROR_IF( m_boxAverageMin.size() > 3, "Box average min must have 3 elements" );
    }

    if(  m_boxAverageMax.size() == 0  )
    {
      m_boxAverageMax.resize( 3 );
      LvArray::tensorOps::copy< 3 >( m_boxAverageMax, m_xGlobalMax );
    }
    else
    {
      GEOS_ERROR_IF( m_boxAverageMax.size() > 3, "Box average min must have 3 elements" );
    }

    GEOS_ERROR_IF( ( m_boxAverageMin[0] > m_boxAverageMax[0] ) || ( m_boxAverageMin[1] > m_boxAverageMax[1] ) || ( m_boxAverageMin[2] > m_boxAverageMax[2] ) , "Box minimums must be less than box maximums");

    computeAndWriteBoxAverage( 0.0, 0.0, particleManager );
  }

  // Resize grid arrays according to number of velocity fields
  int maxLocalGroupNumber = 0; // Maximum contact group number on this partition.
  int maxGlobalGroupNumber; // Maximum contact group number on global domain.
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    arrayView1d< int > const particleGroup = subRegion.getParticleGroup();
    for( int p=0; p<subRegion.size(); p++ )
    {
      if( particleGroup[p] > maxLocalGroupNumber )
      {
        maxLocalGroupNumber = particleGroup[p];
      }
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
  nodeManager.getReference< array2d< real64 > >( viewKeyStruct::massString() ).resize( numNodes, m_numVelocityFields );
  nodeManager.getReference< array2d< real64 > >( viewKeyStruct::damageString() ).resize( numNodes, m_numVelocityFields );
  nodeManager.getReference< array2d< real64 > >( viewKeyStruct::maxDamageString() ).resize( numNodes, m_numVelocityFields );
  nodeManager.getReference< array2d< real64 > >( viewKeyStruct::damageGradientString() ).resize( numNodes, 3 );
  nodeManager.getReference< array3d< real64 > >( viewKeyStruct::velocityString() ).resize( numNodes, m_numVelocityFields, 3 );
  nodeManager.getReference< array3d< real64 > >( viewKeyStruct::momentumString() ).resize( numNodes, m_numVelocityFields, 3 );
  nodeManager.getReference< array3d< real64 > >( viewKeyStruct::accelerationString() ).resize( numNodes, m_numVelocityFields, 3 );
  nodeManager.getReference< array3d< real64 > >( viewKeyStruct::forceInternalString() ).resize( numNodes, m_numVelocityFields, 3 );
  nodeManager.getReference< array3d< real64 > >( viewKeyStruct::forceExternalString() ).resize( numNodes, m_numVelocityFields, 3 );
  nodeManager.getReference< array3d< real64 > >( viewKeyStruct::forceContactString() ).resize( numNodes, m_numVelocityFields, 3 );
  nodeManager.getReference< array3d< real64 > >( viewKeyStruct::surfaceNormalString() ).resize( numNodes, m_numVelocityFields, 3 );
  nodeManager.getReference< array3d< real64 > >( viewKeyStruct::materialPositionString() ).resize( numNodes, m_numVelocityFields, 3 );

  // Load strength scale into constitutive model (for ceramic)
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    // Get constitutive model reference
    string const & solidMaterialName = subRegion.template getReference< string >( viewKeyStruct::solidMaterialNamesString() );
    SolidBase & solidModel = getConstitutiveModel< SolidBase >( subRegion, solidMaterialName );
    if( solidModel.hasWrapper( "strengthScale" ) )
    {
      SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
      arrayView1d< real64 > const particleStrengthScale = subRegion.getParticleStrengthScale();
      arrayView1d< real64 > const constitutiveStrengthScale = solidModel.getReference< array1d< real64 > >( "strengthScale" );
      forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
      {
        localIndex const p = activeParticleIndices[pp];
        constitutiveStrengthScale[p] = particleStrengthScale[p];
      } );
    }
  } );

  // Initialize friction coefficient table
  initializeFrictionCoefficients();
}

real64 SolidMechanicsMPM::explicitStep( real64 const & time_n,
                                        real64 const & dt,
                                        const int cycleNumber,
                                        DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  #define USE_PHYSICS_LOOP

  //#######################################################################################
  GEOS_LOG_RANK_IF( m_debugFlag == 1, "Get spatial partition, get node and particle managers. Resize m_iComm." );
  solverProfiling( "Get spatial partition, get node and particle managers. Resize m_iComm." );
  //#######################################################################################
  SpatialPartition & partition = dynamic_cast< SpatialPartition & >( domain.getPartition() );
  arrayView1d< int const > const periodic = partition.getPeriodic();

  // ***** We assume that there are exactly two mesh bodies, and that one has particles and one does not. *****
  Group & meshBodies = domain.getMeshBodies();

  MeshBody & meshBody1 = meshBodies.getGroup< MeshBody >( 0 );
  MeshBody & meshBody2 = meshBodies.getGroup< MeshBody >( 1 );
  MeshBody & particles = meshBody1.hasParticles() ? meshBody1 : meshBody2;
  MeshBody & grid = !meshBody1.hasParticles() ? meshBody1 : meshBody2;

  ParticleManager & particleManager = particles.getBaseDiscretization().getParticleManager();
  MeshLevel & mesh = grid.getBaseDiscretization();
  NodeManager & nodeManager = mesh.getNodeManager();

  m_iComm.resize( domain.getNeighbors().size() );


  //#######################################################################################
  GEOS_LOG_RANK_IF( m_debugFlag == 1 && cycleNumber == 0 , "At time step zero, perform initialization calculations" );
  solverProfilingIf( "At time step zero, perform initialization calculations", cycleNumber == 0 );
  //#######################################################################################
  if( cycleNumber == 0 )
  {
    initialize( nodeManager, particleManager, partition );
  }

  //#######################################################################################
  GEOS_LOG_RANK_IF( m_debugFlag == 1, "Set grid multi-field labels to avoid a VTK output bug" );
  solverProfiling( "Set grid multi-field labels to avoid a VTK output bug" );
  //#######################################################################################
  // Must be done every time step despite grid fields being registered
  // TODO: Only doing this on the 1st cycle breaks restarts, can we do better? Why aren't labels part of the restart data?
  setGridFieldLabels( nodeManager );


  //#######################################################################################
  GEOS_LOG_RANK_IF( m_debugFlag == 1, "Set grid fields to zero" );
  solverProfiling( "Set grid fields to zero" );
  //#######################################################################################
  initializeGridFields( nodeManager );


  //#######################################################################################
  GEOS_LOG_RANK_IF( m_debugFlag == 1 && m_useEvents == 1, "Check if event has been triggered" );
  solverProfilingIf( "Check if event has been triggered", m_useEvents == 1);
  //#######################################################################################
  if( m_useEvents == 1 )
  {
    triggerEvents( dt,
                   time_n,
                   particleManager,
                   partition );
  }

  //#######################################################################################
  GEOS_LOG_RANK_IF( m_debugFlag == 1, "Update global-to-local map" );
  solverProfiling( "Update global-to-local map" );
  //#######################################################################################
  particleManager.updateMaps();

  //#######################################################################################
  GEOS_LOG_RANK_IF( m_debugFlag == 1 && MpiWrapper::commSize( MPI_COMM_GEOSX ) > 1 && m_needsNeighborList == 1,  "Perform particle ghosting" );
  solverProfilingIf( "Perform particle ghosting", MpiWrapper::commSize( MPI_COMM_GEOSX ) > 1 && m_needsNeighborList == 1 );
  //#######################################################################################
  if( MpiWrapper::commSize( MPI_COMM_GEOSX ) > 1 && m_needsNeighborList == 1 )
  {
    // Move everything into host memory space
    particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
    {
      subRegion.forWrappers( [&]( WrapperBase & wrapper )
      {
        wrapper.move( LvArray::MemorySpace::host, true );
      } );
    } );

    // Perform ghosting
    partition.getGhostParticlesFromNeighboringPartitions( domain, m_iComm, m_neighborRadius );
  }

  //#######################################################################################
  GEOS_LOG_RANK_IF(m_debugFlag == 1, "Get indices of non-ghost particles" );
  solverProfiling( "Get indices of non-ghost particles" );
  //#######################################################################################
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    subRegion.setActiveParticleIndices();
  } );

  //#######################################################################################
  GEOS_LOG_RANK_IF(m_debugFlag == 1 && ( periodic[0] == 1 || periodic[1] == 1 || periodic[2] == 1 ), "Correct ghost particle centers across periodic boundaries" );
  solverProfilingIf( "Correct ghost particle centers across periodic boundaries",  periodic[0] == 1 || periodic[1] == 1 || periodic[2] == 1 );
  //#######################################################################################
  if( periodic[0] == 1 || periodic[1] == 1 || periodic[2] == 1 )
  {
    correctGhostParticleCentersAcrossPeriodicBoundaries(particleManager, partition);
  }

  //#######################################################################################
  GEOS_LOG_RANK_IF(m_debugFlag == 1 && m_needsNeighborList == 1, "Construct neighbor list" );
  solverProfilingIf( "Construct neighbor list", m_needsNeighborList == 1 );
  //#######################################################################################
  if( m_needsNeighborList == 1 )
  { // This optimization compares neighbor list creation time for different
    // bin sizes, and finds an optimum.  It can be non deterministic and
    // is currently disabled to ensure integrated tests pass.
    // Optimize
    // if( cycleNumber == 0 )
    // {
    //   optimizeBinSort( particleManager );
    // }
    // else
    // {
    (void) computeNeighborList( particleManager );
    // }
  }

  //#######################################################################################
  GEOS_LOG_RANK_IF( m_debugFlag == 1 && m_surfaceDetection > 0 && cycleNumber == 0, "Compute surface flags" );
  solverProfilingIf( "Compute surface flags", m_surfaceDetection > 0 && cycleNumber == 0 );
  //#######################################################################################
  if( m_surfaceDetection > 0 && ( cycleNumber == 0 || m_surfaceHealing == true ) )
  {
    m_surfaceHealing = false;
    computeSurfaceFlags( particleManager );
  }

  //#######################################################################################
  GEOS_LOG_RANK_IF( m_debugFlag == 1 && m_damageFieldPartitioning == 1, "Compute damage field gradient" );
  solverProfilingIf( "Compute damage field gradient", m_damageFieldPartitioning == 1 );
  //#######################################################################################
  if( m_damageFieldPartitioning == 1 )
  {
    computeDamageFieldGradient( particleManager );
  }

  //#######################################################################################
  GEOS_LOG_RANK_IF( m_debugFlag == 1, "Update surface flag overload" );
  solverProfiling( "Update surface flag overload" );
  //#######################################################################################
  // We now read in surface flags so don't need to update them based on damage.
  // Keeping this here as it will eventually be used to support other surface
  // flag usages.
  //updateSurfaceFlagOverload( particleManager );

  //#######################################################################################
  GEOS_LOG_RANK_IF( m_debugFlag == 1 && m_prescribedBcTable == 1, "Update BCs based on bcTable" );
  solverProfilingIf( "Update BCs based on bcTable", m_prescribedBcTable == 1 );
  //#######################################################################################
  if( m_prescribedBcTable == 1 )
  {
    boundaryConditionUpdate( dt, time_n );
  }

  //#######################################################################################
  GEOS_LOG_RANK_IF( m_debugFlag == 1 && m_cpdiDomainScaling == 1, "Perform r-vector scaling (CPDI domain scaling)" );
  solverProfilingIf( "Perform r-vector scaling (CPDI domain scaling)", m_cpdiDomainScaling == 1 );
  //#######################################################################################
  if( m_cpdiDomainScaling == 1 )
  {
    cpdiDomainScaling( particleManager );
  }

  //#######################################################################################
  GEOS_LOG_RANK_IF( m_debugFlag == 1, "Resize and populate mapping arrays" );
  solverProfiling( "Resize and populate mapping arrays" );
  //#######################################################################################
  resizeMappingArrays( particleManager );
  populateMappingArrays( particleManager, nodeManager );

  //#######################################################################################
  GEOS_LOG_RANK_IF( m_debugFlag == 1 && m_damageFieldPartitioning == 1, "Project damage field gradient to the grid and then sync" );
  solverProfilingIf( "Project damage field gradient to the grid and then sync", m_damageFieldPartitioning == 1 );
  //#######################################################################################
  if( m_damageFieldPartitioning == 1 )
  {
    projectDamageFieldGradientToGrid( particleManager, nodeManager );
    std::vector< std::string > fieldNames = { viewKeyStruct::damageGradientString() };
    syncGridFields( fieldNames, domain, nodeManager, mesh, MPI_MAX );
  }

  //#######################################################################################
  GEOS_LOG_RANK_IF( m_debugFlag == 1, "Compute particle body forces" );
  solverProfiling( "Compute particle body forces" );
  //#######################################################################################
  computeBodyForce( particleManager);

  //#######################################################################################
  GEOS_LOG_RANK_IF( m_debugFlag == 1 && m_generalizedVortexMMS == 1, "Compute particle body forces" );
  solverProfilingIf( "Compute particle body forces", m_generalizedVortexMMS == 1 );
  //#######################################################################################
  if( m_generalizedVortexMMS == 1 ){
    computeGeneralizedVortexMMSBodyForce( time_n,
                                          particleManager );
  } 

  //#######################################################################################
  GEOS_LOG_RANK_IF( m_debugFlag == 1, "Particle-to-grid interpolation" );
  solverProfiling( "Particle-to-grid interpolation" );
  //#######################################################################################
  particleToGrid( particleManager, nodeManager );

  //#######################################################################################
  GEOS_LOG_RANK_IF( m_debugFlag == 1, "Grid MPI operations" );
  solverProfiling( "Grid MPI operations" );
  //#######################################################################################
  std::vector< std::string > fieldNames1 = { viewKeyStruct::massString(),
                                             viewKeyStruct::momentumString(),
                                             viewKeyStruct::damageString(),
                                             viewKeyStruct::materialPositionString(),
                                             viewKeyStruct::forceInternalString(),
                                             viewKeyStruct::forceExternalString() };
  syncGridFields( fieldNames1, domain, nodeManager, mesh, MPI_SUM );
  std::vector< std::string > fieldNames2 = { viewKeyStruct::maxDamageString() };
  syncGridFields( fieldNames2, domain, nodeManager, mesh, MPI_MAX );

  //#######################################################################################
  GEOS_LOG_RANK_IF( m_debugFlag == 1, "Determine trial momenta and velocities based on acceleration due to internal and external forces, but before contact enforcement" );
  solverProfiling( "Determine trial momenta and velocities based on acceleration due to internal and external forces, but before contact enforcement" );
  //#######################################################################################
  gridTrialUpdate( dt, nodeManager );

  //#######################################################################################
  GEOS_LOG_RANK_IF( m_debugFlag == 1 && m_numVelocityFields > 1, "Contact enforcement" );
  solverProfilingIf( "Contact enforcement", m_numVelocityFields > 1 );
  //#######################################################################################
  if( m_numVelocityFields > 1 )
  {
    enforceContact( dt, domain, particleManager, nodeManager, mesh );
  }

  //#######################################################################################
  GEOS_LOG_RANK_IF( m_debugFlag == 1 && ( m_stressControl[0] == 1 || m_stressControl[1] == 1 || m_stressControl[2] == 1 ), "Interpolate stress table" );
  solverProfilingIf( "Interpolate stress table", m_stressControl[0] == 1 || m_stressControl[1] == 1 || m_stressControl[2] == 1 );
  //#######################################################################################
  if( m_stressControl[0] == 1 || m_stressControl[1] == 1 || m_stressControl[2] == 1 )
  {
    interpolateStressTable( dt, time_n );
  }

  //#######################################################################################
  GEOS_LOG_RANK_IF( m_debugFlag == 1 && ( ( !m_stressControl[0] || !m_stressControl[1] || !m_stressControl[2] ) && ( m_prescribedBoundaryFTable == 1 || m_prescribedFTable == 1 ) ), "Interpolate F table" );
  solverProfilingIf( "Interpolate F table", ( !m_stressControl[0] || !m_stressControl[1] || !m_stressControl[2] ) && ( m_prescribedBoundaryFTable == 1 || m_prescribedFTable == 1 ) );
  //#######################################################################################
  if( ( !m_stressControl[0] || !m_stressControl[1] || !m_stressControl[2] ) && ( m_prescribedBoundaryFTable == 1 || m_prescribedFTable == 1 ) )

  {
    interpolateFTable( dt, time_n );
  }

  //#######################################################################################
  GEOS_LOG_RANK_IF( m_debugFlag == 1, "Apply essential boundary conditions" );
  solverProfiling( "Apply essential boundary conditions" );
  //#######################################################################################
  applyEssentialBCs( dt, time_n, nodeManager );


  // //#######################################################################################
  // solverProfilingIf( "Directional overlap correction", m_directionalOverlapCorrection == 1 );
  // TODO: Figure out syncing on ghost particles and move this to after the F update
  // //#######################################################################################
  // if( m_directionalOverlapCorrection == 1 )
  // {
  //   directionalOverlapCorrection( dt, particleManager );
  // }


  //#######################################################################################
  GEOS_LOG_RANK_IF( m_debugFlag == 1, "Grid-to-particle interpolation" );
  solverProfiling( "Grid-to-particle interpolation" );
  //#######################################################################################
  gridToParticle( dt, particleManager, nodeManager );

  //#######################################################################################
  GEOS_LOG_RANK_IF( m_debugFlag == 1 && m_prescribedFTable == 1, "Update particle positions according to prescribed F Table" );
  solverProfilingIf( "Update particle positions according to prescribed F Table", m_prescribedFTable == 1 );
  //####################################################################################### 
  if( m_prescribedFTable == 1 )
  {
    applySuperimposedVelocityGradient( dt,
                                       particleManager,
                                       partition );
  }

  //#######################################################################################
  GEOS_LOG_RANK_IF( m_debugFlag == 1, "Update deformation gradient" );
  solverProfiling( "Update deformation gradient" );
  //#######################################################################################
  updateDeformationGradient( dt, particleManager );


  //#######################################################################################
  GEOS_LOG_RANK_IF( m_debugFlag == 1, "Update particle geometry (e.g. volume, r-vectors) and density" );
  solverProfiling( "Update particle geometry (e.g. volume, r-vectors) and density" );
  //#######################################################################################
  particleKinematicUpdate( particleManager );


  //#######################################################################################
  GEOS_LOG_RANK_IF( m_debugFlag == 1, "Update constitutive model dependencies" );
  solverProfiling( "Update constitutive model dependencies" );
  //#######################################################################################
  updateConstitutiveModelDependencies( particleManager );


  //#######################################################################################
  GEOS_LOG_RANK_IF( m_debugFlag == 1, "Update stress" );
  solverProfiling( "Update stress" );
  //#######################################################################################
  updateStress( dt, particleManager );


  //#######################################################################################
  GEOS_LOG_RANK_IF( m_debugFlag == 1, "Update solver dependencies" );
  solverProfiling( "Update solver dependencies" );
  //#######################################################################################
  updateSolverDependencies( particleManager );


  //#######################################################################################
  GEOS_LOG_RANK_IF( m_debugFlag == 1 && m_boxAverageHistory == 1, "Compute and write box averages" );
  solverProfilingIf( "Compute and write box averages", m_boxAverageHistory == 1 );
  //#######################################################################################
  if( m_boxAverageHistory == 1 )
  {
    if( time_n + dt >= m_nextBoxAverageWriteTime )
    {
      computeAndWriteBoxAverage( time_n, dt, particleManager );
      m_nextBoxAverageWriteTime += m_boxAverageWriteInterval;
    }
  }

  //#######################################################################################
  GEOS_LOG_RANK_IF( m_debugFlag == 1 && ( m_stressControl[0] == 1 || m_stressControl[1] == 1 || m_stressControl[2] == 1 ), "Stress control" );
  solverProfilingIf("Stress control",  m_stressControl[0] == 1 || m_stressControl[1] == 1 || m_stressControl[2] == 1 );
  //#######################################################################################
  // stress control reads a principal stress table.  we compute the
  // difference between most recent box sum and the prescribed value and increment
  // the domain strain rate accordingly.  This is a simple control loop.
  if( m_stressControl[0] == 1 || m_stressControl[1] == 1 || m_stressControl[2] == 1 )
  {
    stressControl( dt,
                   particleManager );
  }


  //#######################################################################################
  GEOS_LOG_RANK_IF( m_debugFlag == 1, "Calculate stable time step" );
  solverProfiling( "Calculate stable time step" );
  //#######################################################################################
  real64 dtReturn = getStableTimeStep( particleManager );

  //#######################################################################################
  GEOS_LOG_RANK_IF( m_debugFlag == 1, "Update global-to-local map" );
  solverProfiling( "Update global-to-local map" );
  //#######################################################################################
  flagOutOfRangeParticles( particleManager );

  //#######################################################################################
  GEOS_LOG_RANK_IF( m_debugFlag == 1, "Delete bad particles" );
  solverProfiling( "Delete bad particles" );
  //#######################################################################################
  deleteBadParticles( particleManager );


  //#######################################################################################
  GEOS_LOG_RANK_IF( m_debugFlag == 1 && MpiWrapper::commSize( MPI_COMM_GEOSX ) > 1, "Delete bad particles" );
  solverProfilingIf( "Particle repartitioning", MpiWrapper::commSize( MPI_COMM_GEOSX ) > 1 );
  //#######################################################################################
  if( MpiWrapper::commSize( MPI_COMM_GEOSX ) > 1 )
  {
    particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
    {
      subRegion.forWrappers( [&]( WrapperBase & wrapper )
      {
        wrapper.move( LvArray::MemorySpace::host, true );
      } );
      partition.repartitionMasterParticles( subRegion, m_iComm );

      // CC: We need to know which of the particles are master to correct the particle centers across periodic boundaries
      // We could perform this correction on all particles, but that might be slower since we would also be iterating over
      // Ghost particles that get update beginning of next time step anyways
      subRegion.setActiveParticleIndices();
    } );

    //CC: Correct particle centers across periodic boundaries
    if( periodic[0] == 1 || periodic[1] == 1 || periodic[2] == 1 )
    {
        correctParticleCentersAcrossPeriodicBoundaries(particleManager, partition);
    }
  }

  //#######################################################################################
  GEOS_LOG_RANK_IF( m_debugFlag == 1 && ( m_prescribedBoundaryFTable == 1 || m_prescribedFTable == 1 || m_stressControl[0] == 1 || m_stressControl[1] == 1 || m_stressControl[2] == 1 ), "Resize grid based on F-table" );
  solverProfilingIf( "Resize grid based on F-table",  m_prescribedBoundaryFTable == 1 || m_prescribedFTable == 1 || m_stressControl[0] == 1 || m_stressControl[1] == 1 || m_stressControl[2] == 1 );
  //#######################################################################################
  if( m_prescribedBoundaryFTable == 1 || m_prescribedFTable == 1 || m_stressControl[0] == 1 || m_stressControl[1] == 1 || m_stressControl[2] == 1 )
  {
    resizeGrid( partition, nodeManager, dt );
  }


  //#######################################################################################
  GEOS_LOG_RANK_IF( m_debugFlag == 1, "End of explicitStep");
  solverProfiling( "End of explicitStep" );
  //#######################################################################################
  if( m_solverProfiling >= 1 )
  {
    printProfilingResults();
  }


  // Return stable time step
  return dtReturn;
}

void SolidMechanicsMPM::triggerEvents( const real64 dt,
                                       const real64 time_n,
                                       ParticleManager & particleManager,
                                       SpatialPartition & partition )
{
  //Iterate over every MPM Events to check if conditions are met and if so perform event
  m_mpmEventManager->forSubGroups< MPMEventBase >( [&]( MPMEventBase & event )
  {
    real64 eventTime = event.getTime();
    real64 eventInterval = event.getInterval();

    if( ( eventTime - dt / 2 <= time_n && time_n <= eventTime + eventInterval + dt/2 ) && !event.isComplete() )
    {
      if( event.getName() == "MaterialSwap" )
      {
        MaterialSwapMPMEvent & materialSwap = dynamicCast< MaterialSwapMPMEvent & >( event );
        GEOS_LOG_RANK_0("Performing material swap");
        performMaterialSwap( particleManager, materialSwap.getSourceRegion(), materialSwap.getDestinationRegion() );
        event.setIsComplete( 1 );
      }

      if( event.getName() == "Anneal" )
      {
        GEOS_LOG_RANK_0( "Processing anneal event");
        AnnealMPMEvent & anneal = dynamicCast< AnnealMPMEvent & >( event );

        particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
        {
          if( subRegion.getName() == anneal.getTargetRegion() || anneal.getTargetRegion() == "all" )
          {
            // Get constitutive model reference
            string const & solidMaterialName = subRegion.template getReference< string >( viewKeyStruct::solidMaterialNamesString() );
            SolidBase & solidModel = getConstitutiveModel< SolidBase >( subRegion, solidMaterialName );

            GEOS_ERROR_IF( !solidModel.hasWrapper( "oldStress" ), "Cannot anneal constitutive model that does not have oldStress wrapper!");
            arrayView3d< real64 > const constitutiveOldStress = solidModel.getReference< array3d< real64 > >( "oldStress" );
          
            // SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
            forAll< serialPolicy >( constitutiveOldStress.size(0), [=] GEOS_HOST ( localIndex const p )
            {

                real64 knockdown = 0.0; // At the end of the annealing process, strongly enforce zero deviatoric stress.
                if( !( time_n - dt / 2 <= eventTime + eventInterval && time_n + dt / 2 > eventTime + eventInterval ) )
                {
                  // Otherwise, scale it down by a number that gradually goes from 1 to 0
                  knockdown = fmax( 0.0, 1.0 - dt * 20.0 * ( time_n - eventTime ) / ( eventInterval * eventInterval ) ); // Gaussian decay
                }

                // Smoothly knock down the deviatoric stress to simulate annealing
                real64 stress[6] = {0};
                LvArray::tensorOps::copy< 6 >( stress, constitutiveOldStress[p][0]);

                real64 trialP;
                real64 trialQ;
                real64 deviator[6];
                twoInvariant::stressDecomposition( stress,
                                                   trialP,
                                                   trialQ,
                                                   deviator );

                LvArray::tensorOps::copy< 6 >( constitutiveOldStress[p][0], stress );
                twoInvariant::stressRecomposition( trialP,
                                                   knockdown * trialQ,
                                                   deviator,
                                                   stress );
                LvArray::tensorOps::copy< 6 >( constitutiveOldStress[p][0], stress );
            });
          }
        });
      }

      if( event.getName() == "Heal" )
      {
        // HealMPMEvent & heal = dynamicCast< HealMPMEvent & >( event );      
        
        m_surfaceHealing = true;
        particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
        {
          arrayView1d< int > const particleSurfaceFlag = subRegion.getParticleSurfaceFlag();
          
          SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
          forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST ( localIndex const pp )
          {
            localIndex const p = activeParticleIndices[pp];
            particleSurfaceFlag[p] = 0;
          });
        });
        event.setIsComplete( 1 );
      }

      if( event.getName() == "InsertPeriodicContactSurfaces" )
      {
          // InsertPeriodicContactSurfacesMPMEvent & insertPeriodicContactSurfaces = dynamicCast< InsertPeriodicContactSurfacesMPMEvent & >( event );

          real64 xGlobalMin[3];
          real64 xGlobalMax[3];
          real64 hEl[3];
          LvArray::tensorOps::copy< 3 >(xGlobalMin, m_xGlobalMin);
          LvArray::tensorOps::copy< 3 >(xGlobalMax, m_xGlobalMax);
          LvArray::tensorOps::copy< 3 >(hEl, m_hEl);

          // Perhaps I need to find the largest particle size during initialize 
          auto periodic = partition.getPeriodic();

          particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
          {
            arrayView2d< real64 const > const particlePosition = subRegion.getParticleCenter();
            arrayView1d< int > const particleSurfaceFlag = subRegion.getParticleSurfaceFlag();

            SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
            forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST ( localIndex const pp )
            {
              localIndex const p = activeParticleIndices[pp];
              for(int i =0; i < 3; ++i)
              {
                if( periodic[i] && ( ( particlePosition[p][i] - xGlobalMin[i] < hEl[i] ) || (  xGlobalMax[i] - particlePosition[p][i] < hEl[i] ) ) )
                {
                  particleSurfaceFlag[p] = 2;
                  break;
                }
              }
            });
          });
          event.setIsComplete( 1 );
      }

      // if( event.getName() == "CrystalHeal" ){
      //   CrystalHealMPMEvent & crystalHeal = dynamicCast< CrystalHealMPMEvent & >( event );

      //   int healType = crystalHeal.getHealType();
      //   ParticleRegion & targetParticleRegion = particleManager.getRegion< ParticleRegion >( crystalHeal.getTargetRegion() );
        
      //   // Copy particle data from source sub region to destination sub region
      //   auto & targetSubRegions = targetParticleRegion.getSubRegions();
      //   for( int r=0; r < targetSubRegions.size(); ++r)
      //   {
      //     ParticleSubRegion & targetSubRegion = dynamicCast< ParticleSubRegion & >( *targetSubRegions[r] );
    
      //     arrayView2d< real64 > const particleStress = targetSubRegion.getField< fields::mpm::particleStress >();
      //     arrayView1d< real64 > const particleCrystalHealFlag = targetSubRegion.getField< fields::mpm::particleStress >();
      //     arrayView1d< real64 > const particleDamage = subRegion.getParticleDamage();
      //     arrayView3d< real64 > const particleRVectors = subRegion.getParticleRVectors();
      //     arrayView2d< real64 const > const particleDamageGradient = targetSubRegion.getField< fields::mpm::particleDamageGradient >();
      //     arrayView3d< real64 const > const particleDeformationGradient = targetSubRegion.getField< fields::mpm::particleDeformationGradient >();
      //     arrayView1d< real64 > const sourceParticleInitialVolume = sourceSubRegion.getField< fields::mpm::particleInitialVolume >();

      //     // CC: TODO Need to add this for Mike
      //     // arrayView1d< real64 const > const particleReferencePorosity = targetSubRegion.getField< fields::mpm::particleReferencePorosity >();

      //     SortedArrayView< localIndex const > const activeParticleIndices = targetSubRegion.activeParticleIndices();
      //     forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST ( localIndex const pp )
      //     {
      //       localIndex const p = activeParticleIndices[pp];

      //       real64 temp[3] = { 0 };
      //       LvArray::tensorOps::Ri_eq_symAijBj< 3 >( temp, stress, particleDamageGradient[p] );
      //       real64 normalStress = LvArray::tensorOps::AiBi< 3 >( particleDamageGradient[p], temp );

      //       real64 detF = LvArray::tensorOps::determinant< 3 >( particleDeformationGradient[p] );
      //       if( ( healType == 1 || ( healType == 0 && ( normalStress || detF < 1.0 ) ) ) && particleDamage[p] > 0.0 )
      //       {
      //         particleCrystalHealFlag[p] = 1;

      //         if( detF > 1.0 )
      //         {
      //         // particleReferencePorosity[p] = 1.0 - 1.0 / detF;
                
      //           real64 power = m_planeStrain ? 0.5 : 1.0 / 3.0;
      //           real64 scaling = std::powd( detF, power );
      //           LvArray::tensorOps::scale< 3, 3 >( particleDeformationGradient[p], 1/scaling );

      //           LvArray::tensorOps::scale< 3 >( particleRVectors[p][0], scaling );
      //           LvArray::tensorOps::scale< 3 >( particleRVectors[p][1], scaling );

      //           if( !m_planeStrain )
      //           {
      //             LvArray::tensorOps::scale< 3 >( particleRVectors[p][2], scaling );
      //           }

      //           particleInitialVolume[p] *= detF;
      //         }

      //       }
      //     } );

      //   }

      // }

      if( event.getName() == "MachineSample" )
      {
          MachineSampleMPMEvent & machineSample = dynamicCast< MachineSampleMPMEvent & >( event );

          real64 xGlobalMin[3];
          real64 xGlobalMax[3];
          real64 domainExtent[3];
          real64 hEl[3];
          LvArray::tensorOps::copy< 3 >(xGlobalMin, m_xGlobalMin);
          LvArray::tensorOps::copy< 3 >(xGlobalMax, m_xGlobalMax);
          LvArray::tensorOps::copy< 3 >(domainExtent, m_domainExtent);
          LvArray::tensorOps::copy< 3 >(hEl, m_hEl);

          particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
          {
            arrayView2d< real64 const > const particlePosition = subRegion.getParticleCenter();
            arrayView1d< int > const particleIsBad = subRegion.getField< fields::mpm::isBad >();

            SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();

            string sampleType = machineSample.getSampleType();

            real64 gaugeRadius = machineSample.getGaugeRadius();
            real64 gaugeLength = machineSample.getGaugeLength();
            real64 filletRadius = machineSample.getFilletRadius();     
            real64 machiningLength = 2*filletRadius + gaugeLength;
            
            real64 diskRadius = machineSample.getDiskRadius();

            real64 domainCenter[3] = { 0 };
            LvArray::tensorOps::copy< 3 >( domainCenter, xGlobalMin);
            LvArray::tensorOps::add< 3 >( domainCenter, xGlobalMax);
            LvArray::tensorOps::scale< 3 >( domainCenter, 0.5 );
            
            forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST ( localIndex const pp )
            {
              localIndex const p = activeParticleIndices[pp];

              // Dogbone is milled along y direction
              if( sampleType == "dogbone")
              {
                // Particle is inside gauge section
                if( particlePosition[p][1] > xGlobalMin[1] + (domainExtent[1] - machiningLength)/2 && 
                    particlePosition[p][1] < xGlobalMax[1] - (domainExtent[1] - machiningLength)/2 )
                {
                  real64 distSqr = std::pow( particlePosition[p][0] - domainCenter[0], 2 ) + std::pow( particlePosition[p][2] - domainCenter[2], 2 );
                  //Is particle inside the gauge section
                  if( particlePosition[p][1] > xGlobalMin[1] + (domainExtent[1] - gaugeLength)/2 && 
                      particlePosition[p][1] < xGlobalMax[1] - (domainExtent[1] - gaugeLength)/2 )
                  {
                    // Check distance in XZ plane from radius of gauge
                    if( distSqr  > gaugeRadius * gaugeRadius )
                    {
                      particleIsBad[p] = 1;
                    }
                  } 
                  //Check if particle is outside of fillet radius
                  else
                  {
                    real64 y = std::abs( particlePosition[p][1] - domainCenter[1] ) - gaugeLength/2;
                    real64 rr = filletRadius - std::sqrt( ( filletRadius * filletRadius - y * y ) ) + gaugeRadius;
                    if( distSqr > rr * rr)
                    {
                      GEOS_LOG_RANK( p << ": " << std::sqrt(distSqr) << ", Pos: " << particlePosition[p] );
                      particleIsBad[p] = 1;
                    }
                  }
                }
              }

              if( sampleType == "brazilDisk")
              {
                real64 distSqr = std::pow( particlePosition[p][0] - domainCenter[0], 2 ) + std::pow( particlePosition[p][1] - domainCenter[1], 2 )  + std::pow( particlePosition[p][2] - domainCenter[2], 2 );
                  
                if( distSqr > diskRadius * diskRadius )
                {
                  particleIsBad[p] = 1;
                }
              }

              if( sampleType == "cylinder" )
              {
                real64 distSqr = std::pow( particlePosition[p][0] - domainCenter[0], 2 ) + std::pow( particlePosition[p][2] - domainCenter[2], 2 );
                
                if( distSqr  > gaugeRadius * gaugeRadius )
                {
                  particleIsBad[p] = 1;
                }

              }

            });
          });
          event.setIsComplete( 1 );
      }

      if( event.getName() == "FrictionCoefficientSwap" )
      {
        GEOS_LOG_RANK_0("Performing friction coefficient swap");
        FrictionCoefficientSwapMPMEvent & frictionCoefficientSwap = dynamicCast< FrictionCoefficientSwapMPMEvent & >( event );

        m_frictionCoefficient = frictionCoefficientSwap.getFrictionCoefficient();
        m_frictionCoefficientTable = frictionCoefficientSwap.getFrictionCoefficientTable();

        initializeFrictionCoefficients();
        
        event.setIsComplete( 1 );
      }
    }
  } );

}

void SolidMechanicsMPM::performMaterialSwap( ParticleManager & particleManager,
                                             string sourceRegionName,
                                             string destinationRegionName )
{
  // Material swap is performed by copying particle data from source region and constitutive model to destination region and constitutive model
  // The destination constiutive model is initialized as an empty particle subRegion
  // Presently only material swaps for polymers are implemented
  // TODO: come up with a general architecture for performing material swaps between other material models

  ParticleRegion & sourceParticleRegion = particleManager.getRegion< ParticleRegion >( sourceRegionName );
  ParticleRegion & destinationParticleRegion = particleManager.getRegion< ParticleRegion > ( destinationRegionName );

  // Copy particle data from source sub region to destination sub region
  auto & sourceSubRegions = sourceParticleRegion.getSubRegions();
  auto & destinationSubRegions = destinationParticleRegion.getSubRegions();
  for( int r=0; r < sourceSubRegions.size(); ++r)
  {
    ParticleSubRegion & sourceSubRegion = dynamicCast< ParticleSubRegion & >( *sourceSubRegions[r] );
    ParticleSubRegion & destinationSubRegion = dynamicCast< ParticleSubRegion & >( *destinationSubRegions[r] );

    // Copy particle data between sub regions
    destinationSubRegion.copyFromParticleSubRegion( sourceSubRegion );

    // Copy particle fields between sub regions
    arrayView1d< real64 > const sourceParticleMass = sourceSubRegion.getField< fields::mpm::particleMass >();
    arrayView1d< real64 > const sourceParticleInitialVolume = sourceSubRegion.getField< fields::mpm::particleInitialVolume >();
    arrayView1d< real64 > const sourceParticleDensity = sourceSubRegion.getField< fields::mpm::particleDensity >();
    arrayView1d< real64 > sourceParticleOverlap = sourceSubRegion.getField< fields::mpm::particleOverlap >();
    arrayView1d< int > sourceParticleSurfaceFlag = sourceSubRegion.getField< fields::mpm::particleSurfaceFlag >();
    arrayView1d< int > sourceIsBad = sourceSubRegion.getField< fields::mpm::isBad >();

    arrayView2d< real64 > const sourceParticleReferencePosition = sourceSubRegion.getField< fields::mpm::particleReferencePosition >();
    arrayView2d< real64 > const sourceParticleInitialMaterialDirection = sourceSubRegion.getParticleInitialMaterialDirection();
    arrayView2d< real64 > const sourceParticleMaterialDirection = sourceSubRegion.getParticleMaterialDirection();
    arrayView2d< real64 > const sourceParticleStress = sourceSubRegion.getField< fields::mpm::particleStress >();
    arrayView2d< real64 > const sourceParticlePlasticStrain = sourceSubRegion.getField< fields::mpm::particlePlasticStrain >(); 
    
    arrayView3d< real64 > const sourceParticleDeformationGradient = sourceSubRegion.getField< fields::mpm::particleDeformationGradient >();
    arrayView3d< real64 > const sourceParticleFDot = sourceSubRegion.getField< fields::mpm::particleFDot >();
    arrayView3d< real64 > const sourceParticleVelocityGradient = sourceSubRegion.getField< fields::mpm::particleVelocityGradient >();
    arrayView3d< real64 > const sourceParticleInitialRVectors = sourceSubRegion.getField< fields::mpm::particleInitialRVectors >();
    arrayView3d< real64 > const sourceParticleSphF = sourceSubRegion.getField< fields::mpm::particleSphF >();

    arrayView1d< real64 > destinationParticleMass = destinationSubRegion.getField< fields::mpm::particleMass >();
    arrayView1d< real64 > destinationParticleInitialVolume = destinationSubRegion.getField< fields::mpm::particleInitialVolume >();
    arrayView1d< real64 > destinationParticleDensity = destinationSubRegion.getField< fields::mpm::particleDensity >();
    arrayView1d< real64 > destinationParticleOverlap = destinationSubRegion.getField< fields::mpm::particleOverlap >();
    arrayView1d< int > destinationParticleSurfaceFlag = destinationSubRegion.getField< fields::mpm::particleSurfaceFlag >();
    arrayView1d< int > destinationIsBad = destinationSubRegion.getField< fields::mpm::isBad >();

    arrayView2d< real64 > const destinationParticleReferencePosition = destinationSubRegion.getField< fields::mpm::particleReferencePosition >();
    arrayView2d< real64 > const destinationParticleInitialMaterialDirection = destinationSubRegion.getParticleInitialMaterialDirection();
    arrayView2d< real64 > const destinationParticleMaterialDirection = destinationSubRegion.getParticleMaterialDirection();
    arrayView2d< real64 > const destinationParticleStress = destinationSubRegion.getField< fields::mpm::particleStress >();
    arrayView2d< real64 > const destinationParticlePlasticStrain = destinationSubRegion.getField< fields::mpm::particlePlasticStrain >(); 
    
    arrayView3d< real64 > const destinationParticleDeformationGradient = destinationSubRegion.getField< fields::mpm::particleDeformationGradient >();
    arrayView3d< real64 > const destinationParticleFDot = destinationSubRegion.getField< fields::mpm::particleFDot >();
    arrayView3d< real64 > const destinationParticleVelocityGradient = destinationSubRegion.getField< fields::mpm::particleVelocityGradient >();
    arrayView3d< real64 > const destinationParticleInitialRVectors = destinationSubRegion.getField< fields::mpm::particleInitialRVectors >();
    arrayView3d< real64 > const destinationParticleSphF = destinationSubRegion.getField< fields::mpm::particleSphF >();
    
    // destinationParticlePlasticStrain.resize(sourceParticlePlasticStrain.size());
    // destinationParticleDensity.resize(sourceParticleDensity.size());
    // destinationParticleMass.resize(sourceParticleMass.size());
    // destinationParticleDeformationGradient.resize(sourceParticleDeformationGradient.size());
    // destinationParticleFDot.resize(sourceParticleFDot.size());
    // destinationParticleVelocityGradient.resize(sourceParticleVelocityGradient.size());
    // destinationParticleInitialVolume.resize(sourceParticleInitialVolume.size());
    // destinationParticleInitialRVectors.resize(sourceParticleInitialRVectors.size());
    // destinationParticleSphF.resize(sourceParticleSphF.size());
    // destinationParticleReferencePosition.resize(sourceParticleReferencePosition.size());

    for(int pp = 0; pp < sourceSubRegion.activeParticleIndices().size(); pp++){
      localIndex const p = sourceSubRegion.activeParticleIndices()[pp];

      destinationParticleDensity[p] = sourceParticleDensity[p];
      destinationParticleMass[p] = sourceParticleMass[p];
      destinationParticleInitialVolume[p] = sourceParticleInitialVolume[p];
      destinationParticleOverlap[p] = sourceParticleOverlap[p];
      destinationParticleSurfaceFlag[p] = sourceParticleSurfaceFlag[p];
      destinationIsBad[p] = sourceIsBad[p];

      LvArray::tensorOps::copy< 3 >( destinationParticleReferencePosition[p], sourceParticleReferencePosition[p]);
      LvArray::tensorOps::copy< 3 >( destinationParticleInitialMaterialDirection[p], sourceParticleInitialMaterialDirection[p] );
      LvArray::tensorOps::copy< 3 >( destinationParticleMaterialDirection[p], sourceParticleMaterialDirection[p] );
      LvArray::tensorOps::copy< 6 >( destinationParticleStress[p], sourceParticleStress[p]);
      LvArray::tensorOps::copy< 6 >( destinationParticlePlasticStrain[p], sourceParticlePlasticStrain[p]);

      LvArray::tensorOps::copy< 3, 3 >( destinationParticleDeformationGradient[p], sourceParticleDeformationGradient[p]);
      LvArray::tensorOps::copy< 3, 3 >( destinationParticleFDot[p], sourceParticleFDot[p]);
      LvArray::tensorOps::copy< 3, 3 >( destinationParticleVelocityGradient[p], sourceParticleVelocityGradient[p]);
      LvArray::tensorOps::copy< 3, 3 >( destinationParticleInitialRVectors[p], sourceParticleInitialRVectors[p] );
      LvArray::tensorOps::copy< 3, 3 >( destinationParticleSphF[p], sourceParticleSphF[p]);
    }

    // array2d< real64 > sourceParticlePlasticStrain = sourceSubRegion.getField< fields::mpm::particlePlasticStrain >(); 
    // array1d< real64 > sourceParticleDensity = sourceSubRegion.getField< fields::mpm::particleDensity >();
    // array1d< real64 > sourceParticleMass = sourceSubRegion.getField< fields::mpm::particleMass >();
    // array3d< real64 > sourceParticleDeformationGradient = sourceSubRegion.getField< fields::mpm::particleDeformationGradient >();
    // array3d< real64 > sourceParticleFDot = sourceSubRegion.getField< fields::mpm::particleFDot >();
    // array3d< real64 > sourceParticleVelocityGradient = sourceSubRegion.getField< fields::mpm::particleVelocityGradient >();
    // array1d< real64 > sourceParticleInitialVolume = sourceSubRegion.getField< fields::mpm::particleInitialVolume >();
    // array3d< real64 > sourceParticleInitialRVectors = sourceSubRegion.getField< fields::mpm::particleInitialRVectors >();
    // array3d< real64 > sourceParticleSphF = sourceSubRegion.getField< fields::mpm::particleSphF >();
    // array2d< real64 > sourceParticleReferencePosition = sourceSubRegion.getField< fields::mpm::particleReferencePosition >();

    // destinationParticlePlasticStrain = sourceParticlePlasticStrain; 
    // destinationParticleDensity = sourceParticleDensity;
    // destinationParticleMass = sourceParticleMass;
    // destinationParticleDeformationGradient = sourceParticleDeformationGradient;
    // destinationParticleFDot = sourceParticleFDot;
    // destinationParticleVelocityGradient = sourceParticleVelocityGradient;
    // destinationParticleInitialVolume = sourceParticleInitialVolume;
    // destinationParticleInitialRVectors = sourceParticleInitialRVectors;
    // destinationParticleSphF = sourceParticleSphF;
    // destinationParticleReferencePosition = sourceParticleReferencePosition;

    // Copy constitutive data between models
    string const & sourceSolidMaterialName = sourceSubRegion.template getReference< string >( viewKeyStruct::solidMaterialNamesString() );
    string const & destinationSolidMaterialName = destinationSubRegion.template getReference< string >( viewKeyStruct::solidMaterialNamesString() );
    SolidBase & sourceSolidModel = getConstitutiveModel< SolidBase >( sourceSubRegion, sourceSolidMaterialName );
    SolidBase & destinationSolidModel = getConstitutiveModel< SolidBase >( destinationSubRegion, destinationSolidMaterialName );

    SortedArrayView< localIndex const > const activeParticleIndices = sourceSubRegion.activeParticleIndices();

    if( sourceSolidModel.hasWrapper( "plasticStrain" ) && destinationSolidModel.hasWrapper( "plasticStrain" ) )
    {
      arrayView3d< real64 > const sourcePlasticStrain = sourceSolidModel.getReference< array3d< real64 > >( "plasticStrain" );
      arrayView3d< real64 > const destinationPlasticStrain = destinationSolidModel.getReference< array3d< real64 > >( "plasticStrain" );
      forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
      {
        localIndex const p = activeParticleIndices[pp];
        LvArray::tensorOps::copy< 6 >( destinationPlasticStrain[p][0], sourcePlasticStrain[p][0]);
      } );
    }

    if( sourceSolidModel.hasWrapper( "yieldStrength" ) && destinationSolidModel.hasWrapper( "yieldStrength" ) )
    {
      arrayView2d< real64 > const sourceYieldStrength = sourceSolidModel.getReference< array2d< real64 > >( "yieldStrength" );
      arrayView2d< real64 > const destinationYieldStrength = destinationSolidModel.getReference< array2d< real64 > >( "yieldStrength" );
      forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
      {
        localIndex const p = activeParticleIndices[pp];
        destinationYieldStrength[p][0] = sourceYieldStrength[p][0];
      } );
    }

    if( sourceSolidModel.hasWrapper( "damage" ) && destinationSolidModel.hasWrapper( "damage" ) )
    {
      arrayView2d< real64 > const sourcePlasticStrain = sourceSolidModel.getReference< array2d< real64 > >( "damage" );
      arrayView2d< real64 > const destinationPlasticStrain = destinationSolidModel.getReference< array2d< real64 > >( "damage" );
      forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
      {
        localIndex const p = activeParticleIndices[pp];
        destinationPlasticStrain[p][0] = sourcePlasticStrain[p][0];
      } );
    }

    if( sourceSolidModel.hasWrapper( "jacobian" ) && destinationSolidModel.hasWrapper( "jacobian" ) )
    {
      arrayView2d< real64 > const sourceJacobian = sourceSolidModel.getReference< array2d< real64 > >( "jacobian" );
      arrayView2d< real64 > const destinationJacobian = destinationSolidModel.getReference< array2d< real64 > >( "jacobian" );
      forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
      {
        localIndex const p = activeParticleIndices[pp];
        destinationJacobian[p][0] = sourceJacobian[p][0];
      } );
    }

    if( sourceSolidModel.hasWrapper( "newStress" ) && destinationSolidModel.hasWrapper( "newStress" ) )
    {
      arrayView3d< real64 > const sourceNewStress = sourceSolidModel.getReference< array3d< real64 > >( "newStress" );
      arrayView3d< real64 > const destinationNewStress = destinationSolidModel.getReference< array3d< real64 > >( "newStress" );
      forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
      {
        localIndex const p = activeParticleIndices[pp];
        LvArray::tensorOps::copy< 6 >( destinationNewStress[p][0], sourceNewStress[p][0]);
      } );
    }

    if( sourceSolidModel.hasWrapper( "oldStress" ) && destinationSolidModel.hasWrapper( "oldStress" ) )
    {
      arrayView3d< real64 > const sourceOldStress = sourceSolidModel.getReference< array3d< real64 > >( "oldStress" );
      arrayView3d< real64 > const destinationOldStress = destinationSolidModel.getReference< array3d< real64 > >( "oldStress" );
      forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
      {
        localIndex const p = activeParticleIndices[pp];
        LvArray::tensorOps::copy< 6 >( destinationOldStress[p][0], sourceOldStress[p][0]);
      } );
    }

    // Remove particles from subregion since they now reside in destination subregion
    // Need to make a set to pass to erase ( maybe subregion needs a clear that deletes all particles or the like)
    std::set< localIndex > indicesToErase;
    for(int p = 0; p < activeParticleIndices.size(); ++p)
    {
      indicesToErase.insert(activeParticleIndices[p]);
    }
    sourceSubRegion.erase( indicesToErase ); 
    sourceSubRegion.resize( 0 );

    sourceSubRegion.setActiveParticleIndices();
    destinationSubRegion.setActiveParticleIndices();
  }
  // CC: may not need this now
  // // Need to resize regions too
  // destinationParticleRegion.resize(sourceParticleRegion.size());
  // sourceParticleRegion.resize(0);
  
  // // CC: debugging
  // GEOS_LOG_RANK_0("Source region has " << sourceSubRegions.size() << " subregions and " << sourceParticleRegion.getNumberOfParticles() << " particles.");

  // sourceParticleRegion.forParticleSubRegionsIndex( [&]( localIndex const subRegionIndex, 
  //                                                       ParticleSubRegion & subRegion )
  // {
  //   GEOS_LOG_RANK_0("Source subregion index: " << subRegionIndex << ", Num particles: " << subRegion.activeParticleIndices().size() << " from region " << sourceRegionName);
  // });

  // GEOS_LOG_RANK_0("Destination region has " << destinationSubRegions.size() << " subregions  and " << sourceParticleRegion.getNumberOfParticles() << " particles.");
  // destinationParticleRegion.forParticleSubRegionsIndex( [&]( localIndex const subRegionIndex, 
  //                                                            ParticleSubRegion & subRegion )
  // {
  //   GEOS_LOG_RANK_0("Destination subregion index: " << subRegionIndex << ", Num particles: " << subRegion.activeParticleIndices().size() << " from region " << destinationRegionName);
  // });

  // particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  // {
  //   GEOS_LOG_RANK_0(subRegion.size() << ", " << subRegion.activeParticleIndices().size() << ", " << subRegion.getName() );
  // });
}

void SolidMechanicsMPM::syncGridFields( std::vector< std::string > const & fieldNames,
                                        DomainPartition & domain,
                                        NodeManager & nodeManager,
                                        MeshLevel & mesh,
                                        MPI_Op op )
{
  // (0) Bring grid fields to host
  for( auto const & name : fieldNames )
  {
    WrapperBase & wrapper = nodeManager.getWrapperBase( name );
    wrapper.move( LvArray::MemorySpace::host, true );
  }

  // (1) Initialize
  FieldIdentifiers fieldsToBeSynced;
  fieldsToBeSynced.addFields( FieldLocation::Node, fieldNames );
  std::vector< NeighborCommunicator > & neighbors = domain.getNeighbors();

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
  CommunicationTools::getInstance().finalizeUnpack( mesh, neighbors, m_iComm, true, unpackEvents, op ); // needs an extra argument to
                                                                                                        // indicate unpack operation

  // (4) Swap send and receive indices back so we can sync from master to ghost
  for( size_t n=0; n<neighbors.size(); n++ )
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
                                                         arrayView3d< real64 > const & vectorMultiField,
                                                         arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const gridPosition,
                                                         Group & nodeSets )
{
  // This is a helper function for enforcing symmetry BCs on a single face and is meant to be called by other functions, not directly by the
  // solver:
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
    SortedArrayView< localIndex const > const boundaryNodes = m_boundaryNodes[face].toView();
    int const numBoundaryNodes = boundaryNodes.size();
    forAll< serialPolicy >( numBoundaryNodes, [&, vectorMultiField] GEOS_HOST ( localIndex const gg ) // Probably not a big enough loop to
                                                                                                      // warrant parallelization
      {
        int const g = boundaryNodes[gg];
        vectorMultiField[g][fieldIndex][dir0] = 0.0;
      } );

    // Perform field reflection on buffer nodes
    SortedArrayView< localIndex const > const bufferNodes = m_bufferNodes[face].toView();
    int const numBufferNodes = bufferNodes.size();
    forAll< serialPolicy >( numBufferNodes, [&, vectorMultiField, gridPosition] GEOS_HOST ( localIndex const gg ) // Probably not a big
                                                                                                                  // enough loop to warrant
                                                                                                                  // parallelization
      {
        int const g = bufferNodes[gg];
        int ijk[3];
        ijk[dir0] = positiveNormal * ( m_nEl[dir0] - 2 ) + ( 1 - positiveNormal ) * ( 2 );
        ijk[dir1] = std::round( ( gridPosition[g][dir1] - m_xLocalMin[dir1] ) / m_hEl[dir1] );
        ijk[dir2] = std::round( ( gridPosition[g][dir2] - m_xLocalMin[dir2] ) / m_hEl[dir2] );

        localIndex gFrom = m_ijkMap[ijk[0]][ijk[1]][ijk[2]];

        vectorMultiField[g][fieldIndex][dir0] = -vectorMultiField[gFrom][fieldIndex][dir0]; // Negate component aligned with surface normal
        vectorMultiField[g][fieldIndex][dir1] =  vectorMultiField[gFrom][fieldIndex][dir1];
        vectorMultiField[g][fieldIndex][dir2] =  vectorMultiField[gFrom][fieldIndex][dir2];
      } );
  }
}

void SolidMechanicsMPM::enforceGridVectorFieldSymmetryBC( arrayView3d< real64 > const & vectorMultiField,
                                                          arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const gridPosition,
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
                                           NodeManager & nodeManager )
{
  // Get grid fields
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const gridPosition = nodeManager.referencePosition();
  arrayView1d< int const > const gridGhostRank = nodeManager.ghostRank();
  arrayView2d< real64 > const gridMass = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::massString() );
  arrayView3d< real64 > const gridVelocity = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::velocityString() );
  arrayView3d< real64 > const gridAcceleration = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::accelerationString() );

  // Get node sets
  Group & nodeSets = nodeManager.sets();
  array1d< SortedArray< localIndex > > & m_boundaryNodes = nodeSets.getReference< array1d< SortedArray< localIndex > > >( viewKeyStruct::boundaryNodesString() );
  array1d< SortedArray< localIndex > > & m_bufferNodes = nodeSets.getReference< array1d< SortedArray< localIndex > > >( viewKeyStruct::bufferNodesString() );

  // Impose BCs on each face while gathering reaction forces
  real64 localFaceReactions[6] = {0.0};
  for( int face = 0; face < 6; face++ )
  {
    if( m_boundaryConditionTypes[face] == 1 )
    {
      singleFaceVectorFieldSymmetryBC( face, gridVelocity, gridPosition, nodeSets );
      singleFaceVectorFieldSymmetryBC( face, gridAcceleration, gridPosition, nodeSets );
    }
    else if( m_boundaryConditionTypes[face] == 2 && ( m_prescribedBoundaryFTable == 1 || m_stressControl[0] == 1 || m_stressControl[1] == 1 || m_stressControl[2] == 1) )
    {
      for( int fieldIndex = 0; fieldIndex < m_numVelocityFields; fieldIndex++ )
      {
        // Face-associated quantities
        int dir0 = face / 2;           // 0, 0, 1, 1, 2, 2 (x-, x+, y-, y+, z-, z+)
        int dir1 = (dir0 + 1) % 3;     // 1, 1, 2, 2, 0, 0
        int dir2 = (dir0 + 2) % 3;     // 2, 2, 0, 0, 1, 1
        int positiveNormal = face % 2; // even => (-) => 0, odd => (+) => 1

        // Enforce BCs on boundary nodes using F-table
        SortedArrayView< localIndex const > const boundaryNodes = m_boundaryNodes[face].toView();
        int const numBoundaryNodes = boundaryNodes.size();
        forAll< serialPolicy >( numBoundaryNodes, [&, gridPosition, gridVelocity, gridMass] GEOS_HOST ( localIndex const gg ) // Probably
                                                                                                                              // not a big
                                                                                                                              // enough loop
                                                                                                                              // to warrant
                                                                                                                              // parallelization
          {
            int const g = boundaryNodes[gg];
            real64 prescribedVelocity = m_domainL[dir0] * gridPosition[g][dir0];
            real64 accelerationForBC = (prescribedVelocity - gridVelocity[g][fieldIndex][dir0]) / dt; // acceleration needed to satisfy BC
            gridVelocity[g][fieldIndex][dir0] = prescribedVelocity;
            gridAcceleration[g][fieldIndex][dir0] += accelerationForBC;
            if( gridGhostRank[g] <= -1 ) // so we don't double count reactions at partition boundaries
            {
              localFaceReactions[face] += accelerationForBC * gridMass[g][fieldIndex];
            }
          } );

        // Perform field reflection on buffer nodes - accounts for moving boundary effects
        SortedArrayView< localIndex const > const bufferNodes = m_bufferNodes[face].toView();
        int const numBufferNodes = bufferNodes.size();
        // Possibly not a big enough loop to warrant parallelization
        forAll< serialPolicy >( numBufferNodes, [&, gridPosition, gridVelocity, gridAcceleration] GEOS_HOST ( localIndex const gg )
          {
            int const g = bufferNodes[gg];

            // Initialize grid ijk indices
            int ijk[3];
            ijk[dir1] = std::round((gridPosition[g][dir1] - m_xLocalMin[dir1]) / m_hEl[dir1] );
            ijk[dir2] = std::round((gridPosition[g][dir2] - m_xLocalMin[dir2]) / m_hEl[dir2] );

            // Grab the node index that we're copying from
            ijk[dir0] = positiveNormal * (m_nEl[dir0] - 2) + (1 - positiveNormal) * (2);
            localIndex gFrom = m_ijkMap[ijk[0]][ijk[1]][ijk[2]];

            // Grab the associated boundary node index for moving boundary correction
            ijk[dir0] = positiveNormal * (m_nEl[dir0] - 1) + (1 - positiveNormal) * (1);
            localIndex gBoundary = m_ijkMap[ijk[0]][ijk[1]][ijk[2]];

            // Calculate velocity, Negate component aligned with surface normal and correct for moving boundary
            gridVelocity[g][fieldIndex][dir0] = -gridVelocity[gFrom][fieldIndex][dir0] + 2.0 * gridVelocity[gBoundary][fieldIndex][dir0];
            gridVelocity[g][fieldIndex][dir1] = gridVelocity[gFrom][fieldIndex][dir1];
            gridVelocity[g][fieldIndex][dir2] = gridVelocity[gFrom][fieldIndex][dir2];

            // Calculate acceleration, Negate component aligned with surface normal and correct for moving boundary
            gridAcceleration[g][fieldIndex][dir0] = -gridAcceleration[gFrom][fieldIndex][dir0] + 2.0 * gridAcceleration[gBoundary][fieldIndex][dir0];
            gridAcceleration[g][fieldIndex][dir1] = gridAcceleration[gFrom][fieldIndex][dir1];
            gridAcceleration[g][fieldIndex][dir2] = gridAcceleration[gFrom][fieldIndex][dir2];
          } );
      }
    }
  }

  // Reduce reaction forces from all partitions
  real64 globalFaceReactions[6];
  for( int face = 0; face < 6; face++ )
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
  if( MpiWrapper::commRank( MPI_COMM_GEOSX ) == 0 && m_reactionHistory == 1 )
  {
    if(time_n + dt >= m_nextReactionWriteTime )
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
      file << std::setprecision( std::numeric_limits< long double >::digits10 )
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
      m_nextReactionWriteTime += m_reactionWriteInterval;
    }
  }
}

void SolidMechanicsMPM::computeGridSurfaceNormals( ParticleManager & particleManager,
                                                   NodeManager & nodeManager )
{
  // Grid fields
  arrayView3d< real64 > const gridSurfaceNormal = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::surfaceNormalString() );
  arrayView2d< real64 const > const gridDamageGradient = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::damageGradientString() );

  localIndex subRegionIndex = 0;
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    // Particle fields
    arrayView1d< real64 const > const particleVolume = subRegion.getParticleVolume();
    arrayView1d< int const > const particleGroup = subRegion.getParticleGroup();
    arrayView2d< real64 const > const particleDamageGradient = subRegion.getField< fields::mpm::particleDamageGradient >();

    // Get views to mapping arrays
    int const numberOfVerticesPerParticle = subRegion.numberOfVerticesPerParticle();
    arrayView2d< localIndex const > const mappedNodes = m_mappedNodes[subRegionIndex];
    arrayView3d< real64 const > const shapeFunctionGradientValues = m_shapeFunctionGradientValues[subRegionIndex];

    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
    int const numDims = m_numDims;
    int const damageFieldPartitioning = m_damageFieldPartitioning;
    int const numContactGroups = m_numContactGroups;
    forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST ( localIndex const pp ) // Can parallelize with atomics
      {
        localIndex const p = activeParticleIndices[pp];

        // Map to grid
        for( int g = 0; g < 8 * numberOfVerticesPerParticle; g++ )
        {
          localIndex const mappedNode = mappedNodes[pp][g];
          // 0 undamaged or "A" field, 1 for "B" field
          int const nodeFlag = ( damageFieldPartitioning == 1 && LvArray::tensorOps::AiBi< 3 >( gridDamageGradient[mappedNode], particleDamageGradient[p] ) < 0.0 ) ? 1 : 0;
          int const fieldIndex = nodeFlag * numContactGroups + particleGroup[p]; // This ranges from 0 to nMatFields-1
          for( int i=0; i<numDims; i++ )
          {
            gridSurfaceNormal[mappedNode][fieldIndex][i] += shapeFunctionGradientValues[pp][g][i] * particleVolume[p];
          }
        }
      } ); // particle loop

    // Increment subregion index
    subRegionIndex++;
  } ); // subregion loop
}

void SolidMechanicsMPM::normalizeGridSurfaceNormals( arrayView2d< real64 const > const & gridMass,
                                                     arrayView3d< real64 > const & gridSurfaceNormal )
{
  int const numNodes = gridSurfaceNormal.size( 0 );
  int const numVelocityFields = m_numVelocityFields;
  real64 const smallMass = m_smallMass;
  int const planeStrain = m_planeStrain;
  forAll< serialPolicy >( numNodes, [=] GEOS_HOST_DEVICE ( localIndex const g )
  {
    for( localIndex fieldIndex = 0; fieldIndex < numVelocityFields; fieldIndex++ )
    {
      if( gridMass[g][fieldIndex] > smallMass ) // small mass threshold
      {
        arraySlice1d< real64 > const surfaceNormal = gridSurfaceNormal[g][fieldIndex];
        real64 norm = planeStrain == 1 ? sqrt( surfaceNormal[0] * surfaceNormal[0] + surfaceNormal[1] * surfaceNormal[1] ) : LvArray::tensorOps::l2Norm< 3 >( surfaceNormal );
        if( norm > 0.0 ) // TODO: Set a finite threshold?
        {
          LvArray::tensorOps::scale< 3 >( surfaceNormal, 1.0/norm );
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
  } );
}

void SolidMechanicsMPM::initializeFrictionCoefficients()
{
 // Table is only used if the frictionCoefficient in the input file is not specified
  if( static_cast<int>( m_frictionCoefficient ) == -1 )
  {
    GEOS_ERROR_IF( m_frictionCoefficientTable.size(0) == 0, "No frictionCoefficient or frictionCoefficientTable was defined." );

    GEOS_ERROR_IF( m_frictionCoefficientTable.size(0) != m_frictionCoefficientTable.size(1), "frictionCoefficientTable must be square.");
    GEOS_ERROR_IF( m_frictionCoefficientTable.size(0) != m_numContactGroups, "frictionCoefficientTable must have the same number of rows and columns as the number of contact groups.");

    for(int i = 0; i < m_numContactGroups; i++)
    {
      for(int j = i+1; i < m_numContactGroups; j++)
      {
        GEOS_ERROR_IF( std::abs(m_frictionCoefficientTable[i][j] - m_frictionCoefficientTable[j][i]) > DBL_EPSILON, "Off-diagonal friction coefficients must match" );
      }
    }

  } 
  else
  {
    GEOS_ERROR_IF( m_frictionCoefficient < 0.0, "Friction coefficient must be positive.");
    m_frictionCoefficientTable.resize( m_numContactGroups, m_numContactGroups );
    for(int i = 0; i < m_numContactGroups; i++)
    {
      for(int j = 0; j < m_numContactGroups; j++)
      {
        m_frictionCoefficientTable[i][j] = m_frictionCoefficient;
      }
    }
  }
}

void SolidMechanicsMPM::computeContactForces( real64 const dt,
                                              arrayView2d< real64 const > const & gridMass,
                                              arrayView2d< real64 const > const & gridDamage,
                                              arrayView2d< real64 const > const & gridMaxDamage,
                                              arrayView3d< real64 const > const & gridVelocity,
                                              arrayView3d< real64 const > const & gridMomentum,
                                              arrayView3d< real64 const > const & gridSurfaceNormal,
                                              arrayView3d< real64 const > const & gridMaterialPosition,
                                              arrayView3d< real64 > const & gridContactForce )
{
  // Get number of nodes
  int numNodes = gridMass.size( 0 );

  forAll< serialPolicy >( numNodes, [&, gridMass, gridVelocity, gridMomentum, gridSurfaceNormal, gridMaterialPosition, gridContactForce] GEOS_HOST ( localIndex const g )
    {
      // Initialize gridContactForce[g] to zero. TODO: This shouldn't be necessary?
      for( int fieldIndex = 0; fieldIndex < m_numVelocityFields; fieldIndex++ )
      {
        for( int i=0; i<3; i++ )
        {
          gridContactForce[g][fieldIndex][i] = 0.0;
        }
      }

      // Loop over all possible field pairings and enforce contact on each pair.
      // gridContactForce will be gradually updated due to a '+=' in computePairwiseNodalContactForce
      for( localIndex A = 0; A < m_numVelocityFields - 1; A++ )
      {
        for( localIndex B = A + 1; B < m_numVelocityFields; B++ )
        {
          // Make sure both fields in the pair are active
          bool active = ( gridMass[g][A] > m_smallMass ) && ( LvArray::tensorOps::l2NormSquared< 3 >( gridSurfaceNormal[g][A] ) > 1.0e-16 )
                        and
                        ( gridMass[g][B] > m_smallMass ) && ( LvArray::tensorOps::l2NormSquared< 3 >( gridSurfaceNormal[g][B] ) > 1.0e-16 );

          if( active )
          {
            // Evaluate the separability criterion for the contact pair.
            int separable = evaluateSeparabilityCriterion( A,
                                                           B,
                                                           gridDamage[g][A],
                                                           gridDamage[g][B],
                                                           gridMaxDamage[g][A],
                                                           gridMaxDamage[g][B] );

            real64 frictionCoefficient = m_frictionCoefficientTable[A % m_numContactGroups][B % m_numContactGroups];
            computePairwiseNodalContactForce( separable,
                                              dt,
                                              frictionCoefficient,
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
    } );
}

void SolidMechanicsMPM::computePairwiseNodalContactForce( int const & separable,
                                                          real64 const & dt,
                                                          real64 const & frictionCoefficient,
                                                          real64 const & mA,
                                                          real64 const & mB,
                                                          arraySlice1d< real64 const > const vA,
                                                          arraySlice1d< real64 const > const GEOS_UNUSED_PARAM( vB ),
                                                          arraySlice1d< real64 const > const qA,
                                                          arraySlice1d< real64 const > const qB,
                                                          arraySlice1d< real64 const > const (nA),
                                                          arraySlice1d< real64 const > const (nB),
                                                          arraySlice1d< real64 const > const xA, // Position of field A
                                                          arraySlice1d< real64 const > const xB, // Position of field B
                                                          arraySlice1d< real64 > const fA,
                                                          arraySlice1d< real64 > const fB )
{
  // Total mass for the contact pair.
  real64 mAB = mA + mB;

  // Outward normal of field A with respect to field B.
  real64 nAB[3];

  // Use the surface normal for whichever field has more mass.
  // if( mA > mB )
  // {
  //   nAB[0] = nA[0];
  //   nAB[1] = nA[1];
  //   nAB[2] = nA[2];
  // }
  // else
  // {
  //   nAB[0] = -nB[0];
  //   nAB[1] = -nB[1];
  //   nAB[2] = -nB[2];
  // }

  // CC: debug
  // GEOS_LOG_RANK( "Mass-weighted average of the field normals" );

  // Mass-weighted average of the field normals
  for( int i=0; i<3; i++ )
  {
    nAB[i] = nA[i] * mA - nB[i] * mB;
  }

  // Vector pointing from CoM of A field to CoM of B field, stretched by grid spacing
  // for( int i=0; i<3; i++ )
  // {
  //   nAB[i] = (xB[i] - xA[i]) / m_hEl[i];
  // }

  // CC: debug
  // GEOS_LOG_RANK( "Normalize the effective surface normal" );

  // Normalize the effective surface normal
  if( m_planeStrain == 1 )
  {
    nAB[2] = 0.0;
  }
  real64 norm = sqrt( nAB[0] * nAB[0] + nAB[1] * nAB[1] + nAB[2] * nAB[2] );
  nAB[0] /= norm;
  nAB[1] /= norm;
  nAB[2] /= norm;

  // CC: debug
  // GEOS_LOG_RANK("Calculate the contact gap between the fields");

  // Calculate the contact gap between the fields
  real64 gap0;
  if( m_planeStrain == 1 )
  {
      // CC: debug
    // GEOS_LOG_RANK("Plane strain");
    gap0 = ( m_hEl[0]*m_hEl[1] ) / sqrt( m_hEl[1]*m_hEl[1]*nAB[0]*nAB[0] + m_hEl[0]*m_hEl[0]*nAB[1]*nAB[1] );
  }
  else
  {
    // CC: debug
    // GEOS_LOG_RANK("Not plane strain");
    // GEOS_LOG_RANK("m_hEl:" << m_hEl[0] << ", " << m_hEl[1] << ", " << m_hEl[2] );
    // GEOS_LOG_RANK("nAB: " << nAB[0] << ", " << nAB[1] << ", " << nAB[2]);
    // GEOS_LOG_RANK("Numerator: " << (m_hEl[0]*m_hEl[1]*m_hEl[2]));
    // GEOS_LOG_RANK("Determinant: " << m_hEl[2]*m_hEl[2]*nAB[2]*nAB[2]*( m_hEl[1]*m_hEl[1]*nAB[0]*nAB[0] + m_hEl[0]*m_hEl[0]*nAB[1]*nAB[1] ) + m_hEl[0]*m_hEl[0]*m_hEl[1]*m_hEl[1]*( nAB[0]*nAB[0] + nAB[1]*nAB[1] ) );
    // GEOS_LOG_RANK("Denominator: " << sqrt( m_hEl[2]*m_hEl[2]*nAB[2]*nAB[2]*( m_hEl[1]*m_hEl[1]*nAB[0]*nAB[0] + m_hEl[0]*m_hEl[0]*nAB[1]*nAB[1] ) + m_hEl[0]*m_hEl[0]*m_hEl[1]*m_hEl[1]*( nAB[0]*nAB[0] + nAB[1]*nAB[1] ) ));
    
    // CC: debug
    // currently failes with polymer model (suspect compliant materials) when z direction boundary type is 2
    //This needs to be corrected
    // gap0 = (m_hEl[0]*m_hEl[1]*m_hEl[2]) /
    //        sqrt( m_hEl[2]*m_hEl[2]*nAB[2]*nAB[2]*( m_hEl[1]*m_hEl[1]*nAB[0]*nAB[0] + m_hEl[0]*m_hEl[0]*nAB[1]*nAB[1] ) + m_hEl[0]*m_hEl[0]*m_hEl[1]*m_hEl[1]*( nAB[0]*nAB[0] + nAB[1]*nAB[1] ) );
 
    // Ellipse solution
    // Gives element sizes in each direction and reasonable approximations in others
    gap0 = 1/sqrt(std::pow(nAB[0]/m_hEl[0],2) + std::pow(nAB[1]/m_hEl[1],2) + std::pow(nAB[2]/m_hEl[2],2));

  }

  // CC: debug
  // GEOS_LOG_RANK("Fudge factor");

  // TODO: A fudge factor of 0.67 on gap0 makes diagonal surfaces (wrt grid) close better I think, but this is more general
  real64 gap = (xB[0] - xA[0]) * nAB[0] + (xB[1] - xA[1]) * nAB[1] + (xB[2] - xA[2]) * nAB[2] - gap0;

  // CC: debug
  // GEOS_LOG_RANK("Total momentum for the contact pair.");

  // Total momentum for the contact pair.
  real64 qAB[3];
  qAB[0] = qA[0] + qB[0];
  qAB[1] = qA[1] + qB[1];
  qAB[2] = qA[2] + qB[2];

  // CC: debug
  // GEOS_LOG_RANK(" Center-of-mass velocity for the contact pair.");

  // Center-of-mass velocity for the contact pair.
  real64 vAB[3];
  vAB[0] = qAB[0] / mAB;
  vAB[1] = qAB[1] / mAB;
  vAB[2] = qAB[2] / mAB;

  // Compute s1AB and s2AB, to form an orthonormal basis. This uses the method by E. Herbold
  // to ensure consistency between surfaces.
  real64 s1AB[3], s2AB[3]; // Tangential vectors for the contact pair
  computeOrthonormalBasis( nAB, s1AB, s2AB );

  // CC: debug
  // GEOS_LOG_RANK("Compute force decomposition, declare increment in fA from 'this' contact pair");

  // Compute force decomposition, declare increment in fA from "this" contact pair
  real64 fnor =  ( mA / dt ) * ( (vAB[0] - vA[0]) * nAB[0]  + (vAB[1] - vA[1]) * nAB[1]  + (vAB[2] - vA[2]) * nAB[2] ),
         ftan1 = ( mA / dt ) * ( (vAB[0] - vA[0]) * s1AB[0] + (vAB[1] - vA[1]) * s1AB[1] + (vAB[2] - vA[2]) * s1AB[2] ),
         ftan2 = ( mA / dt ) * ( (vAB[0] - vA[0]) * s2AB[0] + (vAB[1] - vA[1]) * s2AB[1] + (vAB[2] - vA[2]) * s2AB[2] );
  real64 dfA[3];

  // CC: debug
  // GEOS_LOG_RANK("Check for separability, and enforce either slip, or no-slip contact, accordingly");

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


    real64 ftan = std::min( frictionCoefficient * std::abs( fnor ), ftanMag ); // This goes to zero when contact=0 due to the std::min
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

void SolidMechanicsMPM::computeOrthonormalBasis( const real64 * e1, // input "normal" unit vector.
                                                 real64 * e2,       // output "tangential" unit vector.
                                                 real64 * e3 )      // output "tangential" unit vector.
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
{// TODO: Find a way to automatically loop over these fields.
 // I feel like a dimension check would work.
  // Generate labels
  std::vector< std::string > labels1( m_numVelocityFields );
  std::generate( labels1.begin(), labels1.end(), [i=0]() mutable { return "velocityField" + std::to_string( i++ ); } );
  string const labels2[] = { "X", "Y", "Z" };

  // Apply labels to scalar multi-fields
  std::vector< std::string > keys2d = { viewKeyStruct::massString(),
                                        viewKeyStruct::damageString(),
                                        viewKeyStruct::maxDamageString()};
  for( auto const & key: keys2d )
  {
    WrapperBase & wrapper = nodeManager.getWrapper< array2d< real64 > >( key );
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
  for( auto const & key: keys3d )
  {
    WrapperBase & wrapper = nodeManager.getWrapper< array3d< real64 > >( key );
    wrapper.setDimLabels( 1, labels1 );
    wrapper.setDimLabels( 2, labels2 );
  }
}

void SolidMechanicsMPM::resizeGrid( SpatialPartition & partition,
                                    NodeManager & nodeManager,
                                    real64 const dt )
{
  // Modify SpatialPartition class members
  partition.updateSizes( m_domainL, dt );

  // Modify SolidMechanicsMPM class members
  real64 ratio[3];
  for( int i=0; i<3; i++ )
  {
    // Incremental stretch
    ratio[i] = 1.0 + m_domainL[i] * dt;

    // Modify SolidMechanicsMPM class members
    m_hEl[i] *= ratio[i];
    m_xLocalMin[i] *= ratio[i];
    m_xLocalMax[i] *= ratio[i];
    m_xLocalMinNoGhost[i] *= ratio[i];
    m_xLocalMaxNoGhost[i] *= ratio[i];
    m_xGlobalMin[i] *= ratio[i];
    m_xGlobalMax[i] *= ratio[i];
    m_partitionExtent[i] *= ratio[i];
    m_domainExtent[i] *= ratio[i];
  }

  // Update nodal positions
  arrayView2d< real64, nodes::REFERENCE_POSITION_USD > const gridPosition = nodeManager.referencePosition();
  forAll< serialPolicy >( gridPosition.size( 0 ), [=] GEOS_HOST_DEVICE ( localIndex const g )
  {
    gridPosition[g][0] *= ratio[0];
    gridPosition[g][1] *= ratio[1];
    gridPosition[g][2] *= ratio[2];
  } );
}

void SolidMechanicsMPM::solverProfiling( std::string label )
{
  if( m_solverProfiling >= 1 )
  {
    MPI_Barrier( MPI_COMM_GEOSX );
    GEOS_LOG_RANK_IF( m_solverProfiling == 2, label );
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
  GEOS_ERROR_IF( solidMaterialName.empty(), GEOS_FMT( "SolidBase model not found on subregion {}", subRegion.getName() ) );
}

void SolidMechanicsMPM::setConstitutiveNames( ParticleSubRegionBase & subRegion ) const
{
  GEOS_UNUSED_VAR( subRegion );
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

  // Declare bin key and bins data structure
  BinKey binKey;
  std::unordered_map< BinKey, std::vector< localIndex >, BinKeyHash > bins;

  // Reverse entries in bins based on even distribution of particles in partition - OPTIONAL
  particleManager.forParticleRegions< ParticleRegion >( [&]( ParticleRegion & region ) // idk why this requires a template argument and the
                                                                                       // subregion loops don't
  {
    binKey.regionIndex = region.getIndexInParent();
    region.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
    {
      binKey.subRegionIndex = subRegion.getIndexInParent();
      int numReserved = std::ceil( subRegion.size() / nbins );

      // Loop over bins
      for( int i=0; i<nxbins; i++ )
      {
        for( int j=0; j<nybins; j++ )
        {
          for( int k=0; k<nzbins; k++ )
          {
            binKey.binIndex = i + j * nxbins + k * nxbins * nybins;
            bins[binKey].reserve( numReserved );
          }
        }
      }
    } );
  } );

  // Populate bins with local particle indices
  particleManager.forParticleRegions< ParticleRegion >( [&]( ParticleRegion & region ) // idk why this requires a template argument and the
                                                                                       // subregion loops don't
  {
    binKey.regionIndex = region.getIndexInParent();
    region.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
    {
      binKey.subRegionIndex = subRegion.getIndexInParent();
      arrayView2d< real64 const > const particlePosition = subRegion.getParticleCenter();
      forAll< serialPolicy >( subRegion.size(), [&, particlePosition] GEOS_HOST ( localIndex const p ) // host only because of bin data
                                                                                                       // structure, need atomics in
                                                                                                       // parallel since multiple particles
                                                                                                       // can write to the same bin
        {
          // Particle bin ijk indices
          int i, j, k;
          i = std::floor( ( particlePosition[p][0] - xmin ) / dx ),
          j = std::floor( ( particlePosition[p][1] - ymin ) / dy ),
          k = std::floor( ( particlePosition[p][2] - zmin ) / dz );

          // Bin number
          binKey.binIndex = i + j * nxbins + k * nxbins * nybins;

          // Add particle to bin
          bins[binKey].push_back( p );
        } );
    } );
  } );

  // Perform neighbor search over appropriate bins
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegionA )
  {
    // Get and initialize neighbor list
    OrderedVariableToManyParticleRelation & neighborList = subRegionA.neighborList();
    neighborList.resize( 0, 0 ); // Clear the existing neighbor list
    neighborList.resize( subRegionA.size() );
    int numToReserve = m_planeStrain == 1 ? 25 : 179; // Assuming square/cubic cells, 2 particles per cell in each dimension, and a neighbor
                                                      // radius equal to the cell diagonal length
    reserveNeighbors( neighborList, numToReserve );

    // Get 'this' particle's location
    arrayView2d< real64 > const xA = subRegionA.getParticleCenter();

    // Find neighbors of 'this' particle
    SortedArrayView< localIndex const > const subRegionAActiveParticleIndices = subRegionA.activeParticleIndices();
    forAll< serialPolicy >( subRegionAActiveParticleIndices.size(), [&, subRegionAActiveParticleIndices, xA] GEOS_HOST ( localIndex const pp )
      {// host only because of bin data structure and invocation of particleManager; not thread-safe for some reason
        // Particle A index
        localIndex a = subRegionAActiveParticleIndices[pp];

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
        binKey.regionIndex = region.getIndexInParent();
        binKey.subRegionIndex = subRegionB.getIndexInParent();

        // Get 'other' particle location
        arrayView2d< real64 > const xB = subRegionB.getParticleCenter();

        // Declare temporary array of local neighbor indices
        std::vector< localIndex > localIndices;
        localIndices.reserve( numToReserve );

        // Loop over bins
        int count = 0; // Count the number of neighbors on this subregion
        for( int iBin=imin; iBin<=imax; iBin++ )
        {
          for( int jBin=jmin; jBin<=jmax; jBin++ )
          {
            for( int kBin=kmin; kBin<=kmax; kBin++ )
            {
              binKey.binIndex = iBin + jBin * nxbins + kBin * nxbins * nybins;
              for( localIndex & b: bins[binKey] )
              {
                real64 xBA[3];
                xBA[0] = xB[b][0] - xA[a][0];
                xBA[1] = xB[b][1] - xA[a][1];
                xBA[2] = xB[b][2] - xA[a][2];
                real64 rSquared = xBA[0] * xBA[0] + xBA[1] * xBA[1] + xBA[2] * xBA[2];
                if( rSquared <= neighborRadiusSquared ) // Would you be my neighbor?
                {
                  count++;
                  localIndices.push_back( b );
                }
              }
            }
          }
        }

        // Populate the temporary arrays of neighbor region and subregion indices
        std::vector< localIndex > regionIndices( count, binKey.regionIndex );
        std::vector< localIndex > subRegionIndices( count, binKey.subRegionIndex );

        // Insert indices into neighbor list
        insertMany( neighborList,
                    a,
                    regionIndices,
                    subRegionIndices,
                    localIndices );
      } );
      } );
  } );

  return( MPI_Wtime() - tStart );
}

void SolidMechanicsMPM::optimizeBinSort( ParticleManager & particleManager )
{
  // Each partition determines its optimal multiplier which results in the minimum time for neighbor list construction
  // The global multiplier is set by a weighted average of each partition's multiplier, with the number of particles
  // on the partition being the weight factor.

  // Start
  int optimalMultiplier = 1;

  // Identify the largest possible multiplier - we can limit this if it's prohibitive to check larger multipliers in large 3D sims
  real64 xL = m_xLocalMaxNoGhost[0] - m_xLocalMinNoGhost[0] + 2 * m_neighborRadius;
  real64 yL = m_xLocalMaxNoGhost[1] - m_xLocalMinNoGhost[1] + 2 * m_neighborRadius;
  real64 zL = m_xLocalMaxNoGhost[2] - m_xLocalMinNoGhost[2] + 2 * m_neighborRadius;
  int maxMultiplier = std::max( std::ceil( xL / m_neighborRadius ), std::max( std::ceil( yL / m_neighborRadius ), std::ceil( zL / m_neighborRadius ) ));
  maxMultiplier = std::max( maxMultiplier, 1 );

  // Identify this partition's optimal multiplier
  real64 minTime = DBL_MAX;
  for( int multiplier=1; multiplier<=maxMultiplier; multiplier++ )
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

real64 SolidMechanicsMPM::kernel( real64 const & r ) // distance from particle to query point.
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

void SolidMechanicsMPM::kernelGradient( arraySlice1d< real64 const > const x,  // query point
                                        std::vector< real64 > & xp,            // particle location
                                        real64 const & r,                      // distance from particle to query point.
                                        real64 * result )
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

real64 SolidMechanicsMPM::computeKernelField( arraySlice1d< real64 const > const x,  // query point
                                              arrayView2d< real64 const > const xp,  // List of neighbor particle locations.
                                              arrayView1d< real64 const > const Vp,  // List of neighbor particle volumes.
                                              arrayView1d< real64 const > const fp ) // scalar field values (e.g. damage) at neighbor
                                                                                     // particles
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
  for( localIndex p = 0; p < Vp.size(); ++p )
  {
    relativePosition[0] = x[0] - xp[p][0];
    relativePosition[1] = x[1] - xp[p][1];
    relativePosition[2] = x[2] - xp[p][2];
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

void SolidMechanicsMPM::computeKernelFieldGradient( arraySlice1d< real64 const > const x,       // query point
                                                    std::vector< std::vector< real64 > > & xp,  // List of neighbor particle locations.
                                                    std::vector< real64 > & Vp,                 // List of neighbor particle volumes.
                                                    std::vector< real64 > & fp,                 // scalar field values (e.g. damage) at
                                                                                                // neighbor particles
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

  for( unsigned int p = 0; p < Vp.size(); ++p )
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

void SolidMechanicsMPM::computeKernelVectorGradient( arraySlice1d< real64 const > const x,       // query point
                                                     std::vector< std::vector< real64 > > & xp,  // List of neighbor particle locations.
                                                     std::vector< real64 > & Vp,                 // List of neighbor particle volumes.
                                                     std::vector< std::vector< real64 > > & fp,                 // vector field values (e.g.
                                                                                                                // velocity) at neighbor
                                                                                                                // particles
                                                     arraySlice2d< real64 > const result )
{
  // Compute the kernel scalar field at a point, for a given list of neighbor particles.
  // The lists xp, fp, and the length np could refer to all the particles in the patch,
  // but generally this function will be evaluated with x equal to some particle center,
  // and xp, fp, will be lists for just the neighbors of the particle.
  // TODO: Modify to also "return" the kernel field value

  // Scalar kernel field values
  real64 kernelVal,
         f[3] = { 0.0 },
         k = 0.0,
         r;

  // Gradient of the scalar field
  real64 relativePosition[3],
         fGrad[3][3] = { { 0.0 } },
         kGrad[3] = { 0.0 },
         kernelGradVal[3];

  for( unsigned int p = 0; p < Vp.size(); ++p )
  {
    for( int i = 0; i < 3; i++ )
    {
      relativePosition[i] = x[i] - xp[p][i];
    }
    r = sqrt( relativePosition[0] * relativePosition[0] + relativePosition[1] * relativePosition[1] + relativePosition[2] * relativePosition[2] );

    kernelVal = kernel( r );
    k += Vp[p] * kernelVal;
    for( int i = 0; i < 3; i++ )
    {
      f[i] += Vp[p] * fp[p][i] * kernelVal;
    }

    kernelGradient( x, xp[p], r, kernelGradVal );
    for( int i = 0; i < 3; i++ )
    {
      kGrad[i] += kernelGradVal[i] * Vp[p];
      for( int j = 0; j < 3; j++ )
      {
        fGrad[i][j] += fp[p][i] * kernelGradVal[j] * Vp[p];
      }
    }
  }

  // Return the normalized kernel field gradient (which eliminates edge effects)
  if( k > 0.0 )
  {
    //kernelField = f/k;
    for( int i = 0; i < 3; i++ )
    {
      for( int j = 0; j < 3; j++ )
      {
        result[i][j] = fGrad[i][j] / k - f[i] * kGrad[j] / (k * k);
      }
    }
  }
  else
  {
    //kernelField = 0.0;
    for( int i = 0; i < 3; i++ )
    {
      for( int j = 0; j < 3; j++ )
      {
        result[i][j] = 0.0;
      }
    }
  }
}

void SolidMechanicsMPM::computeDamageFieldGradient( ParticleManager & particleManager )
{
  // Get accessors for volume, position, damage, surface flag
  ParticleManager::ParticleViewAccessor< arrayView1d< real64 const > > particleVolumeAccessor = particleManager.constructArrayViewAccessor< real64, 1 >( "particleVolume" );
  ParticleManager::ParticleViewAccessor< arrayView2d< real64 const > > particlePositionAccessor = particleManager.constructArrayViewAccessor< real64, 2 >( "particleCenter" );
  ParticleManager::ParticleViewAccessor< arrayView1d< real64 const > > particleDamageAccessor = particleManager.constructArrayViewAccessor< real64, 1 >( "particleDamage" );
  ParticleManager::ParticleViewAccessor< arrayView1d< int const > > particleSurfaceFlagAccessor = particleManager.constructArrayViewAccessor< int, 1 >( "particleSurfaceFlag" );

  // Perform neighbor operations
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    // Get neighbor list
    OrderedVariableToManyParticleRelation & neighborList = subRegion.neighborList();
    arrayView1d< localIndex const > const numNeighborsAll = neighborList.m_numParticles.toViewConst();
    ArrayOfArraysView< localIndex const > const neighborRegions = neighborList.m_toParticleRegion.toViewConst();
    ArrayOfArraysView< localIndex const > const neighborSubRegions = neighborList.m_toParticleSubRegion.toViewConst();
    ArrayOfArraysView< localIndex const > const neighborIndices = neighborList.m_toParticleIndex.toViewConst();

    // Get particle position and damage field gradient
    arrayView2d< real64 const > const particlePosition = subRegion.getParticleCenter();
    arrayView2d< real64 > const particleDamageGradient = subRegion.getField< fields::mpm::particleDamageGradient >();
    
    // Loop over neighbors
    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
    forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST ( localIndex const pp ) // Must be on host since we call a 'this'
                                                                                                // method which uses class variables
      {
        localIndex const p = activeParticleIndices[pp];

        localIndex numNeighbors = numNeighborsAll[p];
        arraySlice1d< localIndex const > const regionIndices = neighborRegions[p];
        arraySlice1d< localIndex const > const subRegionIndices = neighborSubRegions[p];
        arraySlice1d< localIndex const > const particleIndices = neighborIndices[p];

        // Declare and size neighbor data arrays - TODO: switch to std::array? But then we'd need to template computeKernelFieldGradient
        std::vector< real64 > neighborVolumes( numNeighbors );
        std::vector< std::vector< real64 > > neighborPositions;
        neighborPositions.resize( numNeighbors, std::vector< real64 >( 3 ) );
        std::vector< real64 > neighborDamages( numNeighbors );

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
          if( particleSurfaceFlagAccessor[regionIndex][subRegionIndex][particleIndex] == 1 )
          {
            neighborDamages[neighborIndex] = 1.0;
          }
          else
          {
            neighborDamages[neighborIndex] = particleDamageAccessor[regionIndex][subRegionIndex][particleIndex];
          }
        }


        // Call kernel field gradient function
        computeKernelFieldGradient( particlePosition[p],        // input
                                    neighborPositions,          // input
                                    neighborVolumes,            // input
                                    neighborDamages,            // input
                                    particleDamageGradient[p] ); // OUTPUT
      } );
  } );
}

void SolidMechanicsMPM::updateSurfaceFlagOverload( ParticleManager & particleManager )
{
  // Surface flags are overloaded so we can visualize surfaces and damage features simultaneously
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    arrayView1d< int > const particleSurfaceFlag = subRegion.getParticleSurfaceFlag();
    arrayView1d< real64 const > const particleDamage = subRegion.getParticleDamage();
    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
    forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
    {
      localIndex const p = activeParticleIndices[pp];
      if( particleSurfaceFlag[p] != 2 )
      {
        if( particleDamage[p] > 0.0 ) // Activate damage field if any particles in domain have damage.
        {
          particleSurfaceFlag[p] = 1;
        }
        else
        {
          particleSurfaceFlag[p] = 0;
        }
      }
    } );
  } );
}

void SolidMechanicsMPM::projectDamageFieldGradientToGrid( ParticleManager & particleManager,
                                                          NodeManager & nodeManager )
{
  // Grid nodes gain the damage field gradient of the particle mapping to them with the largest damage field gradient

  // Get grid fields
  arrayView2d< real64 > const gridDamageGradient = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::damageGradientString() );
  int subRegionIndex = 0;
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    // Get particle fields
    arrayView2d< real64 const > const particleDamageGradient = subRegion.getField< fields::mpm::particleDamageGradient >();

    // Get nodes this particle maps to
    arrayView2d< localIndex const > const mappedNodes = m_mappedNodes[subRegionIndex];

    // Map to grid
    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
    int const numDims = m_numDims;
    forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST ( localIndex const pp ) // Parallelize with atomics/reduction
      {
        localIndex const p = activeParticleIndices[pp];

        // Map to grid
        for( localIndex const & g: mappedNodes[pp] )
        {
          if( LvArray::tensorOps::l2NormSquared< 3 >( particleDamageGradient[p] ) > LvArray::tensorOps::l2NormSquared< 3 >( gridDamageGradient[g] ) )
          {
            for( int i=0; i<numDims; i++ )
            {
              gridDamageGradient[g][i] = particleDamageGradient[p][i];
            }
          }
        }
      } ); // particle loop

    // Increment subregion index
    subRegionIndex++;
  } ); // subregion loop
}

void SolidMechanicsMPM::updateDeformationGradient( real64 dt,
                                                   ParticleManager & particleManager )
{
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    // Get fields
    arrayView3d< real64 > const particleDeformationGradient = subRegion.getField< fields::mpm::particleDeformationGradient >();
    arrayView3d< real64 > const particleFDot = subRegion.getField< fields::mpm::particleFDot >();
    arrayView3d< real64 const > const particleVelocityGradient = subRegion.getField< fields::mpm::particleVelocityGradient >();

    // Update F
    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
    forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
    {
      localIndex const p = activeParticleIndices[pp];
      LvArray::tensorOps::Rij_eq_AikBkj< 3, 3, 3 >( particleFDot[p], particleVelocityGradient[p], particleDeformationGradient[p] ); // Fdot
                                                                                                                                    // = L.F
      LvArray::tensorOps::scaledAdd< 3, 3 >( particleDeformationGradient[p], particleFDot[p], dt ); // Fnew = Fold + Fdot*dt
    } );
  } );
}

void SolidMechanicsMPM::updateConstitutiveModelDependencies( ParticleManager & particleManager )
{
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    // Get needed particle fields
    arrayView1d< real64 const > const particleVolume = subRegion.getParticleVolume();
    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();

    // Get constitutive model reference
    string const & solidMaterialName = subRegion.template getReference< string >( viewKeyStruct::solidMaterialNamesString() );
    SolidBase & solidModel = getConstitutiveModel< SolidBase >( subRegion, solidMaterialName );

    // Pass whatever data the constitutive models may need
    if( solidModel.hasWrapper( "lengthScale" ) ) // Fragile code because someone could change this key without our knowledge. TODO: Make an
                                                 // integrated test that checks this
    {
      arrayView1d< real64 > const lengthScale = solidModel.getReference< array1d< real64 > >( "lengthScale" );
      forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
      {
        localIndex const p = activeParticleIndices[pp];
        lengthScale[p] = pow( particleVolume[p], 1.0 / 3.0 );
      } );
    }

    if(  solidModel.hasWrapper( "materialDirection" ) )
    {
      // CC: Todo add check for fiber vs plane update to material direction
      arrayView2d< real64 > const particleMaterialDirection = subRegion.getParticleMaterialDirection();
      arrayView2d< real64 > const constitutiveMaterialDirection = solidModel.getReference< array2d< real64 > >( "materialDirection" );
      forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
      {
        localIndex const p = activeParticleIndices[pp];
        LvArray::tensorOps::copy< 3 >(constitutiveMaterialDirection[p], particleMaterialDirection[p]); 
      } );
    }

    if(  solidModel.hasWrapper( "deformationGradient" ) )
    {
      arrayView3d< real64 > const particleDeformationGradient = subRegion.getField< fields::mpm::particleDeformationGradient >();
      arrayView3d< real64 > const constitutiveDeformationGradient = solidModel.getReference< array3d< real64 > >( "deformationGradient" );
      forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
      {
        localIndex const p = activeParticleIndices[pp];
        LvArray::tensorOps::copy< 3, 3 >(constitutiveDeformationGradient[p], particleDeformationGradient[p]); 
      } );
    }

    if(  solidModel.hasWrapper( "velocityGradient" ) )
    {
      arrayView3d< real64 > const particleVelocityGradient = subRegion.getField< fields::mpm::particleVelocityGradient >();
      arrayView3d< real64 > const constitutiveVelocityGradient = solidModel.getReference< array3d< real64 > >( "velocityGradient" );
      forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
      {
        localIndex const p = activeParticleIndices[pp];
        LvArray::tensorOps::copy< 3, 3 >(constitutiveVelocityGradient[p], particleVelocityGradient[p]); 
      } );
    }
  } );
}

void SolidMechanicsMPM::updateStress( real64 dt,
                                      ParticleManager & particleManager )
{
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    // Get constitutive model reference
    string const & solidMaterialName = subRegion.template getReference< string >( viewKeyStruct::solidMaterialNamesString() );
    SolidBase & solid = getConstitutiveModel< SolidBase >( subRegion, solidMaterialName );

    // Get particle kinematic fields that are fed into constitutive model
    arrayView3d< real64 const > const particleDeformationGradient = subRegion.getField< fields::mpm::particleDeformationGradient >();
    arrayView3d< real64 const > const particleFDot = subRegion.getField< fields::mpm::particleFDot >();
    arrayView3d< real64 const > const particleVelocityGradient = subRegion.getField< fields::mpm::particleVelocityGradient >();
    arrayView2d< real64 > const particleStress = subRegion.getField< fields::mpm::particleStress >();

    int hyperelasticUpdate = 0;
    if ( solid.getCatalogName() == "HyperelasticMMS" || solid.getCatalogName() == "Hyperelastic" )
    {
      hyperelasticUpdate = 1;
    }

    // Call constitutive model
    ConstitutivePassThruMPM< SolidBase >::execute( solid, [&] ( auto & castedSolid )
    {
      using SolidType = TYPEOFREF( castedSolid );
      typename SolidType::KernelWrapper constitutiveModelWrapper = castedSolid.createKernelUpdates();
      solidMechanicsMPMKernels::StateUpdateKernel::launch< serialPolicy >( subRegion.activeParticleIndices(),
                                                                           constitutiveModelWrapper,
                                                                           dt,
                                                                           hyperelasticUpdate,
                                                                           particleDeformationGradient,
                                                                           particleFDot,
                                                                           particleVelocityGradient,
                                                                           particleStress );
    } );
  } );
}

void SolidMechanicsMPM::particleKinematicUpdate( ParticleManager & particleManager )
{
  // Update particle volume and density
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    // Get particle fields
    arrayView1d< globalIndex const > const particleID = subRegion.getParticleID();
    arrayView1d< real64 > const particleVolume = subRegion.getParticleVolume();
    arrayView1d< int > const isBad = subRegion.getField< fields::mpm::isBad >();
    arrayView1d< real64 const > const particleInitialVolume = subRegion.getField< fields::mpm::particleInitialVolume >();
    arrayView1d< real64 > const particleDensity = subRegion.getField< fields::mpm::particleDensity >();
    arrayView1d< real64 const > const particleMass = subRegion.getField< fields::mpm::particleMass >();
    arrayView3d< real64 > const particleRVectors = subRegion.getParticleRVectors();
    arrayView3d< real64 const > const particleInitialRVectors = subRegion.getField< fields::mpm::particleInitialRVectors >();
    arrayView3d< real64 const > const particleDeformationGradient = subRegion.getField< fields::mpm::particleDeformationGradient >();
    arrayView2d< real64 const > const particleInitialMaterialDirection = subRegion.getParticleInitialMaterialDirection();
    arrayView2d< real64 > const particleMaterialDirection = subRegion.getParticleMaterialDirection();

    // Get constitutive model reference to check if material represents a fiber
    string const & solidMaterialName = subRegion.template getReference< string >( viewKeyStruct::solidMaterialNamesString() );
    SolidBase & solidModel = getConstitutiveModel< SolidBase >( subRegion, solidMaterialName );
    bool isFiber = solidModel.hasWrapper( "isFiber" ); // dumby variable whose value doesn't matter only that it is defined ( possibly a better way to do this )

    // Update volume and r-vectors
    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
    forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
    {
      localIndex const p = activeParticleIndices[pp];
      real64 detF = LvArray::tensorOps::determinant< 3 >( particleDeformationGradient[p] );

      if( detF <= 0.1 || detF >= 10.0 )
      {
        printf( "Flagging particle with unreasonable Jacobian (J<0.1 or J>10) for deletion! Global particle ID: %lld", particleID[p] );
        isBad[p] = 1; // TODO: Switch to a 'markBad' function which changes isBad and removes this particle from activeParticleIndices. Also
                      // change isBad to deletionFlag.
        particleVolume[p] = particleInitialVolume[p];
        particleDensity[p] = particleMass[p] / particleInitialVolume[p];
        LvArray::tensorOps::copy< 3, 3 >( particleRVectors[p], particleInitialRVectors[p] );
      }
      else
      {
        particleVolume[p] = particleInitialVolume[p] * detF;
        particleDensity[p] = particleMass[p] / particleVolume[p];

        // Update material direction (CC: TODO add check for material direction for plane and fiber, currently plane)
        real64 materialDirection[3];
        real64 deformationGradient[3][3];
        LvArray::tensorOps::copy< 3, 3 >(deformationGradient, particleDeformationGradient[p]);

        if( isFiber )
        {
          LvArray::tensorOps::Ri_eq_AijBj< 3, 3 >( materialDirection, deformationGradient, particleInitialMaterialDirection[p] );
          LvArray::tensorOps::copy< 3 >(particleMaterialDirection[p], materialDirection);
        } 
        else
        {
          real64 deformationGradientCofactor[3][3];
          cofactor(deformationGradient, deformationGradientCofactor);
          LvArray::tensorOps::Ri_eq_AijBj< 3, 3 >( materialDirection, deformationGradientCofactor, particleInitialMaterialDirection[p] );
          LvArray::tensorOps::copy< 3 >(particleMaterialDirection[p], materialDirection);
        } 
      }
    } );
  } );

  // Compute particles R vectors
  computeRVectors( particleManager );
}

void SolidMechanicsMPM::computeAndWriteBoxAverage( const real64 dt,
                                                   const real64 time_n,
                                                   ParticleManager & particleManager )
{
  real64 boxPlasticStrain[6] = { 0.0 };
  real64 boxStress[6] = { 0.0 }; // we sum stress * volume in particles, additive sync, then divide by box volume.
  real64 boxMass = 0.0; // we sum particle mass, additive sync, then divide by box volume
  real64 boxParticleInitialVolume = 0.0;
  real64 boxDamage = 0.0; // we sum damage * initial volume, additive sync, then divide by total initial volume in box

  real64 boxAverageMin[3];
  LvArray::tensorOps::copy< 3 >( boxAverageMin, m_boxAverageMin );
  real64 boxAverageMax[3];
  LvArray::tensorOps::copy< 3 >( boxAverageMax, m_boxAverageMax );

  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    // Get fields
    arrayView2d< real64 const > const particlePosition = subRegion.getParticleCenter();
    arrayView1d< real64 const > const particleMass = subRegion.getField< fields::mpm::particleMass >();
    arrayView1d< real64 const > const particleVolume = subRegion.getParticleVolume();
    arrayView1d< real64 const > const particleInitialVolume = subRegion.getField< fields::mpm::particleInitialVolume >();
    arrayView2d< real64 const > const particleStress = subRegion.getField< fields::mpm::particleStress >();
    arrayView2d< real64 const > const particlePlasticStrain = subRegion.getField< fields::mpm::particlePlasticStrain >();
    arrayView1d< real64 const > const particleDamage = subRegion.getParticleDamage();

    // Accumulate values
    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
    forAll< serialPolicy >( activeParticleIndices.size(), [=, &boxMass, &boxParticleInitialVolume, &boxStress, &boxPlasticStrain, &boxDamage] GEOS_HOST ( localIndex const pp ) // This
                                                                                                                                                                                // can
                                                                                                                                                                                // be
                                                                                                                                                                                // parallelized
                                                                                                                                                                                // via
                                                                                                                                                                                // reduction
    {
      localIndex const p = activeParticleIndices[pp];
      real64 x = particlePosition[p][0];
      real64 y = particlePosition[p][1];
      real64 z = particlePosition[p][2];
      if( x > boxAverageMin[0] && x < boxAverageMax[0] && 
          y > boxAverageMin[1] && y < boxAverageMax[1] && 
          z > boxAverageMin[2] && z < boxAverageMax[2] )
      {
        boxMass += particleMass[p];
        boxParticleInitialVolume += particleInitialVolume[p];
        for( int i=0; i<6; i++ )
        {
          boxStress[i] += particleStress[p][i] * particleVolume[p]; // volume weighted average, will normalize later.
          boxPlasticStrain[i] += particlePlasticStrain[p][i] * particleVolume[p];
        }
        boxDamage += particleDamage[p] * particleInitialVolume[p]; // initial volume weighted average, will normalize later.
      }
    } );
  } );

  // Additive sync: sxx, syy, szz, sxy, syz, sxz, mass, particle volume, damage
  // Check the voigt indexing of stress
  real64 boxSums[15];
  boxSums[0] = boxStress[0];       // sig_xx * volume
  boxSums[1] = boxStress[1];       // sig_yy * volume
  boxSums[2] = boxStress[2];       // sig_zz * volume
  boxSums[3] = boxStress[3];       // sig_yz * volume
  boxSums[4] = boxStress[4];       // sig_xz * volume
  boxSums[5] = boxStress[5];       // sig_xy * volume
  boxSums[6] = boxMass;            // total mass in box
  boxSums[7] = boxParticleInitialVolume > 0.0 ? boxParticleInitialVolume : 1.0;  // total particle initial volume in box; prevent div0 error
  boxSums[8] = boxDamage;          // damage * volume
  boxSums[9] = boxPlasticStrain[0]; // plasticStrain_xx
  boxSums[10] = boxPlasticStrain[1]; // plasticStrain_yy
  boxSums[11] = boxPlasticStrain[2]; // plasticStrain_zz
  boxSums[12] = boxPlasticStrain[3]; // plasticStrain_yz
  boxSums[13] = boxPlasticStrain[4]; // plasticStrain_xz
  boxSums[14] = boxPlasticStrain[5]; // plasticStrain_xy
      
  // Do an MPI sync to total these values and write from proc0 to a file.  Also compute global F
  // so file is directly plottable in excel as CSV or something.
  for( localIndex i = 0; i < 15; i++ )
  {
    real64 localSum = boxSums[i];
    real64 globalSum;
    MPI_Allreduce( &localSum,
                   &globalSum,
                   1,
                   MPI_DOUBLE,
                   MPI_SUM,
                   MPI_COMM_GEOSX );
    boxSums[i] = globalSum;
  }

  int rank;
  MPI_Comm_rank( MPI_COMM_GEOSX, &rank );
  if( rank == 0 )
  {
    // Calculate the box volume
    real64 boxVolume = m_domainExtent[0] * m_domainExtent[1] * m_domainExtent[2];

    //CC: Old GEOS logic for computing box sums 
    // May not need this if an extra cell is not already added as was the case when periodicity was specificed in the xml in old GEOS
    // Now we normalize stress, which is stored as a volume-weighted value:
    // realT boxVolume;
    // if( m_prescribed_boundary_f_table || m_prescribed_f_table )
    // {
    // boxVolume = ( m_xmax_global - m_xmin_global - (!periodic[0])*2.0*m_dx) * ( m_ymax_global - m_ymin_global - (!periodic[1])*2.0*m_dy) * ( m_zmax_global - m_zmin_global - (!periodic[2])*2.0*m_dz);
    //   //boxVolume = ( m_xmax_global - m_xmin_global - 2.0*m_dx) * ( m_ymax_global - m_ymin_global - 2.0*m_dy) * ( m_zmax_global - m_zmin_global - 2.0*m_dz);
    // }
    // else
    // {
    //   boxVolume = ( m_particle_box_xmax - m_particle_box_xmin ) * ( m_particle_box_ymax - m_particle_box_ymin ) * ( m_particle_box_zmax - m_particle_box_zmin ); //Does this need to account for periodic boundaries
    // }

     // Write to file
    std::ofstream file;
    file.open( "boxAverageHistory.csv", std::ios::out | std::ios::app );
    if( file.fail() )
    {
      throw std::ios_base::failure( std::strerror( errno ) );
    }
    //make sure write fails with exception if something is wrong
    file.exceptions( file.exceptions() | std::ios::failbit | std::ifstream::badbit );
    // time | sig_xx | sig_yy | sig_zz | sig_xy | sig_yz | sig_zx | density | damage / total particle volume
    file << time_n + dt
         << ","
         << boxSums[0] / boxVolume
         << ","
         << boxSums[1] / boxVolume
         << ","
         << boxSums[2] / boxVolume
         << ","
         << boxSums[3] / boxVolume
         << ","
         << boxSums[4] / boxVolume
         << ","
         << boxSums[5] / boxVolume
         << ","
         << boxSums[6] / boxVolume
         << ","
         << boxSums[8] / boxSums[7] // We normalize by total particle initial volume because this should equal one if all the material is
                                    // damaged
         << ", "
         << boxSums[9] / boxVolume
         << ", "
         << boxSums[10] / boxVolume
         << ", "
         << boxSums[11] / boxVolume
         << ", "
         << boxSums[12] / boxVolume
         << ", "
         << boxSums[13] / boxVolume
         << ", "
         << boxSums[14] / boxVolume
         << std::endl;
    file.close();
  }
}


void SolidMechanicsMPM::computeBoxStress( ParticleManager & particleManager,
                                          arrayView1d< real64 > currentStress )
{
  real64 boxStress[3] = { 0.0 }; // we sum stress * volume in particles, additive sync, then divide by box volume.

  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    // Get fields
    arrayView1d< real64 > const particleVolume = subRegion.getParticleVolume();
    arrayView2d< real64 > const particleStress = subRegion.getField< fields::mpm::particleStress >();

    // Accumulate values
    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();

    // CC: TODO parallelize via reduction
    forAll< serialPolicy >( activeParticleIndices.size(), [=, &boxStress] GEOS_HOST ( localIndex const pp )
    {
      localIndex const p = activeParticleIndices[pp];
      // Only need diagonal components
      for( int i=0; i < 3; i++ )
      {
        boxStress[i] += particleStress[p][i] * particleVolume[p]; // volume weighted average, will normalize later.
      }
    } );
  } );

  // Additive sync: sxx, syy, szz, sxy, syz, sxz, mass, particle volume, damage
  real64 boxSums[3];
  boxSums[0] = boxStress[0];       // sig_xx * volume
  boxSums[1] = boxStress[1];       // sig_yy * volume
  boxSums[2] = boxStress[2];       // sig_zz * volume

  // Do an MPI sync to total these values and write from proc0 to a file.  Also compute global F
  // so file is directly plottable in excel as CSV or something.
  for( localIndex i = 0; i < 3; i++ )
  {
    real64 localSum = boxSums[i];
    real64 globalSum;
    MPI_Allreduce( &localSum,
                   &globalSum,
                   1,
                   MPI_DOUBLE,
                   MPI_SUM,
                   MPI_COMM_GEOSX );
    boxSums[i] = globalSum;
  }

  real64 boxVolume = m_domainExtent[0] * m_domainExtent[1] * m_domainExtent[2];

  currentStress[0] = boxSums[0] / boxVolume;
  currentStress[1] = boxSums[1] / boxVolume;
  currentStress[2] = boxSums[2] / boxVolume;  
}

void SolidMechanicsMPM::stressControl( real64 dt,
                                       ParticleManager & particleManager )
{
  
  real64 targetStress[3] = {0};
  LvArray::tensorOps::copy< 3 >(targetStress, m_domainStress);

  array1d< real64 > currentStress;
  currentStress.resize( 3 );
	computeBoxStress( particleManager,
                    currentStress );

  // Uses maximum bulk modulus ( lowest effective PID gains ) of all materials
  real64 maximumBulkModulus = 0.0;
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    string const & solidMaterialName = subRegion.template getReference< string >( viewKeyStruct::solidMaterialNamesString() );
    const SolidBase & solidModel = getConstitutiveModel< SolidBase >( subRegion, solidMaterialName );
    
    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();

    array1d< real64 > bulkModulus;
    string constitutiveModelName = solidModel.getCatalogName();
    if( constitutiveModelName == "Hyperelastic" ){
      const Hyperelastic & hyperelastic = dynamic_cast< const Hyperelastic & >( solidModel );
      bulkModulus = hyperelastic.bulkModulus();
    }

    if( constitutiveModelName == "HyperelasticMMS" ){
      const HyperelasticMMS & hyperelasticMMS = dynamic_cast< const HyperelasticMMS & >( solidModel ); 
      arrayView1d< real64 const > const lambda = hyperelasticMMS.lambda();
      arrayView1d< real64 const > const shearModulus = hyperelasticMMS.shearModulus();
      bulkModulus.resize(lambda.size());
      forAll< serialPolicy >( activeParticleIndices.size(), [=, &bulkModulus] GEOS_HOST ( localIndex const pp ) // Could be reduction instead of a loop
      {
        localIndex const p = activeParticleIndices[pp];
        bulkModulus[p] = conversions::lameConstants::toBulkMod( lambda[p], shearModulus[p] );
      } );
    }

    if( constitutiveModelName == "ElasticIsotropic" || constitutiveModelName == "CeramicDamage" || constitutiveModelName == "StrainHardeningPolymer" ){
      const ElasticIsotropic & elasticIsotropic = dynamic_cast< const ElasticIsotropic & >( solidModel );
      bulkModulus = elasticIsotropic.bulkModulus();
    }

    if( constitutiveModelName == "Graphite" || constitutiveModelName == "ElasticTransverseIsotropic" || constitutiveModelName == "ElasticTransverseIsotropicPressureDependent" ){
      const ElasticTransverseIsotropic & elasticTransverseIsotropic = dynamic_cast< const ElasticTransverseIsotropic & >( solidModel );
      bulkModulus = elasticTransverseIsotropic.effectiveBulkModulus();
    }

    forAll< serialPolicy >( activeParticleIndices.size(), [=, &maximumBulkModulus] GEOS_HOST ( localIndex const pp ) // Could be reduction instead of a loop
    {
      localIndex const p = activeParticleIndices[pp];
      maximumBulkModulus = fmax( maximumBulkModulus, bulkModulus[p] );
    } );
  } );

  // Do a reduce of maximumBulkModulus between ranks so stress control calculations are consistent
  real64 maximumBulkModulusThisRank = maximumBulkModulus;
  real64 maximumBulkModulusAllRanks = 0;
  MpiWrapper::allReduce< real64 >( &maximumBulkModulusThisRank,
                                   &maximumBulkModulusAllRanks,
                                   1,
                                   MPI_MAX,
                                   MPI_COMM_GEOSX );
  maximumBulkModulus = maximumBulkModulusAllRanks;

  GEOS_ERROR_IF( maximumBulkModulus <= 0.0, "At least one material must have a positive bulk or effective bulk modulus for stress control!" );

	// This will drive the response towards the desired stress but may be
	// unstable.
	real64 stressControlKp = m_stressControlKp / ( maximumBulkModulus * dt );
	real64 stressControlKd = m_stressControlKd / ( maximumBulkModulus );
	real64 stressControlKi = m_stressControlKi / ( maximumBulkModulus * dt * dt );

	real64 error[3] = {0};
  LvArray::tensorOps::copy< 3 >(error, targetStress);
  LvArray::tensorOps::subtract< 3 >(error, currentStress);

	real64 dedt[3] = {0};
  LvArray::tensorOps::copy< 3 >( dedt, error);
  LvArray::tensorOps::subtract< 3 >( dedt, m_stressControlLastError);
  LvArray::tensorOps::scale< 3 >( dedt, dt);

	real64 stressControlPTerm[3] = {0}; 
  LvArray::tensorOps::copy< 3 >( stressControlPTerm, error);
  LvArray::tensorOps::scale< 3 >( stressControlPTerm, stressControlKp);

  real64 stressControlITerm[3] = {0};
  LvArray::tensorOps::copy< 3 >( stressControlITerm, error);
  LvArray::tensorOps::scale< 3 >( stressControlITerm, stressControlKi*dt);
  LvArray::tensorOps::add< 3 >( m_stressControlITerm, stressControlITerm );

  real64 stressControlDTerm[3] = {0};
  LvArray::tensorOps::copy< 3 >( stressControlDTerm, dedt);
  LvArray::tensorOps::scale< 3 >(stressControlDTerm, stressControlKd);

	real64 L_new[3] = {0}; 
  LvArray::tensorOps::copy< 3 >( L_new, stressControlPTerm );
  LvArray::tensorOps::add< 3 >( L_new, m_stressControlITerm );
  LvArray::tensorOps::add< 3 >( L_new, stressControlDTerm );

	LvArray::tensorOps::copy< 3 >(m_stressControlLastError, error);

	// Limit L so boundary motion doesn't exceed CFL limit
	// boundaryV = strainRate*(m_xmax_global - m_xmin_global) = CFL*m_dx/dt
	real64 Lxx_min, 
         Lxx_max, 
         Lyy_min, 
         Lyy_max, 
         Lzz_min,
         Lzz_max;

	Lxx_max = m_cflFactor * ( m_hEl[0] / dt) / ( m_xGlobalMax[0] - m_xGlobalMin[0]);
	Lxx_min = -1.0 * Lxx_max;
	L_new[0] = std::min( Lxx_max , std::max( L_new[0], Lxx_min ) );

	Lyy_max = m_cflFactor * ( m_hEl[1] / dt) / ( m_xGlobalMax[1] - m_xGlobalMin[1]);
	Lyy_min = -1.0 * Lyy_max;
	L_new[1] = std::min( Lyy_max , std::max( L_new[1], Lyy_min ) );

	Lzz_max = m_cflFactor * ( m_hEl[2] / dt) / ( m_xGlobalMax[2] - m_xGlobalMin[2]);
	Lzz_min = -1.0 * Lzz_max;
	L_new[2] = std::min( Lzz_max , std::max( L_new[2], Lzz_min ) );

	for(int i=0; i < m_numDims; i++ )
	{
		if( m_stressControl[i] == 1 )
		{
			m_domainL[i] = L_new[i];
			m_domainF[i] += m_domainL[i]*m_domainF[i]*dt;
		}
	}

  // // CC: debug
  // GEOS_LOG_RANK_0("error: " << error[0] << ", " << error[1] << ", " << error[2] << 
  //                 ", dedt: " << dedt[0] << ", " << dedt[1] << ", " << dedt[2] <<  
  //                 ", P: " << stressControlPTerm[0] << ", " << stressControlPTerm[1] << ", " << stressControlPTerm[2] << 
  //                 ", I: " << stressControlITerm[0] << ", " << stressControlITerm[1] << ", " << stressControlITerm[2] <<  
  //                 ", D: " << stressControlDTerm[0] << ", " << stressControlDTerm[1] << ", " << stressControlDTerm[2]);
  // GEOS_LOG_RANK_0("Domain L: " << m_domainL << ", Domain F" << m_domainF);
}


void SolidMechanicsMPM::applySuperimposedVelocityGradient( const real64 dt, 
                                                           ParticleManager & particleManager,
                                                           SpatialPartition & partition )
{
  real64 domainL[3] = {0};
  LvArray::tensorOps::copy< 3 >( domainL, m_domainL );
  int const numDims = m_numDims; // CC: do member scalars need to be copied to local variable to be used in a RAJA loops?

  arrayView1d< int const > const periodic = partition.getPeriodic();
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    // Particle fields
    arrayView2d< real64 > const particlePosition = subRegion.getParticleCenter();
    arrayView3d< real64 > const particleVelocityGradient = subRegion.getField< fields::mpm::particleVelocityGradient >();

    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
    forAll< parallelDevicePolicy<> >( activeParticleIndices.size(), [=] GEOS_DEVICE ( localIndex const pp )
       {
        localIndex const p = activeParticleIndices[pp];
        
        for(int i=0; i < numDims; i++)
        {
          if( periodic[i] )
          {
            particleVelocityGradient[p][i][i] += domainL[i];
            particlePosition[p][i] += particlePosition[p][i] * domainL[i] * dt;
          }
        }
      } ); // particle loop
  } ); // subregion loop
}

void SolidMechanicsMPM::initializeGridFields( NodeManager & nodeManager )
{
  int const numNodes = nodeManager.size();
  arrayView2d< real64 > const gridMass = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::massString() );
  arrayView2d< real64 > const gridDamage = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::damageString() );
  arrayView2d< real64 > const gridMaxDamage = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::maxDamageString() );
  arrayView2d< real64 > const gridDamageGradient = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::damageGradientString() );
  arrayView3d< real64 > const gridVelocity = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::velocityString() );
  arrayView3d< real64 > const gridMomentum = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::momentumString() );
  arrayView3d< real64 > const gridAcceleration = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::accelerationString() );
  arrayView3d< real64 > const gridInternalForce = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::forceInternalString() );
  arrayView3d< real64 > const gridExternalForce = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::forceExternalString() );
  arrayView3d< real64 > const gridContactForce = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::forceContactString() );
  arrayView3d< real64 > const gridSurfaceNormal = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::surfaceNormalString() );
  arrayView3d< real64 > const gridMaterialPosition = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::materialPositionString() );

  forAll< serialPolicy >( numNodes, [=] GEOS_HOST ( localIndex const g ) // Switch to .zero()?
    {
      for( int i = 0; i < 3; i++ )
      {
        gridDamageGradient[g][i] = 0.0;
      }
      for( int fieldIndex = 0; fieldIndex < m_numVelocityFields; fieldIndex++ )
      {
        gridMass[g][fieldIndex] = 0.0;
        gridDamage[g][fieldIndex] = 0.0;
        gridMaxDamage[g][fieldIndex] = 0.0;
        for( int i = 0; i < 3; i++ )
        {
          gridVelocity[g][fieldIndex][i] = 0.0;
          gridMomentum[g][fieldIndex][i] = 0.0;
          gridAcceleration[g][fieldIndex][i] = 0.0;
          gridInternalForce[g][fieldIndex][i] = 0.0;
          gridExternalForce[g][fieldIndex][i] = 0.0;
          gridContactForce[g][fieldIndex][i] = 0.0;
          gridSurfaceNormal[g][fieldIndex][i] = 0.0;
          gridMaterialPosition[g][fieldIndex][i] = 0.0;
        }
      }
    } );
}

void SolidMechanicsMPM::boundaryConditionUpdate( real64 dt, real64 time_n )
{
  int bcInterval = 0;

  for( localIndex i = 0; i < m_bcTable.size( 0 ); i++ ) // Naive method for determining what part of BC table we're currently in, can
                                                        // optimize later (TODO)
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


void SolidMechanicsMPM::particleToGrid( ParticleManager & particleManager,
                                        NodeManager & nodeManager )
{
  localIndex subRegionIndex = 0;
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    // Particle fields
    // ParticleType particleType = subRegion.getParticleType(); // CC: unused?
    arrayView2d< real64 const > const particlePosition = subRegion.getParticleCenter();
    arrayView2d< real64 const > const particleVelocity = subRegion.getParticleVelocity();
    arrayView1d< real64 const > const particleMass = subRegion.getField< fields::mpm::particleMass >();
    arrayView1d< real64 const > const particleVolume = subRegion.getParticleVolume();
    arrayView1d< int const > const particleGroup = subRegion.getParticleGroup();
    arrayView1d< int const > const particleSurfaceFlag = subRegion.getParticleSurfaceFlag();

    arrayView2d< real64 const > const particleStress = subRegion.getField< fields::mpm::particleStress >();
    arrayView2d< real64 const > const particleBodyForce = subRegion.getField< fields::mpm::particleBodyForce >();
    arrayView2d< real64 const > const particleDamageGradient = subRegion.getField< fields::mpm::particleDamageGradient >();
    arrayView1d< real64 const > const particleDamage = subRegion.getParticleDamage();

    // Grid fields
    arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const gridPosition = nodeManager.referencePosition();
    arrayView2d< real64 > const gridMass = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::massString() );
    arrayView2d< real64 const > const gridDamageGradient = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::damageGradientString() );
    arrayView3d< real64 > const gridMomentum = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::momentumString() );
    arrayView3d< real64 > const gridInternalForce = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::forceInternalString() );
    arrayView3d< real64 > const gridExternalForce = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::forceExternalString() );
    arrayView3d< real64 > const gridMaterialPosition = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::materialPositionString() );
    arrayView2d< real64 > const gridDamage = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::damageString() );
    arrayView2d< real64 > const gridMaxDamage = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::maxDamageString() );

    // Get views to mapping arrays
    int const numberOfVerticesPerParticle = subRegion.numberOfVerticesPerParticle();
    arrayView2d< localIndex const > const mappedNodes = m_mappedNodes[subRegionIndex];
    arrayView2d< real64 const > const shapeFunctionValues = m_shapeFunctionValues[subRegionIndex];
    arrayView3d< real64 const > const shapeFunctionGradientValues = m_shapeFunctionGradientValues[subRegionIndex];

    // Map to grid
    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
    int const numDims = m_numDims;
    int voigtMap[3][3] = { {0, 5, 4}, {5, 1, 3}, {4, 3, 2} };
    int const damageFieldPartitioning = m_damageFieldPartitioning;

    forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST ( localIndex const pp ) // Can be parallized using atomics -
                                                                                                // remember to pass copies of class
                                                                                                // variables
    {                                                                                          // Grid max damage will require reduction
      localIndex const p = activeParticleIndices[pp];

      for( int g = 0; g < 8 * numberOfVerticesPerParticle; g++ )
      {
        localIndex const mappedNode = mappedNodes[pp][g];
        int const nodeFlag = ( damageFieldPartitioning == 1 && LvArray::tensorOps::AiBi< 3 >( gridDamageGradient[mappedNode], particleDamageGradient[p] ) < 0.0 ) ? 1 : 0; // 0
                                                                                                                                                                            // undamaged
                                                                                                                                                                            // or
                                                                                                                                                                            // "A"
                                                                                                                                                                            // field,
                                                                                                                                                                            // 1
                                                                                                                                                                            // for
                                                                                                                                                                            // "B"
                                                                                                                                                                            // field
        int const fieldIndex = nodeFlag * m_numContactGroups + particleGroup[p]; // This ranges from 0 to nMatFields-1

        gridMass[mappedNode][fieldIndex] += particleMass[p] * shapeFunctionValues[pp][g];
        // TODO: Normalizing by volume might be better
        gridDamage[mappedNode][fieldIndex] += particleMass[p] * ( particleSurfaceFlag[p] == 1 ? 1 : particleDamage[pp] ) * shapeFunctionValues[pp][g];
        gridMaxDamage[mappedNode][fieldIndex] = fmax( gridMaxDamage[mappedNode][fieldIndex], particleSurfaceFlag[p] == 1 ? 1 : particleDamage[pp] );
        for( int i=0; i<numDims; i++ )
        {
          gridMomentum[mappedNode][fieldIndex][i] += particleMass[p] * particleVelocity[p][i] * shapeFunctionValues[pp][g];
          gridExternalForce[mappedNode][fieldIndex][i] += particleBodyForce[p][i] * particleMass[p] * shapeFunctionValues[pp][g];
          
          // TODO: Switch to volume weighting?
          gridMaterialPosition[mappedNode][fieldIndex][i] += particleMass[p] * (particlePosition[p][i] - gridPosition[mappedNode][i]) * shapeFunctionValues[pp][g];
          for( int k=0; k<numDims; k++ )
          {
            int voigt = voigtMap[k][i];
            gridInternalForce[mappedNode][fieldIndex][i] -= particleStress[p][voigt] * shapeFunctionGradientValues[pp][g][k] * particleVolume[p];
          }
        }
      }
    } ); // particle loop

    // Increment subregion index
    subRegionIndex++;
  } ); // subregion loop
}

void SolidMechanicsMPM::gridTrialUpdate( real64 dt,
                                         NodeManager & nodeManager )
{
  // Grid fields
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const gridPosition = nodeManager.referencePosition();
  arrayView2d< real64 const > const & gridMass = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::massString() );
  arrayView3d< real64 > const & gridVelocity = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::velocityString() );
  arrayView3d< real64 > const & gridMomentum = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::momentumString() );
  arrayView3d< real64 > const & gridAcceleration = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::accelerationString() );
  arrayView3d< real64 const > const & gridInternalForce = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::forceInternalString() );
  arrayView3d< real64 const > const & gridExternalForce = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::forceExternalString() );
  arrayView3d< real64 > const & gridMaterialPosition = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::materialPositionString() );
  arrayView2d< real64 > const gridDamage = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::damageString() );

  // Loop over velocity fields
  int const numNodes = nodeManager.size();
  real64 const smallMass = m_smallMass;
  int const numDims = m_numDims;
  for( int fieldIndex=0; fieldIndex<m_numVelocityFields; fieldIndex++ )
  {
    forAll< serialPolicy >( numNodes, [=] GEOS_HOST_DEVICE ( localIndex const g )
    {
      if( gridMass[g][fieldIndex] > smallMass ) // small mass threshold
      {
        gridDamage[g][fieldIndex] /= gridMass[g][fieldIndex];
        for( int i=0; i<numDims; i++ )
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
        gridDamage[g][fieldIndex] = 0.0;
        for( int i=0; i<numDims; i++ )
        {
          gridAcceleration[g][fieldIndex][i] = 0.0;
          gridVelocity[g][fieldIndex][i] = 0.0;
          gridMomentum[g][fieldIndex][i] = 0.0;
          gridMaterialPosition[g][fieldIndex][i] = 0 * gridPosition[g][i]; // TODO: zero? since it's supposed to be relative position?
        }
      }
    } );
  }
}

void SolidMechanicsMPM::enforceContact( real64 dt,
                                        DomainPartition & domain,
                                        ParticleManager & particleManager,
                                        NodeManager & nodeManager,
                                        MeshLevel & mesh )
{
  // Grid fields
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const gridPosition = nodeManager.referencePosition();
  arrayView2d< real64 const > const & gridMass = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::massString() );
  arrayView2d< real64 const > const & gridDamage = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::damageString() );
  arrayView2d< real64 const > const & gridMaxDamage = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::maxDamageString() );
  arrayView3d< real64 > const & gridVelocity = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::velocityString() );
  arrayView3d< real64 > const & gridMomentum = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::momentumString() );
  arrayView3d< real64 > const & gridAcceleration = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::accelerationString() );
  arrayView3d< real64 > const & gridContactForce = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::forceContactString() );
  arrayView3d< real64 > const & gridSurfaceNormal = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::surfaceNormalString() );
  arrayView3d< real64 > const & gridMaterialPosition = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::materialPositionString() );

  // Compute grid surface normals
  computeGridSurfaceNormals( particleManager, nodeManager );

  // Sync surface normals
  syncGridFields( { viewKeyStruct::surfaceNormalString() }, domain, nodeManager, mesh, MPI_SUM );

  // Apply symmetry boundary conditions to surface normals
  enforceGridVectorFieldSymmetryBC( gridSurfaceNormal, gridPosition, nodeManager.sets() );

  // Normalize grid surface normals
  normalizeGridSurfaceNormals( gridMass, gridSurfaceNormal );

  // Apply symmetry boundary conditions to material position
  enforceGridVectorFieldSymmetryBC( gridMaterialPosition, gridPosition, nodeManager.sets() );

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
  int const numNodes = nodeManager.size();
  real64 const smallMass = m_smallMass;
  int const numDims = m_numDims;
  for( int fieldIndex=0; fieldIndex<m_numVelocityFields; fieldIndex++ )
  {
    forAll< serialPolicy >( numNodes, [=] GEOS_HOST_DEVICE ( localIndex const g )
    {
      if( gridMass[g][fieldIndex] > smallMass ) // small mass threshold
      {
        for( int i=0; i<numDims; i++ )
        {
          real64 contactForce = gridContactForce[g][fieldIndex][i];
          gridAcceleration[g][fieldIndex][i] += contactForce / gridMass[g][fieldIndex];
          gridMomentum[g][fieldIndex][i] += contactForce * dt;
          gridVelocity[g][fieldIndex][i] = gridMomentum[g][fieldIndex][i] / gridMass[g][fieldIndex];
        }
      }
    } );
  }
}

void SolidMechanicsMPM::interpolateFTable( real64 dt, real64 time_n )
{
  double Fii_new;
  double Fii_dot;
  int fInterval = 0;
  double timeInterval;
  double timePast;

  if( time_n + dt < m_fTable[m_fTable.size( 0 ) - 1][0] ) // If within F table bounds...
  {
    for( localIndex i = 0; i < m_fTable.size( 0 ); i++ ) // Naive method for determining what part of F table we're currently in, can
                                                         // optimize later (TODO)
    {
      if( time_n + dt / 2 > m_fTable[i][0] )
      {
        fInterval = i;
      }
    }

    timeInterval = m_fTable[fInterval + 1][0] - m_fTable[fInterval][0]; // Time fInterval for current part of F table we're in
    timePast = time_n + dt - m_fTable[fInterval][0]; // Time elapsed since switching intervals in F table

    for( int i = 0; i < m_numDims; i++ )   // Update L and F
    {
      if ( m_stressControl[i] != 1)
      {
        // smooth-step interpolation with cosine, zero endpoint velocity
        if( m_fTableInterpType == 1 )
        {
          Fii_new = m_fTable[fInterval][i + 1] - 0.5 * ( m_fTable[fInterval + 1][i + 1] - m_fTable[fInterval][i + 1] ) * ( cos( 3.141592653589793 * timePast / timeInterval ) - 1.0 );
        }
        // smooth-step interpolation with 5th order polynomial, zero endpoint velocity and acceleration
        else if( m_fTableInterpType == 2 )
        {
          Fii_new = m_fTable[fInterval][i+1] + (m_fTable[fInterval+1][i+1] - m_fTable[fInterval][i+1])*
                    (10.0*pow( timePast/timeInterval, 3 ) - 15.0*pow( timePast/timeInterval, 4 ) + 6.0*pow( timePast/timeInterval, 5 ));
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
  }
  else if( time_n + dt >= m_fTable[m_fTable.size( 0 ) - 1][0] ) // Else (i.e. if we exceed F table upper bound)...
  {
    for( int i = 0; i < m_numDims; i++ )
    {
      if( m_stressControl[i] != 1 )
      {
      Fii_new = m_fTable[m_fTable.size( 0 ) - 1][i + 1]; // Set Fii_new to the last prescribed F table value
      Fii_dot = ( Fii_new - m_domainF[i] ) / dt;
      m_domainL[i] = Fii_dot / Fii_new; // L = Fdot.Finv
      m_domainF[i] = Fii_new;
      }
    }
  }
}

void SolidMechanicsMPM::interpolateStressTable( real64 dt,
                                                real64 time_n )
{
  double stressii_new;
  int stressInterval = 0;
  double timeInterval;
  double timePast;

  if( time_n + dt < m_stressTable[m_stressTable.size( 0 ) - 1][0] ) // If within F table bounds...
  {
    for( localIndex i = 0; i < m_stressTable.size( 0 ); i++ ) // Naive method for determining what part of F table we're currently in, can
                                                         // optimize later (TODO)
    {
      if( time_n + dt / 2 > m_stressTable[i][0] )
      {
        stressInterval = i;
      }
    }

    timeInterval = m_stressTable[stressInterval + 1][0] - m_stressTable[stressInterval][0]; // Time fInterval for current part of F table we're in
    timePast = time_n + dt - m_stressTable[stressInterval][0]; // Time elapsed since switching intervals in F table

    for( int i = 0; i < m_numDims; i++ )   // Update L and F
    {
      // smooth-step interpolation with cosine, zero endpoint velocity
      if( m_stressTableInterpType == 1 )
      {
        stressii_new = m_stressTable[stressInterval][i + 1] - 0.5 * ( m_stressTable[stressInterval + 1][i + 1] - m_stressTable[stressInterval][i + 1] ) * ( cos( 3.141592653589793 * timePast / timeInterval ) - 1.0 );
      }
      // smooth-step interpolation with 5th order polynomial, zero endpoint velocity and acceleration
      else if( m_stressTableInterpType == 2 )
      {
        stressii_new = m_stressTable[stressInterval][i+1] + (m_stressTable[stressInterval+1][i+1] - m_stressTable[stressInterval][i+1])*
                  (10.0*pow( timePast/timeInterval, 3 ) - 15.0*pow( timePast/timeInterval, 4 ) + 6.0*pow( timePast/timeInterval, 5 ));
      }
      // default linear interpolation
      else
      {
        stressii_new = m_stressTable[stressInterval][i + 1] * ( timeInterval - timePast ) / timeInterval + m_stressTable[stressInterval + 1][i + 1] * ( timePast ) / timeInterval;
      }

      m_domainStress[i] = stressii_new;
    }
  }
  else if( time_n + dt >= m_stressTable[m_stressTable.size( 0 ) - 1][0] ) // Else (i.e. if we exceed F table upper bound)...
  {
    for( int i = 0; i < m_numDims; i++ )
    {
      stressii_new = m_stressTable[m_stressTable.size( 0 ) - 1][i + 1]; // Set Fii_new to the last prescribed F table value
      m_domainStress[i] = stressii_new;
    }
  }
}

void SolidMechanicsMPM::gridToParticle( real64 dt,
                                        ParticleManager & particleManager,
                                        NodeManager & nodeManager )
{
      // Grid-to-particle map
      switch(( SolidMechanicsMPM::UpdateMethodOption ) m_updateMethod)
      {
        case SolidMechanicsMPM::UpdateMethodOption::FLIP:
          performFLIPUpdate( dt,
                             particleManager,
                             nodeManager );
          break;
        case SolidMechanicsMPM::UpdateMethodOption::PIC:
          performPICUpdate( dt, 
                            particleManager,
                            nodeManager );
          break;
        case SolidMechanicsMPM::UpdateMethodOption::XPIC:
          performXPICUpdate( dt,
                             particleManager,
                             nodeManager );
          break;
        case SolidMechanicsMPM::UpdateMethodOption::FMPM:
          performFMPMUpdate( dt, 
                             particleManager, 
                             nodeManager );
          break;
        default:
          GEOS_ERROR("SolidMechanicsMPM solver update method not recognized.");
          break;
      }
}

void SolidMechanicsMPM::performFLIPUpdate( real64 dt,
                                           ParticleManager & particleManager,
                                           NodeManager & nodeManager ){
  // Grid fields
  arrayView2d< real64 const > const & gridDamageGradient = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::damageGradientString() );
  arrayView3d< real64 const > const & gridVelocity = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::velocityString() );
  arrayView3d< real64 const > const & gridAcceleration = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::accelerationString() );

  localIndex subRegionIndex = 0;
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {


    // Registered by subregion
    arrayView2d< real64 > const particlePosition = subRegion.getParticleCenter();
    arrayView2d< real64 > const particleVelocity = subRegion.getParticleVelocity();
    arrayView1d< int const > const particleGroup = subRegion.getParticleGroup();

    // Registered by MPM solver
    arrayView3d< real64 > const particleVelocityGradient = subRegion.getField< fields::mpm::particleVelocityGradient >();
    arrayView2d< real64 const > const particleDamageGradient = subRegion.getField< fields::mpm::particleDamageGradient >();

    // Get views to mapping arrays
    int const numberOfVerticesPerParticle = subRegion.numberOfVerticesPerParticle();
    arrayView2d< localIndex const > const mappedNodes = m_mappedNodes[subRegionIndex];
    arrayView2d< real64 const > const shapeFunctionValues = m_shapeFunctionValues[subRegionIndex];
    arrayView3d< real64 const > const shapeFunctionGradientValues = m_shapeFunctionGradientValues[subRegionIndex];

    // Map to particles
    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
    int const numDims = m_numDims;
    int const damageFieldPartitioning = m_damageFieldPartitioning;
    int const numContactGroups = m_numContactGroups;
    forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
    {
      localIndex const p = activeParticleIndices[pp];

      // Initialize velocity gradient for this particle
      // if( m_directionalOverlapCorrection == 0 )
      // {
      for( int i=0; i<3; i++ )
      {
        for( int j=0; j<3; j++ )
        {
          particleVelocityGradient[p][i][j] = 0.0;
        }
      }
      // }
    
      for( int g = 0; g < 8 * numberOfVerticesPerParticle; g++ )
      {
        localIndex const mappedNode = mappedNodes[pp][g];
        int const nodeFlag = ( damageFieldPartitioning == 1 && LvArray::tensorOps::AiBi< 3 >( gridDamageGradient[mappedNode], particleDamageGradient[p] ) < 0.0 ) ? 1 : 0; // 0
                                                                                                                                                                          // undamaged
                                                                                                                                                                          // or
                                                                                                                                                                          // "A"
                                                                                                                                                                          // field,
                                                                                                                                                                          // 1
                                                                                                                                                                          // for
                                                                                                                                                                          // "B"
                                                                                                                                                                          // field
        int const fieldIndex = nodeFlag * numContactGroups + particleGroup[p]; // This ranges from 0 to nMatFields-1
        for( int i=0; i<numDims; i++ )
        {
          // Full step update
          // particlePosition[p][i] += gridVelocity[mappedNode][fieldIndex][i] * shapeFunctionValues[pp][g] * dt;
          
          // Half step update ( equivalent to old geos, but if there is only one cell in the simulation that is anomalous issues with the particle velocity gradient and by extension the stresses)
          // If the simulation size is increased to 2x2x2 then the issue goes away and likewise if the rvectors are scaled to 80% for a single cell simulation
          particlePosition[p][i] += ( gridVelocity[mappedNode][fieldIndex][i] - 0.5 * gridAcceleration[mappedNode][fieldIndex][i] * dt ) * shapeFunctionValues[pp][g] * dt;
          
          particleVelocity[p][i] += gridAcceleration[mappedNode][fieldIndex][i] * dt * shapeFunctionValues[pp][g]; // FLIP
          for( int j=0; j<numDims; j++ )
          {
            particleVelocityGradient[p][i][j] += gridVelocity[mappedNode][fieldIndex][i] * shapeFunctionGradientValues[pp][g][j]; // Technically
                                                                                                                                  // wrong,
                                                                                                                                  // the
                                                                                                                                  // best
                                                                                                                                  // kind of
                                                                                                                                  // wrong
                                                                                                                                  // (end-of-step
                                                                                                                                  // velocity
                                                                                                                                  // with
                                                                                                                                  // beginning-of-step
                                                                                                                                  // gradient)
          }
        }
      }
    } );
    subRegionIndex++;
  } );
}

void SolidMechanicsMPM::performPICUpdate(  real64 dt,
                                           ParticleManager & particleManager,
                                           NodeManager & nodeManager )
{
  // Grid fields
  arrayView2d< real64 const > const & gridDamageGradient = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::damageGradientString() );
  arrayView3d< real64 const > const & gridVelocity = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::velocityString() );
  arrayView3d< real64 const > const & gridAcceleration = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::accelerationString() );

  localIndex subRegionIndex = 0;
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    // Registered by subregion
    arrayView2d< real64 > const particlePosition = subRegion.getParticleCenter();
    arrayView2d< real64 > const particleVelocity = subRegion.getParticleVelocity();
    arrayView1d< int const > const particleGroup = subRegion.getParticleGroup();

    // Registered by MPM solver
    arrayView3d< real64 > const particleVelocityGradient = subRegion.getField< fields::mpm::particleVelocityGradient >();
    arrayView2d< real64 const > const particleDamageGradient = subRegion.getField< fields::mpm::particleDamageGradient >();

    // Get views to mapping arrays
    int const numberOfVerticesPerParticle = subRegion.numberOfVerticesPerParticle();
    arrayView2d< localIndex const > const mappedNodes = m_mappedNodes[subRegionIndex];
    arrayView2d< real64 const > const shapeFunctionValues = m_shapeFunctionValues[subRegionIndex];
    arrayView3d< real64 const > const shapeFunctionGradientValues = m_shapeFunctionGradientValues[subRegionIndex];

    // Map to particles
    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
    int const numDims = m_numDims;
    int const damageFieldPartitioning = m_damageFieldPartitioning;
    int const numContactGroups = m_numContactGroups;
    forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
    {
      localIndex const p = activeParticleIndices[pp];

      // Initialize velocity gradient for this particle
      // if( m_directionalOverlapCorrection == 0 )
      // {
      for( int i=0; i<3; i++ )
      {
        for( int j=0; j<3; j++ )
        {
          particleVelocityGradient[p][i][j] = 0.0;
        }
      }
      // }
    
      for( int g = 0; g < 8 * numberOfVerticesPerParticle; g++ )
      {
        localIndex const mappedNode = mappedNodes[pp][g];
        int const nodeFlag = ( damageFieldPartitioning == 1 && LvArray::tensorOps::AiBi< 3 >( gridDamageGradient[mappedNode], particleDamageGradient[p] ) < 0.0 ) ? 1 : 0; // 0
                                                                                                                                                                          // undamaged
                                                                                                                                                                          // or
                                                                                                                                                                          // "A"
                                                                                                                                                                          // field,
                                                                                                                                                                          // 1
                                                                                                                                                                          // for
                                                                                                                                                                          // "B"
                                                                                                                                                                          // field
        int const fieldIndex = nodeFlag * numContactGroups + particleGroup[p]; // This ranges from 0 to nMatFields-1
        for( int i=0; i<numDims; i++ )
        {
          particlePosition[p][i] += ( gridVelocity[mappedNode][fieldIndex][i] - 0.5 * gridAcceleration[mappedNode][fieldIndex][i] * dt ) * shapeFunctionValues[pp][g] * dt; // CC: position update doesn't seem consistent with old GEOS for FLIP and PIC, need to double check
          particleVelocity[p][i] += gridVelocity[mappedNode][fieldIndex][i] * shapeFunctionValues[pp][g];
          for( int j=0; j<numDims; j++ )
          {
            // CC: what should the velocity gradient update be for PIC?
            particleVelocityGradient[p][i][j] += gridVelocity[mappedNode][fieldIndex][i] * shapeFunctionGradientValues[pp][g][j]; // Technically
                                                                                                                                  // wrong,
                                                                                                                                  // the
                                                                                                                                  // best
                                                                                                                                  // kind of
                                                                                                                                  // wrong
                                                                                                                                  // (end-of-step
                                                                                                                                  // velocity
                                                                                                                                  // with
                                                                                                                                  // beginning-of-step
                                                                                                                                  // gradient)
          }
        }
      }
    } );
    subRegionIndex++;
  } );
}

void SolidMechanicsMPM::performXPICUpdate( real64 dt,
                                           ParticleManager & particleManager,
                                           NodeManager & nodeManager )
{
  // Grid fields
  arrayView2d< real64 const > const & gridMass = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::massString() );
  arrayView2d< real64 const > const & gridDamageGradient = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::damageGradientString() );
  arrayView3d< real64 const > const & gridVelocity = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::velocityString() );
  arrayView3d< real64 const > const & gridAcceleration = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::accelerationString() );

  int numNodes = nodeManager.size();

  localIndex subRegionIndex = 0;
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    // Registered by subregion
    arrayView2d< real64 > const particlePosition = subRegion.getParticleCenter();
    arrayView2d< real64 > const particleVelocity = subRegion.getParticleVelocity();
    arrayView1d< int const > const particleGroup = subRegion.getParticleGroup();

    // Registered by MPM solver
    arrayView1d< real64 > const particleMass = subRegion.getField< fields::mpm::particleMass >();
    arrayView3d< real64 > const particleVelocityGradient = subRegion.getField< fields::mpm::particleVelocityGradient >();
    arrayView2d< real64 const > const particleDamageGradient = subRegion.getField< fields::mpm::particleDamageGradient >();

    // Get views to mapping arrays
    int const numberOfVerticesPerParticle = subRegion.numberOfVerticesPerParticle();
    arrayView2d< localIndex const > const mappedNodes = m_mappedNodes[subRegionIndex];
    arrayView2d< real64 const > const shapeFunctionValues = m_shapeFunctionValues[subRegionIndex];
    // arrayView3d< real64 const > const shapeFunctionGradientValues = m_shapeFunctionGradientValues[subRegionIndex];

    // Map to particles
    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
    int const numDims = m_numDims;
    int const damageFieldPartitioning = m_damageFieldPartitioning;
    int const numContactGroups = m_numContactGroups;
    int const updateOrder = m_updateOrder;

    // For iterative XPIC solve
    array3d< real64 > vStar;
    array3d< real64 > vPlus;
    array3d< real64 > vMinus;

    vStar.resize( numNodes, numContactGroups, numDims );
    vPlus.resize( numNodes, numContactGroups, numDims );
    vMinus.resize( numNodes, numContactGroups, numDims );

    // Zero out vStar for each order iteration
    for( int n=0; n < numNodes; n++)
    {
      for( int cg=0; cg < numContactGroups; cg++)
      {
        for( int i = 0; i < numDims; i++)
        {
          vStar[n][cg][i] = 0.0;
          vMinus[n][cg][i] = gridVelocity[n][cg][i] - gridAcceleration[n][cg][i] * dt;
        }
      }
    }

    // Do XPIC iterations
    for(int r=0; r < updateOrder; ++r)
    {
      // Zero out vPlus for each order iteration
      for( int n=0; n < numNodes; n++)
      {
        for( int cg=0; cg < numContactGroups; cg++)
        {
          for( int i = 0; i < numDims; i++)
          {
            vPlus[n][cg][i] = 0.0;
          }
        }
      }

      forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
      {
        localIndex const p = activeParticleIndices[pp];
      
        for( int gi = 0; gi < 8 * numberOfVerticesPerParticle; gi++ )
        {
          localIndex const mappedNodeI = mappedNodes[pp][gi];
          int const nodeFlagI = ( damageFieldPartitioning == 1 && LvArray::tensorOps::AiBi< 3 >( gridDamageGradient[mappedNodeI], particleDamageGradient[p] ) < 0.0 ) ? 1 : 0;
          int const fieldIndexI = nodeFlagI * numContactGroups + particleGroup[p]; // This ranges from 0 to nMatFields-1
          for(int gj = 0; gj < 8 * numberOfVerticesPerParticle; gj++ )
          {
            localIndex const mappedNodeJ = mappedNodes[pp][gj];
            int const nodeFlagJ = ( damageFieldPartitioning == 1 && LvArray::tensorOps::AiBi< 3 >( gridDamageGradient[mappedNodeJ], particleDamageGradient[p] ) < 0.0 ) ? 1 : 0;
            int const fieldIndexJ = nodeFlagJ * numContactGroups + particleGroup[p]; // This ranges from 0 to nMatFields-1

            for (int i = 0; i < numDims; i++){
              vPlus[mappedNodeI][fieldIndexI][gi] += ( ( updateOrder - r + 1 ) / r ) * 
                                                       ( gridMass[mappedNodeI][fieldIndexI] * shapeFunctionValues[pp][gi] * shapeFunctionValues[pp][gj] / particleMass[p] ) * 
                                                         vMinus[mappedNodeJ][fieldIndexJ][gj];
            }

          }
        }
      } );

      // Update vStar
      for( int n=0; n < numNodes; n++)
      {
        for( int cg=0; cg < numContactGroups; cg++)
        {
          for( int i = 0; i < numDims; i++)
          {
            vStar[n][cg][i] += std::pow(-1, r) * vPlus[n][cg][i];
            vMinus[n][cg][i] = vPlus[n][cg][i];
          }
        }
      }
    } //End of updateOrder iterations

    // Update particles position and velocities now
    forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
    {
      localIndex const p = activeParticleIndices[pp];
    
      for( int g = 0; g < 8 * numberOfVerticesPerParticle; g++ )
      {
        localIndex const mappedNode = mappedNodes[pp][g];
        int const nodeFlag = ( damageFieldPartitioning == 1 && LvArray::tensorOps::AiBi< 3 >( gridDamageGradient[mappedNode], particleDamageGradient[p] ) < 0.0 ) ? 1 : 0;
        int const fieldIndex = nodeFlag * numContactGroups + particleGroup[p]; // This ranges from 0 to nMatFields-1
        for( int i=0; i<numDims; i++ )
        {
          real64 m = updateOrder;
          real64 S = shapeFunctionValues[pp][g];
          real64 gVPlus = gridVelocity[mappedNode][fieldIndex][i];
          real64 gA = gridAcceleration[mappedNode][fieldIndex][i];
          real64 pV = particleVelocity[p][i];

          particlePosition[p][i] += S * gVPlus * dt - ( S * gA + ( -m * S * ( gVPlus - gA * dt ) + pV + m * S * vStar[mappedNode][fieldIndex][i]) ) * dt * dt / 2;
          particleVelocity[p][i] += S * ( m * ( gVPlus - vStar[mappedNode][fieldIndex][i] ) + ( 1 - m ) * gA * dt );
          // CC: What about update to velocity gradient?
        }
      }
    } );
    subRegionIndex++;
  } );
}

void SolidMechanicsMPM::performFMPMUpdate(  real64 dt,
                                            ParticleManager & particleManager,
                                            NodeManager & nodeManager )
{
  // Grid fields
  arrayView2d< real64 const > const & gridMass = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::massString() );
  arrayView2d< real64 const > const & gridDamageGradient = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::damageGradientString() );
  arrayView3d< real64 const > const & gridVelocity = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::velocityString() );

  int numNodes = nodeManager.size();

  localIndex subRegionIndex = 0;
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    // Registered by subregion
    arrayView2d< real64 > const particlePosition = subRegion.getParticleCenter();
    arrayView2d< real64 > const particleVelocity = subRegion.getParticleVelocity();
    arrayView1d< int const > const particleGroup = subRegion.getParticleGroup();

    // Registered by MPM solver
    arrayView1d< real64 > const particleMass = subRegion.getField< fields::mpm::particleMass >();
    arrayView3d< real64 > const particleVelocityGradient = subRegion.getField< fields::mpm::particleVelocityGradient >();
    arrayView2d< real64 const > const particleDamageGradient = subRegion.getField< fields::mpm::particleDamageGradient >();

    // Get views to mapping arrays
    int const numberOfVerticesPerParticle = subRegion.numberOfVerticesPerParticle();
    arrayView2d< localIndex const > const mappedNodes = m_mappedNodes[subRegionIndex];
    arrayView2d< real64 const > const shapeFunctionValues = m_shapeFunctionValues[subRegionIndex];
    // arrayView3d< real64 const > const shapeFunctionGradientValues = m_shapeFunctionGradientValues[subRegionIndex];

    // Map to particles
    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
    int const numDims = m_numDims;
    int const damageFieldPartitioning = m_damageFieldPartitioning;
    int const numContactGroups = m_numContactGroups;
    int const updateOrder = m_updateOrder;

    // For iterative FMPM solve
    array3d< real64 > vStar;
    array3d< real64 > vPlus;
    array3d< real64 > vMinus;

    vStar.resize( numNodes, numContactGroups, numDims );
    vPlus.resize( numNodes, numContactGroups, numDims );
    vMinus.resize( numNodes, numContactGroups, numDims );

    // Zero out vStar for each order iteration
    for( int n=0; n < numNodes; n++)
    {
      for( int cg=0; cg < numContactGroups; cg++)
      {
        for( int i = 0; i < numDims; i++)
        {
          vStar[n][cg][i] = 0.0;
          vMinus[n][cg][i] = gridVelocity[n][cg][i];
        }
      }
    }

    // Do XPIC iterations
    for(int r=0; r < updateOrder; ++r)
    {
      // Zero out vPlus for each order iteration
      for( int n=0; n < numNodes; n++)
      {
        for( int cg=0; cg < numContactGroups; cg++)
        {
          for( int i = 0; i < numDims; i++)
          {
            vPlus[n][cg][i] = 0.0;
          }
        }
      }

      forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
      {
        localIndex const p = activeParticleIndices[pp];
      
        for( int gi = 0; gi < 8 * numberOfVerticesPerParticle; gi++ )
        {
          localIndex const mappedNodeI = mappedNodes[pp][gi];
          int const nodeFlagI = ( damageFieldPartitioning == 1 && LvArray::tensorOps::AiBi< 3 >( gridDamageGradient[mappedNodeI], particleDamageGradient[p] ) < 0.0 ) ? 1 : 0;
          int const fieldIndexI = nodeFlagI * numContactGroups + particleGroup[p]; // This ranges from 0 to nMatFields-1
          for(int gj = 0; gj < 8 * numberOfVerticesPerParticle; gj++ )
          {
            localIndex const mappedNodeJ = mappedNodes[pp][gj];
            int const nodeFlagJ = ( damageFieldPartitioning == 1 && LvArray::tensorOps::AiBi< 3 >( gridDamageGradient[mappedNodeJ], particleDamageGradient[p] ) < 0.0 ) ? 1 : 0;
            int const fieldIndexJ = nodeFlagJ * numContactGroups + particleGroup[p]; // This ranges from 0 to nMatFields-1

            for (int i = 0; i < numDims; i++){
              vPlus[mappedNodeI][fieldIndexI][gi] += ( ( updateOrder - r + 1 ) / r ) * 
                                                       ( gridMass[mappedNodeI][fieldIndexI] * shapeFunctionValues[pp][gi] * shapeFunctionValues[pp][gj] / particleMass[p] ) * 
                                                         vMinus[mappedNodeJ][fieldIndexJ][gj];
            }

          }
        }
      } );

      // Update vStar
      for( int n=0; n < numNodes; n++)
      {
        for( int cg=0; cg < numContactGroups; cg++)
        {
          for( int i = 0; i < numDims; i++)
          {
            vStar[n][cg][i] += std::pow(-1, r) * vPlus[n][cg][i];
            vMinus[n][cg][i] = vPlus[n][cg][i];
          }
        }
      }
    } //End of updateOrder iterations

    // Update particles position and velocities now
    forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
    {
      localIndex const p = activeParticleIndices[pp];
    
      for( int g = 0; g < 8 * numberOfVerticesPerParticle; g++ )
      {
        localIndex const mappedNode = mappedNodes[pp][g];
        int const nodeFlag = ( damageFieldPartitioning == 1 && LvArray::tensorOps::AiBi< 3 >( gridDamageGradient[mappedNode], particleDamageGradient[p] ) < 0.0 ) ? 1 : 0;
        int const fieldIndex = nodeFlag * numContactGroups + particleGroup[p]; // This ranges from 0 to nMatFields-1
        for( int i=0; i<numDims; i++ )
        {
          particlePosition[p][i] += 0.5 * dt * ( particleVelocity[p][g]  + shapeFunctionValues[pp][g] * vStar[mappedNode][fieldIndex][i] );
          particleVelocity[p][i] += shapeFunctionValues[pp][g] * vStar[mappedNode][fieldIndex][i];
          // CC: What about update to velocity gradient?
        }
      }
    } );
    subRegionIndex++;
  } );
}

void SolidMechanicsMPM::updateSolverDependencies( ParticleManager & particleManager )
{
  // CC: should plastic strain updated here?
  // CC: TODO add update for material direction for graphite model
  // Get particle damage values from constitutive model
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    // Get constitutive model reference
    string const & solidMaterialName = subRegion.template getReference< string >( viewKeyStruct::solidMaterialNamesString() );
    SolidBase & solidModel = getConstitutiveModel< SolidBase >( subRegion, solidMaterialName );
    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();

    if( solidModel.hasWrapper( "damage" ) ) // Fragile code because someone could change the damage key without our knowledge. TODO: Make an
                                            // integrated test that checks this
    {
      arrayView1d< real64 > const particleDamage = subRegion.getParticleDamage();
      arrayView2d< real64 const > const constitutiveDamage = solidModel.getReference< array2d< real64 > >( "damage" );
      forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
      {
        localIndex const p = activeParticleIndices[pp];
        if( constitutiveDamage[p][0] > particleDamage[p] ) // Damage can only increase - This will also preserve user-specified damage at
                                                           // initialization
        {
          particleDamage[p] = constitutiveDamage[p][0]; // TODO: Load any pre-damage into the constitutive model at initialization. Or,
                                                        // switch to using VTK input such that we can initialize any field we want without
                                                        // explicitly post-processing the input file.
        }
      } );
    }

    if( solidModel.hasWrapper( "plasticStrain" ) ) // Fragile code because someone could change the damage key without our knowledge. TODO: Make an
                                            // integrated test that checks this
    {
      arrayView2d< real64 > const particlePlasticStrain = subRegion.getField< fields::mpm::particlePlasticStrain >();
      arrayView3d< real64 const > const constitutivePlasticStrain = solidModel.getReference< array3d< real64 > >( "plasticStrain" );
      forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
      {
        localIndex const p = activeParticleIndices[pp];
        LvArray::tensorOps::copy< 6 >( particlePlasticStrain[p], constitutivePlasticStrain[p][0] );
      } );
    }
  } );
}

real64 SolidMechanicsMPM::getStableTimeStep( ParticleManager & particleManager )
{
  real64 wavespeed = 0.0;
  real64 length = m_planeStrain == 1 ? std::fmin( m_hEl[0], m_hEl[1] ) : std::fmin( m_hEl[0], std::fmin( m_hEl[1], m_hEl[2] ) );

  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    arrayView1d< real64 const> const particleDensity = subRegion.getField< fields::mpm::particleDensity >();
    arrayView2d< real64 const > const particleVelocity = subRegion.getParticleVelocity();

    string const & solidMaterialName = subRegion.template getReference< string >( viewKeyStruct::solidMaterialNamesString() );
    SolidBase & constitutiveRelation = getConstitutiveModel< SolidBase >( subRegion, solidMaterialName );

    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();

    // CC: there has to be a better way to do this, have all models calculate their own wavespeed?
    string constitutiveModelName = constitutiveRelation.getCatalogName();
    
    array1d< real64 > bulkModulus;
    array1d< real64 > shearModulus;
    if( constitutiveModelName == "Hyperelastic" ){
      const Hyperelastic & hyperelastic = dynamic_cast< const Hyperelastic & >( constitutiveRelation );
      bulkModulus = hyperelastic.bulkModulus();
      shearModulus = hyperelastic.shearModulus();
    }

    if( constitutiveModelName == "HyperelasticMMS" ){
      const HyperelasticMMS & hyperelasticMMS = dynamic_cast< const HyperelasticMMS & >( constitutiveRelation ); 
      arrayView1d< real64 const > const lambda = hyperelasticMMS.lambda();
      shearModulus = hyperelasticMMS.shearModulus();
      bulkModulus.resize(lambda.size());
      forAll< serialPolicy >( lambda.size(), [=, &bulkModulus] GEOS_HOST ( localIndex const p )
      {
        bulkModulus[p] = conversions::lameConstants::toBulkMod( lambda[p], shearModulus[p] );
      } );
    }

    if( constitutiveModelName == "ElasticIsotropic" || constitutiveModelName == "CeramicDamage" || constitutiveModelName == "StrainHardeningPolymer" ){
      const ElasticIsotropic & elasticIsotropic = dynamic_cast< const ElasticIsotropic & >( constitutiveRelation );
      bulkModulus = elasticIsotropic.bulkModulus();
      shearModulus = elasticIsotropic.shearModulus();
    }

    if( constitutiveModelName == "Graphite" || constitutiveModelName == "ElasticTransverseIsotropic" || constitutiveModelName == "ElasticTransverseIsotropicPressureDependent" ){
      const ElasticTransverseIsotropic & elasticTransverseIsotropic = dynamic_cast< const ElasticTransverseIsotropic & >( constitutiveRelation );
      bulkModulus = elasticTransverseIsotropic.effectiveBulkModulus();
      shearModulus = elasticTransverseIsotropic.effectiveShearModulus();
    }

    forAll< serialPolicy >( activeParticleIndices.size(), [=, &wavespeed] GEOS_HOST ( localIndex const pp ) // would need reduction to parallelize
    {
      localIndex const p = activeParticleIndices[pp];
      wavespeed = fmax( wavespeed, sqrt( ( bulkModulus[p] + (4.0/3.0) * shearModulus[p] ) / particleDensity[p] ) + LvArray::tensorOps::l2Norm< 3 >( particleVelocity[p] ) );
    } );
  } );

  real64 dtReturn = wavespeed > 1.0e-16 ? m_cflFactor * length / wavespeed : DBL_MAX; // This partitions's dt, make it huge if wavespeed=0.0
                                                                                      // (this happens when there are no particles on this
                                                                                      // partition)
  
  // CC: TODO add check to make sure that next time step won't skip a MPM event
  return dtReturn;
}

void SolidMechanicsMPM::deleteBadParticles( ParticleManager & particleManager )
{
  int size = 0;
  // Cases covered:
  // 1.) Particles that map outside the domain (including buffer cells)
  // 2.) Particles with unacceptable Jacobian (<0.1 or >10)
  // 3.) Particles that are removed by the machine sample event
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    // Move everything into the host memory space
    subRegion.forWrappers( [&]( WrapperBase & wrapper )
    {
      wrapper.move( LvArray::MemorySpace::host, true );
    } );

    // Get relevant particle arrays
    arrayView1d< int > const isBad = subRegion.getField< fields::mpm::isBad >();

    // Initialize the set of particles to delete
    std::set< localIndex > indicesToErase;
    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
    forAll< serialPolicy >( activeParticleIndices.size(), [=, &indicesToErase] GEOS_HOST ( localIndex const pp ) // need reduction or
                                                                                                                 // atomics to parallelize,
                                                                                                                 // not sure which
      {
        localIndex const p = activeParticleIndices[pp];
        if( isBad[p] == 1 )
        {
          indicesToErase.insert( p );
        }
      } );
    subRegion.erase( indicesToErase );
    subRegion.setActiveParticleIndices(); // CC: debug
    size = subRegion.activeParticleIndices().size();
  } );

  // particleManager.getRegion< ParticleRegion >( "ParticleRegion1" ).resize( size ) ;
}

void SolidMechanicsMPM::printProfilingResults()
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
  for( unsigned int i = 0; i < numIntervals; i++ )
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
    for( unsigned int i = 0; i < numIntervals; i++ )
    {
      std::cout << " (" << i << ") ";
      std::cout << std::fixed;
      std::cout << std::showpoint;
      std::cout << std::setprecision( 6 ) << timeIntervalsAllRanks[i] / totalStepTimeAllRanks;
      std::cout << ": " << m_profilingLabels[i] << std::endl;
    }
    std::cout << " ** Total step CPU time:  " << totalStepTimeAllRanks << " s **" << std::endl;
    std::cout << "---------------------------------------------" << std::endl;
  }

  // Reset profiling arrays
  m_profilingTimes.resize( 0 );
  m_profilingLabels.resize( 0 );
}

void SolidMechanicsMPM::computeSurfaceFlags( ParticleManager & particleManager )
{
  // Get accessors for volume, position
  ParticleManager::ParticleViewAccessor< arrayView1d< real64 const > > particleVolumeAccessor = particleManager.constructArrayViewAccessor< real64, 1 >( "particleVolume" );
  ParticleManager::ParticleViewAccessor< arrayView2d< real64 const > > particlePositionAccessor = particleManager.constructArrayViewAccessor< real64, 2 >( "particleCenter" );

  // Perform neighbor operations
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    // Get neighbor list
    OrderedVariableToManyParticleRelation & neighborList = subRegion.neighborList();
    arrayView1d< localIndex const > const numNeighborsAll = neighborList.m_numParticles.toViewConst();
    ArrayOfArraysView< localIndex const > const neighborRegions = neighborList.m_toParticleRegion.toViewConst();
    ArrayOfArraysView< localIndex const > const neighborSubRegions = neighborList.m_toParticleSubRegion.toViewConst();
    ArrayOfArraysView< localIndex const > const neighborIndices = neighborList.m_toParticleIndex.toViewConst();

    // Get particle position and surface flags
    arrayView2d< real64 const > const particlePosition = subRegion.getParticleCenter();
    arrayView1d< int > const particleSurfaceFlag = subRegion.getField< fields::mpm::particleSurfaceFlag >();

    // Loop over neighbors
    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
    forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST ( localIndex const pp ) // Must be on host since we call a 'this'
                                                                                                // method which uses class variables
      {
        localIndex const p = activeParticleIndices[pp];

        // Get number of neighbors and accessor indices
        localIndex numNeighbors = numNeighborsAll[p];
        arraySlice1d< localIndex const > const regionIndices = neighborRegions[p];
        arraySlice1d< localIndex const > const subRegionIndices = neighborSubRegions[p];
        arraySlice1d< localIndex const > const particleIndices = neighborIndices[p];

        // Calculate
        real64 totalVolume = 0.0;
        real64 centroid[3] = { 0.0 };
        for( localIndex neighborIndex = 0; neighborIndex < numNeighbors; neighborIndex++ )
        {
          localIndex regionIndex = regionIndices[neighborIndex];
          localIndex subRegionIndex = subRegionIndices[neighborIndex];
          localIndex particleIndex = particleIndices[neighborIndex];
          real64 rSquared = 0.0;
          for( int i=0; i<m_numDims; i++ )
          {
            real64 relativePositionComponent = particlePositionAccessor[regionIndex][subRegionIndex][particleIndex][i] - particlePosition[p][i];
            rSquared += relativePositionComponent * relativePositionComponent;
          }
          real64 kernelVal = kernel( sqrt( rSquared ));
          totalVolume += kernelVal * particleVolumeAccessor[regionIndex][subRegionIndex][particleIndex];
          for( int i=0; i<m_numDims; i++ )
          {
            centroid[i] += kernelVal * particleVolumeAccessor[regionIndex][subRegionIndex][particleIndex] * particlePositionAccessor[regionIndex][subRegionIndex][particleIndex][i];
          }
        }
        real64 positionDifference[3] = { 0.0 };
        for( int i=0; i<m_numDims; i++ )
        {
          centroid[i] /= totalVolume;
          positionDifference[i] = centroid[i] - particlePosition[p][i];
        }
        real64 deltaSquared = 0.0;
        for( int i=0; i<m_numDims; i++ )
        {
          deltaSquared += positionDifference[i] * positionDifference[i];
        }
        // 0.08 was hand-picked based on testing done in Mathematica. Works in 2D and 3D!
        if( particleSurfaceFlag[p] != 1 && deltaSquared >= 0.08 * 0.08 * m_neighborRadius * m_neighborRadius )
        {
          particleSurfaceFlag[p] = 1;
        }
      } );
  } );
}

void SolidMechanicsMPM::computeSphF( ParticleManager & particleManager )
{
  // Get accessors for volume, position, reference position
  ParticleManager::ParticleViewAccessor< arrayView1d< real64 const > > particleVolumeAccessor = particleManager.constructArrayViewAccessor< real64, 1 >( "particleVolume" );
  ParticleManager::ParticleViewAccessor< arrayView2d< real64 const > > particlePositionAccessor = particleManager.constructArrayViewAccessor< real64, 2 >( "particleCenter" );
  ParticleManager::ParticleViewAccessor< arrayView2d< real64 const > > particleReferencePositionAccessor = particleManager.constructArrayViewAccessor< real64, 2 >( "particleReferencePosition" );

  // Perform neighbor operations
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    // Get neighbor list
    OrderedVariableToManyParticleRelation & neighborList = subRegion.neighborList();
    arrayView1d< localIndex const > const numNeighborsAll = neighborList.m_numParticles.toViewConst();
    ArrayOfArraysView< localIndex const > const neighborRegions = neighborList.m_toParticleRegion.toViewConst();
    ArrayOfArraysView< localIndex const > const neighborSubRegions = neighborList.m_toParticleSubRegion.toViewConst();
    ArrayOfArraysView< localIndex const > const neighborIndices = neighborList.m_toParticleIndex.toViewConst();

    // Get particle position and sphF
    arrayView2d< real64 const > const particlePosition = subRegion.getParticleCenter();
    arrayView2d< real64 > const particleReferencePosition = subRegion.getField< fields::mpm::particleReferencePosition >();
    arrayView3d< real64 > const particleSphF = subRegion.getField< fields::mpm::particleSphF >();

    // Loop over neighbors
    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
    forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST ( localIndex const pp ) // I think this must be on host since we
                                                                                                // call a 'this' method which uses class
                                                                                                // variables
      {
        localIndex const p = activeParticleIndices[pp];

        // Get number of neighbors and accessor indices
        localIndex numNeighbors = numNeighborsAll[p];
        arraySlice1d< localIndex const > const regionIndices = neighborRegions[p];
        arraySlice1d< localIndex const > const subRegionIndices = neighborSubRegions[p];
        arraySlice1d< localIndex const > const particleIndices = neighborIndices[p];

        // Declare and size neighbor data arrays - TODO: switch to std::array? But then we'd need to template computeKernelFieldGradient
        std::vector< real64 > neighborVolumes( numNeighbors );
        std::vector< std::vector< real64 > > neighborPositions;
        neighborPositions.resize( numNeighbors, std::vector< real64 >( 3 ) );
        std::vector< std::vector< real64 > > neighborDisplacements;
        neighborDisplacements.resize( numNeighbors, std::vector< real64 >( 3 ) );

        // Populate neighbor data arrays
        for( localIndex neighborIndex = 0; neighborIndex < numNeighbors; neighborIndex++ )
        {
          localIndex regionIndex = regionIndices[neighborIndex];
          localIndex subRegionIndex = subRegionIndices[neighborIndex];
          localIndex particleIndex = particleIndices[neighborIndex];
          real64 r0Squared = 0.0;
          for( int i=0; i<m_numDims; i++ )
          {
            real64 thisComponent = particleReferencePosition[p][i];
            real64 neighborComponent = particleReferencePositionAccessor[regionIndex][subRegionIndex][particleIndex][i];
            r0Squared += ( thisComponent - neighborComponent ) * ( thisComponent - neighborComponent );
          }
          if( r0Squared <= m_neighborRadius * m_neighborRadius ) // We only consider neighbor particles that would have been neighbors in
                                                                 // the reference configuration
          {
            neighborVolumes[neighborIndex] = particleVolumeAccessor[regionIndex][subRegionIndex][particleIndex];
          }
          else
          {
            neighborVolumes[neighborIndex] = 0.0;
          }
          for( int i=0; i<3; i++ )
          {
            neighborPositions[neighborIndex][i] = particlePositionAccessor[regionIndex][subRegionIndex][particleIndex][i];
            neighborDisplacements[neighborIndex][i] = particlePositionAccessor[regionIndex][subRegionIndex][particleIndex][i] -
                                                      particleReferencePositionAccessor[regionIndex][subRegionIndex][particleIndex][i];
          }
        }

        // Call kernel field gradient function
        computeKernelVectorGradient( particlePosition[p],        // input
                                     neighborPositions,          // input
                                     neighborVolumes,            // input
                                     neighborDisplacements,      // input
                                     particleSphF[p] );          // OUTPUT: Spatial displacement gradient

        // Convert spatial displacement gradient to deformation gradient
        real64 temp[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
        LvArray::tensorOps::scaledAdd< 3, 3 >( temp, particleSphF[p], -1.0 );
        LvArray::tensorOps::invert< 3 >( temp );
        LvArray::tensorOps::copy< 3, 3 >( particleSphF[p], temp );
      } );
  } );
}

// void SolidMechanicsMPM::directionalOverlapCorrection( real64 dt, ParticleManager & particleManager )
// {
//   // Get the SPH version of the deformation gradient
//   computeSphF( particleManager );

//   // If we're at a surface, induce corrective deformation normal to the surface
//   particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
//   {
//     arrayView1d< real64 > const particleOverlap = subRegion.getField< fields::mpm::particleOverlap >();
//     particleOverlap.zero();
//     arrayView2d< real64 const > const particleDamageGradient = subRegion.getField< fields::mpm::particleDamageGradient >();
//     arrayView3d< real64 > const particleVelocityGradient = subRegion.getField< fields::mpm::particleVelocityGradient >();
//     particleVelocityGradient.zero();
//     arrayView3d< real64 const > const particleDeformationGradient = subRegion.getField< fields::mpm::particleDeformationGradient >();
//     arrayView3d< real64 const > const particleSphF = subRegion.getField< fields::mpm::particleSphF >();

//     SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
//     forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST ( localIndex const pp )
//     {
//       localIndex const p = activeParticleIndices[pp];
//       if( LvArray::tensorOps::l2NormSquared< 3 >( particleDamageGradient[p] ) > 0.0 )
//       {
//         real64 surfaceNormal[3];
//         LvArray::tensorOps::copy< 3 >( surfaceNormal, particleDamageGradient[p] );
//         if( m_planeStrain == 1 )
//         {
//           surfaceNormal[2] = 0.0; // Just in case
//         }
//         LvArray::tensorOps::normalize< 3 >( surfaceNormal );

//         // Polar decompositions
//         real64 Rp[3][3], Rsph[3][3], Vp[3][3], Vsph[3][3];
//         LvArray::tensorOps::polarDecomposition< 3 >( Rp, particleDeformationGradient[p] );
//         LvArray::tensorOps::polarDecomposition< 3 >( Rsph, particleSphF[p] );
//         LvArray::tensorOps::Rij_eq_AikBjk< 3, 3, 3 >(Vp, particleDeformationGradient[p], Rp );
//         LvArray::tensorOps::Rij_eq_AikBjk< 3, 3, 3 >(Vsph, particleSphF[p], Rsph );

//         // Compute directional stretches
//         real64 stretch[3], normalStretch;
//         real64 stretchSph[3], normalStretchSph;
//         LvArray::tensorOps::Ri_eq_AijBj< 3, 3 >(stretch, Vp, surfaceNormal );
//         normalStretch = LvArray::tensorOps::AiBi< 3 >( stretch, surfaceNormal );
//         LvArray::tensorOps::Ri_eq_AijBj< 3, 3 >(stretchSph, Vsph, surfaceNormal );
//         normalStretchSph = LvArray::tensorOps::AiBi< 3 >( stretchSph, surfaceNormal );

//         // Detect overlap and correct. TODO: Implement a ramp?
//         real64 overlap = normalStretch / normalStretchSph;
//         particleOverlap[p] = overlap;
//         if( overlap > 1.2 && normalStretchSph < 1 )
//         {
//           // Get orthonormal basis to transform V into
//           real64 s1[3], s2[3];
//           computeOrthonormalBasis( surfaceNormal, s1, s2 );
//           real64 Q[3][3], Qtranspose[3][3], targetF[3][3], targetFDot[3][3], FInverse[3][3];
//           for( int i=0; i<3; i++ )
//           {
//             Q[i][0] = surfaceNormal[i];
//             Q[i][1] = s1[i];
//             Q[i][2] = s2[i];
//             Qtranspose[0][i] = surfaceNormal[i];
//             Qtranspose[1][i] = s1[i];
//             Qtranspose[2][i] = s2[i];
//           }

//           // Perform transformation, substitute in the directional stretch predicted by Fsph, transform back, re-assemble F
//           real64 VpPrime[3][3], temp[3][3];
//           LvArray::tensorOps::Rij_eq_AikBkj< 3, 3, 3 >( temp, Qtranspose, Vp );
//           LvArray::tensorOps::Rij_eq_AikBkj< 3, 3, 3 >( VpPrime, temp, Q );
//           VpPrime[0][0] = normalStretchSph;
//           LvArray::tensorOps::Rij_eq_AikBkj< 3, 3, 3 >( temp, Q, VpPrime );
//           LvArray::tensorOps::Rij_eq_AikBkj< 3, 3, 3 >( Vp, temp, Qtranspose );
//           LvArray::tensorOps::Rij_eq_AikBkj< 3, 3, 3 >( targetF, Vp, Rp );

//           // Deduce the L necessary to achieve the new F
//           LvArray::tensorOps::copy< 3, 3 >( targetFDot, targetF );
//           LvArray::tensorOps::scaledAdd< 3, 3 >( targetFDot, particleDeformationGradient[p], -1.0 );
//           LvArray::tensorOps::scale< 3, 3 >( targetFDot, 1.0 / dt );
//           LvArray::tensorOps::invert< 3 >( FInverse, particleDeformationGradient[p] );
//           LvArray::tensorOps::Rij_eq_AikBkj< 3, 3, 3 >( particleVelocityGradient[p], targetFDot, FInverse );
//         }
//       }
//     } );
//   } );

//   // TODO: If we're not at a surface, volumetric overlap correction should be performed after the F update
// }

int SolidMechanicsMPM::evaluateSeparabilityCriterion( localIndex const & A,
                                                      localIndex const & B,
                                                      real64 const & damageA,
                                                      real64 const & damageB,
                                                      real64 const & maxDamageA,
                                                      real64 const & maxDamageB )
// m_treatFullyDamagedAsSingleField makes fields inseparable if damageA = damageB = 1, so we aren't putting
// arbitrary separation planes (and potential surfaces for accumulated overlap that needs corrections)
// between fully damaged materials. There is a potential issue that approaching damaged bodies
// would get contact gaps locked in.
{
  int separable = 0;
  // At least one field is fully damaged and both fields have the minimum separable level of damage.
  // The "a%b" is the "mod(a,b)" command, and indicates whether materials are from same contact group.
  if( ( ( maxDamageA >= 0.9999 || maxDamageB >= 0.9999 ) && ( damageA >= m_separabilityMinDamage && damageB >= m_separabilityMinDamage ) )
      || ( A % m_numContactGroups != B % m_numContactGroups ) )
  {
    if( m_treatFullyDamagedAsSingleField == 1 )
    {
      separable = ( fmin( damageA, damageB ) < 0.9999 ) ? 1 : 0;
    }
    else
    {
      separable = 1;
    }
  }

  return separable;
}

// CC: do I need to modify this to check for periodic boundaries
// All master particles should have centers inside the domain if particleCenters are corrected correctly during repartitioning
// CPDI: Edge case, corner of large or long particle beyond ghost cells but center is still inside domain?
void SolidMechanicsMPM::flagOutOfRangeParticles( ParticleManager & particleManager )
{
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    // Identify global domain bounds
    real64 globalMin[3]; // including buffer cells
    real64 globalMax[3]; // including buffer cells
    for( int i=0; i<3; i++ )
    {
      globalMin[i] = m_xGlobalMin[i] - m_hEl[i];
      globalMax[i] = m_xGlobalMax[i] + m_hEl[i];
    }

    // Get particle fields
    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
    arrayView2d< real64 const > const particlePosition = subRegion.getParticleCenter();
    arrayView1d< int > const isBad = subRegion.getField< fields::mpm::isBad >();
    ParticleType particleType = subRegion.getParticleType();

    // Define tolerance
    real64 tolerance[3];
    for( int i=0; i<3; i++ )
    {
      tolerance[i] = m_hEl[i] * DBL_EPSILON;
    }

    switch( particleType )
    {
      case ParticleType::SinglePoint:
        {
          forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
        {
          localIndex const p = activeParticleIndices[pp];
          for( int i=0; i<3; i++ )
          {
            if( particlePosition[p][i] < globalMin[i] + tolerance[i] || globalMax[i] - tolerance[i] < particlePosition[p][i] )
            {
              isBad[p] = 1;
              break; // TODO: if this doesn't work, just modify "i"
            }
          }
        } );
          break;
        }
      case ParticleType::CPDI:
        {
          arrayView3d< real64 const > const particleRVectors = subRegion.getParticleRVectors();
          int const signs[8][3] = { { 1, 1, 1},
            { 1, 1, -1},
            { 1, -1, 1},
            { 1, -1, -1},
            {-1, 1, 1},
            {-1, 1, -1},
            {-1, -1, 1},
            {-1, -1, -1} };
          forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
        {
          localIndex const p = activeParticleIndices[pp];
          for( int cornerIndex=0; cornerIndex<8; cornerIndex++ )
          {
            for( int i=0; i<3; i++ )
            {
              real64 cornerPositionComponent = particlePosition[p][i] + signs[cornerIndex][0] * particleRVectors[p][0][i] + signs[cornerIndex][1] * particleRVectors[p][1][i] + signs[cornerIndex][2] *
                                               particleRVectors[p][2][i];
              if( cornerPositionComponent < globalMin[i] + tolerance[i] || globalMax[i] - tolerance[i] < cornerPositionComponent )
              {
                isBad[p] = 1;
                break;
              }
            }
            if( isBad[p] == 1 )
            {
              break;
            }
          }
        } );
          break;
        }
      default:
        {
          GEOS_ERROR( "Particle type \"" << particleType << "\" is not yet supported." );
          break;
        }
    }
  } );
}

void SolidMechanicsMPM::computeRVectors( ParticleManager & particleManager )
{
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    if( subRegion.hasRVectors() )
    {
      SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
      arrayView3d< real64 > const particleRVectors = subRegion.getParticleRVectors();
      arrayView3d< real64 const > const particleInitialRVectors = subRegion.getField< fields::mpm::particleInitialRVectors >();
      arrayView3d< real64 const > const particleDeformationGradient = subRegion.getField< fields::mpm::particleDeformationGradient >();
      forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
      {
        localIndex const p = activeParticleIndices[pp];
        
        // for(int i = 0; i < 3; i++)
        // {
        //   LvArray::tensorOps::Ri_eq_AijBj< 3, 3 >( particleRVectors[p][i], particleDeformationGradient[p], particleRVectors[p][i] );
        // }
        for( int i=0; i<3; i++ )
        {
          for( int j=0; j<3; j++ )
          {
            particleRVectors[p][i][j] = particleInitialRVectors[p][i][0] * particleDeformationGradient[p][j][0] +
                                        particleInitialRVectors[p][i][1] * particleDeformationGradient[p][j][1] +
                                        particleInitialRVectors[p][i][2] * particleDeformationGradient[p][j][2];
          }
        }
      } );
    }
  } );
}

void SolidMechanicsMPM::cpdiDomainScaling( ParticleManager & particleManager )
{
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    if( subRegion.getParticleType() == ParticleType::CPDI )
    {
      real64 const lCrit = m_planeStrain == 1 ? 0.49999 * fmin( m_hEl[0], m_hEl[1] ) : 0.49999 * fmin( m_hEl[0], fmin( m_hEl[1], m_hEl[2] ) );
      arrayView3d< real64 > const particleRVectors = subRegion.getParticleRVectors();
      int const planeStrain = m_planeStrain;
      forAll< serialPolicy >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const p )
      {
        arraySlice1d< real64 > const r1 = particleRVectors[p][0];
        arraySlice1d< real64 > const r2 = particleRVectors[p][1];
        arraySlice1d< real64 > const r3 = particleRVectors[p][2];

        if( planeStrain == 1 ) // 2D cpdi domain scaling
        {
          // Initialize l-vectors.  Eq. 8a-d in the CPDI domain scaling paper.
          real64 l[2][3];
          for( int i=0; i<3; i++ )
          {
            l[0][i] = r1[i] + r2[i]; // la
            l[1][i] = r1[i] - r2[i]; // lb
          }

          // scale l-vectors if needed.  Eq. 9 in the CPDI domain scaling paper.
          bool scale = false;
          for( int i = 0; i < 2; i++ )
          {
            real64 lLength = sqrt( l[i][0] * l[i][0] + l[i][1] * l[i][1] + l[i][2] * l[i][2] );
            if( lLength > lCrit )
            {
              l[i][0] *= lCrit / lLength;
              l[i][1] *= lCrit / lLength;
              l[i][2] *= lCrit / lLength;
              scale = true;
            }
          }

          // reconstruct r-vectors.  eq. 11 in the CPDI domain scaling paper.
          if( scale )
          {
            for( int i=0; i<3; i++ )
            {
              r1[i] = 0.5 * (l[0][i] + l[1][i]);
              r2[i] = 0.5 * (l[0][i] - l[1][i]);
            }
          }
        }
        else // 3D cpdi domain scaling
        {
          // Initialize l-vectors.  Eq. 8a-d in the CPDI domain scaling paper.
          real64 l[4][3];
          for( int i=0; i<3; i++ )
          {
            l[0][i] = r1[i] + r2[i] + r3[i]; // la
            l[1][i] = r1[i] - r2[i] + r3[i]; // lb
            l[2][i] = r2[i] - r1[i] + r3[i]; // lc
            l[3][i] = r3[i] - r1[i] - r2[i]; // ld
          }

          // scale l vectors if needed.  Eq. 9 in the CPDI domain scaling paper.
          bool scale = false;
          for( int i = 0; i < 4; i++ )
          {
            real64 lLength = sqrt( l[i][0] * l[i][0] + l[i][1] * l[i][1] + l[i][2] * l[i][2] );
            if( lLength > lCrit )
            {
              l[i][0] *= lCrit / lLength;
              l[i][1] *= lCrit / lLength;
              l[i][2] *= lCrit / lLength;
              scale = true;
            }
          }

          // reconstruct r vectors.  eq. 11 in the CPDI domain scaling paper.
          if( scale )
          {
            for( int i=0; i<3; i++ )
            {
              r1[i] = 0.25 * ( l[0][i] + l[1][i] - l[2][i] - l[3][i] );
              r2[i] = 0.25 * ( l[0][i] - l[1][i] + l[2][i] - l[3][i] );
              r3[i] = 0.25 * ( l[0][i] + l[1][i] + l[2][i] + l[3][i] );
            }
          }
        }
      } );
    }
  } );
}

void SolidMechanicsMPM::resizeMappingArrays( ParticleManager & particleManager )
{
  // Count the number of subregions
  int numberOfSubRegions = 0;
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & GEOS_UNUSED_PARAM( subRegion ) ) // TODO: Check if this always accesses
                                                                                                   // subregions in the same order within a
                                                                                                   // single run
  {
    numberOfSubRegions++;
  } );

  m_mappedNodes.resize( numberOfSubRegions );
  m_shapeFunctionValues.resize( numberOfSubRegions );
  m_shapeFunctionGradientValues.resize( numberOfSubRegions );

  localIndex subRegionIndex = 0;
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    int const numberOfActiveParticles = subRegion.activeParticleIndices().size();
    m_mappedNodes[subRegionIndex].resize( numberOfActiveParticles, 8 * subRegion.numberOfVerticesPerParticle() );
    m_shapeFunctionValues[subRegionIndex].resize( numberOfActiveParticles, 8 * subRegion.numberOfVerticesPerParticle() );
    m_shapeFunctionGradientValues[subRegionIndex].resize( numberOfActiveParticles, 8 * subRegion.numberOfVerticesPerParticle(), 3 );
    subRegionIndex++;
  } );
}


void SolidMechanicsMPM::correctGhostParticleCentersAcrossPeriodicBoundaries(ParticleManager & particleManager,
                                                                         SpatialPartition & partition)
{
  arrayView1d< int const > const periodic = partition.getPeriodic();
  real64 xGlobalMin[3] = {0};
  real64 xGlobalMax[3] = {0};

  for( int i=0; i<3; i++ )
  {
    xGlobalMin[i] = m_xGlobalMin[i];
    xGlobalMax[i] = m_xGlobalMax[i];
  }

  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    
    SortedArrayView< localIndex const > const inactiveParticleIndices = subRegion.inactiveParticleIndices();
    arrayView2d< real64 > const particlePosition = subRegion.getParticleCenter();

    forAll< serialPolicy >( inactiveParticleIndices.size(), [&] GEOS_HOST ( localIndex const pp )
    {
      localIndex const p = inactiveParticleIndices[pp];

      for( int i=0 ; i<3 ; ++i)
      {
        if(periodic[i])
        {
          particlePosition[p][i] = Mod(particlePosition[p][i]-xGlobalMin[i], xGlobalMax[i]-xGlobalMin[i])+xGlobalMin[i];
        }
      }
    });
  });
}


void SolidMechanicsMPM::correctParticleCentersAcrossPeriodicBoundaries(ParticleManager & particleManager,
                                                                       SpatialPartition & partition)
{
  arrayView1d< int const > const periodic = partition.getPeriodic();
  real64 xGlobalMin[3] = {0};
  real64 xGlobalMax[3] = {0};

 for( int i=0; i<3; i++ )
 {
   xGlobalMin[i] = m_xGlobalMin[i];
   xGlobalMax[i] = m_xGlobalMax[i];
 }

  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
    arrayView2d< real64 > particlePosition = subRegion.getParticleCenter();

    forAll< serialPolicy >( activeParticleIndices.size(), [&] GEOS_HOST ( localIndex const pp )
    {
      localIndex const p = activeParticleIndices[pp];

      //Check over every dimension
      for( int i=0 ; i<3 ; ++i) //Should this be outside the RAJA::forAll? Are there RAJA calls better or worse than 1
      {
        if(periodic[i])
        {
          particlePosition[p][i] = Mod(particlePosition[p][i]-xGlobalMin[i], xGlobalMax[i]-xGlobalMin[i])+xGlobalMin[i];
        }
      }
    });
  });
}


real64 SolidMechanicsMPM::Mod(real64 num, real64 denom)
{
  if(isZero(denom))
  {
    return num;
  }
  return num - denom * std::floor(num/denom);
}


int SolidMechanicsMPM::combinations( int n, int k )
{
  return factorial( n )/( factorial( k ) * factorial( n - k ) );
}


int SolidMechanicsMPM::factorial( int n )
{
  int val = 1;
  for(int i = 1; i <= n; i++){
    val *= i;
  }
  return val;
}

 
void SolidMechanicsMPM::populateMappingArrays( ParticleManager & particleManager,
                                               NodeManager & nodeManager )
{
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const gridPosition = nodeManager.referencePosition();
  arrayView3d< int const > const ijkMap = m_ijkMap;
  real64 hEl[3] = {0};
  LvArray::tensorOps::copy< 3 >( hEl, m_hEl );
  real64 xLocalMin[3] = {0};
  LvArray::tensorOps::copy< 3 >( xLocalMin, m_xLocalMin );

  localIndex subRegionIndex = 0;
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    // Get views to mapping arrays
    arrayView2d< localIndex > const mappedNodes = m_mappedNodes[subRegionIndex];
    arrayView2d< real64 > const shapeFunctionValues = m_shapeFunctionValues[subRegionIndex];
    arrayView3d< real64 > const shapeFunctionGradientValues = m_shapeFunctionGradientValues[subRegionIndex];

    // Get particle fields
    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
    arrayView2d< real64 const > const particlePosition = subRegion.getParticleCenter();

    // Populate mapping arrays based on particle type
    switch( subRegion.getParticleType() )
    {
      case ParticleType::SinglePoint:
        {
          forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
        {
          localIndex const p = activeParticleIndices[pp];

          // get IJK associated with particle center
          int centerIJK[3];
          for( int i=0; i<3; i++ )
          {
            centerIJK[i] = floor( ( particlePosition[p][i] - xLocalMin[i] ) / hEl[i] );
          }

          // get node IDs, weights and grad weights
          int node = 0;
          int corner = ijkMap[centerIJK[0]][centerIJK[1]][centerIJK[2]];
          auto corner_x = gridPosition[corner];

          real64 xRel = (particlePosition[p][0] - corner_x[0]) / hEl[0];
          real64 yRel = (particlePosition[p][1] - corner_x[1]) / hEl[1];
          real64 zRel = (particlePosition[p][2] - corner_x[2]) / hEl[2];

          for( int i=0; i<2; i++ )
          {
            real64 xWeight = i*xRel + (1-i)*(1.0-xRel);
            real64 dxWeight = i/hEl[0] - (1-i)/hEl[0];
            for( int j=0; j<2; j++ )
            {
              real64 yWeight = j*yRel + (1-j)*(1.0-yRel);
              real64 dyWeight = j/hEl[1] - (1-j)/hEl[1];
              for( int k=0; k<2; k++ )
              {
                real64 zWeight = k*zRel + (1-k)*(1.0-zRel);
                real64 dzWeight = k/hEl[2] - (1-k)/hEl[2];
                mappedNodes[pp][node] = ijkMap[centerIJK[0]+i][centerIJK[1]+j][centerIJK[2]+k];
                shapeFunctionValues[pp][node] = xWeight * yWeight * zWeight;
                shapeFunctionGradientValues[pp][node][0] = dxWeight * yWeight * zWeight;
                shapeFunctionGradientValues[pp][node][1] = xWeight * dyWeight * zWeight;
                shapeFunctionGradientValues[pp][node][2] = xWeight * yWeight * dzWeight;
                node++;
              }
            }
          }
        } );
          break;
        }
      case ParticleType::CPDI:
        {
          int const signs[8][3] = { { -1, -1, -1 },
            {  1, -1, -1 },
            {  1, 1, -1 },
            { -1, 1, -1 },
            { -1, -1, 1 },
            {  1, -1, 1 },
            {  1, 1, 1 },
            { -1, 1, 1 } };
          arrayView3d< real64 const > const particleRVectors = subRegion.getParticleRVectors();
      
          forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
        {
          localIndex const p = activeParticleIndices[pp];

          real64 alpha[8][8];
          real64 cpdiVolume, oneOverV;
          real64 p_r1[3], p_r2[3], p_r3[3]; // allowing 1-indexed r-vectors to persist to torture future postdocs >:)

          for( int i=0; i<3; i++ )
          {
            p_r1[i] = particleRVectors[p][0][i];
            p_r2[i] = particleRVectors[p][1][i];
            p_r3[i] = particleRVectors[p][2][i];
          }

          // We need this because of CPDI domain scaling
          cpdiVolume = 8.0*(-(p_r1[2]*p_r2[1]*p_r3[0]) + p_r1[1]*p_r2[2]*p_r3[0] + p_r1[2]*p_r2[0]*p_r3[1] - p_r1[0]*p_r2[2]*p_r3[1] - p_r1[1]*p_r2[0]*p_r3[2] + p_r1[0]*p_r2[1]*p_r3[2]);
          oneOverV = 1.0 / cpdiVolume;

          alpha[0][0]=oneOverV*(p_r1[2]*p_r2[1] - p_r1[1]*p_r2[2] - p_r1[2]*p_r3[1] + p_r2[2]*p_r3[1] + p_r1[1]*p_r3[2] - p_r2[1]*p_r3[2]);
          alpha[0][1]=oneOverV*(-(p_r1[2]*p_r2[0]) + p_r1[0]*p_r2[2] + p_r1[2]*p_r3[0] - p_r2[2]*p_r3[0] - p_r1[0]*p_r3[2] + p_r2[0]*p_r3[2]);
          alpha[0][2]=oneOverV*(p_r1[1]*p_r2[0] - p_r1[0]*p_r2[1] - p_r1[1]*p_r3[0] + p_r2[1]*p_r3[0] + p_r1[0]*p_r3[1] - p_r2[0]*p_r3[1]);
          alpha[1][0]=oneOverV*(p_r1[2]*p_r2[1] - p_r1[1]*p_r2[2] - p_r1[2]*p_r3[1] - p_r2[2]*p_r3[1] + p_r1[1]*p_r3[2] + p_r2[1]*p_r3[2]);
          alpha[1][1]=oneOverV*(-(p_r1[2]*p_r2[0]) + p_r1[0]*p_r2[2] + p_r1[2]*p_r3[0] + p_r2[2]*p_r3[0] - p_r1[0]*p_r3[2] - p_r2[0]*p_r3[2]);
          alpha[1][2]=oneOverV*(p_r1[1]*p_r2[0] - p_r1[0]*p_r2[1] - p_r1[1]*p_r3[0] - p_r2[1]*p_r3[0] + p_r1[0]*p_r3[1] + p_r2[0]*p_r3[1]);
          alpha[2][0]=oneOverV*(p_r1[2]*p_r2[1] - p_r1[1]*p_r2[2] + p_r1[2]*p_r3[1] - p_r2[2]*p_r3[1] - p_r1[1]*p_r3[2] + p_r2[1]*p_r3[2]);
          alpha[2][1]=oneOverV*(-(p_r1[2]*p_r2[0]) + p_r1[0]*p_r2[2] - p_r1[2]*p_r3[0] + p_r2[2]*p_r3[0] + p_r1[0]*p_r3[2] - p_r2[0]*p_r3[2]);
          alpha[2][2]=oneOverV*(p_r1[1]*p_r2[0] - p_r1[0]*p_r2[1] + p_r1[1]*p_r3[0] - p_r2[1]*p_r3[0] - p_r1[0]*p_r3[1] + p_r2[0]*p_r3[1]);
          alpha[3][0]=oneOverV*(p_r1[2]*p_r2[1] - p_r1[1]*p_r2[2] + p_r1[2]*p_r3[1] + p_r2[2]*p_r3[1] - p_r1[1]*p_r3[2] - p_r2[1]*p_r3[2]);
          alpha[3][1]=oneOverV*(-(p_r1[2]*p_r2[0]) + p_r1[0]*p_r2[2] - p_r1[2]*p_r3[0] - p_r2[2]*p_r3[0] + p_r1[0]*p_r3[2] + p_r2[0]*p_r3[2]);
          alpha[3][2]=oneOverV*(p_r1[1]*p_r2[0] - p_r1[0]*p_r2[1] + p_r1[1]*p_r3[0] + p_r2[1]*p_r3[0] - p_r1[0]*p_r3[1] - p_r2[0]*p_r3[1]);
          alpha[4][0]=oneOverV*(-(p_r1[2]*p_r2[1]) + p_r1[1]*p_r2[2] - p_r1[2]*p_r3[1] + p_r2[2]*p_r3[1] + p_r1[1]*p_r3[2] - p_r2[1]*p_r3[2]);
          alpha[4][1]=oneOverV*(p_r1[2]*p_r2[0] - p_r1[0]*p_r2[2] + p_r1[2]*p_r3[0] - p_r2[2]*p_r3[0] - p_r1[0]*p_r3[2] + p_r2[0]*p_r3[2]);
          alpha[4][2]=oneOverV*(-(p_r1[1]*p_r2[0]) + p_r1[0]*p_r2[1] - p_r1[1]*p_r3[0] + p_r2[1]*p_r3[0] + p_r1[0]*p_r3[1] - p_r2[0]*p_r3[1]);
          alpha[5][0]=oneOverV*(-(p_r1[2]*p_r2[1]) + p_r1[1]*p_r2[2] - p_r1[2]*p_r3[1] - p_r2[2]*p_r3[1] + p_r1[1]*p_r3[2] + p_r2[1]*p_r3[2]);
          alpha[5][1]=oneOverV*(p_r1[2]*p_r2[0] - p_r1[0]*p_r2[2] + p_r1[2]*p_r3[0] + p_r2[2]*p_r3[0] - p_r1[0]*p_r3[2] - p_r2[0]*p_r3[2]);
          alpha[5][2]=oneOverV*(-(p_r1[1]*p_r2[0]) + p_r1[0]*p_r2[1] - p_r1[1]*p_r3[0] - p_r2[1]*p_r3[0] + p_r1[0]*p_r3[1] + p_r2[0]*p_r3[1]);
          alpha[6][0]=oneOverV*(-(p_r1[2]*p_r2[1]) + p_r1[1]*p_r2[2] + p_r1[2]*p_r3[1] - p_r2[2]*p_r3[1] - p_r1[1]*p_r3[2] + p_r2[1]*p_r3[2]);
          alpha[6][1]=oneOverV*(p_r1[2]*p_r2[0] - p_r1[0]*p_r2[2] - p_r1[2]*p_r3[0] + p_r2[2]*p_r3[0] + p_r1[0]*p_r3[2] - p_r2[0]*p_r3[2]);
          alpha[6][2]=oneOverV*(-(p_r1[1]*p_r2[0]) + p_r1[0]*p_r2[1] + p_r1[1]*p_r3[0] - p_r2[1]*p_r3[0] - p_r1[0]*p_r3[1] + p_r2[0]*p_r3[1]);
          alpha[7][0]=oneOverV*(-(p_r1[2]*p_r2[1]) + p_r1[1]*p_r2[2] + p_r1[2]*p_r3[1] + p_r2[2]*p_r3[1] - p_r1[1]*p_r3[2] - p_r2[1]*p_r3[2]);
          alpha[7][1]=oneOverV*(p_r1[2]*p_r2[0] - p_r1[0]*p_r2[2] - p_r1[2]*p_r3[0] - p_r2[2]*p_r3[0] + p_r1[0]*p_r3[2] + p_r2[0]*p_r3[2]);
          alpha[7][2]=oneOverV*(-(p_r1[1]*p_r2[0]) + p_r1[0]*p_r2[1] + p_r1[1]*p_r3[0] + p_r2[1]*p_r3[0] - p_r1[0]*p_r3[1] - p_r2[0]*p_r3[1]);

          // get IJK associated with each corner
          int cornerIJK[8][3]; // CPDI can map to up to 8 cells
          for( int corner=0; corner<8; corner++ )
          {
            for( int i=0; i<3; i++ )
            {
              real64 cornerPositionComponent = particlePosition[p][i] + 
                                               signs[corner][0] * particleRVectors[p][0][i] + 
                                               signs[corner][1] * particleRVectors[p][1][i] + 
                                               signs[corner][2] * particleRVectors[p][2][i];                             

              cornerIJK[corner][i] = std::floor( ( cornerPositionComponent - xLocalMin[i] ) / hEl[i] ); // TODO: Temporarily store the CPDI
                                                                                                        // corners since they're re-used
                                                                                                        // below?
            }
          }

          // get node IDs associated with each corner from IJK map, along with weights and grad weights
          // *** The order in which we access the IJK map must match the order we evaluate the shape functions! ***
          int node = 0;
          for( int corner=0; corner<8; corner++ )
          {
            int cornerNode = ijkMap[cornerIJK[corner][0]][cornerIJK[corner][1]][cornerIJK[corner][2]];
            // GEOS_LOG_RANK("Particle Corner " << corner << " mapped to corner IJK " << cornerIJK[corner][0] << ", "  << cornerIJK[corner][1] << ", " << cornerIJK[corner][2]);
            auto cornerNodePosition = gridPosition[cornerNode];

            real64 x, y, z;
            x = particlePosition[p][0] + signs[corner][0] * particleRVectors[p][0][0] + signs[corner][1] * particleRVectors[p][1][0] + signs[corner][2] * particleRVectors[p][2][0];
            y = particlePosition[p][1] + signs[corner][0] * particleRVectors[p][0][1] + signs[corner][1] * particleRVectors[p][1][1] + signs[corner][2] * particleRVectors[p][2][1];
            z = particlePosition[p][2] + signs[corner][0] * particleRVectors[p][0][2] + signs[corner][1] * particleRVectors[p][1][2] + signs[corner][2] * particleRVectors[p][2][2];

            real64 xRel = (x - cornerNodePosition[0]) / hEl[0];
            real64 yRel = (y - cornerNodePosition[1]) / hEl[1];
            real64 zRel = (z - cornerNodePosition[2]) / hEl[2];

            for( int i=0; i<2; i++ )
            {
              real64 xWeight = i * xRel + (1 - i) * (1.0 - xRel);
              for( int j=0; j<2; j++ )
              {
                real64 yWeight = j * yRel + (1 - j) * (1.0 - yRel);
                for( int k=0; k<2; k++ )
                {
                  real64 zWeight = k * zRel + (1 - k) * (1.0 - zRel);
                  real64 weight = xWeight * yWeight * zWeight;

                  mappedNodes[pp][node] = ijkMap[cornerIJK[corner][0]+i][cornerIJK[corner][1]+j][cornerIJK[corner][2]+k];
                  shapeFunctionValues[pp][node] = 0.125 * weight;
                  shapeFunctionGradientValues[pp][node][0] = alpha[corner][0] * weight;
                  shapeFunctionGradientValues[pp][node][1] = alpha[corner][1] * weight;
                  shapeFunctionGradientValues[pp][node][2] = alpha[corner][2] * weight;
                  node++;
                }
              }
            }
          }
        } );
          break;
        }

      default:
        {
          GEOS_ERROR( "Particle type \"" << subRegion.getParticleType() << "\" is not yet supported." );
          break;
        }
    }

    // Increment the subRegion index
    subRegionIndex++;
  } );
}

//CC: Either need to return or pass body force variables by reference
inline void GEOS_DEVICE SolidMechanicsMPM::computeGeneralizedVortexMMSBodyForce( real64 const time_n,
                                                                                 ParticleManager & particleManager )
{
  // Method of Manufactured Solutions - generalized vortex.
  // As Described in K. Kamojjala, R Brannon, A Sadeghirad, J. Guilkey "Verification Tests
  // in Solid Mechanics", Engineering with computers, (2015).

  localIndex subRegionIndex = 0;
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    // Particle fields
    // ParticleType particleType = subRegion.getParticleType();
    arrayView2d< real64 const > const particlePosition = subRegion.getParticleCenter();
  
    arrayView2d< real64 > const particleBodyForce = subRegion.getField< fields::mpm::particleBodyForce >();

    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();

    // CC: this feels wrong
    // Constitutive model for body force calculation
    string const & solidMaterialName = subRegion.template getReference< string >( viewKeyStruct::solidMaterialNamesString() );
    SolidBase & constitutiveRelation = getConstitutiveModel< SolidBase >( subRegion, solidMaterialName );

    array1d< real64 > shearModulus;
    string constitutiveModelName = constitutiveRelation.getCatalogName();
    if( constitutiveModelName == "Hyperelastic" ){
      Hyperelastic & hyperelastic = dynamic_cast< Hyperelastic & >( constitutiveRelation );
      shearModulus = hyperelastic.shearModulus();
    }

    if( constitutiveModelName == "HyperelasticMMS" ){
      HyperelasticMMS & hyperelasticMMS = dynamic_cast< HyperelasticMMS & >( constitutiveRelation );
      shearModulus = hyperelasticMMS.shearModulus();
    }

    if( constitutiveModelName == "ElasticIsotropic" || constitutiveModelName == "CeramicDamage" || constitutiveModelName == "StrainHardeningPolymer" ){
      ElasticIsotropic & elasticIsotropic = dynamic_cast< ElasticIsotropic & >( constitutiveRelation );
      shearModulus = elasticIsotropic.shearModulus();
    }

    if( constitutiveModelName == "Graphite" || constitutiveModelName == "ElasticTransverseIsotropic" || constitutiveModelName == "ElasticTransverseIsotropicPressureDependent" ){
      ElasticTransverseIsotropic & elasticTransverseIsotropic = dynamic_cast< ElasticTransverseIsotropic & >( constitutiveRelation );
      shearModulus = elasticTransverseIsotropic.effectiveShearModulus();
    }

    GEOS_ERROR_IF( !constitutiveRelation.hasWrapper( constitutive::SolidBase::viewKeyStruct:: defaultDensityString() ) , "Constitutive model must have particle density for the generalized vortex problem!");
    real64 const initialDensity = constitutiveRelation.getReference< real64 >( constitutive::SolidBase::viewKeyStruct:: defaultDensityString() );

    forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST ( localIndex const pp )                                                                                          
    {                                                                                         
      localIndex const p = activeParticleIndices[pp];

      real64 mu = shearModulus[p],
             rho0 = initialDensity, // CC: WRONG, check if this is supposed to be density
             //mu = 384.615384615384585, // Make consistent with elastic properties
             //lambda = 577,             // Make consistent with elastic properties
             //rho0 = 1000.0,            // Make consistent with initial density.
             A = 1.0, // Amplitude of the deformation.
             ri = 0.75,
             ro = 1.25;

      real64 pi = 3.141592653589793;

      real64 x = particlePosition[p][0],
             y = particlePosition[p][1];

      // Radial and angular coordinates:
      real64 R = sqrt( x * x + y * y );

      // CC: Comment from Mike in old geos?
      // I believe the examples in the paper used the reference radius to compute the body force, which
      // greatly increases accuracy.  I've disabled this option, because it is better test to use
      // the reference in the current configuration.
      //    real64 x0 = p_x0( 0 ),
      //           y0 = p_x0( 1 );
      //    real64 R = sqrt( x0 * x0 + y0 * y0 );

      // In either case use the current angle
      real64 Theta = atan2( y, x );

      if( ( R > ri ) && ( R < ro ) )
      {
        // Evaluate some temporary variables:
        real64 p1 = 4096.0 * R * std::pow( 15.0 - 47.0 * R + 48.0 * R * R - 16.0 * R * R * R, 2 ) * mu * std::pow(
            sin( pi * time_n ), 4 ) / rho0,
            p2 = pi * pi * R * std::pow( 15.0 - 32.0 * R + 16.0 * R * R, 4 ) * pow( sin( 2.0 * pi * time_n ), 2 ),
            p3 = -16.0 * ( -45.0 + 188.0 * R - 240.0 * R * R + 96.0 * R * R * R ),
            p4 = -45.0 + 188.0 * R - 240.0 * R * R + 96.0 * R * R * R,
            p5 = std::pow( 15.0 - 32.0 * R + 16.0 * R * R, 2 );

        // Evaluate radial component of the body force:
        real64 br = p1 - p2;

        // Evaluate circumferential component of the body force:
        real64 bt = ( 2.0 * mu * p3 + 2.0 * cos( 2.0 * pi * time_n ) * ( 16.0 * mu * p4 + pi * pi * R * rho0 * p5 ) ) / rho0;

        // Evaluate the rotation angle
        real64 alpha = A * ( 1.0 - cos( 2.0 * pi * time_n ) ) * ( 1.0 - 32.0 * std::pow( R - 1.0, 2 ) + 256.0 * std::pow(
                        R - 1.0, 4 ) ) / 2.0;

        // Evaluate the deformed angular coordinate
        real64 theta = Theta + alpha;

        // evaluate the cartesian components of the body force:
        particleBodyForce[p][0] += br * cos( theta ) - bt * sin( theta );
        particleBodyForce[p][1] += br * sin( theta ) + bt * cos( theta );
        particleBodyForce[p][2] += 0.0;
      }
      } ); // particle loop
    subRegionIndex++;
  } ); // subregion loop

}

inline void GEOS_DEVICE SolidMechanicsMPM::computeBodyForce( ParticleManager & particleManager )
{
  localIndex subRegionIndex = 0;
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    // Particle fields
    // ParticleType particleType = subRegion.getParticleType();
    // arrayView2d< real64 const > const particlePosition = subRegion.getParticleCenter();
  
    arrayView2d< real64 > const particleBodyForce = subRegion.getField< fields::mpm::particleBodyForce >();

    // Map to grid
    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
    // int const numDims = m_numDims;

    forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST ( localIndex const pp )                                                                                          
    {                                                                                         
      localIndex const p = activeParticleIndices[pp];

      particleBodyForce[p][0] = m_bodyForce[0];
      particleBodyForce[p][1] = m_bodyForce[1];
      particleBodyForce[p][2] = m_bodyForce[2];  

    } ); // particle loop
    subRegionIndex++;
  } ); // subregion loop
}


void SolidMechanicsMPM::cofactor( real64 const (& F)[3][3],
                                  real64 (& Fc)[3][3] )
{
  Fc[0][0] = F[1][1] * F[2][2] - F[1][2] * F[2][1];
  Fc[0][1] = F[1][2] * F[2][0] - F[1][0] * F[2][2];
  Fc[0][2] = F[1][0] * F[2][1] - F[1][1] * F[2][0];
  Fc[1][0] = F[0][2] * F[2][1] - F[0][1] * F[2][2];
  Fc[1][1] = F[0][0] * F[2][2] - F[0][2] * F[2][0];
  Fc[1][2] = F[0][1] * F[2][0] - F[0][0] * F[2][1];
  Fc[2][0] = F[0][1] * F[1][2] - F[0][2] * F[1][1];
  Fc[2][1] = F[0][2] * F[1][0] - F[0][0] * F[1][2];
  Fc[2][2] = F[0][0] * F[1][1] - F[0][1] * F[1][0];
}

REGISTER_CATALOG_ENTRY( SolverBase, SolidMechanicsMPM, string const &, dataRepository::Group * const )
}
