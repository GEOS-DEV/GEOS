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
 * @file MPMBase.cpp
 */

#include "MPMBase.hpp"

#include "dataRepository/InputFlags.hpp"
#include "mainInterface/GeosxState.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp" // includes ParticleManager.hpp
#include "finiteElement/FiniteElementDiscretization.hpp" // necessary for FiniteElementDiscretization but not available through FiniteElementDiscretizationManager.hpp?

#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/contact/ContactBase.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/mpiCommunications/SpatialPartition.hpp"

// #include "finiteElement/FiniteElementDiscretizationManager.hpp"

namespace geos
{

using namespace dataRepository;

MPMBase::MPMBase(const string & name, Group * const parent ):
  SolverBase( name, parent ),
  m_solverProfiling( 0 ),
  m_damageFieldPartitioning( 0 ),
  m_surfaceDetection( 0 ),  
  m_directionalOverlapCorrection( 0 ),
  m_needsNeighborList( 0 ),
  m_numDims( 3 ),
  m_planeStrain( 0 ),  
  m_boundaryConditionTypes()
{
  // GEOS_LOG("    MPMBase CTor ... ");

  registerWrapper( "solverProfiling", &m_solverProfiling ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flag for timing subroutines in the solver" );

  registerWrapper( "damageFieldPartitioning", &m_damageFieldPartitioning ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flag for using the gradient of the particle damage field to partition material into separate velocity fields" );

  registerWrapper( "surfaceDetection", &m_surfaceDetection ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flag for automatic surface detection on the 1st cycle" );

  // registerWrapper( "directionalOverlapCorrection", &m_directionalOverlapCorrection ).
  //   setApplyDefaultValue( 0 ).
  //   setInputFlag( InputFlags::OPTIONAL ).
  //   setDescription( "Flag for mitigating pile-up of particles at contact interfaces" );

  registerWrapper( "needsNeighborList", &m_needsNeighborList ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0 ).
    setDescription( "Flag for whether to construct neighbor list" );

  registerWrapper( "numDims", &m_numDims ).
    setApplyDefaultValue( 3 ).
    setInputFlag( InputFlags::FALSE ).
    setDescription( "The number of active spatial dimensions, 2 for plane strain, 3 otherwise" );

  registerWrapper( "planeStrain", &m_planeStrain ).
    setApplyDefaultValue( false ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flag for performing plane strain calculations" );    

  registerWrapper( "boundaryConditionTypes", &m_boundaryConditionTypes ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Boundary conditions on x-, x+, y-, y+, z- and z+ faces. Options are:\n* " + EnumStrings< BoundaryConditionOption >::concat( "\n* " ) );
}

MPMBase::~MPMBase()
{
  // TODO: Auto-generated destructor stub
}

void MPMBase::initializePreSubGroups()
{
  // GEOS_LOG("    MPMBase::initializePreSubGroups() ... ");
  // MPMBase::initializePreSubGroups();

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
        // solidMaterialName = MPMBase::getConstitutiveName< SolidBase >( subRegion );
        solidMaterialName = getConstitutiveName< constitutive::SolidBase >( subRegion );
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

void MPMBase::updateState( DomainPartition & domain )
{
  // TODO: move from SolideMechanicsMPM.cpp
}

void MPMBase::postProcessInput()
{
  // GEOS_LOG("    MPMBase::postProcessInput() ... ");
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
}

void MPMBase::setConstitutiveNamesCallSuper( ParticleSubRegionBase & subRegion ) const
{
  // GEOS_LOG("    MPMBase::setConstitutiveNamesCallSuper(...) ... ");
  SolverBase::setConstitutiveNamesCallSuper( subRegion );

  subRegion.registerWrapper< string >( viewKeyStruct::solidMaterialNamesString() ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setSizedFromParent( 0 );

  string & solidMaterialName = subRegion.getReference< string >( viewKeyStruct::solidMaterialNamesString() );
  solidMaterialName = MPMBase::getConstitutiveName< constitutive::SolidBase >( subRegion );
  GEOS_ERROR_IF( solidMaterialName.empty(), GEOS_FMT( "SolidBase model not found on subregion {}", subRegion.getName() ) );
}

void MPMBase::setConstitutiveNames( ParticleSubRegionBase & subRegion ) const
{
  // GEOS_LOG("    MPMBase::setConstitutiveNames(...) ... ");
  GEOS_UNUSED_VAR( subRegion );
}

void MPMBase::registerDataOnMesh( Group & meshBodies )
{
  // if( m_solverProfiling == -1 )
  // {
  //   GEOS_LOG("    MPMBase::registerDataOnMesh(...) ... ");
  // }

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
        subRegion.registerField< particleOverlap >( getName() );

        // Double-indexed fields (vectors and symmetric tensors stored in Voigt notation)
        subRegion.registerField< particleStress >( getName() ).setDimLabels( 1, voightLabels ).reference().resizeDimension< 1 >( 6 );
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

} // end namespace geos