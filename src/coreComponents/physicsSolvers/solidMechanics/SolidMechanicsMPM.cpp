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

#include "events/mpmEvents/MPMEvents.hpp"

// CC: debug
#include <sstream>

namespace geos
{

using namespace dataRepository;
using namespace constitutive;

// A helper function to calculate polar decomposition. TODO: Previously this was an LvArray method, hopefully it will be again someday.
GEOS_HOST_DEVICE
void polarDecomposition( real64 (& R)[3][3],
                         real64 const (&matrix)[3][3] );

bool compareFloat( real64 a, real64 b, real64 epsilon){
  return std::fabs( a - b ) < epsilon;
}

SolidMechanicsMPM::SolidMechanicsMPM( const string & name,
                                      Group * const parent ):
  SolverBase( name, parent ),
  m_solverProfiling( 0 ),
  m_timeIntegrationOption( TimeIntegrationOption::ExplicitDynamic ),
  m_updateMethod( UpdateMethodOption::FLIP ),
  m_updateOrder( 2 ),
  m_iComm( CommunicationTools::getInstance().getCommID() ),
  m_prescribedBcTable( 0 ),
  m_boundaryConditionTypes(),
  m_bcTable(),
  m_prescribedFTable( 0 ),
  m_prescribedBoundaryFTable( 0 ),
  m_fTableInterpType( SolidMechanicsMPM::InterpolationOption::Linear ),
  m_fTable(),
  m_domainF(),
  m_domainL(),
  m_globalFaceReactions(),
  m_bodyForce(),
  m_stressControl(),
  m_stressTableInterpType( SolidMechanicsMPM::InterpolationOption::Linear ),
  m_stressControlKp( 0.1 ),
  m_stressControlKi( 0.0 ),
  m_stressControlKd( 0.0 ),
  m_domainStress(),
  m_stressControlLastError(),
  m_stressControlITerm(),
  m_boxAverageHistory( 0 ),
  m_reactionHistory( 0 ),
  m_needsNeighborList( 0 ),
  m_neighborRadius( -1.0 ),
  m_binSizeMultiplier( 1 ),
  m_thinFeatureDFGThreshold( DBL_MAX ),
  m_FSubcycles( 1 ),
  m_LBar( 0 ),
  m_LBarScale( 0.0 ),
  m_exactJIntegration( 0 ),
  m_maxParticleVelocity( 1e6 ), // Floating point exception if this is set to DBL_MAX when squared
  m_maxParticleVelocitySquared( DBL_MAX ),
  m_minParticleJacobian( 0.1 ),
  m_maxParticleJacobian( 10.0 ),
  m_overlapCorrection( 0 ),
  m_overlapThreshold1( 1.00 ),
  m_overlapThreshold2( 1.10 ),
  m_computeSPHJacobian( 0 ),
  m_initialTemperature( 300.0 ),
  m_shockHeating( 0 ),
  m_useArtificialViscosity( 0 ),
  m_artificialViscosityQ0( 0.0 ),
  m_artificialViscosityQ1( 0.0 ),
  m_cpdiDomainScaling( 0 ),
  m_smallMass( DBL_MAX ),
  m_numContactGroups( 0 ),
  m_numContactFlags(),
  m_numVelocityFields(),
  m_separabilityMinDamage( 0.5 ),
  m_treatFullyDamagedAsSingleField( 1 ),
  m_surfaceDetection( 0 ),
  m_damageFieldPartitioning( 0 ),
  m_contactNormalType( 1 ),
  m_contactGapCorrection( 0 ),
  m_resetDefGradForFullyDamagedParticles( 0 ),
  m_plotUnscaledParticles( 0 ),
  // m_directionalOverlapCorrection( 0 ),
  m_frictionCoefficient(),
  m_frictionCoefficientTable(),
  m_planeStrain( 0 ),
  m_numDims( 3 ),
  m_ijkMap(),
  m_useEvents( 0 ),
  m_mpmEventManager( nullptr ),
  m_surfaceHealing( false ),
  m_debugFlag( 0 ),
  m_computeXProfile( 0 ),
  m_xProfileWriteInterval( 0.0 ),
  m_nextXProfileWriteTime( 0.0 ),
  m_xProfileVx0( 0.0 )
{
  GEOS_MARK_FUNCTION;

  registerWrapper( "solverProfiling", &m_solverProfiling ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    // setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Flag for timing subroutines in the solver" );

  registerWrapper( viewKeyStruct::timeIntegrationOptionString(), &m_timeIntegrationOption ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_timeIntegrationOption ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Time integration method. Options are:\n* " + EnumStrings< TimeIntegrationOption >::concat( "\n* " ) );

  registerWrapper( "updateMethod", &m_updateMethod ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_updateMethod ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Update method. Options are:\n* " + EnumStrings< UpdateMethodOption >::concat( "\n* ") );

  registerWrapper( "updateOrder", &m_updateOrder ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_updateOrder ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription("Order for update method, only applies to XPIC and FMPM");

  registerWrapper( "prescribedBcTable", &m_prescribedBcTable ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_prescribedBcTable ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Flag for whether to have time-dependent boundary condition types" );

  registerWrapper( "boxAverageHistory", &m_boxAverageHistory ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_boxAverageHistory ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Flag for whether to output box average history" );

  registerWrapper( "boxAverageWriteInterval", &m_boxAverageWriteInterval ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_boxAverageWriteInterval ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Interval between writing box averages to files" );

  registerWrapper( "boxAverageMin", &m_boxAverageMin).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Minimum corner position of box average" );

  registerWrapper( "boxAverageMax", &m_boxAverageMax).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Maximum corner position of box average" );

  registerWrapper( "nextBoxAverageWriteTime", &m_nextBoxAverageWriteTime ).
    setInputFlag( InputFlags::FALSE ).
    setApplyDefaultValue( 0 ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Next time to write box averages" );

  registerWrapper( "reactionHistory", &m_reactionHistory ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0 ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Flag for whether to output face reaction history" );

  registerWrapper( "reactionWriteInterval", &m_reactionWriteInterval ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0 ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Interval between writing reactions to files" );

  registerWrapper( "nextReactionWriteTime", &m_nextReactionWriteTime ).
    setInputFlag( InputFlags::FALSE ).
    setApplyDefaultValue( 0 ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Next time to write reactions" );

  registerWrapper( "boundaryConditionTypes", &m_boundaryConditionTypes ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Boundary conditions on x-, x+, y-, y+, z- and z+ faces. Options are:\n* " + EnumStrings< BoundaryConditionOption >::concat( "\n* " ) );

  registerWrapper( "bodyForce", &m_bodyForce ).
    setInputFlag(InputFlags::OPTIONAL).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Array that stores uniform body force" );

  registerWrapper( "bcTable", &m_bcTable ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Array that stores time-dependent bc types on x-, x+, y-, y+, z- and z+ faces." );

  registerWrapper( "prescribedFTable", &m_prescribedFTable ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Flag for whether to have time-dependent superimposed velocity gradient for triply periodic simulations" );

  registerWrapper( "prescribedBoundaryFTable", &m_prescribedBoundaryFTable ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Flag for whether to have time-dependent boundary conditions described by a global background grid F" );

  registerWrapper( "fTableInterpType", &m_fTableInterpType ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setApplyDefaultValue( m_fTableInterpType ).
    setDescription( "The type of F table interpolation. Options are 0 (linear), 1 (cosine), 2 (quintic polynomial)." );

  registerWrapper( "fTable", &m_fTable ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Array that stores time-dependent grid-aligned stretches interpreted as a gloabl background grid F read from the XML file." );

  registerWrapper( "stressControl" , &m_stressControl).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Flag for whether stress control using box averages is enabled" );

  registerWrapper( "stressTableInterpType", &m_stressTableInterpType ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setApplyDefaultValue( m_stressTableInterpType ).
    setDescription( "The type of stress table interpolation. Options are 0 (linear), 1 (cosine), 2 (quintic polynomial)." );

  registerWrapper( "stressTable", &m_stressTable ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Array that stores the time-depended grid aligned stresses" );

  registerWrapper( "stressControlKp", &m_stressControlKp ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Proportional gain of stress PID controller" );

  registerWrapper( "stressControlKi", &m_stressControlKi ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Integral gain of stress PID controller" );

   registerWrapper( "stressControlKd", &m_stressControlKd ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Derivative gain of stress PID controller" );

  registerWrapper( "domainStress", &m_domainStress ).
    setInputFlag( InputFlags::FALSE ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Stores current target domain stress as driven by stress control" );

  registerWrapper( "stressControlLastError", &m_stressControlLastError ).
    setInputFlag( InputFlags::FALSE ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Stress control error from previous timestep" );

  registerWrapper( "stressControlITerm", &m_stressControlITerm ).
    setInputFlag( InputFlags::FALSE ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Stress control integral term" );

  registerWrapper( "needsNeighborList", &m_needsNeighborList ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0 ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Flag for whether to construct neighbor list" );

  registerWrapper( "neighborRadius", &m_neighborRadius ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( -1.0 ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Neighbor radius for SPH-type calculations" );

  registerWrapper( "binSizeMultiplier", &m_binSizeMultiplier ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 1 ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Multiplier for setting bin size, used to speed up particle neighbor sorting" );

  registerWrapper( "thinFeatureDFGThreshold", &m_thinFeatureDFGThreshold ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_thinFeatureDFGThreshold ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Threshold to treat relatively thin ( compared to grid spacing ) damaged material to avoid spurious slip surfaces" );

  registerWrapper( "useDamageAsSurfaceFlag", &m_useDamageAsSurfaceFlag ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0 ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Indicates whether particle damage at the beginning of the simulation should be interpreted as a surface flag" );

  registerWrapper( "FSubCycles", &m_FSubcycles ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_FSubcycles ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Number of sub cycles to more accurately integrate the deformation gradient" );

  registerWrapper( "LBar", &m_LBar ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_LBar ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Used to choose the Lbar method. 0 = do nothing, 1 = run trD through the SPH kernel, 2 = straight-up volume average" );

  registerWrapper( "LBarScale", &m_LBarScale ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_LBarScale ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Scales the contribution of volume-averaged trD to particle L, alleviates checkerboarding" );

  registerWrapper( "exactJIntegration", &m_exactJIntegration ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_exactJIntegration ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Will force integration of F to have an exact integral of J." );

  registerWrapper( "maxParticleVelocity", &m_maxParticleVelocity ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDefaultValue( 1e6 ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Velocity above which particles are deleted" );

  registerWrapper( "maxParticleVelocitySquared", &m_maxParticleVelocitySquared).
    setInputFlag( InputFlags::FALSE ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Square of the max particle velocity" );

  registerWrapper( "minParticleJacobian", &m_minParticleJacobian ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDefaultValue( m_minParticleJacobian ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Jacobian value below which particles are deleted" );

  registerWrapper( "maxParticleJacobian", &m_maxParticleJacobian ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDefaultValue( m_maxParticleJacobian ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Jacobian value above which particles are deleted" );

  registerWrapper( "overlapCorrection", &m_overlapCorrection ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDefaultValue( m_overlapCorrection ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Overlap correction flag (0=off, 1=increase normal force to remove gap. (not fully developed), 2=compute the SPH Jacobian and use it to scale particle density to improve the overlap. )" );

  registerWrapper( "overlapThreshold1", &m_overlapThreshold1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDefaultValue( m_overlapThreshold1 ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Overlap correction threshold 1" );

  registerWrapper( "overlapThreshold2", &m_overlapThreshold2 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDefaultValue( m_overlapThreshold2 ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Overlap correction threshold 2" );

  registerWrapper( "computeSPHJacobian", &m_computeSPHJacobian ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDefaultValue( m_computeSPHJacobian ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Overlap correction flag (0=off, 1=increase normal force to remove gap. (not fully developed), 2=compute the SPH Jacobian and use it to scale particle density to improve the overlap. )" );

  registerWrapper( "initialTemperature", &m_initialTemperature ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDefaultValue( m_initialTemperature ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Initial particle temperature" );

  registerWrapper( "shockHeating", &m_shockHeating ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDefaultValue( m_shockHeating ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Initial particle temperature" );

  registerWrapper( "useArtificialViscosity", &m_useArtificialViscosity ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDefaultValue( m_useArtificialViscosity ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Flag to enable artificial viscosity" );

  registerWrapper( "artificialViscosityQ0", &m_artificialViscosityQ0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDefaultValue( m_artificialViscosityQ0 ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Initial particle temperature" );

  registerWrapper( "cpdiDomainScaling", &m_cpdiDomainScaling ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDefaultValue( 0 ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Option for CPDI domain scaling" );

  registerWrapper( "generalizedVortexMMS", &m_generalizedVortexMMS ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDefaultValue( 0 ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Option for CPDI domain scaling" );

  registerWrapper( "smallMass", &m_smallMass ).
    setInputFlag( InputFlags::FALSE ).
    setApplyDefaultValue( DBL_MAX ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "The small mass threshold for ignoring extremely low-mass nodes." );

  registerWrapper( "numContactGroups", &m_numContactGroups ).
    setInputFlag( InputFlags::FALSE ).
    setApplyDefaultValue( m_numContactGroups ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Number of prescribed contact groups" );

  registerWrapper( "numContactFlags", &m_numContactFlags ).
    setInputFlag( InputFlags::FALSE ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Number of contact flags that can appear due to damage" );

  registerWrapper( "numVelocityFields", &m_numVelocityFields ).
    setInputFlag( InputFlags::FALSE ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Number of velocity fields" );

  registerWrapper( "separabilityMinDamage", &m_separabilityMinDamage ).
    setApplyDefaultValue( m_separabilityMinDamage ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Damage threshold for field separability" );

  registerWrapper( "treatFullyDamagedAsSingleField", &m_treatFullyDamagedAsSingleField ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_treatFullyDamagedAsSingleField ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Whether to consolidate fully damaged fields into a single field. Nice for modeling damaged mush." );

  registerWrapper( "surfaceDetection", &m_surfaceDetection ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Flag for automatic surface detection on the 1st cycle" );

  registerWrapper( "damageFieldPartitioning", &m_damageFieldPartitioning ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Flag for using the gradient of the particle damage field to partition material into separate velocity fields" );

  registerWrapper( "contactNormalType", &m_contactNormalType ).
    setApplyDefaultValue( m_contactNormalType ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Flag for contact normal type" );

  registerWrapper( "contactGapCorrection", &m_contactGapCorrection ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_contactGapCorrection ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Flag for mitigating contact gaps" );

  registerWrapper( "resetDefGradForFullyDamagedParticles", &m_resetDefGradForFullyDamagedParticles ).
    setApplyDefaultValue( m_resetDefGradForFullyDamagedParticles ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Flag for resetting deformation gradient of fully damaged particles" );

  registerWrapper( "plotUnscaledParticles", &m_plotUnscaledParticles ).
    setApplyDefaultValue( m_plotUnscaledParticles ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Flag for plotting particles without CPDI domain scaling" );

  // registerWrapper( "directionalOverlapCorrection", &m_directionalOverlapCorrection ).
  //   setApplyDefaultValue( 0 ).
  //   setInputFlag( InputFlags::OPTIONAL ).
  //   setRestartFlags( RestartFlags::NO_WRITE ).
  //   setDescription( "Flag for mitigating pile-up of particles at contact interfaces" );

  registerWrapper( "frictionCoefficient", &m_frictionCoefficient ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Coefficient of friction, currently assumed to be the same everywhere" );

  registerWrapper( "frictionCoefficientTable", &m_frictionCoefficientTable ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Friction coefficient table for different groups" );

  registerWrapper( "planeStrain", &m_planeStrain ).
    setApplyDefaultValue( m_planeStrain ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Flag for performing plane strain calculations" );

  registerWrapper( "numDims", &m_numDims ).
    setApplyDefaultValue( m_numDims ).
    setInputFlag( InputFlags::FALSE ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "The number of active spatial dimensions, 2 for plane strain, 3 otherwise" );

  registerWrapper( "m_ijkMap", &m_ijkMap ).
    setInputFlag( InputFlags::FALSE ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Map from indices in each spatial dimension to local node ID" );

  registerWrapper("useEvents", &m_useEvents).
    setInputFlag( InputFlags::OPTIONAL ).
    setDefaultValue( m_useEvents ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Enable events" );

  registerWrapper( "debugFlag", &m_debugFlag ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDefaultValue( m_debugFlag ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Enables debugging of MPM explicit timestep" );

  registerWrapper( "computeXProfile", &m_computeXProfile ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDefaultValue( m_computeXProfile ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Flag to enable x profiling" );

  registerWrapper( "xProfileWriteInterval", &m_xProfileWriteInterval ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDefaultValue( m_xProfileWriteInterval ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Interval between writing x profiling data" );

  registerWrapper( "nextXProfileWriteTime", &m_nextXProfileWriteTime ).
    setInputFlag( InputFlags::FALSE ).
    setDefaultValue( m_nextXProfileWriteTime ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Next x profile writing time" );

  registerWrapper( "xProfileVx0", &m_xProfileVx0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDefaultValue( m_xProfileVx0 ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Velocity to subtract from x profile" );

  registerWrapper( "elementSize", &m_hEl ).
    setInputFlag( InputFlags::FALSE).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Minimum element size in x, y and z" );
  
  registerWrapper( "localMinimum", &m_xLocalMin ).
    setInputFlag( InputFlags::FALSE).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "local minimum" );
  
  registerWrapper( "localMaximum", &m_xLocalMax ).
    setInputFlag( InputFlags::FALSE).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "local maximum" );

    registerWrapper( "localMinimumNoGhost", &m_xLocalMinNoGhost ).
    setInputFlag( InputFlags::FALSE).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "local minimum without ghost cells" );
  
  registerWrapper( "localMaximumNoGhost", &m_xLocalMaxNoGhost ).
    setInputFlag( InputFlags::FALSE).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "local maximum without ghost cells" );

  registerWrapper( "globalMinimum", &m_xGlobalMin ).
    setInputFlag( InputFlags::FALSE).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "global minimum" );
  
  registerWrapper( "globalMaximum", &m_xGlobalMax ).
    setInputFlag( InputFlags::FALSE).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "global maximum" );

  registerWrapper( "partitionExtent", &m_partitionExtent ).
    setInputFlag( InputFlags::FALSE).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "parititon extent" );

  registerWrapper( "domainExtent", &m_domainExtent ).
    setInputFlag( InputFlags::FALSE).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "domain extent" );

  registerWrapper( "domainF", &m_domainF ).
    setInputFlag( InputFlags::FALSE).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "domain deformation gradient" );

  registerWrapper( "domainL", &m_domainL).
    setInputFlag( InputFlags::FALSE).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "domain L" );

  registerWrapper( "globalFaceReactions", &m_globalFaceReactions).
    setInputFlag( InputFlags::FALSE).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Stores the global face reactions for stress control between timesteps" );

  registerWrapper( "numElements", &m_nEl).
    setInputFlag( InputFlags::FALSE).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "number of elements along partition directions" );
  
  m_mpmEventManager = &registerGroup< MPMEventManager >( groupKeys.mpmEventManager );
}

void SolidMechanicsMPM::postRestartInitialization()
{
  GEOS_MARK_FUNCTION;
  SolverBase::postRestartInitialization();

  // Initialize friction coefficient table
  if(m_numContactGroups > 0 )
  {
    initializeFrictionCoefficients();
  }

  // Initialize box average history file and write its header
  if( m_boxAverageHistory == 1 )
  {
    // Initialize box average extent
    if(  m_boxAverageMin.size() == 0  )
    {
      m_boxAverageMin.resize( 3 );
      LvArray::tensorOps::copy< 3 >( m_boxAverageMin, m_xGlobalMin );
    }
    else
    {
      GEOS_ERROR_IF( m_boxAverageMin.size() != 3, "Box average min must have 3 elements" );
    }

    if(  m_boxAverageMax.size() == 0  )
    {
      m_boxAverageMax.resize( 3 );
      LvArray::tensorOps::copy< 3 >( m_boxAverageMax, m_xGlobalMax );
    }
    else
    {
      GEOS_ERROR_IF( m_boxAverageMax.size() != 3, "Box average min must have 3 elements" );
    }

    GEOS_ERROR_IF( ( m_boxAverageMin[0] > m_boxAverageMax[0] ) || ( m_boxAverageMin[1] > m_boxAverageMax[1] ) || ( m_boxAverageMin[2] > m_boxAverageMax[2] ) , "Box minimums must be less than box maximums");
  }
}

void SolidMechanicsMPM::postProcessInput()
{
  GEOS_MARK_FUNCTION;
  SolverBase::postProcessInput();

  if( m_overlapCorrection == 2 )
  {
    GEOS_LOG_RANK_0_IF( m_computeSPHJacobian != 1, "Warning! overlapCorrection=2 sets computeSPHJacobian=1" );
    m_computeSPHJacobian = 1;
  }

  // Activate neighbor list if necessary
  if( m_damageFieldPartitioning == 1 || m_surfaceDetection == 1 || m_computeSPHJacobian > 0 || m_LBar > 0 /*|| m_directionalOverlapCorrection == 1*/ )
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

    // Read and distribute BC table
  if( m_prescribedBcTable == 1 )
  {
    // Reads the FTable directly from the xml
    int numRows = m_bcTable.size( 0 );

    GEOS_ERROR_IF(numRows == 0, "Prescribed boundary conditions is enabled but no bcTable was specified.");
    for(int i = 0; i < numRows; ++i){
      GEOS_ERROR_IF(m_bcTable[i].size() != 7, "BCtable row " << i+1 << " must have 7 elements.");
      
      GEOS_ERROR_IF(m_bcTable[i][0] < 0, "BCTable times must be positive.");

      GEOS_ERROR_IF( compareFloat( roundf(m_bcTable[i][1]), m_bcTable[i][1], 1e-12 ) || 
                     compareFloat( roundf(m_bcTable[i][2]), m_bcTable[i][1], 1e-12 ) || 
                     compareFloat( roundf(m_bcTable[i][3]), m_bcTable[i][1], 1e-12 ) || 
                     compareFloat( roundf(m_bcTable[i][4]), m_bcTable[i][1], 1e-12 ) || 
                     compareFloat( roundf(m_bcTable[i][5]), m_bcTable[i][1], 1e-12 ) || 
                     compareFloat( roundf(m_bcTable[i][6]), m_bcTable[i][1], 1e-12 ), "Only integer boundary condition types are permitted.");

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
  if( m_domainF.size() == 0 )
  {
    m_domainF.resize( 3 );
    LvArray::tensorOps::fill< 3 >(m_domainF, 1.0);
  }

  if( m_domainL.size() == 0 )
  {
    m_domainL.resize( 3 );
    LvArray::tensorOps::fill< 3 >(m_domainL, 0.0);
  }

  if( m_prescribedBoundaryFTable == 1 && m_prescribedFTable == 1 )
  {
    // Reads the FTable directly from the xml
    int numRows = m_fTable.size( 0 );
    GEOS_ERROR_IF(numRows == 0, "Prescribed boundary deformation is enabled but no F table was specified.");
    for(int i = 0; i < numRows; ++i){
      GEOS_ERROR_IF(m_fTable[i].size() != 4, "F table row " << i+1 << " must have 4 elements.");
      GEOS_ERROR_IF(m_fTable[i][0] < 0, "F table times must be positive.");

      if(i == 0)
      {
        GEOS_ERROR_IF( compareFloat( m_fTable[i][1], 1.0, 1e-12 ) || 
                       compareFloat( m_fTable[i][2], 1.0, 1e-12 ) || 
                       compareFloat( m_fTable[i][3], 1.0, 1e-12 ) , "Deformation of first row of F table must be 1." );
      }
      else
      {
        GEOS_ERROR_IF( m_fTable[i][1] < 0 || 
                       m_fTable[i][2] < 0 || 
                       m_fTable[i][3] < 0, "F table entries must be positive." );
      }
    }
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
  else 
  {
    GEOS_ERROR_IF( m_bodyForce.size() != 3, "Body force input must have size 3." );
  }
  
  // Initialze reactions for each face of box to 0
  if( m_globalFaceReactions.size() == 0 )
  {
    m_globalFaceReactions.resize( 6 );
    LvArray::tensorOps::fill< 6 >( m_globalFaceReactions, 0.0 );
  }

  // Check stress control
  if( m_stressControl.size() == 0 ){
    m_stressControl.resize( 3 );
    LvArray::tensorOps::fill<3>( m_stressControl, 0 );
  }
  else
  {
    GEOS_ERROR_IF( m_stressControl.size() != 3, "Stress control input must have size 3.");
  }

  if( m_domainStress.size() == 0 )
  {
    m_domainStress.resize( 3 );
    LvArray::tensorOps::fill< 3 >( m_domainStress, 0.0 );
  }

  if( m_stressControlLastError.size() == 0 )
  {
    m_stressControlLastError.resize(3);
    LvArray::tensorOps::fill< 3 >( m_stressControlLastError, 0.0 );
  }

  if( m_stressControlITerm.size() == 0 )
  {
    m_stressControlITerm.resize(3);
    LvArray::tensorOps::fill< 3 >( m_stressControlITerm, 0.0 );
  }
    
  if( m_stressControl[0] == 1 || m_stressControl[1] == 1 || m_stressControl[2] == 1 )
  {
    GEOS_ERROR_IF(m_stressTable.size(0) == 0, "Stress table cannot be empty if stress control is enabled");
    GEOS_ERROR_IF(m_stressTable.size(1) == 0, "Stress table must have 4 columns");

    for(int i = 0; i < m_stressTable.size(0) ; i++)
    {
      GEOS_ERROR_IF(m_stressTable[i][0] < 0.0, "Stress table times must be positive");
    }
  }

}

void SolidMechanicsMPM::registerDataOnMesh( Group & meshBodies )
{
  GEOS_MARK_FUNCTION;
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
        subRegion.registerField< particleCrystalHealFlag >( getName() );
        subRegion.registerField< particleMass >( getName() );
        subRegion.registerField< particleWavespeed >( getName() );
        subRegion.registerField< particleHeatCapacity >( getName() );
        subRegion.registerField< particleTemperature >( getName() );
        subRegion.registerField< particleInternalEnergy >( getName() );
        subRegion.registerField< particleKineticEnergy >( getName() );
        subRegion.registerField< particleArtificialViscosity >( getName() );
        subRegion.registerField< particleInitialVolume >( getName() );
        subRegion.registerField< particleDensity >( getName() );
        subRegion.registerField< particleOverlap >( getName() );
        subRegion.registerField< particleSPHJacobian >( getName() );
        subRegion.registerField< particlePorosity >( getName() );
        subRegion.registerField< particleReferencePorosity >( getName() );

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

      nodes.registerWrapper< array3d< real64 > >( viewKeyStruct::displacementString() ).
        setPlotLevel( PlotLevel::LEVEL_0 ).
        setRegisteringObjects( this->getName() ).
        setDescription( "An array that hold the displacement on the nodes for each field");

      nodes.registerWrapper< array3d< real64 > >( viewKeyStruct::particleSurfaceNormalString() ).
        setPlotLevel( PlotLevel::LEVEL_0 ).
        setRegisteringObjects( this->getName() ).
        setDescription( "An array that holds the surface normals interpoalted from particles for each field");

      nodes.registerWrapper< array2d< real64 > >( viewKeyStruct::massString() ).
        setPlotLevel( PlotLevel::LEVEL_0 ).
        setRegisteringObjects( this->getName() ).
        setDescription( "An array that holds the mass on the nodes." );

      nodes.registerWrapper< array2d< real64 > >( viewKeyStruct::materialVolumeString() ).
        setPlotLevel( PlotLevel::LEVEL_0 ).
        setRegisteringObjects( this->getName() ).
        setDescription( "An array that holds the material volume on the nodes." );

      nodes.registerWrapper< array3d< real64 > >( viewKeyStruct::velocityString() ).
        setPlotLevel( PlotLevel::LEVEL_0 ).
        setRegisteringObjects( this->getName() ).
        setDescription( "An array that holds the current velocity on the nodes." );

      nodes.registerWrapper< array3d< real64 > >( viewKeyStruct::dVelocityString() ).
        setPlotLevel( PlotLevel::LEVEL_0 ).
        setRegisteringObjects( this->getName() ).
        setDescription( "An array that holds the change in velocity from enforcing contact." );

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

      nodes.registerWrapper< array3d< real64 > >( viewKeyStruct::normalStressString() ).
        setPlotLevel( PlotLevel::LEVEL_1 ).
        setRegisteringObjects( this->getName() ).
        setDescription( "An array that holds the result of mapping particle normal stresses to the nodes for x profiling." );     

      nodes.registerWrapper< array2d< real64 > >( viewKeyStruct::massWeightedDamageString() ).
        setPlotLevel( PlotLevel::LEVEL_1 ).
        setRegisteringObjects( this->getName() ).
        setDescription( "An array that holds the result of mapping particle mass weighted damage to the nodes for x profiling." );     

      nodes.registerWrapper< array3d< real64 > >( viewKeyStruct::vPlusString() ).
        setPlotLevel( PlotLevel::LEVEL_1 ).
        setRegisteringObjects( this->getName() ).
        setDescription( "An array that holds the result of each XPIC and FMPM order iteration" );   

      nodes.registerWrapper< array3d< real64 > >( viewKeyStruct::dVPlusString() ).
        setPlotLevel( PlotLevel::LEVEL_1 ).
        setRegisteringObjects( this->getName() ).
        setDescription( "An array that holds the result of each XPIC and FMPM order iteration for multifield contact" );   

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
  GEOS_MARK_FUNCTION;
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
  GEOS_MARK_FUNCTION;

  // Initialize neighbor lists
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    OrderedVariableToManyParticleRelation & neighborList = subRegion.neighborList();
    neighborList.setParticleManager( particleManager );
  } );

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
  m_xLocalMin.resize(3);
  LvArray::tensorOps::fill< 3 >(m_xLocalMin, DBL_MAX);
  m_xLocalMax.resize(3);
  LvArray::tensorOps::fill< 3 >(m_xLocalMax, -DBL_MAX);
  for( int g=0; g<numNodes; g++ )
  {
    for( int i=0; i<3; i++ )
    {
      m_xLocalMin[i] = std::fmin( m_xLocalMin[i], gridPosition[g][i] );
      m_xLocalMax[i] = std::fmax( m_xLocalMax[i], gridPosition[g][i] );
    }
  }

  // m_xLocalMin0.resize( 3 );
  // LvArray::tensorOps::copy< 3 >( m_xLocalMin0, m_xLocalMin);
  // m_xLocalMax0.resize( 3 );
  // LvArray::tensorOps::copy< 3 >( m_xLocalMax0, m_xLocalMax);

  m_xLocalMinNoGhost.resize(3);
  m_xLocalMaxNoGhost.resize(3);
  m_partitionExtent.resize(3);
  for( int i=0; i<3; i++ )
  {
    m_xLocalMinNoGhost[i] = partition.getLocalMin()[i];
    m_xLocalMaxNoGhost[i] = partition.getLocalMax()[i];
    m_partitionExtent[i] = m_xLocalMax[i] - m_xLocalMin[i];
  }

  // CC: why not compute element size directly from domain extent and number of 
  // pps across direction?

  // Get element size
  m_hEl.resize(3);
  LvArray::tensorOps::fill< 3 >(m_hEl, DBL_MAX);
  for( int g=0; g<numNodes; g++ )
  {
    for( int i=0; i<3; i++ )
    {
      // By definition, this should always be positive. Furthermore, the 
      // gridPosition should only be those on the local partition
      real64 test = gridPosition[g][i] - m_xLocalMin[i];  

      // We're looking for the smallest nonzero distance from the "min" node. 
      // TODO: Could be vulnerable to a finite precision bug.       
      if( test > 0.0 ) 
      {
        m_hEl[i] = std::fmin( test, m_hEl[i] );
      }
    }
  }

  // m_hEl0.resize( 0 );
  // LvArray::tensorOps::copy< 3 >( m_hEl0, m_hEl);

  // Set SPH neighbor radius if necessary
  if( m_neighborRadius <= 0.0 )
  {
    if( m_planeStrain )
    {
      m_neighborRadius *= -1.0 * sqrt( 
        m_hEl[0] * m_hEl[0] + m_hEl[1] * m_hEl[1] );
    }
    else
    {
      m_neighborRadius *= -1.0 * sqrt( 
        m_hEl[0] * m_hEl[0] + m_hEl[1] * m_hEl[1] + m_hEl[2] * m_hEl[2] );
    }
  }

  // Get global domain extent excluding buffer nodes
  m_xGlobalMin.resize(3);
  m_xGlobalMax.resize(3);
  m_domainExtent.resize(3);
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
  m_nEl.resize( 3 );
  for( int i=0; i<3; i++ )
  {
    m_nEl[i] = std::round( m_partitionExtent[i] / m_hEl[i] );
  }

  // Create element map
  // loop through all nodes, i.e. numElems+1, identify closest integer of 
  // hEl (stepsize) multiple.  m_ijkMap[0][0][0] -> [xmin, ymin, zmin]
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

    // dir0 will be 0, 0, 1, 1, 2, 2 or associated axis of(x-, x+, y-, y+, z-, z+)
    int dir0 = face / 2; 

    // even -> -1 and odd -> +1 to get the sign of the direction i.e. x-, x+ ,...
    int positiveNormal = face % 2; 
    int sign = 2 * positiveNormal - 1;

    // minOrMax is globalMin if normal is positive, else globalMax
    // So that face=0 gives x- and xmin
    real64 minOrMax = positiveNormal == 0 ? m_xGlobalMin[dir0] : m_xGlobalMax[dir0];

    for( localIndex g=0; g<numNodes; g++ )
    {
      real64 positionRelativeToBoundary = gridPosition[g][dir0] - minOrMax;
      // small multiple of element dimension
      // CC::Jay: why is 1e-12 a reasonable choice? either co-located with boundary
      //     or potential buffer node?
      real64 tolerance = 1.0e-12 * m_hEl[dir0]; 

      if( std::fabs( positionRelativeToBoundary ) < tolerance )
      {
        tmpBoundaryNodes.insert( g );
      }

      // basically a dot product with the face normal to ensure that point is 
      // outside of the domain and not co-located with the boundary
      if( sign * positionRelativeToBoundary > 0 && std::fabs( positionRelativeToBoundary ) > tolerance ) 
      {
        tmpBufferNodes.insert( g );
      }
    }

    // CC::Jay: boundary and buffer nodes seem to define background grid elements
    //     and are defined by element faces.  Is that correct?
    m_boundaryNodes[face].insert( tmpBoundaryNodes.begin(), tmpBoundaryNodes.end() );
    m_bufferNodes[face].insert( tmpBufferNodes.begin(), tmpBufferNodes.end() );
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

    arrayView2d< real64 > const particleBodyForce = subRegion.getField< fields::mpm::particleBodyForce >();
    arrayView2d< real64 > const particlePlasticStrain = subRegion.getField< fields::mpm::particlePlasticStrain >(); 
    arrayView1d< real64 > const particleDensity = subRegion.getField< fields::mpm::particleDensity >();
    arrayView1d< real64 > const particleMass = subRegion.getField< fields::mpm::particleMass >();
    arrayView1d< real64 > const particleHeatCapacity = subRegion.getField< fields::mpm::particleHeatCapacity >();
    arrayView1d< real64 > const particleTemperature = subRegion.getField< fields::mpm::particleTemperature >();
    arrayView1d< real64 > const particleInternalEnergy = subRegion.getField< fields::mpm::particleInternalEnergy >();
    arrayView1d< real64 > const particleKineticEnergy = subRegion.getField< fields::mpm::particleKineticEnergy >();
    arrayView1d< real64 > const particleArtificialViscosity = subRegion.getField< fields::mpm::particleArtificialViscosity >();
    arrayView3d< real64 > const particleDeformationGradient = subRegion.getField< fields::mpm::particleDeformationGradient >();
    arrayView3d< real64 > const particleFDot = subRegion.getField< fields::mpm::particleFDot >();
    arrayView3d< real64 > const particleVelocityGradient = subRegion.getField< fields::mpm::particleVelocityGradient >();
    arrayView1d< real64 > const particleInitialVolume = subRegion.getField< fields::mpm::particleInitialVolume >();
    arrayView3d< real64 > const particleInitialRVectors = subRegion.getField< fields::mpm::particleInitialRVectors >();
    arrayView3d< real64 > const particleSphF = subRegion.getField< fields::mpm::particleSphF >();
    arrayView1d< real64 > const particleSPHJacobian = subRegion.getField< fields::mpm::particleSPHJacobian >();
    arrayView2d< real64 > const particleReferencePosition = subRegion.getField< fields::mpm::particleReferencePosition >();
    arrayView1d< int > const particleCrystalHealFlag = subRegion.getField< fields::mpm::particleCrystalHealFlag >();
    arrayView1d< real64 > const particlePorosity = subRegion.getField< fields::mpm::particlePorosity >();
    arrayView1d< real64 > const particleReferencePorosity = subRegion.getField< fields::mpm::particleReferencePorosity >();

    // Set reference position, volume and R-vectors
    for( int p=0; p<subRegion.size(); p++ )
    {
      particleInitialVolume[p] = particleVolume[p];
      particleCrystalHealFlag[p] = 0;
      particleTemperature[p] = m_initialTemperature;
      particleInternalEnergy[p] = 0.0;
      particleKineticEnergy[p] = 0.0;
      particleArtificialViscosity[p] = 0.0;
      particleSPHJacobian[p] = 1.0;
      particleReferencePorosity[p] = 0.0;
      particlePorosity[p] = particleReferencePorosity[p];

      for( int i=0; i<3; i++ )
      {
        particleBodyForce[p][i] = 0;
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
      particleHeatCapacity[p] = 1.0; // CC: TODO Need to get this from constitutive model
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
    m_smallMass = fmin( globalMinMass * 1.0e-12, m_smallMass ); // CC: should the small mass be DBL_MIN instead?

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
      file << "Time, Sxx, Syy, Szz, Syz, Sxz, Sxy, Density, Damage, energy, epxx, epyy, epzz, epyz, epxz, epxy, volume" << std::endl;
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
      GEOS_ERROR_IF( m_boxAverageMin.size() != 3, "Box average min must have 3 elements" );
    }

    if(  m_boxAverageMax.size() == 0  )
    {
      m_boxAverageMax.resize( 3 );
      LvArray::tensorOps::copy< 3 >( m_boxAverageMax, m_xGlobalMax );
    }
    else
    {
      GEOS_ERROR_IF( m_boxAverageMax.size() != 3, "Box average min must have 3 elements" );
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
  nodeManager.getReference< array3d< real64 > >( viewKeyStruct::displacementString() ).resize( numNodes, m_numVelocityFields, 3 );
  nodeManager.getReference< array3d< real64 > >( viewKeyStruct::particleSurfaceNormalString() ).resize( numNodes, m_numVelocityFields, 3 );
  nodeManager.getReference< array2d< real64 > >( viewKeyStruct::massString() ).resize( numNodes, m_numVelocityFields );
  nodeManager.getReference< array2d< real64 > >( viewKeyStruct::materialVolumeString() ).resize( numNodes, m_numVelocityFields );
  nodeManager.getReference< array2d< real64 > >( viewKeyStruct::damageString() ).resize( numNodes, m_numVelocityFields );
  nodeManager.getReference< array2d< real64 > >( viewKeyStruct::maxDamageString() ).resize( numNodes, m_numVelocityFields );
  nodeManager.getReference< array2d< real64 > >( viewKeyStruct::damageGradientString() ).resize( numNodes, 3 );
  nodeManager.getReference< array3d< real64 > >( viewKeyStruct::velocityString() ).resize( numNodes, m_numVelocityFields, 3 );
  nodeManager.getReference< array3d< real64 > >( viewKeyStruct::dVelocityString() ).resize( numNodes, m_numVelocityFields, 3 );
  nodeManager.getReference< array3d< real64 > >( viewKeyStruct::momentumString() ).resize( numNodes, m_numVelocityFields, 3 );
  nodeManager.getReference< array3d< real64 > >( viewKeyStruct::accelerationString() ).resize( numNodes, m_numVelocityFields, 3 );
  nodeManager.getReference< array3d< real64 > >( viewKeyStruct::forceInternalString() ).resize( numNodes, m_numVelocityFields, 3 );
  nodeManager.getReference< array3d< real64 > >( viewKeyStruct::forceExternalString() ).resize( numNodes, m_numVelocityFields, 3 );
  nodeManager.getReference< array3d< real64 > >( viewKeyStruct::forceContactString() ).resize( numNodes, m_numVelocityFields, 3 );
  nodeManager.getReference< array3d< real64 > >( viewKeyStruct::surfaceNormalString() ).resize( numNodes, m_numVelocityFields, 3 );
  nodeManager.getReference< array3d< real64 > >( viewKeyStruct::materialPositionString() ).resize( numNodes, m_numVelocityFields, 3 );

  nodeManager.getReference< array3d< real64 > >( viewKeyStruct::normalStressString() ).resize( numNodes, m_numVelocityFields, 3 );
  nodeManager.getReference< array2d< real64 > >( viewKeyStruct::massWeightedDamageString() ).resize( numNodes, m_numVelocityFields );

  nodeManager.getReference< array3d< real64 > >( viewKeyStruct::vPlusString() ).resize( numNodes, m_numVelocityFields, 3 );
  nodeManager.getReference< array3d< real64 > >( viewKeyStruct::dVPlusString() ).resize( numNodes, m_numVelocityFields, 3 );

  // Precompute the square to save on computation later
  m_maxParticleVelocitySquared = m_maxParticleVelocity * m_maxParticleVelocity;

  initializeConstitutiveModelDependencies( particleManager );

  initializeFrictionCoefficients();

  // Initialize particle wavespeed
  computeWavespeed( particleManager );
}

// Load any particle dependent properties into constitutive model such as strength scale which is determiend per particle not by material model input
void SolidMechanicsMPM::initializeConstitutiveModelDependencies( ParticleManager & particleManager)
{
  GEOS_MARK_FUNCTION;
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
  GEOS_LOG_RANK_IF( m_debugFlag == 1 && m_computeSPHJacobian > 0 , "Compute SPH Jacobian for overlap correction" );
  solverProfilingIf( "Compute SPH Jacobian for overlap correction", m_computeSPHJacobian > 0 );
  //#######################################################################################
  if( m_computeSPHJacobian > 0 )
  {
    computeSPHJacobian( particleManager );
  }


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
  populateMappingArrays( particleManager, nodeManager );//, partition );

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
  particleToGrid( time_n,
                  cycleNumber,
                  particleManager,
                  nodeManager );


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


  //#######################################################################################
  GEOS_LOG_RANK_IF( m_debugFlag == 1, "Compute local averages to generate profile in x-direction: write to file." );
  solverProfiling( "Compute local averages to generate profile in x-direction: write to file." );
  //#######################################################################################
  if( m_computeXProfile == 1 && ( ( m_nextXProfileWriteTime <= time_n ) || ( cycleNumber == 0 ) ) )
  {
    computeXProfile( cycleNumber,
                     time_n,
                     dt,
                     nodeManager,
                     partition );
    m_nextXProfileWriteTime += m_xProfileWriteInterval;
  }

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
  gridToParticle( dt, particleManager, nodeManager, domain, mesh );


  // //#######################################################################################
  // GEOS_LOG_RANK_IF( m_debugFlag == 1, "Compute kinetic energy" );
  // solverProfiling( "Compute kinetic energy" );
  // //#######################################################################################
  // // 1/2*m*v^2 calculation, TODO add micro kinetic energy.
  // computeKineticEnergy( particleManager );


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
  GEOS_LOG_RANK_IF( m_debugFlag == 1 && m_overlapCorrection == 2, "Scale F based on SPH-J" );
  solverProfilingIf( "Scale F based on SPH-J", m_overlapCorrection == 2 );
  //#######################################################################################
  // Scale Jacobian to prevent overdensification if overlap correction type 2 is used.
  if( m_overlapCorrection == 2 )
  {
    overlapCorrection( dt,
                       particleManager );
  }

  //#######################################################################################
  GEOS_LOG_RANK_IF( m_debugFlag == 1, "Update particle geometry (e.g. volume, r-vectors) and density" );
  solverProfiling( "Update particle geometry (e.g. volume, r-vectors) and density" );
  //#######################################################################################
  particleKinematicUpdate( particleManager );

  // //#######################################################################################
  // GEOS_LOG_RANK_IF( m_debugFlag == 1 && m_useArtificialViscosity, "Compute artificial viscosity" );
  // solverProfilingIf( "Compute artificial viscosity", m_useArtificialViscosity );
  // //#######################################################################################
  // if( m_useArtificialViscosity )
  // {
  //   computeArtificialViscosity( particleManager );
  // }

//  //#######################################################################################
//   GEOS_LOG_RANK_IF( m_debugFlag == 1, "Increment internalEnergy with old stress and Fdot" );
//   solverProfiling( "Increment internalEnergy with old stress and Fdot" );
//   //#######################################################################################
//   computeInternalEnergyAndTemperature( dt,
//                                        particleManager );

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

//  //#######################################################################################
//   GEOS_LOG_RANK_IF( m_debugFlag == 1, "Increment internalEnergy with new stress and Fdot" );
//   solverProfiling( "Increment internalEnergy with new stress and Fdot" );
//   //#######################################################################################
//   computeInternalEnergyAndTemperature( dt,
//                                        particleManager );

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
                   particleManager,
                   partition );
  }


  //#######################################################################################
  GEOS_LOG_RANK_IF( m_debugFlag == 1, "Compute particle wavespeeds" );
  solverProfiling( "Compute particle wavespeeds" );
  //#######################################################################################
  computeWavespeed( particleManager );


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
  GEOS_LOG_RANK_IF( m_debugFlag == 1 && MpiWrapper::commSize( MPI_COMM_GEOSX ) > 1, "Particle repartitioning" );
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
  GEOS_LOG_RANK_IF( m_debugFlag == 1 && m_resetDefGradForFullyDamagedParticles == 1, "Set F to scaled value for damaged particles, maintain J" );
  solverProfilingIf( "Set F to scaled value for damaged particles, maintain J",  m_resetDefGradForFullyDamagedParticles == 1 );
  //#######################################################################################
  // Option to set F for fully damaged particles to J^(1/3)*[I], in 3D, so we don't get
  // negative J for super sheared particles with finite precision F Update.  This
  // Should only be used when all materials in the domain have a hypo-elastic
  // deviatoric update.

  if (m_resetDefGradForFullyDamagedParticles == 1)
  {
    resetDeformationGradient( particleManager );
  }

  //#######################################################################################
  GEOS_LOG_RANK_IF( m_debugFlag == 1 && m_plotUnscaledParticles == 1, "Reset CPDI R-vectors for plotting if option is selected to plot unscaled domains" );
  solverProfilingIf( "Reset CPDI R-vectors for plotting if option is selected to plot unscaled domains", m_plotUnscaledParticles == 1 );
  //#######################################################################################
  if (m_plotUnscaledParticles == 1)
  {
    unscaleCPDIVectors( particleManager);
  }

  //#######################################################################################
  GEOS_LOG_RANK_IF( m_debugFlag == 1, "End of explicitStep");
  solverProfiling( "End of explicitStep" );
  //#######################################################################################
  if( m_solverProfiling >= 1 )
  {
    printProfilingResults();
  }

  // // CC: debug
  // GEOS_LOG_RANK( "Domain End: " <<
  //                "xLocalMin=" << m_xLocalMin << ", " <<
  //                "xLocalMax=" << m_xLocalMax << ", " <<
  //                "xGlobalMin=" << m_xGlobalMin << ", " <<
  //                "xGlobalMax=" << m_xGlobalMax << ", " << 
  //                "domainL=" << m_domainL << ", " <<
  //                "domainF=" << m_domainF);

  // // CC: debug
  // GEOS_LOG_RANK( "Step end: LocalMin=" << partition.getLocalMin() << ", " << 
  //                "LocalMax=" << partition.getLocalMax() << ", " << 
  //                "GridMin=" << partition.getGlobalMin() << ", " <<  
  //                "GridMax=" << partition.getGlobalMax() );

  // Return stable time step
  return dtReturn;
}

void SolidMechanicsMPM::triggerEvents( const real64 dt,
                                       const real64 time_n,
                                       ParticleManager & particleManager,
                                       SpatialPartition & partition )
{
  GEOS_MARK_FUNCTION;

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
        AnnealMPMEvent & anneal = dynamicCast< AnnealMPMEvent & >( event );

        particleManager.forParticleRegions< ParticleRegion >( [&]( ParticleRegion & region )
        {
          if( region.getName() == anneal.getTargetRegion() || anneal.getTargetRegion() == "all" )
          {
            auto & targetSubRegions = region.getSubRegions();
            for( int r=0; r < targetSubRegions.size(); ++r)
            {
              ParticleSubRegion & targetSubRegion = dynamicCast< ParticleSubRegion & >( *targetSubRegions[r] );

              // Get constitutive model reference
              string const & solidMaterialName = targetSubRegion.template getReference< string >( viewKeyStruct::solidMaterialNamesString() );
              SolidBase & solidModel = getConstitutiveModel< SolidBase >( targetSubRegion, solidMaterialName );

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
              } );
            }
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

      if( event.getName() == "CrystalHeal" ){
        CrystalHealMPMEvent & crystalHeal = dynamicCast< CrystalHealMPMEvent & >( event );

        int healType = crystalHeal.getHealType();
        particleManager.forParticleRegions< ParticleRegion >( [&]( ParticleRegion & region )
        {
          if( region.getName() == crystalHeal.getTargetRegion() || crystalHeal.getTargetRegion() == "all" )
          {
            // Copy particle data from source sub region to destination sub region
            subGroupMap & targetSubRegions = region.getSubRegions();
            for( int r=0; r < targetSubRegions.size(); ++r)
            {
              ParticleSubRegion & targetSubRegion = dynamicCast< ParticleSubRegion & >( *targetSubRegions[r] );
              
              arrayView1d< real64 > const particleDamage = targetSubRegion.getParticleDamage();
              arrayView1d< int > const particleCrystalHealFlag = targetSubRegion.getField< fields::mpm::particleCrystalHealFlag >();
        
              SortedArrayView< localIndex const > const activeParticleIndices = targetSubRegion.activeParticleIndices();
              string const & solidMaterialName = targetSubRegion.template getReference< string >( viewKeyStruct::solidMaterialNamesString() );
              SolidBase & solidModel = getConstitutiveModel< SolidBase >( targetSubRegion, solidMaterialName );
              if( !crystalHeal.getMarkedParticlesToHeal() )
              {
                arrayView2d< real64 > const particleStress = targetSubRegion.getField< fields::mpm::particleStress >();
                arrayView3d< real64 > const particleRVectors = targetSubRegion.getParticleRVectors();
                arrayView3d< real64 > const particleDeformationGradient = targetSubRegion.getField< fields::mpm::particleDeformationGradient >();
                arrayView1d< real64 > const particleInitialVolume = targetSubRegion.getField< fields::mpm::particleInitialVolume >();

                arrayView2d< real64 const > const particleDamageGradient = targetSubRegion.getField< fields::mpm::particleDamageGradient >();
                
                // CC: TODO Need to add this for Mike
                // arrayView1d< real64 const > const particleReferencePorosity = targetSubRegion.getField< fields::mpm::particleReferencePorosity >();
                
                // Scale constitutive model crackspeed so damage does not evolve while healing
                if( solidModel.hasWrapper( "crackSpeed" ) )
                {
                  real64 & crackSpeed = solidModel.getReference< real64 >( "crackSpeed" );
                  crackSpeed *= 1e-100;
                }

                forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST ( localIndex const pp )
                {
                  localIndex const p = activeParticleIndices[pp];

                  real64 temp[3] = { 0 };
                  // real64 particleDamageGradientNormalized[3] = {0};
                  // LvArray::tensorOps::copy< 3 >( particleDamageGradientNormalized, particleDamageGradient[p] );
                  // norm = LvArray::tensorOps::l2Norm< 3 >( particleDamageGradientNormalized );
                  // // CC: zero stress component in ZZ for plane strain?
                  
                  // // Does the damage gradient need to be normalized since we only care about the sign of the stress magnitude along direction?
                  // real64 normalStress = 0.0
                  // if(norm > 1e-20 ){
                  //   LvArray::tensorOps::normalize< 3 >( particleDamageGradientNormalized ); // CC: debug ( should probably be normalized compute of damage field for particles )
                  //   LvArray::tensorOps::Ri_eq_symAijBj< 3 >( temp, particleStress[p], particleDamageGradientNormalized );
                  //   normalStress = LvArray::tensorOps::AiBi< 3 >( particleDamageGradientNormalized, temp );
                  // }

                  LvArray::tensorOps::Ri_eq_symAijBj< 3 >( temp, particleStress[p], particleDamageGradient[p] );
                  real64 normalStress = LvArray::tensorOps::AiBi< 3 >( particleDamageGradient[p], temp );

                  real64 detF = LvArray::tensorOps::determinant< 3 >( particleDeformationGradient[p] );
                  if( ( healType == 1 || ( healType == 0 && ( normalStress < 0.0 || detF < 1.0 ) ) ) && particleDamage[p] > 0.0 )
                  {
                    particleCrystalHealFlag[p] = 1;

                    if( detF > 1.0 )
                    {
                    // particleReferencePorosity[p] = 1.0 - 1.0 / detF;
                      
                      real64 power = m_planeStrain ? 0.5 : 1.0 / 3.0;
                      real64 scaling = std::pow( detF, power );
                      LvArray::tensorOps::scale< 3, 3 >( particleDeformationGradient[p], 1 / scaling );

                      LvArray::tensorOps::scale< 3 >( particleRVectors[p][0], scaling );
                      LvArray::tensorOps::scale< 3 >( particleRVectors[p][1], scaling );

                      if( !m_planeStrain )
                      {
                        LvArray::tensorOps::scale< 3 >( particleRVectors[p][2], scaling );
                      }

                      particleInitialVolume[p] *= detF;
                    }

                  }
                  else
                  {
                    particleCrystalHealFlag[p] = 0;
                  }
                } );
                crystalHeal.setMarkedParticlesToHeal( 1 );
              }
              else
              {
                // Heal particles by gradually reducing their damage over a user-specified interval.
                // Explicitly force exactly zero damage at the end.
                if( time_n - dt/2 <= eventTime + eventInterval && eventTime + eventInterval < time_n + dt / 2 )
                {
                  forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST ( localIndex const pp )
                  {
                    localIndex const p = activeParticleIndices[pp];
                    if( particleCrystalHealFlag[p] == 1 )
                    {
                      particleDamage[p] = 0.0;
                    }
                  } );
                  event.setIsComplete( 1 );
                }
                else
                {
                  forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST ( localIndex const pp )
                  {
                    localIndex const p = activeParticleIndices[pp];
                    if( particleCrystalHealFlag[p] == 1 )
                    {
                      particleDamage[p] *= fmax( 0.0, 1.0 - dt * 20.0 * ( time_n - eventTime ) / ( eventInterval * eventInterval ) ); // Recursive form
                    }
                  } );
                }
              }
            
              if( event.isComplete() )
              {
                if( solidModel.hasWrapper( "crackSpeed" ) )
                {
                  real64 & crackSpeed = solidModel.getReference< real64 >( "crackSpeed" );
                  crackSpeed *= 1e100;
                }
              }
            }
          }
        } );
      }

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

      if( event.getName() == "BodyForceUpdate" )
      {
        BodyForceUpdateMPMEvent & bodyForceUpdate = dynamicCast< BodyForceUpdateMPMEvent & >( event );

        LvArray::tensorOps::copy< 3 >( m_bodyForce, bodyForceUpdate.getBodyForce() );

        event.setIsComplete( 1 );
      }

      if( event.getName() == "DeformationUpdate" )
      {
        DeformationUpdateMPMEvent & deformationUpdate = dynamicCast< DeformationUpdateMPMEvent & >( event );
        
        int oldPrescribedFTable = m_prescribedFTable;
        int oldPrescribedBoundaryFTable = m_prescribedBoundaryFTable;
        int oldStressControl[3];
        GEOS_UNUSED_VAR( oldPrescribedFTable );
        GEOS_UNUSED_VAR( oldPrescribedBoundaryFTable );
        GEOS_UNUSED_VAR( oldStressControl );
        LvArray::tensorOps::copy< 3 >( oldStressControl, m_stressControl );

        // CC: TODO need to add special handling for turning off and on stress control because the F values will not be known
        m_prescribedFTable = deformationUpdate.getPrescribedFTable();
        m_prescribedBoundaryFTable = deformationUpdate.getPrescribedBoundaryFTable();
        m_stressControl = deformationUpdate.getStressControl();
      }
    }
  } );

}

void SolidMechanicsMPM::performMaterialSwap( ParticleManager & particleManager,
                                             string sourceRegionName,
                                             string destinationRegionName )
{
  GEOS_MARK_FUNCTION;
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
    arrayView1d< real64 > const sourceParticleTemperature = sourceSubRegion.getField< fields::mpm::particleTemperature >();
    arrayView1d< real64 > const sourceParticleInitialVolume = sourceSubRegion.getField< fields::mpm::particleInitialVolume >();
    arrayView1d< real64 > const sourceParticleDensity = sourceSubRegion.getField< fields::mpm::particleDensity >();
    arrayView1d< real64 > const sourceParticleOverlap = sourceSubRegion.getField< fields::mpm::particleOverlap >();
    arrayView1d< int > const sourceParticleSurfaceFlag = sourceSubRegion.getField< fields::mpm::particleSurfaceFlag >();
    arrayView1d< int > const sourceIsBad = sourceSubRegion.getField< fields::mpm::isBad >();
    arrayView1d< int > const sourceParticleCrystalHealFlag = sourceSubRegion.getField< fields::mpm::particleCrystalHealFlag >();
    arrayView1d< real64 > const sourceParticleInternalEnergy = sourceSubRegion.getField< fields::mpm::particleInternalEnergy >();
    arrayView1d< real64 > const sourceParticleArtificialViscosity = sourceSubRegion.getField< fields::mpm::particleArtificialViscosity >();
    arrayView1d< real64 > const sourceParticleHeatCapacity = sourceSubRegion.getField< fields::mpm::particleHeatCapacity >();

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

    arrayView1d< real64 > const sourceParticlePorosity = sourceSubRegion.getField< fields::mpm::particlePorosity >();
    arrayView1d< real64 > const sourceParticleReferencePorosity = sourceSubRegion.getField< fields::mpm::particleReferencePorosity >();

    arrayView1d< real64 > const destinationParticleMass = destinationSubRegion.getField< fields::mpm::particleMass >();
    arrayView1d< real64 > const destinationParticleTemperature = destinationSubRegion.getField< fields::mpm::particleTemperature >();
    arrayView1d< real64 > const destinationParticleInitialVolume = destinationSubRegion.getField< fields::mpm::particleInitialVolume >();
    arrayView1d< real64 > const destinationParticleDensity = destinationSubRegion.getField< fields::mpm::particleDensity >();
    arrayView1d< real64 > const destinationParticleOverlap = destinationSubRegion.getField< fields::mpm::particleOverlap >();
    arrayView1d< int > const destinationParticleSurfaceFlag = destinationSubRegion.getField< fields::mpm::particleSurfaceFlag >();
    arrayView1d< int > const destinationIsBad = destinationSubRegion.getField< fields::mpm::isBad >();
    arrayView1d< int > const destinationParticleCrystalHealFlag = destinationSubRegion.getField< fields::mpm::particleCrystalHealFlag >();
    arrayView1d< real64 > const destinationParticleInternalEnergy = destinationSubRegion.getField< fields::mpm::particleInternalEnergy >();
    arrayView1d< real64 > const destinationParticleArtificialViscosity = destinationSubRegion.getField< fields::mpm::particleArtificialViscosity >();
    arrayView1d< real64 > const destinationParticleHeatCapacity = destinationSubRegion.getField< fields::mpm::particleHeatCapacity >();

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
    
    arrayView1d< real64 > const destinationParticlePorosity = destinationSubRegion.getField< fields::mpm::particlePorosity >();
    arrayView1d< real64 > const destinationParticleReferencePorosity = destinationSubRegion.getField< fields::mpm::particleReferencePorosity >();

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
      destinationParticleTemperature[p] = sourceParticleTemperature[p];
      destinationParticleInitialVolume[p] = sourceParticleInitialVolume[p];
      destinationParticleOverlap[p] = sourceParticleOverlap[p];
      destinationParticleSurfaceFlag[p] = sourceParticleSurfaceFlag[p];
      destinationIsBad[p] = sourceIsBad[p];
      destinationParticleCrystalHealFlag[p] = sourceParticleCrystalHealFlag[p];
      destinationParticleInternalEnergy[p] = sourceParticleInternalEnergy[p];
      destinationParticleArtificialViscosity[p] = sourceParticleArtificialViscosity[p];
      destinationParticleHeatCapacity[p] = sourceParticleHeatCapacity[p];
      destinationParticlePorosity[p] = sourceParticlePorosity[p];
      destinationParticleReferencePorosity[p] = sourceParticleReferencePorosity[p];

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
  GEOS_MARK_FUNCTION;

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
                                                         arrayView3d< real64 > const & dVectorMultiField, //Change in vectorMultiField
                                                         arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const gridPosition,
                                                         Group & nodeSets )
{
  GEOS_MARK_FUNCTION;
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
    forAll< serialPolicy >( numBoundaryNodes, [&, vectorMultiField, dVectorMultiField] GEOS_HOST ( localIndex const gg ) // Probably not a big enough loop to
                                                                                                      // warrant parallelization
      {
        int const g = boundaryNodes[gg];
        dVectorMultiField[g][fieldIndex][dir0] = -vectorMultiField[g][fieldIndex][dir0];
        vectorMultiField[g][fieldIndex][dir0] = 0.0;
      } );

    // Perform field reflection on buffer nodes
    SortedArrayView< localIndex const > const bufferNodes = m_bufferNodes[face].toView();
    int const numBufferNodes = bufferNodes.size();
    forAll< serialPolicy >( numBufferNodes, [&, vectorMultiField, dVectorMultiField, gridPosition] GEOS_HOST ( localIndex const gg ) // Probably not a big
                                                                                                                  // enough loop to warrant
                                                                                                                  // parallelization
      {
        int const g = bufferNodes[gg];
        int ijk[3];
        ijk[dir0] = positiveNormal * ( m_nEl[dir0] - 2 ) + ( 1 - positiveNormal ) * ( 2 );
        ijk[dir1] = std::round( ( gridPosition[g][dir1] - m_xLocalMin[dir1] ) / m_hEl[dir1] );
        ijk[dir2] = std::round( ( gridPosition[g][dir2] - m_xLocalMin[dir2] ) / m_hEl[dir2] );

        localIndex gFrom = m_ijkMap[ijk[0]][ijk[1]][ijk[2]];

        real64 previousVectorMultiField[3] = { 0 };
        previousVectorMultiField[dir0] = vectorMultiField[g][fieldIndex][dir0];
        previousVectorMultiField[dir1] = vectorMultiField[g][fieldIndex][dir1];
        previousVectorMultiField[dir2] = vectorMultiField[g][fieldIndex][dir2];

        vectorMultiField[g][fieldIndex][dir0] = -vectorMultiField[gFrom][fieldIndex][dir0]; // Negate component aligned with surface normal
        vectorMultiField[g][fieldIndex][dir1] =  vectorMultiField[gFrom][fieldIndex][dir1];
        vectorMultiField[g][fieldIndex][dir2] =  vectorMultiField[gFrom][fieldIndex][dir2];

        dVectorMultiField[g][fieldIndex][dir0] += vectorMultiField[g][fieldIndex][dir0] - previousVectorMultiField[dir0];
        dVectorMultiField[g][fieldIndex][dir1] += vectorMultiField[g][fieldIndex][dir1] - previousVectorMultiField[dir1];
        dVectorMultiField[g][fieldIndex][dir2] += vectorMultiField[g][fieldIndex][dir2] - previousVectorMultiField[dir2];
      } );
  }
}

void SolidMechanicsMPM::enforceGridVectorFieldSymmetryBC( arrayView3d< real64 > const & vectorMultiField,
                                                          arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const gridPosition,
                                                          Group & nodeSets )
{
  GEOS_MARK_FUNCTION;
  for( int face=0; face<6; face++ )
  {
    if( m_boundaryConditionTypes[face] == 1 || m_boundaryConditionTypes[face] == 2 )
    {
      array3d< real64 > dVectorMultiField( vectorMultiField.size(0), vectorMultiField.size(1), vectorMultiField.size(2) ); // dumby variable, CC TODO there must be a better solution
      singleFaceVectorFieldSymmetryBC( face, vectorMultiField, dVectorMultiField, gridPosition, nodeSets );
    }
  }
}

void SolidMechanicsMPM::applyEssentialBCs( const real64 dt,
                                           const real64 time_n,
                                           NodeManager & nodeManager )
{
  GEOS_MARK_FUNCTION;

  // Get grid fields
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const gridPosition = nodeManager.referencePosition();
  arrayView1d< int const > const gridGhostRank = nodeManager.ghostRank();
  arrayView2d< real64 > const gridMass = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::massString() );
  arrayView3d< real64 > const gridVelocity = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::velocityString() );
  arrayView3d< real64 > const gridDVelocity = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::dVelocityString() );
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
      // CC: TODO add gridDVelocity update to theses for XPIC
      singleFaceVectorFieldSymmetryBC( face, gridVelocity, gridDVelocity, gridPosition, nodeSets );

      //Dumby paramter, we do not need the change in grid acceleration
      array3d< real64 > gridDAcceleration( gridAcceleration.size(0), gridAcceleration.size(1), gridAcceleration.size(2)); //CC: TODO Probably want to avoid allocating a lot of memory just for the dummy variable
      singleFaceVectorFieldSymmetryBC( face, gridAcceleration, gridDAcceleration, gridPosition, nodeSets );
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
        forAll< serialPolicy >( numBoundaryNodes, [&, gridPosition, gridVelocity, gridDVelocity, gridMass] GEOS_HOST ( localIndex const gg ) // Probably
                                                                                                                              // not a big
                                                                                                                              // enough loop
                                                                                                                              // to warrant
                                                                                                                              // parallelization
          {
            int const g = boundaryNodes[gg];
            real64 prescribedVelocity = m_domainL[dir0] * gridPosition[g][dir0];
            gridDVelocity[g][fieldIndex][dir0] = prescribedVelocity - gridVelocity[g][fieldIndex][dir0]; // CC: TODO double check this, because it overrides the change in velocity that might have been written during enforceContact
            real64 accelerationForBC = gridDVelocity[g][fieldIndex][dir0] / dt; // acceleration needed to satisfy BC
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
        forAll< serialPolicy >( numBufferNodes, [&, gridPosition, gridVelocity, gridDVelocity, gridAcceleration] GEOS_HOST ( localIndex const gg )
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

            //Store previous velocity to compute change in velocity
            real64 gridPreviousVelocity[3] = { 0 };
            gridPreviousVelocity[dir0] = gridVelocity[g][fieldIndex][dir0];
            gridPreviousVelocity[dir1] = gridVelocity[g][fieldIndex][dir1];
            gridPreviousVelocity[dir2] = gridVelocity[g][fieldIndex][dir2];

            // Calculate velocity, Negate component aligned with surface normal and correct for moving boundary
            gridVelocity[g][fieldIndex][dir0] = -gridVelocity[gFrom][fieldIndex][dir0] + 2.0 * gridVelocity[gBoundary][fieldIndex][dir0];
            gridVelocity[g][fieldIndex][dir1] = gridVelocity[gFrom][fieldIndex][dir1];
            gridVelocity[g][fieldIndex][dir2] = gridVelocity[gFrom][fieldIndex][dir2];

            // Compute change in velocity for XPIC calculations
            gridDVelocity[g][fieldIndex][dir0] += gridVelocity[g][fieldIndex][dir0] - gridPreviousVelocity[dir0];
            gridDVelocity[g][fieldIndex][dir1] += gridVelocity[g][fieldIndex][dir1] - gridPreviousVelocity[dir1];
            gridDVelocity[g][fieldIndex][dir2] += gridVelocity[g][fieldIndex][dir2] - gridPreviousVelocity[dir2];

            // Calculate acceleration, Negate component aligned with surface normal and correct for moving boundary
            gridAcceleration[g][fieldIndex][dir0] = -gridAcceleration[gFrom][fieldIndex][dir0] + 2.0 * gridAcceleration[gBoundary][fieldIndex][dir0];
            gridAcceleration[g][fieldIndex][dir1] = gridAcceleration[gFrom][fieldIndex][dir1];
            gridAcceleration[g][fieldIndex][dir2] = gridAcceleration[gFrom][fieldIndex][dir2];
          } );
      }
    }
  }

  // Compute change in grid velocities for XPIC calculations


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
  LvArray::tensorOps::copy< 6 >( m_globalFaceReactions, globalFaceReactions );

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
  GEOS_MARK_FUNCTION;
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
  GEOS_MARK_FUNCTION;
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
  GEOS_MARK_FUNCTION;
  if( m_frictionCoefficientTable.size(0) != 0 )
  {
    GEOS_ERROR_IF( m_frictionCoefficientTable.size(0) != m_frictionCoefficientTable.size(1), "frictionCoefficientTable must be square.");
    GEOS_ERROR_IF( m_frictionCoefficientTable.size(0) != m_numContactGroups, "frictionCoefficientTable must have the same number of rows and columns as the number of contact groups.");

    for(int i = 0; i < m_numContactGroups; i++)
    {
      for(int j = i+1; i < m_numContactGroups; j++)
      {
        GEOS_ERROR_IF( std::abs(m_frictionCoefficientTable[i][j] - m_frictionCoefficientTable[j][i]) > DBL_EPSILON, "Off-diagonal friction coefficients must match" );
      }
    }
    return;
  }
  else
  {
    if( static_cast<int>( m_frictionCoefficient ) == -1 )
    {
      m_frictionCoefficient = 0.0;
    }

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
                                              arrayView2d< real64 const > const & gridMaterialVolume,
                                              arrayView2d< real64 const > const & gridDamage,
                                              arrayView2d< real64 const > const & gridMaxDamage,
                                              arrayView3d< real64 const > const & gridVelocity,
                                              arrayView3d< real64 const > const & gridMomentum,
                                              arrayView3d< real64 const > const & gridSurfaceNormal,
                                              arrayView3d< real64 const > const & gridMaterialPosition,
                                              arrayView3d< real64 > const & gridContactForce )
{
  GEOS_MARK_FUNCTION;
  // Get number of nodes
  int numNodes = gridMass.size( 0 ); // get size from nodeManager instead?

  forAll< serialPolicy >( numNodes, [&, gridMass, gridVelocity, gridMomentum, gridSurfaceNormal, gridMaterialPosition, gridContactForce] GEOS_HOST ( localIndex const g )
    {
      // Initialize gridContactForce[g] to zero. TODO: This shouldn't be necessary?
      // This looks to be zeroed every timestep in initializeGridFields, need to test if removing this breaks anything
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
                        ( gridMass[g][B] > m_smallMass ) && ( LvArray::tensorOps::l2NormSquared< 3 >( gridSurfaceNormal[g][B] ) > 1.0e-16 ); // CC: Should grid surface normal min magnitude be DBL_MIN instead?

            real64 frictionCoefficient = m_frictionCoefficientTable[A % m_numContactGroups][B % m_numContactGroups];

            // CC: debug
            // GEOS_LOG_RANK("Node num: " << g << 
            //               ", A: " << A << 
            //               ", B: " << B << 
            //               ", gA x: " <<  gridMaterialPosition[g][A] <<
            //               ", gA mass: " << gridMass[g][A] <<
            //               ", gA: " << gridSurfaceNormal[g][A] << 
            //               ", gA norm: " << LvArray::tensorOps::l2NormSquared< 3 >( gridSurfaceNormal[g][A] )  << 
            //               ", gB x: " <<  gridMaterialPosition[g][B] <<
            //               ", gB mass: " << gridMass[g][B] <<
            //               ", gB: " << gridSurfaceNormal[g][B] << 
            //               ", gB norm: " << LvArray::tensorOps::l2NormSquared< 3 >( gridSurfaceNormal[g][B] ) );

          if( active )
          {
            // Evaluate the separability criterion for the contact pair.
            int separable = evaluateSeparabilityCriterion( A,
                                                           B,
                                                           gridDamage[g][A],
                                                           gridDamage[g][B],
                                                           gridMaxDamage[g][A],
                                                           gridMaxDamage[g][B],
                                                           gridMaterialPosition[g][A],
                                                           gridMaterialPosition[g][B] );

            computePairwiseNodalContactForce( separable,
                                              dt,
                                              frictionCoefficient,
                                              gridMass[g][A],
                                              gridMass[g][B],
                                              gridMaterialVolume[g][A],
                                              gridMaterialVolume[g][B],
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
                                                          real64 const & VA,
                                                          real64 const & VB,
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
  GEOS_MARK_FUNCTION;
  // Total mass for the contact pair.
  real64 mAB = mA + mB;

  // Outward normal of field A with respect to field B.
  real64 nAB[3];
  
  // contact normal averaging: 0: simple, 1: mass weighted, 2: follow larger mass (useful for platens)
  if( m_contactNormalType == 0 ){
    // no weighting.  this usually isn't very good.
    for( int i=0; i<3; i++ )
    {
      nAB[i] = nA[i] - nB[i];
    }
  }
  else if ( m_contactNormalType == 1 )
  {
    // Mass-weighted average of the field normals
    for( int i=0; i<3; i++ )
    {
      nAB[i] = nA[i] * mA - nB[i] * mB;
    }
  }
  else if ( m_contactNormalType == 2 )
  {
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
  }
  else if ( m_contactNormalType == 3 )
  {
	  // If density is similar use mass weighted average
	  if ( fabs( ( mA / VA ) - ( mB / VB ) ) < 0.1 * ( mA / VA ) )
	  {
      for( int i=0; i<3; i++ )
      {
        nAB[i] = nA[i] * mA - nB[i] * mB;
      }
	  }
	  else
	  { 
      // if one field is significantly more dense,
		  // Use the surface normal for whichever field has more density.
		  // This should be good for corners against flat surfaces.
      if( mA / VA > mB / VB )
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
	  }
  }
  else
  { 
    GEOS_ERROR( "Unknown contact normal type specified!" );
  }

  // Vector pointing from CoM of A field to CoM of B field, stretched by grid spacing
  // for( int i=0; i<3; i++ )
  // {
  //   nAB[i] = (xB[i] - xA[i]) / m_hEl[i];
  // }

  // Normalize the effective surface normal
  if( m_planeStrain == 1 )
  {
    nAB[2] = 0.0;
  }

  real64 norm = sqrt( nAB[0] * nAB[0] + nAB[1] * nAB[1] + nAB[2] * nAB[2] );
 
  // CC: Is this the best solution?
  // Need to prevent crashing from randomly defined normals that might arise when running with damage field gradient partitioning
  // If the normal is small, then assume no separation force and return
  if ( norm < 1e-20 )
  {
      fA[0] = 0.0;
      fA[1] = 0.0;
      fA[2] = 0.0;
      fB[0] = 0.0;
      fB[1] = 0.0;
      fB[2] = 0.0;
      return;
  }

  nAB[0] /= norm;
  nAB[1] /= norm;
  nAB[2] /= norm;


  // Calculate the contact gap between the fields
  real64 gap0;
  if( m_planeStrain == 1 )
  {
    // gap0 = ( m_hEl[0]*m_hEl[1] ) / sqrt( m_hEl[1]*m_hEl[1]*nAB[0]*nAB[0] + m_hEl[0]*m_hEl[0]*nAB[1]*nAB[1] );
    gap0 = 1/sqrt( std::pow(nAB[0]/m_hEl[0],2) + std::pow(nAB[1]/m_hEl[1],2) );
  }
  else
  {
    // CC: debug
    // currently failes with polymer model (suspect compliant materials) when z direction boundary type is 2
    //This needs to be corrected
    // gap0 = (m_hEl[0]*m_hEl[1]*m_hEl[2]) /
    //        sqrt( m_hEl[2]*m_hEl[2]*nAB[2]*nAB[2]*( m_hEl[1]*m_hEl[1]*nAB[0]*nAB[0] + m_hEl[0]*m_hEl[0]*nAB[1]*nAB[1] ) + m_hEl[0]*m_hEl[0]*m_hEl[1]*m_hEl[1]*( nAB[0]*nAB[0] + nAB[1]*nAB[1] ) );
 
    // Ellipse solution
    // Gives element sizes in each direction and reasonable approximations in others
    gap0 = 1/sqrt(std::pow(nAB[0]/m_hEl[0],2) + std::pow(nAB[1]/m_hEl[1],2) + std::pow(nAB[2]/m_hEl[2],2));

  }

  // TODO: A fudge factor of 0.67 on gap0 makes diagonal surfaces (wrt grid) close better I think, but this is more general
  real64 gap = (xB[0] - xA[0]) * nAB[0] + (xB[1] - xA[1]) * nAB[1] + (xB[2] - xA[2]) * nAB[2] - gap0;

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
    real64 contact = 0.0;
    real64 temp[3] = { 0 };
    LvArray::tensorOps::copy< 3 >( temp, vA );
    LvArray::tensorOps::subtract< 3 >( temp, vAB );
    real64 test = LvArray::tensorOps::AiBi< 3 >( temp, nAB);
    
    // real64 test = (vA[0] - vAB[0]) * nAB[0] + (vA[1] - vAB[1]) * nAB[1] + (vA[2] - vAB[2]) * nAB[2];
    if( m_contactGapCorrection == 1 ) // TODO: Implement soft/ramped contact option?
    {
      contact = test > 0.0 && gap < 0.0 ? 1.0 : 0.0;
    }
    else if ( m_contactGapCorrection == 2 )
    {

		  if ( test > 0 )
		  {
        // realT gap0 = normalSpacing*cellSpacing;
			  if (gap <= 0.0)
			  {
				  contact = 1.0;
			  }
			  else if (gap < gap0)
			  {
				  contact = 1.0 - gap/gap0;
			  }
		  }
    }
    else
    {
      contact = test > 0.0 ? 1.0 : 0.0;
    }
  
    // Modify normal contact force
    fnor *= contact;

	  // Additional compressive normal force added to correct overlap, this is added only to the normal force and
	  // doesn't affect the shear calculation (which in FEM cuts down on spurious stress oscillations).
	  real64 fgap = 0.0;
	  if( m_overlapCorrection == 1 )
	  {
      real64 cellSpacing[3];
      LvArray::tensorOps::copy< 3 >( cellSpacing, m_hEl);

      real64 normalSpacing[3];
      normalSpacing[0] = fabs( nAB[0] );
      normalSpacing[1] = fabs( nAB[1] );
      normalSpacing[2] = fabs( nAB[2] );

		  real64 cellVolume = m_hEl[0] * m_hEl[1] * m_hEl[2];
		  real64 cellLength = LvArray::tensorOps::AiBi< 3 >( normalSpacing, cellSpacing );
		  real64 cellArea = cellVolume / cellLength;
		  real64 overlap = m_planeStrain ? ( VA + VB - 0.5*cellVolume ) / cellArea : ( VA + VB - cellVolume ) / cellArea;
		  if( overlap > 0.0 )
		  {
			  fgap = -1.0 * overlap * fmin( mA, mB ) / ( dt * dt ); // This could be -0.5*overlap*mAB/dt**2, but that might be less stable if one mass tiny
		  }
	  }

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

    real64 ftan = std::min( frictionCoefficient * std::abs( fnor ), ftanMag ); // This goes to zero when contact=0 due to the std::min
    dfA[0] = ( fnor + fgap ) * nAB[0] + ftan * sAB[0];
    dfA[1] = ( fnor + fgap ) * nAB[1] + ftan * sAB[1];
    dfA[2] = ( fnor + fgap ) * nAB[2] + ftan * sAB[2];
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
  GEOS_MARK_FUNCTION;
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
  GEOS_MARK_FUNCTION;

 // TODO: Find a way to automatically loop over these fields.
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
  GEOS_MARK_FUNCTION;

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
  GEOS_MARK_FUNCTION;
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
  GEOS_MARK_FUNCTION;

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
  GEOS_MARK_FUNCTION;
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
  GEOS_MARK_FUNCTION;
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
  GEOS_MARK_FUNCTION;
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
  GEOS_MARK_FUNCTION;
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
  GEOS_MARK_FUNCTION;
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
  GEOS_MARK_FUNCTION;
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
  GEOS_MARK_FUNCTION;

  // CC: TODO Why is this done instead of getting these inside the particle region loop like elsewhere?
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

        // CC: should the particle damage gradient be normalized?
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
  GEOS_MARK_FUNCTION;

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
  GEOS_MARK_FUNCTION;

  real64 dtSub = dt / m_FSubcycles;

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

      // UPDATE DEFORMATION GRADIENT L = Fdot * Finv -> F_new = F_old + Fdot*dt = F_old + L*F_old*dt
      real64 oldDeformationGradient[3][3];
      LvArray::tensorOps::copy< 3, 3 >( oldDeformationGradient, particleDeformationGradient[p] );

      // Integrate F with L
      real64 previousDeformationGradient[3][3];
      for( int iter = 0 ; iter < m_FSubcycles ; iter++ )
      {
        LvArray::tensorOps::copy< 3, 3 >( previousDeformationGradient, particleDeformationGradient[p] );

        real64 temp[3][3] = { { 0 } };
        LvArray::tensorOps::copy< 3, 3 >( temp, previousDeformationGradient );

        real64 tempFDot[3][3] = { { 0 } };
        LvArray::tensorOps::Rij_eq_AikBkj< 3, 3, 3 >( tempFDot, particleVelocityGradient[p], temp );
        LvArray::tensorOps::scale< 3, 3 >( tempFDot, dtSub);
        LvArray::tensorOps::add< 3, 3 >( particleDeformationGradient[p], tempFDot );
      }

      if (m_exactJIntegration)
      { // This will improve accuracy with nearly incompressible materials
        // Exact integration of J with tr(L)
        real64 Jold = LvArray::tensorOps::determinant< 3 >( oldDeformationGradient );
        real64 Jnew = Jold * exp( LvArray::tensorOps::trace< 3 >( particleVelocityGradient[p] ) * dt );

        // Modify F to reduce volumetric integration error
        real64 detF = LvArray::tensorOps::determinant< 3 >( particleDeformationGradient[p] );
        real64 scale = 1.0;
        if( detF > 0.0 && Jnew >= 0.0 )
        {
          real64 power = m_planeStrain ? 0.5 : 1.0/3.0;
          scale = pow( Jnew / detF, power );
        }
        LvArray::tensorOps::scale< 3, 3 >( particleDeformationGradient[p], scale );
      }

      real64 temp[3][3];
      LvArray::tensorOps::copy< 3, 3 >( temp, particleDeformationGradient[p] );
      
      real64 minusOldDeformationGradient[3][3];
      LvArray::tensorOps::scaledCopy< 3, 3 >( minusOldDeformationGradient, oldDeformationGradient, -1.0);
      LvArray::tensorOps::add< 3, 3 >( temp, minusOldDeformationGradient );
      LvArray::tensorOps::scale< 3, 3 >( temp, 1 / dt );
      LvArray::tensorOps::copy< 3, 3 >( particleFDot[p], temp );

    } );
  } );
}

void SolidMechanicsMPM::updateConstitutiveModelDependencies( ParticleManager & particleManager )
{
  GEOS_MARK_FUNCTION;

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

    if(  solidModel.hasWrapper( "temperature" ) )
    {
      arrayView1d< real64 > const particleTemperature = subRegion.getField< fields::mpm::particleTemperature >();
      arrayView1d< real64 > const constitutiveTemperature = solidModel.getReference< array1d< real64 > >( "temperature" );
      forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
      {
        localIndex const p = activeParticleIndices[pp];
        constitutiveTemperature[p] = particleTemperature[p]; 
      } );
    }
  } );
}

void SolidMechanicsMPM::updateStress( real64 dt,
                                      ParticleManager & particleManager )
{
  GEOS_MARK_FUNCTION;

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
  GEOS_MARK_FUNCTION;

  // Update particle volume and density
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    // Get particle fields
    arrayView1d< globalIndex const > const particleID = subRegion.getParticleID();
    arrayView1d< real64 > const particleVolume = subRegion.getParticleVolume();
    arrayView2d< real64 const > const particleVelocity = subRegion.getParticleVelocity();
    arrayView1d< int > const isBad = subRegion.getField< fields::mpm::isBad >();
    arrayView1d< real64 const > const particleInitialVolume = subRegion.getField< fields::mpm::particleInitialVolume >();
    arrayView1d< real64 > const particleDensity = subRegion.getField< fields::mpm::particleDensity >();
    arrayView1d< real64 const > const particleMass = subRegion.getField< fields::mpm::particleMass >();
    arrayView3d< real64 > const particleRVectors = subRegion.getParticleRVectors();
    arrayView3d< real64 const > const particleInitialRVectors = subRegion.getField< fields::mpm::particleInitialRVectors >();
    arrayView3d< real64 const > const particleDeformationGradient = subRegion.getField< fields::mpm::particleDeformationGradient >();
    arrayView2d< real64 const > const particleInitialMaterialDirection = subRegion.getParticleInitialMaterialDirection();
    arrayView2d< real64 > const particleMaterialDirection = subRegion.getParticleMaterialDirection();
    arrayView2d< real64 const > const particleInitialSurfaceNormal = subRegion.getParticleInitialSurfaceNormal();
    arrayView2d< real64 > const particleSurfaceNormal = subRegion.getParticleSurfaceNormal();

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

      bool flaggedForDeletion = false; 
      if( detF <= m_minParticleJacobian || detF >= m_maxParticleJacobian )
      {
        // printf( "Flagging particle with unreasonable Jacobian (J<%.2f or J>%.2f) for deletion! Global particle ID: %lld", m_minParticleJacobian, m_maxParticleJacobian, particleID[p] );
        flaggedForDeletion = true;
      }

      real64 particleSpeedSquared = LvArray::tensorOps::l2NormSquared< 3 >( particleVelocity[p] );
      if(particleSpeedSquared > m_maxParticleVelocitySquared )
      {
        // printf( "Flagging particle with unreasonable velocity (v > %.2f) for deletion! Global particle ID: %lld", m_maxParticleVelocity, particleID[p] );
        flaggedForDeletion = true;
      }

      if( !flaggedForDeletion )
      {
        particleVolume[p] = particleInitialVolume[p] * detF;
        particleDensity[p] = particleMass[p] / particleVolume[p];

        real64 deformationGradient[3][3];
        LvArray::tensorOps::copy< 3, 3 >(deformationGradient, particleDeformationGradient[p]);

        // Update particle surface normals
        real64 deformationGradientCofactor[3][3];
        cofactor(deformationGradient, deformationGradientCofactor);
        
        real64 surfaceNormal[3];
        LvArray::tensorOps::Ri_eq_AijBj< 3, 3 >( surfaceNormal, deformationGradient, particleInitialSurfaceNormal[p] );
        LvArray::tensorOps::copy< 3 >( particleSurfaceNormal[p], surfaceNormal );

        real64 materialDirection[3];
        if( isFiber )
        {
          LvArray::tensorOps::Ri_eq_AijBj< 3, 3 >( materialDirection, deformationGradient, particleInitialMaterialDirection[p] );
          LvArray::tensorOps::copy< 3 >(particleMaterialDirection[p], materialDirection);
        } 
        else
        {
          LvArray::tensorOps::Ri_eq_AijBj< 3, 3 >( materialDirection, deformationGradientCofactor, particleInitialMaterialDirection[p] );
          LvArray::tensorOps::copy< 3 >(particleMaterialDirection[p], materialDirection);
        } 
      } 
      else
      {
        isBad[p] = 1; // TODO: Switch to a 'markBad' function which changes isBad and removes this particle from activeParticleIndices. Also
                      // change isBad to deletionFlag.
        particleVolume[p] = particleInitialVolume[p];
        particleDensity[p] = particleMass[p] / particleInitialVolume[p];
        LvArray::tensorOps::copy< 3, 3 >( particleRVectors[p], particleInitialRVectors[p] );
      }
    } );
  
    for( int pp = 0; pp < activeParticleIndices.size(); pp++ )
    {
      localIndex const p = activeParticleIndices[pp];
      if( isBad[p] == 1 ){
        real64 detF = LvArray::tensorOps::determinant< 3 >( particleDeformationGradient[p] );
        GEOS_LOG( "Global particle ID " << particleID[p] << " flagged for deletion due to unreasonable jacobian or velocity! J=" << detF << ", v={" << particleVelocity[p][0] << ", " << particleVelocity[p][1] << ", " << particleVelocity[p][2] << "}"  );
      }
    }
  } );

  // Compute particles R vectors
  computeRVectors( particleManager );
}

void SolidMechanicsMPM::computeAndWriteBoxAverage( const real64 dt,
                                                   const real64 time_n,
                                                   ParticleManager & particleManager )
{
  GEOS_MARK_FUNCTION;

  real64 boxPlasticStrain[6] = { 0.0 };
  real64 boxStress[6] = { 0.0 }; // we sum stress * volume in particles, additive sync, then divide by box volume.
  real64 boxMass = 0.0; // we sum particle mass, additive sync, then divide by box volume
  real64 boxParticleInitialVolume = 0.0;
  real64 boxDamage = 0.0; // we sum damage * initial volume, additive sync, then divide by total initial volume in box
  real64 boxEnergy = 0.0;
  real64 boxMatVolume = 0.0; // Sum volume of all particles

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
    arrayView1d< real64 const > const particleInternalEnergy = subRegion.getField< fields::mpm::particleInternalEnergy >();

    // Accumulate values
    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
    forAll< serialPolicy >( activeParticleIndices.size(), [=, &boxMass, &boxParticleInitialVolume, &boxStress, &boxPlasticStrain, &boxDamage, &boxMatVolume, &boxEnergy] GEOS_HOST ( localIndex const pp ) // This
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
        boxMatVolume += particleVolume[p];
        boxParticleInitialVolume += particleInitialVolume[p];
        boxEnergy += particleInternalEnergy[p] * particleVolume[p];
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
  real64 boxSums[17];
  boxSums[0] = boxStress[0];       // sig_xx * volume
  boxSums[1] = boxStress[1];       // sig_yy * volume
  boxSums[2] = boxStress[2];       // sig_zz * volume
  boxSums[3] = boxStress[3];       // sig_yz * volume
  boxSums[4] = boxStress[4];       // sig_xz * volume
  boxSums[5] = boxStress[5];       // sig_xy * volume
  boxSums[6] = boxMass;            // total mass in box
  boxSums[7] = boxParticleInitialVolume > 0.0 ? boxParticleInitialVolume : 1.0;  // total particle initial volume in box; prevent div0 error
  boxSums[8] = boxDamage;          // damage * volume
  boxSums[9] = boxEnergy;          // energy * volume
  boxSums[10] = boxPlasticStrain[0]; // plasticStrain_xx
  boxSums[11] = boxPlasticStrain[1]; // plasticStrain_yy
  boxSums[12] = boxPlasticStrain[2]; // plasticStrain_zz
  boxSums[13] = boxPlasticStrain[3]; // plasticStrain_yz
  boxSums[14] = boxPlasticStrain[4]; // plasticStrain_xz
  boxSums[15] = boxPlasticStrain[5]; // plasticStrain_xy
  boxSums[16] = boxMatVolume;
      
  // Do an MPI sync to total these values and write from proc0 to a file.  Also compute global F
  // so file is directly plottable in excel as CSV or something.
  for( localIndex i = 0; i < 17; i++ )
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

     // Write to file
    std::ofstream file;
    file.open( "boxAverageHistory.csv", std::ios::out | std::ios::app );
    if( file.fail() )
    {
      throw std::ios_base::failure( std::strerror( errno ) );
    }
    //make sure write fails with exception if something is wrong
    file.exceptions( file.exceptions() | std::ios::failbit | std::ifstream::badbit );
    // time | sig_xx | sig_yy | sig_zz | sig_xy | sig_yz | sig_zx | density | damage | energy | epxx | epyy | epzz | epyz | epxz | epxy | total particle volume
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
         << ", "
         << boxSums[15] / boxVolume
         << ", "
         << boxSums[16]
         << std::endl;
    file.close();
  }
}


void SolidMechanicsMPM::computeBoxMetrics( ParticleManager & particleManager,
                                           arrayView1d< real64 > currentStress,
                                           real64 & boxMaterialVolume )
{
  real64 boxStress[3] = { 0.0 }; // we sum stress * volume in particles, additive sync, then divide by box volume.
  boxMaterialVolume = 0.0;

  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    // Get fields
    arrayView1d< real64 > const particleVolume = subRegion.getParticleVolume();
    arrayView2d< real64 > const particleStress = subRegion.getField< fields::mpm::particleStress >();

    // Accumulate values
    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();

    // CC: TODO parallelize via reduction
    forAll< serialPolicy >( activeParticleIndices.size(), [=, &boxStress, &boxMaterialVolume] GEOS_HOST ( localIndex const pp )
    {
      localIndex const p = activeParticleIndices[pp];

      boxMaterialVolume += particleVolume[p];

      // Only need diagonal components
      for( int i=0; i < 3; i++ )
      {
        boxStress[i] += particleStress[p][i] * particleVolume[p]; // volume weighted average, will normalize later.
      }
    } );
  } );

  // Additive sync: sxx, syy, szz, sxy, syz, sxz, mass, particle volume, damage
  real64 boxSums[4];
  boxSums[0] = boxStress[0];       // sig_xx * volume
  boxSums[1] = boxStress[1];       // sig_yy * volume
  boxSums[2] = boxStress[2];       // sig_zz * volume
  boxSums[3] = boxMaterialVolume;

  // Do an MPI sync to total these values and write from proc0 to a file.  Also compute global F
  // so file is directly plottable in excel as CSV or something.
  for( localIndex i = 0; i < 4; i++ )
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

  boxMaterialVolume = boxSums[3];
}

void SolidMechanicsMPM::stressControl( real64 dt,
                                       ParticleManager & particleManager,
                                       SpatialPartition & GEOS_UNUSED_PARAM( partition ) )
{
  GEOS_MARK_FUNCTION;

  // arrayView1d< int const > const periodic = partition.getPeriodic();

  real64 targetStress[3] = {0};
  LvArray::tensorOps::copy< 3 >(targetStress, m_domainStress);

  array1d< real64 > currentStress;
  currentStress.resize( 3 );
  LvArray::tensorOps::fill< 3 >( currentStress, 0.0 );
  // if( periodic[0] == 1 || periodic[1] == 1 || periodic[2] == 1 )
  // {
  real64 boxMaterialVolume;
  computeBoxMetrics( particleManager,
                     currentStress,
                     boxMaterialVolume );
  // }

  // CC: TODO Still use box stress but enfore Lmax for each direction

  // // Non periodic directions should use boundary reactions instead of box stresses
  // for( int i=0; i<m_numDims; i++)
  // {
  //   if( !periodic[i] )
  //   {
  //     real64 area = 1;
  //     for( int j=0; j< m_numDims; j++)
  //     {
  //       if( j == i )
  //       {
  //         continue;
  //       }
  //       area *= m_domainExtent[j];
  //     }
  //     // x-, x+, y-, y+, z-, z+
  //     currentStress[i] = ( m_globalFaceReactions[2*i + 1] * m_xGlobalMax[i] + m_globalFaceReactions[ 2*i ] * m_xGlobalMin[i] ) / ( area * m_domainExtent[i] );
  //   }
  // }

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

    if( constitutiveModelName == "ElasticTransverseIsotropic" || constitutiveModelName == "ElasticTransverseIsotropicPressureDependent" ){
      const ElasticTransverseIsotropic & elasticTransverseIsotropic = dynamic_cast< const ElasticTransverseIsotropic & >( solidModel );
      bulkModulus = elasticTransverseIsotropic.effectiveBulkModulus();
    }

    if( constitutiveModelName == "Graphite" )
    {
      const Graphite & graphite = dynamic_cast< const Graphite & >( solidModel );
      bulkModulus = graphite.effectiveBulkModulus();
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

  // Could some numerical artifact ever produce a negative relative density? Do I need a check for that or clip is to 0
  real64 relativeDensity = fmin( 1.0, boxMaterialVolume / ( m_domainExtent[0] * m_domainExtent[1] * m_domainExtent[2] ) );

	// This will drive the response towards the desired stress but may be
	// unstable.
  // CC: TODO Use 1- domain porosity to make uncompacted region more stable
	real64 stressControlKp = relativeDensity * m_stressControlKp / ( maximumBulkModulus * dt );
	real64 stressControlKd = relativeDensity * m_stressControlKd / ( maximumBulkModulus );
	real64 stressControlKi = relativeDensity * m_stressControlKi / ( maximumBulkModulus * dt * dt );

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
	// boundaryV = strainRate*(m_xGlobalMax[0] - m_xGlobalMin[0]) = CFL*m_hEl[0]/dt
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
  // GEOS_LOG_RANK_0( "Current stress: " << currentStress << "\n" << 
  //                  "Domain stress: " << m_domainStress << "\n" << 
  //                  "Target stress: " << targetStress[0] << ", " << targetStress[1] << ", " << targetStress[2] << "\n" << 
  //                  "Box Material Volume: " << boxMaterialVolume << "\n" <<
  //                  "Relative density: " << relativeDensity << "\n" << 
  //                  "Gains (kp, ki, kd): " << m_stressControlKp << ", " << m_stressControlKi << ", " << m_stressControlKd << "\n" << 
  //                  "error: " << error[0] << ", " << error[1] << ", " << error[2] << "\n" <<
  //                  "dedt: " << dedt[0] << ", " << dedt[1] << ", " << dedt[2] << "\n" <<  
  //                  "P: " << stressControlPTerm[0] << ", " << stressControlPTerm[1] << ", " << stressControlPTerm[2] << "\n" << 
  //                  "I: " << stressControlITerm[0] << ", " << stressControlITerm[1] << ", " << stressControlITerm[2] << "\n" << 
  //                  "D: " << stressControlDTerm[0] << ", " << stressControlDTerm[1] << ", " << stressControlDTerm[2] << "\n" << 
  //                  "Lnew (x,y,z): " << L_new[0] << ", " << L_new[1] << ", " << L_new[2] << "\n" << 
  //                  "Lmin (x,y,z): " << Lxx_min << ", " << Lyy_min << ", " << Lzz_min << "\n" << 
  //                  "Lmax (x,y,z): " << Lxx_max << ", " << Lyy_max << ", " << Lzz_max << "\n" << 
  //                  "Domain L: " << m_domainL << "\n" << 
  //                  "Domain F: " << m_domainF );
}


void SolidMechanicsMPM::applySuperimposedVelocityGradient( const real64 dt, 
                                                           ParticleManager & particleManager,
                                                           SpatialPartition & partition )
{
  GEOS_MARK_FUNCTION;

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
  GEOS_MARK_FUNCTION;

  int const numNodes = nodeManager.size();
  arrayView3d< real64 > const gridDisplacement = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::displacementString() );
  arrayView3d< real64 > const gridParticleSurfaceNormal = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::particleSurfaceNormalString() );
  arrayView2d< real64 > const gridMass = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::massString() );
  arrayView2d< real64 > const gridMaterialVolume = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::materialVolumeString() );
  arrayView2d< real64 > const gridDamage = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::damageString() );
  arrayView2d< real64 > const gridMaxDamage = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::maxDamageString() );
  arrayView2d< real64 > const gridDamageGradient = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::damageGradientString() );
  arrayView3d< real64 > const gridVelocity = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::velocityString() );
  arrayView3d< real64 > const gridDVelocity = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::dVelocityString() );
  arrayView3d< real64 > const gridMomentum = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::momentumString() );
  arrayView3d< real64 > const gridAcceleration = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::accelerationString() );
  arrayView3d< real64 > const gridInternalForce = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::forceInternalString() );
  arrayView3d< real64 > const gridExternalForce = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::forceExternalString() );
  arrayView3d< real64 > const gridContactForce = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::forceContactString() );
  arrayView3d< real64 > const gridSurfaceNormal = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::surfaceNormalString() );
  arrayView3d< real64 > const gridMaterialPosition = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::materialPositionString() );

  arrayView3d< real64 > const gridNormalStress = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::normalStressString() );
  arrayView2d< real64 > const gridMassWeightedDamage = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::massWeightedDamageString() );

  forAll< serialPolicy >( numNodes, [=] GEOS_HOST ( localIndex const g ) // Switch to .zero()?
    {
      for( int i = 0; i < 3; i++ )
      {
        gridDamageGradient[g][i] = 0.0;
      }
      for( int fieldIndex = 0; fieldIndex < m_numVelocityFields; fieldIndex++ )
      {
        gridMass[g][fieldIndex] = 0.0;
        gridMaterialVolume[g][fieldIndex] = 0.0;
        gridDamage[g][fieldIndex] = 0.0;
        gridMaxDamage[g][fieldIndex] = 0.0;
        gridMassWeightedDamage[g][fieldIndex] = 0.0;
        for( int i = 0; i < 3; i++ )
        {
          gridParticleSurfaceNormal[g][fieldIndex][i] = 0.0;
          gridDisplacement[g][fieldIndex][i] = 0.0;
          gridVelocity[g][fieldIndex][i] = 0.0;
          gridDVelocity[g][fieldIndex][i] = 0.0;
          gridMomentum[g][fieldIndex][i] = 0.0;
          gridAcceleration[g][fieldIndex][i] = 0.0;
          gridInternalForce[g][fieldIndex][i] = 0.0;
          gridExternalForce[g][fieldIndex][i] = 0.0;
          gridContactForce[g][fieldIndex][i] = 0.0;
          gridSurfaceNormal[g][fieldIndex][i] = 0.0;
          gridMaterialPosition[g][fieldIndex][i] = 0.0;
          gridNormalStress[g][fieldIndex][i] = 0.0;
        }
      }
    } );
}

void SolidMechanicsMPM::boundaryConditionUpdate( real64 dt, real64 time_n )
{
  GEOS_MARK_FUNCTION;

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


void SolidMechanicsMPM::particleToGrid( real64 const time_n,
                                        integer const cycleNumber,
                                        ParticleManager & particleManager,
                                        NodeManager & nodeManager )
{
  GEOS_MARK_FUNCTION;

  localIndex subRegionIndex = 0;
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    string const & materialName = subRegion.template getReference< string >( viewKeyStruct::solidMaterialNamesString() );
    // Particle fields
    // ParticleType particleType = subRegion.getParticleType(); // CC: unused?
    arrayView2d< real64 const > const particlePosition = subRegion.getParticleCenter();
    arrayView2d< real64 const > const particleVelocity = subRegion.getParticleVelocity();
    arrayView1d< real64 const > const particleMass = subRegion.getField< fields::mpm::particleMass >();
    arrayView1d< real64 const > const particleVolume = subRegion.getParticleVolume();
    arrayView1d< int const > const particleGroup = subRegion.getParticleGroup();
    arrayView1d< int const > const particleSurfaceFlag = subRegion.getParticleSurfaceFlag();
    arrayView1d< real64 const > const particleArtificialViscosity = subRegion.getField< fields::mpm::particleArtificialViscosity >();
    arrayView2d< real64 const > const particleStress = subRegion.getField< fields::mpm::particleStress >();
    arrayView2d< real64 const > const particleBodyForce = subRegion.getField< fields::mpm::particleBodyForce >();
    arrayView2d< real64 const > const particleDamageGradient = subRegion.getField< fields::mpm::particleDamageGradient >();
    arrayView1d< real64 const > const particleDamage = subRegion.getParticleDamage();

    // Grid fields
    arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const gridPosition = nodeManager.referencePosition();
    arrayView2d< real64 > const gridMass = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::massString() );
    arrayView2d< real64 > const gridMaterialVolume = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::materialVolumeString() );
    arrayView2d< real64 const > const gridDamageGradient = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::damageGradientString() );
    arrayView3d< real64 > const gridMomentum = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::momentumString() );
    arrayView3d< real64 > const gridInternalForce = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::forceInternalString() );
    arrayView3d< real64 > const gridExternalForce = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::forceExternalString() );
    arrayView3d< real64 > const gridMaterialPosition = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::materialPositionString() );
    arrayView2d< real64 > const gridDamage = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::damageString() );
    arrayView2d< real64 > const gridMaxDamage = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::maxDamageString() );
    arrayView2d< real64 > const gridMassWeightedDamage = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::massWeightedDamageString() );
    arrayView3d< real64 > const gridNormalStress = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::normalStressString() );

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

        if( m_computeXProfile == 1 && ( ( m_nextXProfileWriteTime <= time_n ) || ( cycleNumber == 0 ) ) )
        {
          for( int i = 0 ; i < 3 ; i++ )
          { // map the volume-weighted normal stress to nodes for profile output.
            gridNormalStress[mappedNode][fieldIndex][i] += shapeFunctionValues[pp][g] * particleStress[p][i] * particleVolume[p];
          }
          gridMassWeightedDamage[mappedNode][fieldIndex] += shapeFunctionValues[pp][g] * particleDamage[p] * particleMass[p];
        }

        gridMass[mappedNode][fieldIndex] += particleMass[p] * shapeFunctionValues[pp][g];
        gridMaterialVolume[mappedNode][fieldIndex] += shapeFunctionValues[pp][g] * particleVolume[p]; // particleMass[p] / particleDensity[p]; // CC: is there a reason in old geos grid material volume was computed from mass and density and did not use particle volume directly?
        
        // TODO: Normalizing by volume might be better
        gridDamage[mappedNode][fieldIndex] += particleMass[p] * ( particleSurfaceFlag[p] == 1 ? 1 : particleDamage[p] ) * shapeFunctionValues[pp][g];
        gridMaxDamage[mappedNode][fieldIndex] = fmax( gridMaxDamage[mappedNode][fieldIndex], particleSurfaceFlag[p] == 1 ? 1 : particleDamage[p] );
        for( int i=0; i<numDims; i++ )
        {
          gridMomentum[mappedNode][fieldIndex][i] += particleMass[p] * particleVelocity[p][i] * shapeFunctionValues[pp][g];
          gridExternalForce[mappedNode][fieldIndex][i] += particleBodyForce[p][i] * particleMass[p] * shapeFunctionValues[pp][g];
          

          // TODO: Switch to volume weighting?
          gridMaterialPosition[mappedNode][fieldIndex][i] += particleMass[p] * (particlePosition[p][i] - gridPosition[mappedNode][i]) * shapeFunctionValues[pp][g];
          for( int k=0; k < numDims; k++ )
          {
            int voigt = voigtMap[k][i];
            // CC: double check this implementation of artificial viscosity is correct
            gridInternalForce[mappedNode][fieldIndex][i] -= ( particleStress[p][voigt] - particleArtificialViscosity[p] * m_useArtificialViscosity * (k == i) ) 
                                                            * shapeFunctionGradientValues[pp][g][k] * particleVolume[p];
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
  GEOS_MARK_FUNCTION;

  // Grid fields
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const gridPosition = nodeManager.referencePosition();
  arrayView2d< real64 const > const & gridMass = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::massString() );
  arrayView3d< real64 > const & gridVelocity = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::velocityString() );
  arrayView3d< real64 > const & gridDVelocity = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::dVelocityString() );
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
          gridDVelocity[g][fieldIndex][i] = gridAcceleration[g][fieldIndex][i] * dt;
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
          gridDVelocity[g][fieldIndex][i] = 0.0;
          gridVelocity[g][fieldIndex][i] = 0.0;
          gridMomentum[g][fieldIndex][i] = 0.0;
          gridMaterialPosition[g][fieldIndex][i] = 0 * gridPosition[g][i]; // TODO: zero? since it's supposed to be relative position?
        }
      }
    } );
  }

  // // CC: debug 
  // // Write out grid Node mass to check in vtk
  // for( int g = 0; g < numNodes; g++)
  // {
  //   std::stringstream line;
  //   for( int fieldIndex = 0; fieldIndex < m_numVelocityFields; fieldIndex++ )
  //   {
  //     // GEOS_LOG_RANK( "g: " << g << ", field" << fieldIndex );
  //     line << gridMass[g][fieldIndex] << " ";
  //   }
  //   GEOS_LOG_RANK( line.str() );
  // }
}

void SolidMechanicsMPM::enforceContact( real64 dt,
                                        DomainPartition & domain,
                                        ParticleManager & particleManager,
                                        NodeManager & nodeManager,
                                        MeshLevel & mesh )
{
  GEOS_MARK_FUNCTION;

  // Grid fields
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const gridPosition = nodeManager.referencePosition();
  arrayView2d< real64 const > const & gridMass = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::massString() );
  arrayView2d< real64 const > const & gridMaterialVolume = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::materialVolumeString() );
  arrayView2d< real64 const > const & gridDamage = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::damageString() );
  arrayView2d< real64 const > const & gridMaxDamage = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::maxDamageString() );
  arrayView3d< real64 > const & gridVelocity = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::velocityString() );
  arrayView3d< real64 > const & gridDVelocity = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::dVelocityString() );
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
                        gridMaterialVolume,
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
          gridDVelocity[g][fieldIndex][i] = contactForce / gridMass[g][fieldIndex] * dt; // Change in velocity from enforcing contact needed for XPIC update, CC: TODO clean up this code to avoid redundant operations
          gridVelocity[g][fieldIndex][i] = gridMomentum[g][fieldIndex][i] / gridMass[g][fieldIndex];
        }
      }
    } );
  }
}

void SolidMechanicsMPM::interpolateTable( real64 x, 
                                          real64 dx,
                                          array2d< real64 > table,
                                          arrayView1d< real64 > output,
                                          SolidMechanicsMPM::InterpolationOption interpolationType )
{
  int numRows = table.size(0);
  int numColumns = table.size(1);
  int numDims = numColumns - 1;

  real64 xNext = x + dx;

  // If we preceed table use first row
  if( xNext <= table[0][0] )
  {
    for( int i = 0; i < numDims; i++ )
    {
      output[i] = table[0][i + 1];
    }
    return;
  }

  // If we exceed table use last row
  if( xNext >= table[numRows - 1][0] )
  {
    for( int i = 0; i < numDims; i++ )
    {
      output[i] = table[numRows - 1][i + 1];
    }
    return;
  }

  int tableInterval = 0;
  for( localIndex i = 0; i < numRows; i++ )
  {
    if( x + dx / 2 > table[i][0] )
    {
      tableInterval = i;
    }
  }

  real64 timeInterval = table[tableInterval + 1][0] - table[tableInterval][0]; // Time fInterval for current part of F table we're in
  real64 timePast = xNext - table[tableInterval][0]; // Time elapsed since switching intervals in F table
  real64 timeFrac = timePast / timeInterval;

  for( int i = 0; i < numDims; i++ )
  {
      switch( interpolationType )
      {
        case SolidMechanicsMPM::InterpolationOption::Linear:
          // default linear interpolation
          output[i] = table[tableInterval][i + 1] * ( 1.0 - timeFrac ) + table[tableInterval + 1][i + 1] * timeFrac;
          break;
        case SolidMechanicsMPM::InterpolationOption::Cosine:
          // smooth-step interpolation with cosine, zero endpoint velocity
          output[i] = table[tableInterval][i + 1] - 0.5 * ( table[tableInterval + 1][i + 1] - table[tableInterval][i + 1] ) * ( cos( 3.141592653589793 * timeFrac ) - 1.0 );
          break;
        case SolidMechanicsMPM::InterpolationOption::Smoothstep:
          // smooth-step interpolation with 5th order polynomial, zero endpoint velocity and acceleration
          output[i] = table[tableInterval][i+1] + ( table[tableInterval+1][i+1] - table[tableInterval][i+1] ) *
                      ( 10.0 * pow( timeFrac, 3 ) - 15.0 * pow( timeFrac, 4 ) + 6.0 * pow( timeFrac, 5 ) );
          break;
        default:
          GEOS_ERROR( "No interpolation option of that type!" );
          break;
      }
  }
}

void SolidMechanicsMPM::interpolateFTable( real64 dt, real64 time_n )
{
  GEOS_MARK_FUNCTION;

  array1d< real64 > Fii_new(3);
  interpolateTable( time_n, 
                    dt,
                    m_fTable,
                    Fii_new,
                    m_fTableInterpType );

  for( int i = 0; i < m_numDims; i++ )
  {
    if( m_stressControl[i] != 1 )
    {
      real64 Fii_dot = ( Fii_new[i] - m_domainF[i] ) / dt;
      m_domainL[i] = Fii_dot / Fii_new[i]; // L = Fdot.Finv
      m_domainF[i] = Fii_new[i];
    }
  }


  // double Fii_new;
  // double Fii_dot;
  // int fInterval = 0;
  // double timeInterval;
  // double timePast;

  // if( time_n + dt < m_fTable[m_fTable.size( 0 ) - 1][0] ) // If within F table bounds...
  // {
  //   for( localIndex i = 0; i < m_fTable.size( 0 ); i++ ) // Naive method for determining what part of F table we're currently in, can
  //                                                        // optimize later (TODO)
  //   {
  //     if( time_n + dt / 2 > m_fTable[i][0] )
  //     {
  //       fInterval = i;
  //     }
  //   }

  //   timeInterval = m_fTable[fInterval + 1][0] - m_fTable[fInterval][0]; // Time fInterval for current part of F table we're in
  //   timePast = time_n + dt - m_fTable[fInterval][0]; // Time elapsed since switching intervals in F table

  //   for( int i = 0; i < m_numDims; i++ )   // Update L and F
  //   {
  //     if ( m_stressControl[i] != 1)
  //     {
  //       // smooth-step interpolation with cosine, zero endpoint velocity
  //       if( m_fTableInterpType == 1 )
  //       {
  //         Fii_new = m_fTable[fInterval][i + 1] - 0.5 * ( m_fTable[fInterval + 1][i + 1] - m_fTable[fInterval][i + 1] ) * ( cos( 3.141592653589793 * timePast / timeInterval ) - 1.0 );
  //       }
  //       // smooth-step interpolation with 5th order polynomial, zero endpoint velocity and acceleration
  //       else if( m_fTableInterpType == 2 )
  //       {
  //         Fii_new = m_fTable[fInterval][i+1] + (m_fTable[fInterval+1][i+1] - m_fTable[fInterval][i+1])*
  //                   (10.0*pow( timePast/timeInterval, 3 ) - 15.0*pow( timePast/timeInterval, 4 ) + 6.0*pow( timePast/timeInterval, 5 ));
  //       }
  //       // default linear interpolation
  //       else
  //       {
  //         Fii_new = m_fTable[fInterval][i + 1] * ( timeInterval - timePast ) / timeInterval + m_fTable[fInterval + 1][i + 1] * ( timePast ) / timeInterval;
  //       }

  //       // // drive F22 and F33 using F11 to get pure shear
  //       // if(m_fTable_pure_shear && i > 0)
  //       // {
  //       //   Fii_new = m_planeStrain ? 1.0/m_domainF[0] : 1.0/sqrt(m_domainF[0]);
  //       // }

  //       Fii_dot = ( Fii_new - m_domainF[i] ) / dt;
  //       m_domainL[i] = Fii_dot / Fii_new; // L = Fdot.Finv
  //       m_domainF[i] = Fii_new;
  //     }
  //   }
  // }
  // else if( time_n + dt >= m_fTable[m_fTable.size( 0 ) - 1][0] ) // Else (i.e. if we exceed F table upper bound)...
  // {
  //   for( int i = 0; i < m_numDims; i++ )
  //   {
  //     if( m_stressControl[i] != 1 )
  //     {
  //     Fii_new = m_fTable[m_fTable.size( 0 ) - 1][i + 1]; // Set Fii_new to the last prescribed F table value
  //     Fii_dot = ( Fii_new - m_domainF[i] ) / dt;
  //     m_domainL[i] = Fii_dot / Fii_new; // L = Fdot.Finv
  //     m_domainF[i] = Fii_new;
  //     }
  //   }
  // }
}

void SolidMechanicsMPM::interpolateStressTable( real64 dt,
                                                real64 time_n )
{
  GEOS_MARK_FUNCTION;

  interpolateTable( time_n, 
                    dt,
                    m_stressTable,
                    m_domainStress,
                    m_fTableInterpType );

  // double stressii_new;
  // int stressInterval = 0;
  // double timeInterval;
  // double timePast;

  // if( time_n + dt < m_stressTable[m_stressTable.size( 0 ) - 1][0] ) // If within F table bounds...
  // {
  //   for( localIndex i = 0; i < m_stressTable.size( 0 ); i++ ) // Naive method for determining what part of F table we're currently in, can
  //                                                        // optimize later (TODO)
  //   {
  //     if( time_n + dt / 2 > m_stressTable[i][0] )
  //     {
  //       stressInterval = i;
  //     }
  //   }

  //   timeInterval = m_stressTable[stressInterval + 1][0] - m_stressTable[stressInterval][0]; // Time fInterval for current part of F table we're in
  //   timePast = time_n + dt - m_stressTable[stressInterval][0]; // Time elapsed since switching intervals in F table

  //   for( int i = 0; i < m_numDims; i++ )   // Update L and F
  //   {
  //     // smooth-step interpolation with cosine, zero endpoint velocity
  //     if( m_stressTableInterpType == 1 )
  //     {
  //       stressii_new = m_stressTable[stressInterval][i + 1] - 0.5 * ( m_stressTable[stressInterval + 1][i + 1] - m_stressTable[stressInterval][i + 1] ) * ( cos( 3.141592653589793 * timePast / timeInterval ) - 1.0 );
  //     }
  //     // smooth-step interpolation with 5th order polynomial, zero endpoint velocity and acceleration
  //     else if( m_stressTableInterpType == 2 )
  //     {
  //       stressii_new = m_stressTable[stressInterval][i+1] + (m_stressTable[stressInterval+1][i+1] - m_stressTable[stressInterval][i+1])*
  //                 (10.0*pow( timePast/timeInterval, 3 ) - 15.0*pow( timePast/timeInterval, 4 ) + 6.0*pow( timePast/timeInterval, 5 ));
  //     }
  //     // default linear interpolation
  //     else
  //     {
  //       stressii_new = m_stressTable[stressInterval][i + 1] * ( timeInterval - timePast ) / timeInterval + m_stressTable[stressInterval + 1][i + 1] * ( timePast ) / timeInterval;
  //     }

  //     m_domainStress[i] = stressii_new;
  //   }
  // }
  // else if( time_n + dt >= m_stressTable[m_stressTable.size( 0 ) - 1][0] ) // Else (i.e. if we exceed F table upper bound)...
  // {
  //   for( int i = 0; i < m_numDims; i++ )
  //   {
  //     stressii_new = m_stressTable[m_stressTable.size( 0 ) - 1][i + 1]; // Set Fii_new to the last prescribed F table value
  //     m_domainStress[i] = stressii_new;
  //   }
  // }
}

void SolidMechanicsMPM::gridToParticle( real64 dt,
                                        ParticleManager & particleManager,
                                        NodeManager & nodeManager,
                                        DomainPartition & domain,
                                        MeshLevel & mesh )
{
  GEOS_MARK_FUNCTION;
  
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
                         nodeManager,
                         domain,
                         mesh );
      break;
    case SolidMechanicsMPM::UpdateMethodOption::FMPM:
      performFMPMUpdate( dt, 
                         particleManager, 
                         nodeManager,
                         domain,
                         mesh );
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

      for( int i=0; i<3; i++ )
      {
        particleVelocity[p][i] = 0.0;
        for( int j=0; j<3; j++ )
        {
          particleVelocityGradient[p][i][j] = 0.0;
        }
      }
    
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
            particleVelocityGradient[p][i][j] += gridVelocity[mappedNode][fieldIndex][i] * shapeFunctionGradientValues[pp][g][j];
          }
        }
      }
    } );
    subRegionIndex++;
  } );
}

void SolidMechanicsMPM::performXPICUpdate( real64 dt,
                                           ParticleManager & particleManager,
                                           NodeManager & nodeManager,
                                           DomainPartition & domain, 
                                           MeshLevel & mesh )
{
  // Grid fields
  arrayView2d< real64 const > const & gridMass = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::massString() );
  arrayView2d< real64 const > const & gridDamageGradient = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::damageGradientString() );
  // arrayView3d< real64 const > const & gridDVelocity = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::dVelocityString() ); // for multifield contact corrections
  arrayView3d< real64 const > const & gridVelocity = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::velocityString() );
  arrayView3d< real64 const > const & gridAcceleration = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::accelerationString() );
  
  arrayView3d< real64 > const & gridVPlus = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::vPlusString() );
  arrayView3d< real64 > const & gridDVPlus = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::dVPlusString() );

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
    arrayView3d< real64 const > const shapeFunctionGradientValues = m_shapeFunctionGradientValues[subRegionIndex];

    // Map to particles
    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
    int const numDims = m_numDims;
    int const damageFieldPartitioning = m_damageFieldPartitioning;
    int const numContactGroups = m_numContactGroups;
    int const numVelocityFields = m_numVelocityFields;
    int const updateOrder = m_updateOrder;

    // For iterative XPIC solve
    array3d< real64 > vStar( numNodes, numVelocityFields, numDims );
    array3d< real64 > vMinus( numNodes, numVelocityFields, numDims );

    array3d< real64 > dVMinus( numNodes, numVelocityFields, numDims );

    // Zero out vStar for each order iteration
    for( int n=0; n < numNodes; n++)
    {
      for( int cg=0; cg < numVelocityFields; cg++)
      {
        for( int i = 0; i < numDims; i++)
        {
          vStar[n][cg][i] = 0.0;
          vMinus[n][cg][i] = gridVelocity[n][cg][i] - gridAcceleration[n][cg][i] * dt;

          gridDVPlus[n][cg][i] = 0.0; // gridDVelocity[n][cg][i]; // CC: this isn't working, currently disabled by writing 0 to it
          dVMinus[n][cg][i] = gridDVPlus[n][cg][i];
        }
      }
    }

    // Do XPIC iterations
    for(int r=2; r <= updateOrder; ++r)
    {
      // Zero out vPlus for each order iteration
      for( int n=0; n < numNodes; n++)
      {
        for( int cg=0; cg < numVelocityFields; cg++)
        {
          for( int i = 0; i < numDims; i++)
          {
            gridDVPlus[n][cg][i] = 0.0;
          }
        }
      }

      // Seems like LvArrays that aren't views need to be passed by reference to forAll loops (e.g. vPlus)
      forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
      {
        localIndex const p = activeParticleIndices[pp];
      
        for( int gi = 0; gi < 8 * numberOfVerticesPerParticle; gi++ )
        {
          localIndex const mappedNodeI = mappedNodes[pp][gi];
          int const nodeFlagI = ( damageFieldPartitioning == 1 && LvArray::tensorOps::AiBi< 3 >( gridDamageGradient[mappedNodeI], particleDamageGradient[p] ) < 0.0 ) ? 1 : 0;
          int const fieldIndexI = nodeFlagI * numContactGroups + particleGroup[p]; // This ranges from 0 to nMatFields-1
          
          if( gridMass[mappedNodeI][fieldIndexI] > m_smallMass )
          {
            for(int gj = 0; gj < 8 * numberOfVerticesPerParticle; gj++ )
            {
              localIndex const mappedNodeJ = mappedNodes[pp][gj];
              int const nodeFlagJ = ( damageFieldPartitioning == 1 && LvArray::tensorOps::AiBi< 3 >( gridDamageGradient[mappedNodeJ], particleDamageGradient[p] ) < 0.0 ) ? 1 : 0;
              int const fieldIndexJ = nodeFlagJ * numContactGroups + particleGroup[p]; // This ranges from 0 to nMatFields-1

              for (int i = 0; i < numDims; i++)
              {
                gridDVPlus[mappedNodeI][fieldIndexI][i] += ( ( updateOrder - r + 1.0 ) / r ) * 
                                                        ( particleMass[p] * shapeFunctionValues[pp][gi] * shapeFunctionValues[pp][gj] / gridMass[mappedNodeI][fieldIndexI] ) * 
                                                          vMinus[mappedNodeJ][fieldIndexJ][i];
              }
            }

          }
        }
      } );

      // syncGridFields( { viewKeyStruct::vPlusString() }, domain, nodeManager, mesh, MPI_SUM );
      syncGridFields( { viewKeyStruct::vPlusString(), viewKeyStruct::dVPlusString() }, domain, nodeManager, mesh, MPI_SUM ); // Also need to sync dVPlus later for multifield contact correction

      // Update vStar
      for( int n=0; n < numNodes; n++)
      {
        for( int cg=0; cg < numVelocityFields; cg++)
        {
          for( int i = 0; i < numDims; i++)
          {
            vStar[n][cg][i] += std::pow(-1.0, r) * ( gridVPlus[n][cg][i] + ( updateOrder - 1.0 ) / updateOrder * gridDVPlus[n][cg][i] );
            vMinus[n][cg][i] = gridVPlus[n][cg][i];
          }
        }
      }

      // Add if statement? because this doesn't need to be performed for the last updateOrder
      // Zero out vPlus for each order iteration
      for( int n=0; n < numNodes; n++)
      {
        for( int cg=0; cg < numVelocityFields; cg++)
        {
          for( int i = 0; i < numDims; i++)
          {
            gridDVPlus[n][cg][i] = 0.0;
          }
        }
      }

      // Update dVStar
      forAll< serialPolicy >( activeParticleIndices.size(), [&] GEOS_HOST_DEVICE ( localIndex const pp )
      {
        localIndex const p = activeParticleIndices[pp];
      
        for( int gi = 0; gi < 8 * numberOfVerticesPerParticle; gi++ )
        {
          localIndex const mappedNodeI = mappedNodes[pp][gi];
          int const nodeFlagI = ( damageFieldPartitioning == 1 && LvArray::tensorOps::AiBi< 3 >( gridDamageGradient[mappedNodeI], particleDamageGradient[p] ) < 0.0 ) ? 1 : 0;
          int const fieldIndexI = nodeFlagI * numContactGroups + particleGroup[p]; // This ranges from 0 to nMatFields-1
          
          if( gridMass[mappedNodeI][fieldIndexI] > m_smallMass )
          {
            for(int gj = 0; gj < 8 * numberOfVerticesPerParticle; gj++ )
            {
              localIndex const mappedNodeJ = mappedNodes[pp][gj];
              int const nodeFlagJ = ( damageFieldPartitioning == 1 && LvArray::tensorOps::AiBi< 3 >( gridDamageGradient[mappedNodeJ], particleDamageGradient[p] ) < 0.0 ) ? 1 : 0;
              int const fieldIndexJ = nodeFlagJ * numContactGroups + particleGroup[p]; // This ranges from 0 to nMatFields-1

              for (int i = 0; i < numDims; i++)
              {
                gridDVPlus[mappedNodeI][fieldIndexI][i] += ( ( updateOrder - r ) / r ) * 
                                                        ( particleMass[p] * shapeFunctionValues[pp][gi] * shapeFunctionValues[pp][gj] / gridMass[mappedNodeI][fieldIndexI] ) * 
                                                          dVMinus[mappedNodeJ][fieldIndexJ][i];
              }
            }
          }
        }
      } );

      for( int n=0; n < numNodes; n++)
      {
        for( int cg=0; cg < numVelocityFields; cg++)
        {
          for( int i = 0; i < numDims; i++)
          {
            dVMinus[n][cg][i] = gridDVPlus[n][cg][i];
          }
        }
      }

    } //End of updateOrder iterations

    // Update particles position and velocities now
    forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
    {
      localIndex const p = activeParticleIndices[pp];
    
      // Zero velocity gradient
      for( int i=0; i < numDims; i++ )
      {
        particlePosition[p][i] -= particleVelocity[p][i] * dt / 2.0;
        particleVelocity[p][i] = 0.0;
        for( int j=0; j < numDims; j++ )
        {
          particleVelocityGradient[p][i][j] = 0.0;
        }
      }

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

          particlePosition[p][i] += S * gVPlus * dt  - ( S * gA * dt - m * S *( gVPlus - gA * dt ) + m * S * vStar[mappedNode][fieldIndex][i] ) * dt / 2.0;
          // particlePosition[p][i] += 0.5 * S * dt * ( (m + 2.0 ) * gVPlus - ( 1.0 + m ) * gA * dt - m * vStar[mappedNode][fieldIndex][i] );
          // particlePosition[p][i] += S * gVPlus * dt - ( S * gA * dt + ( -m * S * ( gVPlus - gA * dt ) + m * S * vStar[mappedNode][fieldIndex][i]) ) * dt / 2.0;
          particleVelocity[p][i] += S * ( m * ( gVPlus - vStar[mappedNode][fieldIndex][i] ) + ( 1 - m ) * gA * dt );
          // CC: What about update to velocity gradient?
          // Currently copy this from FLIP udpate with change from gridVelocity to vStar
          for( int j=0; j < numDims; j++ )
          {
            particleVelocityGradient[p][i][j] += gridVelocity[mappedNode][fieldIndex][i] * shapeFunctionGradientValues[pp][g][j];
          }
        }
      }
    } );
    subRegionIndex++;
  } );
}

void SolidMechanicsMPM::performFMPMUpdate(  real64 dt,
                                            ParticleManager & particleManager,
                                            NodeManager & nodeManager,
                                            DomainPartition & domain, 
                                            MeshLevel & mesh )
{
  // Grid fields
  arrayView2d< real64 const > const & gridMass = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::massString() );
  arrayView2d< real64 const > const & gridDamageGradient = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::damageGradientString() );
  arrayView3d< real64 const > const & gridVelocity = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::velocityString() );

  arrayView3d< real64 > const & gridVPlus = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::vPlusString() );

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
    arrayView3d< real64 const > const shapeFunctionGradientValues = m_shapeFunctionGradientValues[subRegionIndex];

    // Map to particles
    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
    int const numDims = m_numDims;
    int const damageFieldPartitioning = m_damageFieldPartitioning;
    int const numContactGroups = m_numContactGroups;
    int const numVelocityFields = m_numVelocityFields;
    int const updateOrder = m_updateOrder;

    // Added these as fields to nodeManager to easily sync them in parallelization for each iteration (technically only vPlus needs to be synced )
    // For iterative FMPM solve
    array3d< real64 > vStar( numNodes, numVelocityFields, numDims );
    array3d< real64 > vMinus( numNodes, numVelocityFields, numDims );
    
    // Initialize FMPM variables
    for( int n=0; n < numNodes; n++)
    {
      for( int cg=0; cg < numVelocityFields; cg++)
      {
        for( int i = 0; i < numDims; i++)
        {
          vStar[n][cg][i] = updateOrder * gridVelocity[n][cg][i];
          vMinus[n][cg][i] = vStar[n][cg][i];
        }
      }
    }

    // Perform FMPM order iterations
    for(int r=2; r <= updateOrder; ++r)
    {
      // Zero out vPlus for each order iteration
      for( int n=0; n < numNodes; n++)
      {
        for( int cg=0; cg < numVelocityFields; cg++)
        {
          for( int i = 0; i < numDims; i++)
          {
            gridVPlus[n][cg][i] = 0.0;
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

          if( gridMass[mappedNodeI][fieldIndexI] > m_smallMass )
          {
            real64 Splus = particleMass[p] * shapeFunctionValues[pp][gi] / gridMass[mappedNodeI][fieldIndexI];

            for(int gj = 0; gj < 8 * numberOfVerticesPerParticle; gj++ )
            {
              localIndex const mappedNodeJ = mappedNodes[pp][gj];
              int const nodeFlagJ = ( damageFieldPartitioning == 1 && LvArray::tensorOps::AiBi< 3 >( gridDamageGradient[mappedNodeJ], particleDamageGradient[p] ) < 0.0 ) ? 1 : 0;
              int const fieldIndexJ = nodeFlagJ * numContactGroups + particleGroup[p]; // This ranges from 0 to nMatFields-1

              for (int i = 0; i < numDims; i++)
              {
                gridVPlus[mappedNodeI][fieldIndexI][i] += Splus * shapeFunctionValues[pp][gj] * vMinus[mappedNodeJ][fieldIndexJ][i];;
              }
            }
          }
        }
      } );

      syncGridFields( { viewKeyStruct::vPlusString() }, domain, nodeManager, mesh, MPI_SUM );

      // Update vStar
      real64 orderCoefficient = std::pow(-1.0, 1.0+r) * ( updateOrder - r + 1.0 ) / r;
      for( int n=0; n < numNodes; n++)
      {
        for( int cg=0; cg < numVelocityFields; cg++)
        {
          for( int i = 0; i < numDims; i++)
          {
            vStar[n][cg][i] += orderCoefficient * gridVPlus[n][cg][i];
            vMinus[n][cg][i] = gridVPlus[n][cg][i];
          }
        }
      }
    } //End of updateOrder iterations

    // Update particles positions and velocities now
    forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
    {
      localIndex const p = activeParticleIndices[pp];

      for(int i = 0; i < numDims; i++)
      {
        particlePosition[p][i] += 0.5 * dt * particleVelocity[p][i];
        particleVelocity[p][i] = 0.0;

        for( int j=0; j < numDims; j++ )
        {
          particleVelocityGradient[p][i][j] = 0.0;
        }
      }

      for( int g = 0; g < 8 * numberOfVerticesPerParticle; g++ )
      {
        localIndex const mappedNode = mappedNodes[pp][g];
        int const nodeFlag = ( damageFieldPartitioning == 1 && LvArray::tensorOps::AiBi< 3 >( gridDamageGradient[mappedNode], particleDamageGradient[p] ) < 0.0 ) ? 1 : 0;
        int const fieldIndex = nodeFlag * numContactGroups + particleGroup[p];

        for( int i=0; i < numDims; i++ )
        {
          particlePosition[p][i] += 0.5 * dt * shapeFunctionValues[pp][g] * vStar[mappedNode][fieldIndex][i];
          particleVelocity[p][i] += shapeFunctionValues[pp][g] * vStar[mappedNode][fieldIndex][i];

          for( int j=0; j < numDims; j++ )
          {
            particleVelocityGradient[p][i][j] += vStar[mappedNode][fieldIndex][i] * shapeFunctionGradientValues[pp][g][j];
          }
        }
      }
    } );
    subRegionIndex++;
  } );
}

void SolidMechanicsMPM::updateSolverDependencies( ParticleManager & particleManager )
{
  GEOS_MARK_FUNCTION;

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

    if(  solidModel.hasWrapper( "temperature" ) )
    {
      arrayView1d< real64 > const particleTemperature = subRegion.getField< fields::mpm::particleTemperature >();
      arrayView1d< real64 > const constitutiveTemperature = solidModel.getReference< array1d< real64 > >( "temperature" );
      forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
      {
        localIndex const p = activeParticleIndices[pp];
        particleTemperature[p] = constitutiveTemperature[p]; 
      } );
    }
  } );
}

real64 SolidMechanicsMPM::getStableTimeStep( ParticleManager & particleManager )
{
  GEOS_MARK_FUNCTION;

  real64 maxWavespeed = 0.0;
  real64 length = m_planeStrain == 1 ? std::fmin( m_hEl[0], m_hEl[1] ) : std::fmin( m_hEl[0], std::fmin( m_hEl[1], m_hEl[2] ) );

  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    arrayView1d< real64 const > const particleWavespeed = subRegion.getField< fields::mpm::particleWavespeed >();
    arrayView2d< real64 const > const particleVelocity = subRegion.getParticleVelocity();

    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
    forAll< serialPolicy >( activeParticleIndices.size(), [=, &maxWavespeed] GEOS_HOST ( localIndex const pp ) // would need reduction to parallelize
    {
      localIndex const p = activeParticleIndices[pp];
      maxWavespeed = fmax( maxWavespeed, particleWavespeed[p] + LvArray::tensorOps::l2Norm< 3 >( particleVelocity[p] ) );
    } );
  } );

  real64 dtReturn = maxWavespeed > 1.0e-16 ? m_cflFactor * length / maxWavespeed : DBL_MAX; // This partitions's dt, make it huge if wavespeed=0.0
                                                                                      // (this happens when there are no particles on this
                                                                                      // partition)
  
  // CC: TODO add check to make sure that next time step won't skip a MPM event
  return dtReturn;
}

void SolidMechanicsMPM::deleteBadParticles( ParticleManager & particleManager )
{
  GEOS_MARK_FUNCTION;

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
  GEOS_MARK_FUNCTION;

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
                                                      real64 const & maxDamageB,
                                                      arraySlice1d< real64 const > const xA,
                                                      arraySlice1d< real64 const > const xB )
// m_treatFullyDamagedAsSingleField makes fields inseparable if damageA = damageB = 1, so we aren't putting
// arbitrary separation planes (and potential surfaces for accumulated overlap that needs corrections)
// between fully damaged materials. There is a potential issue that approaching damaged bodies
// would get contact gaps locked in.
{
  int separable = 0;
  // At least one field is fully damaged and both fields have the minimum separable level of damage.
  // The "a%b" is the "mod(a,b)" command, and indicates whether materials are from same contact group.
  if( ( ( maxDamageA >= 0.9999 || maxDamageB >= 0.9999 ) && 
        ( damageA >= m_separabilityMinDamage && damageB >= m_separabilityMinDamage ) ) || 
      ( A % m_numContactGroups != B % m_numContactGroups ) )
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

  // For a thin (relative to grid spacing) strip of material with surface flags or damaged surfaces, the material at opposing surfaces will have opposing
	// damage-field gradients, and this will create a spurious slip surfaces within the material.
	// We can detect this case by seeing if the distance between the position of the two fields is small relative to the grid cell spacing.  This will
	// only happen in the case of a thin strip of material when using DFG, so in this case we set separable = 0
	// The threshold spacing is normalized by the neighbor radius.
	if ( m_thinFeatureDFGThreshold < 10.0 ) // Threshold defaults to 1.e99 so we only do this calculation if a lower value is set.
	{                                        //After testing we may default to a lower value.
		real64 dx[3] = { 0 };
    LvArray::tensorOps::copy< 3 >( dx, xA );
    LvArray::tensorOps::subtract< 3 >( dx, xB );

    real64 productOfSquares = ( dx[0] * dx[0] ) * ( dx[1] * dx[1] ) * ( dx[2] * dx[2] );

    if ( productOfSquares < ( m_thinFeatureDFGThreshold * m_neighborRadius ) * ( m_thinFeatureDFGThreshold * m_neighborRadius ) )
		{
			separable = 0;
		}
	}

  return separable;
}

// CC: do I need to modify this to check for periodic boundaries
// All master particles should have centers inside the domain if particleCenters are corrected correctly during repartitioning
// CPDI: Edge case, corner of large or long particle beyond ghost cells but center is still inside domain?
void SolidMechanicsMPM::flagOutOfRangeParticles( ParticleManager & particleManager )
{
  GEOS_MARK_FUNCTION;

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
  GEOS_MARK_FUNCTION;

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
  GEOS_MARK_FUNCTION;

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
  GEOS_MARK_FUNCTION;

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
  GEOS_MARK_FUNCTION;

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
                                               NodeManager & nodeManager ) //,
                                              //  SpatialPartition & partition )
{
  GEOS_MARK_FUNCTION;

  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const gridPosition = nodeManager.referencePosition();
  arrayView3d< int const > const ijkMap = m_ijkMap;
  real64 hEl[3] = {0};
  LvArray::tensorOps::copy< 3 >( hEl, m_hEl );
  real64 xLocalMin[3] = {0};
  LvArray::tensorOps::copy< 3 >( xLocalMin, m_xLocalMin );
  real64 xLocalMax[3] = {0};
  LvArray::tensorOps::copy< 3 >( xLocalMax, m_xLocalMax);

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

    // GEOS_LOG_RANK( "mappedNodes=" << mappedNodes.size( 0 ) << ", " << mappedNodes.size( 1 ) << ", " << 
    //                "shapeFunctions=" << shapeFunctionValues.size( 0 ) << ", " << shapeFunctionValues.size( 1 ) << ", " << 
    //                "shapeFunctionGradients=" << shapeFunctionGradientValues.size( 0 ) << ", " << shapeFunctionGradientValues.size( 1 ) << ", " << shapeFunctionGradientValues.size( 2 ) << ", " << 
    //                "activeParticleIndices=" << activeParticleIndices.size() << ", " << 
    //                "particlePosition=" << particlePosition.size( 0 ) << ", " << particlePosition.size( 1 ) << ", " << 
    //                "gridPosition=" << gridPosition.size( 0 ) << ", " << gridPosition.size( 1 ) << ", " <<
    //                "ijkMap=" << ijkMap.size( 0 ) << ", " << ijkMap.size( 1 ) << ", " << ijkMap.size( 2 ) );

    // Populate mapping arrays based on particle type
    switch( subRegion.getParticleType() )
    {
      case ParticleType::SinglePoint:
        {
          forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
        {
          localIndex const p = activeParticleIndices[pp];

          computeSinglePointShapeFunctions( gridPosition,
                                            particlePosition[p],
                                            ijkMap,
                                            xLocalMin,
                                            hEl,
                                            mappedNodes[pp],
                                            shapeFunctionValues[pp],
                                            shapeFunctionGradientValues[pp] );
          // // get IJK associated with particle center
          // int centerIJK[3];
          // for( int i=0; i<3; i++ )
          // {
          //   centerIJK[i] = floor( ( particlePosition[p][i] - xLocalMin[i] ) / hEl[i] );
          // }

          // // get node IDs, weights and grad weights
          // int node = 0;
          // int corner = ijkMap[centerIJK[0]][centerIJK[1]][centerIJK[2]];
          // auto corner_x = gridPosition[corner];

          // real64 xRel = (particlePosition[p][0] - corner_x[0]) / hEl[0];
          // real64 yRel = (particlePosition[p][1] - corner_x[1]) / hEl[1];
          // real64 zRel = (particlePosition[p][2] - corner_x[2]) / hEl[2];

          // for( int i=0; i<2; i++ )
          // {
          //   real64 xWeight = i*xRel + (1-i)*(1.0-xRel);
          //   real64 dxWeight = i/hEl[0] - (1-i)/hEl[0];
          //   for( int j=0; j<2; j++ )
          //   {
          //     real64 yWeight = j*yRel + (1-j)*(1.0-yRel);
          //     real64 dyWeight = j/hEl[1] - (1-j)/hEl[1];
          //     for( int k=0; k<2; k++ )
          //     {
          //       real64 zWeight = k*zRel + (1-k)*(1.0-zRel);
          //       real64 dzWeight = k/hEl[2] - (1-k)/hEl[2];
          //       mappedNodes[pp][node] = ijkMap[centerIJK[0]+i][centerIJK[1]+j][centerIJK[2]+k];
          //       shapeFunctionValues[pp][node] = xWeight * yWeight * zWeight;
          //       shapeFunctionGradientValues[pp][node][0] = dxWeight * yWeight * zWeight;
          //       shapeFunctionGradientValues[pp][node][1] = xWeight * dyWeight * zWeight;
          //       shapeFunctionGradientValues[pp][node][2] = xWeight * yWeight * dzWeight;
          //       node++;
          //     }
          //   }
          // }
        } );
          break;
        }
      case ParticleType::CPDI:
        {
          // int const signs[8][3] = { { -1, -1, -1 },
          //   {  1, -1, -1 },
          //   {  1, 1, -1 },
          //   { -1, 1, -1 },
          //   { -1, -1, 1 },
          //   {  1, -1, 1 },
          //   {  1, 1, 1 },
          //   { -1, 1, 1 } };
          arrayView3d< real64 const > const particleRVectors = subRegion.getParticleRVectors();
      
          forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
        {
          localIndex const p = activeParticleIndices[pp];

          computeCPDIShapeFunctions( gridPosition,
                                     particlePosition[p],
                                     particleRVectors[p],
                                     ijkMap,
                                     xLocalMin,
                                     hEl,
                                     mappedNodes[pp],
                                     shapeFunctionValues[pp],
                                     shapeFunctionGradientValues[pp] );

          // real64 alpha[8][8];
          // real64 cpdiVolume, oneOverV;
          // real64 p_r1[3], p_r2[3], p_r3[3]; // allowing 1-indexed r-vectors to persist to torture future postdocs >:)

          // for( int i=0; i<3; i++ )
          // {
          //   p_r1[i] = particleRVectors[p][0][i];
          //   p_r2[i] = particleRVectors[p][1][i];
          //   p_r3[i] = particleRVectors[p][2][i];
          // }

          // // We need this because of CPDI domain scaling
          // cpdiVolume = 8.0*(-(p_r1[2]*p_r2[1]*p_r3[0]) + p_r1[1]*p_r2[2]*p_r3[0] + p_r1[2]*p_r2[0]*p_r3[1] - p_r1[0]*p_r2[2]*p_r3[1] - p_r1[1]*p_r2[0]*p_r3[2] + p_r1[0]*p_r2[1]*p_r3[2]);
          // oneOverV = 1.0 / cpdiVolume;

          // alpha[0][0]=oneOverV*(p_r1[2]*p_r2[1] - p_r1[1]*p_r2[2] - p_r1[2]*p_r3[1] + p_r2[2]*p_r3[1] + p_r1[1]*p_r3[2] - p_r2[1]*p_r3[2]);
          // alpha[0][1]=oneOverV*(-(p_r1[2]*p_r2[0]) + p_r1[0]*p_r2[2] + p_r1[2]*p_r3[0] - p_r2[2]*p_r3[0] - p_r1[0]*p_r3[2] + p_r2[0]*p_r3[2]);
          // alpha[0][2]=oneOverV*(p_r1[1]*p_r2[0] - p_r1[0]*p_r2[1] - p_r1[1]*p_r3[0] + p_r2[1]*p_r3[0] + p_r1[0]*p_r3[1] - p_r2[0]*p_r3[1]);
          // alpha[1][0]=oneOverV*(p_r1[2]*p_r2[1] - p_r1[1]*p_r2[2] - p_r1[2]*p_r3[1] - p_r2[2]*p_r3[1] + p_r1[1]*p_r3[2] + p_r2[1]*p_r3[2]);
          // alpha[1][1]=oneOverV*(-(p_r1[2]*p_r2[0]) + p_r1[0]*p_r2[2] + p_r1[2]*p_r3[0] + p_r2[2]*p_r3[0] - p_r1[0]*p_r3[2] - p_r2[0]*p_r3[2]);
          // alpha[1][2]=oneOverV*(p_r1[1]*p_r2[0] - p_r1[0]*p_r2[1] - p_r1[1]*p_r3[0] - p_r2[1]*p_r3[0] + p_r1[0]*p_r3[1] + p_r2[0]*p_r3[1]);
          // alpha[2][0]=oneOverV*(p_r1[2]*p_r2[1] - p_r1[1]*p_r2[2] + p_r1[2]*p_r3[1] - p_r2[2]*p_r3[1] - p_r1[1]*p_r3[2] + p_r2[1]*p_r3[2]);
          // alpha[2][1]=oneOverV*(-(p_r1[2]*p_r2[0]) + p_r1[0]*p_r2[2] - p_r1[2]*p_r3[0] + p_r2[2]*p_r3[0] + p_r1[0]*p_r3[2] - p_r2[0]*p_r3[2]);
          // alpha[2][2]=oneOverV*(p_r1[1]*p_r2[0] - p_r1[0]*p_r2[1] + p_r1[1]*p_r3[0] - p_r2[1]*p_r3[0] - p_r1[0]*p_r3[1] + p_r2[0]*p_r3[1]);
          // alpha[3][0]=oneOverV*(p_r1[2]*p_r2[1] - p_r1[1]*p_r2[2] + p_r1[2]*p_r3[1] + p_r2[2]*p_r3[1] - p_r1[1]*p_r3[2] - p_r2[1]*p_r3[2]);
          // alpha[3][1]=oneOverV*(-(p_r1[2]*p_r2[0]) + p_r1[0]*p_r2[2] - p_r1[2]*p_r3[0] - p_r2[2]*p_r3[0] + p_r1[0]*p_r3[2] + p_r2[0]*p_r3[2]);
          // alpha[3][2]=oneOverV*(p_r1[1]*p_r2[0] - p_r1[0]*p_r2[1] + p_r1[1]*p_r3[0] + p_r2[1]*p_r3[0] - p_r1[0]*p_r3[1] - p_r2[0]*p_r3[1]);
          // alpha[4][0]=oneOverV*(-(p_r1[2]*p_r2[1]) + p_r1[1]*p_r2[2] - p_r1[2]*p_r3[1] + p_r2[2]*p_r3[1] + p_r1[1]*p_r3[2] - p_r2[1]*p_r3[2]);
          // alpha[4][1]=oneOverV*(p_r1[2]*p_r2[0] - p_r1[0]*p_r2[2] + p_r1[2]*p_r3[0] - p_r2[2]*p_r3[0] - p_r1[0]*p_r3[2] + p_r2[0]*p_r3[2]);
          // alpha[4][2]=oneOverV*(-(p_r1[1]*p_r2[0]) + p_r1[0]*p_r2[1] - p_r1[1]*p_r3[0] + p_r2[1]*p_r3[0] + p_r1[0]*p_r3[1] - p_r2[0]*p_r3[1]);
          // alpha[5][0]=oneOverV*(-(p_r1[2]*p_r2[1]) + p_r1[1]*p_r2[2] - p_r1[2]*p_r3[1] - p_r2[2]*p_r3[1] + p_r1[1]*p_r3[2] + p_r2[1]*p_r3[2]);
          // alpha[5][1]=oneOverV*(p_r1[2]*p_r2[0] - p_r1[0]*p_r2[2] + p_r1[2]*p_r3[0] + p_r2[2]*p_r3[0] - p_r1[0]*p_r3[2] - p_r2[0]*p_r3[2]);
          // alpha[5][2]=oneOverV*(-(p_r1[1]*p_r2[0]) + p_r1[0]*p_r2[1] - p_r1[1]*p_r3[0] - p_r2[1]*p_r3[0] + p_r1[0]*p_r3[1] + p_r2[0]*p_r3[1]);
          // alpha[6][0]=oneOverV*(-(p_r1[2]*p_r2[1]) + p_r1[1]*p_r2[2] + p_r1[2]*p_r3[1] - p_r2[2]*p_r3[1] - p_r1[1]*p_r3[2] + p_r2[1]*p_r3[2]);
          // alpha[6][1]=oneOverV*(p_r1[2]*p_r2[0] - p_r1[0]*p_r2[2] - p_r1[2]*p_r3[0] + p_r2[2]*p_r3[0] + p_r1[0]*p_r3[2] - p_r2[0]*p_r3[2]);
          // alpha[6][2]=oneOverV*(-(p_r1[1]*p_r2[0]) + p_r1[0]*p_r2[1] + p_r1[1]*p_r3[0] - p_r2[1]*p_r3[0] - p_r1[0]*p_r3[1] + p_r2[0]*p_r3[1]);
          // alpha[7][0]=oneOverV*(-(p_r1[2]*p_r2[1]) + p_r1[1]*p_r2[2] + p_r1[2]*p_r3[1] + p_r2[2]*p_r3[1] - p_r1[1]*p_r3[2] - p_r2[1]*p_r3[2]);
          // alpha[7][1]=oneOverV*(p_r1[2]*p_r2[0] - p_r1[0]*p_r2[2] - p_r1[2]*p_r3[0] - p_r2[2]*p_r3[0] + p_r1[0]*p_r3[2] + p_r2[0]*p_r3[2]);
          // alpha[7][2]=oneOverV*(-(p_r1[1]*p_r2[0]) + p_r1[0]*p_r2[1] + p_r1[1]*p_r3[0] + p_r2[1]*p_r3[0] - p_r1[0]*p_r3[1] - p_r2[0]*p_r3[1]);

          // // get IJK associated with each corner
          // int cornerIJK[8][3]; // CPDI can map to up to 8 cells
          // for( int corner=0; corner<8; corner++ )
          // {
          //   for( int i=0; i<3; i++ )
          //   {
          //     real64 cornerPositionComponent = particlePosition[p][i] + 
          //                                      signs[corner][0] * particleRVectors[p][0][i] + 
          //                                      signs[corner][1] * particleRVectors[p][1][i] + 
          //                                      signs[corner][2] * particleRVectors[p][2][i];                             

          //     cornerIJK[corner][i] = std::floor( ( cornerPositionComponent - xLocalMin[i] ) / hEl[i] ); // TODO: Temporarily store the CPDI
          //                                                                                               // corners since they're re-used
          //                                                                                               // below?
          //   }
          // }

          // // for( int corner=0; corner<8; corner++ )
          // // {
          // //   real64 cornerPos[3] ={ 0 };
          // //   for( int w=0; w<3; w++ )
          // //   {
          // //     cornerPos[w] = particlePosition[p][w] + 
          // //                    signs[corner][0] * particleRVectors[p][0][w] + 
          // //                    signs[corner][1] * particleRVectors[p][1][w] + 
          // //                    signs[corner][2] * particleRVectors[p][2][w];                             
          // //   }

          // //   for( int i=0; i<2; i++ )
          // //     {
          // //       for( int j=0; j<2; j++ )
          // //       {
          // //         for( int k=0; k<2; k++ )
          // //         {
          // //           if( cornerIJK[corner][0] + i >= ijkMap.size( 0 ) || cornerIJK[corner][1] + j >= ijkMap.size( 1 ) || cornerIJK[corner][2] + k >= ijkMap.size( 2 ))
          // //           {
          // //             GEOS_LOG_RANK( "Particle: " << p << ", " << 
          // //                            "ijk size={" << ijkMap.size( 0 ) << ", " << ijkMap.size( 1 ) << ", " << ijkMap.size( 2 )  << "}, " << 
          // //                            "cornerIJK=" << cornerIJK[corner][0] << ", " << cornerIJK[corner][1] << ", " << cornerIJK[corner][2] << ", " <<
          // //                            "i=" << i << ", j=" << j << ", k=" << k << ", " << 
          // //                            "particlePos=" << particlePosition[p] << ", " << 
          // //                            "cornerPos={" << cornerPos[0] << ", " << cornerPos[1] << ", " <<  cornerPos[2] << "}" <<
          // //                            "xLocalMin={" << xLocalMin[0] << ", " << xLocalMin[1] << ", " <<  xLocalMin[2] << "}" <<
          // //                            "xLocalMax={" << xLocalMax[0] << ", " << xLocalMax[1] << ", " <<  xLocalMax[2] << "}" << 
          // //                            "InPartition=" << ( partition.isCoordInPartition( particlePosition[p][0], 0 ) && partition.isCoordInPartition( particlePosition[p][1], 1 ) && partition.isCoordInPartition( particlePosition[p][2], 2 ) ) );
          // //                           //  << ", " << 
          // //                           //  "PartitionMin={" << partition.getLocalMin()[0] << ", " << partition.getLocalMin()[1] << ", " << partition.getLocalMin()[2] << "}, " << 
          // //                           //  "PartitionMax={" << partition.getLocalMax()[0] << ", " << partition.getLocalMax()[1] << ", " << partition.getLocalMax()[2] << "}");
          // //           }
          // //         }
          // //       }
          // //     }
          // // }

          // // get node IDs associated with each corner from IJK map, along with weights and grad weights
          // // *** The order in which we access the IJK map must match the order we evaluate the shape functions! ***
          // int node = 0;
          // for( int corner=0; corner<8; corner++ )
          // {
          //   int cornerNode = ijkMap[cornerIJK[corner][0]][cornerIJK[corner][1]][cornerIJK[corner][2]];
          //   // GEOS_LOG_RANK("Particle Corner " << corner << " mapped to corner IJK " << cornerIJK[corner][0] << ", "  << cornerIJK[corner][1] << ", " << cornerIJK[corner][2]);
          //   auto cornerNodePosition = gridPosition[cornerNode];

          //   real64 x, y, z;
          //   x = particlePosition[p][0] + signs[corner][0] * particleRVectors[p][0][0] + signs[corner][1] * particleRVectors[p][1][0] + signs[corner][2] * particleRVectors[p][2][0];
          //   y = particlePosition[p][1] + signs[corner][0] * particleRVectors[p][0][1] + signs[corner][1] * particleRVectors[p][1][1] + signs[corner][2] * particleRVectors[p][2][1];
          //   z = particlePosition[p][2] + signs[corner][0] * particleRVectors[p][0][2] + signs[corner][1] * particleRVectors[p][1][2] + signs[corner][2] * particleRVectors[p][2][2];

          //   real64 xRel = (x - cornerNodePosition[0]) / hEl[0];
          //   real64 yRel = (y - cornerNodePosition[1]) / hEl[1];
          //   real64 zRel = (z - cornerNodePosition[2]) / hEl[2];

          //   for( int i=0; i<2; i++ )
          //   {
          //     real64 xWeight = i * xRel + (1 - i) * (1.0 - xRel);
          //     for( int j=0; j<2; j++ )
          //     {
          //       real64 yWeight = j * yRel + (1 - j) * (1.0 - yRel);
          //       for( int k=0; k<2; k++ )
          //       {
          //         real64 zWeight = k * zRel + (1 - k) * (1.0 - zRel);
          //         real64 weight = xWeight * yWeight * zWeight;

          //         mappedNodes[pp][node] = ijkMap[cornerIJK[corner][0]+i][cornerIJK[corner][1]+j][cornerIJK[corner][2]+k];
          //         shapeFunctionValues[pp][node] = 0.125 * weight;
          //         shapeFunctionGradientValues[pp][node][0] = alpha[corner][0] * weight;
          //         shapeFunctionGradientValues[pp][node][1] = alpha[corner][1] * weight;
          //         shapeFunctionGradientValues[pp][node][2] = alpha[corner][2] * weight;
          //         node++;
          //       }
          //     }
          //   }
          // }
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

void SolidMechanicsMPM::computeSinglePointShapeFunctions( arrayView2d< real64 const > const gridPosition,
                                                          arraySlice1d< real64 const > const particlePosition,
                                                          arrayView3d< int const > const ijkMap,
                                                          real64 const (& xLocalMin)[3],
                                                          real64 const (&hEl)[3],
                                                          arraySlice1d< int > const mappedNodes,
                                                          arraySlice1d< real64 > const shapeFunctionValues,
                                                          arraySlice2d< real64 > const shapeFunctionGradientValues  )
{
  // get IJK associated with particle center
  int centerIJK[3];
  for( int i=0; i<3; i++ )
  {
    centerIJK[i] = floor( ( particlePosition[i] - xLocalMin[i] ) / hEl[i] );
  }

  // get node IDs, weights and grad weights
  int node = 0;
  int corner = ijkMap[centerIJK[0]][centerIJK[1]][centerIJK[2]];
  auto corner_x = gridPosition[corner];

  real64 xRel = (particlePosition[0] - corner_x[0]) / hEl[0];
  real64 yRel = (particlePosition[1] - corner_x[1]) / hEl[1];
  real64 zRel = (particlePosition[2] - corner_x[2]) / hEl[2];

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
        mappedNodes[node] = ijkMap[centerIJK[0]+i][centerIJK[1]+j][centerIJK[2]+k];
        shapeFunctionValues[node] = xWeight * yWeight * zWeight;
        shapeFunctionGradientValues[node][0] = dxWeight * yWeight * zWeight;
        shapeFunctionGradientValues[node][1] = xWeight * dyWeight * zWeight;
        shapeFunctionGradientValues[node][2] = xWeight * yWeight * dzWeight;
        node++;
      }
    }
  }
}


void SolidMechanicsMPM::computeCPDIShapeFunctions( arrayView2d< real64 const > const gridPosition,
                                                   arraySlice1d< real64 const > const particlePosition,
                                                   arraySlice2d< real64 const > const particleRVectors,
                                                   arrayView3d< int const > const ijkMap,
                                                   real64 const (& xLocalMin)[3],
                                                   real64 const (&hEl)[3],
                                                   arraySlice1d< int > const mappedNodes,
                                                   arraySlice1d< real64 > const shapeFunctionValues,
                                                   arraySlice2d< real64 > const shapeFunctionGradientValues )
{
  int const signs[8][3] = { { -1, -1, -1 },
                            {  1, -1, -1 },
                            {  1,  1, -1 },
                            { -1,  1, -1 },
                            { -1, -1,  1 },
                            {  1, -1,  1 },
                            {  1,  1,  1 },
                            { -1,  1,  1 } };
      
  real64 alpha[8][8];
  real64 cpdiVolume, oneOverV;
  real64 p_r1[3], p_r2[3], p_r3[3]; // allowing 1-indexed r-vectors to persist to torture future postdocs >:)

  for( int i=0; i<3; i++ )
  {
    p_r1[i] = particleRVectors[0][i];
    p_r2[i] = particleRVectors[1][i];
    p_r3[i] = particleRVectors[2][i];
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
      real64 cornerPositionComponent = particlePosition[i] + 
                                       signs[corner][0] * particleRVectors[0][i] + 
                                       signs[corner][1] * particleRVectors[1][i] + 
                                       signs[corner][2] * particleRVectors[2][i];                             

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
    x = particlePosition[0] + signs[corner][0] * particleRVectors[0][0] + signs[corner][1] * particleRVectors[1][0] + signs[corner][2] * particleRVectors[2][0];
    y = particlePosition[1] + signs[corner][0] * particleRVectors[0][1] + signs[corner][1] * particleRVectors[1][1] + signs[corner][2] * particleRVectors[2][1];
    z = particlePosition[2] + signs[corner][0] * particleRVectors[0][2] + signs[corner][1] * particleRVectors[1][2] + signs[corner][2] * particleRVectors[2][2];

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

          mappedNodes[node] = ijkMap[cornerIJK[corner][0]+i][cornerIJK[corner][1]+j][cornerIJK[corner][2]+k];
          shapeFunctionValues[node] = 0.125 * weight;
          shapeFunctionGradientValues[node][0] = alpha[corner][0] * weight;
          shapeFunctionGradientValues[node][1] = alpha[corner][1] * weight;
          shapeFunctionGradientValues[node][2] = alpha[corner][2] * weight;
          node++;
        }
      }
    }
  }
}


//CC: Either need to return or pass body force variables by reference
inline void GEOS_DEVICE SolidMechanicsMPM::computeGeneralizedVortexMMSBodyForce( real64 const time_n,
                                                                                 ParticleManager & particleManager )
{
  GEOS_MARK_FUNCTION;

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

    if( constitutiveModelName == "ElasticTransverseIsotropic" || constitutiveModelName == "ElasticTransverseIsotropicPressureDependent" ){
      ElasticTransverseIsotropic & elasticTransverseIsotropic = dynamic_cast< ElasticTransverseIsotropic & >( constitutiveRelation );
      shearModulus = elasticTransverseIsotropic.effectiveShearModulus();
    }

    if( constitutiveModelName == "Graphite" )
    {
      Graphite & graphite = dynamic_cast< Graphite & >( constitutiveRelation );
      shearModulus = graphite.effectiveShearModulus();
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


void SolidMechanicsMPM::computeBodyForce( ParticleManager & particleManager )
{
  GEOS_MARK_FUNCTION;

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

void SolidMechanicsMPM::computeWavespeed( ParticleManager & particleManager )
{
  GEOS_MARK_FUNCTION;

  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    arrayView1d< real64 const> const particleDensity = subRegion.getField< fields::mpm::particleDensity >();
    arrayView1d< real64 > const particleWavespeed = subRegion.getField< fields::mpm::particleWavespeed >();

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

    if( constitutiveModelName == "ElasticTransverseIsotropic" || constitutiveModelName == "ElasticTransverseIsotropicPressureDependent" ){
      const ElasticTransverseIsotropic & elasticTransverseIsotropic = dynamic_cast< const ElasticTransverseIsotropic & >( constitutiveRelation );
      bulkModulus = elasticTransverseIsotropic.effectiveBulkModulus();
      shearModulus = elasticTransverseIsotropic.effectiveShearModulus();
    }

    if( constitutiveModelName == "Graphite" )
    {
      const Graphite & graphite = dynamic_cast< const Graphite & >( constitutiveRelation );
      bulkModulus = graphite.effectiveBulkModulus();
      shearModulus = graphite.effectiveShearModulus();
    }

    forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST ( localIndex const pp ) // would need reduction to parallelize
    {
      localIndex const p = activeParticleIndices[pp];
      particleWavespeed[p] =  sqrt( ( bulkModulus[p] + (4.0/3.0) * shearModulus[p] ) / particleDensity[p] );
    } );
  } );
}

void SolidMechanicsMPM::computeArtificialViscosity( ParticleManager & particleManager )
{
  GEOS_MARK_FUNCTION;

  real64 length = m_planeStrain == 1 ? std::fmin( m_hEl[0], m_hEl[1] ) : std::fmin( m_hEl[0], std::fmin( m_hEl[1], m_hEl[2] ) );

  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    arrayView1d< real64 const > const particleWavespeed = subRegion.getField< fields::mpm::particleWavespeed >();
    arrayView1d< real64 const > const particleDensity = subRegion.getField< fields::mpm::particleDensity >();
    arrayView3d< real64 const > const particleVelocityGradient = subRegion.getField< fields::mpm::particleVelocityGradient >();

    arrayView1d< real64 > const particleArtificialViscosity = subRegion.getField< fields::mpm::particleArtificialViscosity >();

    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
    forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
    {
      localIndex const p = activeParticleIndices[pp];

      // This is computed for du/dx = tr(D) = tr(L);
      real64 trD = LvArray::tensorOps::trace< 3 >( particleVelocityGradient[p] ); // characteristic increment = dx*(du/dx), where du/dx is the characteristic velocity gradient.
      if( trD < 0.0 )
      {
        particleArtificialViscosity[p] = particleDensity[p] * ( m_artificialViscosityQ1 * ( trD * length ) * ( trD * length )
                                         + m_artificialViscosityQ0 * fabs( particleWavespeed[p] * trD * length ) );
      }
      else
      {
        particleArtificialViscosity[p] = 0.0;
      }
    } );
  } );
}


void SolidMechanicsMPM::computeInternalEnergyAndTemperature( const real64 dt,
                                                             ParticleManager & particleManager )
{
  GEOS_MARK_FUNCTION;

  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    // Get constitutive model reference
    // string const & solidMaterialName = subRegion.template getReference< string >( viewKeyStruct::solidMaterialNamesString() );
    // SolidBase & solidModel = getConstitutiveModel< SolidBase >( subRegion, solidMaterialName );
    // if( solidModel.hasWrapper( "strengthScale" ) )
    // {
    // arrayView1d< real64 > const constitutiveStrengthScale = solidModel.getReference< array1d< real64 > >( "strengthScale" );

      arrayView1d< real64 const > const particleDensity = subRegion.getField< fields::mpm::particleDensity >();
      arrayView1d< real64 const > const particleHeatCapacity = subRegion.getField< fields::mpm::particleHeatCapacity >();
      arrayView1d< real64 const > const particleArtificialViscosity = subRegion.getField< fields::mpm::particleArtificialViscosity >();
      arrayView2d< real64 const > const particleStress = subRegion.getField< fields::mpm::particleStress >();
      arrayView3d< real64 const > const particleVelocityGradient = subRegion.getField< fields::mpm::particleVelocityGradient >();
      
      arrayView1d< real64 > const particleTemperature = subRegion.getField< fields::mpm::particleTemperature >();
      arrayView1d< real64 > const particleInternalEnergy = subRegion.getField< fields::mpm::particleInternalEnergy >();

      SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
      forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
      {
        localIndex const p = activeParticleIndices[pp];
        
        //  ========================================================================
        //  Update internal energy with stress power, where Fdot and PK1 stress are a work-
        //  Compute Stress and strain measure The rate of internal energy: udot = (1/rho0)*PK1:Fdot
        //  use L=Fdot*Finv, rho*J=rho0, and PK1=J*sigma(i,k)*Finv(j,k)
        //          udot = (1/rho0)*PK1:Fdot
        //          udot = (1/rho0)*(J*sigma(i,k)*Finv(j,k))*Fdot(i,j)
        //          udot = (1/rho0)*J*sigma(i,k)*L(i,k)
        //          udot = (1/rho)*sigma(i,k)*L(i,k)
        //
        // To compute the temperature it will be necessary to partition the energy into a thermal part
        // and to use some model for the heat capacity.

        real64 energyIncrement = 0.0;

        // 1st half-step internal energy update
        int voigtMap[3][3] = { {0, 5, 4}, {5, 1, 3}, {4, 3, 2} }; // CC: Class member for voigtMap?
        for( int i = 0 ; i < m_numDims ; i++ )
        {
          for( int j = 0 ; j < m_numDims ; j++ )
          {
            energyIncrement += 0.5 * dt * particleStress[p][voigtMap[i][j]] * particleVelocityGradient[p][i][j] / particleDensity[p]; //  General case for increment in internal energy
          }
          energyIncrement -= m_shockHeating * 0.5 * dt * particleArtificialViscosity[p] * particleVelocityGradient[p]( i, i ) / particleDensity[p]; // Artificial viscosity contribution
        }

        particleTemperature[p] += energyIncrement / particleHeatCapacity[p];
        particleInternalEnergy[p] += energyIncrement;

        //      // 1st half-step internal energy update
        //      for( int i = 0 ; i < m_ndim ; i++ )
        //      {
        //        for( int j = 0 ; j < m_ndim ; j++ )
        //        {
        //          p_internalEnergy[pp] += 0.5 * dt * p_stress[pp]( i, j ) * p_L[pp]( i, j ) / p_rho[pp]; //  General case for increment in internal energy
        //        }
        //        p_internalEnergy[pp] -= m_shock_heating * 0.5 * dt * p_q[pp] * p_L[pp]( i, i ) / p_rho[pp]; // Artificial viscosity contribution
        //      }

        // Prevent non-physical negative internal energy from integration error.
        particleInternalEnergy[p] = ( particleInternalEnergy[p] < 0 ) ? 0 : particleInternalEnergy[p];

        // ***************************************************************************
        // TODO: This is the internal energy calculation from Uintah:
        // We need to consider whether to include the artificial viscosity heating.
        //
        //    R2Tensor avgStress = (pStress_new[idx] + pStress[idx])*0.5;
        //    double avgVolume  = (pVolume_deformed[idx]+pVolume[idx])*0.5;
        //
        //    double pSpecificStrainEnergy = (tensorD(0,0)*avgStress(0,0) +
        //                                    tensorD(1,1)*avgStress(1,1) +
        //                                    tensorD(2,2)*avgStress(2,2) +
        //                               2.0*(tensorD(0,1)*avgStress(0,1) +
        //                                    tensorD(0,2)*avgStress(0,2) +
        //                                    tensorD(1,2)*avgStress(1,2)))*
        //                                    avgVolume*delT/pMass[idx];
        //
        //    // Compute rate of change of specific volume
        //    double Vdot = (pVolume_deformed[idx] - pVolume[idx])/(pMass[idx]*delT);
        //
        //    pEnergy_new[idx] = pEnergy[idx] + pSpecificStrainEnergy
        //                                    - p_q[idx]*Vdot*delT*include_AV_heating;
        //
        //    totalStrainEnergy += pSpecificStrainEnergy*pMass[idx];
        // ***************************************************************************

      } );
    // }
  } );

}

void SolidMechanicsMPM::computeSPHJacobian( ParticleManager & particleManager )
{
  GEOS_MARK_FUNCTION;

  // Compute local density based on SPH method (actually compute the local jacobian to account
  // for the possibility of different densities between materials.

  ParticleManager::ParticleViewAccessor< arrayView2d< real64 const > > particlePositionAccessor = particleManager.constructArrayViewAccessor< real64, 2 >( "particleCenter" );
  ParticleManager::ParticleViewAccessor< arrayView1d< real64 const > > particleInitialVolumeAccessor = particleManager.constructArrayViewAccessor< real64, 1 >( "particleInitialVolume" );
  ParticleManager::ParticleViewAccessor< arrayView1d< real64 const > > particleVolumeAccessor = particleManager.constructArrayViewAccessor< real64, 1 >( "particleVolume" );
  ParticleManager::ParticleViewAccessor< arrayView1d< real64 const > > particleMassAccessor = particleManager.constructArrayViewAccessor< real64, 1 >( "particleMass" );

  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    // Get neighbor list
    OrderedVariableToManyParticleRelation & neighborList = subRegion.neighborList();
    arrayView1d< localIndex const > const numNeighborsAll = neighborList.m_numParticles.toViewConst();
    ArrayOfArraysView< localIndex const > const neighborRegions = neighborList.m_toParticleRegion.toViewConst();
    ArrayOfArraysView< localIndex const > const neighborSubRegions = neighborList.m_toParticleSubRegion.toViewConst();
    ArrayOfArraysView< localIndex const > const neighborIndices = neighborList.m_toParticleIndex.toViewConst();

    arrayView2d< real64 const > const particlePosition = subRegion.getParticleCenter();
    arrayView1d< real64 const > const particleVolume = subRegion.getParticleVolume();
    arrayView1d< real64 const > const particleInitialVolume = subRegion.getField< fields::mpm::particleInitialVolume >();
    arrayView1d< real64 const > const particleMass = subRegion.getField< fields::mpm::particleMass >();
    arrayView1d< real64 > const particleSPHJacobian = subRegion.getField< fields::mpm::particleSPHJacobian >();

    ParticleRegion & region = dynamicCast< ParticleRegion & >( subRegion.getParent().getParent() );
    localIndex regionIndexOfSubRegion = region.getIndexInParent();

    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
    forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
    {
      localIndex const p = activeParticleIndices[pp];

      localIndex numNeighbors = numNeighborsAll[p];
      arraySlice1d< localIndex const > const regionIndices = neighborRegions[p];
      arraySlice1d< localIndex const > const subRegionIndices = neighborSubRegions[p];
      arraySlice1d< localIndex const > const particleIndices = neighborIndices[p];

      // Declare and size neighbor data arrays - TODO: switch to std::array? But then we'd need to template computeKernelFieldGradient
      array1d< real64 > neighborVolumes( numNeighbors );
      array2d< real64 > neighborPositions( numNeighbors, 3 );

      // Populate neighbor positions
      for( localIndex neighborIndex = 0; neighborIndex < numNeighbors; neighborIndex++ )
      {
        localIndex regionIndex = regionIndices[neighborIndex];
        localIndex subRegionIndex = subRegionIndices[neighborIndex];
        localIndex particleIndex = particleIndices[neighborIndex];
        neighborPositions[neighborIndex][0] = particlePositionAccessor[regionIndex][subRegionIndex][particleIndex][0];
        neighborPositions[neighborIndex][1] = particlePositionAccessor[regionIndex][subRegionIndex][particleIndex][1];
        neighborPositions[neighborIndex][2] = particlePositionAccessor[regionIndex][subRegionIndex][particleIndex][2];
      }

      // Compare SPH number density to initial particle-wise number density. This gives poor results for
      // non-unform particle spacing and/or local particle refinement.
      if( m_computeSPHJacobian == 1 )
      {
        array1d< real64 > neighborNumberDensities( numNeighbors );

        for( localIndex neighborIndex = 0; neighborIndex < numNeighbors; neighborIndex++ )
        {
          localIndex regionIndex = regionIndices[neighborIndex];
          localIndex subRegionIndex = subRegionIndices[neighborIndex];
          localIndex particleIndex = particleIndices[neighborIndex];
          neighborNumberDensities[neighborIndex] = ( m_planeStrain ? m_hEl[2] : 1.0 ) / particleVolumeAccessor[regionIndex][subRegionIndex][particleIndex];
          neighborVolumes[neighborIndex] = particleVolumeAccessor[regionIndex][subRegionIndex][particleIndex];
        }

        // TODO modify this function to mirror points near symmetry planes (within neighbor radius of the boundary)
        // to obtain good behavior at sym boundaries.

        // Evaluate kernel field
        real64 sphParticleNumberDensity = computeKernelField( particlePosition[p], // query point
                                                              neighborPositions, // List of neighbor particle locations
                                                              neighborVolumes,
                                                              neighborNumberDensities );

        // Determine J
        real64 initialParticleNumberDensity;
        if( m_planeStrain )
        {
          // This is actually the particles per unit area.  This assumes that there is only one particles in
          // the z-direction, which will only be the case if the particleFileWriter script is called
          // with the python flag planeStrain=1 (in addition to setting that GEOS flag)
          initialParticleNumberDensity = m_hEl[2] / particleInitialVolume[p];
        }
        else
        {
          initialParticleNumberDensity = 1.0 / particleInitialVolume[p];
        }
        particleSPHJacobian[p] = initialParticleNumberDensity / fmax( 1.0e-9, sphParticleNumberDensity );
      }
      // Use SPH mass density to determine J. Only neighbors that are the same material as the ppth particle are considered.
      // This approach is more robust against non-uniform spacing/local particle refinement.
      else if ( m_computeSPHJacobian == 2 )
      {
        array1d< real64 > neighborMassDensities( numNeighbors );

        for( localIndex neighborIndex = 0; neighborIndex < numNeighbors; neighborIndex++ )
        {
          localIndex regionIndex = regionIndices[neighborIndex];
          if( regionIndex == regionIndexOfSubRegion ) // same material 
          {
            localIndex subRegionIndex = subRegionIndices[neighborIndex];
            localIndex particleIndex = particleIndices[neighborIndex];
            neighborMassDensities[neighborIndex] = particleMassAccessor[regionIndex][subRegionIndex][particleIndex] / particleVolumeAccessor[regionIndex][subRegionIndex][particleIndex]; // CC: Should we just used particle density directly?
            neighborVolumes[neighborIndex] = particleVolumeAccessor[regionIndex][subRegionIndex][particleIndex];
          }
          else
          {
            neighborMassDensities[neighborIndex] = 0.0;
            neighborVolumes[neighborIndex] = 0.0;
          }
        }

        // Evaluate kernel field
        real64 sphMassDensity = computeKernelField( particlePosition[p], // query point
                                                    neighborPositions,    // List of neighbor particle locations
                                                    neighborVolumes,
                                                    neighborMassDensities );

        particleSPHJacobian[p] = ( particleMass[p] / particleInitialVolume[p] ) / fmax( 1.0e-9, sphMassDensity );
      }
      else
      {
        particleSPHJacobian[p] = particleVolume[p]/particleInitialVolume[p];
      }
    } );
  } );
}


void SolidMechanicsMPM::overlapCorrection( real64 const dt,
                                           ParticleManager & particleManager )
{
  GEOS_MARK_FUNCTION;

  // Compute the grid field gradient as the max of the gradients of the particles mapping to the node
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    arrayView3d< real64 > particleDeformationGradient = subRegion.getField< fields::mpm::particleDeformationGradient >();
    arrayView3d< real64 > particleFDot = subRegion.getField< fields::mpm::particleFDot >();
    arrayView3d< real64 > particleVelocityGradient = subRegion.getField< fields::mpm::particleVelocityGradient >();
    arrayView1d< real64 const > const particleSPHJacobian = subRegion.getField< fields::mpm::particleSPHJacobian >();

    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
    forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
    {
      localIndex const p = activeParticleIndices[pp];

      real64 Fold[3][3] = { { 0 } };
      LvArray::tensorOps::copy< 3, 3 >( Fold, particleDeformationGradient[p] );

      // Update the Jacobian based on particle data.
      real64 J = LvArray::tensorOps::determinant< 3 >( Fold );

      if ( J <= 0.0 )
      {
        return; // CC: For lambda we would return instead of continue?
      }

      // If there is overdensification, the jacobian as computed from the sph kernel will be much less
      // than that from the particle def. grad, so overlap>1 may indicate overdensification. But, since
      // the neighbor radius is small, there is noise in J_sph, so we only want to scale F if the
      // overlap is above a threshold.
      real64 overlap = J / particleSPHJacobian[p];

      // If the densification exceeds 1, slowly transition from the particle J to the non-local SPH J,
      // and scale Fij accordingly.  Perform this so that J=Jsph if J/J_sph exceeds a threshold.
      // Note that at free surfaces, undeformed Jsph will be >1. This is also true at sym boundaries
      // until we can fix the calculation of the kernel field at boundaries.
      if( overlap > m_overlapThreshold1 )
      {
        real64 alpha = 1.0; // Smoothing parameter: Jnew = (1-alpha)*Jold + alpha*Jsph
        if( overlap <= m_overlapThreshold2 )
        {
          // smooth step function
          real64 xi = ( overlap - m_overlapThreshold1 ) / ( m_overlapThreshold2 - m_overlapThreshold1 );
          alpha = 3.0 * xi * xi - 2.0 * xi * xi * xi;
        }

        real64 scale;
        if( m_planeStrain == 1 )
        {
          scale = std::sqrt( 1.0 + alpha * ( 1.0 / overlap - 1.0 ) );
        }
        else
        {
          scale = std::pow( 1.0 + alpha * ( 1.0 / overlap - 1.0 ), 1.0 / 3.0 );
        }

        // put a limit on scale so it can't change the density too much in a single time-step.
        scale = std::fmax( 0.99, scale );

        // TODO modify p_L as well, consistent with the change to p_F in case a hypoelastic constitutive model
        // is used, and so the update to the internal energy is consistent with the change in F.
        //
        // The simple attempt, commented out below, produces extremely high internal energy when
        // the overlap correction is used in a DFG brittle damage problem.

        for( int i = 0 ; i < m_numDims; i++ )
        {
          for( int j = 0 ; j < m_numDims; j++ )
          {
            particleDeformationGradient[p][i][j] *= scale;
            particleFDot[p][i][j] += ( particleDeformationGradient[p][i][j] - Fold[i][j] ) / dt;

            // Update velocity gradient without taking the inverse of the particle deformation gradient
            particleVelocityGradient[p][i][j] /= scale;
          }

          particleVelocityGradient[p][i][i] += ( scale - 1.0 ) / ( scale * dt );
        }
        // Modify p_L as well, consistent with the change to p_F in case a hypoelastic constitutive model
        // is used, and so the update to the internal energy is consistent with the change in F.
        
        // real64 invF[3][3] = { { 0 } };
        // LvArray::tensorOps::invert< 3 >( invF, particleDeformationGradient[p] );
        // LvArray::tensorOps::Rij_eq_AikBkj< 3, 3, 3 >( particleVelocityGradient[p], particleFDot[p], invF );
      }
    } );
  } );

}

void SolidMechanicsMPM::resetDeformationGradient( ParticleManager & particleManager )
{
  // Reset the deformation gradient to be the spherical part.
  // This should only be used for cases where the deviatoric part
  // of the constitutive model is hypoelastic.
  //
  // We also keep to rotation to make the plotting look a bit better.
  //
  // We also use this for the gas constitutive model.
  // TODO: Make this a more general option (per material type)
  
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    string const & solidMaterialName = subRegion.template getReference< string >( viewKeyStruct::solidMaterialNamesString() );
    const SolidBase & solidModel = getConstitutiveModel< SolidBase >( subRegion, solidMaterialName );
    string constitutiveModelName = solidModel.getCatalogName();

    if( constitutiveModelName == "Air" )
    {
      arrayView3d< real64 > const particleDeformationGradient = subRegion.getField< fields::mpm::particleDeformationGradient >();
      arrayView1d< real64 const > const particleDamage = subRegion.getParticleDamage();

      SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
      forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
      {
        localIndex const p = activeParticleIndices[pp];

        if ( particleDamage[p] > 0.9999999 )
        {
          real64 rotation[3][3] = { { 0 } }; 
          real64 deformationGradient[3][3] = { { 0 } };
          LvArray::tensorOps::copy< 3, 3 >( deformationGradient, particleDeformationGradient[p] );
          polarDecomposition( rotation, deformationGradient ); // Start-of-step polar decomposition
          
          real64 J = LvArray::tensorOps::determinant< 3 >( particleDeformationGradient[p]);

          real64 U[3][3] = { { 0 } };
          if ( m_planeStrain == 1 )
          {
            real64 JtoOneHalf = std::sqrt( J );
            U[0][0] = JtoOneHalf;
            U[1][1] = JtoOneHalf;
            U[2][2] = 1.0;
          }
          else
          {
            real64 JtoOneThird = std::pow( J, 1.0 / 3.0 );
            U[0][0] = JtoOneThird;
            U[1][1] = JtoOneThird;
            U[2][2] = JtoOneThird;
          }

          LvArray::tensorOps::Rij_eq_AikBkj< 3, 3, 3 >( particleDeformationGradient[p], rotation, U);
        }

      } );
    }
  } );
}

void SolidMechanicsMPM::unscaleCPDIVectors( ParticleManager & particleManager )
{
  GEOS_MARK_FUNCTION;

  // This undoes CPDI domain scaling. Can be used to plot unscaled domains.
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    // Only perform on CPDI particles
    if( subRegion.getParticleType() == ParticleType::CPDI )
    {
      arrayView3d< real64 const > const particleDeformationGradient = subRegion.getField< fields::mpm::particleDeformationGradient >();
      arrayView3d< real64 const > const particleInitialRVectors = subRegion.getField< fields::mpm::particleInitialRVectors >();
      arrayView3d< real64 > const particleRVectors = subRegion.getParticleRVectors();

      SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
      forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
      {
        localIndex const p = activeParticleIndices[pp];
        // for( int i=0; i < m_numDims; i++)
        // {
          LvArray::tensorOps::Ri_eq_AijBj< 3, 3 >( particleRVectors[p][0], particleDeformationGradient[p], particleInitialRVectors[p][0] );
          LvArray::tensorOps::Ri_eq_AijBj< 3, 3 >( particleRVectors[p][1], particleDeformationGradient[p], particleInitialRVectors[p][1] );
          LvArray::tensorOps::Ri_eq_AijBj< 3, 3 >( particleRVectors[p][2], particleDeformationGradient[p], particleInitialRVectors[p][2] );  
        // }
      } );
    }
  } );
}

void SolidMechanicsMPM::computeKineticEnergy( ParticleManager & particleManager )
{
  GEOS_MARK_FUNCTION;

  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    arrayView1d< real64 > const particleKineticEnergy = subRegion.getField< fields::mpm::particleKineticEnergy >();
    arrayView1d< real64 const > const particleMass = subRegion.getField< fields::mpm::particleMass >();
    arrayView2d< real64 const > const particleVelocity = subRegion.getParticleVelocity();

    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
    forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
    {
      localIndex const p = activeParticleIndices[pp];
      // TODO: add micro kinetic energy term based on L?
      particleKineticEnergy[p] = 0.5 * particleMass[p] * LvArray::tensorOps::l2NormSquared< 3 >( particleVelocity[p] );
    } );
  } );
}

// This function will do sums of nodal data to write a profile to a file of the
// state averaged across the y-z plane at each x-nodal position.  This is very slow
// but is useful for on-the-fly extraction of shock profiles from mesoscale simulations.
// Typically this is not done every time step.
// Sum local arrays to compute x-profile to write to file
void SolidMechanicsMPM::computeXProfile( int const cycleNumber,
                                         real64 const time,
                                         real64 const dt,
                                         NodeManager & nodeManager,
                                         SpatialPartition & partition )
{
  GEOS_MARK_FUNCTION;

  int numXNodesLocal = m_nEl[0] + 1;
  // int numYNodesLocal = m_nEl[0] + 1;
  // int numZNodesLocal = m_nEl[0] + 1;

  int numXNodesGlobal = static_cast< int >( round( ( m_xGlobalMax[0] - m_xGlobalMin[0]) / m_hEl[0] ) ) + 1;
  int numYNodesGlobal = static_cast< int >( round( ( m_xGlobalMax[1] - m_xGlobalMin[1]) / m_hEl[1] ) ) + 1;
  int numZNodesGlobal = static_cast< int >( round( ( m_xGlobalMax[2] - m_xGlobalMin[2]) / m_hEl[2] ) ) + 1;

  array1d< real64 > xProfileDensityLocal( numXNodesGlobal ), 
                    xProfileDensityGlobal( numXNodesGlobal ),
                    xProfileDamageLocal( numXNodesGlobal ), 
                    xProfileDamageGlobal( numXNodesGlobal ),
                    xProfileXStressLocal( numXNodesGlobal ), 
                    xProfileXStressGlobal( numXNodesGlobal ),
                    xProfileYStressLocal( numXNodesGlobal ), 
                    xProfileYStressGlobal( numXNodesGlobal ),
                    xProfileZStressLocal( numXNodesGlobal ), 
                    xProfileZStressGlobal( numXNodesGlobal ),
                    xProfileKineticEnergyLocal( numXNodesGlobal ), 
                    xProfileKineticEnergyGlobal( numXNodesGlobal );

  for( int n=0; n < numXNodesGlobal; n++ )
  {
    xProfileDamageLocal[n] = 0.0;
    xProfileDamageGlobal[n] = 0.0;

    xProfileDensityLocal[n] = 0.0;
    xProfileDensityGlobal[n] = 0.0;

    xProfileXStressLocal[n] = 0.0;
    xProfileXStressGlobal[n] = 0.0;

    xProfileYStressLocal[n] = 0.0;
    xProfileYStressGlobal[n] = 0.0;

    xProfileZStressLocal[n] = 0.0;
    xProfileZStressGlobal[n] = 0.0;

    xProfileKineticEnergyLocal[n] = 0.0;
    xProfileKineticEnergyGlobal[n] = 0.0;
  }

  // array will store density at each x-nodal position, averaged over the domain.
  // Compute the damage, which is stored in the mass-weighted variable
  // g_damage = (sum m_p*D_p)/m_i for each field (which also assumes p_flag=2 equates to D=1)
  // dividing by the summed mass should give mass fractino of damaged material.
  
  // array will store density at each x-nodal position, averaged over the domain.
  // Compute density first by summing mass across all nodes and fields.
  // set profile to 0 before sum (may not be needed).

  // index in the global x-position array of the current node.
  int iGlobal, 
      jGlobal, 
      kGlobal,
      iMin, 
      iMax,
      jMin, 
      jMax,
      kMin,
      kMax;

  real64 smallX = 1.0e-8*(m_xLocalMax[0] - m_xLocalMin[0]);
  real64 smallY = 1.0e-8*(m_xLocalMax[1] - m_xLocalMin[1]);
  real64 smallZ = 1.0e-8*(m_xLocalMax[2] - m_xLocalMin[2]);

  // limits are defined so that a node is counted if imin <= i <=imax
  iMin = fabs(m_xLocalMin[0] - m_xGlobalMin[0]) < smallX ? 0 : round( ( m_xLocalMin[0] - m_xGlobalMin[0] ) / m_hEl[0] ) + 2;
  iMax = fabs(m_xLocalMax[0] - m_xGlobalMax[0]) < smallX ? numXNodesGlobal : round( ( m_xLocalMax[0] - m_xGlobalMin[0] ) / m_hEl[0] ) - 1;

  jMin = fabs(m_xLocalMin[1] - m_xGlobalMin[1]) < smallY ? 0 : round( ( m_xLocalMin[1] - m_xGlobalMin[1] ) / m_hEl[1] ) + 2;
  jMax = fabs(m_xLocalMax[1] - m_xGlobalMax[1]) < smallY ? numYNodesGlobal : round( ( m_xLocalMax[1] - m_xGlobalMin[1] ) / m_hEl[1] ) - 1;

  kMin = fabs(m_xLocalMin[2] - m_xGlobalMin[2]) < smallZ ? 0 : round( ( m_xLocalMin[2] - m_xGlobalMin[2] ) / m_hEl[2] ) + 2;
  kMax = fabs(m_xLocalMax[2] - m_xGlobalMax[2]) < smallZ ? numZNodesGlobal : round( ( m_xLocalMax[2] - m_xGlobalMin[2] ) / m_hEl[2] ) - 1;


  arrayView1d< int const > const periodic = partition.getPeriodic();
  // int numNodes = nodeManager.size();
  arrayView2d< real64, nodes::REFERENCE_POSITION_USD > const & gridPosition = nodeManager.referencePosition();
  arrayView2d< real64 > const gridMass = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::massString() );
  arrayView2d< real64 > const gridMassWeightedDamage = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::massWeightedDamageString() );
  arrayView3d< real64 > const gridNormalStress = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::normalStressString() );
  arrayView3d< real64 > const gridVelocity = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::velocityString() );

  // sum all values to local array.
  for( localIndex g = 0 ; g < numXNodesLocal; g++ )
  {
    iGlobal = round( ( gridPosition[g][0] - m_xGlobalMin[0] ) / m_hEl[0] ); // index in the global x-position array of the current node.
    jGlobal = round( ( gridPosition[g][1] - m_xGlobalMin[1] ) / m_hEl[1] ); // index in the global y-position array of the current node.
    kGlobal = round( ( gridPosition[g][2] - m_xGlobalMin[2] ) / m_hEl[2] ); // index in the global z-position array of the current node.



    for( localIndex fieldIndex = 0; fieldIndex < m_numVelocityFields; ++fieldIndex )
    {
      // Avoid double counting ghost nodes for mass, which is being taken from a synced field..
      if( iGlobal >= iMin && iGlobal <= iMax && 
          jGlobal >= jMin && jGlobal <= (jMax-periodic[1]) && 
          kGlobal >= kMin && kGlobal <= (kMax-periodic[2]) )
      {
        // CC: Moved this inside fieldIndex loop, NEED TO DOUBLE CHECK!
        // mass weighted damage is computed from master particles and never synchronized between patches
        // so a simple sum is adequate

        xProfileDamageLocal[iGlobal] += gridMassWeightedDamage[g][fieldIndex];

        // Normal stress components are mapped to nodes from master particles and never syncronized between patches,
        // so a simple sum should be adequate
        xProfileXStressLocal[iGlobal] += gridNormalStress[g][fieldIndex][0];
        xProfileYStressLocal[iGlobal] += gridNormalStress[g][fieldIndex][1];
        xProfileZStressLocal[iGlobal] += gridNormalStress[g][fieldIndex][2];

        // Kinetic energy is computed relative to some initial velocity ( v0=[vx0,0,0] ) whre vx0 is input.  This allows
        // for micro kinetic energy profiles to be computed for reverse ballistic simulations.
        xProfileKineticEnergyLocal[iGlobal] += 0.5 * gridMass[g][fieldIndex] * (
            ( gridVelocity[g][fieldIndex][0] - m_xProfileVx0 ) * ( gridVelocity[g][fieldIndex][0] - m_xProfileVx0 ) +
            gridVelocity[g][fieldIndex][1] * gridVelocity[g][fieldIndex][1] +
            gridVelocity[g][fieldIndex][2] * gridVelocity[g][fieldIndex][2] );

        xProfileDensityLocal[iGlobal] += gridMass[g][fieldIndex];
      }
    } // end loop over field index
  } // end loop over grid nodes

  // do additive sync on global mass-eighted damage array (not yet divided by mass)
  MpiWrapper::allReduce( xProfileDamageLocal.data(),
                         xProfileDamageGlobal.data(),
                         xProfileDamageLocal.size(),
                         MpiWrapper::getMpiOp( MpiWrapper::Reduction::Sum ),
                         MPI_COMM_GEOSX );

  // do additive sync on global mass array (e.g. mass)
  MpiWrapper::allReduce( xProfileDensityLocal.data(),
                         xProfileDensityGlobal.data(),
                         xProfileDensityLocal.size(),
                         MpiWrapper::getMpiOp( MpiWrapper::Reduction::Sum ),
                         MPI_COMM_GEOSX );

  // do additive sync on global normalStress arrays (e.g. summed stress*volume)
  MpiWrapper::allReduce( xProfileXStressLocal.data(),
                         xProfileXStressGlobal.data(),
                         xProfileXStressLocal.size(),
                         MpiWrapper::getMpiOp( MpiWrapper::Reduction::Sum ),
                         MPI_COMM_GEOSX );

  MpiWrapper::allReduce( xProfileYStressLocal.data(),
                         xProfileYStressGlobal.data(),
                         xProfileYStressLocal.size(),
                         MpiWrapper::getMpiOp( MpiWrapper::Reduction::Sum ),
                         MPI_COMM_GEOSX );

  MpiWrapper::allReduce( xProfileZStressLocal.data(),
                         xProfileZStressGlobal.data(),
                         xProfileZStressLocal.size(),
                         MpiWrapper::getMpiOp( MpiWrapper::Reduction::Sum ),
                         MPI_COMM_GEOSX );

  // do additive sync on global kinetic energy array (summed 0.5*m*V*V)
  MpiWrapper::allReduce( xProfileKineticEnergyLocal.data(),
                         xProfileKineticEnergyGlobal.data(),
                         xProfileKineticEnergyLocal.size(),
                         MpiWrapper::getMpiOp( MpiWrapper::Reduction::Sum ),
                         MPI_COMM_GEOSX );

  // Divide the summed mass-weighted damage by the summed nodal mass (the latter is still stored
  // in the not-yet-normalized-by-volume "density" array.  Divide element-wise by mass, but
  // avoid division by zero for 0-mass nodes.
  //
  // The resulting array is now the mass fraction of damaged material at each x-nodal position.
  for( int i = 0 ; i < numXNodesGlobal; i++ )
  {
    xProfileDamageGlobal[i] = ( xProfileDensityGlobal[i] > m_smallMass ) ? ( xProfileDamageGlobal[i] / xProfileDensityGlobal[i] ) : 0.0;
  }

  // Compute cross section of global domain and divide to convert summed mass into density profile.
  // This is global domain cross section in y-z plane, excluding ghost cells.
  real64 volume = m_hEl[0] * ( m_xGlobalMax[1] - m_xGlobalMin[1] - (2.0 + periodic[1]) * m_hEl[1] ) * ( m_xGlobalMax[2] - m_xGlobalMin[2] - (2.0 + periodic[2]) * m_hEl[2] );

  for( localIndex g = 0 ; g < numXNodesLocal; g++ )
  {
    xProfileDensityGlobal[g] *= 1.0 / volume;
    xProfileXStressGlobal[g] *= 1.0 / volume;
    xProfileYStressGlobal[g] *= 1.0 / volume;
    xProfileZStressGlobal[g] *= 1.0 / volume;
    xProfileKineticEnergyGlobal[g] *= 1.0 / volume;
  }

  // proc 0 writes the profile data.
  int rank = 0;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  if( rank == 0 )
  {
    // ---------------------------------------------
    // Mass Profile
    //
    std::ofstream file;
    //can't enable exception now because of gcc bug that raises ios_base::failure with useless message
    //file.exceptions(file.exceptions() | std::ios::failbit);
    file.open( "x_profile_density.csv", std::ios::out | std::ios::app );
    if( file.fail() )
      throw std::ios_base::failure( std::strerror( errno ) );
    //make sure write fails with exception if something is wrong
    file.exceptions( file.exceptions() | std::ios::failbit | std::ifstream::badbit );

    if( cycleNumber == 0 )
    {
      // Write the header to the file.
      // time, x_0, x_1, ..., x_n

      file << "time";
      for( int i = 0 ; i < ( numXNodesGlobal ) ; i++ )
      {
        file << std::setprecision( std::numeric_limits<long double>::digits10 ) << "," << m_xGlobalMin[0] + i * m_hEl[0];
      }
      file << std::endl;
    }

    // Write the current time and the nodal vaerage for each x position.

    file << std::setprecision( std::numeric_limits<long double>::digits10 ) << time + dt;
    for( int i = 0 ; i < ( numXNodesGlobal ) ; i++ )
    {
      file << std::setprecision( std::numeric_limits<long double>::digits10 ) << "," << xProfileDensityGlobal[i];
    }
    file << std::endl;
    file.close();

    // ---------------------------------------------
    // Damage Profile
    //
    //can't enable exception now because of gcc bug that raises ios_base::failure with useless message
    //file.exceptions(file.exceptions() | std::ios::failbit);
    file.open( "x_profile_damage.csv", std::ios::out | std::ios::app );
    if( file.fail() )
      throw std::ios_base::failure( std::strerror( errno ) );
    //make sure write fails with exception if something is wrong
    file.exceptions( file.exceptions() | std::ios::failbit | std::ifstream::badbit );

    if( cycleNumber == 0 )
    {
      // Write the header to the file.
      // time, x_0, x_1, ..., x_n

      file << "time";
      for( int i = 0 ; i < ( numXNodesGlobal ) ; i++ )
      {
        file << std::setprecision( std::numeric_limits<long double>::digits10 ) << "," << m_xGlobalMin[0] + i * m_hEl[0];
      }
      file << std::endl;
    }

    // Write the current time and the nodal vaerage for each x position.

    file << std::setprecision( std::numeric_limits<long double>::digits10 ) << time + dt;
    for( int i = 0 ; i < ( numXNodesGlobal ) ; i++ )
    {
      file << std::setprecision( std::numeric_limits<long double>::digits10 ) << "," << xProfileDamageGlobal[i];
    }
    file << std::endl;
    file.close();

    // ---------------------------------------------
    // xStress Profile
    //

    //can't enable exception now because of gcc bug that raises ios_base::failure with useless message
    //file.exceptions(file.exceptions() | std::ios::failbit);
    file.open( "x_profile_xStress.csv", std::ios::out | std::ios::app );
    if( file.fail() )
      throw std::ios_base::failure( std::strerror( errno ) );
    //make sure write fails with exception if something is wrong
    file.exceptions( file.exceptions() | std::ios::failbit | std::ifstream::badbit );

    if( cycleNumber == 0 )
    {
      // Write the header to the file.
      // time, x_0, x_1, ..., x_n

      file << "time";
      for( int i = 0 ; i < ( numXNodesGlobal ) ; i++ )
      {
        file << std::setprecision( std::numeric_limits<long double>::digits10 ) << "," << m_xGlobalMin[0] + i * m_hEl[0];
      }
      file << std::endl;
    }

    // Write the current time and the nodal average for each x position.

    file << std::setprecision( std::numeric_limits<long double>::digits10 ) << time + dt;
    for( int i = 0 ; i < ( numXNodesGlobal ) ; i++ )
    {
      file << std::setprecision( std::numeric_limits<long double>::digits10 ) << "," << xProfileXStressGlobal[i];
    }
    file << std::endl;
    file.close();

    // ---------------------------------------------
    // yStress Profile
    //

    //can't enable exception now because of gcc bug that raises ios_base::failure with useless message
    //file.exceptions(file.exceptions() | std::ios::failbit);
    file.open( "x_profile_yStress.csv", std::ios::out | std::ios::app );
    if( file.fail() )
      throw std::ios_base::failure( std::strerror( errno ) );
    //make sure write fails with exception if something is wrong
    file.exceptions( file.exceptions() | std::ios::failbit | std::ifstream::badbit );

    if( cycleNumber == 0 )
    {
      // Write the header to the file.
      // time, x_0, x_1, ..., x_n

      file << "time";
      for( int i = 0 ; i < ( numXNodesGlobal ) ; i++ )
      {
        file << std::setprecision( std::numeric_limits<long double>::digits10 ) << "," << m_xGlobalMin[0] + i * m_hEl[0];
      }
      file << std::endl;
    }

    // Write the current time and the nodal average for each x position.

    file << std::setprecision( std::numeric_limits<long double>::digits10 ) << time + dt;
    for( int i = 0 ; i < numXNodesGlobal ; i++ )
    {
      file << std::setprecision( std::numeric_limits<long double>::digits10 ) << "," << xProfileYStressGlobal[i];
    }
    file << std::endl;
    file.close();

    // ---------------------------------------------
    // zStress Profile
    //

    //can't enable exception now because of gcc bug that raises ios_base::failure with useless message
    //file.exceptions(file.exceptions() | std::ios::failbit);
    file.open( "x_profile_zStress.csv", std::ios::out | std::ios::app );
    if( file.fail() )
      throw std::ios_base::failure( std::strerror( errno ) );
    //make sure write fails with exception if something is wrong
    file.exceptions( file.exceptions() | std::ios::failbit | std::ifstream::badbit );

    if( cycleNumber == 0 )
    {
      // Write the header to the file.
      // time, x_0, x_1, ..., x_n

      file << "time";
      for( int i = 0 ; i < ( numXNodesGlobal ) ; i++ )
      {
        file << std::setprecision( std::numeric_limits<long double>::digits10 ) << "," << m_xGlobalMin[0] + i * m_hEl[0];
      }
      file << std::endl;
    }

    // Write the current time and the nodal average for each x position.

    file << std::setprecision( std::numeric_limits<long double>::digits10 ) << time + dt;
    for( int i = 0 ; i < ( numXNodesGlobal ) ; i++ )
    {
      file << std::setprecision( std::numeric_limits<long double>::digits10 ) << "," << xProfileZStressGlobal[i];
    }
    file << std::endl;
    file.close();

    // ---------------------------------------------
    // KINETICENERGY Profile
    //

    //can't enable exception now because of gcc bug that raises ios_base::failure with useless message
    //file.exceptions(file.exceptions() | std::ios::failbit);
    file.open( "x_profile_kineticEnergy.csv", std::ios::out | std::ios::app );
    if( file.fail() )
      throw std::ios_base::failure( std::strerror( errno ) );
    //make sure write fails with exception if something is wrong
    file.exceptions( file.exceptions() | std::ios::failbit | std::ifstream::badbit );

    if( cycleNumber == 0 )
    {
      // Write the header to the file.
      // time, x_0, x_1, ..., x_n

      file << "time";
      for( int i = 0 ; i < numXNodesGlobal; i++ )
      {
        file << std::setprecision( std::numeric_limits<long double>::digits10 ) << "," << m_xGlobalMin[0] + i * m_hEl[0];
      }
      file << std::endl;
    }

    // Write the current time and the nodal average for each x position.
    file << std::setprecision( std::numeric_limits<long double>::digits10 ) << time + dt;
    for( int i = 0 ; i < numXNodesGlobal; i++ )
    {
      file << std::setprecision( std::numeric_limits<long double>::digits10 ) << "," << xProfileKineticEnergyGlobal[i];
    }
    file << std::endl;
    file.close();
  }
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

// A helper function to calculate polar decomposition. TODO: Previously this was an LvArray method, hopefully it will be again someday.
GEOS_HOST_DEVICE
void polarDecomposition( real64 (& R)[3][3],
                         real64 const (&matrix)[3][3] )
{
  // Initialize
  LvArray::tensorOps::copy< 3, 3 >( R, matrix );
  real64 RInverse[3][3] = { {0} },
         RInverseTranspose[3][3] = { {0} },
         RRTMinusI[3][3] = { {0} };

  // Higham Algorithm
  real64 errorSquared = 1.0;
  real64 tolerance = 10 * LvArray::NumericLimits< real64 >::epsilon;
  int iter = 0;
  while( errorSquared > tolerance * tolerance && iter < 100 )
  {
    iter++;
    errorSquared = 0.0;

    // Average the current R with its inverse tranpose
    LvArray::tensorOps::internal::SquareMatrixOps< 3 >::invert( RInverse, R );
    LvArray::tensorOps::transpose< 3, 3 >( RInverseTranspose, RInverse );
    LvArray::tensorOps::add< 3, 3 >( R, RInverseTranspose );
    LvArray::tensorOps::scale< 3, 3 >( R, 0.5 );

    // Determine how close R is to being orthogonal using L2Norm(R.R^T-I)
    real64 copyR[3][3];
    LvArray::tensorOps::copy< 3, 3 >( copyR, R );
    LvArray::tensorOps::Rij_eq_AikBjk< 3, 3, 3 >( RRTMinusI, R, copyR );
    LvArray::tensorOps::addIdentity< 3 >( RRTMinusI, -1.0 );
    for( std::ptrdiff_t i = 0; i < 3; i++ )
    {
      for( std::ptrdiff_t j = 0; j < 3; j++ )
      {
        errorSquared += RRTMinusI[i][j] * RRTMinusI[i][j];
      }
    }
  }
  if( iter == 100 )
  {
    GEOS_LOG_RANK( "Polar decomposition did not converge in 100 iterations!" );
  }
}

REGISTER_CATALOG_ENTRY( SolverBase, SolidMechanicsMPM, string const &, dataRepository::Group * const )
}
