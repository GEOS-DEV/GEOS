/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file SolidMechanicsMPM.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_SOLIDMECHANICS_MPM_HPP_
#define GEOS_PHYSICSSOLVERS_SOLIDMECHANICS_MPM_HPP_

#include "LvArray/src/tensorOps.hpp"

#include "physicsSolvers/PhysicsSolverBase.hpp"

#include "common/TimingMacros.hpp"
#include "kernels/ExplicitMPM.hpp"

#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "mesh/mpiCommunications/MPI_iCommData.hpp"

#include "events/mpmEvents/MPMEventManager.hpp"
#include "mesh/CohesiveZoneManager.hpp"

#include "MPMSolverFields.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsFields.hpp"
#include "physicsSolvers/solidMechanics/MPMSolverEnums.hpp"

namespace geos
{
namespace axisymmetricSurface
{
struct Support;
}

class SpatialPartition;
class SetInitialTemperatureAndPressureMPMEvent;

/**
 * @brief Explicit material point method solver for deforming solids.
 *
 * Supports point, B-spline, CPDI, and CPTI particle mappings, multi-field contact,
 * cohesive interfaces, and Cartesian or axisymmetric kinematics.
 */
class SolidMechanicsMPM : public PhysicsSolverBase
{
public:
  /**
   * @brief Constructs an MPM solver in the supplied repository group.
   * @param name Solver instance name.
   * @param parent Owning parent group.
   */
  SolidMechanicsMPM( const string & name,
                     Group * const parent );

  SolidMechanicsMPM( SolidMechanicsMPM const & ) = delete;
  SolidMechanicsMPM( SolidMechanicsMPM && ) = default;

  SolidMechanicsMPM & operator=( SolidMechanicsMPM const & ) = delete;
  SolidMechanicsMPM & operator=( SolidMechanicsMPM && ) = delete;

  /// @brief Destroys the solver instance.
  virtual ~SolidMechanicsMPM() override {}

  /// @brief Returns the registered solver catalog key.
  static string catalogName() { return "SolidMechanics_MPM"; }

  /// @copydoc PhysicsSolverBase::getCatalogName()
  string getCatalogName() const override { return catalogName(); }

  /// @brief MPM computes particle interpolation directly, without cell-quadrature gradients.
  bool requiresCellShapeFunctionGradients() const override { return false; }

  /// @brief Creates solver-owned event and cohesive-zone subgroups before input registration.
  virtual void initializePreSubGroups() override;

  /// @brief Registers solver particle, nodal, and constitutive fields on the mesh.
  virtual void registerDataOnMesh( Group & meshBodies ) override final;

  /** @name Time stepping and event interface */
  /** @{ */

  /// @brief Runs one solver time step and returns the time increment selected for continuation.
  virtual
  real64 solverStep( real64 const & time_n,
                     real64 const & dt,
                     integer const cycleNumber,
                     DomainPartition & domain ) override;

  /// @brief Advances the explicit MPM sequence, including mapping, contact, constitutive updates, and output.
  virtual
  real64 explicitStep( real64 const & time_n,
                       real64 const & dt,
                       integer const cycleNumber,
                       DomainPartition & domain ) override;

  /// @brief No separate state update is required outside the explicit step.
  virtual void updateState( DomainPartition & domain ) override final
  {
    // All state changes occur during the explicit step.
    GEOS_UNUSED_VAR( domain );
  };

  /** @} */


  /// @brief Runs the solver when the containing event is triggered.
  virtual bool execute( real64 const time_n,
                        real64 const dt,
                        integer const cycleNumber,
                        integer const eventCounter,
                        real64 const eventProgress,
                        DomainPartition & domain ) override;

  /// @brief Dispatches an explicit element kernel over selected mesh regions.
  template< typename ... PARAMS >
  real64 explicitKernelDispatch( MeshLevel & mesh,
                                 string_array const & targetRegions,
                                 string const & finiteElementName,
                                 real64 const dt,
                                 std::string const & elementListName );

  /// Repository keys for MPM-owned auxiliary groups.
  struct groupKeyStruct : PhysicsSolverBase::groupKeyStruct
  {
    static constexpr char const * cohesiveZoneManagerString() { return "CohesiveZoneManager"; };
    static constexpr char const * mpmEventManagerString() { return "MPMEvents"; };
  };

  /// Registered grid-field and solver-option names used by data repository views.
  struct viewKeyStruct : PhysicsSolverBase::viewKeyStruct
  {
    // Boundary conditions / loading
    static constexpr char const * boundaryNodesString() { return "boundaryNodes"; }
    static constexpr char const * gridSurfaceTensionForceString() { return "gridSurfaceTensionForce"; }

    // Contact
    static constexpr char const * gridContactForceString() { return "gridContactForce"; }
    static constexpr char const * gridInterfaceMechanismString() { return "gridInterfaceMechanism"; }
    static constexpr char const * gridInterfaceNormalForceString() { return "gridInterfaceNormalForce"; }
    static constexpr char const * gridInterfaceTangentialForceString() { return "gridInterfaceTangentialForce"; }
    static constexpr char const * gridInterfaceNormalVelocityJumpString() { return "gridInterfaceNormalVelocityJump"; }
    static constexpr char const * gridInterfaceTangentialVelocityJumpString() { return "gridInterfaceTangentialVelocityJump"; }
    static constexpr char const * gridInterfaceNormalDisplacementJumpString() { return "gridInterfaceNormalDisplacementJump"; }
    static constexpr char const * gridInterfaceTangentialDisplacementJumpString() { return "gridInterfaceTangentialDisplacementJump"; }

    // Cohesive zone
    static constexpr char const * gridCohesiveAreaString() { return "gridCohesiveArea"; }
    static constexpr char const * gridCohesiveFieldFlagString() { return "gridCohesiveFieldFlag"; }
    static constexpr char const * gridCohesiveForceString() { return "gridCohesiveForce"; }
    static constexpr char const * gridCohesiveNodeString() { return "gridCohesiveNode"; }
    static constexpr char const * gridReferenceAreaVectorString() { return "gridReferenceAreaVector"; }
    static constexpr char const * gridReferenceMaterialVolumeString() { return "gridReferenceMaterialVolume"; }
    static constexpr char const * gridReferenceSurfacePositionString() { return "gridReferenceSurfacePosition"; }

    // Damage / DFG / crack tips
    static constexpr char const * gridDamageGradientString() { return "gridDamageGradient"; }
    static constexpr char const * gridDamageString() { return "gridDamage"; }
    static constexpr char const * gridMappingNormalTensorString() { return "gridMappingNormalTensor"; }
    static constexpr char const * gridMaxDamageString() { return "gridMaxDamage"; }
    static constexpr char const * gridSingleFieldStateFractionString() { return "gridSingleFieldStateFraction"; }

    // Integration / update
    static constexpr char const * gridInternalForceString() { return "gridInternalForce"; }
    static constexpr char const * timeIntegrationOptionString() { return "timeIntegrationOption"; }

    // Mapping / shape functions
    static constexpr char const * gridMaxMappedParticleIDString() { return "gridMaxMappedParticleID"; }
    static constexpr char const * gridParticleMappedSurfaceNormalString() { return "gridParticleMappedSurfaceNormal"; }

    // Miscellaneous
    static constexpr char const * bufferNodesString() { return "bufferNodes"; }
    static constexpr char const * cflFactorString() { return "cflFactor"; }
    static constexpr char const * gridAccelerationString() { return "gridAcceleration"; }
    // Shared active field mask computed from post-P2G nodal mass.
    static constexpr char const * gridActiveString() { return "gridActive"; }
    static constexpr char const * gridBackgroundStressString() { return "gridBackgroundStress"; }
    static constexpr char const * gridBasedSurfaceNormalString() { return "gridBasedSurfaceNormal"; }
    static constexpr char const * gridBasedSurfacePositionString() { return "gridBasedSurfacePosition"; }
    static constexpr char const * gridHasBinderString() { return "gridHasBinder"; }
    static constexpr char const * gridCenterOfMassString() { return "gridCenterOfMass"; }
    static constexpr char const * gridCenterOfVolumeString() { return "gridCenterOfVolume"; }
    static constexpr char const * gridDisplacementString() { return "gridDisplacement"; }
    static constexpr char const * gridDVelocityString() { return "gridDVelocity"; }
    static constexpr char const * gridDVPlusString() { return "gridDVPlus"; }
    static constexpr char const * gridExplicitSurfaceNormalString() { return "gridExplicitSurfaceNormal"; }
    static constexpr char const * gridExternalForceString() { return "gridExternalForce"; }
    static constexpr char const * gridMassString() { return "gridMass"; }
    static constexpr char const * gridMaterialVolumeString() { return "gridMaterialVolume"; }
    static constexpr char const * gridSurfaceMaterialVolumeString() { return "gridSurfaceMaterialVolume"; }
    static constexpr char const * gridMomentumString() { return "gridMomentum"; }
    static constexpr char const * gridPrincipalExplicitSurfaceNormalString() { return "gridPrincipalExplicitSurfaceNormal"; }
    static constexpr char const * gridSurfaceAreaString() { return "gridSurfaceArea"; }
    static constexpr char const * gridSurfaceFieldMassString() { return "gridSurfaceFieldMass"; }
    static constexpr char const * gridSurfaceMassString() { return "gridSurfaceMass"; }
    static constexpr char const * gridSurfaceNormalString() { return "gridSurfaceNormal"; }
    static constexpr char const * gridSurfaceNormalWeightNormalizationString() { return "gridSurfaceNormalWeightNormalization"; }
    static constexpr char const * gridSurfaceNormalWeightsString() { return "gridSurfaceNormalWeights"; }
    static constexpr char const * gridSurfacePositionString() { return "gridSurfacePosition"; }
    static constexpr char const * gridUncontactedVelocityString() { return "gridUncontactedVelocity"; }
    static constexpr char const * gridVelocityString() { return "gridVelocity"; }
    static constexpr char const * gridVPlusString() { return "gridVPlus"; }
    static constexpr char const * solidMaterialNamesString() { return "solidMaterialNames"; }

    dataRepository::ViewKey timeIntegrationOption = { timeIntegrationOptionString() };
  } solidMechanicsViewKeys;

  /// Non-owning references to the grid, particle, and partition managers for one explicit step.
  struct ExplicitStepManagers
  {
    SpatialPartition & partition;
    arrayView1d< int const > periodic;
    ParticleManager & particleManager;
    MeshLevel & mesh;
    NodeManager & nodeManager;
  };

  /// @brief Initializes background-grid and particle state for the selected spatial partition.
  void initialize( NodeManager & nodeManager,
                   ParticleManager & particleManager,
                   SpatialPartition & partition );

  /// @brief Selects a particle's nodal velocity field from contact-group and damage-side information.
  static GEOS_FORCE_INLINE
  GEOS_HOST_DEVICE
  localIndex partitionField( int numContactGroups,
                             int damageFieldPartitioning,
                             localIndex particleGroup,
                             arraySlice1d< real64 const > const particleDamageGradient,
                             arraySlice1d< real64 const > const particleSurfaceNormal,
                             arraySlice1d< real64 const > const gridDamageGradient );

  /// @brief Measures directional coherence of two mapped surface-normal tensors for DFG contact.
  static GEOS_FORCE_INLINE
  GEOS_HOST_DEVICE
  real64 computeSurfaceQualityFromMappingNormalTensor( int const & planeStrain,
                                                       arraySlice1d< real64 const > const mappingNormalTensorA,
                                                       arraySlice1d< real64 const > const mappingNormalTensorB,
                                                       real64 const & materialVolumeA,
                                                       real64 const & materialVolumeB );

  /// @brief Applies MPM events scheduled for the current time step.
  void triggerEvents( const real64 dt,
                      const real64 time_n,
                      ParticleManager & particleManager,
                      SpatialPartition & partition );

  /// @brief Initializes particle thermodynamic state from a prescribed event.
  void setInitialTemperatureAndPressure( ParticleManager & particleManager,
                                         SetInitialTemperatureAndPressureMPMEvent const & event );

  /// @brief Updates completion state for events following time advancement.
  void checkEventCompletion( const real64 time_n );

  /// @brief Transfers a selected particle population to another constitutive material.
  void performMaterialSwap( ParticleManager & particleManager,
                            string sourceRegionName,
                            string destinationRegionName );

  /// @brief Updates structured-grid extents and spacing for prescribed domain deformation.
  void resizeGrid( SpatialPartition & partition,
                   NodeManager & nodeManager,
                   real64 const dt );

  /// @brief Reduces selected nodal fields across MPI-shared nodes with the supplied operator.
  void syncGridFields( stdVector< std::string > const & fieldNames,
                       DomainPartition & domain,
                       NodeManager & nodeManager,
                       MeshLevel & mesh,
                       MPI_Op op,
                       bool syncGridOnDevice = false );

  /// @brief Copies finalized nodal field values from owners to ghost replicas.
  void replaceGridFieldsOwnerToGhost( stdVector< std::string > const & fieldNames,
                                      DomainPartition & domain,
                                      NodeManager & nodeManager,
                                      MeshLevel & mesh );

  /// @brief Synchronizes constrained nodal kinematics before the final grid-to-particle gather.
  void synchronizePostBoundaryKinematicFieldsForG2P( DomainPartition & domain,
                                                     NodeManager & nodeManager,
                                                     MeshLevel & mesh );

  /// @brief Builds a shared active-field mask from synchronized nodal mass.
  void computeActiveGridFieldsForExplicitStep( DomainPartition & domain,
                                               NodeManager & nodeManager,
                                               MeshLevel & mesh );

  /// @brief Enforces symmetry constraints on a nodal vector field.
  void enforceGridVectorFieldSymmetryBC( arrayView3d< real64 > const & vectorMultiField,
                                         arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const gridPosition,
                                         Group & nodeSets );

  /// @brief Applies prescribed nodal velocity and displacement boundary conditions.
  void applyEssentialBCs( const real64 dt,
                          const real64 time_n,
                          NodeManager & nodeManager );

  /// @brief Collects non-owning references to the managers used during explicit stepping.
  ExplicitStepManagers getExplicitStepManagers( DomainPartition & domain );

  /// @brief Checks whether any coordinate direction uses a periodic boundary.
  bool hasPeriodicBoundary( arrayView1d< int const > const periodic ) const;

  /// @brief Initializes step-local storage and clears grid accumulators.
  void initializeAndResetExplicitStepState( int const cycleNumber,
                                            NodeManager & nodeManager,
                                            ParticleManager & particleManager,
                                            SpatialPartition & partition );

  /// @brief Prepares particle populations and geometry for the subsequent interpolation pass.
  void prepareParticleTopologyForExplicitStep( real64 const dt,
                                               real64 const time_n,
                                               int const cycleNumber,
                                               DomainPartition & domain,
                                               ParticleManager & particleManager,
                                               SpatialPartition & partition,
                                               arrayView1d< int const > const periodic );

  /// @brief Refreshes neighbor, damage, and contact state before particle-to-grid mapping.
  void buildNeighborhoodAndContactStateForExplicitStep( real64 const dt,
                                                        real64 const time_n,
                                                        int const cycleNumber,
                                                        ParticleManager & particleManager );

  /// @brief Builds valid particle-to-node interpolation arrays for active particles.
  void populateParticleGridMappingForExplicitStep( ParticleManager & particleManager,
                                                  NodeManager & nodeManager );

  /// @brief Updates damage-gradient and surface fields used by contact and cohesive laws.
  void updateDamageAndSurfaceFieldsForExplicitStep( DomainPartition & domain,
                                                    ParticleManager & particleManager,
                                                    NodeManager & nodeManager,
                                                    MeshLevel & mesh );

  /// @brief Assembles particle loads and optional background-stress data.
  void computeParticleLoadsAndBackgroundFieldsForExplicitStep( real64 const time_n,
                                                               ParticleManager & particleManager,
                                                               NodeManager & nodeManager );

  /// @brief Advances trial nodal dynamics and resolves material contact.
  void updateGridDynamicsAndContactForExplicitStep( real64 const dt,
                                                    ParticleManager & particleManager,
                                                    NodeManager & nodeManager );

  /// @brief Constrains the uncontacted FMPM velocity seed without adding material-contact impulse.
  void applyFMPMUncontactedVelocityBoundaryConditions( NodeManager & nodeManager );

  /// @brief Imposes homogeneous constraints on higher-order FMPM velocity increments.
  void applyFMPMVelocityIncrementBoundaryConditions( NodeManager & nodeManager,
                                                    arrayView3d< real64 > const velocityIncrement );

  /// @brief Copies prescribed components between FMPM velocity fields and rebuilds buffer values.
  void copyConstrainedFMPMBoundaryVelocity( NodeManager & nodeManager,
                                            arrayView3d< real64 > const destinationVelocity,
                                            arrayView3d< real64 const > const sourceVelocity );

  /// @brief Applies prescribed domain deformation and grid-boundary motion.
  void applyPrescribedDeformationAndBoundaryConditionsForExplicitStep( real64 const dt,
                                                                       real64 const time_n,
                                                                       int const cycleNumber,
                                                                       ParticleManager & particleManager,
                                                                       NodeManager & nodeManager,
                                                                       SpatialPartition & partition );

  /// @brief Gathers grid motion and updates particle positions and deformation state.
  void updateParticleKinematicsForExplicitStep( real64 const dt,
                                                real64 const time_n,
                                                ParticleManager & particleManager,
                                                SpatialPartition & partition );

  /// @brief Updates constitutive response and optional thermodynamic fields.
  void updateConstitutiveAndThermalStateForExplicitStep( real64 const dt,
                                                         ParticleManager & particleManager );

  /// @brief Writes scheduled diagnostics and computes the next stability-limited time step.
  real64 writeOutputsAndComputeStableTimeStepForExplicitStep( real64 const time_n,
                                                              real64 const dt,
                                                              int const cycleNumber,
                                                              ParticleManager & particleManager );

  /// @brief Updates the grid and removes particles marked invalid or out of bounds.
  void resizeGridAndCleanParticlesForExplicitStep( real64 const dt,
                                                   DomainPartition & domain,
                                                   ParticleManager & particleManager,
                                                   NodeManager & nodeManager,
                                                   SpatialPartition & partition,
                                                   arrayView1d< int const > const periodic );

  /// @brief Records an explicitly labeled stage in the solver performance log.
  void logAndProfile( std::string const & label );

  /// @brief Records a labeled performance stage with particle/grid diagnostics.
  void logAndProfile( std::string const & label,
                      ParticleManager & particleManager,
                      NodeManager & nodeManager );

  /// @brief Logs particle/grid momentum diagnostics at a named step stage.
  void logMomentumSum( std::string const & label,
                       ParticleManager & particleManager,
                       NodeManager & nodeManager );

  /// @brief Logs momentum with an explicitly supplied diagnostic nodal velocity field.
  void logMomentumSum( std::string const & label,
                       ParticleManager & particleManager,
                       NodeManager & nodeManager,
                       arrayView3d< real64 const > const diagnosticVelocity,
                       std::string const & diagnosticVelocityName );

  /// @brief Makes particle wrappers available for host-side diagnostics and output.
  void moveParticleWrappersToHost( ParticleManager & particleManager );

  /// @brief Updates ghost particle data required by local mapping and neighborhoods.
  void performExplicitStepParticleGhosting( DomainPartition & domain,
                                            ParticleManager & particleManager,
                                            SpatialPartition & partition );

  /// @brief Builds the active-particle ordinal map used by compact mapping arrays.
  void setActiveParticleIndices( ParticleManager & particleManager );

  /// @brief Incorporates boundary-condition settings into the contact flags.
  void updateContactFlagsFromBoundaryConditions();

  /// @brief Disables surface reconstruction data on sufficiently damaged particles.
  void disableSurfaceDataOnDamagedParticles( ParticleManager & particleManager );

  /// @brief Disables surface reconstruction data on melted particles.
  void disableSurfaceDataOnMeltedParticles( ParticleManager & particleManager );

  /// @brief Determines whether the current options require a nodal neighbor list.
  void updateNodalNeighborListRequirement();

  /// @brief Clears particle cohesive-force accumulators before traction assembly.
  void zeroParticleCohesiveForces( ParticleManager & particleManager );

  /// @brief Resets surface tags used during cohesive-zone initialization and updates.
  void resetCohesiveSurfaceFlags( ParticleManager & particleManager );

  /// @brief Builds reference cohesive data and evaluates cohesive tractions for the step.
  void updateCohesiveZonesForExplicitStep( real64 const dt,
                                           DomainPartition & domain,
                                           ParticleManager & particleManager,
                                           NodeManager & nodeManager,
                                           MeshLevel & mesh );

  /// @brief Checks interpolation consistency and particle mapping diagnostics.
  void shapeFunctionDiagnostics( ParticleManager & particleManager );

  /// @brief Scatters particle momentum, mass, stresses, and related fields to the grid.
  void performParticleToGridForExplicitStep( real64 const time_n,
                                             integer const cycleNumber,
                                             ParticleManager & particleManager,
                                             NodeManager & nodeManager );

  /// @brief Synchronizes primary nodal accumulators after particle-to-grid mapping.
  void syncGridFieldsForExplicitStep( DomainPartition & domain,
                                      NodeManager & nodeManager,
                                      MeshLevel & mesh );

  /// @brief Synchronizes nodal capillary-force contributions across partitions.
  void syncSurfaceTensionForcesForExplicitStep( DomainPartition & domain,
                                                NodeManager & nodeManager,
                                                MeshLevel & mesh );

  /// @brief Assembles diffuse-interface capillary forces and synchronizes them.
  void computeAndSyncSurfaceTensionForcesForExplicitStep( ParticleManager & particleManager,
                                                          DomainPartition & domain,
                                                          NodeManager & nodeManager,
                                                          MeshLevel & mesh );

  /// @brief Applies nodal symmetry and normalizes mapped grid fields.
  void enforceGridFieldSymmetryAndNormalize( NodeManager & nodeManager );

  /// @brief Updates prescribed particle and domain kinematics for the current time.
  void updatePrescribedKinematicsForExplicitStep( real64 const dt,
                                                  real64 const time_n,
                                                  ParticleManager & particleManager,
                                                  SpatialPartition & partition );

  /// @brief Migrates particles to the spatial partitions owning their updated positions.
  void repartitionParticlesForExplicitStep( DomainPartition & domain,
                                            ParticleManager & particleManager,
                                            SpatialPartition & partition,
                                            arrayView1d< int const > const periodic );

  /// @brief Applies a prescribed background velocity gradient to particle motion.
  void applySuperimposedVelocityGradient( const real64 dt,
                                          ParticleManager & particleManager,
                                          SpatialPartition & partition );

  /// @brief Estimates nodal surface area by direct sampling of a local surface patch.
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  real64 bruteForceNodalAreaIntegration( real64 const (&hEl)[3],
                                         integer const & numSurfaceIntegrationPoints,
                                         real64 const & L,
                                         real64 const & dA,
                                         arraySlice1d< real64 const > const surfaceNormal,
                                         arraySlice1d< real64 const > const surfacePosition );

  /// @brief Integrates nodal surface area over neighboring particle-domain geometry.
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  real64 meshNodalAreaIntegration( real64 const (&hEl)[3],
                                  integer const & numSurfaceIntegrationPoints,
                                  localIndex const & numNeighbors,
                                  arraySlice1d< localIndex const > const neighborRegions,
                                  arraySlice1d< localIndex const > const neighborSubRegions,
                                  arraySlice1d< localIndex const > const neighborIndices,
                                  ParticleManager::ParticleViewConst< arrayView2d< real64 const > > const particleCenter,
                                  ParticleManager::ParticleViewConst< arrayView3d< real64 const > > const particleGeometryCurrent,
                                  arraySlice1d< real64 const > const surfaceNormal,
                                  arraySlice1d< real64 const > const surfacePosition );

  /// @brief Estimates nodal surface area from a neighboring-particle convex hull.
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  real64 convexHullAreaIntegration( real64 const (& hEl)[3],
                                    real64 const (& nodePosition)[3],
                                    integer const & maxNodalNeighbors,
                                    integer const & numSurfaceIntegrationPoints,
                                    localIndex const & numNeighbors,
                                    arraySlice1d< localIndex const > const neighborRegions,
                                    arraySlice1d< localIndex const > const neighborSubRegions,
                                    arraySlice1d< localIndex const > const neighborIndices,
                                    ParticleManager::ParticleViewConst< arrayView2d< real64 const > > const particleCenter,
                                    ParticleManager::ParticleViewConst< arrayView3d< real64 const > > const particleGeometryCurrent,
                                    arraySlice1d< real64 const > const surfaceNormal,
                                    arraySlice1d< real64 const > const surfacePosition );

  /// @brief Evaluates nodal surface areas using the selected integration method.
  void computeNodalAreas( NodeManager & nodeManager );

  /// @brief Normalizes mapped surface normals and surface-location accumulators.
  void normalizeGridSurfaceNormalsAndPositions( NodeManager & nodeManager );

  /// @brief Forms nodal normal weights from nearby particle surface information.
  void computeGridSurfaceNormalWeights( ParticleManager & particleManager,
                                        NodeManager & nodeManager );

  /// @brief Assembles pairwise nodal contact forces, including normal and frictional responses.
  void computeContactForces( real64 const dt,
                             ParticleManager & particleManager,
                             NodeManager & nodeManager );

  /// @brief Records the selected contact or cohesive interface mechanism and jump measures.
  void updateInterfaceMechanismDiagnostics( NodeManager & nodeManager );

  /// @brief Forms the cumulative contact-impulse target for an uncorrected FMPM trial velocity.
  /// The caller converts the change in the cumulative target into an increment.
  void computeFMPMNetContactMomentumTarget( real64 const dt,
                                            ParticleManager & particleManager,
                                            NodeManager & nodeManager,
                                            arrayView3d< real64 const > const vUncorrectedTotal,
                                            arrayView3d< real64 > const contactMomentumTarget );

  /// @brief Initializes global and contact-group-pair friction coefficients.
  void initializeFrictionCoefficients();

  /// @brief Evaluates the impulse form of the two-field nodal contact law for FMPM.
  static GEOS_FORCE_INLINE
  GEOS_HOST_DEVICE
  void computePairwiseNodalContactImpulse( mpm::ContactGapCorrectionOption const & contactGapCorrection,
                                           mpm::OverlapCorrectionOption const & overlapCorrection,
                                           real64 const overlapThreshold1,
                                           real64 const overlapThreshold2,
                                           real64 const maxParticleVelocitySquared,
                                           real64 const (&hEl) [3],
                                           real64 const cellVolume,
                                           int const & planeStrain,
                                           real64 const & smallMass,
                                           int const & useSurfacePositionForContact,
                                           int const & useCohesiveTangentialForces,
                                           bool & separable,
                                           real64 const & dt,
                                           real64 const & frictionCoefficient,
                                           real64 ( &nAB )[3],
                                           real64 const & mA,
                                           real64 const & mB,
                                           real64 const & VA,
                                           real64 const & VB,
                                           arraySlice1d< real64 const > const vA,
                                           arraySlice1d< real64 const > const GEOS_UNUSED_PARAM( vB ),
                                           real64 const & qA0,
                                           real64 const & qA1,
                                           real64 const & qA2,
                                           real64 const & qB0,
                                           real64 const & qB1,
                                           real64 const & qB2,
                                           real64 const spmA,
                                           real64 const spmB,
                                           arraySlice1d< real64 const > const sA,
                                           arraySlice1d< real64 const > const sB,
                                           arraySlice1d< real64 const > const xA,
                                           arraySlice1d< real64 const > const xB,
                                           arraySlice1d< real64 > const jA,
                                           arraySlice1d< real64 > const jB );

  /// @brief Evaluates normal and tangential contact forces between two nodal fields.
  static GEOS_FORCE_INLINE
  GEOS_HOST_DEVICE
  void computePairwiseNodalContactForce( mpm::ContactGapCorrectionOption const & contactGapCorrection,
                                         mpm::OverlapCorrectionOption const & overlapCorrection,
                                         real64 const overlapThreshold1,
                                         real64 const overlapThreshold2,
                                         real64 const maxParticleVelocitySquared,
                                         real64 const (&hEl) [3],
                                         real64 const cellVolume,
                                         int const & planeStrain,
                                         real64 const & smallMass,
                                         int const & useSurfacePositionForContact,
                                         int const & useCohesiveTangentialForces,
                                         bool & separable,
                                         real64 const & dt,
                                         real64 const & frictionCoefficient,
                                         real64 ( &nAB )[3],
                                         real64 const & mA,
                                         real64 const & mB,
                                         real64 const & VA,
                                         real64 const & VB,
                                         arraySlice1d< real64 const > const vA,
                                         arraySlice1d< real64 const > const GEOS_UNUSED_PARAM( vB ),
                                         arraySlice1d< real64 const > const qA,
                                         arraySlice1d< real64 const > const qB,
                                         real64 const spmA,
                                         real64 const spmB,
                                         arraySlice1d< real64 const > const sA,
                                         arraySlice1d< real64 const > const sB,
                                         arraySlice1d< real64 const > const xA,
                                         arraySlice1d< real64 const > const xB,
                                         arraySlice1d< real64 const > const centerOfVolumeA,
                                         arraySlice1d< real64 const > const centerOfVolumeB,
                                         arraySlice1d< real64 > const fA,
                                         arraySlice1d< real64 > const fB );

  /// @brief Builds two tangential unit vectors orthogonal to a supplied normal.
  static GEOS_FORCE_INLINE
  GEOS_HOST_DEVICE
  void computeOrthonormalBasis( const real64 * e1,
                                real64 * e2,
                                real64 * e3 );

  /// @brief Assigns readable labels to nodal multifield data for output.
  void setGridFieldLabels( NodeManager & nodeManager );

  /// @brief Builds particle-neighbor lists using the active neighborhood radius.
  real64 computeNeighborList( ParticleManager & particleManager );

  /// @brief Builds nodal particle-neighborhood relations needed by surface/contact algorithms.
  void generateNodalNeighborList( ParticleManager & particleManager,
                                  NodeManager & nodeManager );

  /// @brief Organizes particles into spatial bins for efficient neighborhood queries.
  void optimizeBinSort( ParticleManager & particleManager );

  /// @brief Evaluates the radial smoothing kernel at a particle-query distance.
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  real64 kernel( real64 const & r );

  /// @brief Computes the inverse radial-kernel relation used by surface estimation.
  GEOS_HOST
  GEOS_FORCE_INLINE
  real64 inverseKernel( real64 const & d );

  /// @brief Evaluates the dimension-aware radial kernel on host or device.
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  real64 kernelDevice( real64 const & r,
                       real64 const & neighborRadius,
                       int const & planeStrain );

  /// @brief Evaluates the spatial gradient of the radial smoothing kernel.
  /// @param x Query point; xp source particle position; r their separation.
  void kernelGradient( arraySlice1d< real64 const > const & x,
                       arraySlice1d< real64 const > const & xp,
                       real64 const & r,
                       real64 ( &result )[3] );

  /// @brief Evaluates the dimension-aware smoothing-kernel gradient on host or device.
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  void kernelGradientDevice( arraySlice1d< real64 const > const & x,
                             real64 const (&xp)[3],
                             real64 const & r,
                             real64 const & neighborRadius,
                             int const & planeStrain,
                             real64 ( &result )[3] );

  /// @brief Interpolates a scalar field from volume-weighted neighboring particles.
  /// Neighbor inputs xp, Vp, and fp contain positions, volumes, and scalar samples.
  GEOS_HOST
  GEOS_FORCE_INLINE
  real64 computeKernelField( arraySlice1d< real64 const > const x,
                             localIndex const numNeighbors,
                             arrayView2d< real64 const > const & xp,
                             arrayView1d< real64 const > const & Vp,
                             arrayView1d< real64 const > const & fp );

  /// @brief Interpolates a particle scalar field through region/subregion neighbor indices.
  template< typename DAMAGE >
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  real64 computeKernelFieldDevice( arraySlice1d< real64 const > const x,
                                   real64 const & neighborRadius,
                                   int const & planeStrain,
                                   localIndex const & numNeighbors,
                                   arraySlice1d< localIndex const > const regionIndices,
                                   arraySlice1d< localIndex const > const subRegionIndices,
                                   arraySlice1d< localIndex const > const particleIndices,
                                   ParticleManager::ParticleView< arrayView1d< real64 const > > particleVolumeView,
                                   ParticleManager::ParticleView< arrayView2d< real64 const > > particlePositionView,
                                   DAMAGE const & damage );

  /// @brief Evaluates the spatial gradient of a smoothed particle scalar field.
  /// Neighbor inputs xp, Vp, and fp contain positions, volumes, and scalar samples.
  GEOS_HOST
  GEOS_FORCE_INLINE
  void computeKernelFieldGradient( arraySlice1d< real64 const > const x,
                                   localIndex const numNeighbors,
                                   arrayView2d< real64 const > const & xp,
                                   arrayView1d< real64 const > const & Vp,
                                   arrayView1d< real64 const > const & fp,
                                   arraySlice1d< real64 > const result );

  /// @brief Evaluates a particle scalar-field gradient using host/device neighbor views.
  template< typename DAMAGE >
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  void computeKernelFieldGradientDevice( arraySlice1d< real64 const > const x,
                                         real64 const & neighborRadius,
                                         int const & planeStrain,
                                         localIndex const & numNeighbors,
                                         arraySlice1d< localIndex const > const regionIndices,
                                         arraySlice1d< localIndex const > const subRegionIndices,
                                         arraySlice1d< localIndex const > const particleIndices,
                                         ParticleManager::ParticleView< arrayView1d< real64 const > > particleVolumeView,
                                         ParticleManager::ParticleView< arrayView2d< real64 const > > particlePositionView,
                                         DAMAGE const & damage,
                                         real64 ( &result )[3] );

  /// @brief Evaluates the divergence of a smoothed particle vector field.
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  void computeKernelFieldDivergenceDevice( arraySlice1d< real64 const > const x,
                                           real64 const & neighborRadius,
                                           int const & planeStrain,
                                           localIndex const & numNeighbors,
                                           arraySlice1d< localIndex const > const regionIndices,
                                           arraySlice1d< localIndex const > const subRegionIndices,
                                           arraySlice1d< localIndex const > const particleIndices,
                                           ParticleManager::ParticleView< arrayView1d< real64 const > > particleVolumeView,
                                           ParticleManager::ParticleView< arrayView2d< real64 const > > particlePositionView,
                                           ParticleManager::ParticleView< arrayView2d< real64 const > > particleVectorView,
                                           real64 & result );

  /// @brief Evaluates the gradient tensor of a smoothed particle vector field.
  /// Neighbor inputs xp, Vp, and fp contain positions, volumes, and vector samples.
  GEOS_HOST
  GEOS_FORCE_INLINE
  void computeKernelVectorGradient( arraySlice1d< real64 const > const x,
                                    localIndex const numNeighbors,
                                    arrayView2d< real64 const > const & xp,
                                    arrayView1d< real64 const > const & Vp,
                                    arrayView2d< real64 const > const & fp,
                                    arraySlice2d< real64 > const result );

  /// @brief Evaluates a particle vector-field gradient using host/device neighbor views.
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  void computeKernelVectorGradientDevice( arraySlice1d< real64 const > const x,
                                          real64 const & neighborRadius,
                                          int const & planeStrain,
                                          localIndex const numNeighbors,
                                          arraySlice1d< localIndex const > const regionIndices,
                                          arraySlice1d< localIndex const > const subRegionIndices,
                                          arraySlice1d< localIndex const > const particleIndices,
                                          ParticleManager::ParticleView< arrayView1d< real64 const > > particleVolumeView,
                                          ParticleManager::ParticleView< arrayView2d< real64 const > > particlePositionView,
                                          ParticleManager::ParticleView< arrayView2d< real64 const > > particleDisplacementView,
                                          arraySlice2d< real64 > const result );

  /// @brief Tests whether a surface classification marks a damage-generated surface.
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  bool markSurfaceAsDamage( integer const & surfaceFlag );

  /// @brief Reconstructs local particle damage gradients from neighboring particles.
  void computeDamageFieldGradient( ParticleManager & particleManager );

  /// @brief Estimates distance from particles to a detected crack tip.
  void computeDistanceToCrackTip( ParticleManager & particleManager );

  /// @brief Estimates particle surface points from local surface information.
  void estimateParticleSurfacePosition( ParticleManager & particleManager );

  /// @brief Updates overloaded damage/interior flags while preserving explicit surface tags.
  void updateSurfaceFlagOverload( ParticleManager & particleManager );

  /// @brief Scatters particle damage gradients to nodal fields for DFG partitioning.
  void projectDamageFieldGradientToGrid( DomainPartition & domain,
                                         ParticleManager & particleManager,
                                         NodeManager & nodeManager,
                                         MeshLevel & mesh );

  /// @brief Assigns effective nodal fields to each mapped particle contribution.
  void computeParticleFieldMappings( ParticleManager & particleManager,
                                     NodeManager & nodeManager );

  /// @brief Advances particle deformation gradients using the reconstructed velocity gradient.
  void updateDeformationGradient( real64 dt,
                                  ParticleManager & particleManager );

  /// @brief Computes prescribed domain deformation rates from the stress-control options.
  void stressControl( const real64 dt,
                      ParticleManager & particleManager,
                      SpatialPartition & partition );

  /// @brief Runs the adaptive stress-seeking feedback controller and its state transitions.
  void adaptiveSeekStressControl( const real64 dt,
                                  real64 const (& targetStress)[3],
                                  arrayView1d< real64 const > const currentStress,
                                  real64 const boxMaterialVolume,
                                  real64 const maximumBulkModulus );

  /// @brief Allocates and initializes particle mechanical and geometry fields.
  void initializeParticleFields( ParticleManager & particleManager );

  /// @brief Copies optional particle strength-scale data into constitutive-model fields.
  void initializeConstitutiveModelDependencies( ParticleManager & particleManager );

  /// @brief Transfers current particle kinematics and state into constitutive wrappers.
  void updateConstitutiveModelDependencies( ParticleManager & particleManager );

  /// @brief Evaluates constitutive stress updates for active particle materials.
  void updateStress( real64 dt,
                     ParticleManager & particleManager );

  /// @brief Updates particle kinematics and geometry from grid-to-particle mapping.
  void particleKinematicUpdate( const real64 dt,
                                ParticleManager & particleManager );

  /// @brief Reduces and writes selected box-averaged particle quantities.
  void computeAndWriteBoxAverage( const real64 dt,
                                  const real64 time_n,
                                  ParticleManager & particleManager );

  /// @brief Writes particle diagnostics and optional history output.
  void writeParticleData( const real64 time_n,
                          ParticleManager & particleManager );

  /// @brief Reduces total material volume and stress measures over the configured box.
  void computeBoxMetrics( ParticleManager & particleManager,
                          arrayView1d< real64 > boxStress,
                          real64 & boxMaterialVolume );

  /// @brief Clears step-local nodal field arrays while retaining persistent grid state.
  void initializeGridFields( NodeManager & nodeManager );

  /// @brief Updates time-dependent prescribed boundary conditions.
  void boundaryConditionUpdate( real64 dt, real64 time_n );

  /// @brief Maps particle surface normals to grid nodes for contact and geometry recovery.
  void projectParticleSurfaceNormalsToGrid( DomainPartition & domain,
                                            ParticleManager & particleManager,
                                            NodeManager & nodeManager,
                                            MeshLevel & mesh );

  /// @brief Constructs reference nodal and particle data for cohesive interfaces.
  void initializeCohesiveReferenceConfiguration( DomainPartition & domain,
                                                 ParticleManager & particleManager,
                                                 NodeManager & nodeManager,
                                                 MeshLevel & mesh );


  /// @brief Removes a vector's component along a specified plane normal.
  void projectToPlane( real64 const (&vector)[3],
                       real64 const (&normal)[3],
                       real64 ( &projection )[3] );

  /// @brief Evaluates cohesive tractions and maps the resulting forces to particles.
  void enforceCohesiveLaw( real64 dt,
                           ParticleManager & particleManager,
                           NodeManager & nodeManager,
                           SpatialPartition & partition );

  /// @brief Constructs optional background-stress data on the grid.
  GEOS_HOST_DEVICE
  void updateGridBackgroundStress( NodeManager & nodeManager );

  /// @brief Assembles nodal mass, momentum, forces, and surface fields from effective particle stencils.
  void particleToGrid( real64 const time_n,
                       integer const cycleNumber,
                       ParticleManager & particleManager,
                       NodeManager & nodeManager );

  /// @brief Assembles diffuse-interface capillary force from mapped phase-color gradients.
  void computePartitionedSurfaceTensionForces( ParticleManager & particleManager,
                                               NodeManager & nodeManager );

  /// @brief Computes trial nodal velocities from mass, momentum, and forces.
  void gridTrialUpdate( real64 dt,
                        NodeManager & nodeManager );

  /// @brief Applies the configured multifield nodal contact response.
  void enforceContact( real64 dt,
                       ParticleManager & particleManager,
                       NodeManager & nodeManager );

  /// @brief Constructs particle and grid surface normals and positions for contact.
  void computeParticleSurfaceNormalsAndPositions( ParticleManager & particleManager,
                                                  NodeManager & nodeManager,
                                                  DomainPartition & domain,
                                                  MeshLevel & mesh );

  /// @brief Projects ring-weighted meridional surface geometry with Q1 interpolation.
  void projectAxisymmetricSurfaceGeometry( ParticleManager & particleManager,
                                           NodeManager & nodeManager,
                                           DomainPartition & domain,
                                           MeshLevel & mesh,
                                           bool const computeNormals );

  /// @brief Reconstructs axisymmetric nodal surfaces from projected geometric fields.
  void reconstructAxisymmetricGridSurfaces( ParticleManager & particleManager,
                                            NodeManager & nodeManager );

  /// @brief Returns the axisymmetric projection support for a grid location.
  axisymmetricSurface::Support axisymmetricSurfaceSupport(
    arraySlice1d< real64 const, nodes::REFERENCE_POSITION_USD - 1 > const position ) const;

  /// @brief Maps a grid node to the corresponding axisymmetric color-field node.
  localIndex axisymmetricSurfaceColorNode(
    localIndex const node,
    arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const position ) const;

  /// @brief Fits pairwise interface normals and positions from particle neighborhoods.
  void computePairwiseLogisticRegressionSurfaceNormalsAndPositions( ParticleManager & particleManager,
                                                                    NodeManager & nodeManager );

  /// @brief Estimates nodal interface positions from mapped material volumes.
  void computeGridSurfacePositionFromVolume( NodeManager & nodeManager );

  /// @brief Integrates the nodal support volume on one side of a planar surface.
  static GEOS_HOST_DEVICE GEOS_FORCE_INLINE
  real64 computeVolumeFromSurfacePosition( real64 const (&n)[3],
                                           real64 const offset,
                                           real64 const (&hEl) [3],
                                           int const numDims,
                                           real64 const cellVolume );

  /// @brief Returns the maximum plane offset supported by a cell for a given normal.
  static GEOS_HOST_DEVICE GEOS_FORCE_INLINE
  real64 computeMaximumSurfacePositionOffset( real64 ( &n )[3],
                                              real64 const (&hEl)[3] );

  static GEOS_HOST_DEVICE

  /// @brief Fits a separating interface to particles assigned to two nodal fields.
  GEOS_FORCE_INLINE
  void logisticRegression( int const & planeStrain,
                           integer const & numContactGroups,
                           integer const & damageFieldPartitioning,
                           integer const & maxLRIterations,
                           real64 const & LRtolerance,
                           real64 const (&hEl)[3],
                           localIndex const & fieldA,
                           localIndex const & fieldB,
                           localIndex const & numNeighboringParticles,
                           arraySlice1d< localIndex const > const neighborRegions,
                           arraySlice1d< localIndex const > const neighborSubRegions,
                           arraySlice1d< localIndex const > const neighborIndices,
                           ParticleManager::ParticleViewConst< arrayView1d< localIndex const > > particleGroup,
                           ParticleManager::ParticleViewConst< arrayView2d< real64 const > > particleDamageGradient,
                           ParticleManager::ParticleViewConst< arrayView2d< real64 const > > particleSurfaceNormal,
                           ParticleManager::ParticleViewConst< arrayView2d< real64 const > > particlePosition,
                           arraySlice1d< real64 const, nodes::REFERENCE_POSITION_USD - 1 > const gridPosition,
                           arraySlice1d< real64 const > const gridDamageGradient,
                           real64 ( &normal0 )[3],
                           real64 ( &normal )[3],
                           real64 ( &surfacePosition )[3] );

  /// @brief Gathers reconstructed nodal surface geometry back to particles.
  void mapSurfaceNormalsAndPositionsToParticles( ParticleManager & particleManager,
                                                 NodeManager & nodeManager );

  /// @brief Interpolates a time-history table and its rate using the selected interpolation method.
  void interpolateTable( real64 x,
                         real64 dx,
                         arrayView2d< real64 const > const & table,
                         arrayView1d< real64 > output,
                         arrayView1d< real64 > outputRate,
                         mpm::InterpolationOption const & interpolationType );

  /// @brief Interpolates one scalar between given bounds with the configured rule.
  void interpolateValueInRange( real64 const & x,
                                real64 const & xmin,
                                real64 const & xmax,
                                real64 const & ymin,
                                real64 const & ymax,
                                real64 & output,
                                mpm::InterpolationOption interpolationType );

  /// @brief Evaluates the prescribed deformation-gradient history.
  void interpolateFTable( real64 dt, real64 time_n );

  /// @brief Evaluates the prescribed stress history.
  void interpolateStressTable( real64 dt, real64 time_n );

  /// @brief Dispatches the selected grid-to-particle update scheme.
  void gridToParticle( real64 dt,
                       ParticleManager & particleManager,
                       NodeManager & nodeManager,
                       DomainPartition & domain,
                       MeshLevel & mesh );

  /// @brief Advances particles using nodal velocity increments and a grid-based position update.
  void performFLIPUpdate( real64 dt,
                          ParticleManager & particleManager,
                          NodeManager & nodeManager );

  /// @brief Advances particle velocity and position directly from nodal values.
  void performPICUpdate( real64 dt,
                         ParticleManager & particleManager,
                         NodeManager & nodeManager );

  /// @brief Performs the selected-order XPIC recursive velocity projection and particle update.
  void performXPICUpdate( real64 dt,
                          ParticleManager & particleManager,
                          NodeManager & nodeManager,
                          DomainPartition & domain,
                          MeshLevel & mesh );

  /// @brief Advances particles with recursive FMPM velocity increments and optional Net contact.
  void performFMPMUpdate( real64 dt,
                          ParticleManager & particleManager,
                          NodeManager & nodeManager,
                          DomainPartition & domain,
                          MeshLevel & mesh );

  /// @brief Applies the change in cumulative contact impulse to an FMPM velocity increment.
  /// Applies (target - prior cumulative impulse) / nodal mass.
  void applyFMPMNetContactCorrection( real64 const dt,
                                      ParticleManager & particleManager,
                                      NodeManager & nodeManager,
                                      arrayView3d< real64 const > const vUncorrectedTotal,
                                      arrayView3d< real64 > const contactMomentumNet,
                                      arrayView3d< real64 > const contactMomentumTarget,
                                      arrayView3d< real64 > const velocityIncrement );

  /// @brief Transfers constitutive response and damage fields back to particle wrappers.
  void updateSolverDependencies( ParticleManager & particleManager );

  /// @brief Computes the stability-limited time increment from the active physical criteria.
  real64 getStableTimeStep( ParticleManager & particleManager );

  /// @brief Removes particles previously flagged for deletion.
  void deleteBadParticles( ParticleManager & particleManager );

  /// @brief Checks particle mapping arrays for invalid grid nodes and weights.
  void validateParticleMappings( ParticleManager & particleManager, NodeManager & nodeManager );

  /// @brief Writes particle mapping data for debugging.
  void dumpMappingDataToFile( ParticleManager & particleManager,
                              int const cycleNumbder,
                              std::string label );

  /// @brief Writes selected grid fields for debugging.
  void dumpGridDataToFile( NodeManager & nodeManager,
                           int const cycleNumber,
                           std::string label );

  /// @brief Writes constitutive-model data for debugging.
  void dumpModelDataToFile( ParticleManager & particleManager,
                            int const cycleNumbder,
                            std::string label );

  /// @brief Writes particle data for debugging.
  void dumpParticleDataToFile( ParticleManager & particleManager,
                               int const cycleNumber,
                               std::string label );

  /// @brief Reports accumulated solver profiling measurements.
  void printProfilingResults();

  /// @brief Classifies particle surfaces based on current geometry and damage.
  void computeSurfaceFlags( ParticleManager & particleManager );

  /// @brief Computes smoothed deformation-gradient information from particle neighbors.
  void computeSphF( ParticleManager & particleManager );

  /// @brief Decides whether two nodal fields may separate under the configured DFG criteria.
  static GEOS_FORCE_INLINE
  GEOS_HOST_DEVICE
  bool evaluateSeparabilityCriterion( int const & numContactGroups,
                                     real64 const & maxSingleFieldStateFractionForSeparability,
                                     real64 const & separabilityMinDamage,
                                     real64 const & thinFeatureDFGThreshold,
                                     real64 const & neighborRadius,
                                     real64 const & surfaceQualityThreshold,
                                     localIndex const & A,
                                     localIndex const & B,
                                     real64 const & damageA,
                                     real64 const & damageB,
                                     real64 const & maxDamageA,
                                     real64 const & maxDamageB,
                                     arraySlice1d< real64 const > const xA,
                                     arraySlice1d< real64 const > const xB,
                                     real64 const & surfaceQuality,
                                     real64 const & singleFieldStateFraction );

  /// @brief Marks particles whose centers or integration domains leave the supported grid.
  void flagOutOfRangeParticles( ParticleManager & particleManager,
                                SpatialPartition & partition );


  /// @brief Limits overly distorted CPDI integration domains according to the selected policy.
  void cpdiDomainScaling( ParticleManager & particleManager );

  /// @brief Resolves automatic subdivision settings from particle and run configuration.
  void resolveParticleSubdivisionDefault( ParticleManager & particleManager );

  /// @brief Splits eligible particle domains and transfers state to the resulting children.
  bool subdivideParticles( ParticleManager & particleManager,
                           localIndex const rankParticleCountForSplitDecision );

  /// @brief Resizes per-subregion mapping and effective-stencil storage.
  void resizeMappingArrays( ParticleManager & particleManager );

  /// @brief Maps active particles and populates shape values and gradients.
  bool populateMappingArraysForActiveParticles( ParticleManager & particleManager,
                                                NodeManager & nodeManager );

  /// @brief Flags invalid mappings and compacts their active-particle row indices.
  bool flagParticlesWithBadMappingArraysAndCompactActiveOrdinalState( ParticleManager & particleManager,
                                                                      stdVector< array1d< int > > & badMappingRows );

  /// @brief Compacts mapping arrays after invalid active rows are removed.
  void compactActiveOrdinalMappingArrays( localIndex const subRegionIndex,
                                          ParticleSubRegion & subRegion,
                                          arrayView1d< int const > const badMappingRows );

  /// @brief Projects simplex CPTI geometry into an admissible integration configuration.
  /// @return Status bits: 1 scaled, 2 integration measure relaxed, 4 invalid.
  template< int DIM = 3 >
  GEOS_HOST_DEVICE int projectCPTIIntegrationGeometry(
    arraySlice2d< real64 const > const geometry, real64 ( &projected )[3][3] ) const;

  /// @brief Prepares a temporary admissible domain for particle interpolation.
  /// @return Status bits: 1 changed, 2 measure relaxed, 4 invalid.
  template< int DIM >
  GEOS_HOST_DEVICE int prepareIntegrationGeometry(
    ParticleType const type, arraySlice1d< real64 const > const center,
    arraySlice2d< real64 const > const geometry, real64 ( &integration )[3][3] ) const;

  /// @brief Builds a particle stencil using its type, current geometry, and active dimension.
  /// Uses an on-demand projected geometry for CPDI/CPTI; physical geometry is unchanged.
  /// @return True only when a valid nodal stencil was constructed.
  template< int DIM, class NODES, class VALUES, class GRADIENTS >
  GEOS_HOST_DEVICE bool computeParticleMapping(
    ParticleType const type, arraySlice1d< real64 const > const center,
    arraySlice2d< real64 const > const geometry, arrayView3d< localIndex const > const map,
    real64 const (&origin)[3], real64 const (&spacing)[3], localIndex const (&cells)[3],
    NODES nodes, VALUES values, GRADIENTS gradients ) const;

  /// @brief Projects a box-type particle domain to a feasible integration configuration.
  /// @return Status bits distinguish scaled/relaxed/invalid geometry.
  template< int DIM = 3 >
  GEOS_HOST_DEVICE int projectIntegrationGeometry(
    arraySlice2d< real64 const > const geometry, real64 ( &projected )[3][3]) const;

  /// @brief Fills nodal indices, shape values, and gradients for one particle.
  template< int DIM >
  GEOS_FORCE_INLINE
  GEOS_HOST_DEVICE
  void mapNodesAndComputeShapeFunctionsForSingleParticle( arrayView3d< localIndex const > const ijkMap,
                                                         real64 const (&xLocalMin)[3],
                                                         real64 const (&hEl)[3],
                                                         localIndex const (&nEl)[3],
                                                         ParticleType particleType,
                                                         arraySlice1d< real64 const > const particlePosition,
                                                         arraySlice2d< real64 const > const particleGeometryCurrent,
                                                         arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const gridPosition,
                                                         localIndex * const mappedNode,
                                                         real64 * const shapeFunctionValues,
                                                         real64 shapeFunctionGradientValues[][3] );

  /// @brief Applies particle-transform events after deformation-gradient updates.
  void transformParticlesForTriggeredEvents( real64 const time_n,
                                             real64 const dt,
                                             ParticleManager & particleManager );

  /// @brief Assembles particle body forces from the configured accelerations and loading.
  void computeBodyForce( ParticleManager & particleManager );

  /// @brief Computes neighborhood-smoothed Jacobian measures for overlap correction.
  void computeSPHJacobian( ParticleManager & particleManager );

  /// @brief Applies the configured SPH-based particle-overlap correction.
  void sphOverlapCorrection( real64 const dt,
                             ParticleManager & particleManager );

  /// @brief Updates particle thermal energy and temperature when enabled.
  void computeInternalEnergyAndTemperature( const real64 dt,
                                            ParticleManager & particleManager );

  /// @brief Evaluates the body-force field for the generalized-vortex verification problem.
  void computeGeneralizedVortexMMSBodyForce( real64 const time_n,
                                             ParticleManager & particleManager );

  /// @brief Shifts periodic ghost particle centers into the local image used by mapping.
  void correctGhostParticleCentersAcrossPeriodicBoundaries( ParticleManager & particleManager,
                                                            SpatialPartition & partition );

  /// @brief Wraps owned particle centers across periodic global boundaries.
  void correctParticleCentersAcrossPeriodicBoundaries( ParticleManager & particleManager,
                                                       SpatialPartition & partition );

  /// @brief Rebases deformed particles and associated geometry when deformation gradients are reset.
  void resetDeformationGradient( ParticleManager & particleManager );


  /// @brief Reduces particle kinetic energy for diagnostics.
  void computeKineticEnergy( ParticleManager & particleManager );

  /// @brief Converts the configured profile axis to a coordinate index.
  int getProfileDirectionIndex() const;

  /// @brief Determines the profile slice count for a chosen coordinate direction.
  int getProfileNumSlices( int const direction ) const;

  /// @brief Tests the cycle/time conditions for writing a profile sample.
  bool shouldWriteProfiles( real64 const outputTime,
                            int const cycleNumber ) const;

  /// @brief Advances the scheduled profile output time.
  void updateNextProfileWriteTime( real64 const outputTime );

  /// @brief Creates profile output files and their column headers.
  void initializeProfileFiles();

  /// @brief Reduces particle/grid quantities into profile slices and writes them.
  void computeAndWriteProfiles( int const cycleNumber,
                                real64 const time,
                                real64 const dt,
                                ParticleManager & particleManager,
                                NodeManager & nodeManager );

  /// @brief Estimates second spatial derivatives of the particle damage field.
  void computeDamageHessian( ParticleManager & particleManager );

  /// @brief Updates cohesive-related surface classifications.
  void czSurfaceFlagUpdate( ParticleManager & particleManager );

  /// @brief Marks candidate cohesive interfaces adjacent to binder material.
  void tagBinderCZSurfaces( ParticleManager & particleManager,
                            NodeManager & nodeManager );

  /// @brief Tests the cycle/time conditions for writing tracer data.
  bool shouldWriteTracers( real64 const outputTime,
                           int const cycleNumber ) const;

  /// @brief Advances the scheduled tracer output time.
  void updateNextTracerWriteTime( real64 const outputTime );

  /// @brief Selects the persistent particle IDs used for tracer histories.
  void initializeTracerParticleIDs( ParticleManager & particleManager );

  /// @brief Creates tracer data files and their column headers.
  void initializeTracerFiles( ParticleManager & particleManager );

  /// @brief Gathers selected tracer quantities across ranks and writes history records.
  void computeAndWriteTracers( int const cycleNumber,
                               real64 const time,
                               real64 const dt,
                               ParticleManager & particleManager );

  /// @brief Returns mathematical modulo, including nonnegative wrap for negative arguments.
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  real64 Mod( real64 num, real64 denom );

  /// @brief Computes the binomial coefficient for integer arguments.
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  int combinations( int n,
                    int k );

  /// @brief Computes the factorial for a nonnegative integer argument.
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  int factorial( int n );

protected:
  /// @brief Extends recursive input parsing for MPM-specific XML sections.
  void processInputFileRecursive( xmlWrapper::xmlDocument & xmlDocument,
                                  xmlWrapper::xmlNode & targetNode ) override;

  /// @brief Continues recursive parsing from an explicitly supplied XML node position.
  void processInputFileRecursive( xmlWrapper::xmlDocument & xmlDocument,
                                  xmlWrapper::xmlNode & targetNode,
                                  xmlWrapper::xmlNodePos const & targetNodePos ) override;

  /// @brief Validates selected solver options after input parsing.
  virtual void postInputInitialization() override final;

  /// @brief Restores solver state and derived storage after a restart.
  virtual void postRestartInitialization() override final;

  /// @brief Registers constitutive names for a particle subregion through the base interface.
  virtual void setConstitutiveNamesCallSuper( ParticleSubRegionBase & subRegion ) const override;

  // Solver configuration, restartable state, and transient mapping data.
  // Maintain this declaration order when editing constructor initialization.
  mpm::AreaIntegrationOption m_areaIntegrationMethod;
  localIndex m_batchSize;
  array2d< real64 > m_bcTable;
  int m_binSizeMultiplier;
  array1d< real64 > m_bodyForce;
  real64 m_boreholeRadius;
  array1d< real64 > m_boreholeStress;
  array1d< int > m_boundaryConditionTypes; // TODO: Surely there's a way to have just one variable here
  array1d< real64 > m_boundaryFaceCoefficientsOfRestitution;
  array1d< real64 > m_boundaryFaceFrictionCoefficients; // Ignored unless face has boundary condition type 3
  int m_boxAverageHistory;
  array1d< real64 > m_boxAverageMax;
  array1d< real64 > m_boxAverageMin;
  int m_boxAverageResizeWithDomain; // 0 fixed sampling box, 1 follow domain stretch.
  real64 m_boxAverageWriteInterval;
  int m_computeInternalEnergyAndTemperature;
  int m_computeNodalArea;
  int m_computeParticleSurfaceNormalsAndPositions;
  int m_computeSPHJacobian;
  array1d< real64 > m_confiningPressureBoxMax;
  array1d< real64 > m_confiningPressureBoxMin;
  real64 m_confiningPressureSphereRadius;
  int m_uniformConfiningPressure;
  array1d< real64 > m_confiningStress;
  // Contact gap, normal, overlap, and field-partition controls.
  mpm::ContactGapCorrectionOption m_contactGapCorrection;
  real64 m_contactNormalExponent;
  mpm::ContactNormalTypeOption m_contactNormalType;
  int m_cpdiDomainScaling;
  mpm::CPDIDomainScalingTypeOption m_cpdiDomainScalingType;
  real64 m_crackTipDetectionThreshold;
  int m_damageFieldPartitioning;
  real64 m_damageHessianSurfaceThreshold;
  int m_computeCZInterfacesFromDamage;
  int m_directionalOverlapCorrection;
  int m_disableSurfaceNormalsAndPositionsOnCPDIScaling; // Suppress surface geometry after severe CPDI scaling.
  int m_disableSurfaceNormalsAndPositionsOnDamage; // Suppress surface geometry after severe damage.
  int m_disableSurfaceNormalsAndPositionsOnMelt; // Suppress surface geometry for melted particles.
  int m_disableSurfaceNormalsAndPositionsOnOversizedDomainReset; // Suppress surface geometry after F reset.
  int m_disableSurfaceTractionsOnDamage; // Suppress surface tractions after severe damage.
  int m_disableSurfaceTractionsOnMelt; // Suppress surface tractions for melted particles.
  array1d< real64 > m_domainExtent; // Global domain extents without buffers.
  array1d< real64 > m_domainF;
  // Step-local domain stretch derived from domainF; not checkpointed.
  array1d< real64 > m_domainIncrementalStretch;
  array1d< real64 > m_domainL;
  mpm::DomainResetTypeOption m_domainResetType;
  integer m_geometryRestartVersion = 0;
  mpm::GeometryResetOption m_geometryReset = mpm::GeometryResetOption::FOnly;
  array1d< real64 > m_domainStress;
  real64 m_domainTemperature;
  real64 m_domainTemperatureRate;
  // Effective/coalesced stencils use active-particle ordinals, not storage indices.
  stdVector< array2d< integer > > m_effectiveMappedFields;
  stdVector< array2d< localIndex > > m_effectiveMappedNodes;
  stdVector< array3d< real64 > > m_effectiveShapeFunctionGradientValues;
  stdVector< array2d< real64 > > m_effectiveShapeFunctionValues;
  int m_enableBoreholePressure;
  int m_enableConfiningPressure;
  int m_enableContact;
  array1d< int > m_enablePrescribedBoundaryTransverseVelocities;
  int m_enableSurfaceTension;
  int m_eventReporting;
  // Deprecated compatibility input; endpoint split always advances scalar volume.
  int m_exactJIntegration;
  real64 m_explicitSurfaceNormalInfluence;
  real64 m_frictionCoefficient;
  array2d< real64 > m_frictionCoefficientTable;
  // Deprecated compatibility input; the fixed-L exponential uses internal scaling and squaring.
  integer m_FSubcycles;
  int m_flagParticlesWithBadMappingArraysForDeletion;
  array2d< real64 > m_fTable;
  mpm::InterpolationOption m_fTableInterpType;
  int m_generalizedVortexMMS;
  array1d< real64 > m_globalFaceReactions;
  mpm::GPUSchemeOption m_gpuScheme;
  int m_hasContact;
  array1d< real64 > m_hEl; // Cartesian background-grid spacing.
  array3d< localIndex > m_ijkMap; // Structured (i,j,k) indices to local node IDs.
  int m_LBar;
  real64 m_LBarScale;
  int m_logMomentum;
  int m_logStartCycle;
  int m_currentMomentumLogCycle;
  real64 m_currentMomentumLogTime;
  real64 m_currentMomentumLogDt;
  bool m_momentumHistoryInitialized;
  real64 m_LRtolerance;
  stdVector< array2d< integer > > m_mappedFields;
  // Raw particle-to-node indices, arranged as [subregion][particle][stencil entry].
  stdVector< array2d< localIndex > > m_mappedNodes;
  integer m_maxLRIterations;
  integer m_maxNodalNeighbors;
  real64 m_maxParticleJacobian;
  real64 m_maxParticleVelocity;
  real64 m_maxParticleVelocitySquared;
  real64 m_minParticleJacobian;
  int m_needsNeighborList;
  int m_needsNodalNeighborList;
  real64 m_neighborRadius;
  array1d< integer > m_nEl; // Cells per direction, including buffer/ghost cells.
  real64 m_nextBoxAverageWriteTime;
  real64 m_nextParticleDataWriteTime;
  real64 m_nextProfileWriteTime;
  real64 m_nextReactionWriteTime;
  real64 m_nextTracerWriteTime;
  OrderedVariableToManyParticleRelation m_nodalNeighborList;
  mpm::NormalsAndPositionsMethodOption m_normalAndPositionMethod;
  localIndex m_numberOfSubRegions;
  integer m_numContactFlags;
  integer m_numContactGroups;
  integer m_numDims;
  stdVector< array1d< localIndex > > m_numEffectiveMappedNodes;
  integer m_numSurfaceIntegrationPoints;
  integer m_numVelocityFields;
  mpm::OverlapCorrectionOption m_overlapCorrection;
  real64 m_overlapThreshold1;
  real64 m_overlapThreshold2;
  int m_overwriteExistingNormalsAndPositions;
  real64 m_particleDataWriteInterval;
  string m_profileDirection;
  int m_profileCycleInterval;
  int m_profileHistory;
  int m_profileNumSlices;
  string_array m_profileVariables;
  real64 m_profileWriteInterval;
  array1d< real64 > m_partitionExtent; // Local extents including buffer/ghost cells.
  int m_planeStrain; // Legacy selector for native 2D, including axisymmetry.
  bool m_axisymmetric = false;
  string m_formulation = "Cartesian3D";
  // Physical 2D thickness, cached from the partition rather than solver input.
  real64 m_planeStrainThickness = 1.0;

  /// @brief Converts extents to physical Cartesian volume, including 2D thickness.
  /// Axisymmetric volume is handled separately by physicalDomainVolume().
  template< class EXTENTS >
  GEOS_HOST_DEVICE
  real64 physicalVolume( EXTENTS const & extent ) const
  {
    return extent[0] * extent[1] *
           ( m_numDims == 2 ? m_planeStrainThickness : extent[2] );
  }

  /// @brief Returns physical global volume; axisymmetry uses a revolved annulus.
  real64 physicalDomainVolume() const
  {
    return m_axisymmetric ?
           3.14159265358979323846 *
           ( m_xGlobalMax[0] * m_xGlobalMax[0] - m_xGlobalMin[0] * m_xGlobalMin[0] ) * m_domainExtent[1] :
           physicalVolume( m_domainExtent );
  }

  int m_plotGridFields;
  string_array m_plottableFields;
  SortedArray< string > m_plottableFieldsSorted;
  real64 m_polymerCZThickness;

  /// @brief Rejects boundary controls incompatible with the symmetry axis.
  void validateAxisymmetricBoundaryControls() const;
  int m_prescribedBcTable;
  int m_prescribedBoundaryFTable;
  array2d< real64 > m_prescribedBoundaryTransverseVelocities; // Two tangential components per face.
  int m_prescribedFTable;
  int m_preventCZInterpenetration;
  stdVector< std::string > m_profilingLabels;
  stdVector< real64 > m_profilingTimes;
  int m_reactionHistory;
  real64 m_reactionWriteInterval;
  int m_resetDefGradForFullyDamagedParticles;
  int m_resetDefGradForMeltedParticles;
  int m_resetDefGradForScaledSurfaceParticles;
  int m_resetDefGradForOversizedParticles;
  mpm::OversizedParticleTreatmentOption m_oversizedParticleTreatment;
  integer m_oversizedParticleResetRankParticleCountThreshold;
  real64 m_defGradResetMaxParticleDomainToGridCellRatio;
  real64 m_separabilityMinDamage;
  real64 m_maxSingleFieldStateFractionForSeparability;
  int m_setDomainTemperature;
  int m_setDomainTemperatureRate;
  int m_shapeFunctionDiagnostics;
  // Raw, pre-coalescence gradients [subregion][particle][entry][component].
  stdVector< array3d< real64 > > m_shapeFunctionGradientValues;
  // Raw, pre-coalescence weights [subregion][particle][entry].
  stdVector< array2d< real64 > > m_shapeFunctionValues;
  real64 m_smallMass;
  int m_solverProfiling;
  // Stress-control inputs and persistent adaptive-seeking history.
  array1d< int > m_stressControl;
  int m_stressControlAdaptiveSeek;
  int m_stressControlAdaptiveInitialized;
  int m_stressControlWaveWarningIssued;
  int m_stressControlSeekWarningIssued;
  array1d< real64 > m_stressControlITerm;
  real64 m_stressControlKd;
  real64 m_stressControlKi;
  real64 m_stressControlKp;
  array1d< real64 > m_stressControlLastError;
  array1d< real64 > m_stressControlFilteredStress;
  array1d< real64 > m_stressControlWindowStartStress;
  array1d< real64 > m_stressControlPreviousDomainL;
  array1d< real64 > m_stressControlAccumulatedStrain;
  array1d< real64 > m_stressControlSeekStrain;
  array1d< real64 > m_stressControlReengageStrain;
  array1d< int > m_stressControlAdaptiveState;
  array2d< real64 > m_stressControlTangent;
  real64 m_stressControlResponseStrain;
  real64 m_stressControlFilterStrain;
  real64 m_stressControlAdaptStrainWindow;
  real64 m_stressControlAdaptGain;
  real64 m_stressControlMaxRateRatio;
  real64 m_stressControlMaxFeedbackRateRatio;
  real64 m_stressControlSeekRateRatio;
  real64 m_stressControlMaxSeekRateRatio;
  real64 m_stressControlMinStrainRate;
  real64 m_stressControlMaxSeekStrain;
  real64 m_stressControlMaxSeekStrainIncrement;
  real64 m_stressControlCommandFilterStrain;
  real64 m_stressControlMaxRateChangeRatio;
  real64 m_stressControlJammingPressureRatio;
  real64 m_stressControlReengageTangentRatio;
  real64 m_stressControlReengageRampStrain;
  real64 m_stressControlWaveWarnRatio;
  real64 m_stressControlWaveCutbackRatio;
  real64 m_stressControlStressFloor;
  real64 m_stressControlTangentFloorRatio;
  real64 m_stressControlTangentCeilingRatio;
  real64 m_stressControlMaxCouplingRatio;
  real64 m_stressControlSolverDampingRatio;
  array2d< real64 > m_stressTable;
  mpm::InterpolationOption m_stressTableInterpType;
  int m_applyVelocityGradientWhenSplittingParticles;
  int m_subdivideParticles; // -1 auto, 0 off, 1 on.
  int m_subdivideParticlesAutomatic; // Resolved auto-subdivision flag.
  localIndex m_subdivisionRankParticleCountAtStepStart;
  int m_surfaceDetection;
  int m_surfaceHealing;
  real64 m_surfaceNormalAndPositionDamageThreshold;
  real64 m_surfaceQualityThreshold; // [0,1] minimum mapped-normal coherence for same-group DFG.
  real64 m_surfaceTensionForceSign;
  real64 m_surfaceTensionGradientTolerance;
  array2d< real64 > m_surfaceTensionPairs;
  real64 m_thinFeatureDFGThreshold;
  // Integration and grid-to-particle update selection.
  mpm::TimeIntegrationOption m_timeIntegrationOption;
  int m_tracerCycleInterval;
  // Particle tracer selection and output scheduling.
  array2d< real64 > m_tracerCoordinates;
  int m_tracerHistory;
  string m_tracerOutputPrefix;
  array1d< globalIndex > m_tracerParticleIDs;
  string_array m_tracerVariables;
  real64 m_tracerWriteInterval;
  real64 m_totalBinderVolume;
  mpm::GridToParticleMappingOption m_gridToParticleMapping;
  mpm::UpdateMethodOption m_updateMethod;
  integer m_updateOrder;
  int m_useCrackTipDetection;
  int m_useEvents; // Enable MPM event processing.
  int m_useInternalForceAsFaceReaction;
  int m_useNodePosForArea;
  int m_useReferenceVectorsForParticleUpdate;
  int m_useSurfacePositionForContact;
  int m_writeParticleData;
  // Global extents exclude buffers; local extents distinguish ghost inclusion.
  array1d< real64 > m_xGlobalMax;
  array1d< real64 > m_xGlobalMin;
  array1d< real64 > m_xLocalMax;
  array1d< real64 > m_xLocalMaxNoGhost;
  array1d< real64 > m_xLocalMin;
  array1d< real64 > m_xLocalMinNoGhost;

private:
  /// Spatial-bin identity, scoped to particle region and subregion.
  struct BinKey
  {
    localIndex regionIndex;
    localIndex subRegionIndex;
    localIndex binIndex;

    bool operator==( BinKey const & other ) const
    {
      return ( regionIndex == other.regionIndex &&
               subRegionIndex == other.subRegionIndex &&
               binIndex == other.binIndex );
    }
  };

  /// Hash adapter for the three-component spatial-bin key.
  struct BinKeyHash
  {
    std::size_t operator()( BinKey const & k ) const
    {
      using std::size_t;
      using std::hash;

      // Mix the region, subregion, and spatial-bin components of the key.
      return ((std::hash< localIndex >()( k.regionIndex )
               ^ (std::hash< localIndex >()( k.subRegionIndex ) << 1)) >> 1)
             ^ (std::hash< localIndex >()( k.binIndex ) << 1);
    }
  };

  /// @brief Assigns constitutive material names to a particle subregion.
  void setParticlesConstitutiveNames( ParticleSubRegionBase & subRegion ) const;
};

} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_SOLIDMECHANICS_MPM_HPP_ */
