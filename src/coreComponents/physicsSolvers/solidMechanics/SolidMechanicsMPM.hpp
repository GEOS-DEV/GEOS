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

class SpatialPartition;
class SetInitialTemperatureAndPressureMPMEvent;

/**
 * @class SolidMechanicsMPM
 *
 * This class implements a material point method solution to the equations of motion.
 */
class SolidMechanicsMPM : public PhysicsSolverBase
{
public:
  /**
   * Constructor
   * @param name The name of the solver instance
   * @param parent the parent group of the solver
   */
  SolidMechanicsMPM( const string & name,
                     Group * const parent );

  SolidMechanicsMPM( SolidMechanicsMPM const & ) = delete;
  SolidMechanicsMPM( SolidMechanicsMPM && ) = default;

  SolidMechanicsMPM & operator=( SolidMechanicsMPM const & ) = delete;
  SolidMechanicsMPM & operator=( SolidMechanicsMPM && ) = delete;

  /**
   * destructor
   */
  virtual ~SolidMechanicsMPM() override {}

  /**
   * @return The string that may be used to generate a new instance from the
   * PhysicsSolverBase::CatalogInterface::CatalogType
   */
  static string catalogName() { return "SolidMechanics_MPM"; }

  /**
   * @copydoc PhysicsSolverBase::getCatalogName()
   */
  string getCatalogName() const override { return catalogName(); }

  virtual void initializePreSubGroups() override;

  virtual void registerDataOnMesh( Group & meshBodies ) override final;

  /**
   * @defgroup Solver Interface Functions
   *
   * These functions provide the primary interface that is required for derived classes
   */
  /**@{*/
  virtual
  real64 solverStep( real64 const & time_n,
                     real64 const & dt,
                     integer const cycleNumber,
                     DomainPartition & domain ) override;

  virtual
  real64 explicitStep( real64 const & time_n,
                       real64 const & dt,
                       integer const cycleNumber,
                       DomainPartition & domain ) override;

  virtual void updateState( DomainPartition & domain ) override final
  {
    // There should be nothing to update
    GEOS_UNUSED_VAR( domain );
  };

  /**@}*/

  /**
   * This method is called when its host event is triggered
   */
  virtual bool execute( real64 const time_n,
                        real64 const dt,
                        integer const cycleNumber,
                        integer const eventCounter,
                        real64 const eventProgress,
                        DomainPartition & domain ) override;

  template< typename ... PARAMS >
  real64 explicitKernelDispatch( MeshLevel & mesh,
                                 string_array const & targetRegions,
                                 string const & finiteElementName,
                                 real64 const dt,
                                 std::string const & elementListName );

  struct groupKeyStruct : PhysicsSolverBase::groupKeyStruct
  {
    static constexpr char const * cohesiveZoneManagerString() { return "CohesiveZoneManager"; };
    static constexpr char const * mpmEventManagerString() { return "MPMEvents"; };
  };

  struct viewKeyStruct : PhysicsSolverBase::viewKeyStruct
  {
    // Boundary conditions / loading
    static constexpr char const * boundaryNodesString() { return "boundaryNodes"; }
    static constexpr char const * gridSurfaceTensionForceString() { return "gridSurfaceTensionForce"; }

    // Contact
    static constexpr char const * gridContactForceString() { return "gridContactForce"; }
    static constexpr char const * gridRigidBodyFieldKeyString() { return "gridRigidBodyFieldKey"; }
    static constexpr char const * gridRigidBodyFieldColorString() { return "gridRigidBodyFieldColor"; }
    static constexpr char const * gridRigidBodyFieldContactGroupString() { return "gridRigidBodyFieldContactGroup"; }
    static constexpr char const * gridWeakInterfaceTraceActiveString() { return "gridWeakInterfaceTraceActive"; }
    static constexpr char const * gridWeakInterfaceTraceContactSuppressedString() { return "gridWeakInterfaceTraceContactSuppressed"; }
    static constexpr char const * gridWeakInterfaceTraceAnchorWeightString() { return "gridWeakInterfaceTraceAnchorWeight"; }
    static constexpr char const * gridWeakInterfaceTraceForceString() { return "gridWeakInterfaceTraceForce"; }
    static constexpr char const * gridWeakInterfaceTracePointString() { return "gridWeakInterfaceTracePoint"; }
    static constexpr char const * gridWeakInterfaceTraceSupportWeightString() { return "gridWeakInterfaceTraceSupportWeight"; }
    static constexpr char const * gridWeakInterfaceTraceSkipReasonString() { return "gridWeakInterfaceTraceSkipReason"; }
    static constexpr char const * gridWeakInterfaceTraceSurfaceJumpString() { return "gridWeakInterfaceTraceSurfaceJump"; }
    static constexpr char const * gridWeakInterfaceTraceVelocityJumpString() { return "gridWeakInterfaceTraceVelocityJump"; }
    static constexpr char const * gridWeakInterfaceTraceVelocityJumpPostString() { return "gridWeakInterfaceTraceVelocityJumpPost"; }
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
    // Per-step active grid-field mask, computed once from post-P2G mass and synchronized.
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

  struct ExplicitStepManagers
  {
    SpatialPartition & partition;
    arrayView1d< int const > periodic;
    ParticleManager & particleManager;
    MeshLevel & mesh;
    NodeManager & nodeManager;
  };

  void initialize( NodeManager & nodeManager,
                   ParticleManager & particleManager,
                   SpatialPartition & partition );

  static GEOS_FORCE_INLINE
  GEOS_HOST_DEVICE
  localIndex partitionField( int numContactGroups,
                             int damageFieldPartitioning,
                             localIndex particleGroup,
                             arraySlice1d< real64 const > const particleDamageGradient,
                             arraySlice1d< real64 const > const particleSurfaceNormal,
                             arraySlice1d< real64 const > const gridDamageGradient );

  static GEOS_FORCE_INLINE
  GEOS_HOST_DEVICE
  real64 computeSurfaceQualityFromMappingNormalTensor( int const & planeStrain,
                                                       arraySlice1d< real64 const > const mappingNormalTensorA,
                                                       arraySlice1d< real64 const > const mappingNormalTensorB,
                                                       real64 const & materialVolumeA,
                                                       real64 const & materialVolumeB );

  void triggerEvents( const real64 dt,
                      const real64 time_n,
                      ParticleManager & particleManager,
                      SpatialPartition & partition );

  void setInitialTemperatureAndPressure( ParticleManager & particleManager,
                                         SetInitialTemperatureAndPressureMPMEvent const & event );

  void checkEventCompletion( const real64 time_n );

  real64 computeRigidBodyJammingStress( real64 const dt ) const;

  real64 computeRigidBodyContactLengthScale() const;

  void performMaterialSwap( ParticleManager & particleManager,
                            string sourceRegionName,
                            string destinationRegionName );

  void resizeGrid( SpatialPartition & partition,
                   NodeManager & nodeManager,
                   real64 const dt );

  void syncGridFields( stdVector< std::string > const & fieldNames,
                       DomainPartition & domain,
                       NodeManager & nodeManager,
                       MeshLevel & mesh,
                       MPI_Op op );

  void replaceGridFieldsOwnerToGhost( stdVector< std::string > const & fieldNames,
                                      DomainPartition & domain,
                                      NodeManager & nodeManager,
                                      MeshLevel & mesh );

  void synchronizePostBoundaryKinematicFieldsForG2P( DomainPartition & domain,
                                                     NodeManager & nodeManager,
                                                     MeshLevel & mesh );

  /**
   * @brief Computes a per-step active grid-field mask from synchronized grid mass.
   *
   * This avoids repeated mass-threshold decisions after P2G. The mask is
   * synchronized with MPI_MAX so all partitions make the same active/inactive
   * decision for shared grid nodes during the rest of the step.
   */
  void computeActiveGridFieldsForExplicitStep( DomainPartition & domain,
                                               NodeManager & nodeManager,
                                               MeshLevel & mesh );

  // void singleFaceVectorFieldSymmetryBC( const int face,
  //                                       arrayView3d< real64 > const & vectorMultiField,
  //                                       arrayView3d< real64 > const & dVectorMultiField,
  //                                       arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const gridPosition,
  //                                       Group & nodeSets );

  void enforceGridVectorFieldSymmetryBC( arrayView3d< real64 > const & vectorMultiField,
                                         arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const gridPosition,
                                         Group & nodeSets );

  void applyEssentialBCs( const real64 dt,
                          const real64 time_n,
                          NodeManager & nodeManager );

  ExplicitStepManagers getExplicitStepManagers( DomainPartition & domain );

  bool hasPeriodicBoundary( arrayView1d< int const > const periodic ) const;

  void initializeAndResetExplicitStepState( int const cycleNumber,
                                            NodeManager & nodeManager,
                                            ParticleManager & particleManager,
                                            SpatialPartition & partition );

  void prepareParticleTopologyForExplicitStep( real64 const dt,
                                               real64 const time_n,
                                               int const cycleNumber,
                                               DomainPartition & domain,
                                               ParticleManager & particleManager,
                                               SpatialPartition & partition,
                                               arrayView1d< int const > const periodic );

  void buildNeighborhoodAndContactStateForExplicitStep( real64 const dt,
                                                        real64 const time_n,
                                                        int const cycleNumber,
                                                        ParticleManager & particleManager );

  void populateParticleGridMappingForExplicitStep( ParticleManager & particleManager,
                                                  NodeManager & nodeManager );

  void updateDamageAndSurfaceFieldsForExplicitStep( DomainPartition & domain,
                                                    ParticleManager & particleManager,
                                                    NodeManager & nodeManager,
                                                    MeshLevel & mesh );

  void computeParticleLoadsAndBackgroundFieldsForExplicitStep( real64 const time_n,
                                                               ParticleManager & particleManager,
                                                               NodeManager & nodeManager );

  void updateGridDynamicsAndContactForExplicitStep( real64 const dt,
                                                    ParticleManager & particleManager,
                                                    NodeManager & nodeManager );

  /**
   * @brief Applies already-computed essential boundary constraints to the FMPM uncontacted seed velocity.
   *
   * The uncontacted seed is saved before material contact, then constrained after the usual essential-boundary
   * pass so Net contact can distinguish prescribed boundary motion from material-contact impulse.
   */
  void applyFMPMUncontactedVelocityBoundaryConditions( NodeManager & nodeManager );

  /**
   * @brief Applies homogeneous moving/symmetry boundary constraints to one FMPM velocity increment.
   *
   * The first-order FMPM seed contains the prescribed boundary velocity. Higher-order increments must therefore
   * have zero constrained components so the accumulated FMPM velocity preserves the prescribed motion.
   */
  void applyFMPMVelocityIncrementBoundaryConditions( NodeManager & nodeManager,
                                                    arrayView3d< real64 > const velocityIncrement );

  /**
   * @brief Copies constrained boundary components from one FMPM grid velocity field to another.
   *
   * Unconstrained components are left unchanged in the destination field. Buffer-node values are rebuilt from the
   * destination interior field with the same moving-grid mirror relation used by the explicit boundary update.
   */
  void copyConstrainedFMPMBoundaryVelocity( NodeManager & nodeManager,
                                            arrayView3d< real64 > const destinationVelocity,
                                            arrayView3d< real64 const > const sourceVelocity );

  void applyPrescribedDeformationAndBoundaryConditionsForExplicitStep( real64 const dt,
                                                                       real64 const time_n,
                                                                       int const cycleNumber,
                                                                       ParticleManager & particleManager,
                                                                       NodeManager & nodeManager,
                                                                       SpatialPartition & partition );

  void updateParticleKinematicsForExplicitStep( real64 const dt,
                                                real64 const time_n,
                                                ParticleManager & particleManager,
                                                SpatialPartition & partition );

  void updateConstitutiveAndThermalStateForExplicitStep( real64 const dt,
                                                         ParticleManager & particleManager );

  real64 writeOutputsAndComputeStableTimeStepForExplicitStep( real64 const time_n,
                                                              real64 const dt,
                                                              int const cycleNumber,
                                                              ParticleManager & particleManager );

  void resizeGridAndCleanParticlesForExplicitStep( real64 const dt,
                                                   DomainPartition & domain,
                                                   ParticleManager & particleManager,
                                                   NodeManager & nodeManager,
                                                   SpatialPartition & partition,
                                                   arrayView1d< int const > const periodic );

  void logAndProfile( std::string const & label );

  void logAndProfile( std::string const & label,
                      ParticleManager & particleManager,
                      NodeManager & nodeManager );

  void logMomentumSum( std::string const & label, // For tagging code location of output
                       ParticleManager & particleManager,
                       NodeManager & nodeManager );

  void logMomentumSum( std::string const & label, // For tagging code location of output
                       ParticleManager & particleManager,
                       NodeManager & nodeManager,
                       arrayView3d< real64 const > const diagnosticVelocity,
                       std::string const & diagnosticVelocityName );

  void moveParticleWrappersToHost( ParticleManager & particleManager );

  void performExplicitStepParticleGhosting( DomainPartition & domain,
                                            ParticleManager & particleManager,
                                            SpatialPartition & partition );

  void setActiveParticleIndices( ParticleManager & particleManager );

  void updateContactFlagsFromBoundaryConditions();

  void disableSurfaceDataOnDamagedParticles( ParticleManager & particleManager );

  void disableSurfaceDataOnMeltedParticles( ParticleManager & particleManager );

  void updateNodalNeighborListRequirement();

  void zeroParticleCohesiveForces( ParticleManager & particleManager );

  void resetCohesiveSurfaceFlags( ParticleManager & particleManager );

  void updateCohesiveZonesForExplicitStep( real64 const dt,
                                           DomainPartition & domain,
                                           ParticleManager & particleManager,
                                           NodeManager & nodeManager,
                                           MeshLevel & mesh );

  void shapeFunctionDiagnostics( ParticleManager & particleManager );

  void performParticleToGridForExplicitStep( real64 const time_n,
                                             integer const cycleNumber,
                                             ParticleManager & particleManager,
                                             NodeManager & nodeManager );

  void syncGridFieldsForExplicitStep( DomainPartition & domain,
                                      NodeManager & nodeManager,
                                      MeshLevel & mesh );

  void syncSurfaceTensionForcesForExplicitStep( DomainPartition & domain,
                                                NodeManager & nodeManager,
                                                MeshLevel & mesh );

  void computeAndSyncSurfaceTensionForcesForExplicitStep( ParticleManager & particleManager,
                                                          DomainPartition & domain,
                                                          NodeManager & nodeManager,
                                                          MeshLevel & mesh );

  void enforceGridFieldSymmetryAndNormalize( NodeManager & nodeManager );

  void updatePrescribedKinematicsForExplicitStep( real64 const dt,
                                                  real64 const time_n,
                                                  ParticleManager & particleManager,
                                                  SpatialPartition & partition );

  void repartitionParticlesForExplicitStep( DomainPartition & domain,
                                            ParticleManager & particleManager,
                                            SpatialPartition & partition,
                                            arrayView1d< int const > const periodic );

  void applySuperimposedVelocityGradient( const real64 dt,
                                          ParticleManager & particleManager,
                                          SpatialPartition & partition );

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  real64 bruteForceNodalAreaIntegration( real64 const (&hEl)[3],
                                         integer const & numSurfaceIntegrationPoints,
                                         real64 const & L,
                                         real64 const & dA,
                                         arraySlice1d< real64 const > const surfaceNormal,
                                         arraySlice1d< real64 const > const surfacePosition );

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  real64 meshNodalAreaIntegration( real64 const (&hEl)[3],
                                  integer const & numSurfaceIntegrationPoints,
                                  localIndex const & numNeighbors,
                                  arraySlice1d< localIndex const > const neighborRegions,
                                  arraySlice1d< localIndex const > const neighborSubRegions,
                                  arraySlice1d< localIndex const > const neighborIndices,
                                  ParticleManager::ParticleViewConst< arrayView2d< real64 const > > const particleCenter,
                                  ParticleManager::ParticleViewConst< arrayView3d< real64 const > > const particleRVectors,
                                  arraySlice1d< real64 const > const surfaceNormal,
                                  arraySlice1d< real64 const > const surfacePosition );

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
                                                     ParticleManager::ParticleViewConst< arrayView3d< real64 const > > const particleRVectors,
                                                     arraySlice1d< real64 const > const surfaceNormal,
                                                     arraySlice1d< real64 const > const surfacePosition );

  void computeNodalAreas( NodeManager & nodeManager );

  void normalizeGridSurfaceNormalsAndPositions( NodeManager & nodeManager );

  void computeGridSurfaceNormalWeights( ParticleManager & particleManager,
                                        NodeManager & nodeManager );

  void computeContactForces( real64 const dt,
                             ParticleManager & particleManager,
                             NodeManager & nodeManager );

  void enforceWeakInterfaceTraceProjection( real64 const dt,
                                            DomainPartition & domain,
                                            NodeManager & nodeManager,
                                            MeshLevel & mesh );

  void zeroWeakInterfaceTraceProjectionFields( NodeManager & nodeManager );

  void computeWeakInterfaceTraceProjectionForces( real64 const dt,
                                                  NodeManager & nodeManager );

  void applyWeakInterfaceTraceProjectionForces( real64 const dt,
                                                NodeManager & nodeManager );

  void computeWeakInterfaceTraceProjectionPostVelocityJump( NodeManager & nodeManager );

  void updateInterfaceMechanismDiagnostics( NodeManager & nodeManager );

  /**
   * @brief Computes the cumulative material-contact impulse target for FMPM Net contact.
   *
   * The target is evaluated from the accumulated uncorrected FMPM velocity. The caller applies only the
   * difference between this target and the previously applied cumulative impulse to the current increment.
   */
  void computeFMPMNetContactMomentumTarget( real64 const dt,
                                            ParticleManager & particleManager,
                                            NodeManager & nodeManager,
                                            arrayView3d< real64 const > const vUncorrectedTotal,
                                            arrayView3d< real64 > const contactMomentumTarget );

  void initializeFrictionCoefficients();

  /**
   * @brief Computes the impulse form of the pairwise nodal contact law.
   *
   * This helper mirrors computePairwiseNodalContactForce, but returns contact impulse directly so FMPM Net
   * contact can accumulate and compare total momentum corrections between FMPM orders.
   */
  static GEOS_FORCE_INLINE
  GEOS_HOST_DEVICE
  void computePairwiseNodalContactImpulse( mpm::ContactGapCorrectionOption const & contactGapCorrection,
                                           mpm::OverlapCorrectionOption const & overlapCorrection,
                                           real64 const overlapThreshold1,
                                           real64 const overlapThreshold2,
                                           real64 const maxParticleVelocitySquared,
                                           real64 const (&hEl) [3],
                                           int const & planeStrain,
                                           real64 const & smallMass,
                                           int const & useSurfacePositionForContact,
                                           int const & useCohesiveTangentialForces,
                                           real64 const rigidBodyPenetrationPenaltyBeta,
                                           real64 & contactPenetration,
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

  static GEOS_FORCE_INLINE
  GEOS_HOST_DEVICE
  void computePairwiseNodalContactForce( mpm::ContactGapCorrectionOption const & contactGapCorrection,
                                         mpm::OverlapCorrectionOption const & overlapCorrection,
                                         real64 const overlapThreshold1,
                                         real64 const overlapThreshold2,
                                         real64 const maxParticleVelocitySquared,
                                         real64 const (&hEl) [3],
                                         int const & planeStrain,
                                         real64 const & smallMass,
                                         int const & useSurfacePositionForContact,
                                         int const & useCohesiveTangentialForces,
                                         real64 const rigidBodyPenetrationPenaltyBeta,
                                         real64 & contactPenetration,
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
                                         arraySlice1d< real64 const > const xA, // Position of field A
                                         arraySlice1d< real64 const > const xB, // Position of field B
                                         arraySlice1d< real64 const > const centerOfVolumeA,
                                         arraySlice1d< real64 const > const centerOfVolumeB,
                                         arraySlice1d< real64 > const fA,
                                         arraySlice1d< real64 > const fB );

  static GEOS_FORCE_INLINE
  GEOS_HOST_DEVICE
  void computeOrthonormalBasis( const real64 * e1,  // input "normal" unit vector.
                                real64 * e2,        // output "tangential" unit vector.
                                real64 * e3 );      // output "tangential" unit vector.

  void setGridFieldLabels( NodeManager & nodeManager );

  real64 computeNeighborList( ParticleManager & particleManager );

  void generateNodalNeighborList( ParticleManager & particleManager,
                                  NodeManager & nodeManager );

  void optimizeBinSort( ParticleManager & particleManager );

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  real64 kernel( real64 const & r ); // distance from particle to query point

  GEOS_HOST
  GEOS_FORCE_INLINE
  real64 inverseKernel( real64 const & d );

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  real64 kernelDevice( real64 const & r,
                       real64 const & neighborRadius,
                       int const & planeStrain ); // distance from particle to query point

  void kernelGradient( arraySlice1d< real64 const > const & x,  // query point
                       arraySlice1d< real64 const > const & xp, // particle location
                       real64 const & r,                        // distance from particle to query point
                       real64 ( &result )[3] );

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  void kernelGradientDevice( arraySlice1d< real64 const > const & x,  // query point
                             real64 const (&xp)[3],  //arraySlice1d< real64 const > const & xp, // particle location
                             real64 const & r,                        // distance from particle to query point
                             real64 const & neighborRadius,
                             int const & planeStrain,
                             real64 ( &result )[3] );

  GEOS_HOST
  GEOS_FORCE_INLINE
  real64 computeKernelField( arraySlice1d< real64 const > const x,    // query point
                             localIndex const numNeighbors,
                             arrayView2d< real64 const > const & xp,    // List of neighbor particle locations.
                             arrayView1d< real64 const > const & Vp,    // List of neighbor particle volumes.
                             arrayView1d< real64 const > const & fp );  // scalar field values (e.g. damage) at neighbor particles

  template< typename DAMAGE >
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  real64 computeKernelFieldDevice( arraySlice1d< real64 const > const x,  // query point
                                   real64 const & neighborRadius,
                                   int const & planeStrain,
                                   localIndex const & numNeighbors,
                                   arraySlice1d< localIndex const > const regionIndices,
                                   arraySlice1d< localIndex const > const subRegionIndices,
                                   arraySlice1d< localIndex const > const particleIndices,
                                   ParticleManager::ParticleView< arrayView1d< real64 const > > particleVolumeView,
                                   ParticleManager::ParticleView< arrayView2d< real64 const > > particlePositionView,
                                   DAMAGE const & damage );

  GEOS_HOST
  GEOS_FORCE_INLINE
  void computeKernelFieldGradient( arraySlice1d< real64 const > const x, // query point
                                   localIndex const numNeighbors,
                                   arrayView2d< real64 const > const & xp,           // List of neighbor particle locations.
                                   arrayView1d< real64 const > const & Vp,           // List of neighbor particle volumes.
                                   arrayView1d< real64 const > const & fp,           // scalar field values (e.g. damage) at neighbor
                                                                                     // particles
                                   arraySlice1d< real64 > const result );

  template< typename DAMAGE >
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  void computeKernelFieldGradientDevice( arraySlice1d< real64 const > const x, // query point
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

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  void computeKernelFieldDivergenceDevice( arraySlice1d< real64 const > const x, // query point
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

  GEOS_HOST
  GEOS_FORCE_INLINE
  void computeKernelVectorGradient( arraySlice1d< real64 const > const x,       // query point
                                    localIndex const numNeighbors,
                                    arrayView2d< real64 const > const & xp,  // List of neighbor particle locations.
                                    arrayView1d< real64 const > const & Vp,                 // List of neighbor particle volumes.
                                    arrayView2d< real64 const > const & fp,  // vector field values (e.g. velocity) at neighbor particles
                                    arraySlice2d< real64 > const result );

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  void computeKernelVectorGradientDevice( arraySlice1d< real64 const > const x,       // query point
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

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  bool markSurfaceAsDamage( integer const & surfaceFlag );

  void computeDamageFieldGradient( ParticleManager & particleManager );

  void computeDistanceToCrackTip( ParticleManager & particleManager );

  void estimateParticleSurfacePosition( ParticleManager & particleManager );

  void updateSurfaceFlagOverload( ParticleManager & particleManager );

  void projectDamageFieldGradientToGrid( DomainPartition & domain,
                                         ParticleManager & particleManager,
                                         NodeManager & nodeManager,
                                         MeshLevel & mesh );

  void computeParticleFieldMappings( DomainPartition & domain,
                                     ParticleManager & particleManager,
                                     NodeManager & nodeManager,
                                     MeshLevel & mesh );

  void computeRigidBodyColorFieldMappings( DomainPartition & domain,
                                           ParticleManager & particleManager,
                                           NodeManager & nodeManager,
                                           MeshLevel & mesh );

  void rigidBodyParticleUpdate( real64 const time_n,
                                real64 const dt,
                                ParticleManager & particleManager,
                                NodeManager & nodeManager,
                                SpatialPartition & partition );

  void resetRigidBodyMPIState();

  void initializeRigidBodyMPIState( ParticleManager & particleManager,
                                    SpatialPartition & partition );

  real64 computeRigidBodyStableTimeStep( real64 localMaximumParticleSpeed,
                                         real64 localMinimumCellLength ) const;

  static GEOS_FORCE_INLINE
  GEOS_HOST_DEVICE
  bool isRigidBodySurfaceFlag( integer const surfaceFlag )
  {
    return surfaceFlag == mpm::toInteger( mpm::SurfaceFlag::Surface ) ||
           surfaceFlag == mpm::toInteger( mpm::SurfaceFlag::Cohesive ) ||
           surfaceFlag == mpm::toInteger( mpm::SurfaceFlag::DamagedCohesive ) ||
           surfaceFlag == mpm::toInteger( mpm::SurfaceFlag::WeakDiscontinuity );
  }

  void updateDeformationGradient( real64 dt,
                                  ParticleManager & particleManager );

  void stressControl( const real64 dt,
                      ParticleManager & particleManager,
                      SpatialPartition & partition );

  void adaptiveSeekStressControl( const real64 dt,
                                  real64 const (& targetStress)[3],
                                  arrayView1d< real64 const > const currentStress,
                                  real64 const boxMaterialVolume,
                                  real64 const maximumBulkModulus );

  void initializeParticleFields( ParticleManager & particleManager );

  void initializeConstitutiveModelDependencies( ParticleManager & particleManager );

  void updateConstitutiveModelDependencies( ParticleManager & particleManager );

  void updateStress( real64 dt,
                     ParticleManager & particleManager );

  void particleKinematicUpdate( const real64 dt,
                                ParticleManager & particleManager );

  void computeAndWriteBoxAverage( const real64 dt,
                                  const real64 time_n,
                                  ParticleManager & particleManager );

  void writeParticleData( const real64 time_n,
                          ParticleManager & particleManager );

  void computeBoxMetrics( ParticleManager & particleManager,
                          arrayView1d< real64 > boxStress,
                          real64 & boxMaterialVolume );

  void initializeGridFields( NodeManager & nodeManager );

  void boundaryConditionUpdate( real64 dt, real64 time_n );

  void projectParticleSurfaceNormalsToGrid( DomainPartition & domain,
                                            ParticleManager & particleManager,
                                            NodeManager & nodeManager,
                                            MeshLevel & mesh );

  void initializeCohesiveReferenceConfiguration( DomainPartition & domain,
                                                 ParticleManager & particleManager,
                                                 NodeManager & nodeManager,
                                                 MeshLevel & mesh );

  bool interiorToParticleProjectedArea( ParticleManager & particleManager,
                                        globalIndex const GEOS_UNUSED_PARAM( gridIndex ),
                                        int const gridFieldIndex,
                                        real64 const (&gridSurfaceNormal)[3],
                                        real64 const (&gridSurfacePoint)[3] );

  void projectToPlane( real64 const (&vector)[3],
                       real64 const (&normal)[3],
                       real64 ( &projection )[3] );

  void enforceCohesiveLaw( real64 dt,
                           ParticleManager & particleManager,
                           NodeManager & nodeManager,
                           SpatialPartition & partition );

  GEOS_HOST_DEVICE
  real64 computeDistanceToParticleSurface( real64 ( &normal )[3],
                                           arraySlice2d< real64 const > const rVectors );

  void updateGridBackgroundStress( NodeManager & nodeManager );

  void particleToGrid( real64 const time_n,
                       integer const cycleNumber,
                       ParticleManager & particleManager,
                       NodeManager & nodeManager );

  // XXX Legacy P2G performance variants were usually worse than atomics.
  // XXX Retained only as inactive reference; delete when the archive is no longer useful.

/**
   * @brief Assembles diffuse-interface surface-tension forces with a weak-form P2G pass.
   *
   * The method uses synchronized gridMaterialVolume as a partitioned phase color field,
   * gathers the color gradient to particles with shape-function gradients, builds a
   * continuum surface stress at each particle, and scatters -V_p tau_s grad N into
   * gridSurfaceTensionForce.
   */
  void computePartitionedSurfaceTensionForces( ParticleManager & particleManager,
                                               NodeManager & nodeManager );

  void gridTrialUpdate( real64 dt,
                        NodeManager & nodeManager );

  void enforceContact( real64 dt,
                       //  DomainPartition & domain,
                       ParticleManager & particleManager,
                       NodeManager & nodeManager
                       //  MeshLevel & mesh
                       );

  void computeParticleSurfaceNormalsAndPositions( ParticleManager & particleManager,
                                                  NodeManager & nodeManager );

  void computePairwiseLogisticRegressionSurfaceNormalsAndPositions( ParticleManager & particleManager,
                                                                    NodeManager & nodeManager );

  void computeGridSurfacePositionFromVolume( NodeManager & nodeManager );

  static GEOS_HOST_DEVICE GEOS_FORCE_INLINE
  real64 computeVolumeFromSurfacePosition( real64 const (&n)[3],
                                           real64 const offset,
                                           real64 const (&hEl) [3] );

  static GEOS_HOST_DEVICE GEOS_FORCE_INLINE
  real64 computeMaximumSurfacePositionOffset( real64 ( &n )[3],
                                              real64 const (&hEl)[3] );

  static GEOS_HOST_DEVICE
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

  void mapSurfaceNormalsAndPositionsToParticles( ParticleManager & particleManager,
                                                 NodeManager & nodeManager );

  void interpolateTable( real64 x,
                         real64 dx,
                         arrayView2d< real64 const > const & table,
                         arrayView1d< real64 > output,
                         arrayView1d< real64 > outputRate,
                         mpm::InterpolationOption const & interpolationType );
  
  void interpolateValueInRange( real64 const & x,
                                real64 const & xmin,
                                real64 const & xmax,
                                real64 const & ymin,
                                real64 const & ymax,
                                real64 & output,
                                mpm::InterpolationOption interpolationType );

  void interpolateFTable( real64 dt, real64 time_n );

  void interpolateStressTable( real64 dt, real64 time_n );

  void gridToParticle( real64 dt,
                       ParticleManager & particleManager,
                       NodeManager & nodeManager,
                       DomainPartition & domain,
                       MeshLevel & mesh );

  void performFLIPUpdate( real64 dt,
                          ParticleManager & particleManager,
                          NodeManager & nodeManager );

  void performPICUpdate( real64 dt,
                         ParticleManager & particleManager,
                         NodeManager & nodeManager );

  void performXPICUpdate( real64 dt,
                          ParticleManager & particleManager,
                          NodeManager & nodeManager,
                          DomainPartition & domain,
                          MeshLevel & mesh );

  /**
   * @brief Performs the revised incremental FMPM particle update.
   *
   * The method builds v^{+(k)} from round-trip grid-particle-grid velocity increments, constrains each
   * higher-order increment on moving/symmetry boundaries, optionally applies FMPM Net material contact, and maps
   * the corrected velocity back to particles for velocity, position, and velocity-gradient updates.
   */
  void performFMPMUpdate( real64 dt,
                          ParticleManager & particleManager,
                          NodeManager & nodeManager,
                          DomainPartition & domain,
                          MeshLevel & mesh );

  /**
   * @brief Applies one FMPM Net material-contact correction to the current velocity increment.
   *
   * The update is Delta v <- Delta v + ( J_target - J_previous ) / M. After the correction, J_previous is replaced
   * by J_target.
   */
  void applyFMPMNetContactCorrection( real64 const dt,
                                      ParticleManager & particleManager,
                                      NodeManager & nodeManager,
                                      arrayView3d< real64 const > const vUncorrectedTotal,
                                      arrayView3d< real64 > const contactMomentumNet,
                                      arrayView3d< real64 > const contactMomentumTarget,
                                      arrayView3d< real64 > const velocityIncrement );

  void updateSolverDependencies( ParticleManager & particleManager );

  real64 getStableTimeStep( ParticleManager & particleManager );

  void deleteBadParticles( ParticleManager & particleManager );

  void printProfilingResults();

  void computeSurfaceFlags( ParticleManager & particleManager );

  void computeSphF( ParticleManager & particleManager );

  static GEOS_FORCE_INLINE
  GEOS_HOST_DEVICE
  bool evaluateSeparabilityCriterion( int const & planeStrain,
                                     int const & numContactGroups,
                                     int const & treatFullyDamagedAsSingleField,
                                     real64 const & separabilityMinDamage,
                                     real64 const & thinFeatureDFGThreshold,
                                     real64 const & neighborRadius,
                                     real64 const & surfaceQualityThreshold,
                                     real64 const (&hEl)[3],
                                     localIndex const & A,
                                     localIndex const & B,
                                     real64 const & damageA,
                                     real64 const & damageB,
                                     real64 const & maxDamageA,
                                     real64 const & maxDamageB,
                                     arraySlice1d< real64 const > const damageGradient,
                                     arraySlice1d< real64 const > const xA,
                                     arraySlice1d< real64 const > const xB,
                                     real64 const & surfaceQuality );

  void flagOutOfRangeParticles( ParticleManager & particleManager,
                                SpatialPartition & partition );

  void computeRVectors( ParticleManager & particleManager );

  void cpdiDomainScaling( ParticleManager & particleManager );

  void subdivideParticles( ParticleManager & particleManager );

  void resizeMappingArrays( ParticleManager & particleManager );

  bool populateMappingArraysForActiveParticles( ParticleManager & particleManager,
                                                NodeManager & nodeManager );

  bool flagParticlesWithBadMappingArraysAndCompactActiveOrdinalState( ParticleManager & particleManager,
                                                                      stdVector< array1d< int > > & badMappingRows );

  void compactActiveOrdinalMappingArrays( localIndex const subRegionIndex,
                                          ParticleSubRegion & subRegion,
                                          arrayView1d< int const > const badMappingRows );

  GEOS_FORCE_INLINE
  GEOS_HOST_DEVICE
  void invalidateShapeFunctionRow( localIndex const numberOfMappedNodesPerParticle,
                                   arraySlice1d< int > const mappedNodes,
                                   arraySlice1d< real64 > const shapeFunctionValues,
                                   arraySlice2d< real64 > const shapeFunctionGradientValues );

  GEOS_FORCE_INLINE
  GEOS_HOST_DEVICE
  bool computeSinglePointShapeFunctions( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const gridPosition,
                                         arraySlice1d< real64 const > const particlePosition,
                                         arrayView3d< int const > const ijkMap,
                                         real64 const (&xLocalMin)[3],
                                         real64 const (&hEl)[3],
                                         localIndex const (&nEl)[3],
                                         arraySlice1d< int > const mappedNodes,
                                         arraySlice1d< real64 > const shapeFunctionValues,
                                         arraySlice2d< real64 > const shapeFunctionGradientValues );

  GEOS_FORCE_INLINE
  GEOS_HOST_DEVICE
  bool computeSinglePointBSplineShapeFunctions( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const gridPosition,
                                                arraySlice1d< real64 const > const particlePosition,
                                                arrayView3d< int const > const ijkMap,
                                                real64 const (&xLocalMin)[3],
                                                real64 const (&hEl)[3],
                                                localIndex const (&nEl)[3],
                                                arraySlice1d< int > const mappedNodes,
                                                arraySlice1d< real64 > const shapeFunctionValues,
                                                arraySlice2d< real64 > const shapeFunctionGradientValues );

  GEOS_FORCE_INLINE
  GEOS_HOST_DEVICE
  bool computeCPDIShapeFunctions( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const gridPosition,
                                  arraySlice1d< real64 const > const particlePosition,
                                  arraySlice2d< real64 const > const particleRVectors,
                                  arrayView3d< int const > const ijkMap,
                                  real64 const (&xLocalMin)[3],
                                  real64 const (&hEl)[3],
                                  localIndex const (&nEl)[3],
                                  arraySlice1d< int > const mappedNodes,
                                  arraySlice1d< real64 > const shapeFunctionValues,
                                  arraySlice2d< real64 > const shapeFunctionGradientValues );

  GEOS_FORCE_INLINE
  GEOS_HOST_DEVICE
  void mapNodesAndComputeShapeFunctionsForSingleParticle( arrayView3d< localIndex const > const ijkMap,
                                         real64 const (&xLocalMin)[3],
                                         real64 const (&hEl)[3],
                                         ParticleType particleType,
                                         arraySlice1d< real64 const > const particlePosition,
                                         arraySlice2d< real64 const > const particleRVectors,
                                         arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const gridPosition,
                                         localIndex * const mappedNode,
                                         real64 * const shapeFunctionValues,
                                         real64 shapeFunctionGradientValues[][3] );

  GEOS_FORCE_INLINE
  GEOS_HOST_DEVICE
  void computeSinglePointParticleShapeFunctions( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const gridPosition,
                                                 arraySlice1d< real64 const > const particlePosition,
                                                 arrayView3d< int const > const ijkMap,
                                                 real64 const (&xLocalMin)[3],
                                                 real64 const (&hEl)[3],
                                                 localIndex * const mappedNodes,
                                                 real64 * const shapeFunctionValues,
                                                 real64 shapeFunctionGradientValues[][3] );

  GEOS_FORCE_INLINE
  GEOS_HOST_DEVICE
  void computeSinglePointBSplineParticleShapeFunctions( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const gridPosition,
                                                        arraySlice1d< real64 const > const particlePosition,
                                                        arrayView3d< int const > const ijkMap,
                                                        real64 const (&xLocalMin)[3],
                                                        real64 const (&hEl)[3],
                                                        localIndex * const mappedNodes,
                                                        real64 * const shapeFunctionValues,
                                                        real64 shapeFunctionGradientValues[][3] );

  GEOS_FORCE_INLINE
  GEOS_HOST_DEVICE
  void computeCPDIParticleShapeFunctions( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const gridPosition,
                                          arraySlice1d< real64 const > const particlePosition,
                                          arraySlice2d< real64 const > const particleRVectors,
                                          arrayView3d< int const > const ijkMap,
                                          real64 const (&xLocalMin)[3],
                                          real64 const (&hEl)[3],
                                          localIndex * const mappedNodes,
                                          real64 * const shapeFunctionValues,
                                          real64 shapeFunctionGradientValues[][3] );

  #ifdef USEOLDSHAPEFUNCTIONFUNCTIONS
  GEOS_FORCE_INLINE
  GEOS_HOST_DEVICE
  void computeSinglePointParticleShapeFunctionsSlower( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const gridPosition,
                                                 arraySlice1d< real64 const > const particlePosition,
                                                 arrayView3d< int const > const ijkMap,
                                                 real64 const (&xLocalMin)[3],
                                                 real64 const (&hEl)[3],
                                                 localIndex * const mappedNodes,
                                                 real64 * const shapeFunctionValues,
                                                 real64 shapeFunctionGradientValues[][3] );

  GEOS_FORCE_INLINE
  GEOS_HOST_DEVICE
  void computeCPDIParticleShapeFunctionsSlower( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const gridPosition,
                                          arraySlice1d< real64 const > const particlePosition,
                                          arraySlice2d< real64 const > const particleRVectors,
                                          arrayView3d< int const > const ijkMap,
                                          real64 const (&xLocalMin)[3],
                                          real64 const (&hEl)[3],
                                          localIndex * const mappedNodes,
                                          real64 * const shapeFunctionValues,
                                          real64 shapeFunctionGradientValues[][3] );
  #endif

  /**
   * @brief Applies TransformParticles events that must run after deformation-gradient updates.
   *
   * @param time_n Beginning-of-step time.
   * @param dt Time-step size.
   * @param particleManager Particle manager containing active particle subregions.
   */
  void transformParticlesForTriggeredEvents( real64 const time_n,
                                             real64 const dt,
                                             ParticleManager & particleManager );

  void computeBodyForce( ParticleManager & particleManager );

  void computeSPHJacobian( ParticleManager & particleManager );

  void sphOverlapCorrection( real64 const dt,
                             ParticleManager & particleManager );

  void computeInternalEnergyAndTemperature( const real64 dt,
                                            ParticleManager & particleManager );

  void computeGeneralizedVortexMMSBodyForce( real64 const time_n,
                                             ParticleManager & particleManager );

  void correctGhostParticleCentersAcrossPeriodicBoundaries( ParticleManager & particleManager,
                                                            SpatialPartition & partition );

  void correctParticleCentersAcrossPeriodicBoundaries( ParticleManager & particleManager,
                                                       SpatialPartition & partition );

  void resetDeformationGradient( ParticleManager & particleManager );

  void unscaleCPDIVectors( ParticleManager & particleManager );

  void computeKineticEnergy( ParticleManager & particleManager );

  int getProfileDirectionIndex() const;

  int getProfileNumSlices( int const direction ) const;

  bool shouldWriteProfiles( real64 const outputTime,
                            int const cycleNumber ) const;

  void updateNextProfileWriteTime( real64 const outputTime );

  void initializeProfileFiles();

  void computeAndWriteProfiles( int const cycleNumber,
                                real64 const time,
                                real64 const dt,
                                ParticleManager & particleManager,
                                NodeManager & nodeManager );

  void computeDamageHessian( ParticleManager & particleManager );

  void czSurfaceFlagUpdate( ParticleManager & particleManager );

  void tagBinderCZSurfaces( ParticleManager & particleManager,
                            NodeManager & nodeManager );
  bool shouldWriteTracers( real64 const outputTime,
                           int const cycleNumber ) const;

  void updateNextTracerWriteTime( real64 const outputTime );

  void initializeTracerParticleIDs( ParticleManager & particleManager );

  void initializeTracerFiles( ParticleManager & particleManager );

  void computeAndWriteTracers( int const cycleNumber,
                               real64 const time,
                               real64 const dt,
                               ParticleManager & particleManager );

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  real64 Mod( real64 num, real64 denom );

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  int combinations( int n,
                    int k );

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  int factorial( int n );

protected:
  void processInputFileRecursive( xmlWrapper::xmlDocument & xmlDocument,
                                  xmlWrapper::xmlNode & targetNode );

  void processInputFileRecursive( xmlWrapper::xmlDocument & xmlDocument,
                                  xmlWrapper::xmlNode & targetNode,
                                  xmlWrapper::xmlNodePos const & targetNodePos );

  virtual void postInputInitialization() override final;

  virtual void postRestartInitialization() override final;

  virtual void setConstitutiveNamesCallSuper( ParticleSubRegionBase & subRegion ) const override;

  // Member fields are ordered alphabetically by member name to match the constructor initializer list.
  mpm::AreaIntegrationOption m_areaIntegrationMethod;
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
  int m_boxAverageResizeWithDomain;  // 0 use constant box_average domain, 1: resize with domainL
  real64 m_boxAverageWriteInterval;
  int m_computeInternalEnergyAndTemperature;
  int m_computeNodalArea;
  int m_computeParticleSurfaceNormalsAndPositions;
  int m_computeSPHJacobian;
  array1d< real64 > m_confiningPressureBoxMax;
  array1d< real64 > m_confiningPressureBoxMin;
  array1d< real64 > m_confiningStress;
  mpm::ContactGapCorrectionOption m_contactGapCorrection;
  real64 m_contactNormalExponent;
  mpm::ContactNormalTypeOption m_contactNormalType;
  int m_cpdiDomainScaling;
  real64 m_crackTipDetectionThreshold;
  int m_damageFieldPartitioning;
  real64 m_damageHessianSurfaceThreshold;
  int m_computeCZInterfacesFromDamage;
  int m_directionalOverlapCorrection;
  int m_disableSurfaceNormalsAndPositionsOnCPDIScaling; // Turns off surface normals and positions for highly deformed particles
  int m_disableSurfaceNormalsAndPositionsOnDamage; // Turns off surface normals and positions for highly damaged particles
  int m_disableSurfaceNormalsAndPositionsOnMelt; // Turns off surface normals and positions for melted particles
  int m_disableSurfaceTractionsOnDamage; // Turns off surface tractions for highly damaged particles
  int m_disableSurfaceTractionsOnMelt; // Turns off surface tractions for melted particles
  array1d< real64 > m_domainExtent;       // Length of each edge of global domain excluding buffer cells
  array1d< real64 > m_domainF;
  array1d< real64 > m_domainL;
  array1d< real64 > m_domainStress;
  real64 m_domainTemperature;
  real64 m_domainTemperatureRate;
  stdVector< array2d< integer > > m_effectiveMappedFields;
  stdVector< array2d< localIndex > > m_effectiveMappedNodes;
  stdVector< array3d< real64 > > m_effectiveShapeFunctionGradientValues;
  stdVector< array2d< real64 > > m_effectiveShapeFunctionValues;
  int m_enableBoreholePressure;
  int m_enableConfiningPressure;
  int m_enableContact;
  int m_enableWeakInterfaceTraceProjection;
  array1d< int > m_enablePrescribedBoundaryTransverseVelocities;
  int m_enableSurfaceTension;
  int m_eventReporting;
  int m_exactJIntegration;
  real64 m_explicitSurfaceNormalInfluence;
  real64 m_frictionCoefficient;
  array2d< real64 > m_frictionCoefficientTable;
  int m_FSubcycles;
  int m_flagParticlesWithBadMappingArraysForDeletion;
  array2d< real64 > m_fTable;
  mpm::InterpolationOption m_fTableInterpType;
  int m_generalizedVortexMMS;
  array1d< real64 > m_globalFaceReactions;
  mpm::GPUSchemeOption m_gpuScheme;
  int m_hasContact;
  array1d< real64 > m_hEl;                // Grid spacing in x-y-z
  array3d< localIndex > m_ijkMap;        // Map from indices in each spatial dimension to local node ID
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
  stdVector< array2d< localIndex > > m_mappedNodes; // mappedNodes[subregion index][particle index][node index]. dims = {# of subregions,
                                                    // # of particles, # of nodes a particle on the subregion maps to}
  int m_maxLRIterations;
  int m_maxNodalNeighbors;
  real64 m_maxParticleJacobian;
  real64 m_maxParticleVelocity;
  real64 m_maxParticleVelocitySquared;
  real64 m_minParticleJacobian;
  int m_needsNeighborList;
  int m_needsNodalNeighborList;
  real64 m_neighborRadius;
  array1d< int > m_nEl;                   // Number of elements in each grid direction including buffer and ghost cells
  real64 m_nextBoxAverageWriteTime;
  real64 m_nextParticleDataWriteTime;
  real64 m_nextProfileWriteTime;
  real64 m_nextReactionWriteTime;
  real64 m_nextRigidBodyHistoryWriteTime;
  real64 m_nextTracerWriteTime;
  OrderedVariableToManyParticleRelation m_nodalNeighborList;
  mpm::NormalsAndPositionsMethodOption m_normalAndPositionMethod;
  localIndex m_numberOfSubRegions;
  int m_numContactFlags;
  int m_numContactGroups;
  int m_numDims;
  stdVector< array1d< localIndex > > m_numEffectiveMappedNodes;
  int m_numSurfaceIntegrationPoints;
  int m_numVelocityFields;
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
  array1d< real64 > m_partitionExtent;    // Length of each edge of partition including buffer and ghost cells
  int m_planeStrain;
  int m_plotGridFields;
  string_array m_plottableFields;
  SortedArray< string > m_plottableFieldsSorted;
  int m_plotUnscaledParticles;
  real64 m_polymerCZThickness;
  int m_prescribedBcTable;
  int m_prescribedBoundaryFTable;
  array2d< real64 > m_prescribedBoundaryTransverseVelocities; // 2 in-plane directions * 6 faces
  int m_prescribedFTable;
  int m_preventCZInterpenetration;
  stdVector< std::string > m_profilingLabels;
  stdVector< real64 > m_profilingTimes;
  int m_reactionHistory;
  real64 m_reactionWriteInterval;
  string m_rigidBodyActiveEventName;
  array1d< globalIndex > m_rigidBodyAnchorParticleIDs;
  real64 m_rigidBodyAngularDamping;
  array1d< integer > m_rigidBodyColors;
  array1d< integer > m_rigidBodyContactGroups;
  real64 m_rigidBodyContactCFL;
  real64 m_rigidBodyContactLengthScale;
  int m_rigidBodyHistory;
  real64 m_rigidBodyHistoryWriteInterval;
  real64 m_rigidBodyKineticEnergy;
  real64 m_rigidBodyLinearDamping;
  int m_rigidBodyMaxGridFields;
  real64 m_rigidBodyMaxForce;
  real64 m_rigidBodyMaxTimeStep;
  int m_rigidBodyMode;
  real64 m_rigidBodyObservedJammingStress;
  real64 m_rigidBodyObservedMaxForce;
  real64 m_rigidBodyObservedMaxPenetration;
  real64 m_rigidBodyPenetrationPenaltyBeta;
  array1d< integer > m_rigidBodyPeriodicDirections;
  int m_rigidBodyRegistryInitialized;
  int m_rigidBodyRegistrySynchronized;
  real64 m_rigidBodyStopKineticEnergy;
  array1d< real64 > m_rigidBodyUnwrappedCenters;
  int m_resetDefGradForFullyDamagedParticles;
  int m_resetDefGradForMeltedParticles;
  int m_resetDefGradForScaledSurfaceParticles;
  real64 m_separabilityMinDamage;
  int m_setDomainTemperature;
  int m_setDomainTemperatureRate;
  int m_shapeFunctionDiagnostics;
  stdVector< array3d< real64 > > m_shapeFunctionGradientValues; // mappedNodes[subregion][particle][nodal shape function gradient
                                                                // value][direction]. dims = {# of subregions, # of
                                                                // particles, # of nodes
                                                                // a particle on the subregion maps to, 3}
  stdVector< array2d< real64 > > m_shapeFunctionValues; // mappedNodes[subregion][particle][nodal shape function value]. dims = {# of
                                                        // subregions, # of particles, # of nodes a particle on the
                                                        // subregion maps to}
  real64 m_smallMass;
  int m_solverProfiling;
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
  int m_subdivideParticles; // Gas particles larger than a grid cell are subdivided
  int m_surfaceDetection;
  int m_surfaceHealing;
  real64 m_surfaceNormalAndPositionDamageThreshold;
  real64 m_surfaceQualityThreshold;  // value [0,1]; same-group DFG surfaces require mapping-normal tensor quality >= this threshold.
  real64 m_surfaceTensionForceSign;
  real64 m_surfaceTensionGradientTolerance;
  array2d< real64 > m_surfaceTensionPairs;
  real64 m_thinFeatureDFGThreshold;
  mpm::TimeIntegrationOption m_timeIntegrationOption;
  int m_tracerCycleInterval;
  array2d< real64 > m_tracerCoordinates;
  int m_tracerHistory;
  string m_tracerOutputPrefix;
  array1d< globalIndex > m_tracerParticleIDs;
  string_array m_tracerVariables;
  real64 m_tracerWriteInterval;
  real64 m_totalBinderVolume;
  int m_treatFullyDamagedAsSingleField;
  mpm::UpdateMethodOption m_updateMethod;
  int m_updateOrder;
  int m_useCrackTipDetection;
  int m_useEvents;                   // Events flag
  int m_useInternalForceAsFaceReaction;
  int m_useNodePosForArea;
  int m_useReferenceVectorsForParticleUpdate;
  int m_useSurfacePositionForContact;
  real64 m_weakInterfaceTraceGapStabilization;
  real64 m_weakInterfaceTraceMinWeight;
  array2d< int > m_weakInterfaceTracePairs;
  int m_weakInterfaceTraceProjectionIterations;
  real64 m_weakInterfaceTraceProjectionScale;
  int m_weakInterfaceTraceSuppressNodalContact;
  int m_writeParticleData;
  array1d< real64 > m_xGlobalMax;         // Maximum global grid coordinate excluding buffer nodes
  array1d< real64 > m_xGlobalMin;         // Minimum global grid coordinate excluding buffer nodes
  array1d< real64 > m_xLocalMax;          // Maximum local grid coordinate including ghost nodes
  array1d< real64 > m_xLocalMaxNoGhost;   // Maximum local grid coordinate EXCLUDING ghost nodes
  array1d< real64 > m_xLocalMin;          // Minimum local grid coordinate including ghost nodes
  array1d< real64 > m_xLocalMinNoGhost;   // Minimum local grid coordinate EXCLUDING ghost nodes

private:
  struct BinKey
  {
    localIndex regionIndex;
    localIndex subRegionIndex;
    localIndex binIndex;

    bool operator==( BinKey const & other ) const
    {
      return (regionIndex == other.regionIndex && subRegionIndex == other.subRegionIndex && binIndex == other.binIndex);
    }
  };

  struct BinKeyHash
  {
    std::size_t operator()( BinKey const & k ) const
    {
      using std::size_t;
      using std::hash;

      // Compute individual hash values for first,
      // second and third and combine them using XOR
      // and bit shifting:
      return ((std::hash< localIndex >()( k.regionIndex )
               ^ (std::hash< localIndex >()( k.subRegionIndex ) << 1)) >> 1)
             ^ (std::hash< localIndex >()( k.binIndex ) << 1);
    }
  };

  void setParticlesConstitutiveNames( ParticleSubRegionBase & subRegion ) const;
};

} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_SOLIDMECHANICS_MPM_HPP_ */
