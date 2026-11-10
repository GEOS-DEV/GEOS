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
 * @file MPMSolverBaseFields.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_SOLIDMECHANICS_MPMSOLVERBASEFIELDS_HPP_
#define GEOS_PHYSICSSOLVERS_SOLIDMECHANICS_MPMSOLVERBASEFIELDS_HPP_

#include "mesh/MeshFields.hpp"
#include "physicsSolvers/solidMechanics/MPMSolverEnums.hpp"

namespace geos
{

/**
 * @brief Registry of mesh field traits used by solver components.
 *
 * Each DECLARE_FIELD fixes the repository key, value type, default, plot level,
 * restart policy, and registered description string. The registered metadata
 * controls field allocation, restart handling and default visualization.
 */
namespace fields
{
namespace mpm
{

// These are field *declarations*. SolidMechanicsMPM::registerDataOnMesh and
// material/particle managers determine which traits are registered and sized
// for an individual particle subregion.

  // Core particle identity, mass and deletion.
  // The first dimension of every trait is the particle-storage index. Scalar fields need no extra
  // dimensions;
  // dimensions of vector, tensor and stencil fields are set at registration.

  /**
   * @brief Current mass assigned to each particle [particle].
   * The mass is used in momentum/force transfers and apportioned among children when a particle is
   * subdivided.
   */
  DECLARE_FIELD( particleMass,
                 "particleMass",
                 array1d< real64 >,
                 0.0,
                 LEVEL_1,
                 WRITE_AND_READ,
                 "An array that stores particle masses." );

  /**
   * @brief Deferred-deletion marker [particle], not a damage variable.
   * The solver can flag a particle during a step and remove its storage in the later cleanup pass.
   */
  DECLARE_FIELD( particleDeleteFlag,
                 "particleDeleteFlag",
                 array1d< int >,
                 0,
                 NOPLOT,
                 WRITE_AND_READ,
                 "An array that remembers particles that should be deleted at the end of the time step." );

  /**
   * @brief Constitutive-model update status [particle][Gauss point].
   * The MPM registration currently allocates one Gauss-point slot per particle; a negative value in slot
   * zero requests particle deletion.
   */
  DECLARE_FIELD( particleConstitutiveUpdateFlag,
                 "constitutiveUpdateFlag",
                 array2d< int >,
                 0,
                 LEVEL_1,
                 WRITE_AND_READ,
                 "An optional per-gauss-point constitutive-model update flag copied back to the solver. "
                 "Negative values at gauss point 0 flag particles for deletion." );

  /**
   * @brief Per-particle state marker for crystal-healing events.
   * Event logic uses this flag independently of the general damage and cohesive-surface classifications.
   */
  DECLARE_FIELD( particleCrystalHealFlag,
                 "particleCrystalHealFlag",
                 array1d< int >,
                 0,
                 LEVEL_1,
                 WRITE_AND_READ,
                 "An array that remembers particles that are undergoing crystal healing for the mpm event." );

  /**
   * @brief Integer material-type index [particle].
   * The solver assumes a single material assignment within each particle region.
   */
  DECLARE_FIELD( particleMaterialType,
                 "particleMaterialType",
                 array1d< int >,
                 0,
                 LEVEL_1,
                 WRITE_AND_READ,
                 "An array that stores index of particle material type "
                 "(assumes single material for each particle region)." );

  // Thermal and energy-related state.
  // Current temperature, its rate, and thermodynamic reference quantities are distinct.
  // Fields used only by optional models may not be registered on every subregion.

  /**
   * @brief Particle acoustic wave speed, used when estimating a stable explicit time step.
   * This is a material-provided speed rather than the particle's transport velocity.
   */
  DECLARE_FIELD( particleWavespeed,
                 "particleWavespeed",
                 array1d< real64 >,
                 0.0,
                 LEVEL_1,
                 WRITE_AND_READ,
                 "An array that stores particle wavespeeds." );

  /**
   * @brief Optional per-particle specific heat / heat-capacity state.
   * Registered when thermal energy updates are enabled; material initialization can supply a tangent
   * specific heat.
   */
  DECLARE_FIELD( particleHeatCapacity,
                 "particleHeatCapacity",
                 array1d< real64 >,
                 0.0,
                 LEVEL_1,
                 WRITE_AND_READ,
                 "An array that stores particle heat capacity." );

  /**
   * @brief Reference temperature associated with the particle's thermal state.
   * This is distinct from the evolving particle temperature and can be reinitialized from constitutive
   * state.
   */
  DECLARE_FIELD( particleReferenceTemperature,
                 "particleReferenceTemperature",
                 array1d< real64 >,
                 0.0,
                 LEVEL_1,
                 WRITE_AND_READ,
                 "An array that stores particle reference temperature." );

  /**
   * @brief Current particle temperature, distinct from particleReferenceTemperature.
   * It evolves during thermal updates and is not the temperature at the constitutive reference state.
   */
  DECLARE_FIELD( particleTemperature,
                 "particleTemperature",
                 array1d< real64 >,
                 0.0,
                 LEVEL_1,
                 WRITE_AND_READ,
                 "An array that stores particle reference temperature." );

  /**
   * @brief Time rate of the evolving particle temperature.
   * The solver may prescribe this rate from a domain thermal loading or initialize it to zero.
   */
  DECLARE_FIELD( particleTemperatureRate,
                 "particleTemperatureRate",
                 array1d< real64 >,
                 0.0,
                 LEVEL_1,
                 WRITE_AND_READ,
                 "An array that stores particle reference temperatureRate." );

  /**
   * @brief Particle specific internal energy (energy per unit mass).
   * The initial thermodynamic state obtains this value from the constitutive model; it is not kinetic
   * energy.
   */
  DECLARE_FIELD( particleInternalEnergy,
                 "particleInternalEnergy",
                 array1d< real64 >,
                 0.0,
                 LEVEL_1,
                 WRITE_AND_READ,
                 "An array that stores particle internal energy." );

  /**
   * @brief Particle translational kinetic energy, accumulated for energy diagnostics.
   * It represents energy of an individual particle rather than specific energy per unit mass.
   */
  DECLARE_FIELD( particleKineticEnergy,
                 "particleKineticEnergy",
                 array1d< real64 >,
                 0.0,
                 LEVEL_1,
                 WRITE_AND_READ,
                 "An array that stores particle kinetic energy." );

  /**
   * @brief Additional scalar pressure with positive-compression sign convention.
   * The MPM force assembly applies this contribution separately from the six stored particleStress
   * components.
   */
  DECLARE_FIELD( particleSupplementalPressure,
                 "particleSupplementalPressure",
                 array1d< real64 >,
                 0.0,
                 LEVEL_1,
                 WRITE_AND_READ,
                 "Positive-compressive scalar pressure supplement used by MPM force assembly; not included in particleStress." );

  /**
   * @brief Constitutive phase/melt indicator stored as a real-valued particle scalar.
   * The surface-data invalidation path treats values greater than 0.5 as melted.
   */
  DECLARE_FIELD( particleMeltFlag,
                 "meltFlag",
                 array1d< real64 >,
                 0.0,
                 LEVEL_1,
                 WRITE_AND_READ,
                 "Melt diagnostic copied from constitutive models that expose a meltFlag wrapper." );

  /**
   * @brief Indicator that the deformation-gradient / particle-domain representation was reset.
   * It identifies a reset performed by the relevant reset pass; it is not the accumulated strain or
   * damage.
   */
  DECLARE_FIELD( particleDomainResetFlag,
                 "particleDomainResetFlag",
                 array1d< int >,
                 0,
                 LEVEL_1,
                 WRITE_AND_READ,
                 "Flag indicating that the particle deformation-gradient/domain representation was reset in the most recent reset pass." );

  // Reference volumes and geometric/kinematic state.
  // The stored deformation tensors retain three physical coordinates in native 2D and axisymmetric
  // formulations.
  // Reference geometry and the separately stored current particle domain have different lifetimes.

  /**
   * @brief Optional neighbor-kernel estimate of the particle Jacobian (volume ratio).
   * Allocated only when the SPH-Jacobian correction is enabled, independently of the stored deformation-
   * gradient determinant.
   */
  DECLARE_FIELD( particleSPHJacobian,
                 "particleSPHJacobian",
                 array1d< real64 >,
                 0.0,
                 LEVEL_1,
                 WRITE_AND_READ,
                 "An array that stores particle SPH computed jacobian." );

  /**
   * @brief Per-particle distance estimate for the optional crack-tip detection path.
   * Its numerical construction belongs to the crack-tip routine; this trait stores the resulting scalar.
   */
  DECLARE_FIELD( particleCrackTipDistance,
                 "particleCrackTipDistance",
                 array1d< real64 >,
                 0.0,
                 LEVEL_1,
                 WRITE_AND_READ,
                 "An array that stores particle Crack Tip Distance." );

  /**
   * @brief Particle material volume in the chosen reference configuration.
   * The current physical volume is carried separately; reference volume is used in material-volume and
   * deformation accounting.
   */
  DECLARE_FIELD( particleReferenceVolume,
                 "particleReferenceVolume",
                 array1d< real64 >,
                 0.0,
                 NOPLOT,
                 WRITE_AND_READ,
                 "ParticleReferenceVolume" );

  /**
   * @brief Particle porosity captured in the reference state.
   * Initialization copies the constitutive particle porosity here; subsequent current porosity is a
   * separate material state.
   */
  DECLARE_FIELD( particleReferencePorosity,
                 "particleReferencePorosity",
                 array1d< real64 >,
                 0.0,
                 NOPLOT,
                 WRITE_AND_READ,
                 "ParticleReferencePorosity" );

  /**
   * @brief Reference affine particle-domain generators, stored as [particle][DIM][DIM].
   * These are full-edge geometry vectors, not the particle-center coordinates or the temporary projected
   * integration geometry.
   */
  DECLARE_FIELD( particleGeometryReference,
                 "particleGeometryReference",
                 array3d< real64 >,
                 0.0,
                 NOPLOT,
                 WRITE_AND_READ,
                 "Reference full-edge particle geometry (DIM by DIM)" );

  /**
   * @brief Current physical deformation gradient F, stored as [particle][3][3].
   * The full tensor is retained even for native 2D; axisymmetric updates use the third diagonal for hoop
   * stretch.
   */
  DECLARE_FIELD( particleDeformationGradient,
                 "particleDeformationGradient",
                 array3d< real64 >,
                 0.0,
                 LEVEL_1,
                 WRITE_AND_READ,
                 "ParticleDeformationGradient" );

  /**
   * @brief Current particle rotation tensor [particle][3][3].
   * Stored separately from F for constitutive kinematics that use a rotational frame.
   */
  DECLARE_FIELD( particleRotation,
                 "particleRotation",
                 array3d< real64 >,
                 0.0,
                 NOPLOT,
                 WRITE_AND_READ,
                 "ParticleRotation" );

  /**
   * @brief Previous particle rotation tensor [particle][3][3].
   * Keeping both rotation states supports incremental transformations between time levels.
   */
  DECLARE_FIELD( particleOldRotation,
                 "particleOldRotation",
                 array3d< real64 >,
                 0.0,
                 NOPLOT,
                 WRITE_AND_READ,
                 "ParticleOldRotation" );

  /**
   * @brief Step-secant rate of deformation gradient, (F_next - F_previous) / dt.
   * The rate is retained for constitutive-model compatibility and is distinct from the spatial velocity
   * gradient.
   */
  DECLARE_FIELD( particleFDot,
                 "particleFDot",
                 array3d< real64 >,
                 0.0,
                 NOPLOT,
                 WRITE_AND_READ,
                 "Step-secant deformation-gradient rate (F_{n+1}-F_n)/dt retained for compatibility." );

  /**
   * @brief Particle spatial velocity gradient L, stored as [particle][3][3].
   * The inactive Cartesian 2D components are handled by the formulation; axisymmetry can supply a hoop
   * component.
   */
  DECLARE_FIELD( particleVelocityGradient,
                 "particleVelocityGradient",
                 array3d< real64 >,
                 0.0,
                 NOPLOT,
                 WRITE_AND_READ,
                 "ParticleVelocityGradient" );

  // Stress, body loading and damage.
  // Vector arrays are [particle][3]; stress and plastic strain use the six-component Voigt convention.
  // Constitutive stress excludes the independently stored supplemental pressure.

  /**
   * @brief Particle stress in six-component Voigt storage [particle][6].
   * The registration labels are XX, YY, ZZ, YZ, XZ, XY; supplemental scalar pressure is stored
   * separately.
   */
  DECLARE_FIELD( particleStress,
                 "particleStress",
                 array2d< real64 >,
                 0.0,
                 LEVEL_1,
                 WRITE_AND_READ,
                 "An array that holds particle stresses in Voigt notation." );

  /**
   * @brief Prescribed particle body acceleration [particle][3].
   * The P2G force assembly multiplies this acceleration by particle mass before adding the applied
   * force.
   */
  DECLARE_FIELD( particleBodyForce,
                 "particleBodyForce",
                 array2d< real64 >,
                 0.0,
                 LEVEL_1,
                 WRITE_AND_READ,
                 "An array that holds particle body forces." );

  /**
   * @brief Particle plastic strain in six-component Voigt storage [particle][6].
   * Uses the same XX, YY, ZZ, YZ, XZ, XY component order as the registered stress array.
   */
  DECLARE_FIELD( particlePlasticStrain,
                 "particlePlasticStrain",
                 array2d< real64 >,
                 0.0,
                 LEVEL_1,
                 WRITE_AND_READ,
                 "An array that holds particle plastic strain in Voigt notation." );

  /**
   * @brief Current particle material density, separate from its mass and volume.
   * Thermodynamic initialization can reset density and mass consistently with the current physical
   * volume.
   */
  DECLARE_FIELD( particleDensity,
                 "particleDensity",
                 array1d< real64 >,
                 0.0,
                 LEVEL_1,
                 WRITE_AND_READ,
                 "An array that holds particle densities." );

  /**
   * @brief Spatial gradient of the particle damage indicator [particle][3].
   * The solver can reconstruct it from neighboring particles for DFG and surface-orientation decisions;
   * the axisymmetric path uses meridional weighting.
   */
  DECLARE_FIELD( particleDamageGradient,
                 "particleDamageGradient",
                 array2d< real64 >,
                 0.0,
                 LEVEL_1,
                 WRITE_AND_READ,
                 "An array that holds particle damage gradients as calculated with an SPH kernel." );

  /**
   * @brief Integer-coded particle surface classification (mpm::SurfaceFlag).
   * Distinguishes interior, damaged, ordinary surface, and cohesive surface states rather than recording
   * a scalar surface area.
   */
  DECLARE_FIELD( particleSurfaceFlag,
                 "particleSurfaceFlag",
                 array1d< integer >,
                 0,
                 LEVEL_1,
                 WRITE_AND_READ,
                 "An array that holds particle surface flags." );

  /**
   * @brief Optional SPH-based deformation-gradient estimate [particle][3][3].
   * The solver reconstructs a neighbor displacement gradient and converts it into a deformation-like
   * tensor for overlap-related diagnostics.
   */
  DECLARE_FIELD( particleSPHF,
                 "particleSPHF",
                 array3d< real64 >,
                 0.0,
                 LEVEL_1,
                 WRITE_AND_READ,
                 "ParticleSPHF" );

  /**
   * @brief Optional per-particle overlap diagnostic associated with directional overlap correction.
   * Registered with the SPH overlap fields when that correction is enabled.
   */
  DECLARE_FIELD( particleOverlap,
                 "particleOverlap",
                 array1d< real64 >,
                 0.0,
                 LEVEL_1,
                 WRITE_AND_READ,
                 "ParticleOverlap" );

  // Particle reference surfaces and material orientation.
  // Reference quantities support surface transport and anisotropic constitutive behavior.
  // Cohesive activation introduces a further reference configuration described below.

  /**
   * @brief Particle center in its reference configuration [particle][3].
   * Retained alongside the evolving center for relative-displacement and neighbor-based deformation
   * calculations.
   */
  DECLARE_FIELD( particleReferencePosition,
                 "particleReferencePosition",
                 array2d< real64 >,
                 0.0,
                 LEVEL_1,
                 WRITE_AND_READ,
                 "ParticleReferencePosition" );

  /**
   * @brief Reference material-orientation basis [particle][3][3].
   * Transport routines combine the reference directions with F or its cofactor for anisotropic material
   * orientations.
   */
  DECLARE_FIELD( particleReferenceMaterialDirection,
                 "particleReferenceMaterialDirection",
                 array3d< real64 >,
                 0.0,
                 LEVEL_1,
                 WRITE_AND_READ,
                 "ParticleReferenceMaterialDirection" );

  /**
   * @brief Reference particle-surface normal [particle][3].
   * Separate from the current normal so surface orientation can be transported or invalidated without
   * discarding the initial direction.
   */
  DECLARE_FIELD( particleReferenceSurfaceNormal,
                 "particleReferenceSurfaceNormal",
                 array2d< real64 >,
                 0.0,
                 LEVEL_1,
                 WRITE_AND_READ,
                 "ParticleReferenceSurfaceNormal" );

  /**
   * @brief Reference surface offset vector relative to the particle center [particle][3].
   * The current surface-position estimate is stored separately and can evolve with deformation.
   */
  DECLARE_FIELD( particleReferenceSurfacePosition,
                 "particleReferenceSurfacePosition",
                 array2d< real64 >,
                 0.0,
                 LEVEL_1,
                 WRITE_AND_READ,
                 "ParticleReferenceSurfacePosition" );

  /**
   * @brief Traction vector associated with the particle's reference surface [particle][3].
   * This reference value is maintained separately from current surface traction and may be zeroed by
   * damage/melt controls.
   */
  DECLARE_FIELD( particleReferenceSurfaceTraction,
                 "particleReferenceSurfaceTraction",
                 array2d< real64 >,
                 0.0,
                 LEVEL_1,
                 WRITE_AND_READ,
                 "ParticleReferenceSurfaceTraction" );

  // Cohesive-zone force and persistent interface mapping.
  // Cohesive references are snapshots taken when an interface is initialized; they are not the
  // ordinary P2G/G2P mapping arrays recomputed each explicit step.

  /**
   * @brief Particle force accumulated from cohesive-interface interactions [particle][3].
   * It is a force contribution, distinct from particleBodyForce, which stores acceleration.
   */
  DECLARE_FIELD( particleCohesiveForce,
                 "particleCohesiveForce",
                 array2d< real64 >,
                 0.0,
                 LEVEL_1,
                 WRITE_AND_READ,
                 "ParticleCohesiveForce" );

  /**
   * @brief Indicator that a particle participates in an initialized cohesive interface.
   * This is separate from its SurfaceFlag classification and is assigned when reference cohesive-node
   * support is found.
   */
  DECLARE_FIELD( particleCohesiveZoneFlag,
                 "particleCohesiveZoneFlag",
                 array1d< int >,
                 0,
                 LEVEL_1,
                 WRITE_AND_READ,
                 "ParticleCohesiveZoneFlag" );

  /**
   * @brief Persistent reference grid-node identifiers for cohesive particles [particle][stencil slot].
   * Stores global grid IDs captured at interface initialization, not the per-step local mapped-node
   * indices.
   */
  DECLARE_FIELD( particleReferenceMappedNodes,
                 "particleReferenceMappedNodes",
                 array2d< globalIndex >,
                 0.0,
                 LEVEL_1,
                 WRITE_AND_READ,
                 "ParticleReferenceMappedNodes" );

  /**
   * @brief Cohesive-node ordinal associated with each saved reference stencil slot.
   * Despite the name and globalIndex type, the stored valid value indexes the cohesive-node list;
   * max(globalIndex) marks non-cohesive nodes.
   */
  DECLARE_FIELD( particleReferenceGlobalNodeIndex,
                 "particleReferenceGlobalNodeIndex",
                 array2d< globalIndex >,
                 0.0,
                 LEVEL_1,
                 WRITE_AND_READ,
                 "Cached cohesive-zone local node index for each particle reference mapped node. "
                 "Non-cohesive mapped nodes are marked with std::numeric_limits< globalIndex >::max()." );

  /**
   * @brief Shape-function weights frozen at cohesive reference initialization [particle][stencil slot].
   * Cohesive updates reuse them after particles migrate or remap to different ordinary grid nodes.
   */
  DECLARE_FIELD( particleReferenceShapeFunctionValues,
                 "particleReferenceShapeFunctionValues",
                 array2d< real64 >,
                 0.0,
                 LEVEL_1,
                 WRITE_AND_READ,
                 "ParticleReferenceShapeFunctionValues" );

  /**
   * @brief Storage for reference stencil shape-function gradients [particle][stencil slot][3].
   * The field is registered with a capacity matching the particle type; the ordinary per-step mapping
   * gradients live in separate solver arrays.
   */
  DECLARE_FIELD( particleReferenceShapeFunctionGradientValues,
                 "particleReferenceShapeFunctionGradientValues",
                 array3d< real64 >,
                 0.0,
                 LEVEL_1,
                 WRITE_AND_READ,
                 "ParticleReferenceShapeFunctionGradientValues" );

  /**
   * @brief Particle-center position at cohesive-interface activation [particle][3].
   * The cohesive relative-deformation calculation measures motion from this activation state, not
   * necessarily the original material reference.
   */
  DECLARE_FIELD( particleCohesiveReferencePosition,
                 "particleCohesiveReferencePosition",
                 array2d< real64 >,
                 0.0,
                 LEVEL_1,
                 WRITE_AND_READ,
                 "particleCohesiveReferencePosition" );

  /**
   * @brief Particle-surface normal frozen at cohesive-interface activation [particle][3].
   * Used with the saved cohesive reference deformation to reconstruct the evolving interface
   * orientation.
   */
  DECLARE_FIELD( particleCohesiveReferenceSurfaceNormal,
                 "particleCohesiveReferenceSurfaceNormal",
                 array2d< real64 >,
                 0.0,
                 LEVEL_1,
                 WRITE_AND_READ,
                 "particleCohesiveReferenceSurfaceNormal" );

  /**
   * @brief Cohesive-interface reference surface offset [particle][3].
   * Its initialization depends on the TypeA, TypeB, or Nodal surface-displacement method.
   */
  DECLARE_FIELD( particleCohesiveReferenceSurfacePosition,
                 "particleCohesiveReferenceSurfacePosition",
                 array2d< real64 >,
                 0.0,
                 LEVEL_1,
                 WRITE_AND_READ,
                 "particleCohesiveReferenceSurfacePosition" );

  /**
   * @brief Physical F at cohesive-interface activation [particle][3][3].
   * Later cohesive kinematics use relative deformation F * inverse(F_at_activation).
   */
  DECLARE_FIELD( particleCohesiveReferenceDeformationGradient,
                 "particleCohesiveReferenceDeformationGradient",
                 array3d< real64 >,
                 0.0,
                 LEVEL_1,
                 WRITE_AND_READ,
                 "particleCohesiveReferenceDeformationGradient" );

  /**
   * @brief Nodal velocity/contact field index for each particle reference-mapping slot.
   * Saved for cohesive force assembly so the interface side can remain consistent with the original
   * field partition.
   */
  DECLARE_FIELD( particleCohesiveFieldMapping,
                 "particleCohesiveFieldMapping",
                 array2d< int >,
                 0,
                 LEVEL_1,
                 WRITE_AND_READ,
                 "particleCohesiveFieldMapping" );

  // Subdivision and domain-scaling work state.
  // These integer arrays are allocated only when their associated optional workflow is active.
  // Unlike constitutive history they can be reset or repopulated during a time step.

  /**
   * @brief Temporary subdivision work field: number of additional children selected for a parent.
   * The flag is cleared when parent and child data have been transferred.
   */
  DECLARE_FIELD( particleSubdivideFlag,
                 "particleSubdivideFlag",
                 array1d< int >,
                 0,
                 LEVEL_1,
                 WRITE_AND_READ,
                 "Temporary number of additional children selected for subdivision; cleared after transfer." );

  /**
   * @brief Temporary child-to-parent particle-storage index for subdivision.
   * A nonnegative value identifies the source parent for state copying; -1 marks no pending copy.
   */
  DECLARE_FIELD( particleCopyFlag,
                 "particleCopyFlag",
                 array1d< int >,
                 -1,
                 LEVEL_1,
                 WRITE_AND_READ,
                 "particleCopyFlag" );

  /**
   * @brief Marker that CPDI domain scaling altered the particle's integration domain.
   * Registered only for enabled domain scaling; downstream routines can use it to adjust surface or
   * reset treatment.
   */
  DECLARE_FIELD( particleDomainScaledFlag,
                 "particleDomainScaledFlag",
                 array1d< int >,
                 -1,
                 LEVEL_1,
                 WRITE_AND_READ,
                 "particleDomainScaledFlag" );

  // Labels and reconstructed surface diagnostics.
  // Output levels control default visualization but do not by themselves define field
  // allocation or restart behavior.

  /**
   * @brief Nonnegative integer color/label associated with each particle.
   * Kept as a separate field from material type, contact-group indexing, and damage/surface
   * classification.
   */
  DECLARE_FIELD( particleColor,
                 "particleColor",
                 array1d< int >,
                 0,
                 LEVEL_0,
                 WRITE_AND_READ,
                 "particleColor" );

  /**
   * @brief Scalar marker for residue/binder identification.
   * Its interpretation is supplied by the component that populates the field; the trait itself fixes
   * only the storage contract.
   */
  DECLARE_FIELD( particleResidueBinderFlag,
                 "particleResidueBinderFlag",
                 array1d< real64 >,
                 0,
                 LEVEL_0,
                 WRITE_AND_READ,
                 "particleResidueBinderFlag" );

  /**
   * @brief Estimated particle-surface offset vector [particle][3].
   * The damage-gradient surface estimator and subsequent subdivision/rebasing logic can update this
   * quantity.
   */
  DECLARE_FIELD( particleEstimatedSurfacePosition,
                 "particleEstimatedSurfacePosition",
                 array2d< real64 >,
                 0,
                 LEVEL_0,
                 WRITE_AND_READ,
                 "particleEstimatedSurfacePosition" );

  /**
   * @brief Second spatial derivative tensor of the reconstructed damage field [particle][3][3].
   * Available for curvature/surface characterization; independent of the first damage gradient.
   */
  DECLARE_FIELD( particleDamageHessian,
                 "particleDamageHessian",
                 array3d< real64 >,
                 0,
                 LEVEL_0,
                 WRITE_AND_READ,
                 "particleDamageHessian" );

  /**
   * @brief Scalar L2 (Frobenius-type) magnitude diagnostic for the particle damage Hessian.
   * Stored separately to support threshold decisions without repeatedly inspecting its tensor
   * components.
   */
  DECLARE_FIELD( particleDamageHessianL2Norm,
                 "particleDamageHessianL2Norm",
                 array1d< real64 >,
                 0,
                 LEVEL_0,
                 WRITE_AND_READ,
                 "particleDamageHessianL2Norm" );

} // namespace mpm
} // namespace fields
} // namespace geos

#endif // GEOS_PHYSICSSOLVERS_SOLIDMECHANICS_MPMSOLVERBASEFIELDS_HPP_
