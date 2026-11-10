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
 * @file MPMSolverEnums.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_SOLIDMECHANICS_MPMSOLVERENUMS_HPP_
#define GEOS_PHYSICSSOLVERS_SOLIDMECHANICS_MPMSOLVERENUMS_HPP_

#include "common/format/EnumStrings.hpp"

namespace geos
{

namespace mpm
{

/**
 * @brief Selects the solver's time-integration mode.
 *
 * These are input/API choices, not a statement that every mode is implemented.
 * The current SolidMechanicsMPM solver step executes ExplicitDynamic and
 * rejects the other selections.
 */
enum class TimeIntegrationOption : integer
{
  QuasiStatic,      //!< Reserved selection for quasi-static integration.
  ImplicitDynamic,  //!< Reserved selection for implicit dynamics.
  ExplicitDynamic   //!< Explicit grid momentum update and particle advection.
};

/**
 * @brief Selects the grid-to-particle velocity update in explicit MPM.
 *
 * Each method shares the same grid/particle state but differs in how particle
 * velocity is reconstructed from nodal motion and prior particle velocity.
 */
enum class UpdateMethodOption : integer
{
  FLIP, //!< Increment particle velocity by the gathered grid velocity change.
  PIC,  //!< Replace particle velocity by the gathered grid velocity.
  XPIC, //!< Apply recursive particle/grid projections to reduce PIC dissipation.
  FMPM  //!< Apply higher-order filtered/incremental MPM velocity updates.
};

/**
 * @brief Selects where G2P obtains its raw particle interpolation stencil.
 *
 * This concerns grid-to-particle evaluation, not the effective/coalesced P2G
 * mapping arrays used when scattering particle quantities to grid nodes.
 */
enum class GridToParticleMappingOption : integer
{
  Precomputed, //!< Read connectivity, values and gradients assembled earlier in the step.
  OnTheFly     //!< Re-evaluate the particle stencil inside the G2P kernel.
};

/**
 * @brief Classifies the prescribed behavior of a physical grid boundary.
 *
 * These values select boundary handling independently of the material-pair
 * contact law used at internal interfaces.
 */
enum struct BoundaryConditionOption : integer
{
  Outflow,  //!< Permit material to pass through the boundary.
  Symmetry, //!< Enforce the prescribed symmetry component(s) of motion.
  Moving,   //!< Prescribe the boundary velocity/motion.
  Contact   //!< Use the wall-contact treatment for this physical boundary.
};

/**
 * @brief Identifies the particle's free-surface, damage and cohesive status.
 *
 * The explicit integer values are stored in particle surface-flag fields;
 * downstream filtering and normal reconstruction distinguish these states.
 */
enum struct SurfaceFlag : integer
{
  Interior = 0,        //!< Particle is not classified as a physical surface.
  FullyDamaged = 1,    //!< Damaged surface state used by DFG reconstruction.
  Surface = 2,         //!< Explicit/free surface particle.
  Cohesive = 3,        //!< Surface carrying a cohesive interface definition.
  DamagedCohesive = 4 //!< Cohesive surface with a damaged interface state.
};

/**
 * @brief Selects the interpolation function between successive table entries.
 *
 * Value and rate use the same interpolation law within each table interval.
 * Cosine has a trigonometric ease-in/ease-out profile; Smoothstep uses a
 * fifth-order polynomial with zero first and second endpoint derivatives.
 */
enum struct InterpolationOption : integer
{
  Linear,    //!< Piecewise linear interpolation with constant interval rate.
  Cosine,    //!< Cosine easing with zero endpoint rate.
  Smoothstep //!< Quintic smoothstep easing with smooth endpoint derivatives.
};

/**
 * @brief Chooses how opposing grid-field normals form one contact normal.
 *
 * The resulting pair normal points from field A toward B. The selected
 * arbitration precedes dimensional projection and normalization.
 */
enum struct ContactNormalTypeOption : integer
{
  Difference,         //!< Difference between the two field normals.
  MassWeighted,       //!< Difference weighted by the respective nodal masses.
  LargerMass,         //!< Use the normal of the field with larger nodal mass.
  Mixed,              //!< Mass-weight if densities agree; otherwise use the denser field.
  Aligned,            //!< Weight field normals by mapped-normal alignment quality.
  LogisticRegression  //!< Fit a separating normal from nearby particle samples.
};

/**
 * @brief Chooses how the contact gap gates the approach correction.
 *
 * All three modes require the pair to be approaching. They differ in whether
 * existing overlap is required or a positive-gap transition is permitted.
 */
enum struct ContactGapCorrectionOption : integer
{
  Simple,   //!< Apply on approach without requiring a negative surface gap.
  Implicit, //!< Apply only on approach with a negative surface gap.
  Softened  //!< Taper the correction across a finite positive-gap interval.
};

/**
 * @brief Chooses how an oversized CPDI box integration domain is limited.
 *
 * These act on the mapping/integration representation, not the physical
 * constitutive stress update. The simplex/CPTI projection is handled
 * separately from box-specific scaling.
 */
enum struct CPDIDomainScalingTypeOption : integer
{
  Homel,    //!< Limit center-to-corner reach, then reconstruct box generators.
  VPHencky  //!< Retain feasible deviatoric Hencky stretch at fixed volume.
};

/**
 * @brief Selects which reference quantities are rebased by a geometry reset.
 */
enum struct GeometryResetOption : integer
{
  FOnly,       //!< Reset the deformation-gradient reference without material-geometry rebasing.
  FAndGeometry //!< Also rebase the stored material-geometry reference.
};

/**
 * @brief Chooses the replacement stretch during a deformation/domain reset.
 *
 * Both methods preserve volume; the selected rule controls the remaining
 * deviatoric stretch while preserving the appropriate rotation.
 */
enum struct DomainResetTypeOption : integer
{
  IsotropicPolar, //!< Preserve polar rotation and volume; use isotropic stretch.
  VPHencky       //!< Preserve the largest feasible fraction of deviatoric Hencky stretch.
};

/**
 * @brief Chooses the response when particle size exceeds the grid-scale limit.
 *
 * Split can fall back to resetting the domain when the per-rank particle
 * count is above the configured subdivision limit.
 */
enum struct OversizedParticleTreatmentOption : integer
{
  None,                     //!< Leave the oversized domain untreated by this policy.
  ResetDeformationGradient, //!< Rebase the deformation-gradient/domain representation.
  Split                     //!< Subdivide when allowed; otherwise use the reset fallback.
};

/**
 * @brief Selects the geometric integration used to estimate nodal surface area.
 */
enum struct AreaIntegrationOption : integer
{
  BruteForce, //!< Integrate sampled surface points around a node.
  Mesh        //!< Integrate reconstructed surface patches using neighboring particles.
};

/**
 * @brief Selects an additional overlap correction for contact/particle motion.
 *
 * These choices supplement the pairwise contact gap condition; they do not
 * change the definition of the underlying particle material geometry.
 */
enum struct OverlapCorrectionOption : integer
{
  Off,         //!< No additional overlap correction.
  NormalForce, //!< Add an overlap-dependent normal contact force.
  SPH,         //!< Use the SPH-based particle overlap correction.
  Volume       //!< Use the volume-based overlap correction.
};

/**
 * @brief Selects the traction-separation constitutive law at a cohesive zone.
 */
enum struct CohesiveLawOption : integer
{
  Uncoupled,  //!< Independent normal and tangential cohesive response.
  NeedlemanXu, //!< Coupled Needleman-Xu traction-separation model.
  Polymer    //!< Polymer-interface cohesive response.
};

/**
 * @brief Chooses how cohesive-interface surface displacement is reconstructed.
 *
 * TypeA and TypeB use a surface offset transported from the cohesive
 * reference configuration; Nodal uses cohesive-node displacement relative
 * to its stored reference position. These are current algorithms, not patch
 * versions.
 */
enum struct CohesiveSurfaceDisplacementUpdateOption : integer
{
  TypeA, //!< Transport the stored reference surface-position vector.
  TypeB, //!< Transport a reference offset derived from the CPDI domain face.
  Nodal  //!< Reconstruct displacement from cohesive-node reference motion.
};

/**
 * @brief Selects a particle-to-grid execution strategy.
 *
 * The current solver dispatch supports Atomics. The other enum labels remain
 * registered but are not executable P2G choices in this implementation.
 */
enum struct GPUSchemeOption : integer
{
  Atomics,        //!< Use atomic accumulation during particle-to-grid scatter.
  NoAtomics,      //!< Registered alternative; not supported by current dispatch.
  RandomMix,      //!< Registered alternative; not supported by current dispatch.
  MinimalAtomics, //!< Registered alternative; not supported by current dispatch.
  Reduction,      //!< Registered alternative; not supported by current dispatch.
  Colors          //!< Registered alternative; not supported by current dispatch.
};

/**
 * @brief Selects reconstruction of nodal surface normals and positions.
 */
enum struct NormalsAndPositionsMethodOption : integer
{
  LogisticRegression,     //!< Fit a local separating surface to neighboring particles.
  DFGAndVolumeIntegration //!< Use mapped DFG normals and geometric support-volume inversion.
};

/**
 * @brief Register the exact input/output spellings of each enum.
 *
 * The strings below are externally visible option values. Their spelling and
 * order (including lowercase entries such as "vpHencky") must remain paired
 * with the enum members above; formatting changes must not alter them.
 */

ENUM_STRINGS( mpm::TimeIntegrationOption,
              "QuasiStatic",
              "ImplicitDynamic",
              "ExplicitDynamic" );

ENUM_STRINGS( mpm::UpdateMethodOption,
              "FLIP",
              "PIC",
              "XPIC",
              "FMPM" );

ENUM_STRINGS( mpm::GridToParticleMappingOption,
              "Precomputed",
              "OnTheFly" );

ENUM_STRINGS( mpm::BoundaryConditionOption,
              "Outflow",
              "Symmetry",
              "Moving",
              "Contact" );

ENUM_STRINGS( mpm::SurfaceFlag,
              "Interior",
              "FullyDamaged",
              "Surface",
              "Cohesive",
              "DamagedCohesive" );

ENUM_STRINGS( mpm::InterpolationOption,
              "Linear",
              "Cosine",
              "Smoothstep" );

ENUM_STRINGS( mpm::ContactNormalTypeOption,
              "Difference",
              "MassWeighted",
              "LargerMass",
              "Mixed",
              "Aligned",
              "LogisticRegression" );

ENUM_STRINGS( mpm::ContactGapCorrectionOption,
              "Simple",
              "Implicit",
              "Softened" );

ENUM_STRINGS( mpm::CPDIDomainScalingTypeOption,
              "homel",
              "vpHencky" );

ENUM_STRINGS( mpm::GeometryResetOption, "FOnly", "FAndGeometry" );

ENUM_STRINGS( mpm::DomainResetTypeOption,
              "isotropicPolar",
              "vpHencky" );

ENUM_STRINGS( mpm::OversizedParticleTreatmentOption,
              "none",
              "resetDeformationGradient",
              "split" );

ENUM_STRINGS( mpm::AreaIntegrationOption,
              "BruteForce",
              "Mesh" );

ENUM_STRINGS( mpm::OverlapCorrectionOption,
              "Off",
              "NormalForce",
              "SPH",
              "Volume" );

ENUM_STRINGS( mpm::CohesiveLawOption,
              "Uncoupled",
              "NeedlemanXu",
              "Polymer" );

ENUM_STRINGS( mpm::CohesiveSurfaceDisplacementUpdateOption,
              "TypeA",
              "TypeB",
              "Nodal" );

ENUM_STRINGS( mpm::GPUSchemeOption,
              "Atomics",
              "NoAtomics",
              "RandomMix",
              "MinimalAtomics",
              "Reduction",
              "Colors" );

ENUM_STRINGS( mpm::NormalsAndPositionsMethodOption,
              "LogisticRegression",
              "DFGAndVolumeIntegration" );

/** @brief Return the stored integer representation of a particle surface flag. */
GEOS_HOST_DEVICE
inline integer toInteger( SurfaceFlag const flag )
{
  return static_cast< integer >( flag );
}

/**
 * @brief Convert a stored integer to SurfaceFlag without range validation.
 *
 * Callers must ensure the value denotes one of the enumerators above.
 */
GEOS_HOST_DEVICE
inline SurfaceFlag toSurfaceFlag( integer const flag )
{
  return static_cast< SurfaceFlag >( flag );
}

} // namespace mpm

} // namespace geos

#endif // GEOS_PHYSICSSOLVERS_SOLIDMECHANICS_MPMSOLVERENUMS_HPP_
