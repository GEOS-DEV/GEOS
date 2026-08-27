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
 * @file MPMEnums.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_SOLIDMECHANICS_MPMSOLVERENUMS_HPP_
#define GEOS_PHYSICSSOLVERS_SOLIDMECHANICS_MPMSOLVERENUMS_HPP_

#include "common/format/EnumStrings.hpp"

namespace geos
{

namespace mpm
{

    
  /**
   * @enum TimeIntegrationOption
   *
   * The options for time integration
   */
  enum class TimeIntegrationOption : integer
  {
    QuasiStatic,      //!< QuasiStatic
    ImplicitDynamic,  //!< ImplicitDynamic
    ExplicitDynamic   //!< ExplicitDynamic
  };

  /**
   * @enum UpdateMethodOption
   *
   * The options for time integration
   */
  enum class UpdateMethodOption : integer
  {
    FLIP,      //!< FLIP
    PIC,       //!< PIC
    XPIC,      //!< XPIC
    FMPM       //!< FMPM
  };

  /**
   * @enum GridToParticleMappingOption
   *
   * Selects how grid-to-particle kernels obtain their node connectivity,
   * shape-function values, and shape-function gradients.
   */
  enum class GridToParticleMappingOption : integer
  {
    Precomputed, //!< Reuse the raw mapping arrays built once per step.
    OnTheFly     //!< Recompute the raw particle mapping inside each G2P kernel.
  };

  /**
   * @enum BoundaryConditionOption
   *
   * The options for essential boundary conditions
   */
  enum struct BoundaryConditionOption : integer
  {
    Outflow,    //!<Outflow
    Symmetry,   //!<Symmetry
    Moving,     //!<Moving
    Contact     //!<Contact
  };

  /**
   * @enum SurfaceFlag
   *
   * The flags associated with different surface types
   */
  enum struct SurfaceFlag : integer
  {
    Interior = 0,
    FullyDamaged = 1,
    Surface = 2,
    Cohesive = 3,
    DamagedCohesive = 4,
    WeakDiscontinuity = 5
  };

  /**
   * @enum InterpolationOption
   *
   * The options for interpolating tables
   */
  enum struct InterpolationOption : integer
  {
    Linear,
    Cosine,
    Smoothstep
  };

  /**
   * @enum ContactNormalTypeOption
   *
   * The options for contact gap correction
   */
  enum struct ContactNormalTypeOption : integer
  {
    Difference,
    MassWeighted,
    LargerMass,
    Mixed,
    Aligned,
    LogisticRegression
  };

  /**
   * @enum ContactGapCorrectionOption
   *
   * The options for contact gap correction
   */
  enum struct ContactGapCorrectionOption : integer
  {
    Simple,
    Implicit,
    Softened
  };

  /**
   * @enum CPDIDomainScalingTypeOption
   *
   * The algorithms available for limiting a CPDI particle integration domain.
   */
  enum struct CPDIDomainScalingTypeOption : integer
  {
    Homel,     //!< Radially clip selected center-to-corner vectors and reconstruct the r-vectors.
    VPHencky  //!< Preserve as much deviatoric Hencky stretch as possible subject to the domain-diagonal limit.
  };

  /**
   * @enum DomainResetTypeOption
   *
   * The algorithms available when the physical deformation gradient is reset.
   */
  enum struct DomainResetTypeOption : integer
  {
    IsotropicPolar, //!< Preserve rotation and volume while replacing the stretch by an isotropic stretch.
    VPHencky       //!< Preserve volume and the largest feasible fraction of the deviatoric Hencky stretch.
  };

  /**
   * @enum OversizedParticleTreatmentOption
   *
   * The actions available when a particle domain exceeds the oversized-domain
   * diagonal threshold.
   */
  enum struct OversizedParticleTreatmentOption : integer
  {
    None,                     //!< Do not apply an oversized-particle-specific treatment.
    ResetDeformationGradient, //!< Reset the deformation gradient/domain representation.
    Split                     //!< Split the particle unless the rank particle-count limit requests reset fallback.
  };

  /**
   * @enum AreaIntegrationOption
   *
   * The options for nodal area integration
   */
  enum struct AreaIntegrationOption : integer
  {
    BruteForce,
    Mesh
  };

  /**
   * @enum OverlapCorrectionOption
   *
   * The options for overlap correction
   */
  enum struct OverlapCorrectionOption : integer
  {
    Off,
    NormalForce,
    SPH,
    Volume
  };

  /**
   * @enum CohesiveLawOption
   *
   * The options for cohesive laws
   */
  enum struct CohesiveLawOption : integer
  {
    Uncoupled,
    NeedlemanXu,
    Polymer
  };

  /**
   * @enum CohesiveSurfaceDisplacementUpdateOption
   *
   * The options for updating the cohesive particle surface displacement.
   */
  enum struct CohesiveSurfaceDisplacementUpdateOption : integer
  {
    TypeA, //!< Surface displacement from the deformed stored surface-position vector, Fig. 4a.
    TypeB, //!< Surface displacement from the deformed CPDI particle-face vector, Fig. 4b.
    Nodal  //!< Bugfix nodal displacement update from the cohesive-node reference position.
  };

  enum struct GPUSchemeOption : integer
  {
    Atomics,
    NoAtomics,
    RandomMix,
    MinimalAtomics,
    Reduction,
    Colors
  };

  enum struct NormalsAndPositionsMethodOption : integer
  {
    LogisticRegression,
    DFGAndVolumeIntegration
  };

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
              "DamagedCohesive",
              "WeakDiscontinuity" );

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

GEOS_HOST_DEVICE
inline integer toInteger( SurfaceFlag const flag )
{
  return static_cast< integer >( flag );
}

GEOS_HOST_DEVICE
inline SurfaceFlag toSurfaceFlag( integer const flag )
{
  return static_cast< SurfaceFlag >( flag );
}

} // namespace mpm

} // namespace geos

#endif // GEOS_PHYSICSSOLVERS_SOLIDMECHANICS_MPMSOLVERENUMS_HPP_
