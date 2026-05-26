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
    Interior,
    FullyDamaged,
    Surface,
    Cohesive,
    DamagedCohesive
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

ENUM_STRINGS( mpm::BoundaryConditionOption,
              "Outflow",
              "Symmetry",
              "Moving",
              "Contact" );

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

} // namespace mpm

} // namespace geos

#endif // GEOS_PHYSICSSOLVERS_SOLIDMECHANICS_MPMSOLVERENUMS_HPP_
