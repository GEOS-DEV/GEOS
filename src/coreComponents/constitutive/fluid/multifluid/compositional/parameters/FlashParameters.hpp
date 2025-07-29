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
 * @file FlashParameters.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_PARAMETERS_FLASHPARAMETERS_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_PARAMETERS_FLASHPARAMETERS_HPP_

#include "ModelParameters.hpp"
#include "common/DataTypes.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

class FlashParameters : public ModelParameters
{
public:
  /**
   * Flash variables
   * These are indices into an array that holds the actual values
   */

  // --- integer parameters ---
  /** The maximum number of iterations used during the flash. */
  static integer constexpr FLASH_MAX_ITERATIONS = 0;

  /** The maximum number of iterations used during the stability test. */
  static integer constexpr STABILITY_MAX_ITERATIONS = 1;

  /* Number of discrete variables */
  static integer constexpr NUM_INTEGER = 2;

  // --- real parameters ---
  /** The tolerance to use to determine convergece of the flash. */
  static integer constexpr FLASH_TOLERANCE = 0;

  /** The limit for negative flash tolerance. */
  static integer constexpr NEGATIVE_FLASH_TOLERANCE = 1;

  /** The value of the tangent plane distance below which a mixture is determined to be unstable. */
  static integer constexpr STABILITY_THRESHOLD = 2;

  /** The tolerance to use to determine convergece to a stationary point in the stability test. */
  static integer constexpr STABILITY_TOLERANCE = 3;

  /* Number of continuous variables */
  static integer constexpr NUM_FLOAT = 4;

public:
  FlashParameters( std::unique_ptr< ModelParameters > parameters );
  ~FlashParameters() override = default;

  static std::unique_ptr< ModelParameters > create( std::unique_ptr< ModelParameters > parameters );

  struct viewKeyStruct
  {
    static constexpr char const * flashToleranceString() { return "flashTolerance"; }
    static constexpr char const * flashMaxIterationsString() { return "flashMaxIterations"; }
    static constexpr char const * stabilityThresholdString() { return "stabilityThreshold"; }
    static constexpr char const * stabilityToleranceString() { return "stabilityTolerance"; }
    static constexpr char const * stabilityMaxIterationsString() { return "stabilityMaxIterations"; }
  };

  array1d< real64 > m_continuousParameters;
  array1d< integer > m_discreteParameters;

protected:
  void registerParametersImpl( MultiFluidBase * fluid ) override;

  void postInputInitializationImpl( MultiFluidBase const * fluid, ComponentProperties const & componentProperties ) override;
};

} // end namespace compositional

} // end namespace constitutive

} // end namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_PARAMETERS_FLASHPARAMETERS_HPP_
