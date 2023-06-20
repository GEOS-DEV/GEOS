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
 * @file NegativeTwoPhaseFlash.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_NEGATIVETWOPHASEFLASH_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_NEGATIVETWOPHASEFLASH_HPP_

#include "common/DataTypes.hpp"

namespace geos
{

namespace constitutive
{

template< typename EOS_TYPE_LIQUID, typename EOS_TYPE_VAPOUR >
struct NegativeTwoPhaseFlash
{
public:
  /// Max number of components alloweeed in the class for now
  static constexpr integer maxNumComps = 5;
  /// Max number of iterations
  static constexpr integer maxIterations = 200;
  /// Epsilon used in the calculations
  static constexpr real64 epsilon = LvArray::NumericLimits< real64 >::epsilon;
  /// Tolerance for checking fugacity ratio convergence
  static constexpr real64 fugacityTolerance = 1.0e-8;

  /**
   * @brief Perform negative two-phase EOS flash
   * @param[in] numComps number of components
   * @param[in] pressure pressure
   * @param[in] temperature temperature
   * @param[in] composition composition of the mixture
   * @param[in] criticalPressure critical pressures
   * @param[in] criticalTemperature critical temperatures
   * @param[in] acentricFactor acentric factors
   * @param[in] binaryInteractionCoefficients binary coefficients (currently not implemented)
   * @param[out] vapourPhaseMoleFraction the calculated vapour (gas) mole fraction
   * @param[out] liquidComposition the calculated liquid phase composition
   * @param[out] vapourComposition the calculated vapour phase composition
   * @return an indicator of success of the flash
   */
  GEOS_HOST_DEVICE
  static bool compute( integer const numComps,
                       real64 const pressure,
                       real64 const temperature,
                       arrayView1d< real64 const > const composition,
                       arrayView1d< real64 const > const criticalPressure,
                       arrayView1d< real64 const > const criticalTemperature,
                       arrayView1d< real64 const > const acentricFactor,
                       real64 const & binaryInteractionCoefficients,
                       real64 & vapourPhaseMoleFraction,
                       arrayView1d< real64 > const liquidComposition,
                       arrayView1d< real64 > const vapourComposition );

private:
  /**
   * @brief Normalise a composition in place to ensure that the components add up to unity
   * @param[in] numComps number of components
   * @param[in/out] composition composition to be normalized
   * @return the sum of the given values
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static real64 normalizeComposition( integer const numComps,
                                      arraySlice1d< real64 > const composition );
};

} // namespace constitutive

} // namespace geos

#include "NegativeTwoPhaseFlash_impl.hpp"

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_NEGATIVETWOPHASEFLASH_HPP_
