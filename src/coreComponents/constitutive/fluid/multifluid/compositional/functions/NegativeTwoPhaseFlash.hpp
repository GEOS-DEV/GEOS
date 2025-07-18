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
 * @file NegativeTwoPhaseFlash.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_NEGATIVETWOPHASEFLASH_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_NEGATIVETWOPHASEFLASH_HPP_

#include "FlashData.hpp"
#include "constitutive/fluid/multifluid/Layouts.hpp"
#include "constitutive/fluid/multifluid/compositional/parameters/ComponentProperties.hpp"

namespace geos
{

namespace constitutive
{
namespace compositional
{

struct NegativeTwoPhaseFlash
{
  using Deriv = constitutive::multifluid::DerivativeOffset;

public:
  /**
   * @brief Perform negative two-phase EOS flash
   * @param[in] numComps number of components
   * @param[in] pressure pressure
   * @param[in] temperature temperature
   * @param[in] composition composition of the mixture
   * @param[in] componentProperties The compositional component properties
   * @param[in] flashData The parameters required for the flash
   * @param[in] continuousFlashParameters List of continuous (float) parameters for flash
   * @param[in] discreteFlashParameters List of discrete (integer) parameters for flash
   * @param[in/out] kValues The phase equilibrium ratios
   * @param[out] vapourPhaseMoleFraction the calculated vapour (gas) mole fraction
   * @param[out] liquidComposition the calculated liquid phase composition
   * @param[out] vapourComposition the calculated vapour phase composition
   * @return an indicator of success of the flash
   */
  template< int USD1, int USD2 >
  GEOS_HOST_DEVICE
  static bool compute( integer const numComps,
                       real64 const pressure,
                       real64 const temperature,
                       arraySlice1d< real64 const > const & composition,
                       ComponentProperties::KernelWrapper const & componentProperties,
                       FlashData const & flashData,
                       arraySlice1d< real64 const > const & continuousFlashParameters,
                       arraySlice1d< integer const > const & discreteFlashParameters,
                       arraySlice2d< real64, USD1 > const & kValues,
                       real64 & vapourPhaseMoleFraction,
                       arraySlice1d< real64, USD2 > const & liquidComposition,
                       arraySlice1d< real64, USD2 > const & vapourComposition );

  /**
   * @brief Calculate derivatives from the two-phase negative flash
   * @param[in] numComps number of components
   * @param[in] pressure pressure
   * @param[in] temperature temperature
   * @param[in] composition composition of the mixture
   * @param[in] componentProperties The compositional component properties
   * @param[in] flashData The parameters required for the flash
   * @param[in] vapourFraction the calculated vapour (gas) mole fraction
   * @param[in] liquidComposition the calculated liquid phase composition
   * @param[in] vapourComposition the calculated vapour phase composition
   * @param[out] vapourFractionDerivs derivatives of the calculated vapour (gas) mole fraction
   * @param[out] liquidCompositionDerivs derivatives of the calculated liquid phase composition
   * @param[out] vapourCompositionDerivs derivatives of the calculated vapour phase composition
   */
  template< integer USD1, integer USD2, integer USD3 >
  GEOS_HOST_DEVICE
  static void computeDerivatives( integer const numComps,
                                  real64 const pressure,
                                  real64 const temperature,
                                  arraySlice1d< real64 const > const & composition,
                                  ComponentProperties::KernelWrapper const & componentProperties,
                                  FlashData const & flashData,
                                  real64 const & vapourFraction,
                                  arraySlice1d< real64 const, USD1 > const & liquidComposition,
                                  arraySlice1d< real64 const, USD1 > const & vapourComposition,
                                  arraySlice1d< real64, USD2 > const & vapourFractionDerivs,
                                  arraySlice2d< real64, USD3 > const & liquidCompositionDerivs,
                                  arraySlice2d< real64, USD3 > const & vapourCompositionDerivs );

private:
  /**
   * @brief Calculate the logarithms of the fugacity ratios
   * @param[in] numComps number of components
   * @param[in] pressure pressure
   * @param[in] temperature temperature
   * @param[in] composition composition of the mixture
   * @param[in] componentProperties The compositional component properties
   * @param[in] flashData The parameters required for the flash
   * @param[in] kValues The k-values
   * @param[in] presentComponents The indices of the present components
   * @param[out] vapourPhaseMoleFraction the calculated vapour (gas) mole fraction
   * @param[out] liquidComposition the calculated liquid phase composition
   * @param[out] vapourComposition the calculated vapour phase composition
   * @param[out] logLiquidFugacity the calculated log fugacity ratios for the liquid phase
   * @param[out] logVapourFugacity the calculated log fugacity ratios for the vapour phase
   * @param[out] fugacityRatios the fugacity rations
   * @return The error
   */
  template< integer USD1, integer USD2 >
  GEOS_HOST_DEVICE
  static real64 computeFugacityRatio(
    integer const numComps,
    real64 const pressure,
    real64 const temperature,
    arraySlice1d< real64 const > const & composition,
    ComponentProperties::KernelWrapper const & componentProperties,
    FlashData const & flashData,
    arraySlice1d< real64 const, USD1 > const & kValues,
    arraySlice1d< integer const > const & presentComponents,
    real64 & vapourPhaseMoleFraction,
    arraySlice1d< real64, USD2 > const & liquidComposition,
    arraySlice1d< real64, USD2 > const & vapourComposition,
    arraySlice1d< real64 > const & logLiquidFugacity,
    arraySlice1d< real64 > const & logVapourFugacity,
    arraySlice1d< real64 > const & fugacityRatios );
};

} // namespace compositional

} // namespace constitutive

} // namespace geos

// Include the implementation
#include "NegativeTwoPhaseFlash_impl.hpp"

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_NEGATIVETWOPHASEFLASH_HPP_
