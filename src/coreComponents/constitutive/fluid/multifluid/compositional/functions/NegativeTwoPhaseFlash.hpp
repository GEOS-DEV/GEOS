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
   * @brief Perform a two-phase EOS flash calculation.
   *
   * @details
   * Executes a negative flash calculation using an equation of state (EOS) model
   * to determine the equilibrium phase split of a multicomponent mixture at a given
   * pressure and temperature. The function computes the vapor phase mole fraction,
   * as well as the compositions of the vapor and liquid phases, based on the provided
   * initial estimates and component properties. Uses successive substitution.
   *
   * @tparam USD1 Unique stride descriptor for the K-values array.
   * @tparam USD2 Unique stride descriptor for the phase composition arrays.
   *
   * @param[in] numComps Number of components in the mixture.
   * @param[in] pressure System pressure [Pa].
   * @param[in] temperature System temperature [K].
   * @param[in] composition Overall composition of the mixture (z_i).
   * @param[in] componentProperties Thermodynamic properties of each component.
   * @param[in] flashData Precomputed data and parameters required for the flash calculation.
   * @param[in] continuousFlashParameters List of continuous (floating-point) parameters for the flash.
   * @param[in] discreteFlashParameters List of discrete (integer) parameters for the flash.
   * @param[in,out] kValues Phase equilibrium ratios (K-values), updated during the flash.
   * @param[out] vapourPhaseMoleFraction Calculated vapor phase mole fraction (V).
   * @param[out] liquidComposition Calculated composition of the liquid phase (x_i).
   * @param[out] vapourComposition Calculated composition of the vapor phase (y_i).
   *
   * @return True if the flash calculation was successful; false otherwise.
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
   * @brief Calculate the derivatives of the flash.
   *
   * @details
   * Implicitly computes the derivatives of the flash problem after the iterative
   * solution has converged. This includes the derivatives of the vapor fraction
   * and the phase compositions with respect to all primary variables
   *
   * @tparam USD1 Unique stride descriptor for phase compositions.
   * @tparam USD2 Unique stride descriptor for vapor fraction derivatives.
   * @tparam USD3 Unique stride descriptor for composition derivatives.
   *
   * @param[in] numComps Number of components in the system.
   * @param[in] totalComposition Overall composition vector (z_i).
   * @param[in] phase1Fraction Mole fraction of phase 1 (typically vapour).
   * @param[in] phase1Composition Composition of phase 1 (y_i).
   * @param[in] phase2Composition Composition of phase 2 (x_i).
   * @param[in] phase1Fugacity Fugacity coefficients in phase 1 (exp(\phi_V,i)).
   * @param[in] phase2Fugacity Fugacity coefficients in phase 2 (exp(\phi_L,i)).
   * @param[in] phase1LogFugacityDerivs Derivatives of log fugacity in phase 1.
   * @param[in] phase2LogFugacityDerivs Derivatives of log fugacity in phase 2.
   * @param[out] phase1FractionDerivs Derivatives of phase1Fraction.
   * @param[out] phase1CompositionDerivs Derivatives of phase1Composition.
   * @param[out] phase2CompositionDerivs Derivatives of phase2Composition.
   */
  template< integer USD1, integer USD2, integer USD3 >
  GEOS_HOST_DEVICE
  static void computeDerivatives(
    integer const numComps,
    arraySlice1d< real64 const > const & totalComposition,
    real64 const & phase1Fraction,
    arraySlice1d< real64 const, USD1 > const & phase1Composition,
    arraySlice1d< real64 const, USD1 > const & phase2Composition,
    arraySlice1d< real64 const > const & phase1Fugacity,
    arraySlice1d< real64 const > const & phase2Fugacity,
    arraySlice2d< real64 const > const & phase1LogFugacityDerivs,
    arraySlice2d< real64 const > const & phase2LogFugacityDerivs,
    arraySlice1d< real64, USD2 > const & phase1FractionDerivs,
    arraySlice2d< real64, USD3 > const & phase1CompositionDerivs,
    arraySlice2d< real64, USD3 > const & phase2CompositionDerivs );

  /**
   * @brief A test for when the result of the negative flash should be truncated immediately
   * @param[in] numComps Number of components in the mixture.
   * @param[in] totalComposition Total overall composition of the system (z_i).
   * @param[in] vapourPhaseMoleFraction Mole fraction of the vapor phase (V).
   * @param[in] liquidComposition Mole fractions in the liquid phase (x_i).
   * @param[in] vapourComposition Mole fractions in the vapor phase (y_i).
   * @return @c true if truncation should be applied
   */
  template< integer USD1, integer USD2 >
  GEOS_FORCE_INLINE
  GEOS_HOST_DEVICE
  static bool truncateCompositions( integer const numComps,
                                    arraySlice1d< real64 const, USD1 > const & totalComposition,
                                    real64 const & vapourPhaseMoleFraction,
                                    arraySlice1d< real64 const, USD2 > const & liquidComposition,
                                    arraySlice1d< real64 const, USD2 > const & vapourComposition );
};

} // namespace compositional

} // namespace constitutive

} // namespace geos

// Include the implementation
#include "NegativeTwoPhaseFlash_impl.hpp"

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_NEGATIVETWOPHASEFLASH_HPP_
