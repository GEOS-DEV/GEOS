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

/**
 * @brief Struct implementing a negative two-phase flash calculation.
 *
 * @details This struct provides methods to perform a two-phase flash calculation
 * for reservoir fluid systems where the mixture is thermodynamically unstable.
 * It uses the algorithm described in Michelsen (1982)
 *
 * The flash calculation is performed in two stages:
 * - Successive Substitution Iteration (SSI): Used for initial convergence.
 * - Newton-Raphson Iteration: Used to refine the solution and accelerate convergence.
 *
 * The algorithm solves the Rachford-Rice equation to determine the vapor phase mole fraction,
 * and iteratively updates K-values and phase compositions to satisfy fugacity equilibrium.
 *
 * Missing components (i.e., those with zero mole fraction in the feed) are excluded from
 * fugacity and Jacobian calculations but handled explicitly in the residual system.
 *
 * Reference: Michelsen, M.L. (1982). The Isothermal Flash Problem. Part II. Phase Split Calculation.
 *            Fluid Phase Equilibria, 9(1), 21–40. https://doi.org/10.1016/0378-3812(82)85002-4
 */
struct NegativeTwoPhaseFlash
{
  using Deriv = constitutive::multifluid::DerivativeOffset;

public:
  /**
   * @brief Perform a negative two-phase flash calculation for unstable mixtures.
   *
   * @details This method solves the isothermal two-phase flash problem for mixtures
   * that are thermodynamically unstable.
   *
   * The algorithm proceeds in two stages:
   * - Successive Substitution Iteration (SSI): Used for initial convergence of K-values.
   * - Newton-Raphson Iteration: Refines the solution using a Jacobian matrix and residual vector.
   *
   * The Jacobian matrix includes partial derivatives of fugacity with respect to composition
   * and vapor fraction. Missing components (zero mole fraction) are handled by setting
   * diagonal Jacobian entries to unity to ensure numerical stability.
   *
   * @tparam USD1 Unit Stride Dimension for the K-values array.
   * @tparam USD2 Unit Stride Dimension for the composition arrays.
   *
   * @param[in] numComps Number of components in the mixture.
   * @param[in] pressure System pressure [Pa].
   * @param[in] temperature System temperature [K].
   * @param[in] composition Overall feed composition (mole fractions).
   * @param[in] componentProperties Thermodynamic properties of each component.
   * @param[in] flashData Flash calculation context including EOS types.
   * @param[in] continuousFlashParameters Continuous parameters for flash calculation.
   * @param[in] discreteFlashParameters Discrete parameters for flash calculation.
   * @param[in,out] kValues Initial and updated K-values (vapor-liquid equilibrium ratios).
   * @param[out] vapourPhaseMoleFraction Computed vapor phase mole fraction.
   * @param[out] liquidComposition Computed liquid phase composition.
   * @param[out] vapourComposition Computed vapor phase composition.
   *
   * @return True if the flash calculation converged; false otherwise.
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
   * @brief Compute derivatives of the vapor fraction and phase compositions for a negative two-phase flash.
   *
   * @details This method computes the sensitivities of the vapor fraction and the liquid/vapor phase compositions
   * with respect to pressure, temperature, and overall composition. It uses implicit differentiation of the
   * equilibrium conditions and mass balance constraints to solve for the derivatives. The Jacobian matrix is
   * constructed from the linearized system of equations derived from fugacity equilibrium and Rachford-Rice relations.
   *
   * @tparam USD1 Unit Stride Dimension for input phase composition arrays.
   * @tparam USD2 Unit Stride Dimension for vapor fraction derivative array.
   * @tparam USD3 Unit Stride Dimension for phase composition derivative arrays.
   *
   * @param[in] numComps Number of components in the mixture.
   * @param[in] pressure System pressure [Pa].
   * @param[in] temperature System temperature [K].
   * @param[in] composition Overall feed composition.
   * @param[in] componentProperties Thermodynamic properties of each component.
   * @param[in] flashData Flash calculation context including EOS types.
   * @param[in] vapourFraction Vapor phase mole fraction.
   * @param[in] liquidComposition Composition of the liquid phase.
   * @param[in] vapourComposition Composition of the vapor phase.
   * @param[out] vapourFractionDerivs Derivatives of vapor fraction with respect to pressure, temperature, and composition.
   * @param[out] liquidCompositionDerivs Derivatives of liquid phase composition.
   * @param[out] vapourCompositionDerivs Derivatives of vapor phase composition.
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

  /**
   * @brief Compute the residual vector for a two-phase flash calculation.
   *
   * @details This method evaluates the phase equilibrium condition by computing
   * the residuals of fugacity balance between liquid and vapor phases for each component.
   * It also includes a mass balance residual for the vapor phase mole fraction.
   * The compositions are calculated using the Rachford-Rice formulation and normalized.
   *
   * The fugacity coefficients are computed for both phases using the specified equation of state.
   * The residuals are defined as:
   *   R_i = x_i * \phi_i^L - y_i * \phi_i^V
   *   R_N = \sum_{i=1}^N (y_i - x_i)
   * where \phi_i^L and \phi_i^V are fugacity coefficients in liquid and vapor phases.
   *
   * @tparam USD1 Unit Stride Dimension for k-values array.
   * @tparam USD2 Unit Stride Dimension for composition arrays.
   * @tparam USD3 Unit Stride Dimension for residual array.
   *
   * @param[in] numComps Number of components in the mixture.
   * @param[in] pressure System pressure [Pa].
   * @param[in] temperature System temperature [K].
   * @param[in] composition Overall feed composition (mole fractions).
   * @param[in] componentProperties Thermodynamic properties of components.
   * @param[in] flashData Flash calculation context including EOS types.
   * @param[in] kValues Equilibrium ratios (K-values) for each component.
   * @param[in] vapourPhaseMoleFraction Mole fraction of the vapor phase.
   * @param[out] liquidComposition Computed composition of the liquid phase.
   * @param[out] vapourComposition Computed composition of the vapor phase.
   * @param[out] logLiquidFugacity Log fugacity coefficients for the liquid phase.
   * @param[out] logVapourFugacity Log fugacity coefficients for the vapor phase.
   * @param[out] residual Residual vector containing fugacity and mass balance errors.
   *
   * @return L2 norm of the residual vector.
   */
  template< int USD1, int USD2, int USD3 >
  GEOS_HOST_DEVICE
  static real64 calculateResidual( integer const numComps,
                                   real64 const pressure,
                                   real64 const temperature,
                                   arraySlice1d< real64 const > const & composition,
                                   ComponentProperties::KernelWrapper const & componentProperties,
                                   FlashData const & flashData,
                                   arraySlice1d< real64 const, USD1 > const & kValues,
                                   real64 const & vapourPhaseMoleFraction,
                                   arraySlice1d< real64, USD2 > const & liquidComposition,
                                   arraySlice1d< real64, USD2 > const & vapourComposition,
                                   arraySlice1d< real64 > const & logLiquidFugacity,
                                   arraySlice1d< real64 > const & logVapourFugacity,
                                   arraySlice1d< real64, USD3 > const & residual );

  /**
   * @brief Calculates the residual vector and Jacobian matrix for a two-phase flash system.
   *
   * @details This method computes the nonlinear residuals and Jacobian matrix for a negative two-phase flash
   * system using fugacity equilibrium and mass balance constraints. It is used in Newton-Raphson iterations
   * to refine the solution of the flash problem. The primary variables in the Jacobian are the logarithm of
   * the K-values (log(K)) and the vapor phase mole fraction. The Jacobian is constructed with respect to these
   * variables. Missing components (i.e., those with zero mole fraction in the feed) are excluded from the
   * fugacity and Jacobian calculations but are handled explicitly by setting diagonal Jacobian entries to unity.
   *
   * The method normalizes phase compositions, computes fugacity coefficients and their derivatives, and
   * assembles the Jacobian matrix using chain rule and normalization sensitivities.
   *
   * @tparam USD1 Unit Stride Dimension for the K-values array.
   * @tparam USD2 Unit Stride Dimension for the composition arrays.
   * @tparam USD3 Unit Stride Dimension for the residual vector.
   * @tparam USD4 Unit Stride Dimension for the Jacobian matrix.
   *
   * @param[in] numComps Number of components in the mixture.
   * @param[in] pressure System pressure [Pa].
   * @param[in] temperature System temperature [K].
   * @param[in] composition Overall feed composition (mole fractions).
   * @param[in] componentProperties Thermodynamic properties of each component.
   * @param[in] flashData Flash calculation context including EOS types.
   * @param[in] kValues Equilibrium ratios (K-values) for each component.
   * @param[in] vapourPhaseMoleFraction Mole fraction of the vapor phase.
   * @param[in] liquidComposition Composition of the liquid phase.
   * @param[in] vapourComposition Composition of the vapor phase.
   * @param[in] logLiquidFugacity Logarithm of fugacity coefficients in the liquid phase.
   * @param[in] logVapourFugacity Logarithm of fugacity coefficients in the vapor phase.
   * @param[in] logLiquidFugacityDerivs Derivatives of log fugacity in the liquid phase.
   * @param[in] logVapourFugacityDerivs Derivatives of log fugacity in the vapor phase.
   * @param[out] residual Residual vector containing fugacity and mass balance errors.
   * @param[out] jacobian Jacobian matrix of the residuals with respect to log(K) and vapor fraction.
   *
   * @return Error norm of the residual vector.
   */
  template< int USD1, int USD2, int USD3, int USD4 >
  GEOS_HOST_DEVICE
  static real64 calculateResidualAndJacobian( integer const numComps,
                                              real64 const pressure,
                                              real64 const temperature,
                                              arraySlice1d< real64 const > const & composition,
                                              ComponentProperties::KernelWrapper const & componentProperties,
                                              FlashData const & flashData,
                                              arraySlice1d< real64 const, USD1 > const & kValues,
                                              real64 const & vapourPhaseMoleFraction,
                                              arraySlice1d< real64, USD2 > const & liquidComposition,
                                              arraySlice1d< real64, USD2 > const & vapourComposition,
                                              arraySlice1d< real64 > const & logLiquidFugacity,
                                              arraySlice1d< real64 > const & logVapourFugacity,
                                              arraySlice2d< real64 > const & logLiquidFugacityDerivs,
                                              arraySlice2d< real64 > const & logVapourFugacityDerivs,
                                              arraySlice1d< real64, USD3 > const & residual,
                                              arraySlice2d< real64, USD4 > const & jacobian );

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
  static void computeDerivatives( integer const numComps,
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
