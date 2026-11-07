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
 * @file StabilityTest.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_STABILITYTEST_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_STABILITYTEST_HPP_

#include "KValueInitialization.hpp"
#include "FugacityCalculator.hpp"
#include "Utilities.hpp"
#include "constitutive/fluid/multifluid/Layouts.hpp"
#include "constitutive/fluid/multifluid/MultiFluidConstants.hpp"
#include "constitutive/fluid/multifluid/compositional/parameters/ComponentProperties.hpp"
#include "constitutive/fluid/multifluid/compositional/parameters/FlashParameters.hpp"
#include "constitutive/fluid/multifluid/compositional/functions/FlashData.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

/**
 * @struct StabilityTest
 * @brief Performs stability analysis for a reservoir fluid.
 *
 * This class implements a stability test to determine if a single-phase fluid mixture is thermodynamically
 * stable or if it will split into multiple phases. The algorithm follows the tangent-plane-distance (TPD)
 * method, as described in the reference by Michelsen (1982).
 *
 * The method involves finding a trial composition ($y_i$) that minimizes the tangent-plane-distance ($TPD$)
 * at a given pressure and temperature. The TPD for a trial phase with composition $y$ relative to a sample
 * phase with composition $z$ is given by:
 * $$ TPD(y) = \sum_{i=1}^{nc} y_i \left[ \ln(y_i) - \ln(z_i) + \ln(\phi}_i(y)) - \ln(\phi_i(z)) \right] $$
 * where $nc$ is the number of components, $z_i$ and $y_i$ are the mole fractions of component $i$ in
 * the sample and trial phases, respectively, and $\phi_i$ and $\hat{\phi}_i$ are the fugacity coefficients
 * of component $i$ in the sample and trial phases, respectively.
 *
 * The stability analysis is performed by solving for the stationary points of the TPD function. The solution
 * is found by a combination of successive substitution iterations (SSI) and Newton-Raphson iterations.
 * SSI is used initially to provide a good starting point for the faster, but more sensitive, Newton-Raphson
 * method. The mixture is considered unstable if the minimum TPD is less than zero, indicating that a more
 * stable phase exists.
 *
 * Reference: Michelsen, M.L. (1982). The Isothermal Flash Problem. Part I. Stability.
 *            Fluid Phase Equilibria, 9(1), 1–19. https://doi.org/10.1016/0378-3812(82)85001-2
 */
struct StabilityTest
{
private:
  static constexpr integer maxNumComps = MultiFluidConstants::MAX_NUM_COMPONENTS;
  static constexpr integer maxNumDofs = MultiFluidConstants::MAX_NUM_COMPONENTS + 2;
  using Deriv = constitutive::multifluid::DerivativeOffset;

public:
  /**
   * @brief Perform a two-phase stability test using the Tangent Plane Distance (TPD) method.
   *
   * The stability analysis is performed by solving for the stationary points of the TPD function. The solution
   * is found by a combination of successive substitution iterations (SSI) and Newton-Raphson iterations.
   * SSI is used initially to provide a good starting point for the faster, but more sensitive, Newton-Raphson
   * method. The mixture is considered unstable if the minimum TPD is less than zero, indicating that a more
   * stable phase exists.
   *
   * This implementation supports asymmetric equations of state, where the EOS used for the
   * incipient phase may differ from that of the original mixture. This is handled by testing
   * both configurations independently.
   *
   * @tparam USD1 Unit Stride Dimension for the composition array.
   * @tparam USD2 Unit Stride Dimension for the K-values array.
   *
   * @param[in] numComps Number of components in the mixture.
   * @param[in] pressure Pressure of the system [Pa].
   * @param[in] temperature Temperature of the system [K].
   * @param[in] composition Mole fractions of each component in the mixture (must be normalized).
   * @param[in] componentProperties Thermodynamic properties of each component.
   * @param[in] flashData Flash calculation context, including EOS types for liquid and vapor phases.
   * @param[in] kValues Initial equilibrium ratios (K-values) used to generate trial compositions.
   * @param[in] continuousFlashParameters Continuous parameters (e.g., tolerances, thresholds) for the flash calculation.
   * @param[in] discreteFlashParameters Discrete parameters (e.g., max iterations) for the flash calculation.
   * @param[out] unstableMixture Flag indicating whether the mixture is unstable.
   * @param[out] incipientEquationOfState EOS type of the incipient phase (liquid or vapor).
   * @param[out] incipientComposition Composition of the incipient phase at the stationary point furthest from the trivial solution.
   *
   * @return True if the mixture is unstable or all trial compositions converge to stationary points; false otherwise.
   *
   */
  template< integer USD1, integer USD2 >
  GEOS_HOST_DEVICE
  static bool compute( integer const numComps,
                       real64 const pressure,
                       real64 const temperature,
                       arraySlice1d< real64 const, USD1 > const & composition,
                       ComponentProperties::KernelWrapper const & componentProperties,
                       FlashData const & flashData,
                       arraySlice1d< real64 const, USD2 > const & kValues,
                       arraySlice1d< real64 const > const & continuousFlashParameters,
                       arraySlice1d< integer const > const & discreteFlashParameters,
                       bool & unstableMixture,
                       EquationOfStateType & incipientEquationOfState,
                       arraySlice1d< real64 > const & incipientComposition );

  /**
   * @brief Compute derivatives of the incipient phase composition with respect to pressure, temperature, and overall composition.
   *
   * @details This method calculates the sensitivity of the incipient phase composition resulting from a phase stability test.
   * It computes the partial derivatives of the incipient composition with respect to pressure, temperature, and each component
   * of the overall mixture composition.
   *
   * The implementation solves a linear system derived from the stationarity condition of the Tangent Plane Distance (TPD)
   * function. It uses fugacity derivatives from both the sample and incipient phases, and enforces the constraint that the
   * incipient composition must sum to unity.
   *
   * @tparam USD1 Unit Stride Dimension for the composition array.
   * @tparam USD2 Unit Stride Dimension for the derivative arrays.
   *
   * @param[in] numComps Number of components in the mixture.
   * @param[in] pressure Pressure of the system [Pa].
   * @param[in] temperature Temperature of the system [K].
   * @param[in] composition Mole fractions of each component in the original mixture.
   * @param[in] componentProperties Thermodynamic properties of each component.
   * @param[in] sampleEquationOfState EOS used for the original mixture.
   * @param[in] incipientEquationOfState EOS used for the incipient phase.
   * @param[in] flashData Flash calculation context and EOS configuration.
   * @param[in] incipientComposition Composition of the incipient phase.
   * @param[out] incipientCompositionDerivs Derivatives of the incipient composition with respect to pressure, temperature, and composition.
   * @param[out] compositionDerivs Workspace for storing derivatives of the sample composition.
   *
   */
  template< integer USD1, integer USD2 >
  GEOS_HOST_DEVICE
  static void computeDerivatives( integer const numComps,
                                  real64 const pressure,
                                  real64 const temperature,
                                  arraySlice1d< real64 const, USD1 > const & composition,
                                  ComponentProperties::KernelWrapper const & componentProperties,
                                  EquationOfStateType const & sampleEquationOfState,
                                  EquationOfStateType const & incipientEquationOfState,
                                  FlashData const & flashData,
                                  arraySlice1d< real64 const > const & incipientComposition,
                                  arraySlice2d< real64, USD2 > const & incipientCompositionDerivs,
                                  arraySlice2d< real64, USD2 > const & compositionDerivs );

  /**
   * @brief Calculates the residual vector and Tangent Plane Distance (TPD) for a trial composition.
   *
   * @details This method evaluates the deviation of a trial composition from stationarity
   * in the context of phase stability analysis. It computes the residual vector used in
   * iterative solvers (SSI/Newton) and the tangent plane distance (TPD) which indicates
   * whether the trial composition has lower Gibbs energy than the original mixture.
   *
   * The trial composition is provided in logarithmic form and exponentiated internally.
   * Missing components (not present in presentComponents) are ignored in the residual
   * and TPD calculation. The composition is normalized after exponentiation.
   *
   * @tparam USD1 Unit Stride Dimension for logTestComposition array.
   * @tparam USD2 Unit Stride Dimension for testComposition array.
   * @tparam USD3 Unit Stride Dimension for residual array.
   *
   * @param[in] numComps Number of components in the mixture.
   * @param[in] presentComponents Indices of components present in the mixture.
   * @param[in] pressure System pressure [Pa].
   * @param[in] temperature System temperature [K].
   * @param[in] componentProperties Thermodynamic properties of each component.
   * @param[in] equationOfState Equation of state used for fugacity calculations.
   * @param[in] flashData Flash calculation context.
   * @param[in] logTestComposition Logarithm of trial composition.
   * @param[out] tangentPlaneDistance Computed tangent plane distance (TPD).
   * @param[out] testComposition Trial composition (exponentiated and normalized).
   * @param[out] testLogFugacity Log fugacity coefficients of the trial composition.
   * @param[out] residual Residual vector for iterative solver.
   *
   * @return Error norm of the residual vector.
   *
   */
  template< int USD1, int USD2, int USD3 >
  GEOS_HOST_DEVICE
  static real64 calculateResidual( integer const numComps,
                                   arraySlice1d< integer const > const & presentComponents,
                                   real64 const pressure,
                                   real64 const temperature,
                                   ComponentProperties::KernelWrapper const & componentProperties,
                                   EquationOfStateType const equationOfState,
                                   FlashData const & flashData,
                                   arraySlice1d< real64 const, USD1 > const & logTestComposition,
                                   real64 & tangentPlaneDistance,
                                   arraySlice1d< real64, USD2 > const & testComposition,
                                   arraySlice1d< real64 > const & testLogFugacity,
                                   arraySlice1d< real64, USD3 > const & residual );

  /**
   * @brief Calculates the residual vector and Jacobian matrix for a trial composition in phase stability analysis.
   *
   * @details This method evaluates the stationarity of a trial composition using the Tangent Plane Distance (TPD) criterion.
   * It computes the residual vector and the Jacobian matrix required for Newton-Raphson iterations.
   * The fugacity and its derivatives are calculated using the specified equation of state.
   * Missing components (not present in the trial composition) are ignored in the residual and Jacobian calculations.
   *
   * @tparam USD1 Unit Stride Dimension for logTestComposition array.
   * @tparam USD2 Unit Stride Dimension for testComposition array.
   * @tparam USD3 Unit Stride Dimension for residual array.
   * @tparam USD4 Unit Stride Dimension for jacobian array.
   *
   * @param[in] numComps Number of components in the mixture.
   * @param[in] presentComponents Indices of components present in the trial composition.
   * @param[in] pressure Pressure of the system [Pa].
   * @param[in] temperature Temperature of the system [K].
   * @param[in] componentProperties Thermodynamic properties of each component.
   * @param[in] equationOfState Equation of state used for fugacity calculations.
   * @param[in] flashData Flash calculation context and parameters.
   * @param[in] logTestComposition Logarithm of the trial composition.
   * @param[out] tangentPlaneDistance Computed TPD value for the trial composition.
   * @param[out] testComposition Normalized trial composition.
   * @param[out] testFugacity Log fugacity coefficients of the trial composition.
   * @param[out] logTestFugacityDerivs Derivatives of log fugacity with respect to composition.
   * @param[out] residual Residual vector used in Newton-Raphson iterations.
   * @param[out] jacobian Jacobian matrix of the residual vector.
   *
   * @return Error norm of the residual vector.
   */
  template< int USD1, int USD2, int USD3, int USD4 >
  GEOS_HOST_DEVICE
  static real64 calculateResidualAndJacobian( integer const numComps,
                                              arraySlice1d< integer const > const & presentComponents,
                                              real64 const pressure,
                                              real64 const temperature,
                                              ComponentProperties::KernelWrapper const & componentProperties,
                                              EquationOfStateType const equationOfState,
                                              FlashData const & flashData,
                                              arraySlice1d< real64 const, USD1 > const & logTestComposition,
                                              real64 & tangentPlaneDistance,
                                              arraySlice1d< real64, USD2 > const & testComposition,
                                              arraySlice1d< real64 > const & testFugacity,
                                              arraySlice2d< real64 > const & logTestFugacityDerivs,
                                              arraySlice1d< real64, USD3 > const & residual,
                                              arraySlice2d< real64, USD4 > const & jacobian );
};

} // namespace compositional

} // namespace constitutive

} // namespace geos

// Include the implementation
#include "StabilityTest_impl.hpp"

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_STABILITYTEST_HPP_
