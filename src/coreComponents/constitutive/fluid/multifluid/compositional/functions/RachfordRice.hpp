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
 * @file RachfordRice.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_RACHFORDRICE_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_RACHFORDRICE_HPP_

#include "common/DataTypes.hpp"
#include "constitutive/fluid/multifluid/MultiFluidConstants.hpp"
#include "constitutive/fluid/multifluid/Layouts.hpp"

namespace geos
{

namespace constitutive
{

struct RachfordRice
{
public:
  /// Tolerance of the SSI loop
  static constexpr real64 SSITolerance = MultiFluidConstants::SSITolerance;
  /// Tolerance of the Newton loop
  static constexpr real64 newtonTolerance = MultiFluidConstants::newtonTolerance;
  /// Max number of SSI iterations
  static constexpr integer maxSSIIterations = MultiFluidConstants::maxSSIIterations;
  /// Max number of Newton iterations
  static constexpr integer maxNewtonIterations = MultiFluidConstants::maxNewtonIterations;
  /// Epsilon used in the calculations
  static constexpr real64 epsilon = LvArray::NumericLimits< real64 >::epsilon;
  /// Used to avoid exploring boundary of feasibility region
  static constexpr real64 threshold = 1.0 - 1.0e-4;

  /**
   * @brief Function solving the Rachford-Rice equation
   * @input[in] kValues the array of K-values
   * @input[in] feed the component fractions
   * @input[in] presentComponentIds the indices of components with a non-zero fractions
   * @return the gas mole fraction
   **/
  template< integer USD1, integer USD2 >
  GEOS_HOST_DEVICE
  real64
  static
  solve( arraySlice1d< real64 const, USD2 > const & kValues,
         arraySlice1d< real64 const, USD1 > const & feed,
         arraySlice1d< integer const > const & presentComponentIds )
  {
    real64 gasPhaseMoleFraction = 0;

    // min and max Kvalues for non-zero composition
    real64 minK, maxK;
    findKValueRange( kValues, presentComponentIds, minK, maxK );

    // check for trivial solutions.
    // this corresponds to bad Kvalues
    if( maxK < 1.0 )
    {
      return 0.0;
    }
    if( minK > 1.0 )
    {
      return 1.0;
    }

    // start the solver loop

    // step 1: find solution window
    real64 xMin = 1.0 / ( 1 - maxK );
    real64 xMax = 1.0 / ( 1 - minK );
    real64 const sqrtEpsilon = sqrt( epsilon );
    xMin += sqrtEpsilon * ( LvArray::math::abs( xMin ) + sqrtEpsilon );
    xMax -= sqrtEpsilon * ( LvArray::math::abs( xMax ) + sqrtEpsilon );

    real64 currentError = 1 / epsilon;

    // step 2: start the SSI loop
    // Evaluate at the bounds
    real64 funcXMin = evaluate( kValues, feed, presentComponentIds, xMin );
    real64 funcXMax = evaluate( kValues, feed, presentComponentIds, xMax );
    real64 funcXMid = 0.0;

    // If the bound values are the same sign then we have a trivial solution
    if( 0.0 < funcXMin * funcXMax )
    {
      gasPhaseMoleFraction = (0.0 < funcXMin) ? 1.0 : 0.0;
      return gasPhaseMoleFraction;
    }

    integer SSIIteration = 0;

    while( ( currentError > SSITolerance ) && ( SSIIteration < maxSSIIterations ) )
    {
      real64 const xMid = 0.5 * ( xMin + xMax );
      funcXMid = evaluate( kValues, feed, presentComponentIds, xMid );

      if( 0.0 < funcXMax * funcXMid )
      {
        xMax = xMid;
        funcXMax = funcXMid;
      }
      else if( 0.0 < funcXMin * funcXMid )
      {
        xMin = xMid;
        funcXMin = funcXMid;
      }

      currentError = LvArray::math::min( LvArray::math::abs( funcXMid ),
                                         LvArray::math::abs( xMax - xMin ) );
      SSIIteration++;

      // TODO: add warning if max number of SSI iterations is reached
    }

    gasPhaseMoleFraction = 0.5 * ( xMax + xMin );

    // step 3: start the Newton loop
    integer newtonIteration = 0;
    real64 newtonValue = gasPhaseMoleFraction;
    real64 funcNewton = evaluate( kValues, feed, presentComponentIds, newtonValue );

    while( ( currentError > newtonTolerance ) && ( newtonIteration < maxNewtonIterations ) )
    {
      real64 deltaNewton = -funcNewton / evaluateDerivative( kValues, feed, presentComponentIds, newtonValue );

      // test if we are stepping out of the [xMin;xMax] interval
      if( newtonValue + deltaNewton < xMin )
      {
        deltaNewton = 0.5 * ( xMin - newtonValue );
      }
      else if( newtonValue + deltaNewton > xMax )
      {
        deltaNewton = 0.5 * ( xMax - newtonValue );
      }

      newtonValue = newtonValue + deltaNewton;

      funcNewton = evaluate( kValues, feed, presentComponentIds, newtonValue );

      currentError = LvArray::math::min( LvArray::math::abs( funcNewton ),
                                         LvArray::math::abs( deltaNewton ) );
      newtonIteration++;

      // TODO: add warning if max number of Newton iterations is reached
    }
    return gasPhaseMoleFraction = newtonValue;
  }

private:

  /**
   * @brief Function evaluating the Rachford-Rice function
   * @input[in] kValues the array fo K-values
   * @input[in] feed the component fractions
   * @input[in] presentComponentIds the indices of components with a non-zero fractions
   * @input[in] x the value at which the Rachford-Rice function is evaluated
   * @return the value of the Rachford-Rice function at x
   **/
  template< integer USD1, integer USD2 >
  GEOS_HOST_DEVICE
  real64
  static
  evaluate( arraySlice1d< real64 const, USD2 > const & kValues,
            arraySlice1d< real64 const, USD1 > const & feed,
            arraySlice1d< integer const > const & presentComponentIds,
            real64 const & x )
  {
    real64 value = 0.0;
    for( integer i = 0; i < presentComponentIds.size(); ++i )
    {
      integer const ic = presentComponentIds[i];
      real64 const k = ( kValues[ic] - 1.0 );
      value += feed[ic] * k / ( 1.0 + x * k );
    }
    return value;
  }

  /**
   * @brief Function evaluating the derivative of the Rachford-Rice function
   * @input[in] kValues the array fo K-values
   * @input[in] feed the component fractions
   * @input[in] presentComponentIds the indices of components with a non-zero fractions
   * @input[in] x the value at which the derivative of the Rachford-Rice function is evaluated
   * @return the value of the derivative of the Rachford-Rice function at x
   **/
  template< integer USD1, integer USD2 >
  GEOS_HOST_DEVICE
  real64
  static
  evaluateDerivative( arraySlice1d< real64 const, USD2 > const & kValues,
                      arraySlice1d< real64 const, USD1 > const & feed,
                      arraySlice1d< integer const > const & presentComponentIds,
                      real64 const & x )
  {
    real64 value = 0.0;
    for( integer i = 0; i < presentComponentIds.size(); ++i )
    {
      integer const ic = presentComponentIds[i];
      real64 const k = ( kValues[ic] - 1.0 );
      real64 const r = k / ( 1.0 + x * k );
      value -= feed[ic] * r * r;
    }
    return value;
  }

  /**
   * @brief Calculate the minimum and maximum k-value
   * @input[in] kValues the array fo K-values
   * @input[in] presentComponentIds the indices of components with a non-zero fractions
   * @input[out] minK the minimum k-value for non-zero components
   * @input[out] maxK the maximum k-value for non-zero components
   **/
  template< integer USD >
  GEOS_FORCE_INLINE
  GEOS_HOST_DEVICE
  void
  static
  findKValueRange( arraySlice1d< real64 const, USD > const & kValues,
                   arraySlice1d< integer const > const & presentComponentIds,
                   real64 & minK,
                   real64 & maxK )
  {
    minK = 1.0 / epsilon;
    maxK = 0.0;
    for( integer const ic : presentComponentIds )
    {
      if( kValues[ic] > maxK )
      {
        maxK = kValues[ic];
      }
      if( kValues[ic] < minK )
      {
        minK = kValues[ic];
      }
    }
  }

  /**
   * @brief Computes the objective function, gradient vector, and Hessian matrix.
   *
   * @details Evaluates Michelsen's objective function reformulated for the negative flash
   * region as described in Leibovici and Nichita (2008), Section 3, Equations 10-14.
   * Paper: "A new look at multiphase Rachford-Rice equations for negative flashes"
   * DOI: 10.1016/j.fluid.2008.03.006
   * The function F is convex and its unbounded minimization yields the phase fractions.
   * Reference phase is Phase 1 (implicitly index 1).
   *
   * @param kValues Array slice of equilibrium constants (K-values).
   * @param feed Array slice of feed composition (mole fractions).
   * @param presentComponentIds Array slice of indices for components present in the mixture.
   * @param v2 Mole fraction of the second phase.
   * @param v3 Mole fraction of the third phase.
   * @param F Computed objective function value (Equation 10).
   * @param g Computed gradient vector of size 2 (Equation 13).
   * @param H Computed Hessian matrix of size 4 (explicitly computed 2x2 matrix).
   * @return The L2 norm of the error
   */
  template< integer USD1, integer USD2 >
  GEOS_HOST_DEVICE
  static real64 computeObjective(
    arraySlice2d< real64 const, USD2 > const & kValues,
    arraySlice1d< real64 const, USD1 > const & feed,
    arraySlice1d< integer const > const & presentComponentIds,
    real64 const v2, real64 const v3,
    real64 & F, real64 g[2], real64 H[4] )
  {
    F = 0.0;
    g[0] = 0.0;
    g[1] = 0.0;
    H[0] = 0.0;
    H[1] = 0.0;
    H[2] = 0.0;
    H[3] = 0.0;

    for( integer const & ic : presentComponentIds )
    {
      real64 const z = feed[ic];

      real64 const k2 = kValues( 0, ic ) - 1.0;
      real64 const k3 = kValues( 1, ic ) - 1.0;

      // E_i = 1 + (K_2i - 1)v_2 + (K_3i - 1)v_3 (Equation 12)
      real64 const e = 1.0 + k2 * v2 + k3 * v3;

      // Objective function summation (Equation 10)
      F -= z * LvArray::math::log( e );

      // Gradients (Equation 13)
      real64 const term1 = z / e;
      g[0] -= term1 * k2;
      g[1] -= term1 * k3;

      // Hessian entries (symmetric 2x2 matrix computed explicitly)
      real64 const term2 = term1 / e;
      H[0] += term2 * k2 * k2;
      H[1] += term2 * k2 * k3;
      H[2] += term2 * k3 * k2;
      H[3] += term2 * k3 * k3;
    }
    return LvArray::math::sqrt( g[0]*g[0] + g[1]*g[1] );
  }

public:
  /**
   * @brief Solves the 3-phase Rachford-Rice equations allowing negative phase fractions.
   *
   * @details Implements the constrained optimization approach by Michelsen,
   * extended to the negative flash region by Leibovici and Nichita (2008).
   * Paper: "A new look at multiphase Rachford-Rice equations for negative flashes"
   * DOI: 10.1016/j.fluid.2008.03.006
   * The problem is formulated as an unbounded minimization with linear constraints
   * (hyperplanes E_i > 0) to ensure continuous differentiability at phase boundaries.
   * The inner loop uses a bounded Newton-Raphson scheme with Armijo backtracking line search.
   *
   * @param kValues Array slice of equilibrium constants (K-values).
   * @param feed Array slice of feed composition (mole fractions).
   * @param presentComponentIds Array slice of indices for components present in the mixture.
   * @param v2 Output mole fraction of the second phase. Serves as initial guess if feasible.
   * @param v3 Output mole fraction of the third phase. Serves as initial guess if feasible.
   * @return bool True if the minimization converges successfully within tolerances, false otherwise.
   */
  template< integer USD1, integer USD2 >
  GEOS_HOST_DEVICE
  static bool solve(
    arraySlice2d< real64 const, USD2 > const & kValues,
    arraySlice1d< real64 const, USD1 > const & feed,
    arraySlice1d< integer const > const & presentComponentIds,
    real64 & v2,
    real64 & v3 )
  {
    // Quick exit for trivial solutions
    // Checks if the mixture universally prefers a single phase
    // This avoids singularities in the Hessian matrix during Newton-Raphson.

    real64 maxK2 = -LvArray::NumericLimits< real64 >::max;
    real64 minK2 = LvArray::NumericLimits< real64 >::max;
    real64 maxK3 = -LvArray::NumericLimits< real64 >::max;
    real64 minK3 = LvArray::NumericLimits< real64 >::max;
    real64 maxK2OverK3 = -LvArray::NumericLimits< real64 >::max;
    real64 minK2OverK3 = LvArray::NumericLimits< real64 >::max;

    for( integer const & ic : presentComponentIds )
    {
      real64 const k2 = kValues( 0, ic );
      real64 const k3 = kValues( 1, ic );

      maxK2 = LvArray::math::max( maxK2, k2 );
      minK2 = LvArray::math::min( minK2, k2 );
      maxK3 = LvArray::math::max( maxK3, k3 );
      minK3 = LvArray::math::min( minK3, k3 );

      real64 const k2OverK3 = k2 / k3;
      maxK2OverK3 = LvArray::math::max( maxK2OverK3, k2OverK3 );
      minK2OverK3 = LvArray::math::min( minK2OverK3, k2OverK3 );
    }

    if( maxK2 < 1.0 && maxK3 < 1.0 )
    {
      v2 = 0.0;
      v3 = 0.0;
      return true;
    }
    if( minK2 > 1.0 && minK2OverK3 > 1.0 )
    {
      v2 = 1.0;
      v3 = 0.0;
      return true;
    }
    if( minK3 > 1.0 && maxK2OverK3 < 1.0 )
    {
      v2 = 0.0;
      v3 = 1.0;
      return true;
    }

    // Initial feasibility check
    // Evaluate if the provided initial conditions (v2, v3) are feasible (E_i > 0).
    bool feasible = true;
    for( integer const & ic : presentComponentIds )
    {
      real64 const k2 = kValues( 0, ic ) - 1.0;
      real64 const k3 = kValues( 1, ic ) - 1.0;
      real64 const e = 1.0 + k2 * v2 + k3 * v3;
      if( e <= 0.0 )
      {
        feasible = false;
        break;
      }
    }

    // Fallback inherently feasible initialization (1 / np) for 3 phases
    if( !feasible )
    {
      v2 = 1.0 / 3.0;
      v3 = 1.0 / 3.0;
    }

    // Bound-constrained Newton-Raphson Optimization
    real64 F = 0.0;
    real64 g[2]{};
    real64 H[4]{};

    for( integer iterations = 0; iterations < maxNewtonIterations; ++iterations )
    {
      real64 const error = computeObjective( kValues, feed, presentComponentIds, v2, v3, F, g, H );
      if( error < newtonTolerance )
      {
        return true;
      }

      // Add small regularization to strictly enforce positive definiteness
      // in cases where phases might artificially merge
      H[0] += epsilon;
      H[3] += epsilon;

      // Invert explicit 2x2 Hessian to compute Newton step (H * dv = -g)
      real64 const determinant = H[0] * H[3] - H[1] * H[2];
      if( determinant < epsilon )
      {
        return false;
      }

      real64 const dv2 = (-g[0] * H[3] + g[1] * H[1]) / determinant;
      real64 const dv3 = (-g[1] * H[0] + g[0] * H[2]) / determinant;

      // If the step size is too small we will declare convergence
      real64 const stepSize = LvArray::math::sqrt( dv2*dv2 + dv3*dv3 );
      if( stepSize < newtonTolerance )
      {
        return true;
      }

      // Line Search Step 1: Find distance to hyperplane boundaries.
      // The feasible domain is limited by hyperplanes E_i = 0 (Leibovici and Nichita 2008, Section 3).
      // We compute the maximum step size alphaBound that keeps E_i > 0 for all components.
      // The step along the search direction dv is restricted such that we do not cross the hyperplanes.
      // For each component, if the step decreases E_i (i.e., de < 0), we find the intersection.
      real64 alphaBound = 1.0;
      for( integer const & ic : presentComponentIds )
      {
        real64 const k2 = kValues( 0, ic );
        real64 const k3 = kValues( 1, ic );

        real64 const de = (k2 - 1.0) * dv2 + (k3 - 1.0) * dv3;
        real64 const e = 1.0 + (k2 - 1.0) * v2 + (k3 - 1.0) * v3;
        if( epsilon < LvArray::math::abs( de ) )
        {
          real64 const alphaI = -e / de;
          if( 0.0 < alphaI && alphaI < alphaBound )
          {
            // We scale the boundary by 1.0 - eps to prevent E_i from becoming exactly 0.
            alphaBound = threshold * alphaI;
          }
        }
      }

      // Limit the maximum step length to remain strictly inside the domain.
      real64 alpha = alphaBound;

      // Line Search Step 2: Backtracking (Armijo Rule).
      // Ensures sufficient decrease in the objective function to guarantee global convergence.
      // We start with the bounded step size 'alpha' and iteratively shrink it by half.
      // The Armijo condition checks if the new objective value FNew is sufficiently smaller
      // than the current value F plus a fraction (c1) of the directional derivative.
      real64 const c1 = 1.0e-4;
      real64 const directionalDerivative = g[0] * dv2 + g[1] * dv3;
      bool stepAccepted = false;

      for( integer lineSearchIteration = 0; lineSearchIteration < 20; ++lineSearchIteration )
      {
        real64 const v2New = v2 + alpha * dv2;
        real64 const v3New = v3 + alpha * dv3;

        real64 FNew = 0.0;

        for( integer const & ic : presentComponentIds )
        {
          real64 const k2 = kValues( 0, ic ) - 1.0;
          real64 const k3 = kValues( 1, ic ) - 1.0;
          real64 const e = 1.0 + k2 * v2New + k3 * v3New;
          FNew -= feed[ic] * LvArray::math::log( e );
        }

        // Armijo condition check: F(x + alpha * p) <= F(x) + c1 * alpha * grad(F)^T * p
        if( FNew <= F + c1 * alpha * directionalDerivative )
        {
          v2 = v2New;
          v3 = v3New;
          stepAccepted = true;
          break;
        }

        // Step too large; shrink
        alpha *= 0.5;
      }

      if( !stepAccepted )
      {
        // Step size became too small (likely converged to machine precision limit)
        return false;
      }
    }

    return false;
  }
};

} // namespace constitutive

} // namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_RACHFORDRICE_HPP_
