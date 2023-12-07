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
  using Deriv = multifluid::DerivativeOffset;
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

  /**
   * @brief Function solving the Rachford-Rice equation
   * @input[in] kValues the array of K-values
   * @input[in] feed the component fractions
   * @input[in] presentComponentIds the indices of components with a non-zero fractions
   * @return the gas mole fraction
   **/
  GEOS_HOST_DEVICE
  real64
  static
  solve( arraySlice1d< real64 const > const kValues,
         arraySlice1d< real64 const > const feed,
         arraySlice1d< integer const > const presentComponentIds )
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
    real64 funcXMin = 0.0;
    real64 funcXMid = 0.0;
    real64 funcXMax = 0.0;
    bool recomputeMin = true;
    bool recomputeMax = true;
    integer SSIIteration = 0;

    while( ( currentError > SSITolerance ) && ( SSIIteration < maxSSIIterations ) )
    {
      real64 const xMid = 0.5 * ( xMin + xMax );
      if( recomputeMin )
      {
        funcXMin = evaluate( kValues, feed, presentComponentIds, xMin );
      }
      if( recomputeMax )
      {
        funcXMax = evaluate( kValues, feed, presentComponentIds, xMax );
      }
      funcXMid = evaluate( kValues, feed, presentComponentIds, xMid );
      if( ( funcXMin < 0 ) && ( funcXMax < 0 ) )
      {
        return gasPhaseMoleFraction = 0.0;
      }
      else if( ( funcXMin > 0 ) && ( funcXMax > 0 ) )
      {
        return gasPhaseMoleFraction = 1.0;
      }
      else if( funcXMin * funcXMid < 0.0 )
      {
        xMax = xMid;
        recomputeMax = true;
        recomputeMin = false;
      }
      else if( funcXMax * funcXMid < 0.0 )
      {
        xMin = xMid;
        recomputeMax = false;
        recomputeMin = true;
      }

      currentError = LvArray::math::min( LvArray::math::abs( funcXMax - funcXMin ),
                                         LvArray::math::abs( xMax - xMin ) );
      SSIIteration++;

      // TODO: add warning if max number of SSI iterations is reached
    }

    gasPhaseMoleFraction = 0.5 * ( xMax + xMin );

    // step 3: start the Newton loop
    integer newtonIteration = 0;
    real64 newtonValue = gasPhaseMoleFraction;

    while( ( currentError > newtonTolerance ) && ( newtonIteration < maxNewtonIterations ) )
    {
      real64 const deltaNewton = -evaluate( kValues, feed, presentComponentIds, newtonValue )
                                 / evaluateDerivative( kValues, feed, presentComponentIds, newtonValue );
      currentError = LvArray::math::abs( deltaNewton ) / LvArray::math::abs( newtonValue );

      // test if we are stepping out of the [xMin;xMax] interval
      if( newtonValue + deltaNewton < xMin )
      {
        newtonValue = 0.5 * ( newtonValue + xMin );
      }
      else if( newtonValue + deltaNewton > xMax )
      {
        newtonValue = 0.5 * ( newtonValue + xMax );
      }
      else
      {
        newtonValue = newtonValue + deltaNewton;
      }
      newtonIteration++;

      // TODO: add warning if max number of Newton iterations is reached
    }
    return gasPhaseMoleFraction = newtonValue;
  }

  /**
   * @brief Compute derivatives after solving the Rachford-Rice equation
   * @input[in] kValues the array of K-values
   * @input[in] kValueDerivs derivatives of the K-values
   * @input[in] feed the component fractions
   * @input[in] presentComponentIds the indices of components with a non-zero fractions
   * @input[in] vapourFraction the calculated vapour fraction
   * @input[out] vapourFractionDerivs derivatives of the vapour fraction
   * @return the gas mole fraction
   **/
  GEOS_HOST_DEVICE
  void
  static
  computeDerivatives( arraySlice1d< real64 const > const kValues,
                      arraySlice2d< real64 const > const kValueDerivs,
                      arraySlice1d< real64 const > const feed,
                      arraySlice1d< integer const > const presentComponentIds,
                      real64 const vapourFraction,
                      arraySlice1d< real64 > const vapourFractionDerivs )
  {
    integer const numComps = feed.size();
    LvArray::forValuesInSlice( vapourFractionDerivs, []( real64 & val ){ val = 0.0; } );

    // Look for trivial solutions
    real64 minK, maxK;
    findKValueRange( kValues, presentComponentIds, minK, maxK );
    if( maxK < 1.0 || 1.0 < minK )
    {
      return;
    }

    for( integer const kc : {Deriv::dP, Deriv::dT} )
    {
      real64 numerator = 0.0;
      real64 denominator = 0.0;
      for( integer const ic : presentComponentIds )
      {
        real64 const k = ( kValues[ic] - 1.0 );
        real64 const dk = 1.0 / ( 1.0 + vapourFraction * k );
        numerator += feed[ic] * kValueDerivs( ic, kc ) * dk * (1.0 - k * vapourFraction * dk);
        denominator += feed[ic] * k * k * dk * dk;
      }
      vapourFractionDerivs[kc] = numerator / denominator;
    }

    for( integer jc = 0; jc < numComps; ++jc )
    {
      integer const kc = Deriv::dC + jc;

      real64 numerator = 0.0;
      real64 denominator = 0.0;

      for( integer const ic : presentComponentIds )
      {
        real64 const k = ( kValues[ic] - 1.0 );
        real64 const dk = 1.0 / ( 1.0 + vapourFraction * k );
        numerator += feed[ic] * kValueDerivs( ic, kc ) * dk * (1.0 - k * vapourFraction * dk);
        denominator += feed[ic] * k * k * dk * dk;
      }
      real64 const k = ( kValues[jc] - 1.0 );
      real64 const dk = 1.0 / ( 1.0 + vapourFraction * k );
      numerator += k * dk;
      vapourFractionDerivs[kc] = numerator / denominator;
    }
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
  GEOS_HOST_DEVICE
  real64
  static
  evaluate( arraySlice1d< real64 const > const kValues,
            arraySlice1d< real64 const > const feed,
            arraySlice1d< integer const > const presentComponentIds,
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
  GEOS_HOST_DEVICE
  real64
  static
  evaluateDerivative( arraySlice1d< real64 const > const kValues,
                      arraySlice1d< real64 const > const feed,
                      arraySlice1d< integer const > const presentComponentIds,
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
  GEOS_FORCE_INLINE
  GEOS_HOST_DEVICE
  void
  static
  findKValueRange( arraySlice1d< real64 const > const kValues,
                   arraySlice1d< integer const > const presentComponentIds,
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

};

} // namespace constitutive

} // namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_RACHFORDRICE_HPP_
