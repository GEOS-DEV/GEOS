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

#include "RachfordRice.hpp"

namespace geos
{

namespace constitutive
{

GEOS_HOST_DEVICE
real64
RachfordRice::solve( arraySlice1d< real64 const > const kValues,
                     arraySlice1d< real64 const > const feed,
                     arraySlice1d< integer const > const presentComponentIds )
{
  real64 gasPhaseMoleFraction = 0;

  // min and max Kvalues for non-zero composition
  real64 maxK = 0.0;
  real64 minK = 1 / epsilon;

  for( integer i = 0; i < presentComponentIds.size(); ++i )
  {
    integer const ic = presentComponentIds[i];
    if( kValues[ic] > maxK )
    {
      maxK = kValues[ic];
    }
    if( kValues[ic] < minK )
    {
      minK = kValues[ic];
    }
  }

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
    else if( ( funcXMin > 1 ) && ( funcXMax > 1 ) )
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

GEOS_HOST_DEVICE
real64
RachfordRice::evaluate( arraySlice1d< real64 const > const kValues,
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

GEOS_HOST_DEVICE
real64
RachfordRice::evaluateDerivative( arraySlice1d< real64 const > const kValues,
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

} // namespace constitutive

} // namespace geos
