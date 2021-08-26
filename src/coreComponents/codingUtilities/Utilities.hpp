/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file Utilities.hpp
 */

#ifndef GEOSX_CODINGUTILITIES_UTILITIES_H_
#define GEOSX_CODINGUTILITIES_UTILITIES_H_

#include "common/DataTypes.hpp"
#include "LvArray/src/limits.hpp"

namespace geosx
{

/**
 * @brief Compare two real values with a tolerance.
 * @tparam type of real value
 * @param val1 first value
 * @param val2 second value
 * @param relTol relative tolerance for comparison
 * @return @p true if @p val1 and @p val2 are within tolerance of each other
 */
template< typename T >
GEOSX_FORCE_INLINE GEOSX_HOST_DEVICE constexpr
bool isEqual( T const val1, T const val2, T const relTol = 0.0 )
{
  T const absTol = ( relTol > 10 * LvArray::NumericLimits< T >::epsilon ) ? relTol * ( fabs( val1 ) + fabs( val2 ) ) * 0.5 : 0.0;
  return ( val2 - absTol ) <= val1 && val1 <= ( val2 + absTol );
}

/**
 * @brief Test if a real value is (almost) zero.
 * @tparam T type of real value
 * @param val the value to test
 * @param tol absolute tolerance for comparison
 * @return @p true if @p val is within @p tol of zero
 */
template< typename T >
GEOSX_FORCE_INLINE GEOSX_HOST_DEVICE constexpr
bool isZero( T const val, T const tol = LvArray::NumericLimits< T >::epsilon )
{
  return -tol <= val && val <= tol;
}

template< typename T >
GEOSX_FORCE_INLINE GEOSX_HOST_DEVICE constexpr
bool isOdd( T x )
{
  static_assert( std::is_integral< T >::value, "Not meaningful for non-integral types" );
  return (x & 1);
}

template< typename T >
GEOSX_FORCE_INLINE GEOSX_HOST_DEVICE constexpr
bool isEven( T x )
{
  return !isOdd( x );
}

template< typename T1, typename T2, typename SORTED >
T2 & stlMapLookup( mapBase< T1, T2, SORTED > & Map, const T1 & key )
{
  typename mapBase< T1, T2, SORTED >::iterator MapIter = Map.find( key );
  GEOSX_ERROR_IF( MapIter==Map.end(), "Key not found: " << key );
  return MapIter->second;
}


template< typename T1, typename T2, typename SORTED >
const T2 & stlMapLookup( const mapBase< T1, T2, SORTED > & Map, const T1 & key )
{
  return (stlMapLookup( const_cast< mapBase< T1, T2, SORTED > & >(Map), key ));
}

template< typename T1, typename T2, typename SORTED, typename LAMBDA >
bool executeOnMapValue( mapBase< T1, T2, SORTED > const & Map, const T1 & key, LAMBDA && lambda )
{
  bool rval = false;
  typename mapBase< T1, T2, SORTED >::const_iterator MapIter = Map.find( key );
  if( MapIter!=Map.end() )
  {
    rval = true;
    lambda( MapIter->second );
  }

  return rval;
}

template< typename T_KEY, typename T_VALUE, typename SORTED >
T_VALUE softMapLookup( mapBase< T_KEY, T_VALUE, SORTED > const & theMap,
                       T_KEY const & key,
                       T_VALUE const failValue )
{
  T_VALUE rvalue;
  typename mapBase< T_KEY, T_VALUE, SORTED >::const_iterator iter = theMap.find( key );
  if( iter==theMap.end() )
  {
    rvalue = failValue;
  }
  else
  {
    rvalue = iter->second;
  }
  return rvalue;
}

// The code below should work with any subscriptable vector/matrix types

template< typename VEC1, typename VEC2 >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void copy( localIndex const N, VEC1 const & v1, VEC2 const & v2 )
{
  for( localIndex i = 0; i < N; ++i )
  {
    v2[i] = v1[i];
  }
}

template< typename MATRIX, typename VEC1, typename VEC2 >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void applyChainRule( localIndex const N,
                     MATRIX const & dy_dx,
                     VEC1 const & df_dy,
                     VEC2 && df_dx )
{
  // this could use some dense linear algebra
  for( localIndex i = 0; i < N; ++i )
  {
    df_dx[i] = 0.0;
    for( localIndex j = 0; j < N; ++j )
    {
      df_dx[i] += df_dy[j] * dy_dx[j][i];
    }
  }
}

template< typename MATRIX, typename VEC1, typename VEC2 >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void applyChainRuleInPlace( localIndex const N,
                            MATRIX const & dy_dx,
                            VEC1 && df_dxy,
                            VEC2 && work )
{
  applyChainRule( N, dy_dx, df_dxy, work );
  copy( N, work, df_dxy );
}

// TODO: move to the interpolation folder

/**
 * @brief Perform linear extrapolation of function f(x) at x3 using the known values of f on interval between x1 and x2
 * @tparam T the value type
 * @param[in] x1 first bound of the interval
 * @param[in] x2 second bound of the interval
 * @param[in] f1 value of function f(x) evaluated at x1
 * @param[in] f2 value of function f(x) evaluated at x2
 * @param[in] x3 value at which we want to extrapolate linearly function f(x)
 * @return the value of function f(x) extrapolated at x3
 */
template< typename T >
T linearExtrapolation( T const x1,
                       T const x2,
                       T const f1,
                       T const f2,
                       T const x3 )
{
  return ( f2 - f1 ) / ( x2 - x1 ) * ( x3 - x2 ) + f2;
}

/**
 * @brief Perform log extrapolation of function f(x) at x3 using the known values of f on interval between x1 and x2
 * @tparam T the value type
 * @param[in] x1 first bound of the interval
 * @param[in] x2 second bound of the interval
 * @param[in] f1 value of function f(x) evaluated at x1
 * @param[in] f2 value of function f(x) evaluated at x2
 * @param[in] x3 value at which we want to extrapolate logarithmically function f(x)
 * @return the value of function f(x) extrapolated at x3
 */
template< typename T >
T logExtrapolation( T x1,
                    T x2,
                    T f1,
                    T f2,
                    T x3 )
{
  T const lnf3 = ( log( f2 ) - log( f1 ) ) / ( log( x2 ) - log( x1 ) ) * ( log( x3 ) - log( x2 ) ) + log( f2 );
  return exp( lnf3 );
}

/**
 * @brief Perform linear interpolation of function f(x) at x3 on interval [x1,x2]
 * @tparam T the value type
 * @param[in] dx1 length of interval [x1,x3]: x3-x1
 * @param[in] dx2 length of interval [x3,x2]: x2-x3
 * @param[in] f1 value of function f(x) evaluated at x1
 * @param[in] f2 value of function f(x) evaluated at x2
 * @return the interpolated value of function f(x) at x3
 */
template< typename T >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
T linearInterpolation( T const dx1,
                       T const dx2,
                       T const f1,
                       T const f2 )
{
  T const alpha = dx1 / ( dx1 + dx2 );
  return alpha * f2 + ( 1.0 - alpha ) * f1;
}

/**
 * @brief Perform linear interpolation of function f(x) at x3 on interval [x1,x2]
 * @tparam T the value type
 * @param[in] dx1 length of interval [x1,x3]: x3-x1
 * @param[in] dx2 length of interval [x3,x2]: x2-x3
 * @param[in] f1 value of function f(x) evaluated at x1
 * @param[in] f2 value of function f(x) evaluated at x2
 * @param[out] f the interpolated value of function f(x) at x3
 * @param[out] df_dx the derivative (wrt x) of the interpolated value of function f(x) at x3
 */
template< typename T >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void linearInterpolation( T const dx1,
                          T const dx2,
                          T const f1,
                          T const f2,
                          T & f,
                          T & df_dx )
{
  T const alpha = dx1 / ( dx1 + dx2 );
  T const dAlpha_dx = 1.0 / ( dx1 + dx2 );
  f = alpha * f2 + ( 1.0 - alpha ) * f1;
  df_dx = dAlpha_dx * ( f2 - f1 );
}

/**
 * @brief Perform linear interpolation of function f(x) for all the values of a given array
 * @tparam T the value type
 * @param[in] xIn the values at which we know the value of f(x)
 * @param[in] fIn the known values of function f(x) at the values in xIn
 * @param[in] xOut the values where we want to interpolate f(x)
 * @param[out] fOut the interpolated values of f(x) at xOut
 */
template< typename T >
void linearInterpolation( arrayView1d< T const > const & xIn,
                          arrayView1d< T const > const & fIn,
                          arrayView1d< T const > const & xOut,
                          arrayView1d< T > const & fOut )
{
  GEOSX_ASSERT_EQ_MSG( xIn.size(), fIn.size(), "linearInterpolation: size mismatch in the xIn and fIn input vectors" );
  GEOSX_ASSERT_EQ_MSG( xOut.size(), fOut.size(), "linearInterpolation: size mismatch in the xOut and fOut input vectors" );

  for( localIndex i = 0; i < xOut.size(); ++i )
  {
    integer const idx = LvArray::sortedArrayManipulation::find( xIn.begin(),
                                                                xIn.size(),
                                                                xOut[i] );
    integer const iUp  = LvArray::math::min( LvArray::math::max( idx, 1 ),
                                             LvArray::integerConversion< integer >( xIn.size()-1 ) );
    integer const iLow = iUp-1;
    fOut[i] = linearInterpolation( xOut[i] - xIn[iLow], xIn[iUp] - xOut[i], fIn[iLow], fIn[iUp] );
  }
}

} // namespace geosx

#endif /* GEOSX_CODINGUTILITIES_UTILITIES_H_ */
