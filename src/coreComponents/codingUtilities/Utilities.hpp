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

// TODO: make sure this is a good place
// TODO: check/fix the implementation of these functions

template< typename T >
void findSurroundingIndex( array1d< T > const & x,
                           T xval,
                           integer & iminus,
                           integer & iplus )
{
  GEOSX_THROW_IF( x.empty(), "Interpolation table is empty", InputError );
  GEOSX_THROW_IF( x[0] > xval, "Input x value is out of range, extrapolation not allowed", InputError );

  // search for interval
  for( iplus = 1; iplus < x.size() - 1 && x[iplus] < xval; ++iplus )
  { }
  iminus = iplus - 1;
}


// TODO: replace with something from LvArray
template< typename T >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void findSurroundingIndex( arrayView1d< T const > const & x,
                           T xval,
                           integer & iminus,
                           integer & iplus )
{
  // search for interval
  for( iplus = 1; iplus < x.size() - 1 && x[iplus] < xval; ++iplus )
  { }
  iminus = iplus - 1;
}

template< typename T >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void findSurroundingIndex( arraySlice1d< const T > const & x,
                           T xval,
                           integer & iminus,
                           integer & iplus )
{
  // search for interval
  for( iplus = 1; iplus < x.size() - 1 && x[iplus] < xval; ++iplus )
  { }
  iminus = iplus - 1;
}

template< typename T >
void interpolation1( array1d< T > const & xin,
                     array1d< T > const & yin,
                     array1d< T > const & xout,
                     array1d< T > & yout )
{
  GEOSX_THROW_IF( xin.size() != yin.size(), "interpolation1: size mismatch!", InputError );
  GEOSX_THROW_IF( xout.size() != yout.size(), "interpolation1: size mismatch!", InputError );
  for( localIndex n = 0; n != xout.size(); ++n )
  {
    T const x = xout[n];
    GEOSX_THROW_IF( xin[0] > x, "interpolation1: input x out of range, cannot extrapolate below the limits", InputError );

    // search for interval
    integer i_minus, i_plus;
    findSurroundingIndex( xin, x, i_minus, i_plus );

    if( i_minus == i_plus )
    {
      yout[n] = yin[i_minus];
    }
    else
    {
      yout[n] = yin[i_minus] * ( xin[i_plus] - x ) / ( xin[i_plus] - xin[i_minus] )
                + yin[i_plus] * ( ( x - xin[i_minus] ) / ( xin[i_plus] - xin[i_minus] ) );
      GEOSX_THROW_IF( yout[n] < 0.0, "interpolation1: negative value", InputError );
    }
  }
}

template< typename T >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
T linearInterpolation( T dminus,
                       T dplus,
                       T xminus,
                       T xplus )
{
  T const f = dminus / ( dminus + dplus );
  return f * xplus + ( 1 - f ) * xminus;
}

template< typename T >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void linearInterpolation( T dminus,
                          T dplus,
                          T xminus,
                          T xplus,
                          T & x,
                          T & dx )
{
  T const f = dminus / ( dminus + dplus );
  T const df = 1./( dminus + dplus );
  x = f * xplus + ( 1 - f ) * xminus;
  dx = df*(xplus - xminus);
}

template< typename T >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
T linearInterpolation( T x1,
                       T y1,
                       T x2,
                       T y2,
                       T x3 )
{
  return linearInterpolation( x3 - x1, x2 - x3, y1, y2 );
}

template< typename T >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void linearInterpolation( T x1,
                          T y1,
                          T x2,
                          T y2,
                          T x3,
                          T & x,
                          T & dx )
{
  linearInterpolation( x3 - x1, x2 - x3, y1, y2, x, dx );
}

template< typename T >
T linearExtrapolation( T x1,
                       T y1,
                       T x2,
                       T y2,
                       T x3 )
{
  return ( y2 - y1 ) / ( x2 - x1 ) * ( x3 - x2 ) + y2;
}

template< typename T >
T logExtrapolation( T x1,
                    T y1,
                    T x2,
                    T y2,
                    T x3 )
{
  T const lny = ( log( y2 ) - log( y1 ) ) / ( log( x2 ) - log( x1 ) ) * ( log( x3 ) - log( x2 ) ) + log( y2 );
  return exp( lny );
}


} // namespace geosx

#endif /* GEOSX_CODINGUTILITIES_UTILITIES_H_ */
