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

namespace geosx
{

/**
 * @brief GPU-friendly analogue of std::numeric_limits
 * @tparam T type of numeric value
 */
template< typename T >
struct NumericTraits
{
  static constexpr T min = std::numeric_limits< T >::lowest();
  static constexpr T max = std::numeric_limits< T >::max();
  static constexpr T eps = std::numeric_limits< T >::epsilon();
};

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
  T const absTol = ( relTol > 10 * NumericTraits< T >::eps ) ? relTol * ( fabs( val1 ) + fabs( val2 ) ) * 0.5 : 0.0;
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
bool isZero( T const val, T const tol = NumericTraits< T >::eps )
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
T2 & stlMapLookup( mapBase< T1, T2, SORTED > & map, const T1 & key )
{
  typename mapBase< T1, T2, SORTED >::iterator MapIter = map.find( key );
  GEOSX_ERROR_IF( MapIter==map.end(), "Key not found: " << key );
  return MapIter->second;
}


template< typename T1, typename T2, typename SORTED >
const T2 & stlMapLookup( const mapBase< T1, T2, SORTED > & map, const T1 & key )
{
  return (stlMapLookup( const_cast< mapBase< T1, T2, SORTED > & >(map), key ));
}

template< typename T1, typename T2, typename SORTED, typename LAMBDA >
bool executeOnMapValue( mapBase< T1, T2, SORTED > const & map, const T1 & key, LAMBDA && lambda )
{
  bool rval = false;
  typename mapBase< T1, T2, SORTED >::const_iterator MapIter = map.find( key );
  if( MapIter!=map.end() )
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
void copy( localIndex const n, VEC1 const & v1, VEC2 const & v2 )
{
  for( localIndex i = 0; i < n; ++i )
  {
    v2[i] = v1[i];
  }
}

template< typename MATRIX, typename VEC1, typename VEC2 >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void applyChainRule( localIndex const n,
                     MATRIX const & dyDx,
                     VEC1 const & dfDy,
                     VEC2 && dfDx )
{
  // this could use some dense linear algebra
  for( localIndex i = 0; i < n; ++i )
  {
    dfDx[i] = 0.0;
    for( localIndex j = 0; j < n; ++j )
    {
      dfDx[i] += dfDy[j] * dyDx[j][i];
    }
  }
}

template< typename MATRIX, typename VEC1, typename VEC2 >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void applyChainRuleInPlace( localIndex const n,
                            MATRIX const & dyDx,
                            VEC1 && dfDxy,
                            VEC2 && work )
{
  applyChainRule( n, dyDx, dfDxy, work );
  copy( n, work, dfDxy );
}

} // namespace geosx

#endif /* GEOSX_CODINGUTILITIES_UTILITIES_H_ */
