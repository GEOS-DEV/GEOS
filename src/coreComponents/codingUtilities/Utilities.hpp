/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
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

inline constexpr bool isEven(int x) {
  return !(x&1);
}

inline constexpr bool isOdd(int x) {
  return (x&1);
}

template< typename T1, typename T2, typename SORTED >
T2& stlMapLookup( mapBase<T1,T2, SORTED>& Map, const T1& key )
{
  typename mapBase<T1, T2, SORTED>::iterator MapIter = Map.find( key );
  GEOSX_ERROR_IF(MapIter==Map.end(), "Key not found: " << key);
  return MapIter->second;
}


template< typename T1, typename T2, typename SORTED >
const T2& stlMapLookup( const mapBase<T1,T2, SORTED>& Map, const T1& key)
{
  return (stlMapLookup( const_cast<mapBase<T1,T2, SORTED>&>(Map), key ));
}

template< typename T1, typename T2, typename SORTED, typename LAMBDA >
bool executeOnMapValue( mapBase<T1,T2, SORTED> const & Map, const T1& key, LAMBDA&& lambda )
{
  bool rval = false;
  typename mapBase<T1,T2, SORTED>::const_iterator MapIter = Map.find( key );
  if( MapIter!=Map.end() )
  {
    rval = true;
    lambda(MapIter->second);
  }

  return rval;
}

/**
 * @param val1
 * @param val2
 * @param tolfac
 * @return
 */
inline bool isEqual( const realT& val1, const realT& val2, const realT& tolfac=0.0 )
{
  realT tol = 0.0;
  if( tolfac > 1.0e-15 )
    tol = fabs(tolfac) * (fabs(val1)+fabs(val2))*0.5;
  return val1<=(val2+tol) && val1>=(val2-tol);
}

inline bool isZero( const realT& val, const realT& tol=std::numeric_limits<realT>::epsilon() )
{
  if( val<=tol && val>=-tol )
  {
    return true;
  }
  else
  {
    return false;
  }
}

template< typename T_KEY, typename T_VALUE, typename SORTED >
T_VALUE softMapLookup( mapBase<T_KEY,T_VALUE, SORTED> const & theMap,
                       T_KEY const & key,
                       T_VALUE const failValue )
{
  T_VALUE rvalue;
  typename mapBase<T_KEY,T_VALUE, SORTED>::const_iterator iter = theMap.find(key);
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

template<typename VEC1, typename VEC2>
inline void copy( localIndex N, VEC1 const & v1, VEC2 const & v2 )
{
  for (localIndex i = 0; i < N; ++i)
    v2[i] = v1[i];
}

template<typename MATRIX, typename VEC1, typename VEC2>
inline void applyChainRule( localIndex N,
                            MATRIX const & dy_dx,
                            VEC1 const & df_dy,
                            VEC2 const & df_dx )
{
  // this could use some dense linear algebra
  for (localIndex i = 0; i < N; ++i)
  {
    df_dx[i] = 0.0;
    for (localIndex j = 0; j < N; ++j)
    {
      df_dx[i] += df_dy[j] * dy_dx[j][i];
    }
  }
}

template<typename MATRIX, typename VEC1, typename VEC2>
inline void applyChainRuleInPlace( localIndex N,
                                   MATRIX const & dy_dx,
                                   VEC1 const & df_dxy,
                                   VEC2 & work )
{
  applyChainRule( N, dy_dx, df_dxy, work );
  copy( N, work, df_dxy );
}

} // namespace geosx

#endif /* GEOSX_CODINGUTILITIES_UTILITIES_H_ */
