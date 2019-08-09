/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * @file Utilities.hpp
 */

#ifndef UTILITIES_H_
#define UTILITIES_H_

#include "common/DataTypes.hpp"

namespace geosx
{

template< typename T >
inline void CopyGlobalToLocal(arraySlice1d<localIndex const> const & globalToLocalRelation,
                              arraySlice1d< T const> const& globalField1,
                              arraySlice1d< T const> const& globalField2,
                              arraySlice1d< T > & localField1,
                              arraySlice1d< T > & localField2,
                              localIndex N)
{
  for( localIndex a=0 ; a<N ; ++a )
  {
    localField1[a] = globalField1[ globalToLocalRelation[a] ];
    localField2[a] = globalField2[ globalToLocalRelation[a] ];
  }
}

template< localIndex N, typename T >
inline void CopyGlobalToLocal(arraySlice1d<localIndex const > const & globalToLocalRelation,
                              arraySlice1d< T const > const & globalField1,
                              arraySlice1d< T const > const & globalField2,
                              T * const restrict localField1,
                              T * const restrict localField2 )
{
  for( localIndex a=0 ; a<N ; ++a )
  {
    localField1[a] = globalField1[ globalToLocalRelation[a] ];
    localField2[a] = globalField2[ globalToLocalRelation[a] ];
  }
}

template< localIndex N, typename T >
inline void CopyGlobalToLocal(arraySlice1d<localIndex> const & globalToLocalRelation,
                              arraySlice1d< T > const & globalField1,
                              arraySlice1d< T > const & globalField2,
                              arraySlice1d< T > const & globalField3,
                              T * const restrict localField1,
                              T * const restrict localField2,
                              T * const restrict localField3 )
{
  for( localIndex a=0 ; a<N ; ++a )
  {
    localField1[a] = globalField1[ globalToLocalRelation[a] ];
    localField2[a] = globalField2[ globalToLocalRelation[a] ];
    localField3[a] = globalField3[ globalToLocalRelation[a] ];
  }
}

template< int N, typename T >
inline void CopyGlobalToLocal(arraySlice1d<localIndex const> const & globalToLocalRelation,
                              arraySlice1d< T const > const & globalField1,
                              arraySlice1d< T const > const & globalField2,
                              arraySlice1d< T const > const & globalField3,
                              arraySlice1d< T const > const & globalField4,
                              T * const localField1,
                              T * const localField2,
                              T * const localField3,
                              T * const localField4 )
{
  for( localIndex a=0 ; a<N ; ++a )
  {
    localField1[a] = globalField1[ globalToLocalRelation[a] ];
    localField2[a] = globalField2[ globalToLocalRelation[a] ];
    localField3[a] = globalField3[ globalToLocalRelation[a] ];
    localField4[a] = globalField4[ globalToLocalRelation[a] ];
  }
}

inline bool isEven(int x) {
  return !(x&1);
}

inline bool isOdd(int x) {
  return (x&1);
}

template< typename T1, typename T2, typename SORTED >
T2& stlMapLookup( mapBase<T1,T2, SORTED>& Map, const T1& key )
{
  typename mapBase<T1, T2, SORTED>::iterator MapIter = Map.find( key );
  GEOS_ERROR_IF(MapIter==Map.end(), "Key not found: " << key);
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


inline realT Power(const realT val, const realT exponent)
{
  if(isEqual(exponent, 0.5))
    return sqrt(val);
  else if(isEqual(exponent, 1.5))
    return val*sqrt(val);
  else if(isEqual(exponent, 2.0))
    return val*val;
  else
    return pow(val, exponent);
}

/**
 * @return cpu usage
 *
 * This function uses the rusage structure to query elapsed system time and user
 * time, and returns
 * the result.
 */
inline realT getcputime(void)
{
  struct timeval tim;
  struct rusage ru;
  getrusage(RUSAGE_SELF, &ru);

  tim=ru.ru_utime;
  realT t= tim.tv_sec + tim.tv_usec / 1.0e6;

  tim=ru.ru_stime;
  t+= tim.tv_sec + tim.tv_usec / 1.0e6;
  return t;
}

template< typename TYPE >
inline void Intersection( const set<TYPE>& set1, const set<TYPE>& set2, set<TYPE>& intersection )
{
  intersection.clear();
  typename set<TYPE>::const_iterator iter_1 = set1.begin();
  typename set<TYPE>::const_iterator iter_2 = set2.begin();

  while( iter_1!=set1.end() && iter_2!=set2.end() )
  {
    if( *iter_1 == *iter_2 )
    {
      intersection.insert(*iter_1);
      ++iter_1;
      ++iter_2;
    }
    else if( (*iter_1)<(*iter_2) )
    {
      ++iter_1;
    }
    else if( (*iter_1)>(*iter_2) )
    {
      ++iter_2;
    }
  }
}

template< typename TYPE >
inline void Intersection( const set<TYPE>& input, const array1d<TYPE>& arr, set<TYPE>& intersection )
{
  intersection.clear();

  for( typename array1d<TYPE>::const_iterator iter_arr=arr.begin() ; iter_arr!=arr.end() ; ++iter_arr )
  {
    if( input.count( *iter_arr ) == 1 )
    {
      intersection.insert(*iter_arr);
    }
  }
}

template< typename TYPE >
inline void Intersection( const set<TYPE>& input, const array1d<TYPE>& arr, array1d<TYPE>& intersection )
{
  intersection.clear();

  for( typename set<TYPE>::const_iterator iter_arr=arr.begin() ; iter_arr!=arr.end() ; ++iter_arr )
  {
    if( input.count( *iter_arr ) == 1 )
    {
      intersection.push_back(*iter_arr);
    }
  }
}

inline
real64_array linspace(realT start, realT stop, int count){
  real64_array rv(count);
  rv = start;

  if(count > 1)
  {

    realT dL = (stop-start)/double(count-1);
    realT sum = start;
    for(int i = 1 ; i < count-1 ; ++i)
    {
      sum += dL;
      rv[i] =  sum;
    }

    rv[count -1] = stop;
  }
  return rv;
}


inline
real64_array logspace(realT start, realT stop, int count){
  real64_array rv(count);
  rv = start;

  if(count > 1)
  {

    realT dL = std::pow(stop/start,1.0/double(count-1));
    realT prod = start;
    for(int i = 1 ; i < count-1 ; ++i)
    {
      prod *= dL;
      rv[i] =  prod;
    }

    rv[count -1] = stop;
  }
  return rv;
}



inline double GetOrder( double number, const unsigned int digits = 1 )
{

  const int exp = static_cast<int>(std::log10(number));
  const double magnitude = std::pow(10,exp);

  return ( digits>0 ? round( ( number / magnitude ) * std::pow(10,digits-1) ) / std::pow(10,digits-1) : 1 ) * magnitude;
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

#endif /* UTILITIES_H_ */
