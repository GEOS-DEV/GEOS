/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
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

/*
 * Utilities.h
 *
 *  Created on: Sep 13, 2010
 *      Author: settgast1
 */

#ifndef UTILITIES_H_
#define UTILITIES_H_

#include "common/DataTypes.hpp"
#include <limits>
#include <sys/resource.h>
#include <map>
#include <algorithm>
#include "RAJA/RAJA.hpp"
#include "RAJA/util/defines.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"


#ifdef USE_ATK
#include "slic/slic.hpp"
#endif

//#include "../legacy/Common/GPException.h"

/////////////////////////////////////////////////
// Forward declaration of templated functions

namespace geosx
{
template< class T >
void PushFieldForwardInTime( const realT& dt,
                             const array1d< T >& dfield,
                             array1d< T >& field );

template< class T >
void IntegrateFieldInTime( const realT& dt,
                           const array1d< T >& field,
                           array1d< T >& Ifield );

template< class T >
inline void IntegrateField( const realT& dt,
                            const T& dfdt,
                            T& df );

template< typename T >
void SetConstPointer( T* const& pointer,  T*  newpointer );

template < typename T1, typename T2 >
void ClearStlMapValues( std::map<T1,T2>& Map );

template <class ElementClass, class SetClass>
bool isMember(const ElementClass& x, const SetClass& aSetOrMap);

template< typename T1, typename T2 >
const T2* stlMapLookupPointer( const std::map<T1,T2>& Map, const T1& key );

template< typename T1, typename T2 >
T2* stlMapLookupPointer( std::map<T1,T2>& Map, const T1& key );


template< typename T1, typename T2 >
const T2& stlMapLookup( const std::map<T1,T2>& Map, const T1& key, const std::string& message="" );

template< typename T1, typename T2 >
T2& stlMapLookup( std::map<T1,T2>& Map, const T1& key, const std::string& message="" );


real64_array logspace(realT start, realT stop, int count=100);
real64_array linspace(realT start, realT stop, int count=100);

//#pragma GCC diagnostic push
//#pragma GCC diagnostic ignored "-Wsign-compare"
//
//
//template< typename RTYPE, typename T >
//typename std::enable_if< std::is_unsigned<T>::value && std::is_signed<RTYPE>::value, RTYPE >::type
//integer_conversion( T input )
//{
//  static_assert( std::numeric_limits<T>::is_integer, "input is not an integer type" );
//  static_assert( std::numeric_limits<RTYPE>::is_integer, "requested conversion is not an integer type" );
//
//  if( input > std::numeric_limits<RTYPE>::max()  )
//  {
//    abort();
//  }
//  return static_cast<RTYPE>(input);
//}
//
//template< typename RTYPE, typename T >
//typename std::enable_if< std::is_signed<T>::value && std::is_unsigned<RTYPE>::value, RTYPE >::type
//integer_conversion( T input )
//{
//  static_assert( std::numeric_limits<T>::is_integer, "input is not an integer type" );
//  static_assert( std::numeric_limits<RTYPE>::is_integer, "requested conversion is not an integer type" );
//
//  if( input > std::numeric_limits<RTYPE>::max() ||
//      input < 0 )
//  {
//    abort();
//  }
//  return static_cast<RTYPE>(input);
//}
//
//
//template< typename RTYPE, typename T >
//typename std::enable_if< ( std::is_signed<T>::value && std::is_signed<RTYPE>::value ) ||
//                         ( std::is_unsigned<T>::value && std::is_unsigned<RTYPE>::value ), RTYPE >::type
//integer_conversion( T input )
//{
//  static_assert( std::numeric_limits<T>::is_integer, "input is not an integer type" );
//  static_assert( std::numeric_limits<RTYPE>::is_integer, "requested conversion is not an integer type" );
//
//  if( input > std::numeric_limits<RTYPE>::max() ||
//      input < std::numeric_limits<RTYPE>::lowest() )
//  {
//    abort();
//  }
//  return static_cast<RTYPE>(input);
//}
//
//
//
//#pragma GCC diagnostic pop


/////////////////////////////////////////////////



template< class T >
inline void PushFieldForwardInTime( const realT& dt,
                                    const array1d< T >& dfield,
                                    array1d< T >& field )
{
  T dfieldDt;

  const int N = field.size();
  for( int a=0 ; a<N ; ++a )
  {
    dfieldDt = dfield(a);
    dfieldDt *= dt;

    field(a) += dfieldDt;
  }
}


template< class T >
inline void IntegrateFieldInTime( const realT& dt,
                                  const array1d< T >& field,
                                  array1d< T >& Ifield )
{
  const int N = field.size();
  for( int a=0 ; a<N ; ++a )
  {
    Ifield(a) = field(a);
    Ifield(a) *= dt;
  }
}


template< class T >
inline void IntegrateField( const realT& dt,
                            const T& dfdt,
                            T& df )
{
  df = dfdt;
  df *= dt;
}



template< typename T >
inline void CopyGlobalToLocal(const localIndex_array& globalToLocalRelation,
                              const array1d< T >& globalField,
                              array1d< T >& localField)
{
  const typename array1d<T>::size_type N = globalToLocalRelation.size();

  for( typename array1d<T>::size_type a=0 ; a<N ; ++a )
  {
    localField[a] = globalField[ globalToLocalRelation[a] ];
  }
}



template< typename T >
inline void CopyGlobalToLocal(const localIndex* __restrict__ const globalToLocalRelation,
                              const array1d< T >& globalField,
                              array1d< T >& localField)
{
  const typename array1d<T>::size_type N = localField.size();

  for( typename array1d<T>::size_type a=0 ; a<N ; ++a )
  {
    localField[a] = globalField[ globalToLocalRelation[a] ];
  }
}

template< typename T >
inline void CopyGlobalToLocal(const localIndex* __restrict__ const globalToLocalRelation,
                              const array1d< T >& globalField1,
                              const array1d< T >& globalField2,
                              array1d< T >& localField1,
                              array1d< T >& localField2 )
{
  const typename array1d<T>::size_type N = localField1.size();

  for( typename array1d<T>::size_type a=0 ; a<N ; ++a )
  {
    localField1[a] = globalField1[ globalToLocalRelation[a] ];
    localField2[a] = globalField2[ globalToLocalRelation[a] ];
  }
}

template< typename T >
inline void CopyGlobalToLocal(const localIndex* __restrict__ const globalToLocalRelation,
                              const array1d< T >& globalField1,
                              const array1d< T >& globalField2,
                              T * __restrict__ const localField1,
                              T * __restrict__ const localField2,
                              localIndex N)
{
  for( localIndex a=0 ; a<N ; ++a )
  {
    localField1[a] = globalField1[ globalToLocalRelation[a] ];
    localField2[a] = globalField2[ globalToLocalRelation[a] ];
  }
}

template< typename T >
inline void CopyGlobalToLocal(const localIndex* __restrict__ const globalToLocalRelation,
                              T const * __restrict__ const globalField1,
                              T const * __restrict__ const globalField2,
                              T * __restrict__ const localField1,
                              T * __restrict__ const localField2,
                              localIndex N)
{
  for( localIndex a=0 ; a<N ; ++a )
  {
    localField1[a] = globalField1[ globalToLocalRelation[a] ];
    localField2[a] = globalField2[ globalToLocalRelation[a] ];
  }
}

template< typename T >
inline void CopyGlobalToLocal(const localIndex* __restrict__ const globalToLocalRelation,
                              const array1d< T >& globalField1,
                              const array1d< T >& globalField2,
                              const array1d< T >& globalField3,
                              array1d< T >& localField1,
                              array1d< T >& localField2,
                              array1d< T >& localField3 )
{
  const typename array1d<T>::size_type N = localField1.size();

  for( typename array1d<T>::size_type a=0 ; a<N ; ++a )
  {
    localField1[a] = globalField1[ globalToLocalRelation[a] ];
    localField2[a] = globalField2[ globalToLocalRelation[a] ];
    localField3[a] = globalField3[ globalToLocalRelation[a] ];
  }
}

template< typename T >
inline void CopyGlobalToLocal(const localIndex* __restrict__ const globalToLocalRelation,
                              const array1d< T >& globalField1,
                              const array1d< T >& globalField2,
                              const array1d< T >& globalField3,
                              const array1d< T >& globalField4,
                              array1d< T >& localField1,
                              array1d< T >& localField2,
                              array1d< T >& localField3,
                              array1d< T >& localField4 )
{
  const typename array1d<T>::size_type N = localField1.size();

  for( typename array1d<T>::size_type a=0 ; a<N ; ++a )
  {
    localField1[a] = globalField1[ globalToLocalRelation[a] ];
    localField2[a] = globalField2[ globalToLocalRelation[a] ];
    localField3[a] = globalField3[ globalToLocalRelation[a] ];
    localField4[a] = globalField4[ globalToLocalRelation[a] ];
  }
}


template< typename T >
inline void CopyGlobalToLocal(const localIndex* __restrict__ const globalToLocalRelation,
                              T const * __restrict__ const globalField1,
                              T const * __restrict__ const globalField2,
                              T const * __restrict__ const globalField3,
                              T const * __restrict__ const globalField4,
                              T * __restrict__ const localField1,
                              T * __restrict__ const localField2,
                              T * __restrict__ const localField3,
                              T * __restrict__ const localField4,
                              localIndex N)
{
  for( localIndex a=0 ; a<N ; ++a )
  {
    localField1[a] = globalField1[ globalToLocalRelation[a] ];
    localField2[a] = globalField2[ globalToLocalRelation[a] ];
    localField3[a] = globalField3[ globalToLocalRelation[a] ];
    localField4[a] = globalField4[ globalToLocalRelation[a] ];
  }
}

template< typename T >
inline void AddLocalToGlobal( const localIndex* __restrict__ const globalToLocalRelation,
                              const array1d< T >& localField,
                              array1d< T >& globalField)
{
  const typename array1d<T>::size_type N = localField.size();

  for( typename array1d<T>::size_type a=0 ; a<N ; ++a )
  {
    globalField[ globalToLocalRelation[a] ] += localField[a];
  }
}

template< typename T >
inline void AddLocalToGlobal( const localIndex* __restrict__ const globalToLocalRelation,
                              T const * __restrict__ const localField,
                              array1d< T >& globalField,
                              localIndex const N )
{
  for( localIndex a=0 ; a<N ; ++a )
  {
    globalField[ globalToLocalRelation[a] ] += localField[a];
  }
}


template<typename T, typename atomicPol=atomicPolicy>
inline void AddLocalToGlobal( const localIndex* __restrict__ const globalToLocalRelation,
                                    T const * __restrict__ const localField,
                                    T * __restrict__ const globalField,
                                    localIndex const N )
{
  for( typename array1d<T>::size_type a=0 ; a<N ; ++a )
  {
    geosx::raja::atomicAdd<atomicPol>(&globalField[globalToLocalRelation[a]],localField[a]);
  }
}

template<typename atomicPol=atomicPolicy>
inline void AddLocalToGlobal( const localIndex * __restrict__ const globalToLocalRelation,
                                    R1Tensor const * __restrict__ const localField,
                                    R1Tensor * __restrict__ const globalField,
                                    localIndex const N )
{

  for( localIndex a=0 ; a<N ; ++a )
  {
    real64 * __restrict__ const gData = globalField[globalToLocalRelation[a]].Data();
    real64 const * __restrict__ const lData = localField[a].Data();
    for( localIndex i=0 ; i<3 ; ++i )
    {
      geosx::raja::atomicAdd<atomicPol>( &gData[i], lData[i] );
    }
  }
}

//01-22-2018 - Hack, we will have to fix. 
#ifdef USE_OPENMP
template<>
inline void AddLocalToGlobal<R1Tensor,RAJA::atomic::omp_atomic>( const localIndex* __restrict__ const globalToLocalRelation,
                                         R1Tensor const * __restrict__ const localField,
                                         R1Tensor * __restrict__ const globalField,
                                         localIndex const N )
{

  for( typename array1d<R1Tensor>::size_type a=0 ; a<N ; ++a )
    {
      double * const lhs = globalField[ globalToLocalRelation[a] ].Data();
      double const * const rhs = localField[a].Data();
      for( int i=0; i<3; ++i )
        {          
          geosx::raja::atomicAdd<RAJA::atomic::omp_atomic>(&lhs[i],rhs[i]);
        }
    }
}
#endif

template< typename T >
inline void AddLocalToGlobal( const localIndex* __restrict__ const globalToLocalRelation,
                              const array1d< T >& localField1,
                              const array1d< T >& localField2,
                              array1d< T >& globalField1,
                              array1d< T >& globalField2 )
{
  const typename array1d<T>::size_type  N = localField1.size();

  for( typename array1d<T>::size_type a=0 ; a<N ; ++a )
  {
    globalField1[ globalToLocalRelation[a] ] += localField1[a];
    globalField2[ globalToLocalRelation[a] ] += localField2[a];
  }
}



template< typename T >
inline bool listsHaveEqualPermutations( const T* const list1, const T* const list2, const localIndex n )
{
  localIndex count = 0;
  bool rval = false;
  for( localIndex i1=0 ; i1<n ; ++i1 )
    for( localIndex i2=0 ; i2<n ; ++i2 )
    {
      if( list1[i1] == list2[i2] )
      {
        ++count;
        break;
      }
    }

  if( count == n )
    rval = true;

  return rval;

}

inline bool isEven(int x) {
  return !(x&1);
}

inline bool isOdd(int x) {
  return (x&1);
}


/// find if object is member of vector
template <class ElementClass>
bool isMember(const ElementClass& x, const std::vector<ElementClass>& aVec) {return ( std::find(aVec.begin(), aVec.end(), x) !=  aVec.end() );}
inline bool isMember(localIndex x, const localIndex_array& aVec) {return ( std::find(aVec.begin(), aVec.end(), x) !=  aVec.end() );}


/// find if object is member of set or map
template <class ElementClass, class SetClass>
bool isMember(const ElementClass& x, const SetClass& aSetOrMap) {return ( aSetOrMap.find(x) != aSetOrMap.end() );}



/// permutation tensor
inline int eijk(int i,int j,int k){
  return ((i-j)*(j-k)*(k-i))/2;
}


template< typename T >
void SetConstPointer( T* const& pointer,  T*  newpointer )
{
  T** temp = const_cast< T**>(&pointer);
  *temp = newpointer;

  return;
}



template< typename T1, typename T2 >
T2* stlMapLookupPointer( std::map<T1,T2>& Map, const T1& key )
{
  T2* rval = NULL;
  typename std::map<T1,T2>::iterator MapIter = Map.find( key );
  if( MapIter!=Map.end()  )
  {
    rval = &(MapIter->second);
  }


  return rval;
}


template< typename T1, typename T2 >
const T2* stlMapLookupPointer( const std::map<T1,T2>& Map, const T1& key )
{
  const T2* rval = NULL;
  typename std::map<T1,T2>::const_iterator MapIter = Map.find( key );
  if( MapIter!=Map.end()  )
  {
    rval = &(MapIter->second);
  }


  return rval;
}


template< typename T1, typename T2 >
T2& stlMapLookup( std::map<T1,T2>& Map, const T1& key, const std::string& message )
{
  typename std::map<T1,T2>::iterator MapIter = Map.find( key );
  if( MapIter==Map.end()  )
  {
    std::cout<<std::endl;
    std::stringstream st;
    st << "Error in stlMapLookup. Key not found in map! key: " << key << " message: " << message <<"\n";
//    throw GPException(st.str().c_str());
#ifdef USE_ATK
    SLIC_ERROR(st.str());
#endif

  }

  return MapIter->second;
}


template< typename T1, typename T2 >
const T2& stlMapLookup( const std::map<T1,T2>& Map, const T1& key, const std::string& message )
{
  return (stlMapLookup( const_cast<std::map<T1,T2>&>(Map), key, message ));
}


/*
 * Initialize map values.
 * std::map<int,int> aMap = CreateStlMap<int, int >(1,2)(3,4)(5,6);
 */
template <typename K, typename V>
class CreateStlMap
{
private:
  std::map<K, V> m_map;
public:
  CreateStlMap(const K& key, const V& value){ m_map[key] = value; }

  CreateStlMap<K, V>& operator()(const K& key, const V& value)
  {
    m_map[key] = value;
    return *this;
  }

  operator std::map<K, V>() { return m_map; }
};


template < typename T1, typename T2 >
void ClearStlMapValues( std::map<T1,T2>& Map )
{
  for( typename std::map<T1,T2>::iterator iter=Map.begin() ; iter!=Map.end() ; ++iter )
  {
    iter->second.clear();
  }
}

/*
 * Initialize vector values.
 * Usage:
 *   std::vector<int> aVect;
 *   aVect += 1,1,2,3,4;
 */
template <class T> class vector_inserter
{
public:
  std::vector<T>& v;
  vector_inserter(std::vector<T>& vv): v(vv){}
  vector_inserter& operator,(const T& val){v.push_back(val); return *this;}
};
template <class T> vector_inserter<T>& operator+=(std::vector<T>& v,const T& x);

template <class T> vector_inserter<T>& operator+=(std::vector<T>& v,const T& x){
  return vector_inserter<T>(v),x;
}



/*
 * Initialize string vector values with chars
 * Usage:
 *   string_array aVect;
 *   aVect += "The","quick","brown","fox";
 */
class svector_inserter
{
public:
  string_array& v;
  svector_inserter(string_array& vv): v(vv){}
  svector_inserter& operator,(const char* val){v.push_back(std::string(val)); return *this;}
};
svector_inserter& operator+=(string_array& v,const std::string& x);

inline
svector_inserter& operator+=(string_array& v,const char* x){
  return svector_inserter(v),x;
}



/// Static Assert
/// Check satement at compile time
template <bool b>
struct gp_static_assert {};
// Specialization with member function
template <>
struct gp_static_assert<true>
{
  static void is_valid() {}
};
// use: gp_static_assert<TEST>is_valid();



inline bool isGTE0( const int i )
{ return (i>=0); }



/**
 * @author Randy Settgast
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
 * @author Randy Settgast
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


template< typename T_KEY, typename T_VALUE >
T_VALUE softMapLookup( map<T_KEY,T_VALUE> const & theMap,
                       T_KEY const & key,
                       T_VALUE const failValue )
{
  T_VALUE rvalue;
  typename map<T_KEY,T_VALUE>::const_iterator iter = theMap.find(key);
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

}

#endif /* UTILITIES_H_ */
