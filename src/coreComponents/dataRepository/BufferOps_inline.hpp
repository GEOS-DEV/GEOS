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



#ifndef DATAREPOSITORY_BUFFEROPS_INLINE_H_
#define DATAREPOSITORY_BUFFEROPS_INLINE_H_

#include "common/DataTypes.hpp"
#include "codingUtilities/Utilities.hpp"
#include "codingUtilities/static_if.hpp"
#include "codingUtilities/GeosxTraits.hpp"
#include "IntegerConversion.hpp"

#include <type_traits>

namespace geosx
{

namespace bufferOps
{

/**
 *
 * @param buffer
 * @param var
 * @return
 */
template< bool DO_PACKING, typename T >
typename std::enable_if< std::is_trivial<T>::value, localIndex >::type
Pack( char*&  buffer, T const & var )
{
  localIndex const sizeOfPackedChars = sizeof(T);
  static_if( DO_PACKING )
  {
    memcpy( buffer, &var, sizeOfPackedChars );
    buffer += sizeOfPackedChars;
  }
  end_static_if
  return sizeOfPackedChars;
}

/**
 * @param buffer
 * @param var
 * @return size (in bytes) of unpacked data
 */
template< typename T>
typename std::enable_if< std::is_trivial<T>::value, localIndex >::type
Unpack( char const *& buffer, T & var )
{
  localIndex const sizeOfUnpackedChars = sizeof(T);
  memcpy( &var, buffer, sizeOfUnpackedChars );
  buffer += sizeOfUnpackedChars;
  return sizeOfUnpackedChars;
}


template< bool DO_PACKING, typename T, typename INDEX_TYPE >
typename std::enable_if< std::is_trivial<T>::value, localIndex >::type
Pack( char*&  buffer, T const * const restrict var, INDEX_TYPE const length )
{
  localIndex sizeOfPackedChars = Pack<DO_PACKING>( buffer, length );

  sizeOfPackedChars += length * sizeof(T);
  static_if( DO_PACKING )
  {
    memcpy( buffer, var, length * sizeof(T) );
    buffer += length * sizeof(T);
  }
  end_static_if

  return sizeOfPackedChars;
}


template< typename T, typename INDEX_TYPE >
typename std::enable_if< std::is_trivial<T>::value, localIndex >::type
Unpack( char const *& buffer, T * const restrict var, INDEX_TYPE const expectedLength )
{
  INDEX_TYPE length;
  localIndex sizeOfUnpackedChars = Unpack( buffer, length );

  GEOS_ASSERT_MSG( length == expectedLength, "expectedLength != length: " <<
                   expectedLength << " != " << length );

  memcpy( var, buffer, length * sizeof(T) );
  sizeOfUnpackedChars += length * sizeof(T);
  buffer += length * sizeof(T);

  return sizeOfUnpackedChars;
}


template< bool DO_PACKING >
localIndex Pack( char*& buffer, const std::string& var )
{
  const string::size_type varSize = var.size();
  localIndex sizeOfPackedChars = Pack<DO_PACKING>( buffer, varSize );

  static_if( DO_PACKING )
  {
    memcpy( buffer, var.data(), varSize );
    buffer += varSize;
  }
  end_static_if
  
  sizeOfPackedChars += varSize;
  return sizeOfPackedChars;
}


inline localIndex Unpack( char const *& buffer, string& var )
{
  string::size_type stringsize = 0;
  localIndex sizeOfUnpackedChars = Unpack( buffer, stringsize );

  var.resize(stringsize);
  var.assign( buffer, stringsize );
  buffer += stringsize;
  sizeOfUnpackedChars += stringsize;

  return sizeOfUnpackedChars;
}


template< bool DO_PACKING, typename T >
typename std::enable_if< traits::is_tensorT<T>::value, localIndex >::type
Pack( char*& buffer, T const & var )
{
  localIndex sizeOfPackedChars = 0;
  sizeOfPackedChars += Pack<DO_PACKING>( buffer, var.Data(), var.Length());
  return sizeOfPackedChars;
}


template< typename T >
typename std::enable_if< traits::is_tensorT<T>::value, localIndex >::type
Unpack( char const *& buffer, T & var )
{
  localIndex sizeOfUnpackedChars = 0;
  real64 * const pVar = var.Data();
  int const length = var.Length();
  sizeOfUnpackedChars += Unpack( buffer, pVar, length );
  return sizeOfUnpackedChars;
}


template< bool DO_PACKING, typename T, typename INDEX_TYPE >
typename std::enable_if< !std::is_trivial<T>::value, localIndex >::type
Pack( char*&  buffer, T const * const restrict var, INDEX_TYPE const length )
{
  localIndex sizeOfPackedChars = Pack<DO_PACKING>( buffer, length );

  for( INDEX_TYPE a=0 ; a<length ; ++a )
  {
    sizeOfPackedChars += Pack<DO_PACKING>( buffer, var[ a ] );
  }

  return sizeOfPackedChars;
}


template< typename T, typename INDEX_TYPE >
typename std::enable_if< !std::is_trivial<T>::value, localIndex >::type
Unpack( char const *& buffer, T * const restrict var, INDEX_TYPE const expectedLength )
{
  INDEX_TYPE length;
  localIndex sizeOfUnpackedChars = Unpack( buffer, length );

  GEOS_ASSERT_MSG( length == expectedLength, "expectedLength != length: " <<
                   expectedLength << " != " << length );

  for( INDEX_TYPE a=0 ; a<length ; ++a )
  {
    sizeOfUnpackedChars += Unpack( buffer, var[a] );
  }

  return sizeOfUnpackedChars;
}


template< bool DO_PACKING, typename T, typename INDEX_TYPE >
localIndex
Pack( char*&  buffer,
      T const * const restrict var,
      arraySlice1d<INDEX_TYPE const> const & indices,
      INDEX_TYPE const length )
{
  localIndex sizeOfPackedChars = Pack<DO_PACKING>( buffer, length );

  for( INDEX_TYPE a=0 ; a<length ; ++a )
  {
    sizeOfPackedChars += Pack<DO_PACKING>( buffer, var[ indices[a] ] );
  }

  return sizeOfPackedChars;
}


template< typename T, typename INDEX_TYPE >
localIndex
Unpack( char const *& buffer,
        T * const restrict var,
        arraySlice1d<INDEX_TYPE const> const & indices,
        INDEX_TYPE & length )
{
  localIndex sizeOfUnpackedChars = Unpack( buffer, length );

  for( INDEX_TYPE a=0 ; a<length ; ++a )
  {
    sizeOfUnpackedChars += Unpack( buffer, var[ indices[a] ] );
  }

  return sizeOfUnpackedChars;
}

#ifdef GEOSX_USE_ARRAY_BOUNDS_CHECK


template< bool DO_PACKING, typename T, typename INDEX_TYPE >
typename std::enable_if< std::is_trivial<T>::value, localIndex >::type
Pack( char*& buffer,
      arraySlice1d<T> const & var,
      INDEX_TYPE const length )
{
  localIndex sizeOfPackedChars = Pack<DO_PACKING>( buffer, length );
  sizeOfPackedChars += length * sizeof(T);

  static_if( DO_PACKING )
  {
    T * const restrict buffer_T = static_cast< T * >( buffer );
    for ( INDEX_TYPE i = 0; i < length; ++i )
    {
      buffer_T[ i ] = var[ i ];
    }

    buffer += length * sizeof(T);
  }
  end_static_if

  return sizeOfPackedChars;
}


template< bool DO_PACKING, typename T, typename INDEX_TYPE >
typename std::enable_if< !std::is_trivial<T>::value, localIndex >::type
Pack( char*& buffer,
      arraySlice1d<T> const & var,
      INDEX_TYPE const length )
{
  localIndex sizeOfPackedChars = Pack<DO_PACKING>( buffer, length );

  for( INDEX_TYPE a=0 ; a<length ; ++a )
  {
    sizeOfPackedChars += Pack<DO_PACKING>( buffer, var[ a ] );
  }

  return sizeOfPackedChars;
}


template< typename T, typename INDEX_TYPE >
typename std::enable_if< std::is_trivial<T>::value, localIndex >::type
Unpack( char const *& buffer,
        arraySlice1d<T> & var,
        INDEX_TYPE const expectedLength )
{
  INDEX_TYPE length;
  localIndex sizeOfUnpackedChars = Unpack( buffer, length );

  GEOS_ASSERT_MSG( length == expectedLength, "expectedLength != length: " <<
                   expectedLength << " != " << length );

  T const * const restrict buffer_T = reinterpret_cast< T const * >( buffer );
  for ( INDEX_TYPE i = 0; i < length; ++i )
  {
    var[ i ] = buffer_T[ i ];
  }

  buffer += length * sizeof(T);

  return sizeOfUnpackedChars;
}


template< typename T, typename INDEX_TYPE >
typename std::enable_if< !std::is_trivial<T>::value, localIndex >::type
Unpack( char const *& buffer,
        arraySlice1d<T> & var,
        INDEX_TYPE const expectedLength )
{
  INDEX_TYPE length;
  localIndex sizeOfUnpackedChars = Unpack( buffer, length );

  GEOS_ASSERT_MSG( length == expectedLength, "expectedLength != length: " <<
                   expectedLength << " != " << length );

  for( INDEX_TYPE a=0 ; a<length ; ++a )
  {
    sizeOfUnpackedChars += Unpack( buffer, var[a] );
  }

  return sizeOfUnpackedChars;
}


template< bool DO_PACKING, typename T, typename INDEX_TYPE >
localIndex
Pack( char*&  buffer,
      arraySlice1d<T> const & var,
      arraySlice1d<INDEX_TYPE> const & indices,
      INDEX_TYPE const length )
{
  localIndex sizeOfPackedChars = Pack<DO_PACKING>( buffer, length );

  for( INDEX_TYPE a=0 ; a<length ; ++a )
  {
    sizeOfPackedChars += Pack<DO_PACKING>( buffer, var[ indices[a] ] );
  }

  return sizeOfPackedChars;
}


template< typename T, typename INDEX_TYPE >
localIndex
Unpack( char const *& buffer,
        arraySlice1d<T> & var,
        arraySlice1d<INDEX_TYPE> const & indices,
        INDEX_TYPE & length )
{
  localIndex sizeOfUnpackedChars = 0;

  sizeOfUnpackedChars += Unpack( buffer, length );

  for( INDEX_TYPE a=0 ; a<length ; ++a )
  {
    sizeOfUnpackedChars += Unpack( buffer, var[ indices[a] ] );
  }

  return sizeOfUnpackedChars;
}

#endif /* GEOSX_USE_ARRAY_BOUNDS_CHECK */


template< bool DO_PACKING, typename T >
localIndex Pack( char *& buffer, set<T> const & var )
{
  const localIndex length = integer_conversion<localIndex>(var.size());
  localIndex sizeOfPackedChars = Pack<DO_PACKING>( buffer, length );;

  for( typename set<T>::const_iterator i=var.begin() ; i!=var.end() ; ++i )
  {
    sizeOfPackedChars += Pack<DO_PACKING>(buffer, *i);
  }

  return sizeOfPackedChars;
}


template< typename T>
localIndex Unpack( char const *& buffer, set<T> & var )
{
  var.clear();

  localIndex set_length;
  localIndex sizeOfUnpackedChars = Unpack( buffer, set_length );

  for( localIndex a=0 ; a<set_length ; ++a )
  {
    T temp;
    sizeOfUnpackedChars += Unpack( buffer, temp );
    var.insert( temp );
  }

  return sizeOfUnpackedChars;
}

template< bool DO_PACKING >
localIndex Pack( char *& buffer,
                 set<localIndex> const & var,
                 set<globalIndex> const & unmappedGlobalIndices,
                 arraySlice1d<globalIndex const> const & localToGlobal )
{
  const localIndex length = integer_conversion<localIndex>(var.size()+unmappedGlobalIndices.size());
  localIndex sizeOfPackedChars = Pack<DO_PACKING>( buffer, length );

  for( typename set<localIndex>::const_iterator i=var.begin() ; i!=var.end() ; ++i )
  {
    sizeOfPackedChars += Pack<DO_PACKING>( buffer, localToGlobal[*i]);
  }

  for( typename set<globalIndex>::const_iterator i=unmappedGlobalIndices.begin() ;
      i!=unmappedGlobalIndices.end() ; ++i )
  {
    sizeOfPackedChars += Pack<DO_PACKING>( buffer, *i);
  }


  return sizeOfPackedChars;
}

inline
localIndex Unpack( char const *& buffer,
                   set<localIndex> & var,
                   set<globalIndex> & unmappedGlobalIndices,
                   map<globalIndex,localIndex> const & globalToLocalMap,
                   bool const clearFlag )
{
  if( clearFlag )
  {
    var.clear();
  }
  localIndex set_length;
  localIndex sizeOfUnpackedChars = Unpack( buffer, set_length );

  for( localIndex a=0 ; a<set_length ; ++a )
  {
    globalIndex temp;
    sizeOfUnpackedChars += Unpack( buffer, temp );
    map<globalIndex,localIndex>::const_iterator iter = globalToLocalMap.find(temp);
    if( iter==globalToLocalMap.end() )
    {
      unmappedGlobalIndices.insert(temp);
    }
    else
    {
      var.insert( iter->second );
    }
  }

  return sizeOfUnpackedChars;
}




template< bool DO_PACKING >
localIndex Pack( char *& buffer,
                 set<localIndex> const & var,
                 arrayView1d<localIndex const> const & packList,
                 set<globalIndex> const & unmappedGlobalIndices,
                 arraySlice1d<globalIndex const> const & localToGlobal )
{

  localIndex length = 0;
  array1d<localIndex> temp(var.size());

  for( auto a : packList )
  {
    if( var.count(a) )
    {
      temp[length] = a;
      ++length;
    }
  }
  temp.resize(length);
  localIndex sizeOfPackedChars = Pack<DO_PACKING>( buffer, length );


  for( localIndex a=0 ; a<length ; ++a )
  {
    sizeOfPackedChars += Pack<DO_PACKING>( buffer, localToGlobal[temp[a]]);
  }

  for( typename set<globalIndex>::const_iterator i=unmappedGlobalIndices.begin() ;
      i!=unmappedGlobalIndices.end() ; ++i )
  {
    sizeOfPackedChars += Pack<DO_PACKING>( buffer, *i);
  }

  return sizeOfPackedChars;
}








template< bool DO_PACKING, typename T, int NDIM, typename INDEX_TYPE >
typename std::enable_if< bufferOps::is_packable_array< LvArray::ArrayView<T,NDIM,INDEX_TYPE> >::value, localIndex >::type
Pack( char*& buffer,
      LvArray::ArrayView<T,NDIM,INDEX_TYPE> const & var )
{
  localIndex sizeOfPackedChars = Pack<DO_PACKING>( buffer, var.dims(), NDIM );
  sizeOfPackedChars += Pack<DO_PACKING>( buffer, var.strides(), NDIM );

  const localIndex length = var.size();
  localIndex sizeOfPackedArrayChars = length * sizeof(T);

  sizeOfPackedChars += Pack<DO_PACKING>( buffer, var.data(), length );

  return sizeOfPackedChars;
}


template< typename T, int NDIM, typename INDEX_TYPE >
typename std::enable_if< bufferOps::is_packable_array< LvArray::Array<T,NDIM,INDEX_TYPE> >::value, localIndex >::type
Unpack( char const *& buffer, LvArray::Array<T,NDIM,INDEX_TYPE> & var )
{
  INDEX_TYPE dims[NDIM];
  localIndex sizeOfUnpackedChars = Unpack( buffer, dims, NDIM );

  var.resize( NDIM, dims );

  INDEX_TYPE strides[NDIM];
  sizeOfUnpackedChars += Unpack( buffer, strides, NDIM );

  INDEX_TYPE const * const existingStrides = var.strides();
  for( int i=0 ; i<NDIM ; ++i )
  {
    GEOS_ASSERT_MSG( strides[i] == existingStrides[i], "Strides are inconsistent: " <<
                     strides[i] << " != " << existingStrides[i] );
  }

  sizeOfUnpackedChars += Unpack( buffer, var.data(), var.size() );

  return sizeOfUnpackedChars;
}


template< bool DO_PACKING, typename T, int NDIM, typename T_indices, typename INDEX_TYPE >
typename std::enable_if< bufferOps::is_packable_array< LvArray::ArrayView<T,NDIM,INDEX_TYPE> >::value, localIndex >::type
Pack( char*& buffer,
      LvArray::ArrayView<T,NDIM,INDEX_TYPE> const & var,
      const T_indices& indices )
{
  localIndex sizeOfPackedChars = Pack<DO_PACKING>( buffer, var.strides(), NDIM );

  for( localIndex a=0 ; a<indices.size() ; ++a )
  {
    sizeOfPackedChars += Pack<DO_PACKING>( buffer, var.data(indices[a]), var.strides()[0] );
  }
  return sizeOfPackedChars;
}


template< typename T, int NDIM, typename T_indices, typename INDEX_TYPE >
localIndex
Unpack( char const *& buffer,
        LvArray::ArrayView<T,NDIM,INDEX_TYPE> & var,
        const T_indices& indices )
{
  INDEX_TYPE strides[NDIM];
  localIndex sizeOfUnpackedChars = Unpack( buffer, strides, NDIM );

  for( int dim=0 ; dim<NDIM ; ++dim )
  {
    GEOS_ASSERT_MSG( strides[dim] == var.strides()[dim], "Strides are inconsistent: " <<
                     strides[dim] << " != " << var.strides()[dim] );
  }

  for( localIndex a=0 ; a<indices.size() ; ++a )
  {
    sizeOfUnpackedChars += Unpack( buffer, var.data( indices[a] ), strides[0] );
  }

  return sizeOfUnpackedChars;
}


template< bool DO_PACKING >
localIndex
Pack( char*& buffer,
      arraySlice1d<localIndex const> const & var,
      arraySlice1d<globalIndex const> const & unmappedGlobalIndices,
      localIndex const length,
      arraySlice1d<globalIndex const> const & localToGlobalMap )
{
  localIndex sizeOfPackedChars = Pack<DO_PACKING>( buffer, length );
  sizeOfPackedChars += length*sizeof(globalIndex);

  static_if( DO_PACKING )
  {
    globalIndex * const buffer_GI = reinterpret_cast<globalIndex * const>(buffer);
    for( localIndex a=0 ; a<length ; ++a )
    {
      if( var[a] != unmappedLocalIndexValue )
      {
        buffer_GI[ a ] = localToGlobalMap[var[a]];
      }
      else
      {
        buffer_GI[ a ] = unmappedGlobalIndices[a];
      }
    }

    buffer += length * sizeof(globalIndex);
  }
  end_static_if

  return sizeOfPackedChars;
}


inline
localIndex
Unpack( char const *& buffer,
        localIndex_array & var,
        array1d<globalIndex> & unmappedGlobalIndices,
        map<globalIndex,localIndex> const & globalToLocalMap )
{
  localIndex length;
  localIndex sizeOfUnpackedChars = Unpack( buffer, length );
  var.resize(length);
  unmappedGlobalIndices.resize(length);
  unmappedGlobalIndices = unmappedLocalIndexValue;

  bool unpackedGlobalFlag = false;
  for( localIndex a=0 ; a<length ; ++a )
  {
    globalIndex unpackedGlobalIndex;
    sizeOfUnpackedChars += Unpack( buffer, unpackedGlobalIndex );

    map<globalIndex,localIndex>::const_iterator
    iter = globalToLocalMap.find(unpackedGlobalIndex);
    if( iter == globalToLocalMap.end() )
    {
      var[a] = unmappedLocalIndexValue;
      unmappedGlobalIndices[a] = unpackedGlobalIndex;
      unpackedGlobalFlag = true;
    }
    else
    {
      var[a] = iter->second;
    }
  }
  if( !unpackedGlobalFlag )
  {
    unmappedGlobalIndices.clear();
  }

  return sizeOfUnpackedChars;
}


inline
localIndex
Unpack( char const *& buffer,
        arraySlice1d<localIndex> & var,
        array1d<globalIndex> & unmappedGlobalIndices,
        localIndex const expectedLength,
        map<globalIndex,localIndex> const & globalToLocalMap )
{
  localIndex sizeOfUnpackedChars = 0;

  localIndex length;
  sizeOfUnpackedChars += Unpack( buffer, length );

  GEOS_ASSERT_MSG( length == expectedLength, "expectedLength != length: " <<
                   expectedLength << " != " << length );

  unmappedGlobalIndices.resize(length);
  unmappedGlobalIndices = unmappedLocalIndexValue;

  bool unpackedGlobalFlag = false;
  for( localIndex a=0 ; a<length ; ++a )
  {
    globalIndex unpackedGlobalIndex;
    sizeOfUnpackedChars += Unpack( buffer, unpackedGlobalIndex );

    map<globalIndex,localIndex>::const_iterator
    iter = globalToLocalMap.find(unpackedGlobalIndex);
    if( iter == globalToLocalMap.end() )
    {
      var[a] = unmappedLocalIndexValue;
      unmappedGlobalIndices[a] = unpackedGlobalIndex;
      unpackedGlobalFlag = true;
    }
    else
    {
      var[a] = iter->second;
    }
  }
  if( !unpackedGlobalFlag )
  {
    unmappedGlobalIndices.clear();
  }

  return sizeOfUnpackedChars;
}


template< bool DO_PACKING >
localIndex
Pack( char*& buffer,
      arrayView1d<localIndex_array const> const & var,
      map< localIndex, array1d<globalIndex> > const & unmappedGlobalIndices,
      arrayView1d<localIndex const> const & indices,
      arrayView1d<globalIndex const> const & localToGlobalMap,
      arrayView1d<globalIndex const> const & relatedObjectLocalToGlobalMap )
{
  localIndex sizeOfPackedChars=0;

  sizeOfPackedChars += Pack<DO_PACKING>( buffer, indices.size() );
  for( localIndex a=0 ; a<indices.size() ; ++a )
  {
    localIndex const li = indices[a];
    sizeOfPackedChars += Pack<DO_PACKING>( buffer, localToGlobalMap[li] );

    map< localIndex, array1d<globalIndex> >::const_iterator
    iterUnmappedGI = unmappedGlobalIndices.find(li);

    array1d<globalIndex> junk;
    array1d<globalIndex> const & unmappedGI = iterUnmappedGI==unmappedGlobalIndices.end() ?
                                              junk :
                                              iterUnmappedGI->second;


    sizeOfPackedChars += Pack<DO_PACKING>( buffer,
                                           var[li],
                                           unmappedGI,
                                           var[li].size(),
                                           relatedObjectLocalToGlobalMap );
  }

  return sizeOfPackedChars;
}


inline
localIndex
Unpack( char const *& buffer,
        arrayView1d<localIndex_array> & var,
        array1d<localIndex> & indices,
        map<localIndex, array1d<globalIndex> > & unmappedGlobalIndices,
        map<globalIndex,localIndex> const & globalToLocalMap,
        map<globalIndex,localIndex> const & relatedObjectGlobalToLocalMap )
{
  localIndex numIndicesUnpacked;
  localIndex const sizeOfIndicesPassedIn = indices.size();

  localIndex sizeOfUnpackedChars = Unpack( buffer, numIndicesUnpacked );

  GEOS_ERROR_IF( sizeOfIndicesPassedIn!=0 && numIndicesUnpacked!=indices.size(),
                 "number of unpacked indices("<<numIndicesUnpacked<<") does not equal size of "
                 "indices passed into Unpack function("<<sizeOfIndicesPassedIn );

  indices.resize(numIndicesUnpacked);

  for( localIndex a=0 ; a<indices.size() ; ++a )
  {
    globalIndex gi;
    sizeOfUnpackedChars += Unpack( buffer, gi );

    localIndex & li = indices[a];
    if( sizeOfIndicesPassedIn > 0 )
    {
      GEOS_ERROR_IF( li!=globalToLocalMap.at(gi),
                     "global index "<<gi<<" unpacked from buffer does not equal the lookup "
                     <<li<<" for localIndex "<<li<<" on this rank");
    }
    else
    {
      li = globalToLocalMap.at(gi);
    }

    array1d<globalIndex> unmappedIndices;
    sizeOfUnpackedChars += Unpack( buffer,
                                   var[li],
                                   unmappedIndices,
                                   relatedObjectGlobalToLocalMap );
    if( unmappedIndices.size() > 0 )
    {
      unmappedGlobalIndices[li] = unmappedIndices;
    }
  }
  return sizeOfUnpackedChars;
}


template< bool DO_PACKING >
localIndex
Pack( char*& buffer,
      arrayView1d< set<localIndex> const > const & var,
      map< localIndex, set<globalIndex> > const & unmappedGlobalIndices,
      arrayView1d<localIndex const> const & indices,
      arrayView1d<globalIndex const> const & localToGlobalMap,
      arrayView1d<globalIndex const> const & relatedObjectLocalToGlobalMap )
{
  localIndex sizeOfPackedChars = Pack<DO_PACKING>( buffer, indices.size() );

  for( localIndex a=0 ; a<indices.size() ; ++a )
  {
    localIndex li = indices[a];
    sizeOfPackedChars += Pack<DO_PACKING>( buffer, localToGlobalMap[li] );

    map< localIndex, set<globalIndex> >::const_iterator
    iterUnmappedGI = unmappedGlobalIndices.find(li);

    set<globalIndex> junk;
    set<globalIndex> const & unmappedGI = iterUnmappedGI==unmappedGlobalIndices.end() ?
                                          junk :
                                          iterUnmappedGI->second;

    sizeOfPackedChars += Pack<DO_PACKING>( buffer,
                                           var[li],
                                           unmappedGI,
                                           relatedObjectLocalToGlobalMap );
  }

  return sizeOfPackedChars;
}


inline
localIndex
Unpack( char const *& buffer,
        arrayView1d< set<localIndex> > & var,
        localIndex_array & indices,
        map<localIndex, set<globalIndex> > & unmappedGlobalIndices,
        map<globalIndex,localIndex> const & globalToLocalMap,
        map<globalIndex,localIndex> const & relatedObjectGlobalToLocalMap,
        bool const clearFlag )
{
  localIndex sizeOfUnpackedChars=0;
  localIndex const sizeOfIndicesPassedIn = indices.size();

  localIndex numIndicesUnpacked;
  sizeOfUnpackedChars += Unpack( buffer, numIndicesUnpacked );
  GEOS_ERROR_IF( sizeOfIndicesPassedIn!=0 && numIndicesUnpacked!=indices.size(),
                 "number of unpacked indices("<<numIndicesUnpacked<<") does not equal size of "
                 "indices passed into Unpack function("<<sizeOfIndicesPassedIn );

  indices.resize(numIndicesUnpacked);

  for( localIndex a=0 ; a<numIndicesUnpacked ; ++a )
  {

    globalIndex gi;
    sizeOfUnpackedChars += Unpack( buffer, gi );

    localIndex & li = indices[a];
    if( sizeOfIndicesPassedIn > 0 )
    {
      GEOS_ERROR_IF( li!=globalToLocalMap.at(gi),
                     "global index "<<gi<<" unpacked from buffer does equal the lookup "
                     <<li<<" for localIndex "<<li<<" on this rank");
    }
    else
    {
      li = globalToLocalMap.at(gi);
    }

    set<globalIndex> unmappedIndices;
    sizeOfUnpackedChars += Unpack( buffer,
                                   var[li],
                                   unmappedIndices,
                                   relatedObjectGlobalToLocalMap,
                                   clearFlag );

    unmappedGlobalIndices[li].insert( unmappedIndices.begin(),unmappedIndices.end() );
  }
  return sizeOfUnpackedChars;
}


template< bool DO_PACKING >
localIndex
Pack( char*& buffer,
      arrayView2d<localIndex> const & var,
      arrayView1d<localIndex> const & indices,
      arraySlice1d<globalIndex> const & localToGlobalMap )
{
  localIndex sizeOfPackedChars = Pack<DO_PACKING>( buffer, indices.size() );
  for( localIndex a=0 ; a<indices.size() ; ++a )
  {
    localIndex li = indices[a];
    sizeOfPackedChars += Pack<DO_PACKING>( buffer, localToGlobalMap[li] );

    sizeOfPackedChars += Pack<DO_PACKING>( buffer, var[li], var.size(1) );
  }

  return sizeOfPackedChars;
}


inline
localIndex
Unpack( char const *& buffer,
        arrayView2d<localIndex> & var,
        array1d<localIndex> & indices,
        map<globalIndex,localIndex> const & globalToLocalMap )
{
  localIndex sizeOfUnpackedChars = 0;
  localIndex const sizeOfIndicesPassedIn = indices.size();

  localIndex numIndicesUnpacked;
  sizeOfUnpackedChars += Unpack( buffer, numIndicesUnpacked );

  GEOS_ERROR_IF( sizeOfIndicesPassedIn!=0 && numIndicesUnpacked!=indices.size(),
                 "number of unpacked indices("<<numIndicesUnpacked<<") does not equal size of "
                 "indices passed into Unpack function("<<sizeOfIndicesPassedIn );

  indices.resize(numIndicesUnpacked);

  for( localIndex a=0 ; a<numIndicesUnpacked ; ++a )
  {
    globalIndex gi;
    sizeOfUnpackedChars += Unpack( buffer, gi );

    localIndex & li = indices[a];
    if( sizeOfIndicesPassedIn > 0 )
    {
      GEOS_ERROR_IF( li!=globalToLocalMap.at(gi),
                     "global index "<<gi<<" unpacked from buffer does equal the lookup "
                     <<li<<" for localIndex "<<li<<" on this rank");
    }
    else
    {
      li = globalToLocalMap.at(gi);
    }

    localIndex * const varSlice = var[li];
    sizeOfUnpackedChars += Unpack( buffer, varSlice, var.size(1) );
  }


  return sizeOfUnpackedChars;
}


template< bool DO_PACKING >
localIndex
Pack( char*& buffer,
      arrayView2d<localIndex const> const & var,
      map< localIndex, array1d<globalIndex> > const & unmappedGlobalIndices,
      arrayView1d<localIndex const> const & indices,
      arraySlice1d<globalIndex const> const & localToGlobalMap,
      arraySlice1d<globalIndex const> const & relatedObjectLocalToGlobalMap )
{
  localIndex sizeOfPackedChars = 0;

  sizeOfPackedChars += Pack<DO_PACKING>( buffer, indices.size() );
  for( localIndex a=0 ; a<indices.size() ; ++a )
  {
    localIndex li = indices[a];
    sizeOfPackedChars += Pack<DO_PACKING>( buffer, localToGlobalMap[li] );

    map< localIndex, array1d<globalIndex> >::const_iterator
    iterUnmappedGI = unmappedGlobalIndices.find(li);

    array1d<globalIndex> junk;
    array1d<globalIndex> const & unmappedGI = iterUnmappedGI==unmappedGlobalIndices.end() ?
                                              junk :
                                              iterUnmappedGI->second;

    sizeOfPackedChars += Pack<DO_PACKING>( buffer,
                                           var[li],
                                           unmappedGI,
                                           var.size(1),
                                           relatedObjectLocalToGlobalMap );
  }

  return sizeOfPackedChars;
}


inline
localIndex
Unpack( char const *& buffer,
        arrayView2d<localIndex> & var,
        localIndex_array & indices,
        map< localIndex, array1d<globalIndex> > & unmappedGlobalIndices,
        map<globalIndex,localIndex> const & globalToLocalMap,
        map<globalIndex,localIndex> const & relatedObjectGlobalToLocalMap )
{
  localIndex sizeOfUnpackedChars = 0;
  localIndex const sizeOfIndicesPassedIn = indices.size();

  localIndex numIndicesUnpacked;
  sizeOfUnpackedChars += Unpack( buffer, numIndicesUnpacked );
  GEOS_ERROR_IF( sizeOfIndicesPassedIn!=0 && numIndicesUnpacked!=indices.size(),
                 "number of unpacked indices("<<numIndicesUnpacked<<") does not equal size of "
                 "indices passed into Unpack function("<<sizeOfIndicesPassedIn );

  indices.resize(numIndicesUnpacked);

  for( localIndex a=0 ; a<numIndicesUnpacked ; ++a )
  {
    globalIndex gi;
    sizeOfUnpackedChars += Unpack( buffer, gi );

    localIndex & li = indices[a];
    if( sizeOfIndicesPassedIn > 0 )
    {
      GEOS_ERROR_IF( li!=globalToLocalMap.at(gi),
                     "global index "<<gi<<" unpacked from buffer does equal the lookup "
                     <<li<<" for localIndex "<<li<<" on this rank");
    }
    else
    {
      li = globalToLocalMap.at(gi);
    }

    arraySlice1d<localIndex> varSlice = var[li];
    array1d<globalIndex> unmappedIndices;

    sizeOfUnpackedChars += Unpack( buffer,
                                   varSlice,
                                   unmappedIndices,
                                   var.size(1),
                                   relatedObjectGlobalToLocalMap );

    if( unmappedIndices.size()>0 )
    {
      unmappedGlobalIndices[li] = unmappedIndices;
    }

  }

  return sizeOfUnpackedChars;
}


template< bool DO_PACKING, typename T_KEY, typename T_VAL >
typename std::enable_if< bufferOps::is_packable_map< map< T_KEY, T_VAL > >::value, localIndex >::type
Pack( char*& buffer, map<T_KEY, T_VAL> const & var )
{
  const typename std::map<T_KEY,T_VAL>::size_type length = var.size();
  localIndex sizeOfPackedChars = Pack<DO_PACKING>( buffer, length );

  for( typename std::map<T_KEY,T_VAL>::const_iterator i=var.begin() ; i!=var.end() ; ++i )
  {
    sizeOfPackedChars += Pack<DO_PACKING>( buffer, i->first );
    sizeOfPackedChars += Pack<DO_PACKING>( buffer, i->second );
  }

  return sizeOfPackedChars;
}


template< typename T_KEY, typename T_VAL >
typename std::enable_if< bufferOps::is_packable_map< map<T_KEY,T_VAL> >::value, localIndex >::type
Unpack( char const *& buffer, std::map<T_KEY,T_VAL>& map )
{
  map.clear();
  typename std::map<T_KEY,T_VAL>::size_type map_length;

  localIndex sizeOfUnpackedChars = Unpack( buffer, map_length );

  for( typename std::map<T_KEY,T_VAL>::size_type a=0 ; a<map_length ; ++a )
  {
    T_KEY key;
    T_VAL value;
    sizeOfUnpackedChars += Unpack( buffer, key );
    sizeOfUnpackedChars += Unpack( buffer, value );

    map[key] = std::move(value);
  }

  return sizeOfUnpackedChars;
}


template< bool DO_PACKING, typename T_KEY, typename T_VAL, typename T_INDICES >
typename std::enable_if< bufferOps::is_packable_map< map<T_KEY, T_VAL> >::value && bufferOps::is_packable_by_index<T_VAL>::value, localIndex >::type
Pack( char*& buffer, std::map<T_KEY,T_VAL> const & var, T_INDICES const & packIndices )
{
  const typename std::map<T_KEY,T_VAL>::size_type length = var.size();
  localIndex sizeOfPackedChars = Pack<DO_PACKING>( buffer, length );

  for( typename std::map<T_KEY,T_VAL>::const_iterator i=var.begin() ; i!=var.end() ; ++i )
  {
    sizeOfPackedChars += Pack<DO_PACKING>( buffer, i->first );
    sizeOfPackedChars += Pack<DO_PACKING>( buffer, i->second, packIndices );
  }

  return sizeOfPackedChars;
}


template< typename T_KEY, typename T_VAL, typename T_INDICES >
typename std::enable_if< bufferOps::is_packable_map< map<T_KEY,T_VAL> >::value && bufferOps::is_packable_by_index<T_VAL>::value, localIndex >::type
Unpack( char const *& buffer, std::map<T_KEY,T_VAL>& map, T_INDICES const & unpackIndices )
{
  map.clear();
  typename std::map<T_KEY,T_VAL>::size_type map_length;

  localIndex sizeOfUnpackedChars = Unpack( buffer, map_length );

  for( typename std::map<T_KEY,T_VAL>::size_type a=0 ; a<map_length ; ++a )
  {
    T_KEY key;
    T_VAL value;
    sizeOfUnpackedChars += Unpack( buffer, key );
    sizeOfUnpackedChars += Unpack( buffer, value, unpackIndices );

    map[key] = value;
  }

  return sizeOfUnpackedChars;
}


template< bool DO_PACKING, typename T_FIRST, typename T_SECOND >
localIndex
Pack( char*& buffer, std::pair< T_FIRST, T_SECOND > const & var )
{
  localIndex sizeOfPackedChars = Pack<DO_PACKING>(buffer, var.first);
  sizeOfPackedChars += Pack<DO_PACKING>(buffer, var.second);
  return sizeOfPackedChars;
}


template< typename T_FIRST, typename T_SECOND >
localIndex
Unpack( char const *& buffer, std::pair< T_FIRST, T_SECOND >& var )
{
  localIndex sizeOfUnpackedChars = Unpack(buffer, var.first);
  sizeOfUnpackedChars += Unpack(buffer, var.second);
  return sizeOfUnpackedChars;
}

} /* namespace bufferOps */
} /* namespace geosx */

#endif /* DATAREPOSITORY_BUFFEROPS_INLINE_H_ */
