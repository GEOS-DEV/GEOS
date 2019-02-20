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



#ifndef DATAREPOSITORY_BUFFEROPS_H_
#define DATAREPOSITORY_BUFFEROPS_H_

#include "common/DataTypes.hpp"
#include "codingUtilities/Utilities.hpp"
#include "codingUtilities/static_if.hpp"
#include "codingUtilities/GeosxTraits.hpp"
#include "IntegerConversion.hpp"

#include <type_traits>

namespace geosx
{

/* Forward declaration of InterObjectRelation */
template < typename T >
class InterObjectRelation;

namespace bufferOps
{

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename T >
typename std::enable_if< std::is_trivial<T>::value, localIndex >::type
Pack( char*&  buffer,
      T const & var );

//------------------------------------------------------------------------------
template< typename T>
typename std::enable_if< std::is_trivial<T>::value, localIndex >::type
Unpack( char const *& buffer,
        T & var );

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename T, typename INDEX_TYPE >
typename std::enable_if< std::is_trivial<T>::value, localIndex >::type
Pack( char*&  buffer,
      T const * const restrict var,
      INDEX_TYPE const length );

//------------------------------------------------------------------------------
template< typename T, typename INDEX_TYPE >
typename std::enable_if< std::is_trivial<T>::value, localIndex >::type
Unpack( char const *& buffer,
        T * const restrict var,
        INDEX_TYPE const expectedLength );

//------------------------------------------------------------------------------
template< bool DO_PACKING >
localIndex
Pack( char*& buffer,
      const std::string& var );

//------------------------------------------------------------------------------
inline
localIndex
Unpack( char const *& buffer,
        string& var );

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename T >
typename std::enable_if< traits::is_tensorT<T>::value, localIndex >::type
Pack( char*& buffer,
      T const & var );

//------------------------------------------------------------------------------
template< typename T >
typename std::enable_if< traits::is_tensorT<T>::value, localIndex >::type
Unpack( char const *& buffer,
        T & var );

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename T, typename INDEX_TYPE >
typename std::enable_if< !std::is_trivial<T>::value, localIndex >::type
Pack( char*&  buffer,
      T const * const restrict var,
      INDEX_TYPE const length );

//------------------------------------------------------------------------------
template< typename T, typename INDEX_TYPE >
typename std::enable_if< !std::is_trivial<T>::value, localIndex >::type
Unpack( char const *& buffer,
        T * const restrict var,
        INDEX_TYPE const expectedLength );

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename T, typename INDEX_TYPE >
localIndex
Pack( char*&  buffer,
      T const * const restrict var,
      arraySlice1d<INDEX_TYPE const> const & indices,
      INDEX_TYPE const length );

//------------------------------------------------------------------------------
template< typename T, typename INDEX_TYPE >
localIndex
Unpack( char const *& buffer,
        T * const restrict var,
        arraySlice1d<INDEX_TYPE const> const & indices,
        INDEX_TYPE & length );

#ifdef GEOSX_USE_ARRAY_BOUNDS_CHECK

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename T, typename INDEX_TYPE >
typename std::enable_if< std::is_trivial<T>::value, localIndex >::type
Pack( char*& buffer,
      arraySlice1d<T> const & var,
      INDEX_TYPE const length );

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename T, typename INDEX_TYPE >
typename std::enable_if< !std::is_trivial<T>::value, localIndex >::type
Pack( char*& buffer,
      arraySlice1d<T> const & var,
      INDEX_TYPE const length );

//------------------------------------------------------------------------------
template< typename T, typename INDEX_TYPE >
typename std::enable_if< std::is_trivial<T>::value, localIndex >::type
Unpack( char const *& buffer,
        arraySlice1d<T> & var,
        INDEX_TYPE const expectedLength );

//------------------------------------------------------------------------------
template< typename T, typename INDEX_TYPE >
typename std::enable_if< !std::is_trivial<T>::value, localIndex >::type
Unpack( char const *& buffer,
        arraySlice1d<T> & var,
        INDEX_TYPE const expectedLength );

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename T, typename INDEX_TYPE >
localIndex
Pack( char*&  buffer,
      arraySlice1d<T> const & var,
      arraySlice1d<INDEX_TYPE> const & indices,
      INDEX_TYPE const length );

//------------------------------------------------------------------------------
template< typename T, typename INDEX_TYPE >
localIndex
Unpack( char const *& buffer,
        arraySlice1d<T> & var,
        arraySlice1d<INDEX_TYPE> const & indices,
        INDEX_TYPE & length );

#endif /* GEOSX_USE_ARRAY_BOUNDS_CHECK */


//------------------------------------------------------------------------------
template< bool DO_PACKING, typename T >
localIndex
Pack( char *& buffer, set<T> const & var );

//------------------------------------------------------------------------------
template< typename T>
localIndex 
Unpack( char const *& buffer, set<T> & var );

//------------------------------------------------------------------------------
template< bool DO_PACKING >
localIndex 
Pack( char *& buffer,
      set<localIndex> const & var,
      set<globalIndex> const & unmappedGlobalIndices,
      arraySlice1d<globalIndex const> const & localToGlobal );

//------------------------------------------------------------------------------
inline
localIndex
Unpack( char const *& buffer,
        set<localIndex> & var,
        set<globalIndex> & unmappedGlobalIndices,
        map<globalIndex,localIndex> const & globalToLocalMap,
        bool const clearFlag );

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename T, int NDIM, typename INDEX_TYPE >
typename std::enable_if< bufferOps::is_packable_array< LvArray::ArrayView<T,NDIM,INDEX_TYPE> >::value, localIndex >::type
Pack( char*& buffer,
      LvArray::ArrayView<T,NDIM,INDEX_TYPE> const & var );

//------------------------------------------------------------------------------
template< typename T, int NDIM, typename INDEX_TYPE >
typename std::enable_if< bufferOps::is_packable_array< LvArray::Array<T,NDIM,INDEX_TYPE> >::value, localIndex >::type
Unpack( char const *& buffer, LvArray::Array<T,NDIM,INDEX_TYPE> & var );

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename T, int NDIM, typename T_indices, typename INDEX_TYPE >
typename std::enable_if< bufferOps::is_packable_array< LvArray::ArrayView<T,NDIM,INDEX_TYPE> >::value, localIndex >::type
Pack( char*& buffer,
      LvArray::ArrayView<T,NDIM,INDEX_TYPE> const & var,
      const T_indices& indices );

//------------------------------------------------------------------------------
template< typename T, int NDIM, typename T_indices, typename INDEX_TYPE >
localIndex
Unpack( char const *& buffer,
        LvArray::ArrayView<T,NDIM,INDEX_TYPE> & var,
        const T_indices& indices );

//------------------------------------------------------------------------------
template< bool DO_PACKING >
localIndex
Pack( char*& buffer,
      arraySlice1d<localIndex const> const & var,
      arraySlice1d<globalIndex const> const & unmappedGlobalIndices,
      localIndex const length,
      arraySlice1d<globalIndex const> const & localToGlobalMap );

//------------------------------------------------------------------------------
inline
localIndex
Unpack( char const *& buffer,
        localIndex_array & var,
        array1d<globalIndex> & unmappedGlobalIndices,
        map<globalIndex,localIndex> const & globalToLocalMap );

//------------------------------------------------------------------------------
inline
localIndex
Unpack( char const *& buffer,
        arraySlice1d<localIndex> & var,
        array1d<globalIndex> & unmappedGlobalIndices,
        localIndex const expectedLength,
        map<globalIndex,localIndex> const & globalToLocalMap );

//------------------------------------------------------------------------------
template< bool DO_PACKING >
localIndex
Pack( char*& buffer,
      arrayView1d<localIndex const> const & var,
      arrayView1d<localIndex const> const & indices,
      arrayView1d<globalIndex const> const & localToGlobalMap,
      arrayView1d<globalIndex const> const & relatedObjectLocalToGlobalMap );

//------------------------------------------------------------------------------
inline
localIndex
Unpack( char const *& buffer,
        arrayView1d<localIndex> & var,
        array1d<localIndex> const & indices,
        map<globalIndex,localIndex> const & globalToLocalMap,
        map<globalIndex,localIndex> const & relatedObjectGlobalToLocalMap );


//------------------------------------------------------------------------------
template< bool DO_PACKING >
localIndex
Pack( char*& buffer,
      arrayView1d<localIndex_array const> const & var,
      map< localIndex, array1d<globalIndex> > const & unmappedGlobalIndices,
      arrayView1d<localIndex const> const & indices,
      arrayView1d<globalIndex const> const & localToGlobalMap,
      arrayView1d<globalIndex const> const & relatedObjectLocalToGlobalMap );

//------------------------------------------------------------------------------
inline
localIndex
Unpack( char const *& buffer,
        arrayView1d<localIndex_array> & var,
        array1d<localIndex> & indices,
        map<localIndex, array1d<globalIndex> > & unmappedGlobalIndices,
        map<globalIndex,localIndex> const & globalToLocalMap,
        map<globalIndex,localIndex> const & relatedObjectGlobalToLocalMap );

//------------------------------------------------------------------------------
template< bool DO_PACKING >
localIndex
Pack( char*& buffer,
      arrayView1d< set<localIndex> const > const & var,
      map< localIndex, set<globalIndex> > const & unmappedGlobalIndices,
      arrayView1d<localIndex const> const & indices,
      arrayView1d<globalIndex const> const & localToGlobalMap,
      arrayView1d<globalIndex const> const & relatedObjectLocalToGlobalMap );

//------------------------------------------------------------------------------
inline
localIndex
Unpack( char const *& buffer,
        arrayView1d< set<localIndex> > & var,
        localIndex_array & indices,
        map<localIndex, set<globalIndex> > & unmappedGlobalIndices,
        map<globalIndex,localIndex> const & globalToLocalMap,
        map<globalIndex,localIndex> const & relatedObjectGlobalToLocalMap,
        bool const clearFlag );

//------------------------------------------------------------------------------
template< bool DO_PACKING >
localIndex
Pack( char*& buffer,
      arrayView2d<localIndex> const & var,
      arrayView1d<localIndex> const & indices,
      arraySlice1d<globalIndex> const & localToGlobalMap );

//------------------------------------------------------------------------------
inline
localIndex
Unpack( char const *& buffer,
        arrayView2d<localIndex> & var,
        array1d<localIndex> & indices,
        map<globalIndex,localIndex> const & globalToLocalMap );

//------------------------------------------------------------------------------
template< bool DO_PACKING >
localIndex
Pack( char*& buffer,
      arrayView2d<localIndex const> const & var,
      map< localIndex, array1d<globalIndex> > const & unmappedGlobalIndices,
      arrayView1d<localIndex const> const & indices,
      arraySlice1d<globalIndex const> const & localToGlobalMap,
      arraySlice1d<globalIndex const> const & relatedObjectLocalToGlobalMap );

//------------------------------------------------------------------------------
inline
localIndex
Unpack( char const *& buffer,
        arrayView2d<localIndex> & var,
        localIndex_array & indices,
        map< localIndex, array1d<globalIndex> > & unmappedGlobalIndices,
        map<globalIndex,localIndex> const & globalToLocalMap,
        map<globalIndex,localIndex> const & relatedObjectGlobalToLocalMap );

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename T_KEY, typename T_VAL >
typename std::enable_if< bufferOps::is_packable_map< map< T_KEY, T_VAL > >::value, localIndex >::type
Pack( char*& buffer, map<T_KEY, T_VAL> const & var );

//------------------------------------------------------------------------------
template< typename T_KEY, typename T_VAL >
typename std::enable_if< bufferOps::is_packable_map< map<T_KEY,T_VAL> >::value, localIndex >::type
Unpack( char const *& buffer, std::map<T_KEY,T_VAL>& map );

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename T_KEY, typename T_VAL, typename T_INDICES >
typename std::enable_if< bufferOps::is_packable_map< map<T_KEY, T_VAL> >::value && bufferOps::is_packable_by_index<T_VAL>::value, localIndex >::type
Pack( char*& buffer, std::map<T_KEY,T_VAL> const & var, T_INDICES const & packIndices );


//------------------------------------------------------------------------------
template< typename T_KEY, typename T_VAL, typename T_INDICES >
typename std::enable_if< bufferOps::is_packable_map< map<T_KEY,T_VAL> >::value && bufferOps::is_packable_by_index<T_VAL>::value, localIndex >::type
Unpack( char const *& buffer, std::map<T_KEY,T_VAL>& map, T_INDICES const & unpackIndices );

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename T_FIRST, typename T_SECOND >
localIndex
Pack( char*& buffer, std::pair< T_FIRST, T_SECOND > const & var );

//------------------------------------------------------------------------------
template< typename T_FIRST, typename T_SECOND >
localIndex
Unpack( char const *& buffer, std::pair< T_FIRST, T_SECOND >& var );

//------------------------------------------------------------------------------
template < bool DO_PACKING, typename T >
localIndex
Pack( char*& buffer, InterObjectRelation<T> const & var )
{
  return Pack< DO_PACKING >(buffer, static_cast<T const &>(var));
}

//------------------------------------------------------------------------------
template < typename T >
localIndex
Unpack( char const *& buffer, InterObjectRelation<T> & var )
{
  return Unpack(buffer, static_cast<T&>(var));
}

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename T >
typename std::enable_if< !bufferOps::is_packable<T>::value, localIndex >::type
Pack( char*&  buffer, T const & var )
{
  GEOS_ERROR("Trying to pack data type ("<<typeid(T).name()<<") but type is not packable.");
  return 0;
}

//------------------------------------------------------------------------------
template< typename T >
typename std::enable_if< !bufferOps::is_packable<T>::value, localIndex >::type
Unpack( char const *& buffer, T & var )
{
  GEOS_ERROR("Trying to unpack data type ("<<typeid(T).name()<<") but type is not packable.");
  return 0;
}

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename T, typename T_INDICES >
typename std::enable_if< !bufferOps::is_packable_by_index<T>::value, localIndex >::type
Pack( char*&  buffer, T const & var, T_INDICES const& )
{
  GEOS_ERROR("Trying to pack data type ("<<typeid(T).name()<<") but type is not packable by index.");
  return 0;
}

//------------------------------------------------------------------------------
template< typename T, typename T_INDICES >
typename std::enable_if< !bufferOps::is_packable_by_index<T>::value, localIndex >::type
Unpack( char const *& buffer, T & var, T_INDICES const& )
{
  GEOS_ERROR("Trying to unpack data type ("<<typeid(T).name()<<") but type is not packable by index.");
  return 0;
}

//------------------------------------------------------------------------------
template< typename... VARPACK >
localIndex
PackSize( VARPACK const &&...  pack )
{
  char* junk = nullptr;
  return Pack<false>( junk, pack... );
}

//------------------------------------------------------------------------------
template< typename... VARPACK >
localIndex
PackSize( VARPACK &&...  pack )
{
  char* junk = nullptr;
  return Pack<false>( junk, pack... );
}

} /* namespace bufferOps */
} /* namespace geosx */

#include "BufferOps_inline.hpp"

#endif /* DATAREPOSITORY_BUFFEROPS_H_ */
