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
typename std::enable_if< std::is_trivial< T >::value, localIndex >::type
Pack( buffer_unit_type * & buffer, T const & var )
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
template< typename T >
typename std::enable_if< std::is_trivial< T >::value, localIndex >::type
Unpack( buffer_unit_type const * & buffer, T & var )
{
  localIndex const sizeOfUnpackedChars = sizeof(T);
  memcpy( &var, buffer, sizeOfUnpackedChars );
  buffer += sizeOfUnpackedChars;
  return sizeOfUnpackedChars;
}


template< bool DO_PACKING, typename T, typename INDEX_TYPE >
typename std::enable_if< std::is_trivial< T >::value, localIndex >::type
Pack( buffer_unit_type * & buffer, T const * const restrict var, INDEX_TYPE const length )
{
  localIndex sizeOfPackedChars = Pack< DO_PACKING >( buffer, length );

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
typename std::enable_if< std::is_trivial< T >::value, localIndex >::type
Unpack( buffer_unit_type const * & buffer, T * const restrict var, INDEX_TYPE const expectedLength )
{
  INDEX_TYPE length;
  localIndex sizeOfUnpackedChars = Unpack( buffer, length );

  GEOS_ASSERT_MSG( length == expectedLength, "expectedLength != length: " <<
                   expectedLength << " != " << length );
  GEOSX_DEBUG_VAR( expectedLength );

  memcpy( var, buffer, length * sizeof(T) );
  sizeOfUnpackedChars += length * sizeof(T);
  buffer += length * sizeof(T);

  return sizeOfUnpackedChars;
}


template< bool DO_PACKING >
localIndex Pack( buffer_unit_type * & buffer, const std::string & var )
{
  const string::size_type varSize = var.size();
  localIndex sizeOfPackedChars = Pack< DO_PACKING >( buffer, varSize );

  static_if( DO_PACKING )
  {
    memcpy( buffer, var.data(), varSize );
    buffer += varSize;
  }
  end_static_if

    sizeOfPackedChars += varSize;
  return sizeOfPackedChars;
}


inline localIndex Unpack( buffer_unit_type const * & buffer, string & var )
{
  string::size_type stringsize = 0;
  localIndex sizeOfUnpackedChars = Unpack( buffer, stringsize );

  var.resize( stringsize );
  var.assign( reinterpret_cast< char const * >( buffer ), stringsize );
  buffer += stringsize;
  sizeOfUnpackedChars += stringsize;

  return sizeOfUnpackedChars;
}


template< bool DO_PACKING, typename T >
typename std::enable_if< traits::is_tensorT< T >, localIndex >::type
Pack( buffer_unit_type * & buffer, T const & var )
{
  localIndex sizeOfPackedChars = 0;
  sizeOfPackedChars += Pack< DO_PACKING >( buffer, var.Data(), var.Length());
  return sizeOfPackedChars;
}


template< typename T >
typename std::enable_if< traits::is_tensorT< T >, localIndex >::type
Unpack( buffer_unit_type const * & buffer, T & var )
{
  localIndex sizeOfUnpackedChars = 0;
  real64 * const pVar = var.Data();
  int const length = var.Length();
  sizeOfUnpackedChars += Unpack( buffer, pVar, length );
  return sizeOfUnpackedChars;
}


template< bool DO_PACKING, typename T, typename INDEX_TYPE >
typename std::enable_if< !std::is_trivial< T >::value, localIndex >::type
Pack( buffer_unit_type * & buffer, T const * const restrict var, INDEX_TYPE const length )
{
  localIndex sizeOfPackedChars = Pack< DO_PACKING >( buffer, length );

  for( INDEX_TYPE a=0 ; a<length ; ++a )
  {
    sizeOfPackedChars += Pack< DO_PACKING >( buffer, var[ a ] );
  }

  return sizeOfPackedChars;
}


template< typename T, typename INDEX_TYPE >
typename std::enable_if< !std::is_trivial< T >::value, localIndex >::type
Unpack( buffer_unit_type const * & buffer, T * const restrict var, INDEX_TYPE const expectedLength )
{
  INDEX_TYPE length;
  localIndex sizeOfUnpackedChars = Unpack( buffer, length );

  GEOS_ASSERT_MSG( length == expectedLength, "expectedLength != length: " <<
                   expectedLength << " != " << length );
  GEOSX_DEBUG_VAR( expectedLength );

  for( INDEX_TYPE a=0 ; a<length ; ++a )
  {
    sizeOfUnpackedChars += Unpack( buffer, var[a] );
  }

  return sizeOfUnpackedChars;
}


template< bool DO_PACKING, typename T, typename INDEX_TYPE, int UNIT_STRIDE_DIM >
localIndex
Pack( buffer_unit_type * & buffer,
      T const * const restrict var,
      arraySlice1d< INDEX_TYPE const, UNIT_STRIDE_DIM > const & indices,
      INDEX_TYPE const length )
{
  localIndex sizeOfPackedChars = Pack< DO_PACKING >( buffer, length );

  for( INDEX_TYPE a=0 ; a<length ; ++a )
  {
    sizeOfPackedChars += Pack< DO_PACKING >( buffer, var[ indices[a] ] );
  }

  return sizeOfPackedChars;
}


template< typename T, typename INDEX_TYPE >
localIndex
Unpack( buffer_unit_type const * & buffer,
        T * const restrict var,
        arraySlice1d< INDEX_TYPE const > const & indices,
        INDEX_TYPE & length )
{
  localIndex sizeOfUnpackedChars = Unpack( buffer, length );

  for( INDEX_TYPE a=0 ; a<length ; ++a )
  {
    sizeOfUnpackedChars += Unpack( buffer, var[ indices[a] ] );
  }

  return sizeOfUnpackedChars;
}


template< bool DO_PACKING, typename T, typename INDEX_TYPE >
typename std::enable_if< std::is_trivial< T >::value, localIndex >::type
Pack( buffer_unit_type * & buffer,
      arraySlice1d< T const > const & var,
      INDEX_TYPE const length )
{
  localIndex sizeOfPackedChars = Pack< DO_PACKING >( buffer, length );
  sizeOfPackedChars += length * sizeof(T);

  static_if( DO_PACKING )
  {
    T * const restrict buffer_T = reinterpret_cast< T * >( buffer );
    for( INDEX_TYPE i = 0 ; i < length ; ++i )
    {
      buffer_T[ i ] = var[ i ];
    }

    buffer += length * sizeof(T);
  }
  end_static_if

  return sizeOfPackedChars;
}


template< bool DO_PACKING, typename T, typename INDEX_TYPE, int UNIT_STRIDE_DIM >
typename std::enable_if< !std::is_trivial< T >::value, localIndex >::type
Pack( buffer_unit_type * & buffer,
      arraySlice1d< T const, UNIT_STRIDE_DIM > const & var,
      INDEX_TYPE const length )
{
  localIndex sizeOfPackedChars = Pack< DO_PACKING >( buffer, length );

  for( INDEX_TYPE a=0 ; a<length ; ++a )
  {
    sizeOfPackedChars += Pack< DO_PACKING >( buffer, var[ a ] );
  }

  return sizeOfPackedChars;
}


template< typename T, typename INDEX_TYPE, int UNIT_STRIDE_DIM >
typename std::enable_if< std::is_trivial< T >::value, localIndex >::type
Unpack( buffer_unit_type const * & buffer,
        arraySlice1d< T, UNIT_STRIDE_DIM > const & var,
        INDEX_TYPE const expectedLength )
{
  INDEX_TYPE length;
  localIndex sizeOfUnpackedChars = Unpack( buffer, length );

  GEOSX_DEBUG_VAR( expectedLength );
  GEOS_ASSERT_MSG( length == expectedLength, "expectedLength != length: " <<
                   expectedLength << " != " << length );

  T const * const restrict buffer_T = reinterpret_cast< T const * >( buffer );
  for( INDEX_TYPE i = 0 ; i < length ; ++i )
  {
    var[ i ] = buffer_T[ i ];
  }

  buffer += length * sizeof(T);
  sizeOfUnpackedChars += length * sizeof(T);

  return sizeOfUnpackedChars;
}


template< typename T, typename INDEX_TYPE, int UNIT_STRIDE_DIM >
typename std::enable_if< !std::is_trivial< T >::value, localIndex >::type
Unpack( buffer_unit_type const * & buffer,
        arraySlice1d< T, UNIT_STRIDE_DIM > const & var,
        INDEX_TYPE const expectedLength )
{
  INDEX_TYPE length;
  localIndex sizeOfUnpackedChars = Unpack( buffer, length );

  GEOSX_DEBUG_VAR( expectedLength );
  GEOS_ASSERT_MSG( length == expectedLength, "expectedLength != length: " <<
                   expectedLength << " != " << length );

  for( INDEX_TYPE a=0 ; a<length ; ++a )
  {
    sizeOfUnpackedChars += Unpack( buffer, var[a] );
  }

  return sizeOfUnpackedChars;
}


template< bool DO_PACKING, typename T, typename INDEX_TYPE, int UNIT_STRIDE_DIM0, int UNIT_STRIDE_DIM1 >
localIndex
Pack( buffer_unit_type * & buffer,
      arraySlice1d< T const, UNIT_STRIDE_DIM0 > const & var,
      arraySlice1d< INDEX_TYPE const, UNIT_STRIDE_DIM1 > const & indices,
      INDEX_TYPE const length )
{
  localIndex sizeOfPackedChars = Pack< DO_PACKING >( buffer, length );

  for( INDEX_TYPE a=0 ; a<length ; ++a )
  {
    sizeOfPackedChars += Pack< DO_PACKING >( buffer, var[ indices[a] ] );
  }

  return sizeOfPackedChars;
}


template< typename T, typename INDEX_TYPE, int UNIT_STRIDE_DIM0, int UNIT_STRIDE_DIM1 >
localIndex
Unpack( buffer_unit_type const * & buffer,
        arraySlice1d< T, UNIT_STRIDE_DIM0 > const & var,
        arraySlice1d< INDEX_TYPE const, UNIT_STRIDE_DIM1 > const & indices,
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


template< bool DO_PACKING, typename T >
localIndex Pack( buffer_unit_type * & buffer, set< T > const & var )
{
  const localIndex length = integer_conversion< localIndex >( var.size());
  localIndex sizeOfPackedChars = Pack< DO_PACKING >( buffer, length );

  for( typename set< T >::const_iterator i=var.begin() ; i!=var.end() ; ++i )
  {
    sizeOfPackedChars += Pack< DO_PACKING >( buffer, *i );
  }

  return sizeOfPackedChars;
}


template< typename T >
localIndex Unpack( buffer_unit_type const * & buffer, set< T > & var )
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

template< bool DO_PACKING, int UNIT_STRIDE_DIM >
localIndex Pack( buffer_unit_type * & buffer,
                 set< localIndex > const & var,
                 set< globalIndex > const & unmappedGlobalIndices,
                 arraySlice1d< globalIndex const, UNIT_STRIDE_DIM > const & localToGlobal )
{
  const localIndex length = integer_conversion< localIndex >( var.size()+unmappedGlobalIndices.size());
  localIndex sizeOfPackedChars = Pack< DO_PACKING >( buffer, length );

  for( typename set< localIndex >::const_iterator i=var.begin() ; i!=var.end() ; ++i )
  {
    sizeOfPackedChars += Pack< DO_PACKING >( buffer, localToGlobal[*i] );
  }

  for( typename set< globalIndex >::const_iterator i=unmappedGlobalIndices.begin() ;
       i!=unmappedGlobalIndices.end() ; ++i )
  {
    sizeOfPackedChars += Pack< DO_PACKING >( buffer, *i );
  }


  return sizeOfPackedChars;
}

template< typename SORTED >
inline
localIndex Unpack( buffer_unit_type const * & buffer,
                   set< localIndex > & var,
                   set< globalIndex > & unmappedGlobalIndices,
                   mapBase< globalIndex, localIndex, SORTED > const & globalToLocalMap,
                   bool const clearExistingSet )
{
  if( clearExistingSet )
  {
    var.clear();
  }
  localIndex set_length;
  localIndex sizeOfUnpackedChars = Unpack( buffer, set_length );

  for( localIndex a=0 ; a<set_length ; ++a )
  {
    globalIndex temp;
    sizeOfUnpackedChars += Unpack( buffer, temp );
    typename mapBase< globalIndex, localIndex, SORTED >::const_iterator iter = globalToLocalMap.find( temp );
    if( iter==globalToLocalMap.end() )
    {
      unmappedGlobalIndices.insert( temp );
    }
    else
    {
      var.insert( iter->second );
    }
  }

  return sizeOfUnpackedChars;
}

template< typename SORTED >
inline
localIndex Unpack( buffer_unit_type const * & buffer,
                   ArrayOfSets< localIndex > & var,
                   localIndex const setIndex,
                   set< globalIndex > & unmappedGlobalIndices,
                   mapBase< globalIndex, localIndex, SORTED > const & globalToLocalMap,
                   bool const clearExistingSet )
{
  if( clearExistingSet )
  {
    var.clearSet( setIndex );
  }
  localIndex set_length;
  localIndex sizeOfUnpackedChars = Unpack( buffer, set_length );

  for( localIndex a=0 ; a<set_length ; ++a )
  {
    globalIndex temp;
    sizeOfUnpackedChars += Unpack( buffer, temp );
    typename mapBase< globalIndex, localIndex, SORTED >::const_iterator iter = globalToLocalMap.find( temp );
    if( iter==globalToLocalMap.end() )
    {
      unmappedGlobalIndices.insert( temp );
    }
    else
    {
      var.insertIntoSet( setIndex, iter->second );
    }
  }

  return sizeOfUnpackedChars;
}



template< bool DO_PACKING >
localIndex Pack( buffer_unit_type * & buffer,
                 set< localIndex > const & var,
                 arrayView1d< localIndex const > const & packList,
                 set< globalIndex > const & unmappedGlobalIndices,
                 arraySlice1d< globalIndex const > const & localToGlobal )
{

  localIndex length = 0;
  array1d< localIndex > temp( var.size());

  for( auto a : packList )
  {
    if( var.count( a ) )
    {
      temp[length] = a;
      ++length;
    }
  }
  temp.resize( length );
  localIndex sizeOfPackedChars = Pack< DO_PACKING >( buffer, length );


  for( localIndex a=0 ; a<length ; ++a )
  {
    sizeOfPackedChars += Pack< DO_PACKING >( buffer, localToGlobal[temp[a]] );
  }

  for( typename set< globalIndex >::const_iterator i=unmappedGlobalIndices.begin() ;
       i!=unmappedGlobalIndices.end() ; ++i )
  {
    sizeOfPackedChars += Pack< DO_PACKING >( buffer, *i );
  }

  return sizeOfPackedChars;
}


template< bool DO_PACKING, typename T, int NDIM, int UNIT_STRIDE_DIM, typename INDEX_TYPE >
typename std::enable_if< is_packable< T >, localIndex >::type
Pack( buffer_unit_type * & buffer,
      LvArray::ArrayView< T, NDIM, UNIT_STRIDE_DIM, INDEX_TYPE > const & var )
{
  localIndex sizeOfPackedChars = Pack< DO_PACKING >( buffer, var.dims(), NDIM );
  sizeOfPackedChars += Pack< DO_PACKING >( buffer, var.strides(), NDIM );

  const localIndex length = var.size();
  T const * const data = var.data();
  sizeOfPackedChars += Pack< DO_PACKING >( buffer, data, length );

  return sizeOfPackedChars;
}


template< typename T, int NDIM, typename PERMUTATION, typename INDEX_TYPE >
typename std::enable_if< is_packable< T >, localIndex >::type
Unpack( buffer_unit_type const * & buffer, LvArray::Array< T, NDIM, PERMUTATION, INDEX_TYPE > & var )
{
  INDEX_TYPE dims[NDIM];
  localIndex sizeOfUnpackedChars = Unpack( buffer, dims, NDIM );

  var.resize( NDIM, dims );

  INDEX_TYPE strides[NDIM];
  sizeOfUnpackedChars += Unpack( buffer, strides, NDIM );

  for( int i=0 ; i<NDIM ; ++i )
  {
    GEOS_ASSERT_MSG( strides[i] == var.strides()[i], "Strides are inconsistent: " <<
                     strides[i] << " != " << var.strides()[i] );
  }

  sizeOfUnpackedChars += Unpack( buffer, var.data(), var.size() );

  return sizeOfUnpackedChars;
}


template< bool DO_PACKING, typename T, int NDIM, int UNIT_STRIDE_DIM, typename T_indices, typename INDEX_TYPE >
typename std::enable_if< is_packable< T >, localIndex >::type
Pack( buffer_unit_type * & buffer,
      LvArray::ArrayView< T, NDIM, UNIT_STRIDE_DIM, INDEX_TYPE > const & var,
      const T_indices & indices )
{
  localIndex sizeOfPackedChars = Pack< DO_PACKING >( buffer, var.strides(), NDIM );

  for( localIndex a=0 ; a<indices.size() ; ++a )
  {
    T const * const data = var.data( indices[a] );
    sizeOfPackedChars += Pack< DO_PACKING >( buffer, data, var.strides()[0] );
  }
  return sizeOfPackedChars;
}


template< typename T, int NDIM, int UNIT_STRIDE_DIM, typename T_indices, typename INDEX_TYPE >
localIndex
Unpack( buffer_unit_type const * & buffer,
        LvArray::ArrayView< T, NDIM, UNIT_STRIDE_DIM, INDEX_TYPE > & var,
        const T_indices & indices )
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


template< bool DO_PACKING, typename T, typename INDEX_TYPE >
localIndex
Pack( buffer_unit_type * & buffer,
      LvArray::ArrayOfArrays< T, INDEX_TYPE > const & var )
{
  localIndex sizeOfPackedChars = 0;

  sizeOfPackedChars += Pack< DO_PACKING >( buffer, var.size() );
  for( localIndex a=0 ; a<var.size() ; ++a )
  {
    sizeOfPackedChars += Pack< DO_PACKING >( buffer, var.sizeOfArray( a ) );
    T const * const data = var[a];
    sizeOfPackedChars += Pack< DO_PACKING >( buffer, data, var.sizeOfArray( a ) );
  }
  return sizeOfPackedChars;
}


template< typename T, typename INDEX_TYPE >
localIndex
Unpack( buffer_unit_type const * & buffer,
        LvArray::ArrayOfArrays< T, INDEX_TYPE > & var )
{
  localIndex sizeOfUnpackedChars = 0;
  localIndex numOfArrays;
  sizeOfUnpackedChars += Unpack( buffer, numOfArrays );
  var.resize( numOfArrays );

  for( localIndex a=0 ; a<numOfArrays ; ++a )
  {
    localIndex sizeOfArray;
    sizeOfUnpackedChars += Unpack( buffer, sizeOfArray );
    var.resizeArray( a, sizeOfArray );
    sizeOfUnpackedChars += Unpack( buffer, var[a], sizeOfArray );
  }

  return sizeOfUnpackedChars;
}

template< bool DO_PACKING, typename T, typename INDEX_TYPE >
localIndex
Pack( buffer_unit_type * & buffer,
      LvArray::ArrayOfSets< T, INDEX_TYPE > const & var )
{
  localIndex sizeOfPackedChars = 0;

  sizeOfPackedChars += Pack< DO_PACKING >( buffer, var.size() );
  for( localIndex a=0 ; a<var.size() ; ++a )
  {
    sizeOfPackedChars += Pack< DO_PACKING >( buffer, var.sizeOfSet( a ) );
    T const * const data = var[a];
    sizeOfPackedChars += Pack< DO_PACKING >( buffer, data, var.sizeOfSet( a ) );
  }
  return sizeOfPackedChars;
}

template< typename T, typename INDEX_TYPE >
localIndex
Unpack( buffer_unit_type const * & buffer,
        LvArray::ArrayOfSets< T, INDEX_TYPE > & var )
{
  LvArray::ArrayOfArrays< T, INDEX_TYPE > varAsArray;
  localIndex sizeOfUnpackedChars = Unpack( buffer, varAsArray );

  var.stealFrom( std::move( varAsArray ), sortedArrayManipulation::SORTED_UNIQUE );

  return sizeOfUnpackedChars;
}

template< bool DO_PACKING, typename T, typename INDEX_TYPE, typename T_indices >
localIndex
Pack( buffer_unit_type * & buffer,
      LvArray::ArrayOfArrays< T, INDEX_TYPE > const & var,
      T_indices const & indices )
{
  localIndex sizeOfPackedChars = 0;

  sizeOfPackedChars += Pack< DO_PACKING >( buffer, indices.size() );
  for( localIndex a=0 ; a<indices.size() ; ++a )
  {
    sizeOfPackedChars += Pack< DO_PACKING >( buffer, var.sizeOfArray( indices[a] ) );
    T const * const data = var[indices[a]];
    sizeOfPackedChars += Pack< DO_PACKING >( buffer, data, var.sizeOfArray( indices[a] ) );
  }
  return sizeOfPackedChars;
}


template< typename T, typename INDEX_TYPE, typename T_indices >
localIndex
Unpack( buffer_unit_type const * & buffer,
        LvArray::ArrayOfArrays< T, INDEX_TYPE > & var,
        T_indices const & indices )
{
  localIndex sizeOfUnpackedChars = 0;
  localIndex numUnpackedIndices;
  sizeOfUnpackedChars += Unpack( buffer, numUnpackedIndices );

  GEOS_ERROR_IF( numUnpackedIndices!=indices.size(), "number of unpacked indices does not equal expected number" );

  for( localIndex a=0 ; a<indices.size() ; ++a )
  {
    localIndex sizeOfSubArray;
    sizeOfUnpackedChars += Unpack( buffer, sizeOfSubArray );
    var.resizeArray( indices[a], sizeOfSubArray );
    sizeOfUnpackedChars += Unpack( buffer, var[indices[a]], sizeOfSubArray );
  }


  return sizeOfUnpackedChars;
}


template< bool DO_PACKING, int UNIT_STRIDE_DIM >
localIndex
Pack( buffer_unit_type * & buffer,
      arraySlice1d< localIndex const, UNIT_STRIDE_DIM > const & var,
      globalIndex const * const unmappedGlobalIndices,
      localIndex const length,
      arraySlice1d< globalIndex const > const & localToGlobalMap )
{
  localIndex sizeOfPackedChars = Pack< DO_PACKING >( buffer, length );
  sizeOfPackedChars += length*sizeof(globalIndex);

  static_if( DO_PACKING )
  {
    globalIndex * const buffer_GI = reinterpret_cast< globalIndex * >(buffer);
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

template< typename SORTED >
inline
localIndex
Unpack( buffer_unit_type const * & buffer,
        localIndex_array & var,
        array1d< globalIndex > & unmappedGlobalIndices,
        mapBase< globalIndex, localIndex, SORTED > const & globalToLocalMap )
{
  localIndex length;
  localIndex sizeOfUnpackedChars = Unpack( buffer, length );
  var.resize( length );
  unmappedGlobalIndices.resize( length );
  unmappedGlobalIndices = unmappedLocalIndexValue;

  bool unpackedGlobalFlag = false;
  for( localIndex a=0 ; a<length ; ++a )
  {
    globalIndex unpackedGlobalIndex;
    sizeOfUnpackedChars += Unpack( buffer, unpackedGlobalIndex );

    typename mapBase< globalIndex, localIndex, SORTED >::const_iterator
      iter = globalToLocalMap.find( unpackedGlobalIndex );
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

template< typename SORTED >
inline
localIndex
Unpack( buffer_unit_type const * & buffer,
        ArrayOfArrays< localIndex > & var,
        localIndex const subArrayIndex,
        array1d< globalIndex > & unmappedGlobalIndices,
        mapBase< globalIndex, localIndex, SORTED > const & globalToLocalMap )
{
  localIndex length;
  localIndex sizeOfUnpackedChars = Unpack( buffer, length );
  var.resizeArray( subArrayIndex, length );
  unmappedGlobalIndices.resize( length );
  unmappedGlobalIndices = unmappedLocalIndexValue;

  bool unpackedGlobalFlag = false;
  for( localIndex a=0 ; a<length ; ++a )
  {
    globalIndex unpackedGlobalIndex;
    sizeOfUnpackedChars += Unpack( buffer, unpackedGlobalIndex );

    typename mapBase< globalIndex, localIndex, SORTED >::const_iterator
      iter = globalToLocalMap.find( unpackedGlobalIndex );
    if( iter == globalToLocalMap.end() )
    {
      var( subArrayIndex, a ) = unmappedLocalIndexValue;
      unmappedGlobalIndices[a] = unpackedGlobalIndex;
      unpackedGlobalFlag = true;
    }
    else
    {
      var( subArrayIndex, a ) = iter->second;
    }
  }
  if( !unpackedGlobalFlag )
  {
    unmappedGlobalIndices.clear();
  }

  return sizeOfUnpackedChars;
}

template< typename SORTED, int UNIT_STRIDE_DIM >
inline
localIndex
Unpack( buffer_unit_type const * & buffer,
        arraySlice1d< localIndex, UNIT_STRIDE_DIM > const & var,
        array1d< globalIndex > & unmappedGlobalIndices,
        localIndex const expectedLength,
        mapBase< globalIndex, localIndex, SORTED > const & globalToLocalMap )
{
  localIndex sizeOfUnpackedChars = 0;

  localIndex length;
  sizeOfUnpackedChars += Unpack( buffer, length );

  GEOS_ASSERT_EQ( length, expectedLength );
  GEOSX_DEBUG_VAR( expectedLength );

  unmappedGlobalIndices.resize( length );
  unmappedGlobalIndices = unmappedLocalIndexValue;

  bool unpackedGlobalFlag = false;
  for( localIndex a=0 ; a<length ; ++a )
  {
    globalIndex unpackedGlobalIndex;
    sizeOfUnpackedChars += Unpack( buffer, unpackedGlobalIndex );

    typename mapBase< globalIndex, localIndex, SORTED >::const_iterator
      iter = globalToLocalMap.find( unpackedGlobalIndex );
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
Pack( buffer_unit_type * & buffer,
      arrayView1d< localIndex const > const & var,
      arrayView1d< localIndex const > const & indices,
      arrayView1d< globalIndex const > const & localToGlobalMap,
      arrayView1d< globalIndex const > const & relatedObjectLocalToGlobalMap )
{
  localIndex sizeOfPackedChars=0;

  sizeOfPackedChars += Pack< DO_PACKING >( buffer, indices.size() );
  for( localIndex a=0 ; a<indices.size() ; ++a )
  {
    localIndex const li = indices[a];
    sizeOfPackedChars += Pack< DO_PACKING >( buffer, localToGlobalMap[li] );
    if( var[li] != -1 )
    {
      sizeOfPackedChars += Pack< DO_PACKING >( buffer,
                                               relatedObjectLocalToGlobalMap[var[li]] );
    }
    else
    {
      sizeOfPackedChars += Pack< DO_PACKING >( buffer,
                                               localIndex( -1 ) );
    }
  }

  return sizeOfPackedChars;
}

template< typename SORTED0, typename SORTED1 >
inline
localIndex
Unpack( buffer_unit_type const * & buffer,
        arrayView1d< localIndex > & var,
        array1d< localIndex > const & indices,
        mapBase< globalIndex, localIndex, SORTED0 > const & globalToLocalMap,
        mapBase< globalIndex, localIndex, SORTED1 > const & relatedObjectGlobalToLocalMap )
{
  localIndex numIndicesUnpacked;
  localIndex const sizeOfIndicesPassedIn = indices.size();

  localIndex sizeOfUnpackedChars = Unpack( buffer, numIndicesUnpacked );

  GEOS_ERROR_IF( sizeOfIndicesPassedIn!=0 && numIndicesUnpacked!=indices.size(),
                 "number of unpacked indices("<<numIndicesUnpacked<<") does not equal size of "
                                                                    "indices passed into Unpack function("<<sizeOfIndicesPassedIn );

  for( localIndex a=0 ; a<indices.size() ; ++a )
  {
    globalIndex gi;
    sizeOfUnpackedChars += Unpack( buffer, gi );
    localIndex & li = indices[a];
    if( sizeOfIndicesPassedIn > 0 )
    {
      GEOS_ERROR_IF( li!=globalToLocalMap.at( gi ),
                     "global index "<<gi<<" unpacked from buffer does not equal the lookup "
                                    <<li<<" for localIndex "<<li<<" on this rank" );
    }
    else
    {
      li = globalToLocalMap.at( gi );
    }

    globalIndex mappedGlobalIndex;
    sizeOfUnpackedChars += Unpack( buffer, mappedGlobalIndex );
    if( mappedGlobalIndex != -1 )
    {
      var[li] = relatedObjectGlobalToLocalMap.at( mappedGlobalIndex );
    }
    else
    {
      var[li] = -1;
    }

  }
  return sizeOfUnpackedChars;
}

template< bool DO_PACKING, typename SORTED >
localIndex
Pack( buffer_unit_type * & buffer,
      arrayView1d< localIndex_array const > const & var,
      mapBase< localIndex, array1d< globalIndex >, SORTED > const & unmappedGlobalIndices,
      arrayView1d< localIndex const > const & indices,
      arrayView1d< globalIndex const > const & localToGlobalMap,
      arrayView1d< globalIndex const > const & relatedObjectLocalToGlobalMap )
{
  localIndex sizeOfPackedChars=0;

  sizeOfPackedChars += Pack< DO_PACKING >( buffer, indices.size() );
  for( localIndex a=0 ; a<indices.size() ; ++a )
  {
    localIndex const li = indices[a];
    sizeOfPackedChars += Pack< DO_PACKING >( buffer, localToGlobalMap[li] );

    typename mapBase< localIndex, array1d< globalIndex >, SORTED >::const_iterator
      iterUnmappedGI = unmappedGlobalIndices.find( li );

    array1d< globalIndex > junk;
    array1d< globalIndex > const & unmappedGI = iterUnmappedGI==unmappedGlobalIndices.end() ?
                                                junk :
                                                iterUnmappedGI->second;

    sizeOfPackedChars += Pack< DO_PACKING >( buffer,
                                             var[li].toSliceConst(),
                                             unmappedGI.data(),
                                             var[li].size(),
                                             relatedObjectLocalToGlobalMap );
  }

  return sizeOfPackedChars;
}

template< typename SORTED0, typename SORTED1, typename SORTED2 >
inline
localIndex
Unpack( buffer_unit_type const * & buffer,
        arrayView1d< localIndex_array > & var,
        array1d< localIndex > & indices,
        mapBase< localIndex, array1d< globalIndex >, SORTED0 > & unmappedGlobalIndices,
        mapBase< globalIndex, localIndex, SORTED1 > const & globalToLocalMap,
        mapBase< globalIndex, localIndex, SORTED2 > const & relatedObjectGlobalToLocalMap )
{
  localIndex numIndicesUnpacked;
  localIndex const sizeOfIndicesPassedIn = indices.size();

  localIndex sizeOfUnpackedChars = Unpack( buffer, numIndicesUnpacked );

  GEOS_ERROR_IF( sizeOfIndicesPassedIn!=0 && numIndicesUnpacked!=indices.size(),
                 "number of unpacked indices("<<numIndicesUnpacked<<") does not equal size of "
                                                                    "indices passed into Unpack function("<<sizeOfIndicesPassedIn );

  indices.resize( numIndicesUnpacked );

  for( localIndex a=0 ; a<indices.size() ; ++a )
  {
    globalIndex gi;
    sizeOfUnpackedChars += Unpack( buffer, gi );

    localIndex & li = indices[a];
    if( sizeOfIndicesPassedIn > 0 )
    {
      GEOS_ERROR_IF( li!=globalToLocalMap.at( gi ),
                     "global index "<<gi<<" unpacked from buffer does not equal the lookup "
                                    <<li<<" for localIndex "<<li<<" on this rank" );
    }
    else
    {
      li = globalToLocalMap.at( gi );
    }

    array1d< globalIndex > unmappedIndices;
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

template< bool DO_PACKING, typename SORTED >
localIndex
Pack( buffer_unit_type * & buffer,
      ArrayOfArraysView< localIndex const > const & var,
      mapBase< localIndex, array1d< globalIndex >, SORTED > const & unmappedGlobalIndices,
      arrayView1d< localIndex const > const & indices,
      arrayView1d< globalIndex const > const & localToGlobalMap,
      arrayView1d< globalIndex const > const & relatedObjectLocalToGlobalMap )
{
  localIndex sizeOfPackedChars=0;

  sizeOfPackedChars += Pack< DO_PACKING >( buffer, indices.size() );
  for( localIndex a=0 ; a<indices.size() ; ++a )
  {
    localIndex const li = indices[a];
    sizeOfPackedChars += Pack< DO_PACKING >( buffer, localToGlobalMap[li] );

    typename mapBase< localIndex, array1d< globalIndex >, SORTED >::const_iterator
      iterUnmappedGI = unmappedGlobalIndices.find( li );

    array1d< globalIndex > junk;
    array1d< globalIndex > const & unmappedGI = iterUnmappedGI==unmappedGlobalIndices.end() ?
                                                junk :
                                                iterUnmappedGI->second;

    sizeOfPackedChars += Pack< DO_PACKING >( buffer,
                                             var[li],
                                             unmappedGI.data(),
                                             var.sizeOfArray( li ),
                                             relatedObjectLocalToGlobalMap );
  }

  return sizeOfPackedChars;
}

template< bool DO_PACKING, typename SORTED >
localIndex
Pack( buffer_unit_type * & buffer,
      ArrayOfArraysView< localIndex const > const & var,
      mapBase< localIndex, SortedArray< globalIndex >, SORTED > const & unmappedGlobalIndices,
      arrayView1d< localIndex const > const & indices,
      arrayView1d< globalIndex const > const & localToGlobalMap,
      arrayView1d< globalIndex const > const & relatedObjectLocalToGlobalMap )
{
  localIndex sizeOfPackedChars=0;

  sizeOfPackedChars += Pack< DO_PACKING >( buffer, indices.size() );
  for( localIndex a=0 ; a<indices.size() ; ++a )
  {
    localIndex const li = indices[a];
    sizeOfPackedChars += Pack< DO_PACKING >( buffer, localToGlobalMap[li] );

    typename mapBase< localIndex, SortedArray< globalIndex >, SORTED >::const_iterator
      iterUnmappedGI = unmappedGlobalIndices.find( li );

    SortedArray< globalIndex > junk;
    SortedArray< globalIndex > const & unmappedGI = iterUnmappedGI==unmappedGlobalIndices.end() ?
                                                    junk :
                                                    iterUnmappedGI->second;

    sizeOfPackedChars += Pack< DO_PACKING >( buffer,
                                             var[li],
                                             unmappedGI.values(),
                                             var.sizeOfArray( li ),
                                             relatedObjectLocalToGlobalMap );
  }

  return sizeOfPackedChars;
}

template< typename SORTED0, typename SORTED1, typename SORTED2 >
inline
localIndex
Unpack( buffer_unit_type const * & buffer,
        ArrayOfArrays< localIndex > & var,
        array1d< localIndex > & indices,
        mapBase< localIndex, array1d< globalIndex >, SORTED0 > & unmappedGlobalIndices,
        mapBase< globalIndex, localIndex, SORTED1 > const & globalToLocalMap,
        mapBase< globalIndex, localIndex, SORTED2 > const & relatedObjectGlobalToLocalMap )
{
  localIndex numIndicesUnpacked;
  localIndex const sizeOfIndicesPassedIn = indices.size();

  localIndex sizeOfUnpackedChars = Unpack( buffer, numIndicesUnpacked );

  GEOS_ERROR_IF( sizeOfIndicesPassedIn!=0 && numIndicesUnpacked!=indices.size(),
                 "number of unpacked indices("<<numIndicesUnpacked<<") does not equal size of "
                                                                    "indices passed into Unpack function("<<sizeOfIndicesPassedIn );

  indices.resize( numIndicesUnpacked );

  for( localIndex a=0 ; a<indices.size() ; ++a )
  {
    globalIndex gi;
    sizeOfUnpackedChars += Unpack( buffer, gi );

    localIndex & li = indices[a];
    if( sizeOfIndicesPassedIn > 0 )
    {
      GEOS_ERROR_IF( li!=globalToLocalMap.at( gi ),
                     "global index "<<gi<<" unpacked from buffer does not equal the lookup "
                                    <<li<<" for localIndex "<<li<<" on this rank" );
    }
    else
    {
      li = globalToLocalMap.at( gi );
    }

    array1d< globalIndex > unmappedIndices;
    sizeOfUnpackedChars += Unpack( buffer,
                                   var,
                                   li,
                                   unmappedIndices,
                                   relatedObjectGlobalToLocalMap );

    if( unmappedIndices.size() > 0 )
    {
      unmappedGlobalIndices[li] = unmappedIndices;
    }
  }
  return sizeOfUnpackedChars;
}

template< bool DO_PACKING, typename SORTED >
localIndex
Pack( buffer_unit_type * & buffer,
      arrayView1d< set< localIndex > const > const & var,
      mapBase< localIndex, set< globalIndex >, SORTED > const & unmappedGlobalIndices,
      arrayView1d< localIndex const > const & indices,
      arrayView1d< globalIndex const > const & localToGlobalMap,
      arrayView1d< globalIndex const > const & relatedObjectLocalToGlobalMap )
{
  localIndex sizeOfPackedChars = Pack< DO_PACKING >( buffer, indices.size() );

  for( localIndex a=0 ; a<indices.size() ; ++a )
  {
    localIndex li = indices[a];
    sizeOfPackedChars += Pack< DO_PACKING >( buffer, localToGlobalMap[li] );

    typename mapBase< localIndex, set< globalIndex >, SORTED >::const_iterator
      iterUnmappedGI = unmappedGlobalIndices.find( li );

    set< globalIndex > junk;
    set< globalIndex > const & unmappedGI = iterUnmappedGI==unmappedGlobalIndices.end() ?
                                            junk :
                                            iterUnmappedGI->second;

    sizeOfPackedChars += Pack< DO_PACKING >( buffer,
                                             var[li],
                                             unmappedGI,
                                             relatedObjectLocalToGlobalMap );
  }

  return sizeOfPackedChars;
}


template< typename SORTED0, typename SORTED1, typename SORTED2 >
inline
localIndex
Unpack( buffer_unit_type const * & buffer,
        arrayView1d< set< localIndex > > & var,
        localIndex_array & indices,
        mapBase< localIndex, set< globalIndex >, SORTED0 > & unmappedGlobalIndices,
        mapBase< globalIndex, localIndex, SORTED1 > const & globalToLocalMap,
        mapBase< globalIndex, localIndex, SORTED2 > const & relatedObjectGlobalToLocalMap,
        bool const clearFlag )
{
  localIndex sizeOfUnpackedChars=0;
  localIndex const sizeOfIndicesPassedIn = indices.size();

  localIndex numIndicesUnpacked;
  sizeOfUnpackedChars += Unpack( buffer, numIndicesUnpacked );
  GEOS_ERROR_IF( sizeOfIndicesPassedIn!=0 && numIndicesUnpacked!=indices.size(),
                 "number of unpacked indices("<<numIndicesUnpacked<<") does not equal size of "
                                                                    "indices passed into Unpack function("<<sizeOfIndicesPassedIn );

  indices.resize( numIndicesUnpacked );

  for( localIndex a=0 ; a<numIndicesUnpacked ; ++a )
  {

    globalIndex gi;
    sizeOfUnpackedChars += Unpack( buffer, gi );

    localIndex & li = indices[a];
    if( sizeOfIndicesPassedIn > 0 )
    {
      GEOS_ERROR_IF( li!=globalToLocalMap.at( gi ),
                     "global index "<<gi<<" unpacked from buffer does equal the lookup "
                                    <<li<<" for localIndex "<<li<<" on this rank" );
    }
    else
    {
      li = globalToLocalMap.at( gi );
    }

    set< globalIndex > unmappedIndices;
    sizeOfUnpackedChars += Unpack( buffer,
                                   var[li],
                                   unmappedIndices,
                                   relatedObjectGlobalToLocalMap,
                                   clearFlag );

    unmappedGlobalIndices[li].insert( unmappedIndices.values(), unmappedIndices.size() );
  }
  return sizeOfUnpackedChars;
}

template< typename SORTED0, typename SORTED1, typename SORTED2 >
inline
localIndex
Unpack( buffer_unit_type const * & buffer,
        ArrayOfSets< localIndex > & var,
        localIndex_array & indices,
        mapBase< localIndex, SortedArray< globalIndex >, SORTED0 > & unmappedGlobalIndices,
        mapBase< globalIndex, localIndex, SORTED1 > const & globalToLocalMap,
        mapBase< globalIndex, localIndex, SORTED2 > const & relatedObjectGlobalToLocalMap,
        bool const clearFlag )
{
  localIndex sizeOfUnpackedChars=0;
  localIndex const sizeOfIndicesPassedIn = indices.size();

  localIndex numIndicesUnpacked;
  sizeOfUnpackedChars += Unpack( buffer, numIndicesUnpacked );
  GEOS_ERROR_IF( sizeOfIndicesPassedIn!=0 && numIndicesUnpacked!=indices.size(),
                 "number of unpacked indices("<<numIndicesUnpacked<<") does not equal size of "
                                                                    "indices passed into Unpack function("<<sizeOfIndicesPassedIn );

  indices.resize( numIndicesUnpacked );

  for( localIndex a=0 ; a<numIndicesUnpacked ; ++a )
  {

    globalIndex gi;
    sizeOfUnpackedChars += Unpack( buffer, gi );

    localIndex & li = indices[a];
    if( sizeOfIndicesPassedIn > 0 )
    {
      GEOS_ERROR_IF( li!=globalToLocalMap.at( gi ),
                     "global index "<<gi<<" unpacked from buffer does equal the lookup "
                                    <<li<<" for localIndex "<<li<<" on this rank" );
    }
    else
    {
      li = globalToLocalMap.at( gi );
    }

    SortedArray< globalIndex > unmappedIndices;
    sizeOfUnpackedChars += Unpack( buffer,
                                   var,
                                   li,
                                   unmappedIndices,
                                   relatedObjectGlobalToLocalMap,
                                   clearFlag );

    unmappedGlobalIndices[li].insert( unmappedIndices.values(), unmappedIndices.size() );
  }
  return sizeOfUnpackedChars;
}


template< bool DO_PACKING, int UNIT_STRIDE_DIM0, int UNIT_STRIDE_DIM1 >
localIndex
Pack( buffer_unit_type * & buffer,
      arrayView2d< localIndex const, UNIT_STRIDE_DIM0 > const & var,
      arrayView1d< localIndex > const & indices,
      arraySlice1d< globalIndex const, UNIT_STRIDE_DIM1 > const & localToGlobalMap )
{
  localIndex sizeOfPackedChars = Pack< DO_PACKING >( buffer, indices.size() );
  for( localIndex a=0 ; a<indices.size() ; ++a )
  {
    localIndex li = indices[a];
    sizeOfPackedChars += Pack< DO_PACKING >( buffer, localToGlobalMap[li] );

    sizeOfPackedChars += Pack< DO_PACKING >( buffer, var[li], var.size( 1 ) );
  }

  return sizeOfPackedChars;
}

template< typename SORTED, int UNIT_STRIDE_DIM >
inline
localIndex
Unpack( buffer_unit_type const * & buffer,
        arrayView2d< localIndex, UNIT_STRIDE_DIM > const & var,
        array1d< localIndex > & indices,
        mapBase< globalIndex, localIndex, SORTED > const & globalToLocalMap )
{
  localIndex sizeOfUnpackedChars = 0;
  localIndex const sizeOfIndicesPassedIn = indices.size();

  localIndex numIndicesUnpacked;
  sizeOfUnpackedChars += Unpack( buffer, numIndicesUnpacked );

  GEOS_ERROR_IF( sizeOfIndicesPassedIn!=0 && numIndicesUnpacked!=indices.size(),
                 "number of unpacked indices("<<numIndicesUnpacked<<") does not equal size of "
                                                                    "indices passed into Unpack function("<<sizeOfIndicesPassedIn );

  indices.resize( numIndicesUnpacked );

  for( localIndex a=0 ; a<numIndicesUnpacked ; ++a )
  {
    globalIndex gi;
    sizeOfUnpackedChars += Unpack( buffer, gi );

    localIndex & li = indices[a];
    if( sizeOfIndicesPassedIn > 0 )
    {
      GEOS_ERROR_IF( li!=globalToLocalMap.at( gi ),
                     "global index "<<gi<<" unpacked from buffer does equal the lookup "
                                    <<li<<" for localIndex "<<li<<" on this rank" );
    }
    else
    {
      li = globalToLocalMap.at( gi );
    }

    localIndex * const varSlice = var[li];
    sizeOfUnpackedChars += Unpack( buffer, varSlice, var.size( 1 ) );
  }


  return sizeOfUnpackedChars;
}


template< bool DO_PACKING, typename SORTED, int UNIT_STRIDE_DIM0 >
localIndex
Pack( buffer_unit_type * & buffer,
      arrayView2d< localIndex const, UNIT_STRIDE_DIM0 > const & var,
      mapBase< localIndex, array1d< globalIndex >, SORTED > const & unmappedGlobalIndices,
      arrayView1d< localIndex const > const & indices,
      arraySlice1d< globalIndex const > const & localToGlobalMap,
      arraySlice1d< globalIndex const > const & relatedObjectLocalToGlobalMap )
{
  localIndex sizeOfPackedChars = 0;

  sizeOfPackedChars += Pack< DO_PACKING >( buffer, indices.size() );
  for( localIndex a=0 ; a<indices.size() ; ++a )
  {
    localIndex li = indices[a];
    sizeOfPackedChars += Pack< DO_PACKING >( buffer, localToGlobalMap[li] );

    typename mapBase< localIndex, array1d< globalIndex >, SORTED >::const_iterator
      iterUnmappedGI = unmappedGlobalIndices.find( li );

    array1d< globalIndex > junk;
    array1d< globalIndex > const & unmappedGI = iterUnmappedGI==unmappedGlobalIndices.end() ?
                                                junk :
                                                iterUnmappedGI->second;

    sizeOfPackedChars += Pack< DO_PACKING >( buffer,
                                             var[li],
                                             unmappedGI.data(),
                                             var.size( 1 ),
                                             relatedObjectLocalToGlobalMap );
  }

  return sizeOfPackedChars;
}


template< typename SORTED0, typename SORTED1, typename SORTED2, int UNIT_STRIDE_DIM >
inline
localIndex
Unpack( buffer_unit_type const * & buffer,
        arrayView2d< localIndex, UNIT_STRIDE_DIM > const & var,
        localIndex_array & indices,
        mapBase< localIndex, array1d< globalIndex >, SORTED0 > & unmappedGlobalIndices,
        mapBase< globalIndex, localIndex, SORTED1 > const & globalToLocalMap,
        mapBase< globalIndex, localIndex, SORTED2 > const & relatedObjectGlobalToLocalMap )
{
  localIndex sizeOfUnpackedChars = 0;
  localIndex const sizeOfIndicesPassedIn = indices.size();

  localIndex numIndicesUnpacked;
  sizeOfUnpackedChars += Unpack( buffer, numIndicesUnpacked );
  GEOS_ERROR_IF( sizeOfIndicesPassedIn!=0 && numIndicesUnpacked!=indices.size(),
                 "number of unpacked indices("<<numIndicesUnpacked<<") does not equal size of "
                                                                    "indices passed into Unpack function("<<sizeOfIndicesPassedIn );

  indices.resize( numIndicesUnpacked );

  for( localIndex a=0 ; a<numIndicesUnpacked ; ++a )
  {
    globalIndex gi;
    sizeOfUnpackedChars += Unpack( buffer, gi );

    localIndex & li = indices[a];
    if( sizeOfIndicesPassedIn > 0 )
    {
      GEOS_ERROR_IF( li!=globalToLocalMap.at( gi ),
                     "global index "<<gi<<" unpacked from buffer does equal the lookup "
                                    <<li<<" for localIndex "<<li<<" on this rank" );
    }
    else
    {
      li = globalToLocalMap.at( gi );
    }

    arraySlice1d< localIndex, UNIT_STRIDE_DIM - 1 > varSlice = var[li];
    array1d< globalIndex > unmappedIndices;

    sizeOfUnpackedChars += Unpack( buffer,
                                   varSlice,
                                   unmappedIndices,
                                   var.size( 1 ),
                                   relatedObjectGlobalToLocalMap );

    if( unmappedIndices.size()>0 )
    {
      unmappedGlobalIndices[li] = unmappedIndices;
    }

  }

  return sizeOfUnpackedChars;
}


template< bool DO_PACKING, typename MAP_TYPE >
typename std::enable_if< is_packable_map< MAP_TYPE >, localIndex >::type
Pack( buffer_unit_type * & buffer, MAP_TYPE const & var )
{
  const typename MAP_TYPE::size_type length = var.size();
  localIndex sizeOfPackedChars = Pack< DO_PACKING >( buffer, length );

  for( typename MAP_TYPE::const_iterator i=var.begin() ; i!=var.end() ; ++i )
  {
    sizeOfPackedChars += Pack< DO_PACKING >( buffer, i->first );
    sizeOfPackedChars += Pack< DO_PACKING >( buffer, i->second );
  }

  return sizeOfPackedChars;
}


template< typename MAP_TYPE >
typename std::enable_if< is_packable_map< MAP_TYPE >, localIndex >::type
Unpack( buffer_unit_type const * & buffer, MAP_TYPE & map )
{
  map.clear();
  typename MAP_TYPE::size_type map_length;

  localIndex sizeOfUnpackedChars = Unpack( buffer, map_length );

  for( typename MAP_TYPE::size_type a=0 ; a<map_length ; ++a )
  {
    typename MAP_TYPE::key_type key;
    typename MAP_TYPE::mapped_type value;
    sizeOfUnpackedChars += Unpack( buffer, key );
    sizeOfUnpackedChars += Unpack( buffer, value );

    map[ key ] = std::move( value );
  }

  return sizeOfUnpackedChars;
}


template< bool DO_PACKING, typename MAP_TYPE, typename T_INDICES >
typename std::enable_if< is_map_packable_by_index< MAP_TYPE >, localIndex >::type
Pack( buffer_unit_type * & buffer, MAP_TYPE const & var, T_INDICES const & packIndices )
{
  typename MAP_TYPE::size_type const length = var.size();
  localIndex sizeOfPackedChars = Pack< DO_PACKING >( buffer, length );

  for( typename MAP_TYPE::const_iterator i=var.begin() ; i!=var.end() ; ++i )
  {
    sizeOfPackedChars += Pack< DO_PACKING >( buffer, i->first );
    sizeOfPackedChars += Pack< DO_PACKING >( buffer, i->second, packIndices );
  }

  return sizeOfPackedChars;
}


template< typename MAP_TYPE, typename T_INDICES >
typename std::enable_if< is_map_packable_by_index< MAP_TYPE >, localIndex >::type
Unpack( buffer_unit_type const * & buffer, MAP_TYPE & map, T_INDICES const & unpackIndices )
{
  map.clear();
  typename MAP_TYPE::size_type map_length;

  localIndex sizeOfUnpackedChars = Unpack( buffer, map_length );

  for( typename MAP_TYPE::size_type a=0 ; a<map_length ; ++a )
  {
    typename MAP_TYPE::key_type key;
    typename MAP_TYPE::mapped_type value;
    sizeOfUnpackedChars += Unpack( buffer, key );
    sizeOfUnpackedChars += Unpack( buffer, value, unpackIndices );

    map[key] = value;
  }

  return sizeOfUnpackedChars;
}


template< bool DO_PACKING, typename T_FIRST, typename T_SECOND >
localIndex
Pack( buffer_unit_type * & buffer, std::pair< T_FIRST, T_SECOND > const & var )
{
  localIndex sizeOfPackedChars = Pack< DO_PACKING >( buffer, var.first );
  sizeOfPackedChars += Pack< DO_PACKING >( buffer, var.second );
  return sizeOfPackedChars;
}


template< typename T_FIRST, typename T_SECOND >
localIndex
Unpack( buffer_unit_type const * & buffer, std::pair< T_FIRST, T_SECOND > & var )
{
  localIndex sizeOfUnpackedChars = Unpack( buffer, var.first );
  sizeOfUnpackedChars += Unpack( buffer, var.second );
  return sizeOfUnpackedChars;
}

} /* namespace bufferOps */
} /* namespace geosx */
