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
 * @file WrapperHelpers.hpp
 */

#ifndef GEOSX_DATAREPOSITORY_WRAPPERHELPERS_HPP_
#define GEOSX_DATAREPOSITORY_WRAPPERHELPERS_HPP_

/**
 * @brief Enables verbose logging of restart output
 */
#define RESTART_TYPE_LOGGING 0

// Source includes
#include "BufferOps.hpp"
#include "BufferOpsDevice.hpp"
#include "DefaultValue.hpp"
#include "ConduitRestart.hpp"
#include "common/DataTypes.hpp"
#include "common/GeosxMacros.hpp"
#include "codingUtilities/traits.hpp"

// TPL includes
#include <conduit.hpp>

// System includes
#include <cstring>

#if RESTART_TYPE_LOGGING
#include <unordered_set>
#endif

namespace geosx
{
namespace dataRepository
{
namespace wrapperHelpers
{

inline void logOutputType( std::string const & typeString, std::string const & msg )
{
#if RESTART_TYPE_LOGGING
  static std::unordered_set< std::string > m_types;

  if( !m_types.count( typeString ) )
  {
    m_types.insert( typeString );
    GEOSX_LOG( msg << typeString );
  }
#else
  GEOSX_DEBUG_VAR( typeString );
  GEOSX_DEBUG_VAR( msg );
#endif
}

template< typename T >
inline std::enable_if_t< traits::has_size_method< T >, localIndex >
size( T const & value )
{
  return integer_conversion< localIndex >( value.size() );
}

template< typename T >
inline std::enable_if_t< !traits::has_size_method< T >, localIndex >
size( T const & GEOSX_UNUSED_ARG( value ) )
{
  return 1;
}

template< typename T >
inline std::enable_if_t< traits::is_string< T >, char * >
dataPtr( T & var )
{
  return const_cast< char * >( var.data() );
}

template< typename T >
inline std::enable_if_t< !traits::is_string< T > && traits::has_data_method< T >, typename traits::Pointer< T > >
dataPtr( T & value )
{
  return value.data();
}

template< typename T >
inline std::enable_if_t< !traits::is_string< T > && !traits::has_data_method< T >, typename traits::Pointer< T > >
dataPtr( T & value )
{
  return &value;
}

template< class T >
inline T const *
dataPtr( SortedArray< T > & value )
{
  return value.values();
}

template< class T >
inline typename traits::ConstPointer< T >
dataPtr( T const & value )
{
  return dataPtr( const_cast< T & >( value ) );
}

template< typename T >
inline std::enable_if_t< traits::has_resize_method< T > >
resize( T & value, localIndex const newSize )
{
  value.resize( newSize );
}

template< typename T >
inline std::enable_if_t< !traits::has_resize_method< T > >
resize( T & GEOSX_UNUSED_ARG( value ),
        localIndex const GEOSX_UNUSED_ARG( newSize ) )
{}

template< typename T >
inline std::enable_if_t< traits::has_resize_default_method< T, typename DefaultValue< T >::value_type > >
resizeDefault( T & value, localIndex const newSize, DefaultValue< T > const & defaultValue )
{
  value.resizeDefault( newSize, defaultValue.value );
}

template< typename T >
inline std::enable_if_t< !traits::has_resize_default_method< T, typename DefaultValue< T >::value_type > >
resizeDefault( T & value, localIndex const newSize, DefaultValue< T > const & GEOSX_UNUSED_ARG( defaultValue ) )
{
  resize( value, newSize );
}

template< typename T >
inline std::enable_if_t< traits::has_resize_dimensions_method< T > >
resizeDimensions( T & value, int num_dims, localIndex const * const dims )
{
  value.resize( num_dims, dims );
}

template< typename T >
inline std::enable_if_t< !traits::has_resize_dimensions_method< T > >
resizeDimensions( T & value, int num_dims, localIndex const * const dims )
{
  if( num_dims != 1 )
  {
    GEOSX_ERROR( "Data is only 1D" );
    return;
  }
  resize( value, dims[ 0 ] );
}

template< typename T >
inline std::enable_if_t< traits::has_alias_value_type< T >, localIndex >
byteSize( T const & value )
{
  return size( value ) * sizeof( typename T::value_type );
}

template< typename T >
inline std::enable_if_t< !traits::has_alias_value_type< T >, localIndex >
byteSize( T const & value )
{
  return size( value ) * sizeof( T );
}

template< typename T >
inline std::enable_if_t< traits::has_alias_value_type< T >, localIndex >
numElementsFromByteSize( localIndex const byteSize )
{
  GEOSX_ERROR_IF_NE( byteSize % sizeof( typename T::value_type ), 0 );
  return byteSize / sizeof( typename T::value_type );
}


template< typename T >
inline std::enable_if_t< !traits::has_alias_value_type< T >, localIndex >
numElementsFromByteSize( localIndex const byteSize )
{
  GEOSX_ERROR_IF_NE( byteSize % sizeof( T ), 0 );
  return byteSize / sizeof( T );
}


// This is for an object that needs to be packed.
template< typename T >
std::enable_if_t< !bufferOps::can_memcpy< typename traits::Pointer< T > > >
pushDataToConduitNode( T const & var, conduit::Node & node )
{
  logOutputType( cxx_utilities::demangle< decltype( var ) >(), "Packing for output: " );

  // Get the number of bytes in the packed object.
  localIndex const byteSize = bufferOps::PackSize( var );

  // Create a conduit data type that describes the array.
  conduit::DataType const dtype( conduitTypeInfo< buffer_unit_type >::id, byteSize );

  // Allocate the array in the "__values__" child.
  conduit::Node & valuesNode = node[ "__values__" ];
  valuesNode.set( dtype );

  // Get the pointer to the array and pack the object into it.
  buffer_unit_type * buffer = valuesNode.value();
  bufferOps::Pack< true >( buffer, var );
}

// This is for an object that needs to be packed.
template< typename T >
std::enable_if_t< !bufferOps::can_memcpy< typename traits::Pointer< T > > >
pullDataFromConduitNode( T & var, conduit::Node const & node )
{
  conduit::Node const & valuesNode = node.fetch_child( "__values__" );

  // Get the number of bytes in the array and a pointer to the array.
  localIndex const byteSize = valuesNode.dtype().number_of_elements();
  buffer_unit_type const * buffer = valuesNode.value();

  // Unpack the object from the array.
  localIndex const bytesRead = bufferOps::Unpack( buffer, var );
  GEOSX_ERROR_IF_NE( bytesRead, byteSize );
}

// This is for an std::string since the type of char is different on different platforms :(.
inline
void
pushDataToConduitNode( std::string const & var, conduit::Node & node )
{
  logOutputType( cxx_utilities::demangle< decltype( var ) >(), "Output via external pointer: " );

  constexpr int conduitTypeID = conduitTypeInfo< signed char >::id;
  conduit::DataType const dtype( conduitTypeID, var.size() );

  signed char * const ptr = const_cast< signed char * >( reinterpret_cast< signed char const * >( var.data() ) );
  node[ "__values__" ].set_external( dtype, ptr );
}

// This is for an object that doesn't need to be packed but isn't an LvArray.
template< typename T >
std::enable_if_t< bufferOps::can_memcpy< typename traits::Pointer< T > > >
pushDataToConduitNode( T const & var, conduit::Node & node )
{
  logOutputType( cxx_utilities::demangle< decltype( var ) >(), "Output via external pointer: " );

  constexpr int conduitTypeID = conduitTypeInfo< typename traits::Pointer< T > >::id;
  constexpr int sizeofConduitType = conduitTypeInfo< typename traits::Pointer< T > >::sizeOfConduitType;
  localIndex const numBytes = byteSize( var );
  conduit::DataType const dtype( conduitTypeID, numBytes / sizeofConduitType );

  void * const ptr = const_cast< void * >( static_cast< void const * >( dataPtr( var ) ) );
  node[ "__values__" ].set_external( dtype, ptr );
}

// This is for an object that doesn't need to be packed but isn't an LvArray or a SortedArray.
template< typename T >
std::enable_if_t< bufferOps::can_memcpy< typename traits::Pointer< T > > >
pullDataFromConduitNode( T & var, conduit::Node const & node )
{
  conduit::Node const & valuesNode = node.fetch_child( "__values__" );

  localIndex const byteSize = integer_conversion< localIndex >( valuesNode.dtype().strided_bytes() );
  localIndex const numElements = numElementsFromByteSize< T >( byteSize );

  resize( var, numElements );

  std::memcpy( dataPtr( var ), valuesNode.data_ptr(), byteSize );
}

// This is for a SortedArray that doesn't need to be packed.
template< typename T >
std::enable_if_t< bufferOps::can_memcpy< T > >
pullDataFromConduitNode( SortedArray< T > & var, conduit::Node const & node )
{
  conduit::Node const & valuesNode = node.fetch_child( "__values__" );

  localIndex const byteSize = integer_conversion< localIndex >( valuesNode.dtype().strided_bytes() );
  localIndex const numElements = numElementsFromByteSize< T >( byteSize );

  T const * const values = reinterpret_cast< T const * >( valuesNode.data_ptr() );
  var.insertSorted( values, numElements );
}


// This is an LvArray that doesn't need to be packed.
template< typename T, int NDIM, typename PERMUTATION >
std::enable_if_t< bufferOps::can_memcpy< T > >
pushDataToConduitNode( Array< T, NDIM, PERMUTATION > const & var,
                       conduit::Node & node )
{
  logOutputType( cxx_utilities::demangle< decltype( var ) >(), "Output array via external pointer: " );

  // Push the data into conduit
  constexpr int conduitTypeID = conduitTypeInfo< T >::id;
  constexpr int sizeofConduitType = conduitTypeInfo< T >::sizeOfConduitType;
  conduit::DataType const dtype( conduitTypeID, var.size() * sizeof( T ) / sizeofConduitType );
  void * const ptr = const_cast< void * >( static_cast< void const * >( var.data() ) );
  node[ "__values__" ].set_external( dtype, ptr );

  // Create a copy of the dimensions
  localIndex temp[ NDIM + 1 ];
  for( int i = 0 ; i < NDIM ; ++i )
  {
    temp[ i ] = var.size( i );
  }

  // If T is something like a Tensor than there is an extra implicit dimension.
  constexpr int const implicitDimensionLength = conduitTypeInfo< T >::numConduitValues;
  constexpr bool const hasImplicitDimension = implicitDimensionLength != 1;
  constexpr int totalNumDimensions = NDIM + hasImplicitDimension;
  if( hasImplicitDimension )
  {
    temp[ NDIM ] = implicitDimensionLength;
  }

  // push the dimensions into the node
  conduit::DataType const dimensionType( conduitTypeInfo< localIndex >::id, totalNumDimensions );
  node[ "__dimensions__" ].set( dimensionType, temp );

  // Create a copy of the permutation
  constexpr std::array< camp::idx_t, NDIM > const perm = RAJA::as_array< PERMUTATION >::get();
  for( int i = 0 ; i < NDIM ; ++i )
  {
    temp[ i ] = perm[ i ];
  }

  if( hasImplicitDimension )
  {
    temp[ NDIM ] = NDIM;
  }

  node[ "__permutation__" ].set( dimensionType, temp );
}

// This is an LvArray that doesn't need to be packed.
template< typename T, int NDIM, typename PERMUTATION >
std::enable_if_t< bufferOps::can_memcpy< T > >
pullDataFromConduitNode( Array< T, NDIM, PERMUTATION > & var,
                         conduit::Node const & node )
{
  // Get the number of dimensions written out, accounting for an implicit dimension and the permutation.
  constexpr int const implicitDimensionLength = conduitTypeInfo< T >::numConduitValues;
  constexpr bool const hasImplicitDimension = implicitDimensionLength != 1;
  constexpr int totalNumDimensions = NDIM + hasImplicitDimension;

  // Check that the permutations match.
  conduit::Node const & permutationNode = node.fetch_child( "__permutation__" );
  GEOSX_ERROR_IF_NE( permutationNode.dtype().number_of_elements(), totalNumDimensions );

  constexpr std::array< camp::idx_t, NDIM > const perm = RAJA::as_array< PERMUTATION >::get();
  camp::idx_t const * const permFromConduit = permutationNode.value();
  for( int i = 0 ; i < NDIM ; ++i )
  {
    GEOSX_ERROR_IF_NE_MSG( permFromConduit[ i ], perm[ i ],
                           "The permutation of the data in conduit and the provided Array don't match." );
  }

  if( hasImplicitDimension )
  {
    GEOSX_ERROR_IF_NE_MSG( permFromConduit[ NDIM ], NDIM,
                           "The permutation of the data in conduit and the provided Array don't match." );
  }

  // Now pull out the dimensions and resize the array.
  conduit::Node const & dimensionNode = node.fetch_child( "__dimensions__" );
  GEOSX_ERROR_IF_NE( dimensionNode.dtype().number_of_elements(), totalNumDimensions );
  localIndex const * const dims = dimensionNode.value();

  if( hasImplicitDimension )
  {
    GEOSX_ERROR_IF_NE( dims[ NDIM ], implicitDimensionLength );
  }

  var.resize( NDIM, dims );

  // Finally memcpy
  conduit::Node const & valuesNode = node.fetch_child( "__values__" );
  localIndex numBytesFromArray =  var.size() * sizeof( T );
  GEOSX_ERROR_IF_NE( numBytesFromArray, valuesNode.dtype().strided_bytes() );
  std::memcpy( var.data(), valuesNode.data_ptr(), numBytesFromArray );
}

template< typename T >
void pushDataToConduitNode( InterObjectRelation< T > const & var,
                            conduit::Node & node )
{
  return pushDataToConduitNode( var.Base(), node );
}

template< typename T >
void pullDataFromConduitNode( InterObjectRelation< T > & var,
                              conduit::Node const & node )
{
  return pullDataFromConduitNode( var.Base(), node );
}

template< bool DO_PACKING, typename T, typename IDX >
inline std::enable_if_t< bufferOps::is_packable_by_index< T >, localIndex >
PackByIndex( buffer_unit_type * & buffer, T & var, IDX & idx )
{
  return bufferOps::PackByIndex< DO_PACKING >( buffer, var, idx );
}

template< bool DO_PACKING, typename T, typename IDX >
inline std::enable_if_t< !bufferOps::is_packable_by_index< T >, localIndex >
PackByIndex( buffer_unit_type * &, T &, IDX & )
{
  return 0;
}

template< typename T, typename IDX >
inline std::enable_if_t< bufferOps::is_packable_by_index< T >, localIndex >
UnpackByIndex( buffer_unit_type const * & buffer, T & var, IDX & idx )
{
  return bufferOps::UnpackByIndex( buffer, var, idx );
}

template< typename T, typename IDX >
inline std::enable_if_t< !bufferOps::is_packable_by_index< T >, localIndex >
UnpackByIndex( buffer_unit_type const * &, T &, IDX & )
{
  return 0;
}


template< bool DO_PACKING, typename T >
inline std::enable_if_t< bufferOps::is_container< T > && bufferOps::can_memcpy< T >, localIndex >
PackDevice( buffer_unit_type * & buffer, T & var )
{
  return bufferOps::PackDevice< parallelDevicePolicy< >, DO_PACKING >( buffer, var );
}

template< bool DO_PACKING, typename T >
inline std::enable_if_t< !bufferOps::is_container< T > || !bufferOps::can_memcpy< T >, localIndex >
PackDevice( buffer_unit_type * &, T & )
{
  GEOSX_ERROR( "Trying to pack data type ("<<typeid(T).name()<<") on device but type is not packable on device." );
  return 0;
}

template< bool DO_PACKING, typename T, typename IDX >
inline std::enable_if_t< bufferOps::is_container< T >, localIndex >
PackByIndexDevice( buffer_unit_type * & buffer, T & var, IDX & idx )
{
  return bufferOps::PackByIndexDevice< parallelDevicePolicy< >, DO_PACKING >( buffer, var, idx );
}

template< bool DO_PACKING, typename T, typename IDX >
inline std::enable_if_t< !bufferOps::is_container< T >, localIndex >
PackByIndexDevice( buffer_unit_type * &, T &, IDX & )
{
  GEOSX_ERROR( "Trying to pack data type ("<<typeid(T).name()<<") on device but type is not packable by index." );
  return 0;
}

template< typename T >
inline std::enable_if_t< bufferOps::is_container< T >, localIndex >
UnpackDevice( buffer_unit_type const * & buffer, T & var )
{
  return bufferOps::UnpackDevice< parallelDevicePolicy< > >( buffer, var );
}

template< typename T >
inline std::enable_if_t< !bufferOps::is_container< T >, localIndex >
UnpackDevice( buffer_unit_type const * &, T & )
{
  return 0;
}

template< typename T, typename IDX >
inline std::enable_if_t< bufferOps::is_container< T >, localIndex >
UnpackByIndexDevice( buffer_unit_type const * & buffer, T & var, IDX & idx )
{
  return bufferOps::UnpackByIndexDevice< parallelDevicePolicy< > >( buffer, var, idx );
}

template< typename T, typename IDX >
inline std::enable_if_t< !bufferOps::is_container< T >, localIndex >
UnpackByIndexDevice( buffer_unit_type const * &, T &, IDX & )
{
  return 0;
}

} // namespace WrapperHelpers
} // namespace dataRepository
} // namespace geosx

#undef RESTART_TYPE_LOGGING

#endif // GEOSX_DATAREPOSITORY_WRAPPERHELPERS_HPP_
