/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file HDF5Utilities.cpp
 */

#include "HDF5Utilities.hpp"
#include "HDF5Utilities_impl.hpp"
#include "common/format/Format.hpp"
#include "common/logger/Logger.hpp"

#include <variant>

#include <algorithm>

#include <hdf5.h>

#define GEOS_HDF5_CHECK_ERROR_WITH_THRESHOLD( call, context, type, threshold ) \
  do {                                                                         \
    type __geos_hdf_internal_result__ = -1;                                    \
    H5E_BEGIN_TRY {                                                            \
      __geos_hdf_internal_result__ = (call);                                   \
    } H5E_END_TRY;                                                             \
    if( __geos_hdf_internal_result__ < (threshold) ) {                         \
      H5Eclear2( H5E_DEFAULT );                                                \
      GEOS_THROW( GEOS_FMT( "Error in call to:\n"                              \
                            "{}\n"                                             \
                            "({})",                                            \
                            #call, context ),                                  \
                  InputError );                                                \
    }                                                                          \
  } while (false)

#define GEOS_HDF5_CHECK_WITH_THRESHOLD_AND_ASSIGN( var, call, context, type, threshold ) \
  do {                                                                                   \
    type __geos_hdf_internal_result__ = -1;                                              \
    H5E_BEGIN_TRY {                                                                      \
      __geos_hdf_internal_result__ = (call);                                             \
    } H5E_END_TRY;                                                                       \
    if( __geos_hdf_internal_result__ < (threshold) ) {                                   \
      H5Eclear2( H5E_DEFAULT );                                                          \
      GEOS_THROW( GEOS_FMT( "Error in call to:\n"                                        \
                            "{}\n"                                                       \
                            "({})",                                                      \
                            #call, context ),                                            \
                  InputError );                                                          \
    }                                                                                    \
    var = __geos_hdf_internal_result__;                                                  \
  } while (false)

#define GEOS_HDF5_CHECK_ERROR( call, context ) \
  GEOS_HDF5_CHECK_ERROR_WITH_THRESHOLD( call, context, herr_t, 0 )

#define GEOS_HDF5_CHECK_AND_ASSIGN_HID( var, call, context ) \
  GEOS_HDF5_CHECK_WITH_THRESHOLD_AND_ASSIGN( var, call, context, hid_t, 0 )

#define GEOS_HDF5_CHECK_AND_ASSIGN_INT( var, call, context ) \
  GEOS_HDF5_CHECK_WITH_THRESHOLD_AND_ASSIGN( var, call, context, int, 0 )

namespace geos
{

namespace hdf5Utils
{

namespace
{
enum class HDF5NativeDataType
{
  UInt8,
  Int,
  Float,
  Double,
  Unknown
};

inline HDF5NativeDataType resolveDataType( hid_t typeId )
{
  if( H5Tequal( typeId, H5T_NATIVE_UINT8 ))
  {
    return HDF5NativeDataType::UInt8;
  }
  if( H5Tequal( typeId, H5T_NATIVE_INT ))
  {
    return HDF5NativeDataType::Int;
  }
  if( H5Tequal( typeId, H5T_NATIVE_FLOAT ))
  {
    return HDF5NativeDataType::Float;
  }
  if( H5Tequal( typeId, H5T_NATIVE_DOUBLE ))
  {
    return HDF5NativeDataType::Double;
  }
  return HDF5NativeDataType::Unknown;
}

template< typename T >
array1d< T > readTypedData( hid_t datasetId,
                            hid_t spaceId,
                            hid_t nativeType,
                            hsize_t total_elements,
                            string const & datasetName,
                            string const & filename )
{
  string contextCheckMessage{ GEOS_FMT( "Dataset {} in {}", datasetName, filename ) };
  array1d< T > buffer( total_elements );

  GEOS_HDF5_CHECK_ERROR( H5Dread( datasetId, nativeType, H5S_ALL, spaceId, H5P_DEFAULT, buffer.data()),
                         contextCheckMessage );

  return buffer;
}

using TypedArray1d = std::variant<
  array1d< uint8_t >,
  array1d< localIndex >,
  array1d< real32 >,
  array1d< real64 > >;

static TypedArray1d read1D( SerialHDF5File const & file,
                            string const & datasetName,
                            int const expectedDims )
{
  // Create a DatasetHandle to manage resources
  DatasetHandle handle( file, datasetName, expectedDims );

  // Compute the total number of elements
  hsize_t total_elements = 1;
  for( const auto & dim : handle.dims )
  {
    total_elements *= dim;
  }

  // Determine the type and read the data
  switch( resolveDataType( handle.typeId ) )
  {
    case HDF5NativeDataType::UInt8:
    {
      return readTypedData< uint8_t >( handle.datasetId, handle.spaceId, H5T_NATIVE_UINT8, total_elements, datasetName, file.getFilename() );
    }
    case HDF5NativeDataType::Int:
    {
      return readTypedData< localIndex >( handle.datasetId, handle.spaceId, H5T_NATIVE_INT, total_elements, datasetName, file.getFilename() );
    }
    case HDF5NativeDataType::Float:
    {
      return readTypedData< real32 >( handle.datasetId, handle.spaceId, H5T_NATIVE_FLOAT, total_elements, datasetName, file.getFilename() );
    }
    case HDF5NativeDataType::Double:
    {
      return readTypedData< real64 >( handle.datasetId, handle.spaceId, H5T_NATIVE_DOUBLE, total_elements, datasetName, file.getFilename() );
    }
    default:
    {
      GEOS_THROW( GEOS_FMT( "Unsupported dataset type for dataset {} in {}", datasetName, file.getFilename()), InputError );
    }
  }
}



template< typename SOURCE_T, typename TARGET_T >
array1d< TARGET_T > staticCastArray( array1d< SOURCE_T > const & source )
{
  array1d< TARGET_T > casted( source.size() );
  std::transform( source.begin(), source.end(), casted.begin(),
                  []( auto value ) { return static_cast< TARGET_T >( value ); } );
  return casted;
}

} // end anonymous namespace

static_assert( H5_SIZEOF_UINT8_T == sizeof( uint8_t ),
               "H5T_NATIVE_UINT8 and uint8_t must have the same size" );
static_assert( H5_SIZEOF_INT == sizeof( localIndex ),
               "H5T_NATIVE_INT and geos::integer must have the same size" );
static_assert( H5_SIZEOF_FLOAT == sizeof( real32 ),
               "H5T_NATIVE_FLOAT and geos::real32 must have the same size" );
static_assert( H5_SIZEOF_DOUBLE == sizeof( real64 ),
               "H5T_NATIVE_DOUBLE and geos::real64 must have the same size" );

//SerialHDF5File Implementation
SerialHDF5File::SerialHDF5File( const string & filename ): m_filename( filename )
{
  openFile();
}

SerialHDF5File::~SerialHDF5File()
{
  closeFile();
}

SerialHDF5File::SerialHDF5File( SerialHDF5File && other ) noexcept: m_fileId( other.m_fileId ), m_filename( std::move( other.m_filename ))
{
  other.m_fileId = -1;
}

SerialHDF5File & SerialHDF5File::operator=( SerialHDF5File && other ) noexcept
{
  if( this != &other )
  {
    closeFile();
    m_fileId = other.m_fileId;
    m_filename = std::move( other.m_filename );
    other.m_fileId = -1;
  }
  return *this;
}

void SerialHDF5File::openFile()
{
  closeFile();

  GEOS_HDF5_CHECK_AND_ASSIGN_HID( m_fileId,
                                  H5Fopen( m_filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT ),
                                  GEOS_FMT( "Filename: {}", m_filename ) );
}

void SerialHDF5File::closeFile()
{
  if( m_fileId >= 0 )
  {
    GEOS_HDF5_CHECK_ERROR( H5Fclose( m_fileId ), GEOS_FMT( "Filename: {}", m_filename ) );
    m_fileId = -1;
  }
}

// DatasetHandle Implementation

DatasetHandle::DatasetHandle( SerialHDF5File const & hdf5File, string const & datasetName, int const expectedDims )
  : m_datasetName( datasetName )
{
  string contextCheckMessage{ GEOS_FMT( "Dataset {} in {}", datasetName, hdf5File.getFilename() ) };

  if( datasetExists( hdf5File.getFileId(), datasetName )  )
  {
    GEOS_HDF5_CHECK_AND_ASSIGN_HID( datasetId,
                                    H5Dopen2( hdf5File.getFileId(), datasetName.c_str(), H5P_DEFAULT ),
                                    contextCheckMessage );

    GEOS_HDF5_CHECK_AND_ASSIGN_HID( spaceId,
                                    H5Dget_space( datasetId ),
                                    contextCheckMessage );

    GEOS_HDF5_CHECK_AND_ASSIGN_HID( typeId,
                                    H5Dget_type( datasetId ),
                                    contextCheckMessage );

    int ndims{};
    GEOS_HDF5_CHECK_AND_ASSIGN_INT( ndims,
                                    H5Sget_simple_extent_ndims( spaceId ),
                                    contextCheckMessage );

    // Validate dimensions
    GEOS_THROW_IF( ndims != expectedDims,
                   GEOS_FMT( "Incosistent number of dimensions for dataset {} in {}", datasetName, hdf5File.getFilename() ),
                   InputError );

    dims.resize( ndims );
    GEOS_HDF5_CHECK_ERROR( H5Sget_simple_extent_dims( spaceId, dims.data(), nullptr ),
                           contextCheckMessage );
  }
  else
  {
    GEOS_THROW( GEOS_FMT( "Dataset {} doesn't exist in {}", datasetName, hdf5File.getFilename() ),
                InputError );
  }
}

DatasetHandle::~DatasetHandle()
{
  if( typeId >= 0 )
  {
    H5E_BEGIN_TRY
    {
      H5Tclose( typeId );
    }
    H5E_END_TRY
  }
  if( spaceId >= 0 )
  {
    H5E_BEGIN_TRY
    {
      H5Sclose( spaceId );
    }
    H5E_END_TRY
  }
  if( datasetId >= 0 )
  {
    H5E_BEGIN_TRY
    {
      H5Dclose( datasetId );
    }
    H5E_END_TRY
  }
}

DatasetHandle::DatasetHandle( DatasetHandle && other ) noexcept
  : datasetId( other.datasetId ), spaceId( other.spaceId ), typeId( other.typeId ), dims( std::move( other.dims ))
{
  other.datasetId = -1;
  other.spaceId = -1;
  other.typeId = -1;
}

DatasetHandle & DatasetHandle::operator=( DatasetHandle && other ) noexcept
{
  if( this != &other )
  {
    // Close existing resources
    if( typeId >= 0 )
    {
      H5E_BEGIN_TRY
      {
        H5Tclose( typeId );
      }
      H5E_END_TRY
    }
    if( spaceId >= 0 )
    {
      H5E_BEGIN_TRY
      {
        H5Sclose( spaceId );
      }
      H5E_END_TRY
    }
    if( datasetId >= 0 )
    {
      H5E_BEGIN_TRY
      {
        H5Dclose( datasetId );
      }
      H5E_END_TRY
    }

    // Transfer ownership
    datasetId = other.datasetId;
    spaceId = other.spaceId;
    typeId = other.typeId;
    dims = std::move( other.dims );

    other.datasetId = -1;
    other.spaceId = -1;
    other.typeId = -1;
  }
  return *this;
}

// Check if a dataset exists
bool DatasetHandle::datasetExists( hid_t const & fileId, string const & datasetName )
{
  htri_t status{-1};
  H5E_BEGIN_TRY
  {
    status = H5Oexists_by_name( fileId, datasetName.c_str(), H5P_DEFAULT );
  }
  H5E_END_TRY
  return status > 0 ? true : false;
}

// SerialHDF5Reader Implementation
SerialHDF5Reader::SerialHDF5Reader( const string & filename )
  : m_file( filename ) {}

// Specialization for array1d<uint8_t>
template<>
array1d< uint8_t > SerialHDF5Reader::read1DAs< uint8_t >( const std::string & datasetName, const int expectedDims ) const
{
  TypedArray1d result = read1D( m_file, datasetName, expectedDims );

  if( std::holds_alternative< array1d< uint8_t > >( result ) )
  {
    return std::get< array1d< uint8_t > >( result );
  }

  GEOS_THROW( GEOS_FMT( "Dataset {} in {} does not contain a compatible type for array1d<uint8_t> based on a promotion-only casting policy", datasetName, m_file.getFilename() ),
              InputError );
}

// Specialization for array1d<localIndex>
template<>
array1d< localIndex > SerialHDF5Reader::read1DAs< localIndex >( const std::string & datasetName, const int expectedDims ) const
{
  TypedArray1d result = read1D( m_file, datasetName, expectedDims );

  if( std::holds_alternative< array1d< localIndex > >( result ) )
  {
    return std::get< array1d< localIndex > >( result );
  }
  else if( std::holds_alternative< array1d< uint8_t > >( result ) )
  {
    return staticCastArray< uint8_t, localIndex >( std::get< array1d< uint8_t > >( result ) );
  }

  GEOS_THROW( GEOS_FMT( "Dataset {} in {} does not contain a compatible type for array1d<locaIndex> based on a promotion-only casting policy", datasetName, m_file.getFilename() ),
              InputError );
}

// Specialization for array1d<real32>
template<>
array1d< real32 > SerialHDF5Reader::read1DAs< real32 >( const std::string & datasetName, const int expectedDims ) const
{
  TypedArray1d result = read1D( m_file, datasetName, expectedDims );

  if( std::holds_alternative< array1d< real32 > >( result ) )
  {
    return std::get< array1d< real32 > >( result );
  }
  else if( std::holds_alternative< array1d< uint8_t > >( result ) )
  {
    return staticCastArray< uint8_t, real32 >( std::get< array1d< uint8_t > >( result ) );
  }
  else if( std::holds_alternative< array1d< localIndex > >( result ) )
  {
    return staticCastArray< localIndex, real32 >( std::get< array1d< localIndex > >( result ) );
  }

  GEOS_THROW( GEOS_FMT( "Dataset {} in {} does not contain a compatible type for array1d<real32> based on a promotion-only casting policy", datasetName, m_file.getFilename() ),
              InputError );
}

// Specialization for array1d<real64>
template<>
array1d< real64 > SerialHDF5Reader::read1DAs< real64 >( const std::string & datasetName, const int expectedDims ) const
{
  TypedArray1d result = read1D( m_file, datasetName, expectedDims );

  if( std::holds_alternative< array1d< real64 > >( result ) )
  {
    return std::get< array1d< real64 > >( result );
  }
  else if( std::holds_alternative< array1d< uint8_t > >( result ) )
  {
    return staticCastArray< uint8_t, real64 >( std::get< array1d< uint8_t > >( result ) );
  }
  else if( std::holds_alternative< array1d< localIndex > >( result ) )
  {
    return staticCastArray< localIndex, real64 >( std::get< array1d< localIndex > >( result ) );
  }
  else if( std::holds_alternative< array1d< real32 > >( result ) )
  {
    return staticCastArray< real32, real64 >( std::get< array1d< real32 > >( result ) );
  }

  GEOS_THROW( GEOS_FMT( "Dataset {} in {} does not contain a compatible type for array1d<real64> based on a promotion-only casting policy", datasetName, m_file.getFilename() ),
              InputError );
}


} // end of namespace hdf5Utils

} // end of namespace geos
