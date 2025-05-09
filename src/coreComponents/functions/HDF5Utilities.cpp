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
#include "common/format/Format.hpp"
#include "common/logger/Logger.hpp"

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
      GEOS_THROW( GEOS_FMT( "Error in call to:\n"                             \
                            "{}\n"                                            \
                            "({})",                                           \
                            #call, context ),                                 \
                  InputError );                                              \
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
      GEOS_THROW( GEOS_FMT( "Error in call to:\n"                                       \
                            "{}\n"                                                      \
                            "({})",                                                     \
                            #call, context ),                                           \
                  InputError );                                                        \
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
  std::string contextCheckMessage{ GEOS_FMT( "Dataset {} in {}", datasetName, hdf5File.getFilename() ) };

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
SerialHDF5Reader::SerialHDF5Reader( const std::string & filename )
  : m_file( filename ) {}

TypedArray1d SerialHDF5Reader::read1D( const std::string & datasetName, const int expectedDims ) const
{
  // Create a DatasetHandle to manage resources
  DatasetHandle handle( m_file, datasetName, expectedDims );

  // Compute the total number of elements
  hsize_t total_elements = 1;
  for( const auto & dim : handle.dims )
    total_elements *= dim;

  // Determine the type and read the data
  if( H5Tequal( handle.typeId, H5T_NATIVE_INT ))
  {
    array1d< globalIndex > buffer( total_elements );
    if( H5Dread( handle.datasetId, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer.data()) < 0 )
      throw std::runtime_error( "Failed to read dataset: " + datasetName );
    return buffer;
  }
  else if( H5Tequal( handle.typeId, H5T_NATIVE_FLOAT ))
  {
    array1d< real32 > buffer( total_elements );
    if( H5Dread( handle.datasetId, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer.data()) < 0 )
      throw std::runtime_error( "Failed to read dataset: " + datasetName );
    return buffer;
  }
  else if( H5Tequal( handle.typeId, H5T_NATIVE_DOUBLE ))
  {
    array1d< real64 > buffer( total_elements );
    if( H5Dread( handle.datasetId, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer.data()) < 0 )
      throw std::runtime_error( "Failed to read dataset: " + datasetName );
    return buffer;
  }
  else
  {
    throw std::runtime_error( "Unsupported dataset type for dataset: " + datasetName );
  }
}



} // end of namespace hdf5Utils

} // end of namespace geos
