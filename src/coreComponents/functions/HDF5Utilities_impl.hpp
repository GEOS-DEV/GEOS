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
 * @file HDF5Utilities_impl.hpp
 */

#ifndef GEOS_FUNCTIONS_HDF5UTILITIES_IMPL_HPP_
#define GEOS_FUNCTIONS_HDF5UTILITIES_IMPL_HPP_

#include "common/DataTypes.hpp"
#include "common/format/Format.hpp"
#include "common/logger/Logger.hpp"

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


template< typename T >
array1d< T > readTypedData( hid_t datasetId,
                            hid_t nativeType,
                            hsize_t total_elements,
                            string const & datasetName,
                            string const & filename )
{
  string contextCheckMessage{ GEOS_FMT( "Dataset {} in {}", datasetName, filename ) };
  array1d< T > buffer( total_elements );

  GEOS_HDF5_CHECK_ERROR( H5Dread( datasetId, nativeType, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer.data()),
                         contextCheckMessage );

  return buffer;
}

template< typename SOURCE_T, typename TARGET_T >
array1d< TARGET_T > staticCastArray( array1d< SOURCE_T > const & source )
{
  array1d< TARGET_T > casted( source.size() );
  std::transform( source.begin(), source.end(), casted.begin(),
                  []( auto value ) { return static_cast< TARGET_T >( value ); } );
  return casted;
}

// TypedArray1d read1D( string const & filename, string const & datasetName, int const expectedDims ) const;

} // end of namespace hdf5Utils

}  /* namespace geos */

 #endif /* GEOS_FUNCTIONS_HDF5UTILITIES_IMPL_HPP_ */
