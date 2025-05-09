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
 * @file HDF5Utilities.hpp
 */

 #ifndef GEOS_FUNCTIONS_HDF5UTILITIES_HPP_
#define GEOS_FUNCTIONS_HDF5UTILITIES_HPP_

#include "common/DataTypes.hpp"
#include <algorithm>
#include <hdf5.h>
#include <variant>

namespace geos
{

namespace hdf5Utils
{

static_assert( sizeof( H5T_NATIVE_INT ) == sizeof( globalIndex ),
               "H5T_NATIVE_INT and geos::integer must have the same size" );
static_assert( sizeof( H5T_NATIVE_FLOAT ) == sizeof( real64 ),
               "H5T_NATIVE_FLOAT and geos::real32 must have the same size" );
static_assert( sizeof( H5T_NATIVE_DOUBLE ) == sizeof( real64 ),
               "H5T_NATIVE_DOUBLE and geos::real64 must have the same size" );

using TypedArray1d = std::variant<
  array1d< globalIndex >,
  array1d< real32 >,
  array1d< real64 > >;


class SerialHDF5File
{
public:
  explicit SerialHDF5File( const string & filename );
  ~SerialHDF5File();

  SerialHDF5File( const SerialHDF5File & ) = delete;
  SerialHDF5File & operator=( const SerialHDF5File & ) = delete;

  // Allow move semantics
  SerialHDF5File( SerialHDF5File && other ) noexcept;
  SerialHDF5File & operator=( SerialHDF5File && other ) noexcept;

  hid_t const & getFileId() const
  { return m_fileId; }

  string const & getFilename() const
  { return m_filename; }

private:
  void openFile();
  void closeFile();

  hid_t m_fileId{-1};
  string m_filename;
};

struct DatasetHandle
{
  DatasetHandle( SerialHDF5File const & hdf5File, string const & datasetName, int const expectedDims );
  ~DatasetHandle();

  DatasetHandle( const DatasetHandle & ) = delete;
  DatasetHandle & operator=( const DatasetHandle & ) = delete;

  DatasetHandle( DatasetHandle && other ) noexcept;
  DatasetHandle & operator=( DatasetHandle && other ) noexcept;

  bool datasetExists( hid_t const & fileId, string const & datasetName );

  string m_datasetName;
  hid_t datasetId = -1;
  hid_t spaceId = -1;
  hid_t typeId = -1;
  array1d< hsize_t > dims;   // Store dataset dimensions
};

class SerialHDF5Reader
{
public:
  explicit SerialHDF5Reader( const std::string & filename );
  TypedArray1d read1D( const std::string & datasetName, const int expectedDims ) const;

private:
  SerialHDF5File m_file;
};


} // end of namespace hdf5Utils


}  /* namespace geos */

 #endif /* GEOS_FUNCTIONS_HDF%UTILITIES_HPP_ */
