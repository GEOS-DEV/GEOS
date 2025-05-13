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

#include <variant>

// Forward declarations for HDF5 types
typedef int64_t hid_t;    // HDF5 identifier type
typedef unsigned long long hsize_t; // HDF5 size type

namespace geos
{

namespace hdf5Utils
{
using TypedArray1d = std::variant<
  array1d< globalIndex >,
  array1d< real32 >,
  array1d< real64 > >;

template< typename SOURCE_T, typename TARGET_T >
array1d< TARGET_T > staticCastArray( array1d< SOURCE_T > const & source )
{
  array1d< TARGET_T > casted( source.size() );
  std::transform( source.begin(), source.end(), casted.begin(),
                  []( auto value ) { return static_cast< TARGET_T >( value ); } );
  return casted;
}

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

  template< typename T >
  array1d< T > read1DAs( const std::string & datasetName, const int expectedDims ) const;

private:
  SerialHDF5File m_file;
};


// Specializations for supported types
template<>
array1d<globalIndex> SerialHDF5Reader::read1DAs< globalIndex >( const std::string & datasetName, const int expectedDims ) const;

template<>
array1d<real32> SerialHDF5Reader::read1DAs< real32 >( const std::string & datasetName, const int expectedDims ) const;

template<>
array1d<real64> SerialHDF5Reader::read1DAs< real64 >( const std::string & datasetName, const int expectedDims ) const;

} // end of namespace hdf5Utils


}  /* namespace geos */

 #endif /* GEOS_FUNCTIONS_HDF%UTILITIES_HPP_ */
