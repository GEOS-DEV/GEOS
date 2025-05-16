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


// Forward declarations for HDF5 types
typedef int64_t hid_t;    // HDF5 identifier type

namespace geos
{

namespace hdf5Utils
{
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

class SerialHDF5Reader
{
public:
  explicit SerialHDF5Reader( string const & filename );

  template< typename T >
  array1d< T > read1DAs( string const & datasetName, int const expectedDims ) const
  {
    static_assert( std::is_same_v< T, uint8_t > ||
                   std::is_same_v< T, localIndex > ||
                   std::is_same_v< T, real32 > ||
                   std::is_same_v< T, real64 >,
                   "Unsupported template type in read1DAs" );
    return {};
  }

private:
  SerialHDF5File m_file;
};


// Specializations for supported types
template<>
array1d< uint8_t > SerialHDF5Reader::read1DAs< uint8_t >( const std::string & datasetName, const int expectedDims ) const;

template<>
array1d< localIndex > SerialHDF5Reader::read1DAs< localIndex >( const std::string & datasetName, const int expectedDims ) const;

template<>
array1d< real32 > SerialHDF5Reader::read1DAs< real32 >( const std::string & datasetName, const int expectedDims ) const;

template<>
array1d< real64 > SerialHDF5Reader::read1DAs< real64 >( const std::string & datasetName, const int expectedDims ) const;

} // end of namespace hdf5Utils

}  /* namespace geos */

 #endif /* GEOS_FUNCTIONS_HDF5UTILITIES_HPP_ */
