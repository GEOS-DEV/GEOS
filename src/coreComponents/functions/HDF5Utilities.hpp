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
#include <H5Ipublic.h>

namespace geos
{

namespace hdf5Utils
{

/**
 * @class SerialHDF5File
 * @brief Manages the lifetime of a serial HDF5 file.
 */
class SerialHDF5File
{
public:
  /**
   * @brief Construct and open an HDF5 file.
   * @param[in] filename The name of the file to open.
   */
  explicit SerialHDF5File( const string & filename );

  /**
   * @brief Destructor. Closes the HDF5 file if open.
   */
  ~SerialHDF5File();

  // Deleted copy constructor and assignment operator
  SerialHDF5File( const SerialHDF5File & ) = delete;
  SerialHDF5File & operator=( const SerialHDF5File & ) = delete;

  /**
   * @brief Move constructor.
   * @param[in] other The file to move from.
   */
  SerialHDF5File( SerialHDF5File && other ) noexcept;

  /**
   * @brief Move assignment operator.
   * @param[in] other The file to move from.
   * @return Reference to this object.
   */
  SerialHDF5File & operator=( SerialHDF5File && other ) noexcept;

  /**
   * @brief Get the HDF5 file identifier.
   * @return The HDF5 file id.
   */
  hid_t const & getFileId() const { return m_fileId; }

  /**
   * @brief Get the filename associated with this file.
   * @return The filename.
   */
  string const & getFilename() const { return m_filename; }

private:
  /**
   * @brief Open the HDF5 file.
   */
  void openFile();

  /**
   * @brief Close the HDF5 file.
   */
  void closeFile();

  /// The HDF5 file identifier.
  hid_t m_fileId{-1};
  /// The filename.
  string m_filename;
};

/**
 * @class SerialHDF5Reader
 * @brief Provides read access to datasets in a serial HDF5 file.
 */
class SerialHDF5Reader
{
public:
  /**
   * @brief Construct a reader for the given HDF5 file.
   * @param[in] filename The name of the file to read.
   */
  explicit SerialHDF5Reader( string const & filename );

  /**
   * @brief Read an N-dimensional dataset and return it as a flat 1D array in Fortran (column-major) order.
   *
   * This function reads a dataset of arbitrary rank (specified by @p expectedDims),
   * checks that the rank matches the expectation, and returns the data flattened
   * into a one-dimensional array in Fortran array order.
   *
   * Only specific types are supported: uint8_t, localIndex, real32, and real64.
   * If an unsupported type is used, a compile-time error will occur.
   *
   * @tparam T The element type to return (must be one of: uint8_t, localIndex, real32, real64).
   * @param[in] datasetName The name of the dataset in the HDF5 file.
   * @param[in] expectedDims The expected number of dimensions of the dataset.
   * @return The contents of the dataset, flattened to a 1D array in Fortran order.
   */
  template< typename T >
  array1d< T > readAsFortranFlatArray( string const & datasetName,
                                       int const expectedDims ) const = delete;

private:
  /// The underlying HDF5 file.
  SerialHDF5File m_file;
};

//! @name Specializations for supported types
///@{

/**
 * @brief Specialization for reading a dataset as array1d<uint8_t>.
 * @param[in] datasetName The name of the dataset in the HDF5 file.
 * @param[in] expectedDims The expected number of dimensions of the dataset.
 * @return The contents of the dataset, flattened to a 1D array in Fortran order.
 */
template<>
array1d< uint8_t > SerialHDF5Reader::readAsFortranFlatArray< uint8_t >( const std::string & datasetName, const int expectedDims ) const;

/**
 * @brief Specialization for reading a dataset as array1d<localIndex>.
 * @param[in] datasetName The name of the dataset in the HDF5 file.
 * @param[in] expectedDims The expected number of dimensions of the dataset.
 * @return The contents of the dataset, flattened to a 1D array in Fortran order.
 */
template<>
array1d< localIndex > SerialHDF5Reader::readAsFortranFlatArray< localIndex >( const std::string & datasetName, const int expectedDims ) const;

/**
 * @brief Specialization for reading a dataset as array1d<real32>.
 * @param[in] datasetName The name of the dataset in the HDF5 file.
 * @param[in] expectedDims The expected number of dimensions of the dataset.
 * @return The contents of the dataset, flattened to a 1D array in Fortran order.
 */
template<>
array1d< real32 > SerialHDF5Reader::readAsFortranFlatArray< real32 >( const std::string & datasetName, const int expectedDims ) const;

/**
 * @brief Specialization for reading a dataset as array1d<real64>.
 * @param[in] datasetName The name of the dataset in the HDF5 file.
 * @param[in] expectedDims The expected number of dimensions of the dataset.
 * @return The contents of the dataset, flattened to a 1D array in Fortran order.
 */
template<>
array1d< real64 > SerialHDF5Reader::readAsFortranFlatArray< real64 >( const std::string & datasetName, const int expectedDims ) const;

///@}

} // end of namespace hdf5Utils

} // end of namespace geos

#endif /* GEOS_FUNCTIONS_HDF5UTILITIES_HPP_ */
