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
 * @brief Implementation details and helper structures for HDF5 file utilities in GEOS.
 */

#ifndef GEOS_FUNCTIONS_HDF5UTILITIES_IMPL_HPP_
#define GEOS_FUNCTIONS_HDF5UTILITIES_IMPL_HPP_

#include "common/DataTypes.hpp"
#include <H5Ipublic.h>
#include <H5public.h>

namespace geos
{

namespace hdf5Utils
{

/**
 * @struct DatasetHandle
 * @brief Helper struct to manage the lifetime and metadata of an HDF5 dataset.
 */
struct DatasetHandle
{
  /**
   * @brief Construct a handle for a dataset in the given file.
   * @param[in] hdf5File The HDF5 file object.
   * @param[in] datasetName The name of the dataset to open.
   * @param[in] expectedDims The expected dimensionality of the dataspace for the dataset.
   */
  DatasetHandle( SerialHDF5File const & hdf5File, string const & datasetName, int const expectedDims );

  /**
   * @brief Destructor. Closes all open HDF5 handles (datatype, dataspace, and dataset) if valid.
   */
  ~DatasetHandle();

  // Deleted copy constructor and assignment operator
  DatasetHandle( const DatasetHandle & ) = delete;
  DatasetHandle & operator=( const DatasetHandle & ) = delete;

  /**
   * @brief Move constructor.
   * @param[in] other The handle to move from.
   */
  DatasetHandle( DatasetHandle && other ) noexcept;

  /**
   * @brief Move assignment operator.
   * @param[in] other The handle to move from.
   * @return Reference to this object.
   */
  DatasetHandle & operator=( DatasetHandle && other ) noexcept;

  /**
   * @brief Check if a dataset exists in the file.
   * @param[in] fileId The HDF5 file identifier.
   * @param[in] datasetName The name of the dataset.
   * @return True if the dataset exists, false otherwise.
   */
  bool datasetExists( hid_t const & fileId, string const & datasetName );

  /// The name of the dataset.
  string m_datasetName;
  /// The HDF5 dataset identifier.
  hid_t datasetId = -1;
  /// The HDF5 dataspace identifier.
  hid_t spaceId = -1;
  /// The HDF5 datatype identifier.
  hid_t typeId = -1;
  /// The dimensions of the dataspace.
  array1d< hsize_t > dims;
};

} // end of namespace hdf5Utils

}  /* namespace geos */

#endif /* GEOS_FUNCTIONS_HDF5UTILITIES_IMPL_HPP_ */
