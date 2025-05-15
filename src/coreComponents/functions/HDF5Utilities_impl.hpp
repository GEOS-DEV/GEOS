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

// Forward declarations for HDF5 types
typedef int64_t hid_t;    // HDF5 identifier type
typedef unsigned long long hsize_t; // HDF5 size type

namespace geos
{

namespace hdf5Utils
{

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

} // end of namespace hdf5Utils

}  /* namespace geos */

 #endif /* GEOS_FUNCTIONS_HDF5UTILITIES_IMPL_HPP_ */
