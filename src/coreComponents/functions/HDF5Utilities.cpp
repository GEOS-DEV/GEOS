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
 
#define GEOS_HDF5_CHECK_ERROR(call, context)                         \
  do {                                                                            \
    herr_t __geos_lai_internal_ierr__ = -1;                                       \
    H5E_BEGIN_TRY {                                                               \
      __geos_lai_internal_ierr__ = (call);                                        \
    } H5E_END_TRY;                                                                \
    if (__geos_lai_internal_ierr__ < 0) {                                        \
      H5Eclear2(H5E_DEFAULT);                                                     \
      throw std::runtime_error(std::string("Error in call to ") +                \
                               #call + " (" + context + ")");                    \
    }                                                                             \
  } while (false)

  #define GEOS_HDF5_CHECK_AND_ASSIGN_ID(var, call, context)                       \
  do {                                                                            \
    hid_t __geos_lai_internal_result__ = H5I_INVALID_HID;                          \
    H5E_BEGIN_TRY {                                                               \
      __geos_lai_internal_result__ = (call);                                      \
    } H5E_END_TRY;                                                                \
    if (__geos_lai_internal_result__ < 0) {                                       \
      H5Eclear2(H5E_DEFAULT);                                                     \
      throw std::runtime_error(std::string("Error in call to ") +                \
                               #call + " (" + context + ")");                    \
    }                                                                             \
    var = __geos_lai_internal_result__;                                           \
  } while (false)

  #define GEOS_HDF5_CHECK_AND_ASSIGN_INT(var, call, context)                       \
  do {                                                                            \
    int __geos_lai_internal_result__ = -1;                          \
    H5E_BEGIN_TRY {                                                               \
      __geos_lai_internal_result__ = (call);                                      \
    } H5E_END_TRY;                                                                \
    if (__geos_lai_internal_result__ < 0) {                                       \
      H5Eclear2(H5E_DEFAULT);                                                     \
      throw std::runtime_error(std::string("Error in call to ") +                \
                               #call + " (" + context + ")");                    \
    }                                                                             \
    var = __geos_lai_internal_result__;                                           \
  } while (false)


 namespace geos {
 
namespace hdf5Utils
{
 
      //SerialHDF5File Implementation
      SerialHDF5File::SerialHDF5File(const string &filename) : m_filename(filename)
      {
        openFile();
      }
  
      // Destructor: Close the HDF5 file
      SerialHDF5File::~SerialHDF5File()
      {
        closeFile();
      }
  
      SerialHDF5File::SerialHDF5File(SerialHDF5File &&other) noexcept : m_fileId(other.m_fileId), m_filename(std::move(other.m_filename))
      {
        other.m_fileId = -1;
      }
  
      SerialHDF5File &SerialHDF5File::operator=(SerialHDF5File &&other) noexcept
      {
        if (this != &other)
        {
          closeFile();
          m_fileId = other.m_fileId;
          m_filename = std::move(other.m_filename);
          other.m_fileId = -1;
        }
        return *this;
      }
  
      hid_t const &SerialHDF5File::getFileId() const
      {
        return m_fileId;
    }
  
      // Access the filename
      string const &SerialHDF5File::getFilename() const
      {
        return
      m_filename; }
  
      void SerialHDF5File::openFile()
      {
        closeFile(); // Ensure any previously opened file is closed
        H5E_BEGIN_TRY
        {
          m_fileId = H5Fopen(m_filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
        }
        H5E_END_TRY
        GEOS_THROW_IF( m_fileId < 0,
                       GEOS_FMT( "Cannot open HDF5 file {}", getFilename() ),
                       InputError );
  
      }
  
      void SerialHDF5File::closeFile()
      {
        if (m_fileId >= 0)
        {
          herr_t status = -1;
          H5E_BEGIN_TRY
          {
            status = H5Fclose(m_fileId);
            // GEOS_HDF5_CHECK_ERROR( H5Fclose(m_fileId) ); //TODO add error checking as in linear algebra
          }
          H5E_END_TRY
          GEOS_THROW_IF( status < 0,
                         GEOS_FMT( "Cannot close HDF5 file {}", getFilename() ),
                         InputError );
          m_fileId = -1;
        }
      }
  
      // DatasetHandle Implementation
      
      DatasetHandle::DatasetHandle(SerialHDF5File const &hdf5File, string const &datasetName, int const expectedDims)
        : m_datasetName(datasetName)
      {
        GEOS_THROW_IF( !datasetExists( hdf5File.getFileId(), datasetName ),
                       GEOS_FMT( "Dataset {} doesn't exist in {}", datasetName, hdf5File.getFilename() ),
                       InputError );
  
        H5E_BEGIN_TRY
        {
          datasetId = H5Dopen2( hdf5File.getFileId(), datasetName.c_str(), H5P_DEFAULT);
        }
        H5E_END_TRY
        GEOS_THROW_IF( datasetId < 0,
                       GEOS_FMT( "Dataset {} cannot be opened in {}", datasetName, hdf5File.getFilename() ),
                       InputError );
  
        // Get the dataspace
        H5E_BEGIN_TRY
        {
          spaceId = H5Dget_space(datasetId);
        }
        H5E_END_TRY
        GEOS_THROW_IF( spaceId < 0,
                       GEOS_FMT( "Cannot get the dataspace for dataset {} in {}", datasetName, hdf5File.getFilename() ),
                       InputError );
  
        // Get the datatype
        H5E_BEGIN_TRY
        {
          typeId = H5Dget_type(datasetId);
        }
        H5E_END_TRY
        GEOS_THROW_IF( typeId < 0,
                       GEOS_FMT( "Cannot get the datatype for dataset {} in {}", datasetName, hdf5File.getFilename() ),
                       InputError );
  
        // Get the number of dimensions
        int ndims;
        H5E_BEGIN_TRY
        {
          ndims = H5Sget_simple_extent_ndims(spaceId);
        }
        H5E_END_TRY
        GEOS_THROW_IF( ndims < 0,
                       GEOS_FMT( "Cannot get the number of dimensions for dataset {} in {}", datasetName, hdf5File.getFilename() ),
                       InputError );
  
        // Validate dimensions if expectedDims is provided
        GEOS_THROW_IF( ndims != expectedDims,
                       GEOS_FMT( "Cannot get the number of dimensions for dataset {} in {}", datasetName, hdf5File.getFilename() ),
                       InputError );
  
        // Get the dimensions
        dims.resize(ndims);
        herr_t status;
        H5E_BEGIN_TRY
        {
          status = H5Sget_simple_extent_dims(spaceId, dims.data(), nullptr);
        }
        H5E_END_TRY
        GEOS_THROW_IF( status < 0,
                       GEOS_FMT( "Cannot get the dimensions for dataset {} in {}", datasetName, hdf5File.getFilename() ),
                       InputError );
      }
  
      // Destructor: Ensure all resources are closed
      DatasetHandle::~DatasetHandle()
      {
        if (typeId >= 0)
          H5Tclose(typeId);
        if (spaceId >= 0)
          H5Sclose(spaceId);
        if (datasetId >= 0)
          H5Dclose(datasetId);
      }
  
      // Allow move semantics
      DatasetHandle::DatasetHandle(DatasetHandle &&other) noexcept
          : datasetId(other.datasetId), spaceId(other.spaceId), typeId(other.typeId), dims(std::move(other.dims))
      {
        other.datasetId = -1;
        other.spaceId = -1;
        other.typeId = -1;
      }
  
      DatasetHandle &DatasetHandle::operator=(DatasetHandle &&other) noexcept
      {
        if (this != &other)
        {
          // Close existing resources
          if (typeId >= 0)
          {
            H5E_BEGIN_TRY
            {
              H5Tclose(typeId);
            }
            H5E_END_TRY
          }   
          if (spaceId >= 0)
          {
            H5E_BEGIN_TRY
            {
              H5Sclose(spaceId);
            }
            H5E_END_TRY
          }
          if (datasetId >= 0)
          {
            H5E_BEGIN_TRY
            {
              H5Dclose(datasetId);
            }
            H5E_END_TRY
          }
  
          // Transfer ownership
          datasetId = other.datasetId;
          spaceId = other.spaceId;
          typeId = other.typeId;
          dims = std::move(other.dims);
  
          other.datasetId = -1;
          other.spaceId = -1;
          other.typeId = -1;
        }
        return *this;
      }
  
      // Check if a dataset exists
      bool DatasetHandle::datasetExists(hid_t const &fileId, string const &datasetName)
      {
          herr_t err;
  
          H5E_BEGIN_TRY
          {
              err = H5Oexists_by_name(fileId, datasetName.c_str(), H5P_DEFAULT);
          }
          H5E_END_TRY
          return err > 0 ? true : false;
      }    
  
  
  
    // SerialHDF5Reader Implementation
      SerialHDF5Reader::SerialHDF5Reader(const std::string &filename)
          : m_file(filename) {}
  
      TypedArray1d SerialHDF5Reader::read1D(const std::string &datasetName, const int expectedDims) const
      {
        // Create a DatasetHandle to manage resources
        DatasetHandle handle(m_file, datasetName, expectedDims);
  
        // Compute the total number of elements
        hsize_t total_elements = 1;
        for (const auto &dim : handle.dims)
          total_elements *= dim;
  
        // Determine the type and read the data
        if (H5Tequal(handle.typeId, H5T_NATIVE_INT))
        {
          array1d<globalIndex> buffer(total_elements);
          if (H5Dread(handle.datasetId, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer.data()) < 0)
            throw std::runtime_error("Failed to read dataset: " + datasetName);
          return buffer;
        }
        else if (H5Tequal(handle.typeId, H5T_NATIVE_FLOAT))
        {
          array1d<real32> buffer(total_elements);
          if (H5Dread(handle.datasetId, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer.data()) < 0)
            throw std::runtime_error("Failed to read dataset: " + datasetName);
          return buffer;
        }
        else if (H5Tequal(handle.typeId, H5T_NATIVE_DOUBLE))
        {
          array1d<real64> buffer(total_elements);
          if (H5Dread(handle.datasetId, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer.data()) < 0)
            throw std::runtime_error("Failed to read dataset: " + datasetName);
          return buffer;
        }
        else
        {
          throw std::runtime_error("Unsupported dataset type for dataset: " + datasetName);
        }
      }
  

 
} // end of namespace hdf5Utils
 
} // end of namespace geos
 