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
 * @file TableFunction.cpp
 */

#include "common/format/table/TableData.hpp"
#include "common/format/table/TableFormatter.hpp"
#include "common/format/table/TableLayout.hpp"
#include "TableFunction.hpp"
#include "LogLevelsInfo.hpp"
#include "codingUtilities/Parsing.hpp"
#include "common/DataTypes.hpp"
#include "common/MpiWrapper.hpp"

#include "HDF5Utilities.hpp"

#include <algorithm>

#include <hdf5.h>

namespace geos
{

// namespace hdf5Utils
// {

// static_assert( sizeof( H5T_NATIVE_INT ) == sizeof( globalIndex ),
//              "H5T_NATIVE_INT and geos::integer must have the same size" );
// static_assert( sizeof( H5T_NATIVE_FLOAT ) == sizeof( real64 ),
//              "H5T_NATIVE_FLOAT and geos::real32 must have the same size" );
// static_assert( sizeof( H5T_NATIVE_DOUBLE ) == sizeof( real64 ),
//              "H5T_NATIVE_DOUBLE and geos::real64 must have the same size" );

// using TypedArray1d = std::variant<
// array1d< globalIndex >,
// array1d< float >,
// array1d< double > >;


//   class SerialHDF5File
//   {
//   public:
//     // Constructor: Open the HDF5 file
//     explicit SerialHDF5File(const string &filename) : m_filename(filename)
//     {
//       openFile();
//     }

//     // Destructor: Close the HDF5 file
//     ~SerialHDF5File()
//     {
//       closeFile();
//     }

//     // Disable copy constructor and assignment operator
//     SerialHDF5File(const SerialHDF5File &) = delete;
//     SerialHDF5File &operator=(const SerialHDF5File &) = delete;

//     // Allow move semantics
//     SerialHDF5File(SerialHDF5File &&other) noexcept : m_fileId(other.m_fileId), m_filename(std::move(other.m_filename))
//     {
//       other.m_fileId = -1;
//     }

//     SerialHDF5File &operator=(SerialHDF5File &&other) noexcept
//     {
//       if (this != &other)
//       {
//         closeFile();
//         m_fileId = other.m_fileId;
//         m_filename = std::move(other.m_filename);
//         other.m_fileId = -1;
//       }
//       return *this;
//     }

//     // Access the file ID
//     hid_t const &getFileId() const { return m_fileId; }

//     // Access the filename
//     string const &getFilename() const { return m_filename; }

//   private:
//     void openFile()
//     {
//       closeFile(); // Ensure any previously opened file is closed
//       H5E_BEGIN_TRY
//       {
//         m_fileId = H5Fopen(m_filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
//       }
//       H5E_END_TRY
//       GEOS_THROW_IF( m_fileId < 0,
//                      GEOS_FMT( "Cannot open HDF5 file {}", getFilename() ),
//                      InputError );
//     }

//     void closeFile()
//     {
//       if (m_fileId >= 0)
//       {
//         herr_t status = -1;
//         H5E_BEGIN_TRY
//         {
//           status = H5Fclose(m_fileId);
//         }
//         H5E_END_TRY
//         GEOS_THROW_IF( status < 0,
//                        GEOS_FMT( "Cannot close HDF5 file {}", getFilename() ),
//                        InputError );
//         m_fileId = -1;
//       }
//     }

//     hid_t m_fileId{-1}; // HDF5 file handle
//     string m_filename;  // HDF5 file name
//   };


//   struct DatasetHandle
//   {
//     string m_datasetName;
//     hid_t datasetId = -1;
//     hid_t spaceId = -1;
//     hid_t typeId = -1;
//     array1d< hsize_t > dims; // Store dataset dimensions

//     // Constructor: Open dataset, dataspace, and datatype, and validate dimensions
//     DatasetHandle(SerialHDF5File const &hdf5File, string const &datasetName, int const expectedDims)
//     {
//       m_datasetName = datasetName;

//       // Check dataset existance
//       GEOS_THROW_IF( !datasetExists( hdf5File.getFileId(), datasetName ),
//                      GEOS_FMT( "Dataset {} doesn't exist in {}", datasetName, hdf5File.getFilename() ),
//                      InputError );

//       // Open the dataset
//       H5E_BEGIN_TRY
//       {
//         datasetId = H5Dopen2( hdf5File.getFileId(), datasetName.c_str(), H5P_DEFAULT);
//       }
//       H5E_END_TRY
//       GEOS_THROW_IF( datasetId < 0,
//                      GEOS_FMT( "Dataset {} cannot be opened in {}", datasetName, hdf5File.getFilename() ),
//                      InputError );

//       // Get the dataspace
//       H5E_BEGIN_TRY
//       {
//         spaceId = H5Dget_space(datasetId);
//       }
//       H5E_END_TRY
//       GEOS_THROW_IF( spaceId < 0,
//                      GEOS_FMT( "Cannot get the dataspace for dataset {} in {}", datasetName, hdf5File.getFilename() ),
//                      InputError );

//       // Get the datatype
//       H5E_BEGIN_TRY
//       {
//         typeId = H5Dget_type(datasetId);
//       }
//       H5E_END_TRY
//       GEOS_THROW_IF( typeId < 0,
//                      GEOS_FMT( "Cannot get the datatype for dataset {} in {}", datasetName, hdf5File.getFilename() ),
//                      InputError );

//       // Get the number of dimensions
//       int ndims;
//       H5E_BEGIN_TRY
//       {
//         ndims = H5Sget_simple_extent_ndims(spaceId);
//       }
//       H5E_END_TRY
//       GEOS_THROW_IF( ndims < 0,
//                      GEOS_FMT( "Cannot get the number of dimensions for dataset {} in {}", datasetName, hdf5File.getFilename() ),
//                      InputError );

//       // Validate dimensions if expectedDims is provided
//       GEOS_THROW_IF( ndims != expectedDims,
//                      GEOS_FMT( "Cannot get the number of dimensions for dataset {} in {}", datasetName, hdf5File.getFilename() ),
//                      InputError );

//       // Get the dimensions
//       dims.resize(ndims);
//       herr_t status;
//       H5E_BEGIN_TRY
//       {
//         status = H5Sget_simple_extent_dims(spaceId, dims.data(), nullptr);
//       }
//       H5E_END_TRY
//       GEOS_THROW_IF( status < 0,
//                      GEOS_FMT( "Cannot get the dimensions for dataset {} in {}", datasetName, hdf5File.getFilename() ),
//                      InputError );
//     }

//     // Destructor: Ensure all resources are closed
//     ~DatasetHandle()
//     {
//       if (typeId >= 0)
//         H5Tclose(typeId);
//       if (spaceId >= 0)
//         H5Sclose(spaceId);
//       if (datasetId >= 0)
//         H5Dclose(datasetId);
//     }

//     // Disable copy semantics
//     DatasetHandle(const DatasetHandle &) = delete;
//     DatasetHandle &operator=(const DatasetHandle &) = delete;

//     // Allow move semantics
//     DatasetHandle(DatasetHandle &&other) noexcept
//         : datasetId(other.datasetId), spaceId(other.spaceId), typeId(other.typeId), dims(std::move(other.dims))
//     {
//       other.datasetId = -1;
//       other.spaceId = -1;
//       other.typeId = -1;
//     }

//     DatasetHandle &operator=(DatasetHandle &&other) noexcept
//     {
//       if (this != &other)
//       {
//         // Close existing resources
//         if (typeId >= 0)
//         {
//           H5E_BEGIN_TRY
//           {
//             H5Tclose(typeId);
//           }
//           H5E_END_TRY
//         }
//         if (spaceId >= 0)
//         {
//           H5E_BEGIN_TRY
//           {
//             H5Sclose(spaceId);
//           }
//           H5E_END_TRY
//         }
//         if (datasetId >= 0)
//         {
//           H5E_BEGIN_TRY
//           {
//             H5Dclose(datasetId);
//           }
//           H5E_END_TRY
//         }

//         // Transfer ownership
//         datasetId = other.datasetId;
//         spaceId = other.spaceId;
//         typeId = other.typeId;
//         dims = std::move(other.dims);

//         other.datasetId = -1;
//         other.spaceId = -1;
//         other.typeId = -1;
//       }
//       return *this;
//     }

//     // Check if a dataset exists
//     static bool datasetExists(hid_t const &fileId, string const &datasetName)
//     {
//         herr_t err;

//         H5E_BEGIN_TRY
//         {
//             err = H5Oexists_by_name(fileId, datasetName.c_str(), H5P_DEFAULT);
//         }
//         H5E_END_TRY
//         return err > 0 ? true : false;
//     }
//   };



//   class SerialHDF5Reader
//   {
//   public:
//     explicit SerialHDF5Reader(const std::string &filename)
//         : m_file(filename) {}

//     TypedArray1d read1D(const std::string &datasetName, const int expectedDims) const
//     {
//       // Create a DatasetHandle to manage resources
//       DatasetHandle handle(m_file, datasetName, expectedDims);

//       // Compute the total number of elements
//       hsize_t total_elements = 1;
//       for (const auto &dim : handle.dims)
//         total_elements *= dim;

//       // Determine the type and read the data
//       if (H5Tequal(handle.typeId, H5T_NATIVE_INT))
//       {
//         array1d<globalIndex> buffer(total_elements);
//         if (H5Dread(handle.datasetId, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer.data()) < 0)
//           throw std::runtime_error("Failed to read dataset: " + datasetName);
//         return buffer;
//       }
//       else if (H5Tequal(handle.typeId, H5T_NATIVE_FLOAT))
//       {
//         array1d<real32> buffer(total_elements);
//         if (H5Dread(handle.datasetId, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer.data()) < 0)
//           throw std::runtime_error("Failed to read dataset: " + datasetName);
//         return buffer;
//       }
//       else if (H5Tequal(handle.typeId, H5T_NATIVE_DOUBLE))
//       {
//         array1d<real64> buffer(total_elements);
//         if (H5Dread(handle.datasetId, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer.data()) < 0)
//           throw std::runtime_error("Failed to read dataset: " + datasetName);
//         return buffer;
//       }
//       else
//       {
//         throw std::runtime_error("Unsupported dataset type for dataset: " + datasetName);
//       }
//     }

//   private:
//     SerialHDF5File m_file;
//   };

// struct serialHDF5File
// {
//   hid_t m_fileId;
//   string m_filename;

//   serialHDF5File( string const & name ):
//     m_fileId( -1 ),
//     m_filename( name )
//   {}
// };

// void openHDF5File( serialHDF5File & hdfFile )
// {
//   H5E_BEGIN_TRY
//   {
//     hdfFile.m_fileId = H5Fopen( hdfFile.m_filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );
//   }
//   H5E_END_TRY

//   GEOS_THROW_IF( hdfFile.m_fileId < 0,
//                  GEOS_FMT( "hdf5 file {} cannot be opened", hdfFile.m_filename ),
//                  InputError );
// }

// void closeHDF5File( serialHDF5File & hdfFile )
// {
//   herr_t err;

//   H5E_BEGIN_TRY
//   {
//     err = H5Fclose( hdfFile.m_fileId );
//   }
//   H5E_END_TRY

//   GEOS_THROW_IF( err < 0,
//                  GEOS_FMT( "hdf5 file {} cannot be closed", hdfFile.m_filename ),
//                  InputError );
// }

// bool datasetExists( hid_t const & fileId, string const & datasetName )
// {
//   herr_t err;

//   H5E_BEGIN_TRY
//   {
//     err = H5Oexists_by_name( fileId, datasetName.c_str(), H5P_DEFAULT );
//   }
//   H5E_END_TRY
//   return err > 0 ? true : false;
// }

// array1d< real64 > readDataset( serialHDF5File & hdfFile, string const & datasetName, int const expectedNumDims )
// {
//   // Check dataset existance
//   GEOS_THROW_IF( !datasetExists( hdfFile.m_fileId, datasetName ),
//                  GEOS_FMT( "Dataset {} doesn't exist in {}", datasetName, hdfFile.m_filename ),
//                  InputError );

//   // Open the dataset
//   hid_t dataset_id = H5Dopen( hdfFile.m_fileId, datasetName.c_str(), H5P_DEFAULT );
//   GEOS_ERROR_IF( dataset_id < 0, GEOS_FMT( "Error opening dataset {} in {}", datasetName, hdfFile.m_filename ) );

//   // Get the dataspace of the dataset
//   hid_t dataspace_id = H5Dget_space( dataset_id );

//   // Get the dataset dimensions
//   int ndims = H5Sget_simple_extent_ndims( dataspace_id );
//   if( ndims != expectedNumDims )
//   {
//     H5Sclose( dataspace_id );
//     H5Dclose( dataset_id );
//     H5Fclose( hdfFile.m_fileId );
//     GEOS_ERROR( GEOS_FMT( "Dataset {} is not {}D", datasetName, std::to_string( expectedNumDims ) ) );
//   }

//   array1d< hsize_t >  dims( ndims );
//   H5Sget_simple_extent_dims( dataspace_id, dims.data(), nullptr );

//   // Compute the total number of entries
//   hsize_t numDatasetEntries = 1;
//   for( auto const & dim : dims )
//   {
//     numDatasetEntries *= dim;
//   }

//   // Get the dataset type
//   hid_t dataset_type = H5Dget_type( dataset_id );

//   // Read the dataset
//   array1d< real64 > datasetValues( numDatasetEntries );
//   herr_t err{};

//   if( (H5Tequal( dataset_type, H5T_NATIVE_FLOAT ) ) )
//   {
//     array1d< real32 > tmp( numDatasetEntries );
//     H5E_BEGIN_TRY
//     {
//       err = H5Dread( dataset_id, H5T_NATIVE_FLOAT, dataspace_id, H5S_ALL, H5P_DEFAULT, tmp.data() );
//     }
//     H5E_END_TRY
//     GEOS_ERROR_IF( err < 0, GEOS_FMT( "{}: error reading single-precision dataset", datasetName ) );
//     for( int i = 0; i < tmp.size(); ++i )
//     {
//       datasetValues[i] = static_cast< double >( tmp[i] );
//     }
//   }
//   else if( (H5Tequal( dataset_type, H5T_NATIVE_DOUBLE ) ) )
//   {
//     H5E_BEGIN_TRY
//     {
//       err = H5Dread( dataset_id, H5T_NATIVE_DOUBLE, dataspace_id, H5S_ALL, H5P_DEFAULT, datasetValues.data() );
//     }
//     H5E_END_TRY
//     GEOS_ERROR_IF( err < 0, GEOS_FMT( "{}: error reading double-precision dataset", datasetName ) );
//   }
//   else
//   {
//     GEOS_ERROR( GEOS_FMT( "Dataset {} is not of type float", datasetName ) );
//   }

//   GEOS_ERROR_IF( err < 0, GEOS_FMT( "{}: error reading dataset", datasetName ) );

//   H5Tclose( dataset_type );
//   H5Sclose( dataspace_id );
//   H5Dclose( dataset_id );

//   return datasetValues;
// }

// } // end of namespace hdf5Utils

using namespace dataRepository;

TableFunction::TableFunction( const string & name,
                              Group * const parent ):
  FunctionBase( name, parent ),
  m_interpolationMethod( InterpolationType::Linear ),
  m_valueUnit( units::Unknown ),
  m_kernelWrapper( createKernelWrapper() )
{
  registerWrapper( viewKeyStruct::coordinatesString(), &m_tableCoordinates1D ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Coordinates inputs for 1D tables" );

  registerWrapper( viewKeyStruct::valuesString(), &m_values ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Values for 1D tables" );

  registerWrapper( viewKeyStruct::coordinateFilesString(), &m_coordinateFiles ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "List of coordinate file names for ND Table" );

  registerWrapper( viewKeyStruct::voxelFileString(), &m_voxelFile ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Voxel file name for ND Table" );

  registerWrapper( viewKeyStruct::hdf5FileString(), &m_hdf5FileName ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "HDF5 file name for ND Table" );

  registerWrapper( viewKeyStruct::hdf5CoordinateDatasetNamesString(), &m_hdf5CoordinateDatasetNames ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "List of coordinate dataset names in HDF5 file" );

  registerWrapper( viewKeyStruct::hdf5TableDatasetNameString(), &m_hdf5TableDatasetName ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Table dataset name in HDF5 file" );

  registerWrapper( viewKeyStruct::interpolationString(), &m_interpolationMethod ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Interpolation method. Valid options:\n* " + EnumStrings< InterpolationType >::concat( "\n* " ) ).
    setApplyDefaultValue( m_interpolationMethod );

  registerWrapper( viewKeyStruct::writeCSVFlagString(), &m_writeCSV ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "When set to 1, write the table into a CSV file" );

  addLogLevel< logInfo::TableDataOutput >();
}

TableFunction::TableInputType TableFunction::determineTableInputType() const
{
  // Determine which input types are provided
  bool has1DTableInputs = !m_tableCoordinates1D.empty() && !m_values.empty();
  bool hasNDTableInputs = !m_coordinateFiles.empty() && !m_voxelFile.empty();
  bool hasHDF5TableInputs = !m_hdf5FileName.empty() && !m_hdf5CoordinateDatasetNames.empty() && !m_hdf5TableDatasetName.empty();

  // Ensure mutual exclusivity of the three options
  int const inputTypeCount = static_cast< int >(has1DTableInputs) + static_cast< int >(hasNDTableInputs) + static_cast< int >(hasHDF5TableInputs);
  GEOS_THROW_IF( inputTypeCount > 1,
                 GEOS_FMT( "{} {}: Multiple table input types are provided. Only one of the following is allowed:\n"
                           "1. 1D table inputs (coordinates and values)\n"
                           "2. ND table inputs (coordinate files and voxel file)\n"
                           "3. HDF5 table inputs (HDF5 file, coordinate dataset names, and table dataset name).",
                           catalogName(), getDataContext()),
                 InputError );

  // Return the determined input type
  if( has1DTableInputs )
  {
    return TableInputType::OneD;
  }
  if( hasNDTableInputs )
  {
    return TableInputType::ND;
  }
  if( hasHDF5TableInputs )
  {
    return TableInputType::HDF5;
  }

  return TableInputType::None;
}

void TableFunction::readFile( string const & filename, array1d< real64 > & target )
{
  auto const skipped = []( char const c ){ return std::isspace( c ) || c == ','; };
  try
  {
    parseFile( filename, target, skipped );
  }
  catch( std::runtime_error const & e )
  {
    GEOS_THROW( GEOS_FMT( "{} {}: {}", catalogName(), getDataContext(), e.what() ), InputError );
  }
}

void TableFunction::setInterpolationMethod( InterpolationType const method )
{
  m_interpolationMethod = method;
  reInitializeFunction();
}

void TableFunction::setTableCoordinates( array1d< real64_array > const & coordinates,
                                         stdVector< units::Unit > const & dimUnits )
{
  m_dimUnits = dimUnits;
  m_coordinates.resize( 0 );
  for( localIndex i = 0; i < coordinates.size(); ++i )
  {
    for( localIndex j = 1; j < coordinates[i].size(); ++j )
    {
      GEOS_THROW_IF( coordinates[i][j] - coordinates[i][j-1] <= 0,
                     GEOS_FMT( "{} {}: coordinates must be strictly increasing, but axis {} is not",
                               catalogName(), getDataContext(), i ),
                     InputError );
    }
    m_coordinates.appendArray( coordinates[i].begin(), coordinates[i].end() );
  }
  reInitializeFunction();
}

void TableFunction::setTableValues( real64_array values, units::Unit unit )
{
  m_values = std::move( values );
  m_valueUnit = unit;
  reInitializeFunction();
}

void TableFunction::initializeFunction()
{
  if( m_coordinates.size() > 0 )
  {
    // This function appears to be already initialized
    // Apparently, this can be called multiple times during unit tests?
  }
  else
  {
    TableFunction::TableInputType inputType = determineTableInputType();

    // Read in data
    switch( inputType )
    {
      case TableInputType::OneD:
      {
        m_coordinates.appendArray( m_tableCoordinates1D.begin(), m_tableCoordinates1D.end());
        GEOS_THROW_IF_NE_MSG( m_tableCoordinates1D.size(), m_values.size(),
                              GEOS_FMT( "{} {}: 1D table function coordinates and values must have the same length",
                                        catalogName(), getDataContext()),
                              InputError );
        break;
      }

      case TableInputType::ND:
      {
        array1d< real64 > tmp;
        localIndex numValues = 1;
        for( localIndex ii = 0; ii < m_coordinateFiles.size(); ++ii )
        {
          tmp.clear();
          readFile( m_coordinateFiles[ii], tmp );
          m_coordinates.appendArray( tmp.begin(), tmp.end());
          numValues *= tmp.size();
        }
        m_values.reserve( numValues );
        readFile( m_voxelFile, m_values );
        break;
      }

      case TableInputType::HDF5:
      {
        hdf5Utils::SerialHDF5Reader reader( m_hdf5FileName );

        // Read table coordinates
        int tableDim = 0;
        for( auto const & coordName : m_hdf5CoordinateDatasetNames )
        {
          hdf5Utils::TypedArray1d tmp = reader.read1D( coordName, 1 );

          if( std::holds_alternative< array1d< real64 > >( tmp ))
          {
            m_coordinates.appendArray( std::get< array1d< real64 > >( tmp ).begin(), std::get< array1d< real64 > >( tmp ).end());
          }
          else if( std::holds_alternative< array1d< float > >( tmp ))
          {
            const auto & floatArray = std::get< array1d< globalIndex > >( tmp );
            m_coordinates.reserve( floatArray.size() );
            std::transform( floatArray.begin(), floatArray.end(), m_coordinates[tableDim].begin(), []( real32 val )
            { return static_cast< real64 >(val); } );
          }
          else if( std::holds_alternative< array1d< globalIndex > >( tmp ))
          {
            const auto & intArray = std::get< array1d< globalIndex > >( tmp );
            m_coordinates.reserve( intArray.size());
            std::transform( intArray.begin(), intArray.end(), m_coordinates[tableDim].begin(), []( globalIndex val )
            { return static_cast< real64 >(val); } );
          }
          else
          {
            // If tmp holds an unsupported type, throw an error
            throw std::runtime_error( "Unexpected data type for voxel dataset: " + coordName );
          }
          tableDim += 1;
        }

        // Read table dataset
        hdf5Utils::TypedArray1d tmp = reader.read1D( m_hdf5TableDatasetName, tableDim );
        if( std::holds_alternative< array1d< real64 > >( tmp ))
        {
          m_values = std::move( std::get< array1d< real64 > >( tmp ) );
        }
        // else if (std::holds_alternative<std::vector<float>>(tmp))
        // {
        //   // If tmp holds std::vector<real32>, cast to std::vector<double> and push
        //   const auto &floatVec = std::get<array1d<real32>>(tmp);
        //   m_values.resize(floatVec.size());
        //   std::transform(floatVec.begin(), floatVec.end(), m_values.begin(), [](float val)
        //                  { return static_cast<double>(val); });
        // }
        else
        {
          // If tmp holds an unsupported type, throw an error
          throw std::runtime_error( "Unexpected data type for voxel dataset: " + m_hdf5TableDatasetName );
        }

        break;
      }

      case TableInputType::None:
      default:
      {
        GEOS_THROW( GEOS_FMT( "{} {}: No valid table input type is provided.",
                              catalogName(), getDataContext()),
                    InputError );
        break;
      }
    }
  }

  reInitializeFunction();
}

void TableFunction::reInitializeFunction()
{
  // Setup index increment (assume data is in Fortran array order)
  localIndex increment = 1;
  for( localIndex ii = 0; ii < m_coordinates.size(); ++ii )
  {
    increment *= m_coordinates.sizeOfArray( ii );
    for( localIndex j = 1; j < m_coordinates[ii].size(); ++j )
    {
      GEOS_THROW_IF( m_coordinates[ii][j] - m_coordinates[ii][j-1] <= 0,
                     GEOS_FMT( "{} {}: coordinates must be strictly increasing, but axis {} is not",
                               catalogName(), getDataContext(), ii ),
                     InputError );
    }
  }
  if( m_coordinates.size() > 0 && !m_values.empty() ) // coordinates and values have been set
  {
    GEOS_THROW_IF_NE_MSG( increment, m_values.size(),
                          GEOS_FMT( "{} {}: number of values does not match total number of table coordinates",
                                    catalogName(), getDataContext() ),
                          InputError );
  }

  // Create the kernel wrapper
  m_kernelWrapper = createKernelWrapper();
}

void TableFunction::checkCoord( real64 const coord, localIndex const dim ) const
{
  GEOS_THROW_IF( dim >= m_coordinates.size() || dim < 0,
                 GEOS_FMT( "{}: The {} dimension ( no. {} ) doesn't exist in the table.",
                           getDataContext(), units::getDescription( getDimUnit( dim ) ), dim ),
                 SimulationError );
  real64 const lowerBound = m_coordinates[dim][0];
  real64 const upperBound = m_coordinates[dim][m_coordinates.sizeOfArray( dim ) - 1];
  GEOS_THROW_IF( coord > upperBound || coord < lowerBound,
                 GEOS_FMT( "{}: Requested {} is out of the table bounds ( lower bound: {} -> upper bound: {} ).",
                           getDataContext(),
                           units::formatValue( coord, getDimUnit( dim ) ),
                           units::formatValue( lowerBound, getDimUnit( dim ) ),
                           units::formatValue( upperBound, getDimUnit( dim ) ) ),
                 SimulationError );
}

TableFunction::KernelWrapper TableFunction::createKernelWrapper() const
{
  return { m_interpolationMethod,
           m_coordinates.toViewConst(),
           m_values.toViewConst() };
}

real64 TableFunction::evaluate( real64 const * const input ) const
{
  return m_kernelWrapper.compute( input );
}

real64 TableFunction::getCoord( real64 const * const input, localIndex const dim, InterpolationType interpolationMethod ) const
{
  GEOS_ASSERT_MSG( interpolationMethod != InterpolationType::Linear,
                   GEOS_FMT( "{}: TableFunction::getCoord should not be called with {} interpolation method",
                             getDataContext(), EnumStrings< InterpolationType >::toString( interpolationMethod )));
  GEOS_ASSERT( dim >= 0 && dim < m_coordinates.size() );
  return m_kernelWrapper.getCoord( input, dim, interpolationMethod );
}

TableFunction::KernelWrapper::KernelWrapper( InterpolationType const interpolationMethod,
                                             ArrayOfArraysView< real64 const > const & coordinates,
                                             arrayView1d< real64 const > const & values )
  :
  m_interpolationMethod( interpolationMethod ),
  m_coordinates( coordinates ),
  m_values( values )
{}


string TableFunction::getTableDescription() const
{
  std::ostringstream description;
  auto streamArrayDescription = [&]( string_view name,
                                     auto const & array,
                                     units::Unit const unit,
                                     string_view path )
  {
    description << GEOS_FMT( "- {}", string( name ) );
    if( unit != units::Unknown )
      description << GEOS_FMT( " in {} units", units::getDescription( unit ) );
    if( !path.empty() )
      description << GEOS_FMT( " from file: {}", path );
    description << '\n';

    auto const [minValue, maxValue] = std::minmax_element( array.begin(), array.end());
    description << GEOS_FMT( "  * {} values, from {} [{}] to {} [{}].",
                             array.size(),
                             *minValue, units::getSymbol( unit ),
                             *maxValue, units::getSymbol( unit ) );
  };

  for( integer dimId = 0; dimId < numDimensions(); ++dimId )
  {
    string const coordFilePath = dimId < m_coordinateFiles.size() ? m_coordinateFiles[dimId].relativeFilePath() : "";
    streamArrayDescription( "Coordinates " + getCoordsDescription( dimId, true ),
                            m_coordinates[dimId].toSliceConst(),
                            getDimUnit( dimId ),
                            coordFilePath );
    description << '\n';
  }
  streamArrayDescription( "Values",
                          m_values.toViewConst(),
                          getValueUnit(),
                          m_voxelFile.relativeFilePath() );

  return description.str();
}

string TableFunction::getCoordsDescription( integer dimId, bool shortUnitsToVariables ) const
{
  integer const numDims = numDimensions();
  units::Unit const dimCoordsUnit = getDimUnit( dimId );
  if( dimCoordsUnit != units::Unknown )
  {
    return string( shortUnitsToVariables ?
                   units::getVariableSymbol( dimCoordsUnit ) :
                   units::getDescription( dimCoordsUnit ) );
  }
  else
  {
    static constexpr string_view table2DGenericAxes[] = {"x", "y"};
    return numDims<=2 ? string( table2DGenericAxes[dimId] ) : GEOS_FMT( "Coord_{}", dimId );
  }
}

string TableFunction::getValuesDescription() const
{
  return string( getValueUnit() != units::Unknown ?
                 units::getDescription( getValueUnit() ) :
                 "Value" );
}

/**
 * @brief Retrieve all data values
 * @param table The table which contains the formatted data:
 *              Each row contains the coordinates followed by the value.
 * @param numDimensions Numbers of axes in the table
 * @param coordinates The tables axis values
 * @param values The table values to be retrived
 */
void collectTableValues( TableData & table,
                         integer const numDimensions,
                         ArrayOfArraysView< real64 const > const coordinates,
                         arrayView1d< real64 const > const values )
{
  // prepare dividers
  stdVector< integer > div( numDimensions );
  div[0] = 1;
  for( integer d = 1; d < numDimensions; d++ )
  {
    div[d] = div[d-1] * coordinates[d-1].size();
  }
  // loop through all the values
  stdVector< integer > coordsIdx( numDimensions );
  stdVector< TableData::CellData > rowData( numDimensions + 1,
                                            { CellType::Value, string() } );
  for( integer v = 0; v < values.size(); v++ )
  {
    // find coords indices
    integer r = v;
    for( integer d = numDimensions-1; d >= 0; d-- )
    {
      coordsIdx[d] = r / div[d];
      r = r % div[d];
    }
    // finally print out in right order
    for( integer d = 0; d < numDimensions; d++ )
    {
      rowData[d].value = GEOS_FMT( "{}", coordinates[d][coordsIdx[d]] );
    }
    rowData.back().value = GEOS_FMT( "{}", values[v] );
    table.addRow( rowData );
  }
}

bool isTableTooLargeForLog( TableFunction const & table )
{
  static constexpr integer maxWidth = 250;
  static constexpr integer maxRows = 500;
  // for now, we only estimate the table width from approximations
  static constexpr integer meanColumnWidth = 16;
  integer const numDims = table.numDimensions();
  integer const columnCount = numDims != 2 ? numDims + 1 : table.getCoordinates()[0].size();
  integer const columnSepWidth = numDims != 2 ? 5 : 3;
  integer tableApproxWidth = columnCount * (columnSepWidth + meanColumnWidth);
  integer tableRowsCount = numDims != 2 ? table.getValues().size() : table.getCoordinates()[1].size();

  return tableApproxWidth > maxWidth || tableRowsCount > maxRows;
}

void TableFunction::outputTableData( OutputOptions const outputOpts ) const
{
  // we only output from rank 0
  if( MpiWrapper::commRank() != 0 )
    return;

  bool const logOutputFailed = outputOpts.writeInLog && isTableTooLargeForLog( *this );
  if( outputOpts.writeInLog && !logOutputFailed )
  { // log output
    TableLayout tableLayout( getName(), {} );
    if( outputOpts.writeCSV )
    {
      tableLayout.addColumn( GEOS_FMT( "- CSV Generated to:\n  {}/{}.csv", getOutputDirectory(), getName() ) );
    }
    TableTextFormatter textFormatter( tableLayout );
    GEOS_LOG( textFormatter.toString( *this ));
  }

  if( outputOpts.writeCSV || logOutputFailed )
  {
    { // csv output
      std::ofstream logStream( joinPath( FunctionBase::getOutputDirectory(), getName() + ".csv" ) );
      TableCSVFormatter csvFormatter;
      logStream << csvFormatter.toString( *this );
    }

    if( !outputOpts.writeInLog )
    { // mini-table in log to notice user where csv has been output (if only csv output is enabled)
      // only one column which serve as "title" (centered), next, stats & texts are designed to be left-aligned
      TableLayout const tableLayout( { TableLayout::Column().
                                         setName( getName() ).
                                         setHeaderAlignment( TableLayout::Alignment::center ).
                                         setValuesAlignment( TableLayout::Alignment::left ) } );
      TableTextFormatter const tableLog( tableLayout );

      TableData tableData;
      tableData.addRow( getTableDescription());
      if( logOutputFailed )
      {
        tableData.addSeparator();
        tableData.addRow( " / \\ The table was too heavy for log output.\n"
                          "/ ! \\ To visualize the table, please refer to the generated csv." );
      }
      tableData.addSeparator();
      tableData.addRow( GEOS_FMT( "CSV Generated to:\n{}/{}.csv", getOutputDirectory(), getName() ) );
      GEOS_LOG( tableLog.toString( tableData ) );
    }
  }
}

void TableFunction::initializePostSubGroups()
{
  // Output user defined tables (not generated PVT tables)
  outputTableData( OutputOptions{
      m_writeCSV != 0,   // writeCSV
      isLogLevelActive< logInfo::TableDataOutput >( getLogLevel() )   // writeInLog
    } );
}

template<>
string TableCSVFormatter::toString< TableFunction >( TableFunction const & tableFunction ) const
{
  ArrayOfArraysView< real64 const > const coordinates = tableFunction.getCoordinates();
  arrayView1d< real64 const > const values = tableFunction.getValues();
  TableLayout tableLayout;

  integer const numDimensions = LvArray::integerConversion< integer >( coordinates.size() );
  if( numDimensions != 2 )
  {
    for( integer d = 0; d < numDimensions; d++ )
    {
      tableLayout.addColumn( units::getDescription( tableFunction.getDimUnit( d ) ) );
    }
    tableLayout.addColumn( units::getDescription( tableFunction.getValueUnit() ) );

    TableData tableData;
    collectTableValues( tableData, numDimensions, coordinates, values );

    TableCSVFormatter const csvFormat( tableLayout );
    return csvFormat.toString( tableData );
  }
  else
  {
    TableData2D tableData2D;
    TableData2D::TableDataHolder const tableConverted =
      tableData2D.convertTable2D( coordinates,
                                  tableFunction.getCoordsDescription( 0, false ),
                                  tableFunction.getCoordsDescription( 1, false ),
                                  values,
                                  false,
                                  tableFunction.getValuesDescription() );

    tableLayout.addColumns( tableConverted.headerNames );

    TableCSVFormatter const csvFormat( tableLayout );
    return csvFormat.toString( tableConverted.tableData );
  }
}

template<>
string TableTextFormatter::toString< TableFunction >( TableFunction const & tableFunction ) const
{
  ArrayOfArraysView< real64 const > coordinates = tableFunction.getCoordinates();
  arrayView1d< real64 const > const values = tableFunction.getValues();
  integer const numDimensions = LvArray::integerConversion< integer >( coordinates.size() );

  string_view tableTitle = m_tableLayout.getTitleStr();
  TableLayout::Column parentColumn = TableLayout::Column().
                                       setName( tableFunction.getTableDescription() ).
                                       setHeaderAlignment( TableLayout::Alignment::left );
  if( !m_tableLayout.getColumns().empty() )
  { // concatainating description parent column if existing
    parentColumn.m_header.setText( GEOS_FMT( "{}\n{}",
                                             parentColumn.m_header.getText(),
                                             m_tableLayout.getColumns()[0].m_header.getText() ) );
  }

  string logOutput;
  {
    if( numDimensions != 2 )
    {
      for( int i = 0; i < numDimensions; ++i )
      {
        bool const shortenDescription = tableFunction.getDimUnit( i ) == units::Unknown;
        parentColumn.addSubColumn( tableFunction.getCoordsDescription( i, shortenDescription ) );
      }
      parentColumn.addSubColumn( tableFunction.getValuesDescription() );

      TableLayout const tableLayout( tableTitle, { parentColumn } );
      TableTextFormatter const logTable( tableLayout );

      TableData tableData;
      collectTableValues( tableData, numDimensions, coordinates, values );
      logOutput = logTable.toString( tableData );
    }
    else if( numDimensions == 2 )
    {
      TableData2D tableData2D;
      TableData2D::TableDataHolder tableConverted;
      tableConverted = tableData2D.convertTable2D( coordinates,
                                                   tableFunction.getCoordsDescription( 1, true ),
                                                   tableFunction.getCoordsDescription( 0, true ),
                                                   values,
                                                   true,
                                                   tableFunction.getValuesDescription() );

      parentColumn.addSubColumns( tableConverted.headerNames );
      TableLayout const tableLayout = TableLayout( tableTitle, { parentColumn } ).
                                        setMargin( TableLayout::MarginValue::small );
      TableTextFormatter const table2DLog( tableLayout );
      logOutput =  table2DLog.toString( tableConverted.tableData );
    }
  }
  return logOutput;
}

REGISTER_CATALOG_ENTRY( FunctionBase, TableFunction, string const &, Group * const )

} // end of namespace geos
