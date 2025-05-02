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

#include <algorithm>

#include <hdf5.h>

namespace geos
{

namespace hdf5Utils
{
hid_t openHDFFile( const string & filename )
{
  hid_t file_id;

  H5E_BEGIN_TRY
  {
    file_id = H5Fopen( filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );
  }
  H5E_END_TRY

  GEOS_ERROR_IF( file_id < 0, GEOS_FMT( "{}: error opening file", filename ) );
  return file_id;
}

void closeHDFFile( hid_t const file_id, const string & filename )
{
  herr_t err;

  H5E_BEGIN_TRY
  {
    err = H5Fclose( file_id );
  }
  H5E_END_TRY

  GEOS_ERROR_IF( err < 0, GEOS_FMT( "{}: error closing file", filename ) );
}

bool datasetExists( hid_t const & file_id, string const & datasetName )
{
  herr_t err;

  H5E_BEGIN_TRY
  {
    err = H5Oexists_by_name( file_id, datasetName.c_str(), H5P_DEFAULT );
  }
  H5E_END_TRY
  return err > 0 ? true : false;
}

array1d< real64 > readDataset( hid_t file_id, string const & filename, string const & datasetName, int const expectedNumDims )
{
  // Check dataset existance
  GEOS_ERROR_IF( !datasetExists( file_id, datasetName ), GEOS_FMT( "Dataset {} doesn't exist in {}", datasetName, filename ) );

  // Open the dataset
  hid_t dataset_id = H5Dopen( file_id, datasetName.c_str(), H5P_DEFAULT );
  GEOS_ERROR_IF( dataset_id < 0, GEOS_FMT( "Error opening dataset {} in {}", datasetName, filename ) );

  // Get the dataspace of the dataset
  hid_t dataspace_id = H5Dget_space( dataset_id );

  // Get the dataset dimensions
  int ndims = H5Sget_simple_extent_ndims( dataspace_id );
  if( ndims != expectedNumDims )
  {
    H5Sclose( dataspace_id );
    H5Dclose( dataset_id );
    H5Fclose( file_id );
    GEOS_ERROR( GEOS_FMT( "Dataset {} is not {}D", datasetName, std::to_string( expectedNumDims ) ) );
  }

  array1d< hsize_t >  dims( ndims );
  H5Sget_simple_extent_dims( dataspace_id, dims.data(), nullptr );

  // Compute the total number of entries
  hsize_t total_entries = 1;
  for( auto const & dim : dims )
  {
    total_entries *= dim;
  }

  // Get the dataset type
  hid_t dataset_type = H5Dget_type(dataset_id);

  // Read the dataset
  array1d< real64 > datasetValues( total_entries );
  herr_t err{};

  if ( (H5Tequal( dataset_type, H5T_NATIVE_FLOAT) ) )
  {
    array1d< real32 > tmp( total_entries );
    H5E_BEGIN_TRY
    {
      err = H5Dread( dataset_id, H5T_NATIVE_FLOAT, dataspace_id, H5S_ALL, H5P_DEFAULT, tmp.data() );
    }
    H5E_END_TRY
    GEOS_ERROR_IF( err < 0, GEOS_FMT( "{}: error reading single-precision dataset", datasetName ) );
    for (int i = 0; i < tmp.size(); ++i)
    {
      datasetValues[i] = static_cast<double>( tmp[i] ) ;
    }
  }
  else if ( (H5Tequal( dataset_type, H5T_NATIVE_DOUBLE) ) )
  {
    H5E_BEGIN_TRY
    {
      err = H5Dread( dataset_id, H5T_NATIVE_DOUBLE, dataspace_id, H5S_ALL, H5P_DEFAULT, datasetValues.data() );
    }
    H5E_END_TRY
    GEOS_ERROR_IF( err < 0, GEOS_FMT( "{}: error reading double-precision dataset", datasetName ) );
  }
  else
  {
    GEOS_ERROR( GEOS_FMT("Dataset {} is not of type float", datasetName) );
  }

  GEOS_ERROR_IF( err < 0, GEOS_FMT( "{}: error reading dataset", datasetName ) );

  H5Tclose( dataset_type );
  H5Sclose( dataspace_id );
  H5Dclose( dataset_id );

  return datasetValues;
}

} // end of namespace hdf5Utils

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

  registerWrapper( viewKeyStruct::hdf5FileString(), &m_hdf5File ).
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
                                         std::vector< units::Unit > const & dimUnits )
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
  // Read in data
  if( m_coordinates.size() > 0 )
  {
    // This function appears to be already initialized
    // Apparently, this can be called multiple times during unit tests?
  }
  else if( !m_hdf5CoordinateDatasetNames.empty() )
  {
    hid_t file_id = hdf5Utils::openHDFFile( m_hdf5File );

    // Read table coordinates
    int tableDim = 0;
    for( auto const & coordName : m_hdf5CoordinateDatasetNames )
    {
      array1d< real64 > tmpNew = hdf5Utils::readDataset( file_id, m_hdf5File, coordName, 1 );
      m_coordinates.appendArray( tmpNew.begin(), tmpNew.end() );
      tableDim += 1;
    }

    // Read table dataset
    m_values = hdf5Utils::readDataset( file_id, m_hdf5File, m_hdf5TableDatasetName, tableDim );

    hdf5Utils::closeHDFFile( file_id, m_hdf5File );
  }
  else if( m_coordinateFiles.empty() )
  {
    // 1D Table
    m_coordinates.appendArray( m_tableCoordinates1D.begin(), m_tableCoordinates1D.end() );
    GEOS_THROW_IF_NE_MSG( m_tableCoordinates1D.size(), m_values.size(),
                          GEOS_FMT( "{} {}: 1D table function coordinates and values must have the same length",
                                    catalogName(), getDataContext() ),
                          InputError );
  }
  else
  {
    array1d< real64 > tmp;
    localIndex numValues = 1;
    for( localIndex ii = 0; ii < m_coordinateFiles.size(); ++ii )
    {
      tmp.clear();
      readFile( m_coordinateFiles[ii], tmp );
      m_coordinates.appendArray( tmp.begin(), tmp.end() );
      numValues *= tmp.size();
    }
    // ND Table
    m_values.reserve( numValues );
    readFile( m_voxelFile, m_values );
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
  std::vector< integer > div( numDimensions );
  div[0] = 1;
  for( integer d = 1; d < numDimensions; d++ )
  {
    div[d] = div[d-1] * coordinates[d-1].size();
  }
  // loop through all the values
  std::vector< integer > coordsIdx( numDimensions );
  std::vector< TableData::CellData > rowData( numDimensions + 1,
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
