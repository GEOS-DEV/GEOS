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

#include "TableFunction.hpp"
#include "LogLevelsInfo.hpp"
#include "codingUtilities/Parsing.hpp"
#include "common/DataTypes.hpp"
#include "common/MpiWrapper.hpp"

#include <algorithm>

namespace geos
{

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
  size_t const numDims = size_t( numDimensions() );
  std::vector< string > labels;
  labels.reserve( numDims + 1 );

  if( !m_coordinateFiles.empty() )
  {
    for( size_t i = 0; i < numDims; ++i )
      labels.push_back( GEOS_FMT( "Coordinates {} from file: {}",
                                  getCoordsDescription( i ), m_coordinateFiles[i].relativeFilePath() ) );
  }
  else
  {
    for( size_t i = 0; i < numDims; ++i )
      labels.push_back( GEOS_FMT( "{}, the {}", getCoordsDescription( i ), units::getDescription( getDimUnit( i ) ) ) );
  }

  labels.push_back( !m_voxelFile.empty() ?
                    GEOS_FMT( "Values from file: {}", m_voxelFile.relativeFilePath() ) :
                    GEOS_FMT( "Values of {}", units::getDescription( getValueUnit() ) ) );

  return stringutilities::join( labels, "\n" );
}

string TableFunction::getCoordsDescription( integer dimId ) const
{
  integer const numDims = numDimensions();
  units::Unit const dimCoordsUnit = getDimUnit( dimId );
  if( dimCoordsUnit != units::Unknown && dimCoordsUnit != units::Dimensionless )
  {
    return string( units::getSymbol( dimCoordsUnit ) );
  }
  else if( numDims==1 )
  {
    return "coordinates";
  }
  else if( numDims<=3 )
  {
    static const string labels[] = {"X", "Y", "Z"};
    return labels[dimId];
  }
  else
  {
    return GEOS_FMT( "coord[{}]", dimId );
  }
}

string TableFunction::getValuesDescription() const
{
  units::Unit const valuesUnits = getValueUnit();
  return valuesUnits != units::Unknown ? GEOS_FMT( "values in\n{}", units::getDescription( valuesUnits ) ) :
         "values";
}

/**
 * @brief Retrieve all data headers from a table function
 * @param formatterStream The stream who contains the csv table string
 * @param tableFunction The table function to be process
 * @param numDimensions Numbers of axes in the table
 */
void collectHeaders( std::ostringstream & formatterStream,
                     TableFunction const & tableFunction,
                     integer const numDimensions )
{
  for( integer d = 0; d < numDimensions; d++ )
  {
    formatterStream << tableFunction.getCoordsDescription( d ) << ",";
  }
  formatterStream << tableFunction.getValuesDescription() << "\n";
}

/**
 * @brief Retrieve all data values
 * @param formatterStream  The stream who contains the csv table string
 * @param numDimensions Numbers of axes in the table
 * @param coordinates The tables axis values
 * @param values The table values to be retrived
 */
void collectValues( std::ostringstream & formatterStream,
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
  for( integer v = 0; v < values.size(); v++ )
  {
    // find coords indices
    std::vector< integer > idx( numDimensions );
    integer r = v;
    for( integer d = numDimensions-1; d >= 0; d-- )
    {
      idx[d] = r / div[d];
      r = r % div[d];
    }
    // finally print out in right order
    for( integer d = 0; d < numDimensions; d++ )
    {
      arraySlice1d< real64 const > const coords = coordinates[d];
      formatterStream << coords[idx[d]] << ",";
    }
    formatterStream << values[v] << "\n";
  }
}

void TableFunction::outputTableData( OutputOptions const outputOpts ) const
{
  // we only output from rank 0
  if( MpiWrapper::commRank() != 0 )
    return;

  if( outputOpts.writeInLog &&  this->numDimensions() <= 2 )
  {
    TableTextFormatter textFormatter;
    GEOS_LOG( textFormatter.toString( *this ));
  }
  if( outputOpts.writeCSV || ( outputOpts.writeInLog && this->numDimensions() >= 3 ) )
  {
    string const filename = this->getName();
    std::ofstream logStream( joinPath( FunctionBase::getOutputDirectory(), filename + ".csv" ) );
    GEOS_LOG( GEOS_FMT( "CSV Generated to {}/{}.csv \n",
                        FunctionBase::getOutputDirectory(),
                        filename ));
    TableCSVFormatter csvFormatter;
    logStream << csvFormatter.toString( *this );
  }
}

void TableFunction::initializePostSubGroups()
{
  // Output user defined tables (not generated PVT tables)
  outputTableData( OutputOptions{
      m_writeCSV != 0, // writeCSV
      isLogLevelActive< logInfo::TableDataOutput >( getLogLevel() ) // writeInLog
    } );
}

template<>
string TableCSVFormatter::toString< TableFunction >( TableFunction const & tableFunction ) const
{
  ArrayOfArraysView< real64 const > const coordinates = tableFunction.getCoordinates();
  arrayView1d< real64 const > const values = tableFunction.getValues();
  std::ostringstream formatterStream;

  integer const numDimensions = LvArray::integerConversion< integer >( coordinates.size() );
  if( numDimensions != 2 )
  {
    GEOS_LOG( GEOS_FMT( "[START TABLE {}] CSVing 1/ND table", tableFunction.getName()));
    collectHeaders( formatterStream, tableFunction, numDimensions );
    collectValues( formatterStream, numDimensions, coordinates, values );
  }
  else
  {
    GEOS_LOG( GEOS_FMT( "[START TABLE {}] CSVing 2D table", tableFunction.getName()));
    TableData2D tableData2D;
    TableData2D::TableDataHolder tableConverted;
    tableConverted = tableData2D.convertTable2D( coordinates,
                                                 tableFunction.getCoordsDescription( 0 ),
                                                 tableFunction.getCoordsDescription( 1 ),
                                                 values,
                                                 false,
                                                 tableFunction.getValuesDescription() );

    TableLayout const tableLayout( "", tableConverted.headerNames );

    TableCSVFormatter csvFormat( tableLayout );
    formatterStream << csvFormat.headerToString() << csvFormat.dataToString( tableConverted.tableData );
  }
  GEOS_LOG( GEOS_FMT( "[END TABLE {}]", tableFunction.getName()));
  return formatterStream.str();
}

template<>
string TableTextFormatter::toString< TableFunction >( TableFunction const & tableFunction ) const
{
  static constexpr integer maxRows = 500;
  ArrayOfArraysView< real64 const > coordinates = tableFunction.getCoordinates();
  arrayView1d< real64 const > const values = tableFunction.getValues();
  integer const numDimensions = LvArray::integerConversion< integer >( coordinates.size() );
  string const tableTitle = GEOS_FMT( "{}\n{}", tableFunction.getName(), tableFunction.getTableDescription() );
  string logOutput;

  if( numDimensions == 1 && coordinates[0].size() < maxRows )
  {
    GEOS_LOG( GEOS_FMT( "[START TABLE {}] drawing 1D table", tableFunction.getName()));
    TableData tableData;
    arraySlice1d< real64 const > const coords = coordinates[0];
    for( integer idx = 0; idx < values.size(); idx++ )
    {
      tableData.addRow( coords[idx], values[idx] );
    }
    TableLayout const tableLayout( tableTitle, {
        string( tableFunction.getCoordsDescription( 0 ) ),
        string( tableFunction.getValuesDescription() )
      } );
    TableTextFormatter const logTable( tableLayout );
    logOutput = logTable.toString( tableData );
  }
  else if( numDimensions == 2 && ( coordinates[0].size() * coordinates[1].size() ) < maxRows )
  {
    GEOS_LOG( GEOS_FMT( "[START TABLE {}] drawing 2D table", tableFunction.getName()));
    for (int i=0;i<=1;++i) {
    GEOS_LOG( "Column major mode = "<<i );
    TableData2D tableData2D;
    TableData2D::TableDataHolder tableConverted;
    tableConverted = tableData2D.convertTable2D( coordinates,
                                                tableFunction.getCoordsDescription( 1 ),
                                                tableFunction.getCoordsDescription( 0 ),
                                                values,
                                                bool(i),
                                                tableFunction.getValuesDescription() );

    TableLayout const tableLayout = TableLayout( tableTitle, tableConverted.headerNames ).
                                      setMargin( TableLayout::MarginValue::small );
    TableTextFormatter const table2DLog( tableLayout );
    logOutput =  table2DLog.toString( tableConverted.tableData );
    GEOS_LOG(logOutput);//test
    }
  }
  else
  {
    GEOS_LOG( GEOS_FMT( "[START TABLE {}] drawing out of context table", tableFunction.getName()));
    string const tooLongOutputMsg = GEOS_FMT( "The table is too heavy for log output.\n"
                                              "To visualize the table, please refer to the generated csv.",
                                              maxRows );
    TableLayout const tableLayoutInfos( tableTitle, {tooLongOutputMsg} );
    TableTextFormatter const tableLog( tableLayoutInfos );
    logOutput = tableLog.toString();
  }
  GEOS_LOG( GEOS_FMT( "[END TABLE {}]", tableFunction.getName()));
  return logOutput;
}

REGISTER_CATALOG_ENTRY( FunctionBase, TableFunction, string const &, Group * const )

} // end of namespace geos
