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
 * @file ConstitutiveDriver.cpp
 */

#include "ConstitutiveDriver.hpp"
#include "LogLevelsInfo.hpp"

#include "constitutive/ConstitutiveManager.hpp"

#include "common/format/table/TableFormatter.hpp"
#include "common/format/StringUtilities.hpp"
#include "common/MpiWrapper.hpp"

namespace geos
{

using namespace dataRepository;
using namespace constitutive;

ConstitutiveDriver::ConstitutiveDriver( const string & name,
                                        Group * const parent ):
  TaskBase( name, parent )
{
  registerWrapper( viewKeyStruct::numStepsString(), &m_numSteps ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Number of steps to take" );

  registerWrapper( viewKeyStruct::outputString(), &m_outputFile ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Output file" );

  registerWrapper( viewKeyStruct::precisionString(), &m_precision ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_precision ).
    setDescription( "The precision to use to data out to files" );

  addLogLevel< logInfo::LogOutput >();

  registerWrapper( viewKeyStruct::outputDataString(), &m_table ).
    setInputFlag( InputFlags::FALSE ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Time-history of data" );
}

void ConstitutiveDriver::postInputInitialization()
{
  bool const is_root = MpiWrapper::commRank() == 0;

  GEOS_THROW_IF_LT_MSG( m_numSteps, 1,
                        "Number of steps must be at least 1",
                        InputError );

  GEOS_WARNING_IF( ( m_precision < minPrecision) && is_root,
                   GEOS_FMT( "{}: option should be between {} and {}. A value of {} will be used.",
                             getWrapperDataContext( viewKeyStruct::precisionString() ),
                             minPrecision, maxPrecision, minPrecision ));

  GEOS_WARNING_IF( ( maxPrecision < m_precision) && is_root,
                   GEOS_FMT( "{}: option should be between {} and {}. A value of {} will be used.",
                             getWrapperDataContext( viewKeyStruct::precisionString() ),
                             minPrecision, maxPrecision, maxPrecision ) );
}

bool ConstitutiveDriver::execute( real64 const GEOS_UNUSED_PARAM( time_n ),
                                  real64 const GEOS_UNUSED_PARAM( dt ),
                                  integer const GEOS_UNUSED_PARAM( cycleNumber ),
                                  integer const GEOS_UNUSED_PARAM( eventCounter ),
                                  real64 const GEOS_UNUSED_PARAM( eventProgress ),
                                  DomainPartition & GEOS_UNUSED_PARAM( domain ) )
{
  // This code only makes sense in serial
  GEOS_THROW_IF( MpiWrapper::commRank() > 0, "Constitutive Driver should only be run in serial", InputError, getDataContext() );

  bool executeResult = execute();

  // move table back to host for output
  m_table.move( hostMemorySpace );

  validateResults();

  outputResults();

  return executeResult;
}

bool ConstitutiveDriver::validateResults()
{
  return true;
}

void ConstitutiveDriver::outputResults() const
{
  if( m_outputFile.empty())
  {
    outputToConsole();
  }
  else
  {
    outputToFile();
  }
}

void ConstitutiveDriver::outputToFile() const
{
  string const outputDir = TaskBase::getOutputDirectory();
  string const outputPath = joinPath( outputDir, m_outputFile );
  FILE * fp = fopen( outputPath.c_str(), "w" );

  string_array columnNames;
  getColumnNames( columnNames );
  integer const numColumns = static_cast< integer >( columnNames.size() );
  for( integer col = 0; col < numColumns; ++col )
  {
    fprintf( fp, "# column %d = %s\n", col+1, columnNames[col].c_str() );
  }

  integer const precision = LvArray::math::max( LvArray::math::min( m_precision, maxPrecision ), minPrecision );
  string const format = GEOS_FMT( "%{}.{}e ", precision+7, precision );

  for( integer step = 0; step <= m_numSteps; ++step )
  {
    for( integer col = 0; col < numColumns; ++col )
    {
      fprintf( fp, format.c_str(), m_table( step, col ) );
    }
    fprintf( fp, "\n" );
  }

  fclose( fp );
}

void ConstitutiveDriver::outputToConsole() const
{
  string_array columnNames;
  getColumnNames( columnNames );
  integer const numColumns = static_cast< integer >( columnNames.size() );

  integer const precision = LvArray::math::max( LvArray::math::min( m_precision, maxPrecision ), minPrecision );
  string const format = GEOS_FMT( "{{:{}.{}e}}", precision+7, precision );

  stdVector< TableLayout::Column > columns;
  stdMap< std::string, size_t > uniqueColumnNames;
  auto const getColumn = [&uniqueColumnNames, &columns]( string const & name ) -> TableLayout::Column *
  {
    size_t location = columns.size();
    auto [it, inserted] = uniqueColumnNames.try_emplace( name, location );
    if( inserted )
    {
      columns.emplace_back( name );
    }
    else
    {
      location = it->second;
    }
    return &(columns[location]);
  };

  for( integer col = 0; col < numColumns; ++col )
  {
    stdVector< string > titles = stringutilities::tokenize( columnNames[col], ",", false );
    size_t const titleCount = titles.size();

    TableLayout::Column * currentColumn = getColumn( titles[0] );
    for( size_t ti = 1; ti < titleCount; ti++ )
    {
      auto & subColumns = currentColumn->m_subColumns;
      auto colIndex = std::find_if( subColumns.begin(), subColumns.end(), [&]( TableLayout::Column const & column ) {
        return column.getName() == titles[ti];
      } );
      if( colIndex == subColumns.end())
      {
        currentColumn->addSubColumn( TableLayout::Column( titles[ti] ));
        currentColumn = &(subColumns.back());
      }
      else
      {
        currentColumn = &(*colIndex);
      }
    }
  }

  TableData tableData;

  stdVector< TableData::CellData > tableRow( numColumns );
  for( integer step = 0; step <= m_numSteps; ++step )
  {
    for( integer col = 0; col < numColumns; ++col )
    {
      tableRow[col].value = GEOS_FMT( format, m_table( step, col ) );
    }
    tableData.addRow( tableRow );
  }

  string const tableTitle = GEOS_FMT( "Output for {}", getName());
  TableLayout const tableLayout( tableTitle, columns );
  TableTextFormatter const tableText( tableLayout );

  std::cout << tableText.toString( tableData ) << std::endl;
}

void ConstitutiveDriver::allocateTable( integer const numColumns,
                                        real64 const minTime,
                                        real64 const maxTime )
{
  // resize data table to fit number of timesteps:
  m_table.resize( m_numSteps+1, numColumns );
  real64 const dt = (maxTime-minTime) / m_numSteps;

  // set time columns
  for( integer step = 0; step <= m_numSteps; ++step )
  {
    m_table( step, TIME ) = minTime + step*dt;
  }
}

ConstitutiveManager & ConstitutiveDriver::getConstitutiveManager()
{
  return this->getGroupByPath< ConstitutiveManager >( "/Problem/domain/Constitutive" );
}

ConstitutiveManager const & ConstitutiveDriver::getConstitutiveManager() const
{
  return this->getGroupByPath< ConstitutiveManager >( "/Problem/domain/Constitutive" );
}

} /* namespace geos */
