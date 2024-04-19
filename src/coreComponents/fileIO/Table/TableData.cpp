/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file TableData.cpp
 */

#include "TableData.hpp"

namespace geos
{

void TableData::addRow( std::vector< string > const & row )
{
  if( m_rows.size() != 0 && row.size() != m_rows[m_rows.size() - 1].size() )
  {
    string msg = "Remarks : some cells may be missing";
    if( std::find( m_errorsMsg.begin(), m_errorsMsg.end(), msg ) == m_errorsMsg.end())
    {
      m_errorsMsg.push_back( msg );
    }
  }
  m_rows.push_back( row );
}

void TableData::clear()
{
  m_rows.clear();
}

std::vector< std::vector< string > > const & TableData::getTableDataRows() const
{
  return m_rows;
}

std::vector< string > const & TableData::getErrorMsgs() const
{
  return m_errorsMsg;
}

void TableData::addErrorMsgs( string const & msg )
{
  std::vector< string > splitHeaderParts;
  std::istringstream ss( msg );
  string splitErrors;

  while( std::getline( ss, splitErrors, '\n' ))
  {
    m_errorsMsg.push_back( splitErrors );
  }
}

TableData2D::Conversion1D TableData2D::buildTableData( string_view targetUnit,
                                                       string_view rowFmt,
                                                       string_view columnFmt ) const
{
  TableData2D::Conversion1D tableData1D;
  std::vector< real64 > headerValues;
  std::vector< size_t > rowsLength;

  tableData1D.headerNames.push_back( string( targetUnit ) );

  for( auto const & columnValue : m_columnValues )
  {
    tableData1D.headerNames.push_back( GEOS_FMT( columnFmt, columnValue ) );
    headerValues.push_back( columnValue );
  }

  // insert row value and row cell values
  for( auto const & [rowValue, rowMap] : m_data )
  {
    std::vector< string > currentRowValues;
    currentRowValues.reserve( rowMap.size() );
    currentRowValues.push_back( GEOS_FMT( rowFmt, rowValue ) );

    integer idxColumn = 0;
    for( auto const & [columnValue, cellValue] : rowMap )
    {
      //check if column numbers in the evaluated row is consistent
      while( columnValue > headerValues[idxColumn] )
      {
        currentRowValues.push_back( "" );
        ++idxColumn;
      }
      currentRowValues.push_back( GEOS_FMT( "{}", cellValue ) );
      ++idxColumn;
    }

    tableData1D.tableData.addRow( std::move( currentRowValues ) );
    rowsLength.push_back( currentRowValues.size() );
  }

  return tableData1D;
}
}
