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
    if( std::find( errorsMsg.begin(), errorsMsg.end(), msg ) == errorsMsg.end())
    {
      errorsMsg.push_back( msg );
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
  return errorsMsg;
}

void TableData::addErrorMsgs( string const & msg )
{
  std::vector< string > splitHeaderParts;
  std::istringstream ss( msg );
  string splitErrors;

  while( getline( ss, splitErrors, '\n' ))
  {
    errorsMsg.push_back( splitErrors );
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

  for( auto const & columnValue : columnValues )
  {
    tableData1D.headerNames.push_back( GEOS_FMT( columnFmt, columnValue ) );
    headerValues.push_back( columnValue );
  }

  // insert row value and row cell values
  for( auto const & [rowValue, rowMap] : m_data )
  {
    integer i = 0;
    std::vector< string > currentRowValues;
    integer idxColumn = 0;

    currentRowValues.push_back( GEOS_FMT( rowFmt, rowValue ) );
    for( auto const & [columnValue, cellValue] : rowMap )
    {
      //check if column numbers in the evaluated row is consistent
      if( columnValue  > headerValues[i] )
      {
        while( columnValue > headerValues[idxColumn] )
        {
          currentRowValues.push_back( "" );
          ++idxColumn;
        }
      }
      currentRowValues.push_back( GEOS_FMT( "{}", cellValue ) );
      ++idxColumn;
      ++i;
    }

    tableData1D.tableData.addRow( currentRowValues );
    rowsLength.push_back( currentRowValues.size());
  }

  return tableData1D;
}
}
