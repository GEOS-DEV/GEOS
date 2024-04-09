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

TableData2D::Conversion1D TableData2D::buildTableData( string_view targetUnit,
                                                       string_view rowFmt,
                                                       string_view columnFmt ) const
{
  TableData2D::Conversion1D tableData1D;
  std::vector< real64 > headerValues;
  std::vector< size_t > rowsLength;

  tableData1D.headerNames.push_back( string( targetUnit ) );

  // looping over first line to fill columnNames
  for( auto const & [ columnValue, cellValue] : m_data.begin()->second )
  {
    tableData1D.headerNames.push_back( GEOS_FMT( columnFmt, columnValue ) );
    headerValues.push_back( columnValue );
  }

  rowsLength.reserve( headerValues.size());

  // insert row value and row cell values
  bool flag = true;
  for( auto const & [rowValue, rowMap] : m_data )
  {
    integer i = 0;
    std::vector< string > currentRowValues;
    currentRowValues.push_back( GEOS_FMT( rowFmt, rowValue ) );
    for( auto const & [columnValue, cellValue] : rowMap )
    {
      //check is if the column are offset
      if( std::abs( columnValue - headerValues[i] ) > 0.0 )
      {
        flag = false;
      }
      currentRowValues.push_back( GEOS_FMT( "{}", cellValue ) );
      ++i;
    }
    tableData1D.tableData.addRow( currentRowValues );
    rowsLength.push_back( currentRowValues.size());
  }

  if( std::adjacent_find( rowsLength.begin(), rowsLength.end(), std::not_equal_to<>() ) != rowsLength.end() )
  {
    flag = false;
    GEOS_WARNING( "Cell(s) are missing in row" );
  }

  if( !flag )
  {
    tableData1D.isConsistent = flag;
    GEOS_WARNING( "Table isn't consistent" );
  }

  return tableData1D;
}
}
