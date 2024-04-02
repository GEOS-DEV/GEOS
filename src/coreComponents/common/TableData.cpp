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

#include "common/TableData.hpp"

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

std::set< real64 > const & TableData2D::getColumns() const
{
  return m_columns;
}
std::set< real64 > const & TableData2D::getRows() const
{
  return m_rows;
}

TableData TableData2D::buildTableData( std::vector< string > & columnNames,
                                       string_view rowFmt, string_view columnFmt ) const
{
  TableData tableDataToBeBuilt;
  std::set< real64 > columnValues;

  { // insert row value and row cell values
    std::vector< string > currentRowValues;
    RowType currentRow = m_data.begin()->first.first;
    currentRowValues.push_back( GEOS_FMT( rowFmt, currentRow ) );
    for( auto const & [rowColumnPair, cellValue] : m_data )
    {
      if( rowColumnPair.first == currentRow )
      {
        currentRowValues.push_back( GEOS_FMT( "{}", cellValue ) );
        columnValues.insert( rowColumnPair.second );
      }
      else // we are changing line
      {
        tableDataToBeBuilt.addRow( currentRowValues );
        currentRowValues.clear();
        currentRowValues.push_back( GEOS_FMT( rowFmt, currentRow ) );
        firstRow = false;
      }
    }
  }

  // fill columnNames
  std::transform( columnValues.begin(), columnValues.end(),
                  std::back_inserter( columnValues, columnValues.begin() ),
                  [&] ( real64 const columnValue ) {
      return GEOS_FMT( columnFmt, columnValue );
    } );

  return tableDataToBeBuilt;
}

}
