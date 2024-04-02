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

  // looping over first line to fill columnNames
  columnNames.clear();
  for( auto const & [columnValue, GEOS_UNUSED_VAR( cellValue )] : m_data.begin()->second )
  {
    columnNames.push_back( GEOS_FMT( columnFmt, columnValue ) );
    ++columnCount;
  }

  // insert row value and row cell values
  for( auto const & [rowValue, rowMap] : m_data )
  {
    std::vector< string > currentRowValues;
    currentRowValues.push_back( GEOS_FMT( rowFmt, rowValue ) );
    for( auto const & [columnValue, cellValue] : m_data )
    {
      currentRowValues.push_back( GEOS_FMT( "{}", cellValue ) );
    }
    tableDataToBeBuilt.addRow( currentRowValues );
  }

  return tableDataToBeBuilt;
}

}
