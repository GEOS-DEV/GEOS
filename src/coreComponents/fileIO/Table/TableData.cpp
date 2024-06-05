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
  m_errorsMsg.clear();
}

std::vector< std::vector< string > > const & TableData::getTableDataRows() const
{
  return m_rows;
}

std::vector< string > const & TableData::getErrorMsgs() const
{
  return m_errorsMsg;
}

void TableData2D::collect2DData( arraySlice1d< real64 const > const rowAxis,
                                 arraySlice1d< real64 const > const columnAxis,
                                 arrayView1d< real64 const > values )
{
  arraySlice1d< real64 const > const coordsX = rowAxis;
  arraySlice1d< real64 const > const coordsY = columnAxis;
  integer const nX = coordsX.size();
  integer const nY = coordsY.size();

  for( integer i = 0; i < nX; i++ )
  {
    for( integer y = 0; y < nY; y++ )
    {
      addCell( rowAxis[i], columnAxis[y], values[ y*nX + i ] );
    }
  }
}

TableData2D::Conversion1D TableData2D::convert2DData( units::Unit valueUnit,
                                                      string_view rowUnitDescription,
                                                      string_view columnUnitDescription )
{
  string const rowFmt = GEOS_FMT( "{} = {{}}", rowUnitDescription );
  string const columnFmt = GEOS_FMT( "{} = {{}}", columnUnitDescription );
  return buildTableData( string( units::getDescription( valueUnit )),
                         rowFmt,
                         columnFmt );
}

size_t TableData2D::getNbRows() const
{
  return m_data.size();
}

TableData2D::Conversion1D TableData2D::buildTableData( string_view targetUnit,
                                                       string_view rowFmt,
                                                       string_view columnFmt ) const
{
  TableData2D::Conversion1D tableData1D;
  std::vector< size_t > rowsLength;

  tableData1D.headerNames.push_back( string( targetUnit ) );

  for( auto const & columnValue : m_columnValues )
  {
    tableData1D.headerNames.push_back( GEOS_FMT( columnFmt, columnValue ) );
  }

  // insert row value and row cell values
  for( auto const & [rowValue, rowMap] : m_data )
  {
    std::vector< string > currentRowValues;
    currentRowValues.reserve( rowMap.size() );
    currentRowValues.push_back( GEOS_FMT( rowFmt, rowValue ) );

    std::set< real64 >::const_iterator columnIt = m_columnValues.begin();
    for( auto const & [columnValue, cellValue] : rowMap )
    {
      // if a column value(s) is/are missing, insert empty entry(ies)
      while( columnValue > *( columnIt++ ) && columnIt != m_columnValues.end() )
      {
        currentRowValues.push_back( "" );
      }
      currentRowValues.push_back( GEOS_FMT( "{}", cellValue ) );
    }

    tableData1D.tableData.addRow( std::move( currentRowValues ) );
    rowsLength.push_back( currentRowValues.size() );
  }

  return tableData1D;
}
}
