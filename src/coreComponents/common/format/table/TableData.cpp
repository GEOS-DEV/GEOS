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
 * @file TableData.cpp
 */

#include "TableData.hpp"
#include "common/logger/Logger.hpp"

namespace geos
{

void TableData::addRow( std::vector< TableData::CellData > const & row )
{
  m_rows.push_back( row );
}

void TableData::addSeparator()
{
  if( m_rows.empty())
  {
    m_errors->addError( "Warning : Bad use of a Tabledata::addSeparator(). Make sure you have added values in TableData" );
  }
  else
  {
    integer rowSize = m_rows.front().size();
    m_rows.emplace_back( std::vector< TableData::CellData >( rowSize, { CellType::Separator, "-" } ));
  }


}

void TableData::clear()
{
  m_rows.clear();
  getErrorsList().clear();
}

std::vector< std::vector< TableData::CellData > > const & TableData::getTableDataRows() const
{
  return m_rows;
}

void TableData2D::collectTableValues( arraySlice1d< real64 const > dim0AxisCoordinates,
                                      arraySlice1d< real64 const > dim1AxisCoordinates,
                                      arrayView1d< real64 const > values,
                                      bool columnMajorInputValues )
{
  arraySlice1d< real64 const > rowAxisCoordinates = columnMajorInputValues ? dim1AxisCoordinates : dim0AxisCoordinates;
  arraySlice1d< real64 const > columAxisCoordinates = columnMajorInputValues ? dim0AxisCoordinates : dim1AxisCoordinates;
  integer const nCol = columAxisCoordinates.size();
  integer const nRow = rowAxisCoordinates.size();

  array1d< real64 > wellConstructValues( values.size() );

  if( nRow * nCol != values.size())
  {
    wellConstructValues = values;
    if( nRow * nCol > values.size())
    {
      m_errors->addError( GEOS_FMT( "Warning : The number of values is insufficient for the number of columns * rows:\n"
                                    "Expected {} values, Found {} values", nRow * nCol, values.size() ) );
      wellConstructValues.resizeDefault( nRow * nCol, 0 );
    }
    else
    {
      m_errors->addError( GEOS_FMT( "Warning : Too much data for the number of columns * rows:\n"
                                    "Expected {} values, Found {} values\n"
                                    "Data may be missaligned", nRow * nCol, values.size() ) );
    }
  }

  for( integer y = 0; y < nRow; y++ )
  {
    for( integer x = 0; x < nCol; x++ )
    {
      addCell( rowAxisCoordinates[y], columAxisCoordinates[x], wellConstructValues[ x + y*nCol ] );
    }
  }
}

TableData2D::TableDataHolder TableData2D::convertTable2D( ArrayOfArraysView< real64 const > const coordinates,
                                                          string_view rowAxisDescription,
                                                          string_view columnAxisDescription,
                                                          arrayView1d< real64 const > const values,
                                                          bool columnMajorValues,
                                                          string_view valueDescription )
{
  string const rowFmt = GEOS_FMT( "{} = {{}}", rowAxisDescription );
  string const columnFmt = GEOS_FMT( "{} = {{}}", columnAxisDescription );
  collectTableValues( coordinates[0], coordinates[1], values, columnMajorValues );
  return buildTableData( valueDescription, rowFmt, columnFmt );
}

TableData2D::TableDataHolder TableData2D::buildTableData( string_view targetUnit,
                                                          string_view rowFmt,
                                                          string_view columnFmt ) const
{
  TableData2D::TableDataHolder tableData1D;

  tableData1D.headerNames.push_back( string( targetUnit ) );

  for( auto const & columnValue : m_columnValues )
  {
    tableData1D.headerNames.push_back( GEOS_FMT( columnFmt, columnValue ) );
  }

  for( auto const & error : m_errors->errorText )
  {
    tableData1D.tableData.getErrorsList().addError( error );
  }

  // insert row value and row cell values
  for( auto const & [rowValue, rowMap] : m_data )
  {
    std::vector< TableData::CellData > currentRowValues;
    currentRowValues.reserve( rowMap.size() );
    currentRowValues.push_back( {CellType::Value, GEOS_FMT( rowFmt, rowValue )} );

    std::set< real64 >::const_iterator columnIt = m_columnValues.begin();
    for( auto const & [columnValue, cellValue] : rowMap )
    {
      // if a column value(s) is/are missing, insert empty entry(ies)
      while( columnValue > *( columnIt++ ) && columnIt != m_columnValues.end() )
      {
        currentRowValues.push_back( {CellType::Value, ""} );
      }

      currentRowValues.push_back( {CellType::Value, GEOS_FMT( "{}", cellValue )} );
    }

    tableData1D.tableData.addRow( currentRowValues );
  }

  return tableData1D;
}
}
