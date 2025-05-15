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

TableData::TableData():
  m_errors( std::make_unique< TableErrorListing >() )
{}

TableData::TableData( TableData const & other ):
  m_rows( other.m_rows ),
  m_errors( std::make_unique< TableErrorListing >( *other.m_errors ) )
{}

TableData::TableData( TableData && other ):
  m_rows( std::move( other.m_rows )),
  m_errors( std::move( other.m_errors ))
{}

TableData & TableData::operator=( TableData && other )
{
  if( this != &other )
  {
    m_rows = std::move( other.m_rows );
    m_errors = std::move( other.m_errors );
  }
  return *this;
}

TableData & TableData::operator=( TableData const & other )
{
  if( this != &other )
  {
    m_rows = other.m_rows;
    *m_errors = *other.m_errors;
  }
  return *this;
}

bool TableData::operator<( TableData const & other ) const
{
  return m_rows < other.m_rows;
}



void TableData::addRow( stdVector< TableData::CellData > const & row )
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
    m_rows.emplace_back( stdVector< TableData::CellData >( rowSize, { CellType::Separator, "-" } ));
  }

}

void TableData::clear()
{
  m_rows.clear();
  getErrorsList().clear();
}

stdVector< stdVector< TableData::CellData > > const & TableData::getTableDataRows() const
{
  return m_rows;
}

void TableData2D::collectTableValues( arrayView1d< real64 const > dim0AxisCoordinates,
                                      arrayView1d< real64 const > dim1AxisCoordinates,
                                      arrayView1d< real64 const > values,
                                      bool columnMajorInputValues )
{
  arrayView1d< real64 const > rowAxisCoordinates = columnMajorInputValues ? dim1AxisCoordinates : dim0AxisCoordinates;
  arrayView1d< real64 const > columAxisCoordinates = columnMajorInputValues ? dim0AxisCoordinates : dim1AxisCoordinates;
  integer const nCol = columAxisCoordinates.size();
  integer const nRow = rowAxisCoordinates.size();

  array1d< real64 > wellFormedValues( values.size() );

  if( nRow * nCol != values.size())
  {
    wellFormedValues = values;
    if( nRow * nCol > values.size())
    {
      m_errors->addError( GEOS_FMT( "Warning : Not enough data for the number of columns & rows:\n"
                                    "Expected {} values, Found {} values", nRow * nCol, values.size() ) );
      wellFormedValues.resizeDefault( nRow * nCol, 0 );
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
      addCell( rowAxisCoordinates[y], columAxisCoordinates[x], wellFormedValues[ x + y*nCol ] );
    }
  }
}

TableData2D::TableDataHolder TableData2D::convertTable2D( arrayView1d< real64 const > coordX, arrayView1d< real64 const > coordY,
                                                          string_view rowAxisDescription,
                                                          string_view columnAxisDescription,
                                                          arrayView1d< real64 const > const values,
                                                          bool columnMajorValues,
                                                          string_view valueDescription )
{
  string const rowFmt = GEOS_FMT( "{} = {{}}", rowAxisDescription );
  string const columnFmt = GEOS_FMT( "{} = {{}}", columnAxisDescription );
  collectTableValues( coordX, coordY, values, columnMajorValues );
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

  for( auto errorIt = m_errors->begin(); errorIt != m_errors->end(); errorIt++ )
  {
    tableData1D.tableData.getErrorsList().addError( *errorIt );
  }

  // insert row value and row cell values
  for( auto const & [rowValue, rowMap] : m_data )
  {
    stdVector< TableData::CellData > currentRowValues;
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
