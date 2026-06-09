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
  // Compare row by row
  for( size_t i = 0; i < std::min( m_rows.size(), other.m_rows.size() ); ++i )
  {
    // Compare cells in current row
    for( size_t j = 0; j < std::min( m_rows[i].size(), other.m_rows[i].size() ); ++j )
    {
      if( m_rows[i][j].value < other.m_rows[i][j].value )
        return true;
      if( m_rows[i][j].value > other.m_rows[i][j].value )
        return false;
    }

    // If all compared cells are equal, the shorter row is considered less
    if( m_rows[i].size() != other.m_rows[i].size() )
      return m_rows[i].size() < other.m_rows[i].size();
  }

  // If all compared rows are equal, the table with fewer rows is considered less
  return m_rows.size() < other.m_rows.size();
}

bool TableData::operator==( TableData const & other ) const
{
  if( m_rows.size() != other.m_rows.size() )
    return false;

  for( size_t i = 0; i < m_rows.size(); ++i )
  {
    if( m_rows[i].size() != other.m_rows[i].size() )
      return false;

    for( size_t j = 0; j < m_rows[i].size(); ++j )
    {
      if( m_rows[i][j].value != other.m_rows[i][j].value )
        return false;
    }
  }

  return true;
}

void TableData::CellData::serialize( stdVector< buffer_unit_type > & out ) const
{
  basicSerialization::serializePrimitive( type, out );
  basicSerialization::serializeString( value, out );
}


size_t TableData::CellData::getSerializedSize() const
{
  return basicSerialization::sizeOfPrimitive( type ) + basicSerialization::sizeOfString( value );
}

size_t TableData::getSerializedSize() const
{
  size_t totalSize  =0;

  if( m_rows.empty())
    return totalSize;

  for( auto & row : m_rows )
  {
    size_t rowSize = 0;
    for( auto & cell : row )
    {
      rowSize += cell.getSerializedSize();
    }
    totalSize += sizeof(size_t) + rowSize;
  }
  return totalSize;
}

void TableData::serialize( stdVector< buffer_unit_type > & serializedTableData ) const
{
  if( m_rows.empty())
    return;

  for( auto & row : m_rows )
  {
    { // pack row size;
      size_t rowSize = 0;
      for( auto const & cell : row )
        rowSize += cell.getSerializedSize();
      basicSerialization::serializePrimitive( rowSize, serializedTableData );
    }

    { // pack cells
      for( auto const & cell : row )
      {
        cell.serialize( serializedTableData );
      }
    }
  }
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
  wellFormedValues = values;
  if( values.size() < nRow * nCol )
  {
    m_errors->addError( GEOS_FMT( "Warning: Not enough for the number of columns & rows:\n"
                                  "  - Expected {} values ({} columns x {} rows),\n  - Found {} values",
                                  nRow * nCol, nCol, nRow, values.size() ) );
    wellFormedValues.resizeDefault( nRow * nCol, 0 );
  }
  else if( values.size() > nRow * nCol )
  {
    m_errors->addError( GEOS_FMT( "Warning: Too much data for the number of columns & rows:\n"
                                  "  - Expected {} values ({} columns x {} rows),\n  - Found {} values."
                                  " Data may be misaligned",
                                  nRow * nCol, nCol, nRow, values.size() ) );
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

  for( auto const & error : *m_errors )
  {
    tableData1D.tableData.getErrorsList().addError( error );
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

void basicSerialization::serializeString ( string const & data, stdVector< buffer_unit_type > & out )
{
  basicSerialization::serializePrimitive( data.size(), out );
  auto * begin = data.data();
  auto * end = begin + data.size();
  out.insert( out.end(), begin, end );
}

void basicSerialization::deserializeString( string & str, buffer_unit_type const * & ptr, buffer_unit_type const * end )
{
  string::size_type strSize = 0;
  basicSerialization::deserializePrimitive( strSize, ptr, end );
  if( static_cast< long >(strSize) > std::distance( ptr, end ) )
  {
    throw std::runtime_error( "buffer overflow reading string" );
  }
  str.assign( ptr, ptr + strSize );
  ptr += str.size();
}

bool tableDataSorting::positiveNumberStringComp( string_view s1, string_view s2 )
{
  auto split = []( string_view s, string & intPart, string & decPart )
  {
    size_t dotPos = s.find( '.' );
    if( dotPos == string::npos )
    {
      intPart = s;
      decPart = "";
    }
    else
    {
      intPart = s.substr( 0, dotPos );
      decPart = s.substr( dotPos + 1 );
    }
  };

  string s1Int, s1Dec, s2Int, s2Dec;
  split( s1, s1Int, s1Dec );
  split( s2, s2Int, s2Dec );

  if( s1Int.length() != s2Int.length())
    return s1Int.length() < s2Int.length();

  if( s1Int != s2Int )
    return s1Int < s2Int;

  size_t minLen = std::min( s1Dec.length(), s2Dec.length());
  for( size_t i = 0; i < minLen; ++i )
  {
    if( s1Dec[i] != s2Dec[i] )
      return s1Dec[i] < s2Dec[i];
  }


  return false;
}

}
