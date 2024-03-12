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
 * @file TableFormatter.cpp
 */

#include "codingUtilities/TableFormatter.hpp"
namespace geos
{

TableFormatter::TableFormatter( TableLayout tableLayout )
{
  m_tableLayout = tableLayout;
}

void TableFormatter::fillTableColumnsFromRows( std::vector< TableLayout::Column > & columns,
                                               std::vector< std::vector< string > > const & rows )
{
  for( size_t idxRow = 0; idxRow < rows.size(); idxRow++ )
  {
    for( size_t idxColumn = 0; idxColumn < columns.size(); idxColumn++ )
    {
      columns[idxColumn].columnValues.push_back( rows[idxRow][idxColumn] );
    }
  }
}

///////////////////////////////////////////////////////////////////////
////// CSV Formatter implementation
///////////////////////////////////////////////////////////////////////

TableCSVFormatter::TableCSVFormatter( TableLayout tableLayout )
{
  m_tableLayout = tableLayout;
}

string TableCSVFormatter::headerToString()
{
  string headerValues = "";
  string separator = ",";

  for( std::size_t idxColumn = 0; idxColumn < m_tableLayout.m_columns.size(); ++idxColumn )
  {
    headerValues += m_tableLayout.m_columns[idxColumn].parameter.columnName;
    if( idxColumn < m_tableLayout.m_columns.size() - 1 )
    {
      headerValues += separator;
    }
  }
  headerValues += "\n";
  return headerValues;
}

string TableCSVFormatter::dataToString( TableData2D & tableData )
{

  string data;

  std::vector< std::vector< string > > rowsValues = tableData.getTableDataRows();

  tableData.buildRows( rowsValues );

  for( std::vector< string > const & row : rowsValues )
  {
    for( size_t idxRow = 0; idxRow < row.size(); idxRow++ )
    {
      data += row[idxRow];
      if( idxRow < row.size() - 1 )
      {
        data += ",";
      }
    }
    data += "\n";
  }
  return data;
}

string TableCSVFormatter::dataToString( TableData & tableData )
{
  std::vector< std::vector< string > > rowsValues = tableData.getTableDataRows();
  std::vector< TableLayout::Column > columns = m_tableLayout.m_columns;
  string data;

  fillTableColumnsFromRows( columns, rowsValues );

  for( std::vector< string > const & row : rowsValues )
  {
    for( size_t idxRow = 0; idxRow < row.size(); idxRow++ )
    {
      data += row[idxRow];
      if( idxRow < row.size() - 1 )
      {
        data += ",";
      }
    }
    data += "\n";
  }
  return data;
}

///////////////////////////////////////////////////////////////////////
////// Log Formatter implementation
///////////////////////////////////////////////////////////////////////

/**
 * @brief Build a value cell given an alignment and spaces from "|"
 *
 * @param alignment
 * @param value
 * @param spaces
 * @return A cell value
 */
string buildValueCell( TableLayout::Alignment const alignment, string_view value, integer const spaces )
{
  switch( alignment )
  {
    case TableLayout::right:   return GEOS_FMT( "{:>{}}", value, spaces );
    case TableLayout::left:    return GEOS_FMT( "{:<{}}", value, spaces );
    case TableLayout::middle:  return GEOS_FMT( "{:^{}}", value, spaces );
    default:             return GEOS_FMT( "{:<{}}", value, spaces );
  }
}

TableTextFormatter::TableTextFormatter( TableLayout tableLayout ):
  TableFormatter( tableLayout )
{}

string TableTextFormatter::ToString( TableData & tableData )
{
  return constructAndReturnTable( tableData.getTableDataRows() );
}

string TableTextFormatter::ToString( TableData2D & tableData )
{
  std::vector< std::vector< string > > tableRows = tableData.getTableDataRows();
  tableData.buildRows( tableRows );
  return constructAndReturnTable( tableRows );
}

string TableTextFormatter::layoutToString()
{
  string titleRows;
  string topSeparator;
  string sectionSeparator;

  size_t largestHeaderVectorSize = 0;
  std::vector< std::vector< string > > splitHeader;
  std::vector< TableLayout::Column > columns = m_tableLayout.m_columns;

  string tableTitle = string( m_tableLayout.getTitle());

  parseAndStoreHeaderSections( columns, largestHeaderVectorSize, splitHeader );
  adjustHeaderSizesAndStore( columns, largestHeaderVectorSize, splitHeader );

  findAndSetMaxStringSize( columns, 0 );
  computeAndBuildSeparator( columns, topSeparator, sectionSeparator );

  if( !tableTitle.empty())
  {
    buildTitleRow( titleRows, topSeparator, sectionSeparator );
  }

  string tableRows = GEOS_FMT( "{}\n", sectionSeparator );
  buildSectionRows( columns, sectionSeparator, tableRows, largestHeaderVectorSize, TableLayout::Section::header );

  string tableOutput = titleRows + tableRows + '\n';

  return tableOutput;
}

void TableTextFormatter::parseAndStoreHeaderSections( std::vector< TableLayout::Column > & columns,
                                                      size_t & largestHeaderVectorSize,
                                                      std::vector< std::vector< string > > & splitHeader )
{
  for( size_t columnParamIdx = 0; columnParamIdx< columns.size(); columnParamIdx++ )
  {
    std::vector< string > splitHeaderParts;
    std::istringstream ss( columns[columnParamIdx].parameter.columnName );
    string subColumnNames;

    while( getline( ss, subColumnNames, '\n' ))
    {
      splitHeaderParts.push_back( subColumnNames );
    }

    size_t const cellSize = splitHeaderParts.size();
    largestHeaderVectorSize = std::max( largestHeaderVectorSize, cellSize );

    splitHeader.push_back( splitHeaderParts );
  }
}

void TableTextFormatter::adjustHeaderSizesAndStore( std::vector< TableLayout::Column > & columns,
                                                    size_t largestHeaderVectorSize,
                                                    std::vector< std::vector< string > > & splitHeader )
{
  for( size_t columnParamIdx = 0; columnParamIdx < columns.size(); columnParamIdx++ )
  {
    if( splitHeader[columnParamIdx].size() < largestHeaderVectorSize )
    {
      integer const whiteRowToAdd = largestHeaderVectorSize - splitHeader[columnParamIdx].size();
      splitHeader[columnParamIdx].insert( splitHeader[columnParamIdx].end(), whiteRowToAdd, " " );
    }
    columns[columnParamIdx].parameter.splitColumnName = splitHeader[columnParamIdx];
  }
}

void TableTextFormatter::findAndSetMaxStringSize( std::vector< TableLayout::Column > & columns,
                                                  size_t nbRows )
{
  string maxStringSize = "";
  for( size_t idxColumn  = 0; idxColumn <  columns.size(); idxColumn++ )
  {
    auto it = std::max_element( columns[idxColumn].parameter.splitColumnName.begin(),
                                columns[idxColumn].parameter.splitColumnName.end(),
                                []( const auto & a, const auto & b ) {
      return a.size() < b.size();
    } );

    maxStringSize = *it;
    for( size_t idxRow = 0; idxRow <  nbRows; idxRow++ )
    {
      string cell = columns[idxColumn].columnValues[idxRow];

      if( maxStringSize.length() < cell.length())
      {
        maxStringSize = cell;
      }
    }

    columns[idxColumn].m_maxStringSize = maxStringSize;
  }
}

void TableTextFormatter::computeAndSetMaxStringSize( std::vector< TableLayout::Column > & columns,
                                                     string::size_type sectionlineLength,
                                                     string::size_type titleLineLength )
{
  integer extraLinesPerColumn;
  integer extraLines;
  integer newStringSize;

  extraLines = titleLineLength - sectionlineLength;
  extraLinesPerColumn = std::ceil( extraLines / columns.size() );

  for( std::size_t idxColumn = 0; idxColumn < columns.size(); ++idxColumn )
  {
    newStringSize = extraLinesPerColumn + columns[idxColumn].m_maxStringSize.size();
    if( idxColumn == columns.size() - 1 ||  columns.size() == 1 )
    {
      columns[idxColumn].m_maxStringSize = GEOS_FMT( "{:>{}}",
                                                     columns[idxColumn].m_maxStringSize,
                                                     newStringSize + m_tableLayout.getColumnMargin() );
    }
    else
    {
      columns[idxColumn].m_maxStringSize = GEOS_FMT( "{:>{}}",
                                                     columns[idxColumn].m_maxStringSize,
                                                     newStringSize );
    }
  }
}

void TableTextFormatter::computeAndBuildSeparator( std::vector< TableLayout::Column > & columns,
                                                   string & topSeparator,
                                                   string & sectionSeparator )
{
  integer columnMargin = m_tableLayout.getColumnMargin();
  integer borderMargin = m_tableLayout.getBorderMargin();
  integer marginTitle = m_tableLayout.getMarginTitle();
  string tableTitle = string( m_tableLayout.getTitle() );

  string::size_type sectionlineLength = 0;
  string::size_type titleLineLength = tableTitle.length() + ( marginTitle * 2 );
  integer nbSpaceBetweenColumn = ( ( columns.size() - 1 ) *  columnMargin ) + (borderMargin * 2);

  if( !tableTitle.empty())
  {
    tableTitle = GEOS_FMT( "{:^{}}", tableTitle, titleLineLength );
  }

  for( std::size_t i = 0; i < columns.size(); ++i )
  {
    sectionlineLength += columns[i].m_maxStringSize.length();
  }

  sectionlineLength += nbSpaceBetweenColumn;
  if( sectionlineLength < titleLineLength )
  {
    computeAndSetMaxStringSize( columns, sectionlineLength, titleLineLength );
  }
  if( columns.size() == 1 )
  {
    sectionSeparator +=  GEOS_FMT( "+{:-<{}}+",
                                   "",
                                   ( columns[0].m_maxStringSize.length() + (borderMargin - 1) + columnMargin ));
  }
  else
  {
    for( std::size_t idxColumn = 0; idxColumn < columns.size(); ++idxColumn )
    {
      integer const cellSize = columns[idxColumn].m_maxStringSize.length();
      if( idxColumn == 0 )
      {
        sectionSeparator +=  GEOS_FMT( "+{:-<{}}", "", ( cellSize + borderMargin ));
      }
      else if( idxColumn == (columns.size() - 1))
      {
        sectionSeparator += GEOS_FMT( "{:-^{}}", "+", columnMargin );
        sectionSeparator += GEOS_FMT( "{:->{}}", "+", ( cellSize + borderMargin + 1 ) );
      }
      else
      {
        sectionSeparator += GEOS_FMT( "{:-^{}}", "+", columnMargin );
        sectionSeparator += GEOS_FMT( "{:->{}}", "", cellSize );
      }
    }
  }
  topSeparator = GEOS_FMT( "+{:-<{}}+", "", sectionSeparator.size() - 2 );// -2 for ++
}

void TableTextFormatter::buildTitleRow( string & titleRows, string_view topSeparator, string_view sectionSeparator )
{
  titleRows = GEOS_FMT( "\n{}\n|", topSeparator );
  titleRows +=  buildValueCell( TableLayout::Alignment::middle,
                                m_tableLayout.getTitle(),
                                (sectionSeparator.length() - 2) // -2 for ||
                                );
  titleRows += GEOS_FMT( "{}\n", "|" );
}

void TableTextFormatter::buildSectionRows( std::vector< TableLayout::Column > & columns,
                                           string_view sectionSeparator,
                                           string & tableRows,
                                           integer const nbRows,
                                           TableLayout::Section const section )
{
  integer columnMargin = m_tableLayout.getColumnMargin();
  integer borderMargin = m_tableLayout.getBorderMargin();

  for( integer idxRow = 0; idxRow< nbRows; idxRow++ )
  {
    tableRows += GEOS_FMT( "{:<{}}", "|", 1 +  borderMargin );
    for( std::size_t idxColumn = 0; idxColumn < columns.size(); ++idxColumn )
    {
      string cell;

      if( section == TableLayout::Section::header )
      {
        cell = columns[idxColumn].parameter.splitColumnName[idxRow];
      }
      else
      {
        cell = columns[idxColumn].columnValues[idxRow];
      }
      integer const cellSize = columns[idxColumn].m_maxStringSize.length();
      tableRows += buildValueCell( columns[idxColumn].parameter.alignment,
                                   cell,
                                   cellSize );

      if( idxColumn < columns.size() - 1 )
      {
        tableRows += GEOS_FMT( "{:^{}}", "|", columnMargin );
      }

    }
    if( columns.size() == 1 )
    {
      tableRows +=  GEOS_FMT( "{:>{}}\n", "|", columnMargin );
    }
    else
    {
      tableRows += GEOS_FMT( "{:>{}}\n", "|", borderMargin + 1 );
    }

  }
  if( nbRows != 0 )
  {
    tableRows += GEOS_FMT( "{}\n", sectionSeparator );
  }
}

string TableTextFormatter::constructAndReturnTable( std::vector< std::vector< string > > & rowsValues )
{
  string titleRows;
  string topSeparator;
  string sectionSeparator;

  size_t largestHeaderVectorSize = 0;
  std::vector< std::vector< string > > splitHeader;
  std::vector< TableLayout::Column > columns = m_tableLayout.m_columns;

  string tableTitle = string( m_tableLayout.getTitle());

  fillTableColumnsFromRows( columns, rowsValues );

  parseAndStoreHeaderSections( columns, largestHeaderVectorSize, splitHeader );
  adjustHeaderSizesAndStore( columns, largestHeaderVectorSize, splitHeader );

  findAndSetMaxStringSize( columns, rowsValues.size());
  computeAndBuildSeparator( columns, topSeparator, sectionSeparator );

  if( !tableTitle.empty())
  {
    buildTitleRow( titleRows, topSeparator, sectionSeparator );
  }

  string tableRows = GEOS_FMT( "{}\n", sectionSeparator );
  buildSectionRows( columns, sectionSeparator, tableRows, largestHeaderVectorSize, TableLayout::Section::header );
  buildSectionRows( columns, sectionSeparator, tableRows, rowsValues.size(), TableLayout::Section::values );

  string tableOutput = titleRows + tableRows + '\n';

  return tableOutput;
}

}
