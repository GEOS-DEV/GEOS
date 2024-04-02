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

#include "common/TableFormatter.hpp"
namespace geos
{

TableFormatter::TableFormatter( TableLayout const & tableLayout ):
  m_tableLayout( tableLayout )
{}

void TableFormatter::fillTableColumnsFromRows( std::vector< TableLayout::Column > & columns,
                                               std::vector< std::vector< string > > const & rows ) const
{
  //TODO : reserve
  for( size_t idxRow = 0; idxRow < rows.size(); idxRow++ )
  {
    //TODO : reserve
    //TODO : if rows[idxRow].size()!=columns.size() ERROR/THROW/WARNING/ignore
    for( size_t idxColumn = 0; idxColumn < columns.size(); idxColumn++ )
    {
      columns[idxColumn].columnValues.push_back( rows[idxRow][idxColumn] );
    }
  }
}

///////////////////////////////////////////////////////////////////////
////// CSV Formatter implementation
///////////////////////////////////////////////////////////////////////

TableCSVFormatter::TableCSVFormatter( TableLayout const & tableLayout ):
  TableFormatter( tableLayout )
{
  m_tableLayout = tableLayout;
}

string TableCSVFormatter::headerToString() const
{
  std::stringstream oss;
  static constexpr string_view separator = ",";

  for( std::size_t idxColumn = 0; idxColumn < m_tableLayout.getColumns().size(); ++idxColumn )
  {
    oss << m_tableLayout.getColumns()[idxColumn].parameter.columnName;
    if( idxColumn < m_tableLayout.getColumns().size() - 1 )
    {
      oss << separator;
    }
  }
  oss << "\n";
  return oss.str();
}

string TableCSVFormatter::dataToString( TableData const & tableData ) const
{
  std::vector< std::vector< string > > const rowsValues = tableData.getTableDataRows();
  std::ostringstream oss;

  for( const auto & row : rowsValues )
  {
    //TODO : if row.size()!=tableLayout.columnsCount ERROR/THROW/WARNING/ignore
    for( size_t idxColumn = 0; idxColumn < row.size(); ++idxColumn )
    {
      oss << row[idxColumn];
      if( idxColumn < row.size() - 1 )
      {
        oss << ",";
      }
    }
    oss << "\n";
  }
  return oss.str();
}

string TableCSVFormatter::toString( TableData const & tableData ) const
{
  return headerToString() + dataToString( tableData );
}

///////////////////////////////////////////////////////////////////////
////// Log Formatter implementation
///////////////////////////////////////////////////////////////////////

/**
 * @brief Build a value cell given an alignment and spaces from "|"
 *
 * @param alignment The aligment of cell value
 * @param value The cell value
 * @param spaces The number of spaces in the cell
 * @return A formated cell
 */
string buildValueCell( TableLayout::Alignment const alignment, string_view value, integer const spaces )
{
  switch( alignment )
  {
    case TableLayout::right:   return GEOS_FMT( "{:>{}}", value, spaces );
    case TableLayout::left:    return GEOS_FMT( "{:<{}}", value, spaces );
    case TableLayout::center:  return GEOS_FMT( "{:^{}}", value, spaces );
    default:                   return GEOS_FMT( "{:>{}}", value, spaces );
  }
}

TableTextFormatter::TableTextFormatter( TableLayout const & tableLayout ):
  TableFormatter( tableLayout )
{}

string TableTextFormatter::toString( TableData const & tableData ) const
{
  std::ostringstream tableOutput;
  string sectionSeparator;
  std::vector< TableLayout::Column > columns = m_tableLayout.getColumns();
  integer const nbRows = tableData.getTableDataRows().size();

  fillTableColumnsFromRows( columns, tableData.getTableDataRows() );

  layoutToString( tableOutput, columns, nbRows, sectionSeparator );
  buildSectionRows( columns, sectionSeparator, tableOutput, nbRows, TableLayout::Section::values );
  tableOutput << '\n';

  return tableOutput.str();
}

string TableTextFormatter::layoutToString() const
{
  std::ostringstream tableOutput;
  std::vector< TableLayout::Column > columns = m_tableLayout.getColumns();
  string sectionSeparator;

  layoutToString( tableOutput, columns, 0, sectionSeparator );

  return tableOutput.str();
}

void TableTextFormatter::layoutToString( std::ostringstream & tableOutput,
                                         std::vector< TableLayout::Column > & columns,
                                         integer const nbRow,
                                         string & sectionSeparator ) const
{
  string titleRows;
  string topSeparator;
  size_t largestHeaderVectorSize = 0;
  std::vector< std::vector< string > > splitHeader;
  string const tableTitle = string( m_tableLayout.getTitle());

  parseAndStoreHeaderSections( columns, largestHeaderVectorSize, splitHeader );
  adjustHeaderSizesAndStore( columns, largestHeaderVectorSize, splitHeader );

  findAndSetMaxStringSize( columns, nbRow );
  computeAndBuildSeparator( columns, topSeparator, sectionSeparator );

  if( !tableTitle.empty())
  {
    buildTitleRow( titleRows, topSeparator, sectionSeparator );
    tableOutput << titleRows;
  }

  tableOutput << sectionSeparator + '\n';
  buildSectionRows( columns, sectionSeparator, tableOutput, largestHeaderVectorSize, TableLayout::Section::header );
}

void TableTextFormatter::parseAndStoreHeaderSections( std::vector< TableLayout::Column > const & columns,
                                                      size_t & largestHeaderVectorSize,
                                                      std::vector< std::vector< string > > & splitHeader ) const
{
  for( auto const & column : columns )
  {
    std::vector< string > splitHeaderParts;
    std::istringstream ss( column.parameter.columnName );
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
                                                    size_t const & largestHeaderVectorSize,
                                                    std::vector< std::vector< string > > & splitHeader ) const
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
                                                  size_t const & nbRows ) const
{
  string maxStringSize = "";
  for( auto & column : columns )
  {
    auto it = std::max_element( column.parameter.splitColumnName.begin(),
                                column.parameter.splitColumnName.end(),
                                []( const auto & a, const auto & b ) {
      return a.size() < b.size();
    } );

    maxStringSize = *it;
    for( size_t idxRow = 0; idxRow <  nbRows; idxRow++ )
    {
      string cell = column.columnValues[idxRow];

      if( maxStringSize.length() < cell.length())
      {
        maxStringSize = cell;
      }
    }

    column.m_maxStringSize = maxStringSize;
  }
}

void TableTextFormatter::computeAndSetMaxStringSize( std::vector< TableLayout::Column > & columns,
                                                     string::size_type const sectionlineLength,
                                                     string::size_type const titleLineLength ) const
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
                                                   string & sectionSeparator ) const
{
  integer const columnMargin = m_tableLayout.getColumnMargin();
  integer const borderMargin = m_tableLayout.getBorderMargin();
  integer const marginTitle = m_tableLayout.getMarginTitle();
  string tableTitle = string( m_tableLayout.getTitle() );

  string::size_type sectionlineLength = 0;
  string::size_type const titleLineLength = tableTitle.length() + ( marginTitle * 2 );
  integer const nbSpaceBetweenColumn = ( ( columns.size() - 1 ) *  columnMargin ) + (borderMargin * 2);

  if( !tableTitle.empty())
  {
    tableTitle = GEOS_FMT( "{:^{}}", tableTitle, titleLineLength );
  }

  for( auto const & column : columns )
  {
    sectionlineLength += column.m_maxStringSize.length();
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
  topSeparator = GEOS_FMT( "+{:-<{}}+", "", sectionSeparator.size() - 2 );  // -2 for ++
}

void TableTextFormatter::buildTitleRow( string & titleRows,
                                        string_view topSeparator,
                                        string_view sectionSeparator ) const
{
  titleRows = GEOS_FMT( "\n{}\n|", topSeparator );
  titleRows +=  buildValueCell( TableLayout::Alignment::center,
                                m_tableLayout.getTitle(),
                                (sectionSeparator.length() - 2)   // -2 for ||
                                );
  titleRows += GEOS_FMT( "{}\n", "|" );
}

void TableTextFormatter::buildSectionRows( std::vector< TableLayout::Column > & columns,
                                           string_view sectionSeparator,
                                           std::ostringstream & tableRows,
                                           integer const nbRows,
                                           TableLayout::Section const section ) const
{
  integer const columnMargin = m_tableLayout.getColumnMargin();
  integer const borderMargin = m_tableLayout.getBorderMargin();

  for( integer idxRow = 0; idxRow< nbRows; idxRow++ )
  {
    tableRows << GEOS_FMT( "{:<{}}", "|", 1 +  borderMargin );
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
      tableRows << buildValueCell( columns[idxColumn].parameter.alignment,
                                   cell,
                                   cellSize );

      if( idxColumn < columns.size() - 1 )
      {
        tableRows << GEOS_FMT( "{:^{}}", "|", columnMargin );
      }

    }

    if( columns.size() == 1 )
    {
      tableRows <<  GEOS_FMT( "{:>{}}\n", "|", columnMargin );
    }
    else
    {
      tableRows << GEOS_FMT( "{:>{}}\n", "|", borderMargin + 1 );
    }

  }
  if( nbRows != 0 )
  {
    tableRows << GEOS_FMT( "{}\n", sectionSeparator );
  }
}


}
