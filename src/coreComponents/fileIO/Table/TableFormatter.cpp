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

#include <numeric>
#include "TableFormatter.hpp"
namespace geos
{

TableFormatter::TableFormatter( TableLayout const & tableLayout ):
  m_tableLayout( tableLayout )
{}

void TableFormatter::fillTableColumnsFromRows( std::vector< TableLayout::Column > & columns,
                                               std::vector< std::vector< string > > const & rows,
                                               std::vector< string > & msgTableError ) const
{
  bool isConsistent = true;
  for( size_t idxRow = 0; idxRow < rows.size(); idxRow++ )
  {
    for( size_t idxColumn = 0; idxColumn < columns.size(); idxColumn++ )
    {
      // Case of a hidden column during initialization
      if( m_tableLayout.getColumns()[idxColumn].m_parameter.enabled )
      {
        columns[idxColumn].m_columnValues.push_back( rows[idxRow][idxColumn] );
      }
    }

    if( rows[idxRow].size()!=columns.size() )
    {
      isConsistent = false;
    }
  }

  if( !isConsistent )
  {
    msgTableError.push_back( "The number of columns displayed on the table does not match to the columns that have been initialized in TableLayout" );
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
    oss << m_tableLayout.getColumns()[idxColumn].m_parameter.columnName;
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
 * @brief Build  cell given an alignment, a value and spaces
 * @param alignment The aligment of cell value
 * @param value The cell value
 * @param spaces The number of spaces in the cell
 * @return A formated cell
 */
string buildCell( TableLayout::Alignment const alignment, string_view value, integer const spaces )
{
  switch( alignment )
  {
    case TableLayout::right:   return GEOS_FMT( "{:>{}}", value, spaces );
    case TableLayout::left:    return GEOS_FMT( "{:<{}}", value, spaces );
    case TableLayout::center:  return GEOS_FMT( "{:^{}}", value, spaces );
    default:                   return GEOS_FMT( "{:>{}}", value, spaces );
  }
}

/**
 * @brief Detect columns who are not displayed from TableLayout and therefore modify columns / tableDataRows vectors
 * @param columns Vector built in TableLayout containing all columns with their parameters
 * @param tableDataRows Vector built in TableData containing all rows values
 */
void formatColumnsFromLayout( std::vector< TableLayout::Column > & columns, std::vector< std::vector< string > > & tableDataRows )
{
  integer idxColumn = 0;
  for( auto & column : columns )
  {
    if( !column.m_parameter.enabled )
    {
      columns.erase( columns.begin() + idxColumn );
      for( auto & row : tableDataRows )
      {
        row.erase( row.begin() + idxColumn );
      }
    }
    ++idxColumn;
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
  std::vector< std::vector< string > > tableDataRows = tableData.getTableDataRows();
  std::vector< string > msgTableError = tableData.getErrorMsgConversion();
  integer const nbRows = tableDataRows.size();

  formatColumnsFromLayout( columns, tableDataRows );
  fillTableColumnsFromRows( columns, tableDataRows, msgTableError );

  tableOutput << '\n';
  layoutToString( tableOutput, columns, msgTableError, sectionSeparator );
  buildSectionRows( columns, sectionSeparator, tableOutput, nbRows, TableLayout::Section::values );
  tableOutput << '\n';

  return tableOutput.str();
}

void TableTextFormatter::layoutToString( std::ostringstream & tableOutput,
                                         std::vector< TableLayout::Column > & columns,
                                         std::vector< string > & msgTableError,
                                         string & sectionSeparator ) const
{
  string topSeparator;
  size_t largestHeaderVectorSize = 0;
  std::vector< std::vector< string > > splitHeader;
  string const tableTitle = string( m_tableLayout.getTitle());

  parseAndStoreHeaderSections( columns, largestHeaderVectorSize, splitHeader );
  adjustHeaderSizesAndStore( columns, largestHeaderVectorSize, splitHeader );

  findAndSetMaxStringSize( columns );
  computeTableMaxLineLength( columns, msgTableError );
  buildTableSeparators( columns, topSeparator, sectionSeparator );

  addTopRow( tableOutput, msgTableError, topSeparator, sectionSeparator );
  addTopRow( tableOutput, tableTitle, topSeparator, sectionSeparator );

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
    std::istringstream ss( column.m_parameter.columnName );
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
  for( size_t columnParamIdx = 0; columnParamIdx < columns.size(); ++columnParamIdx )
  {
    if( splitHeader[columnParamIdx].size() < largestHeaderVectorSize )
    {
      integer const whiteRowToAdd = largestHeaderVectorSize - splitHeader[columnParamIdx].size();
      splitHeader[columnParamIdx].insert( splitHeader[columnParamIdx].end(), whiteRowToAdd, " " );
    }
    columns[columnParamIdx].m_parameter.splitColumnName = splitHeader[columnParamIdx];
  }
}

void TableTextFormatter::findAndSetMaxStringSize( std::vector< TableLayout::Column > & columns ) const
{
  string maxStringSize;
  for( auto & column : columns )
  {
    auto it = std::max_element( column.m_parameter.splitColumnName.begin(),
                                column.m_parameter.splitColumnName.end(),
                                []( const auto & a, const auto & b ) {
      return a.size() < b.size();
    } );
    maxStringSize = *it;

    for( size_t idxRow = 0; idxRow <  column.m_columnValues.size(); ++idxRow )
    {
      string cell = column.m_columnValues[idxRow];

      if( maxStringSize.length() < cell.length())
      {
        maxStringSize = cell;
      }
    }

    column.m_maxStringSize = maxStringSize;
  }
}

void TableTextFormatter::recalculateMaxStringSize( std::vector< TableLayout::Column > & columns,
                                                   integer const extraLines ) const
{
  integer const extraLinesPerColumn = std::ceil( extraLines / columns.size() );

  for( std::size_t idxColumn = 0; idxColumn < columns.size(); ++idxColumn )
  {
    integer newStringSize = extraLinesPerColumn + columns[idxColumn].m_maxStringSize.size();
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

void TableTextFormatter::computeTableMaxLineLength( std::vector< TableLayout::Column > & columns,
                                                    std::vector< string > & msgTableError ) const
{
  integer const columnMargin = m_tableLayout.getColumnMargin();
  integer const borderMargin = m_tableLayout.getBorderMargin();
  integer const marginTitle = m_tableLayout.getMarginTitle();
  string const tableTitle = string( m_tableLayout.getTitle() );

  string::size_type sectionlineLength = 0;
  string::size_type maxTopLineLength =  tableTitle.length() + ( marginTitle * 2 );
  string::size_type msgTableErrorLength = marginTitle * 2;
  integer const nbSpaceBetweenColumn = ( ( columns.size() - 1 ) *  columnMargin ) + (borderMargin * 2);

  if( msgTableError.size() != 0 )
  {
    auto it = std::max_element( msgTableError.begin(), msgTableError.end(),
                                []( const auto & a, const auto & b ) {
      return a.size() < b.size();
    } );
    string maxStringSize = *it;

    msgTableErrorLength += maxStringSize.size() + 1; // for \n set later

    if( maxTopLineLength  < msgTableErrorLength )
    {
      maxTopLineLength = msgTableErrorLength;
    }
  }

  for( auto const & column : columns )
  {
    sectionlineLength += column.m_maxStringSize.length();
  }
  sectionlineLength += nbSpaceBetweenColumn;

  if( sectionlineLength < maxTopLineLength )
  {
    integer const extraLines = maxTopLineLength - sectionlineLength;
    recalculateMaxStringSize( columns, extraLines );
  }
}

void TableTextFormatter::buildTableSeparators( std::vector< TableLayout::Column > & columns,
                                               string & topSeparator,
                                               string & sectionSeparator ) const
{
  integer const columnMargin = m_tableLayout.getColumnMargin();
  integer const borderMargin = m_tableLayout.getBorderMargin();

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

void TableTextFormatter::addTopRow( std::ostringstream & tableOutput,
                                    string const & msg,
                                    string_view topSeparator,
                                    string_view sectionSeparator ) const
{
  if( !msg.empty())
  {
    buildTopRow( tableOutput, msg, topSeparator, sectionSeparator );
  }
}

void TableTextFormatter::addTopRow( std::ostringstream & tableOutput,
                                    std::vector< string > const & msg,
                                    string_view topSeparator,
                                    string_view sectionSeparator ) const
{
  if( msg.size() != 0 )
  {
    string strConcat;
    for( const std::string & str : msg )
    {
      if( !strConcat.empty())
      {
        strConcat += '\n';
      }
      strConcat += str;
    }
    buildTopRow( tableOutput, strConcat, topSeparator, sectionSeparator );
  }
}


void TableTextFormatter::buildTopRow( std::ostringstream & tableOutput,
                                      string const & msg,
                                      string_view topSeparator,
                                      string_view sectionSeparator ) const
{
  size_t nbLine = std::count_if( msg.begin(), msg.end(), []( char c ){return c =='\n';} ) + 1;
  std::vector< string > messages;
  std::istringstream ss( msg );
  string subMsg;

  while( getline( ss, subMsg, '\n' ))
  {
    messages.push_back( subMsg );
  }

  tableOutput << GEOS_FMT( "{}\n", topSeparator );
  for( size_t idxLine = 0; idxLine< nbLine; ++idxLine )
  {
    tableOutput << GEOS_FMT( "{}", "|" );
    tableOutput << buildCell( TableLayout::Alignment::center,
                              messages[idxLine],
                              (sectionSeparator.length() - 2)  // -2 for ||
                              );
    tableOutput << GEOS_FMT( "{}\n", "|" );
  }
}

void TableTextFormatter::buildSectionRows( std::vector< TableLayout::Column > const & columns,
                                           string_view sectionSeparator,
                                           std::ostringstream & tableRows,
                                           integer const nbRows,
                                           TableLayout::Section const section ) const
{
  integer const columnMargin = m_tableLayout.getColumnMargin();
  integer const borderMargin = m_tableLayout.getBorderMargin();

  for( integer idxRow = 0; idxRow< nbRows; ++idxRow )
  {
    tableRows << GEOS_FMT( "{:<{}}", "|", 1 +  borderMargin );
    for( std::size_t idxColumn = 0; idxColumn < columns.size(); ++idxColumn )
    {
      string cell;

      if( section == TableLayout::Section::header )
      {
        cell = columns[idxColumn].m_parameter.splitColumnName[idxRow];
      }
      else
      {
        cell = columns[idxColumn].m_columnValues[idxRow];
      }
      integer const cellSize = columns[idxColumn].m_maxStringSize.length();
      tableRows << buildCell( columns[idxColumn].m_parameter.alignment,
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
