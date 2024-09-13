/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */


//TODO REAJUSTER LAPPEL DES FONCTIONS
/**
 * @file TableFormatter.cpp
 */

#include "TableFormatter.hpp"
#include <numeric>
#include "common/format/StringUtilities.hpp"
#include "TableFormatter.hpp"

namespace geos
{

TableFormatter::TableFormatter( TableLayout const & tableLayout ):
  m_tableLayout( tableLayout )
{}

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
    oss << stringutilities::join( row.cbegin(), row.cend(), "," ) << "\n";
  }

  return oss.str();
}

template<>
string TableCSVFormatter::toString< TableData >( TableData const & tableData ) const
{
  return headerToString() + dataToString( tableData );
}

///////////////////////////////////////////////////////////////////////
////// Log Formatter implementation
///////////////////////////////////////////////////////////////////////

void transpose( std::vector< std::vector< std::string > > & dest, std::vector< std::vector< std::string > > const & source )
{
  for( size_t idxRow = 0; idxRow < source.size(); ++idxRow )
  {
    for( size_t idxCol = 0; idxCol < source[idxRow].size(); ++idxCol )
    {
      dest[idxCol][idxRow] = source[idxRow][idxCol];
    }
  }
}

/**
 * @brief Build  cell given an alignment, a value and spaces
 * @param alignment The aligment of cell value
 * @param value The cell value
 * @param spaces The number of spaces in the cell
 * @return A formated cell
 */
string buildCell( TableLayout::Alignment const alignment, string_view value, size_t const spaces )
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
void updateVisibleColumns( std::vector< TableLayout::Column > & columns,
                           std::vector< std::vector< string > > & tableDataRows )
{
  integer idxColumn = 0;
  for( auto iterColumn = columns.begin(); iterColumn!=columns.end(); )
  {
    if( !iterColumn->m_parameter.enabled )
    {
      iterColumn = columns.erase( iterColumn );
      for( auto & row : tableDataRows )
      {
        row.erase( row.begin() + idxColumn );
      }
    }
    else
    {
      ++iterColumn;
      ++idxColumn;
    }
  }
}

TableTextFormatter::TableTextFormatter( TableLayout const & tableLayout ):
  TableFormatter( tableLayout )
{}

string TableTextFormatter::layoutToString() const
{
  std::ostringstream tableOutput;
  std::vector< TableLayout::Column > columns = m_tableLayout.getColumns();
  string const tableTitle = string( m_tableLayout.getTitle());
  size_t nbHeaderRows = 0;
  TableData tableData;
  std::string sectionSeparatingLine;
  std::string topSeparator;

  prepareAndBuildTable( columns, tableData, nbHeaderRows, sectionSeparatingLine, topSeparator );

  tableOutput << '\n'; //TODO Fonction
  outputTitleRow( tableOutput, topSeparator );
  tableOutput << GEOS_FMT( "{}\n", sectionSeparatingLine );

  outputHeaderSectionRows( columns, tableOutput, nbHeaderRows, sectionSeparatingLine );

  return tableOutput.str();
}

template<>
string TableTextFormatter::toString< TableData >( TableData const & tableData ) const
{
  std::ostringstream tableOutput;
  std::vector< TableLayout::Column > columns         = m_tableLayout.getColumns();
  size_t nbHeaderRows = 0;
  std::string sectionSeparatingLine;
  std::string topSeparator;

  prepareAndBuildTable( columns, tableData, nbHeaderRows, sectionSeparatingLine, topSeparator );

  outputTable( tableOutput, columns, tableData, nbHeaderRows, sectionSeparatingLine, topSeparator );

  return tableOutput.str();
}

void TableTextFormatter::prepareAndBuildTable( std::vector< TableLayout::Column > & columns,
                                               TableData const & tableData,
                                               size_t & nbHeaderRows,
                                               std::string & sectionSeparatingLine,
                                               std::string & topSeparator ) const
{
  std::vector< std::vector< string > > tableDataRows = tableData.getTableDataRows();
  std::vector< std::vector< string > > splitHeaders;

  if( !tableData.getTableDataRows().empty())
  {
    updateVisibleColumns( columns, tableDataRows );
    populateColumnsFromTableData( columns, tableDataRows, false );
  }

  splitAndMergeColumnHeaders( columns, nbHeaderRows, splitHeaders );

  for( auto & column : columns )
  {
    findAndSetLongestColumnString( column, column.m_maxStringSize, 0 );
  }
  computeTableWidth( columns );
  buildTableSeparators( columns, sectionSeparatingLine, topSeparator );
}

void TableTextFormatter::outputTable( std::ostringstream & tableOutput,
                                      std::vector< TableLayout::Column > & columns,
                                      TableData const & tableData,
                                      size_t & nbHeaderRows,
                                      std::string const & sectionSeparatingLine,
                                      std::string const & topSeparator ) const
{
  integer const nbValuesRows = tableData.getTableDataRows().size();

  tableOutput << '\n'; //TODO Fonction
  outputTitleRow( tableOutput, topSeparator );
  tableOutput << GEOS_FMT( "{}\n", sectionSeparatingLine );

  outputHeaderSectionRows( columns, tableOutput, nbHeaderRows, sectionSeparatingLine );

  outputValuesSectionRows( columns, tableOutput, nbValuesRows, sectionSeparatingLine );

  tableOutput << '\n';
}

void TableTextFormatter::populateColumnsFromTableData( std::vector< TableLayout::Column > & columns,
                                                       std::vector< std::vector< std::string > > const & tableDataRows,
                                                       bool isSubColumn ) const
{
  size_t currentColumn = 0;
  std::vector< std::vector< std::string > > tColumnsValues( tableDataRows[0].size(), std::vector< std::string >( tableDataRows.size()));
  if( !isSubColumn )
  {
    transpose( tColumnsValues, tableDataRows );
  }

  for( auto & column : columns )
  {
    if( column.subColumn.empty())
    {
      column.m_columnValues = !isSubColumn ?
                              tColumnsValues[currentColumn++] : tableDataRows[currentColumn++];
    }
    else
    {
      std::vector< std::vector< std::string > > subColumnValues( tColumnsValues.begin() + currentColumn,
                                                                 tColumnsValues.begin() + currentColumn + column.m_parameter.subColumns.size());
      populateColumnsFromTableData( column.subColumn, subColumnValues, true );

      std::vector< std::vector< std::string > > tSubColumnValues( subColumnValues[0].size(), std::vector< std::string >( subColumnValues.size()) );
      transpose( tSubColumnValues, subColumnValues );
      for( const auto & columnValues : tSubColumnValues )
      { // add all subcolumn values  in parent column
        column.m_columnValues.insert( column.m_columnValues.end(), columnValues.begin(), columnValues.end() );
      }
      currentColumn += subColumnValues.size();
    }
  }
}

void TableTextFormatter::splitAndMergeColumnHeaders( std::vector< TableLayout::Column > & columns,
                                                     size_t & nbHeaderRows,
                                                     std::vector< std::vector< string > > & splitHeaders ) const
{

  splitHeaders.reserve( columns.size() );

  for( auto & column : columns )
  {
    std::vector< string > splitHeaderParts;
    std::istringstream ss( column.m_parameter.columnName );
    string subColumnNames;

    while( getline( ss, subColumnNames, '\n' ))
    {
      splitHeaderParts.push_back( subColumnNames );
    }

    splitHeaders.push_back( std::move( splitHeaderParts ) );

    if( !column.subColumn.empty())
    {
      std::vector< std::vector< string > > splitSubColHeaders;
      size_t nbHeaderSubColRows = 0;
      splitAndMergeColumnHeaders( column.subColumn, nbHeaderSubColRows, splitSubColHeaders );
    }
  }

  nbHeaderRows = std::max_element( splitHeaders.begin(), splitHeaders.end(),
                                   []( auto const & v1, auto const & v2 ) { return v1.size() < v2.size(); } )->size();

  for( auto & headerParts : splitHeaders )
  {
    if( headerParts.size() < nbHeaderRows )
    {
      headerParts.resize( nbHeaderRows, " " );
    }
    columns[&headerParts - &splitHeaders[0]].m_parameter.splitColumnNames = headerParts;
  }

}

void TableTextFormatter::findAndSetLongestColumnString( TableLayout::Column & column, std::vector< std::string > & maxStringSize, integer const idxMaxString ) const
{
  string maxStringColumn;
  { // header case
    auto const maxStringSizeHeader = *std::max_element( column.m_parameter.splitColumnNames.begin(),
                                                        column.m_parameter.splitColumnNames.end(),
                                                        []( const auto & a, const auto & b ) {return a.size() < b.size();} );
    maxStringColumn = maxStringSizeHeader;
    maxStringSize.push_back( maxStringSizeHeader );
  }

  {  // values case
    if( column.subColumn.empty() )
    {
      auto const maxStringSizeCell = *std::max_element( column.m_columnValues.begin(),
                                                        column.m_columnValues.end(),
                                                        []( const auto & a, const auto & b ) {return a.size() < b.size();} );
      if( maxStringColumn.length() < maxStringSizeCell.length())
      {
        maxStringColumn = maxStringSizeCell;
      }
    }
  }

  { // Update max string size if necessary
    if( maxStringSize[idxMaxString].length() < maxStringColumn.length() )
    {
      maxStringSize[idxMaxString] = maxStringColumn;
    }
  }

  { // subcolumn values case
    if( !column.subColumn.empty() )
    {
      column.m_maxStringSize.clear();
      for( size_t idxSubColumn = 0; idxSubColumn < column.subColumn.size(); ++idxSubColumn )
      {
        findAndSetLongestColumnString( column.subColumn[idxSubColumn], column.m_maxStringSize, idxSubColumn );
      }
    }
  }

  if( column.m_maxStringSize.empty() )
  {
    column.m_maxStringSize.push_back( maxStringColumn );
  }
}

void TableTextFormatter::computeTableWidth( std::vector< TableLayout::Column > & columns ) const
{
  integer const columnMargin = m_tableLayout.getColumnMargin();
  integer const borderMargin = m_tableLayout.getBorderMargin();
  string const tableTitle = string( m_tableLayout.getTitle() );

  string::size_type const sectionLengthWithSpacing = ( ( columns.size() - 1 ) *  columnMargin ) + (borderMargin * 2);
  string::size_type sectionlineLength = sectionLengthWithSpacing;
  string const spaces =  std::string( m_tableLayout.getColumnMargin(), ' ' );

  { // Compute total length of all columns with margins
    sectionlineLength += std::accumulate( columns.begin(), columns.end(), 0,
                                          [&]( auto sum, auto & column ) -> auto
    {
      std::string sumOfString = stringutilities::join( column.m_maxStringSize, spaces );
      return static_cast< decltype(sum) >(sum + sumOfString.length());
    } );
  }

  string::size_type maxTopLineLength =  tableTitle.length();
  maxTopLineLength = std::max( {maxTopLineLength, sectionlineLength} );
}


void TableTextFormatter::buildTableSeparators( std::vector< TableLayout::Column > const & columns,
                                               std::string & sectionSeparatingLine,
                                               std::string & topSeparator ) const
{
  integer const columnMargin = m_tableLayout.getColumnMargin();
  integer const borderMargin = m_tableLayout.getBorderMargin();

  std::vector< string > maxStringsPerColumn;
  for( auto const & column : columns )
  {
    std::for_each( column.m_maxStringSize.begin(), column.m_maxStringSize.end(), [&] ( std::string maxString ) {
      maxStringsPerColumn.push_back( string( maxString.length(), m_horizontalLine ) );
    } );
  }

  std::string const patternBetweenColumns = GEOS_FMT( "{:-^{}}", m_horizontalLine, columnMargin );
  std::string const leftBorder = GEOS_FMT( "{}{:-<{}}", m_horizontalLine, "", borderMargin - 1 );
  std::string const rightBorder = GEOS_FMT( "{}{:-<{}}", m_horizontalLine, "", borderMargin - 1 );

  // Join all parts to form the section separating line
  std::string const columnJoin = stringutilities::join( maxStringsPerColumn, patternBetweenColumns );
  std::ostringstream oss;
  oss << leftBorder << columnJoin << rightBorder;
  sectionSeparatingLine = oss.str();

  integer const topSeparatorLength = sectionSeparatingLine.size() - 2; // Adjust for border characters
  topSeparator = GEOS_FMT( "{}{:-<{}}{}", m_horizontalLine, "", topSeparatorLength, m_horizontalLine );
}

void TableTextFormatter::outputSubSection( std::vector< TableLayout::Column > const & columns,
                                           std::ostringstream & tableOutput,
                                           integer idxRow ) const
{
  integer const columnMargin = m_tableLayout.getColumnMargin();
  for( size_t idxCol = 0; idxCol< columns.size(); ++idxCol )
  {
    tableOutput << buildCell( columns[idxCol].m_parameter.alignment,
                              columns[idxCol].m_columnValues[idxRow],
                              columns[idxCol].m_maxStringSize[0].length() );
    tableOutput << GEOS_FMT( "{:^{}}", m_verticalLine, columnMargin );
  }
}

void TableTextFormatter::outputValuesSectionRows( std::vector< TableLayout::Column > const & columns,
                                                  std::ostringstream & tableOutput,
                                                  size_t const nbRows,
                                                  std::string const & sectionSeparatingLine ) const
{
  integer const columnMargin = m_tableLayout.getColumnMargin();
  integer const borderMargin = m_tableLayout.getBorderMargin();
  string const spaces =  std::string( columnMargin, ' ' );

  for( size_t idxRow = 0; idxRow < nbRows; ++idxRow )
  {
    // Append the left border
    tableOutput << GEOS_FMT( "{:<{}}", m_verticalLine, borderMargin );

    for( std::size_t idxColumn = 0; idxColumn < columns.size(); ++idxColumn )
    {
      auto const & column = columns[idxColumn];

      if( !column.subColumn.empty())
      {
        outputSubSection( column.subColumn, tableOutput, idxRow );
      }
      else
      {
        string cell = column.m_columnValues.at( idxRow );
        std::string cellSize = stringutilities::join( column.m_maxStringSize, spaces );
        tableOutput << buildCell( column.m_parameter.alignment, cell, cellSize.length());

        // Add space between column
        if( idxColumn < columns.size() - 1 )
        {
          tableOutput << GEOS_FMT( "{:^{}}", m_verticalLine, columnMargin );
        }
      }

    }
    // Append right border with line return
    tableOutput << GEOS_FMT( "{:>{}}\n", m_verticalLine, borderMargin );
  }

  if( nbRows != 0 )
  {
    tableOutput << GEOS_FMT( "{}\n", sectionSeparatingLine );
  }
}

void TableTextFormatter::outputTitleRow( std::ostringstream & tableOutput,
                                         std::string const & topSeparator ) const
{
  string const tableTitle = string( m_tableLayout.getTitle());
  if( !tableTitle.empty() )
  {
    tableOutput << GEOS_FMT( "{}\n", topSeparator );
    tableOutput << GEOS_FMT( "{:<{}}", m_verticalLine, m_tableLayout.getBorderMargin());
    tableOutput << buildCell( TableLayout::Alignment::left, tableTitle, (topSeparator.length() - 6));
    tableOutput << GEOS_FMT( "{:>{}}\n", m_verticalLine, m_tableLayout.getBorderMargin() );
  }
}

void TableTextFormatter::outputHeaderSectionRows( std::vector< TableLayout::Column > const & columns,
                                                  std::ostringstream & tableOutput,
                                                  size_t const nbRows,
                                                  std::string const & sectionSeparatingLine ) const
{
  integer const columnMargin = m_tableLayout.getColumnMargin();
  integer const borderMargin = m_tableLayout.getBorderMargin();
  string const spaces =  std::string( columnMargin, ' ' );
  bool containSubColumn = false;

  for( size_t idxRow = 0; idxRow < nbRows; ++idxRow )
  {
    // Append the left border
    tableOutput << GEOS_FMT( "{:<{}}", m_verticalLine, borderMargin );

    for( std::size_t idxColumn = 0; idxColumn < columns.size(); ++idxColumn )
    {
      auto const & column = columns[idxColumn];
      string cell = column.m_parameter.splitColumnNames.at( idxRow );
      std::string cellSize =  stringutilities::join( column.m_maxStringSize, spaces );
      tableOutput << buildCell( column.m_parameter.alignment, cell, cellSize.length());

      // Add space between columns
      if( idxColumn < columns.size() - 1 )
      {
        tableOutput << GEOS_FMT( "{:^{}}", m_verticalLine, columnMargin );
      }

      if( !column.subColumn.empty())
      {
        containSubColumn = true;
      }
    }
    // Append right border with line return
    tableOutput << GEOS_FMT( "{:>{}}\n", m_verticalLine, borderMargin );
  }

  if( nbRows != 0 )
  {
    tableOutput << GEOS_FMT( "{}\n", sectionSeparatingLine );
  }

  // Check and build subrow header
  if( containSubColumn )
  {
    std::vector< TableLayout::Column > rowSubColumns;

    for( auto const & column : columns )
    {
      if( column.subColumn.empty())
      {
        rowSubColumns.push_back( {TableLayout::ColumnParam{""}, {}, column.m_maxStringSize, {}} );
      }
      else
      {
        rowSubColumns.insert( rowSubColumns.end(), column.subColumn.begin(), column.subColumn.end());
      }
    }
    outputHeaderSectionRows( rowSubColumns, tableOutput, 1, sectionSeparatingLine );
  }
}


}
