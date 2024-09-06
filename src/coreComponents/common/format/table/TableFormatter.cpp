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

void transpose( std::vector< std::vector< std::string > > & dest, std::vector< std::vector< std::string > > & source )
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
void formatColumnsFromLayout( std::vector< TableLayout::Column > & columns,
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

void TableTextFormatter::fillTableColumnsFromRows( std::vector< TableLayout::Column > & columns,
                                                   std::vector< std::vector< string > > const & columnsValues ) const
{
  size_t realNbColumn = 0;
  size_t columnWithSubColumns = 0;

  // Compute index column because we will iterate through columnsValues
  for( auto & column : columns )
  {
    if( column.subColumn.empty())
    {
      ++realNbColumn;
    }
    else
    {
      ++columnWithSubColumns;
      realNbColumn += columnWithSubColumns - column.m_parameter.subColumns.size();
    }
  }

  std::vector< std::vector< std::string > > subColumnValues;
  auto currentColumn = 0;
  for( size_t idxRow = 0; idxRow < realNbColumn; idxRow++ )
  {
    auto & column = columns[idxRow];
    if( column.subColumn.size() > 0 )
    {
      auto startSubColumn = columnsValues.begin() + currentColumn;
      auto endSubColumn = columnsValues.begin() + currentColumn + column.m_parameter.subColumns.size();
      subColumnValues = std::vector< std::vector< std::string > >( startSubColumn, endSubColumn );
      fillTableColumnsFromRows( column.subColumn, subColumnValues );
      currentColumn += column.subColumn.size();
    }
    else
    {
      column.m_columnValues = columnsValues[currentColumn++];
    }

  }

  for( const auto & row : columnsValues )
  {
    for( const auto & value : row )
    {
      std::cout << value << " ";
    }
    std::cout << std::endl;
  }
  std::cout << "End fillTableColumnsFromRows" << std::endl;
}

string TableTextFormatter::layoutToString() const
{
  std::ostringstream tableOutput;
  string sectionSeparatingLine;
  std::vector< TableLayout::Column > columns = m_tableLayout.getColumns();

  outputLayout( tableOutput, columns, {}, sectionSeparatingLine );
  return tableOutput.str();
}

template<>
string TableTextFormatter::toString< TableData >( TableData const & tableData ) const
{
  std::ostringstream tableOutput;
  string sectionSeparatingLine;
  std::vector< TableLayout::Column > columns         = m_tableLayout.getColumns();
  std::vector< std::vector< string > > tableDataRows = tableData.getTableDataRows();
  std::vector< string > const & msgTableError        = tableData.getErrorMsgs();
  integer const nbValuesRows                         = tableDataRows.size();
  formatColumnsFromLayout( columns, tableDataRows );

  std::vector< std::vector< std::string > > columnsValues( tableDataRows[0].size(), std::vector< std::string >( tableDataRows.size()));
  transpose( columnsValues, tableDataRows );

  fillTableColumnsFromRows( columns, columnsValues );

  outputLayout( tableOutput, columns, msgTableError, sectionSeparatingLine );
  outputSectionRows( columns, sectionSeparatingLine, tableOutput, nbValuesRows, 0, TableLayout::Section::values );
  tableOutput << '\n';

  return tableOutput.str();
}

void TableTextFormatter::outputLayout( std::ostringstream & tableOutput,
                                       std::vector< TableLayout::Column > & columns,
                                       std::vector< string > const & msgTableError,
                                       string & sectionSeparatingLine ) const
{
  string topSeparator;
  size_t nbHeaderRows = 0;
  std::vector< std::vector< string > > splitHeaders;
  string const tableTitle = string( m_tableLayout.getTitle());

  splitAndSetColumnNames( columns, nbHeaderRows, splitHeaders );
  std::cout << "split finished " << std::endl;
  std::string maxStringSize = " ";
  findAndSetMaxStringSize( columns, maxStringSize );
  std::cout << "find finished " << std::endl;

  computeTableWidth( columns, msgTableError );
  std::cout << "computeTableWidth finished " << std::endl;
  buildTableSeparators( columns, topSeparator, sectionSeparatingLine );
  std::cout << "buildTableSeparators finished " << std::endl;

  tableOutput << '\n';
  outputTopRows( tableOutput, {tableTitle}, topSeparator, TableLayout::Alignment::center );
  outputTopRows( tableOutput, msgTableError, topSeparator, TableLayout::Alignment::left );
  tableOutput << GEOS_FMT( "{}\n", sectionSeparatingLine );

  outputSectionRows( columns, sectionSeparatingLine, tableOutput, nbHeaderRows, 0, TableLayout::Section::header );
  std::cout << "outputSectionRows finished " << tableOutput.str() <<std::endl;
}

void TableTextFormatter::splitAndSetColumnNames( std::vector< TableLayout::Column > & columns,
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
      splitAndSetColumnNames( column.subColumn, nbHeaderSubColRows, splitSubColHeaders );
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
    columns[&headerParts - &splitHeaders[0]].m_parameter.splitColumnNameLines = headerParts;
  }

}

void TableTextFormatter::findAndSetMaxStringSize( std::vector< TableLayout::Column > & columns, std::string & maxStringSize ) const
{
  std::vector< std::string > maxStringsColumn;
  maxStringsColumn.reserve( columns.size());

  // max header columns length 
  for( auto & column : columns )
  {
    auto const maxStringSizeHeader = *std::max_element( column.m_parameter.splitColumnNameLines.begin(),
                                                        column.m_parameter.splitColumnNameLines.end(),
                                                        []( const auto & a, const auto & b ) {return a.size() < b.size();} );
    column.m_maxStringSize = maxStringSizeHeader;
    maxStringsColumn.push_back( column.m_parameter.columnName );
  }

  for( auto & column : columns )
  {

    for( auto const & cell : column.m_columnValues )
    {
      if( column.m_maxStringSize.length() < cell.length())
      {
        column.m_maxStringSize = cell;
      }
    }

    if( !column.subColumn.empty())
    {
      findAndSetMaxStringSize( column.subColumn, column.m_maxStringSize );
    }
    std::cout << "verif max string size before " << column.m_maxStringSize << std::endl;
    if( column.m_maxStringSize.length() > maxStringSize.length())
    {
      maxStringSize = std::accumulate( maxStringsColumn.begin(), maxStringsColumn.end(), std::string());
    }

    std::cout << "verif max string size after " << column.m_maxStringSize << std::endl;

  }
}

void TableTextFormatter::increaseColumnsSize( std::vector< TableLayout::Column > & columns,
                                              integer const extraCharacters ) const
{
  integer const extraCharactersPerColumn = std::ceil( extraCharacters / columns.size() );
  for( std::size_t idxColumn = 0; idxColumn < columns.size(); ++idxColumn )
  {
    integer newMaxStringSize = extraCharactersPerColumn + columns[idxColumn].m_maxStringSize.size();
    if( idxColumn == columns.size() - 1 )
    {
      newMaxStringSize += m_tableLayout.getColumnMargin();
    }

    columns[idxColumn].m_maxStringSize = GEOS_FMT( "{:>{}}",
                                                   columns[idxColumn].m_maxStringSize,
                                                   newMaxStringSize );
  }
}

void computeTableSectionLength( std::vector< TableLayout::Column > & columns, string::size_type & sectionlineLength, string_view maxColumnStringSize )
{
  sectionlineLength += std::accumulate( columns.begin(), columns.end(), 0,
                                        [&]( auto sum, auto & column ) -> auto
  {
    if( !column.subColumn.empty())
    {
      computeTableSectionLength( column.subColumn, sectionlineLength, column.m_maxStringSize );
    }

    if( !column.m_maxStringSize.compare( maxColumnStringSize ))
    {
      return sum;
    }
    std::cout << " alors je calcul" << column.m_maxStringSize << " egal " << sum <<  std::endl;
    return static_cast< decltype(sum) >(sum + column.m_maxStringSize.length());
  } );
  std::cout << " fincalcul" << sectionlineLength <<  std::endl;
}

void TableTextFormatter::computeTableWidth( std::vector< TableLayout::Column > & columns,
                                            std::vector< string > const & msgTableError ) const
{
  integer const columnMargin = m_tableLayout.getColumnMargin();
  integer const borderMargin = m_tableLayout.getBorderMargin();
  string const tableTitle = string( m_tableLayout.getTitle() );

  // compute table error msg length
  string::size_type msgTableErrorLength = borderMargin;
  {
    if( !msgTableError.empty() )
    {
      auto maxStringSize = *(std::max_element( msgTableError.begin(), msgTableError.end(),
                                               []( const auto & a, const auto & b ) {
        return a.size() < b.size();
      } ));

      msgTableErrorLength += maxStringSize.size() + 1; // max string size + line return at the end
    }
  }

  // compute table section length
  string::size_type sectionLengthWithSpacing = ( ( columns.size() - 1 ) *  columnMargin ) + (borderMargin * 2);
  string::size_type sectionlineLength = sectionLengthWithSpacing;
  computeTableSectionLength( columns, sectionlineLength, " " );

  string::size_type maxTopLineLength =  tableTitle.length();
  maxTopLineLength = std::max( maxTopLineLength, msgTableErrorLength );
  if( sectionlineLength < maxTopLineLength )
  {
    integer const extraCharacters = maxTopLineLength - sectionlineLength;
    increaseColumnsSize( columns, extraCharacters );
  }
}


void TableTextFormatter::buildTableSeparators( std::vector< TableLayout::Column > const & columns,
                                               string & topSeparator,
                                               string & sectionSeparatingLine ) const
{
  {   // section separator line
    integer const columnMargin = m_tableLayout.getColumnMargin();
    integer const borderMargin = m_tableLayout.getBorderMargin();

    std::vector< string > maxStringPerColumn;
    for( auto const & column : columns )
    {
      maxStringPerColumn.push_back( string( column.m_maxStringSize.length(), m_horizontalLine ) );
    }

    string const patternBetweenColumns = GEOS_FMT( "{:-^{}}", m_horizontalLine, columnMargin );

    std::string const leftBorder = GEOS_FMT( "{}{:-<{}}", m_horizontalLine, "", borderMargin );
    std::string const rightBorder = GEOS_FMT( "{}{:-<{}}", m_horizontalLine, "", borderMargin );
    std::string const columnJoin = stringutilities::join( maxStringPerColumn, patternBetweenColumns );

    std::ostringstream oss;
    oss << leftBorder << columnJoin << rightBorder;
    sectionSeparatingLine = oss.str();
  }

  {   // top line separator
      // -2 because we can have '+' to the extremity of top separator
    integer const topSeparatorLength = sectionSeparatingLine.size() - 2;
    topSeparator = GEOS_FMT( "{}{:-<{}}{}", m_horizontalLine, "", topSeparatorLength, m_horizontalLine );
  }
}

void TableTextFormatter::outputTopRows( std::ostringstream & tableOutput,
                                        std::vector< string > const & msg,
                                        string_view topSeparator,
                                        TableLayout::Alignment alignment ) const
{
  if( msg.size() != 0 && msg[0] != "" )
  {
    tableOutput << GEOS_FMT( "{}\n", topSeparator );
    for( std::string const & str : msg )
    {
      tableOutput << m_verticalLine << string( m_tableLayout.getBorderMargin(), ' ' );
      tableOutput << buildCell( alignment, str, (topSeparator.length() - 6));
      tableOutput << string( m_tableLayout.getBorderMargin(), ' ' ) << "|\n";
    }
  }
}

void TableTextFormatter::outputSectionRows( std::vector< TableLayout::Column > const & columns,
                                            string_view sectionSeparatingLine,
                                            std::ostringstream & tableOutput,
                                            integer const nbRows,
                                            size_t const nbSubColumns,
                                            TableLayout::Section const section ) const
{
  integer const columnMargin = m_tableLayout.getColumnMargin();
  integer const borderMargin = m_tableLayout.getBorderMargin();
  size_t const nbColumns = nbSubColumns + columns.size();
  std::cout << nbRows << std::endl;
  std::cout << nbColumns << std::endl;
  std::cout << "testeuh"<< std::endl;
  for( auto idxRow = 0; idxRow< nbRows; ++idxRow )
  {
    // Append the left border
    tableOutput << GEOS_FMT( "{:<{}}", m_verticalLine, 1 +  borderMargin );
    for( std::size_t idxColumn = 0; idxColumn < nbColumns; ++idxColumn )
    {
      auto & column = columns[idxColumn];
      TableLayout::Column const currentColumn = column;
      auto const & columnContent = section == TableLayout::Section::header ?
                                   column.m_parameter.splitColumnNameLines :
                                   column.m_columnValues;
      string cell = columnContent.at( idxRow );
      integer const cellSize = currentColumn.m_maxStringSize.length();

      tableOutput << buildCell( currentColumn.m_parameter.alignment, cell, cellSize );

      // Add space between column
      if( idxColumn < columns.size() - 1 )
      {
        tableOutput << GEOS_FMT( "{:^{}}", m_verticalLine, columnMargin );
      }
    }

    // Append right border with line return
    tableOutput << GEOS_FMT( "{:>{}}\n", m_verticalLine, borderMargin + 1 );

    // if( !column.subColumn.empty() )
    // {
    //   outputSectionRows( column.subColumn, sectionSeparatingLine, tableOutput, nbRows, columns.size(), section );
    // }
  }

  if( nbRows != 0 )
  {
    tableOutput << GEOS_FMT( "{}\n", sectionSeparatingLine );
  }


}

}
