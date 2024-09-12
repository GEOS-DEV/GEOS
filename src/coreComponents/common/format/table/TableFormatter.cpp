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
      std::cout << " nb  sous col "<< column.m_parameter.subColumns.size() << std::endl;
      std::vector< std::vector< std::string > > subColumnValues( tColumnsValues.begin() + currentColumn,
                                                                 tColumnsValues.begin() + currentColumn + column.m_parameter.subColumns.size());
      fillTableColumnsFromRows( column.subColumn, subColumnValues, true );

      std::vector< std::vector< std::string > > tSubColumnValues( subColumnValues[0].size(), std::vector< std::string >( subColumnValues.size()) );
      transpose( tSubColumnValues, subColumnValues );
      for( const auto & columnValues : tSubColumnValues )
      {
        column.m_columnValues.insert( column.m_columnValues.end(), columnValues.begin(), columnValues.end() );
      }
      currentColumn += subColumnValues.size();
    }
  }
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
  fillTableColumnsFromRows( columns, tableDataRows, false );
  outputLayout( tableOutput, columns, msgTableError, sectionSeparatingLine );
  outputSectionRows( columns, sectionSeparatingLine, tableOutput, nbValuesRows, TableLayout::Section::values );
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
  std::string maxStringSize = " ";

  splitAndSetColumnNames( columns, nbHeaderRows, splitHeaders );
  for( auto & column : columns )
  {
    findAndSetMaxStringSize( column, column.m_maxStringSize, 0 );
  }

  computeTableWidth( columns, msgTableError );
  buildTableSeparators( columns, topSeparator, sectionSeparatingLine );

  tableOutput << '\n';
  outputTopRows( tableOutput, {tableTitle}, topSeparator, TableLayout::Alignment::center );
  outputTopRows( tableOutput, msgTableError, topSeparator, TableLayout::Alignment::left );
  tableOutput << GEOS_FMT( "{}\n", sectionSeparatingLine );
  outputSectionRows( columns, sectionSeparatingLine, tableOutput, nbHeaderRows, TableLayout::Section::header );
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

void TableTextFormatter::findAndSetMaxStringSize( TableLayout::Column & column, std::vector< std::string > & maxStringSize, integer const idxMaxString ) const
{
  string maxStringColumn = "";

  {// header case
    auto const maxStringSizeHeader = *std::max_element( column.m_parameter.splitColumnNameLines.begin(),
                                                        column.m_parameter.splitColumnNameLines.end(),
                                                        []( const auto & a, const auto & b ) {return a.size() < b.size();} );
    maxStringColumn = maxStringSizeHeader;
    maxStringSize.push_back( maxStringSizeHeader );
  }

  std::cout << "- debug maxString 1- "<<  column.m_parameter.columnName << std::endl;
  for( auto const & value : maxStringSize )
  {
    std::cout << value << " ";
  }
  std::cout << std::endl;

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

  { // affectation
    if( maxStringSize[idxMaxString].length() < maxStringColumn.length() )
    {
      maxStringSize[idxMaxString] = maxStringColumn ;
    }
  }

  {  // subcolumn values case
    if( !column.subColumn.empty() )
    {
      integer idxSubColumn = 0;
      column.m_maxStringSize.clear();
      for( auto & subC : column.subColumn )
      {
      
        findAndSetMaxStringSize( subC, column.m_maxStringSize, idxSubColumn );

        std::cout << "- debug maxString - "<<  column.m_parameter.columnName << std::endl;
        for( auto const & value : column.m_maxStringSize )
        {
          std::cout << value << " ";
        }

        //TODO Condition pour check Nodes et tout les cols OU ALORS DEPLACER AU DEBUT
        std::cout << std::endl;

        idxSubColumn++;
      }
    }
  }

    std::cout << "- debug maxString 2- "<<  column.m_parameter.columnName << std::endl;   //TODO DOU VIENS LE NODE JE VOIS
  for( auto const & value : column.m_maxStringSize )
  {
    std::cout << value << " ";
  }
  std::cout << std::endl;

  if( column.m_maxStringSize.empty() )
  {
    column.m_maxStringSize.push_back( maxStringColumn );
  }
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

      msgTableErrorLength += maxStringSize.size() + 1;   // max string size + line return at the end
    }
  }

  string::size_type sectionLengthWithSpacing = ( ( columns.size() - 1 ) *  columnMargin ) + (borderMargin * 2);
  string::size_type sectionlineLength = sectionLengthWithSpacing;
  { // compute table section length TODO A FONCTION
    integer const spacesToAdd = m_tableLayout.getColumnMargin();
    bool isFirstIteration = true;
    sectionlineLength += std::accumulate( columns.begin(), columns.end(), 0,
                                          [&]( auto sum, auto & column ) -> auto
    {
      std::string maxStringSize = std::accumulate( column.m_maxStringSize.begin(), column.m_maxStringSize.end(), std::string(),
                                                   [spacesToAdd, &isFirstIteration]( std::string acc, const std::string & curr ) {
        if( isFirstIteration )
        {
          isFirstIteration = false;
          return curr;
        }
        else
        {
          return acc  +  std::string( spacesToAdd, ' ' ) + curr;
        }
      } );

      return static_cast< decltype(sum) >(sum + maxStringSize.length());
    } );
  }

  string::size_type maxTopLineLength =  tableTitle.length();
  maxTopLineLength = std::max( maxTopLineLength, msgTableErrorLength );
  maxTopLineLength = std::max( maxTopLineLength, sectionlineLength );
  std::cout << "- max length - " << sectionlineLength << " " << maxTopLineLength <<  std::endl;

  //refaire le calcul du max string par colonnes ici


}


void TableTextFormatter::buildTableSeparators( std::vector< TableLayout::Column > const & columns,
                                               string & topSeparator,
                                               string & sectionSeparatingLine ) const
{
  {   // section separator line
    integer const columnMargin = m_tableLayout.getColumnMargin();
    integer const borderMargin = m_tableLayout.getBorderMargin();

    std::vector< string > maxStringPerColumn;
    for( auto const & column : columns )  //TODO Check ca
    {
      std::for_each( column.m_maxStringSize.begin(), column.m_maxStringSize.end(), [&] ( std::string maxString ) {
        maxStringPerColumn.push_back( string( maxString.length(), m_horizontalLine ) );
      } );
    }

    std::cout << "seprator" << std::endl;
    for( auto const & seprator : maxStringPerColumn )
    {
      std::cout << seprator;
    }
    std::cout << std::endl;

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

void TableTextFormatter::outputSectionRows( std::vector< TableLayout::Column > const & columns,
                                            string_view sectionSeparatingLine,
                                            std::ostringstream & tableOutput,
                                            size_t const nbRows,
                                            TableLayout::Section const section ) const
{
  integer const columnMargin = m_tableLayout.getColumnMargin();
  integer const borderMargin = m_tableLayout.getBorderMargin();
  bool containSubColumn = false;
  for( size_t idxRow = 0; idxRow< nbRows; ++idxRow )
  {
    // Append the left border
    tableOutput << GEOS_FMT( "{:<{}}", m_verticalLine, 1 +  borderMargin );

    for( std::size_t idxColumn = 0; idxColumn < columns.size(); ++idxColumn )
    {
      auto const & column = columns[idxColumn];

      if( section == TableLayout::Section::values && !column.subColumn.empty())
      {
        outputSubSection( column.subColumn, tableOutput, idxRow );
      }
      else
      {
        TableLayout::Column const currentColumn = column;
        auto const & columnContent = section == TableLayout::Section::header ?
                                     column.m_parameter.splitColumnNameLines :
                                     column.m_columnValues;
        string cell = columnContent.at( idxRow );
        bool isFirstIteration = true;
        integer const cellSize =  std::accumulate( column.m_maxStringSize.begin(), column.m_maxStringSize.end(), 0,
                                                   [&isFirstIteration, &columnMargin]( int sum, const std::string & str ) {
          if( isFirstIteration )
          {
            isFirstIteration = false;
            return sum + str.length();
          }
          else
          {
            return sum +  columnMargin + str.length();
          }
        } );
        tableOutput << buildCell( currentColumn.m_parameter.alignment, cell, cellSize );

        // Add space between column
        if( idxColumn < columns.size() - 1 )
        {
          tableOutput << GEOS_FMT( "{:^{}}", m_verticalLine, columnMargin );
        }
      }

      if( !column.subColumn.empty() )
      {
        containSubColumn = true;
      }
    }
    // Append right border with line return
    tableOutput << GEOS_FMT( "{:>{}}\n", m_verticalLine, borderMargin + 1 );
  }

  if( nbRows != 0 )
  {
    tableOutput << GEOS_FMT( "{}\n", sectionSeparatingLine );
  }

  //check and build sub row header
  if( section == TableLayout::Section::header && containSubColumn )
  {
    std::vector< TableLayout::Column > rowSubColumns;

    for( auto const & column : columns )
    {
      if( column.subColumn.empty() )
      {
        TableLayout::Column emptyColumn = {TableLayout::ColumnParam{""}, {}, column.m_maxStringSize, {}};
        rowSubColumns.push_back( emptyColumn );
      }
      else
      {
        for( auto const & subColumn : column.subColumn )
        {
          rowSubColumns.push_back( subColumn );
        }
      }
    }
    outputSectionRows( rowSubColumns, sectionSeparatingLine, tableOutput, 1, TableLayout::Section::header );
  }

}

}
