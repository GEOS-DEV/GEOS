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
  m_columns = tableLayout.m_columns;
}

void TableFormatter::fillTableColumnsFromRows( std::vector< std::vector< string > > const & rows )
{
  std::vector< Column > & columns = tableLayout.m_columns();

  for( size_t idxRow = 0; idxRow < rows.size(); idxRow++ )
  {
    for( size_t idxColumn = 0; idxColumn < columns.size(); idxColumn++ )
    {
      columns[idxColumn].columnValues.push_back( rows[idxRow][idxColumn] );
    }
  }

}

//----------------------------------------------//

TableCSVFormatter::TableCSVFormatter( TableLayout tableLayout )
{
  m_columns = tableLayout.m_columns;

}

string TableCSVFormatter::headerToString()
{
  string headerValues = "";

  for( std::size_t idxColumn = 0; idxColumn < m_columns.size(); ++idxColumn )
  {
    headerValues += m_columns[idxColumn].parameter.headerName[idxRow];
  }

  return headerValues;
}

string TableCSVFormatter::dataToString( TableData2D tableData )
{
  std::vector< std::vector< string > > rowsValues = tableData.m_rows;
  string data;

  tableData.buildRows();

  for( std::vector< string > row : rowsValues )
  {
    for( string value : row )
    {
      data += value + ",";
    }
    data += "\n";
  }
  return data;
}

string TableCSVFormatter::dataToString( TableData tableData )
{
  std::vector< std::vector< string > > rowsValues = tableData.m_rows;
  string data;

  fillTableColumnsFromRows( rowsValues );

  for( std::vector< string > row : rowsValues )
  {
    for( string value : row )
    {
      data += value + ",";
    }
    data += "\n";
  }
  return data;
}

//----------------------------------------------//

TableTextFormatter::TableTextFormatter( TableLayout tableLayout ):
  TableFormatter( tableLayout )
{}

string TableTextFormatter::ToString( TableData2D & tableData )
{
  tableData.buildRows();
  getTableBuilt( tableData.m_rows );
}

string TableTextFormatter::ToString( TableData & tableData )
{
  return getTableBuilt( tableData.m_rows );
}

string TableTextFormatter::getTableBuilt( std::vector< std::vector< string > > rowsValues )
{
  string titleRows;
  string topSeparator;
  string sectionSeparator;

  size_t largestHeaderVectorSize = 0;
  std::vector< std::vector< string > > splitHeader;
  fillTableColumnsFromRows( rowsValues );

  parseAndStoreHeaderSections( largestHeaderVectorSize, splitHeader );
  adjustHeaderSizesAndStore( largestHeaderVectorSize, splitHeader );

  findAndSetMaxStringSize();
  computeAndBuildSeparator( topSeparator, sectionSeparator );

  if( !tableTitle.empty())
  {
    buildTitleRow( titleRows, topSeparator, sectionSeparator );
  }

  string tableRows += GEOS_FMT( "{}\n", sectionSeparator );
  buildSectionRows( sectionSeparator, tableRows, largestHeaderVectorSize, TableLayout::Section::header );
  buildSectionRows( sectionSeparator, tableRows, rowsValues.size(), TableLayout::Section::values );

  string const tableOutput = titleRows + tableRows + '\n';

  return tableOutput;
}

void TableTextFormatter::parseAndStoreHeaderSections( size_t & largestHeaderVectorSize,
                                                      std::vector< std::vector< string > > & splitHeader )
{
  for( size_t columnParamIdx = 0; columnParamIdx< m_columns.size(); columnParamIdx++ )
  {
    std::vector< string > splitHeaderParts;

    std::istringstream ss( tableLayout.m_columns[columnParamIdx].parameter.headerName[0] );
    string subHeaderName;

    while( getline( ss, subHeaderName, '\n' ))
    {
      splitHeaderParts.push_back( subHeaderName );
    }

    size_t const cellSize = splitHeaderParts.size();
    largestHeaderVectorSize = std::max( largestHeaderVectorSize, cellSize );

    splitHeader.push_back( splitHeaderParts );
  }
}

void TableTextFormatter::adjustHeaderSizesAndStore( size_t largestHeaderVectorSize,
                                                    std::vector< std::vector< string > > & splitHeader )
{
  for( size_t columnParamIdx = 0; columnParamIdx < m_columns.size(); columnParamIdx++ )
  {
    if( splitHeader[columnParamIdx].size() < largestHeaderVectorSize )
    {
      integer const whiteRowToAdd = largestHeaderVectorSize - splitHeader[columnParamIdx].size();
      splitHeader[columnParamIdx].insert( splitHeader[columnParamIdx].end(), whiteRowToAdd, " " );
    }
    m_columns[columnParamIdx].parameter.headerName = splitHeader[columnParamIdx];
  }
}

void TableTextFormatter::findAndSetMaxStringSize()
{
  string maxStringSize = "";
  for( size_t idxColumn  = 0; idxColumn <  m_columns.size(); idxColumn++ )
  {
    auto it = std::max_element( m_columns[idxColumn].parameter.headerName.begin(),
                                m_columns[idxColumn].parameter.headerName.end(),
                                []( const auto & a, const auto & b ) {
      return a.size() < b.size();
    } );

    maxStringSize = *it;

    for( size_t idxRow = 0; idxRow <  rowsValues.size(); idxRow++ )
    {
      string cell =  m_columns[idxColumn].columnValues[idxRow];
      if( maxStringSize.length() < cell.length())
      {
        maxStringSize = cell;
      }
    }
    m_columns[idxColumn].m_maxStringSize = maxStringSize;
  }
}

void TableTextFormatter::computeAndSetMaxStringSize( string::size_type sectionlineLength,
                                                     string::size_type titleLineLength )
{
  integer extraLinesPerColumn;
  integer extraLines;
  integer newStringSize;

  extraLines = titleLineLength - sectionlineLength;
  extraLinesPerColumn = std::ceil( extraLines / m_columns.size() );

  for( std::size_t idxColumn = 0; idxColumn < m_columns.size(); ++idxColumn )
  {
    newStringSize = extraLinesPerColumn + m_columns[idxColumn].m_maxStringSize.size();
    if( idxColumn == m_columns.size() - 1 ||  m_columns.size() == 1 )
    {
      m_columns[idxColumn].m_maxStringSize = GEOS_FMT( "{:>{}}",
                                                       m_columns[idxColumn].m_maxStringSize,
                                                       newStringSize + columnMargin );
    }
    else
    {
      m_columns[idxColumn].m_maxStringSize = GEOS_FMT( "{:>{}}",
                                                       m_columns[idxColumn].m_maxStringSize,
                                                       newStringSize );
    }
  }
}

void TableTextFormatter::computeAndBuildSeparator( string & topSeparator, string & sectionSeparator )
{
  string::size_type sectionlineLength = 0;
  string::size_type titleLineLength = tableTitle.length() + ( marginTitle * 2 );
  integer nbSpaceBetweenColumn = ( ( m_columns.size() - 1 ) *  columnMargin ) + (borderMargin * 2);
  if( !tableTitle.empty())
  {
    tableTitle = GEOS_FMT( "{:^{}}", tableTitle, titleLineLength );
  }

  for( std::size_t i = 0; i < m_columns.size(); ++i )
  {
    sectionlineLength += m_columns[i].m_maxStringSize.length();
  }

  sectionlineLength += nbSpaceBetweenColumn;
  if( sectionlineLength < titleLineLength )
  {
    computeAndSetMaxStringSize( sectionlineLength, titleLineLength );
  }
  if( m_columns.size() == 1 )
  {
    sectionSeparator +=  GEOS_FMT( "+{:-<{}}+",
                                   "",
                                   ( m_columns[0].m_maxStringSize.length() + (borderMargin - 1) + columnMargin ));
  }
  else
  {
    for( std::size_t idxColumn = 0; idxColumn < m_columns.size(); ++idxColumn )
    {
      integer const cellSize = m_columns[idxColumn].m_maxStringSize.length();
      if( idxColumn == 0 )
      {
        sectionSeparator +=  GEOS_FMT( "+{:-<{}}", "", ( cellSize + borderMargin ));
      }
      else if( idxColumn == (m_columns.size() - 1))
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
  titleRows +=  buildValueCell( Alignment::middle,
                                tableTitle,
                                (sectionSeparator.length() - 2) // -2 for ||
                                );
  titleRows += GEOS_FMT( "{}\n", "|" );
}

void TableTextFormatter::buildSectionRows( string_view sectionSeparator,
                                           string & tableRows,
                                           integer const nbRows,
                                           TableLayout::Section const section )
{
  for( integer idxRow = 0; idxRow< nbRows; idxRow++ )
  {
    tableRows += GEOS_FMT( "{:<{}}", "|", 1 +  borderMargin );
    for( std::size_t idxColumn = 0; idxColumn < m_columns.size(); ++idxColumn )
    {
      string cell;

      if( section == Section::header )
      {
        cell = m_columns[idxColumn].parameter.headerName[idxRow];
      }
      else
      {
        cell = m_columns[idxColumn].columnValues[idxRow];
      }
      integer const cellSize = m_columns[idxColumn].m_maxStringSize.length();
      tableRows += buildValueCell( m_columns[idxColumn].parameter.alignment,
                                   cell,
                                   cellSize );

      if( idxColumn < m_columns.size() - 1 )
      {
        tableRows += GEOS_FMT( "{:^{}}", "|", columnMargin );
      }

    }
    if( m_columns.size() == 1 )
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

}
