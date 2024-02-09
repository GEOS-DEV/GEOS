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
 * @file Table.cpp
 */

#include "Table.hpp"

namespace geos
{

std::string const cellAlignment( Table::Alignment const & a, std::string const & value, int spaces )
{
  switch( a )
  {
    case Table::right:   return GEOS_FMT( "{:>{}}", value, spaces );
    case Table::left:   return GEOS_FMT( "{:<{}}", value, spaces );
    case Table::middle: return GEOS_FMT( "{:^{}}", value, spaces );
    default:      return GEOS_FMT( "{:<{}}", value, spaces );
  }
}

std::string Table::getStringSection( Section section ) const
{
  switch( section )
  {
    case Section::header: return "header";
    case Section::values: return "values";
    default: return "values";
  }
}


Table::Table( std::vector< std::string > const & headers ):
  maxRowHeader( 0 )
{
  for( size_t idx = 0; idx< headers.size(); idx++ )
  {
    m_columns.push_back( {Table::ColumnParam{{headers[idx]}, Alignment::middle, true}, {}, ""} );
  }
}

Table::Table( std::vector< ColumnParam > const & columnParameter )
{
  for( size_t idx = 0; idx< columnParameter.size(); idx++ )
  {
    if( columnParameter[idx].enabled )
    {
      m_columns.push_back( {columnParameter[idx], {}, ""} );
    }
  }
  maxRowHeader = 0;
}

void Table::splitHeadersStringAndStore()
{
  for( size_t columnParamIdx = 0; columnParamIdx< m_columns.size(); columnParamIdx++ )
  {
    std::vector< std::string > splitHeaderParts;
    std::stringstream ss( m_columns[columnParamIdx].parameter.headerName[0] );
    std::string subHeaderName;

    while( getline( ss, subHeaderName, '\n' ))
    {
      splitHeaderParts.push_back( subHeaderName );
    }

    size_t cellSize = splitHeaderParts.size();
    maxRowHeader = std::max( maxRowHeader, cellSize );

    m_splitHeader.push_back( splitHeaderParts );
  }
}

void Table::addSpaceToSplitHeaderAndStore()
{
  for( size_t columnParamIdx = 0; columnParamIdx < m_columns.size(); columnParamIdx++ )
  {
    if( m_splitHeader[columnParamIdx].size() < maxRowHeader )
    {
      integer whiteRowToAdd = maxRowHeader -  m_splitHeader[columnParamIdx].size();
      m_splitHeader[columnParamIdx].insert( m_splitHeader[columnParamIdx].end(), whiteRowToAdd, " " );
    }
    m_columns[columnParamIdx].parameter.headerName = m_splitHeader[columnParamIdx];
  }
}

void Table::setTitle( std::string const & title_ )
{
  title = title_;
}

std::string const & Table::getTitle()
{
  return title;
}

void Table::findMaxStringSize()
{
  std::string maxStringSize = "";
  for( size_t idxColumn  = 0; idxColumn <  m_columns.size(); idxColumn++ )
  {
    auto it = std::max_element( m_columns[idxColumn].parameter.headerName.begin(),
                                m_columns[idxColumn].parameter.headerName.end(),
                                []( const auto & a, const auto & b ) {
      return a.size() < b.size();
    } );

    maxStringSize = *it;

    for( size_t idxRow = 0; idxRow <  m_cellsRows.size(); idxRow++ )
    {
      std::string cell =  m_columns[idxColumn].columnValues[idxRow];
      if( maxStringSize.length() < cell.length())
      {
        maxStringSize = cell;
      }
    }
    m_columns[idxColumn].m_maxStringSize = maxStringSize;
  }
}

void Table::computeAndSetMaxStringSize( string::size_type const & sectionlineLength,
                                        string::size_type const & titleLineLength )
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

void Table::computeAndBuildLines()
{
  string::size_type sectionlineLength = 0;
  string::size_type titleLineLength = title.length() + ( marginTitle * 2 );
  integer nbSpaceBetweenColumn = ( ( m_columns.size() - 1 ) *  columnMargin ) + (borderMargin * 2);

  if( !title.empty())
  {
    title = GEOS_FMT( "{:^{}}", title, titleLineLength );
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
      integer cellSize = m_columns[idxColumn].m_maxStringSize.length();

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

void Table::buildTitleRow()
{
  titleRow = GEOS_FMT( "\n{}\n|", topSeparator );
  titleRow +=  cellAlignment( Alignment::middle,
                              title,
                              (sectionSeparator.length() - 2) // -2 for ||
                              );
  titleRow += GEOS_FMT( "{}\n", "|" );
}

void Table::buildSectionRows( integer const & nbRows, Section const & sectionName )
{

  for( integer idxRow = 0; idxRow< nbRows; idxRow++ )
  {
    rows += GEOS_FMT( "{:<{}}", "|", 1 +  borderMargin );
    for( std::size_t idxColumn = 0; idxColumn < m_columns.size(); ++idxColumn )
    {
      std::string cell;

      if( getStringSection( sectionName ) == "header" )
      {
        cell = m_columns[idxColumn].parameter.headerName[idxRow];
      }
      else
      {
        cell = m_columns[idxColumn].columnValues[idxRow];
      }
      integer cellSize = m_columns[idxColumn].m_maxStringSize.length();
      rows += cellAlignment( m_columns[idxColumn].parameter.alignment,
                             cell,
                             cellSize );

      if( idxColumn < m_columns.size() - 1 )
      {
        rows += GEOS_FMT( "{:^{}}", "|", columnMargin );
      }

    }
    if( m_columns.size() == 1 )
    {
      rows +=  GEOS_FMT( "{:>{}}\n", "|", columnMargin );
    }
    else
    {
      rows += GEOS_FMT( "{:>{}}\n", "|", borderMargin + 1 );
    }

  }
  if( nbRows != 0 )
  {
    rows += GEOS_FMT( "{}\n", sectionSeparator );
  }
}

void Table::fillColumnsValuesFromCellsRows()
{
  for( size_t idxRow = 0; idxRow < m_cellsRows.size(); idxRow++ )
  {
    for( size_t idxColumn = 0; idxColumn < m_columns.size(); idxColumn++ )
    {
      m_columns[idxColumn].columnValues.push_back( m_cellsRows[idxRow][idxColumn] );
    }
  }
}

void Table::draw( std::ostringstream & oss )
{
  oss.clear();
  std::string tableOutput;

  fillColumnsValuesFromCellsRows();

  splitHeadersStringAndStore();
  addSpaceToSplitHeaderAndStore();

  findMaxStringSize();
  computeAndBuildLines();

  if( !title.empty())
  {
    buildTitleRow();
  }

  rows += GEOS_FMT( "{}\n", sectionSeparator );
  buildSectionRows( maxRowHeader, Section::header );
  buildSectionRows( m_cellsRows.size(), Section::values );

  tableOutput = titleRow + rows;

  std::cout << tableOutput;

  oss << tableOutput;
}

}
