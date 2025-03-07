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
 * @file TableData.hpp
 */
#include "TableLayout.hpp"
#include <numeric>

namespace geos
{

void TableLayout::addToColumns( std::vector< string > const & columnNames )
{
  for( auto const & m_header : columnNames )
  {
    addToColumns( m_header );
  }
}

void TableLayout::addToColumns( string_view m_header )
{
  TableLayout::Column column = TableLayout::Column().setName( m_header );
  m_tableColumnsData.push_back( column );
}

void TableLayout::addToColumns( TableLayout::Column const & column )
{
  m_tableColumnsData.push_back( column );
}

TableLayout & TableLayout::setTitle( string_view title )
{
  m_tableTitle = CellLayout( CellType::Header, title, Alignment::center );
  return *this;
}

TableLayout & TableLayout::enableLineBreak( bool value )
{
  m_wrapLine = value;
  return *this;
}

TableLayout & TableLayout::setMargin( MarginValue marginValue )
{
  m_marginValue = marginValue;
  m_borderMargin = marginValue;
  m_columnMargin = integer( marginValue ) * 2 + 1;

  return *this;
}

bool TableLayout::isLineBreakEnabled() const
{ return m_wrapLine; }

size_t TableLayout::getMaxDepth() const
{
  size_t depthMax = 1;
  size_t currDepth = 1;
  for( auto const & column : m_tableColumnsData )
  {
    currDepth = 1;
    TableLayout::Column const * currColumn = &column;
    while( !currColumn->m_subColumns.empty())
    {
      currColumn = &currColumn->m_subColumns[0];
      currDepth++;
    }
    depthMax = std::max( currDepth, depthMax );
  }
  return depthMax;
}

void divideCell( std::vector< string > & lines, string_view value )
{
  std::istringstream strStream{ string( value ) };
  std::string line;
  lines.clear();
  while( getline( strStream, line, '\n' ))
  {
    lines.push_back( line );
  }

  if( line.empty())
  {
    lines.push_back( "" );
  }
}

size_t computeCellWidth( std::vector< string > const & lines )
{
  size_t cellWidth = 0;;
  if( !lines.empty())
  {
    cellWidth = std::max_element( lines.begin(), lines.end(), []( const auto & a, const auto & b )
    {
      return a.length() < b.length();
    } )->length();
  }
  return cellWidth;
}

TableLayout::CellLayout::CellLayout():
  m_lines( {""} ),
  m_cellType( CellType::Header ),
  m_alignment( TableLayout::Alignment::center ),
  m_cellWidth( 0 )
{}

TableLayout::CellLayout::CellLayout( CellType const cellType ):
  m_lines( {""} ),
  m_cellType( cellType ),
  m_alignment( TableLayout::Alignment::center ),
  m_cellWidth( 0 )
{}

TableLayout::CellLayout::CellLayout( CellType type, string_view cellValue, TableLayout::Alignment alignment ):
  m_cellType( type ),
  m_alignment( alignment )
{
  divideCell( m_lines, cellValue );
  m_cellWidth = computeCellWidth( m_lines );
}

TableLayout::Column::Column():
  m_parent( nullptr ), m_next( nullptr )
{
  m_header.m_lines = {};
  m_header.m_cellType  = CellType::Header;
  m_header.m_alignment = Alignment::center;
}

TableLayout::Column::Column( TableLayout::CellLayout cell ):
  m_parent( nullptr ), m_next( nullptr )
{ m_header = cell; }


TableLayout::Column & TableLayout::Column::setName( string_view name )
{
  m_header.m_lines.push_back( std::string( name ) );
  divideCell( m_header.m_lines, m_header.m_lines[0] );
  m_header.m_cellWidth = computeCellWidth( m_header.m_lines );
  m_header.m_cellType = CellType::Header;
  return *this;
}

TableLayout::Column & TableLayout::Column::setVisibility( CellType celltype )
{
  m_header.m_cellType = celltype;
  return *this;
}

/**
 * @brief Creates a vector of sub-columns from a list of names.
 * @tparam CONTAINER Container type of sub-column names (e.g. std::vector<std::string>).
 * @param names Sub-column names list.
 * @param alignment alignment of the sub-columns to create.
 * @return A vector of TableLayout::Column, ready to use for TableLayout::Column::m_subColumns.
 */
template< typename CONTAINER >
static std::vector< TableLayout::Column > makeSubColumnsFromStrings( CONTAINER const & names,
                                                                     TableLayout::ColumnAlignement const alignment )
{
  std::vector< TableLayout::Column > subColumns;
  subColumns.reserve( names.size());
  for( auto const & name : names )
  {
    TableLayout::Column col{ TableLayout::CellLayout( CellType::Header, name, alignment.headerAlignment ) };
    col.m_alignment = alignment;
    subColumns.emplace_back( col );
  }
  return subColumns;
}

TableLayout::Column & TableLayout::Column::addSubColumns( std::initializer_list< string > subColNames )
{
  m_subColumns = makeSubColumnsFromStrings( subColNames, m_alignment );
  return *this;
}

TableLayout::Column & TableLayout::Column::addSubColumns( std::vector< string > const & subColNames )
{
  m_subColumns = makeSubColumnsFromStrings( subColNames, m_alignment );
  return *this;
}

TableLayout::Column & TableLayout::Column::addSubColumns( std::initializer_list< TableLayout::Column > subCol )
{
  m_subColumns = subCol;
  return *this;
}

TableLayout::Column & TableLayout::Column::addSubColumns( string const & subColName )
{
  TableLayout::CellLayout cell{CellType::Header, subColName, TableLayout::Alignment::center};
  TableLayout::Column col{cell};
  col.m_alignment = m_alignment;
  this->m_subColumns.push_back( col );
  return *this;
}

TableLayout::Column & TableLayout::Column::setHeaderAlignment( Alignment headerAlignment )
{
  m_alignment.headerAlignment = headerAlignment;
  m_header.m_alignment = headerAlignment;
  std::vector< Column > & currColumns = m_subColumns;
  if( !currColumns.empty() )
  {
    for( auto & subColumn : currColumns )
    {
      subColumn.setHeaderAlignment( headerAlignment );
    }
  }
  return *this;
}

TableLayout::Column & TableLayout::Column::setValuesAlignment( Alignment valueAlignment )
{
  m_alignment.valueAlignment = valueAlignment;
  std::vector< Column > & currColumns = m_subColumns;
  if( !currColumns.empty() )
  {
    for( auto & subColumn : m_subColumns )
    {
      subColumn.setValuesAlignment( valueAlignment );
    }
  }
  return *this;
}

TableLayout::DeepFirstIterator & TableLayout::DeepFirstIterator::operator++()
{
  if( m_currentColumn->getNext() != nullptr )
  {
    m_currentColumn = m_currentColumn->getNext();
    while( m_currentColumn->hasChild() )
    {
      m_currentLayer++;
      m_currentColumn = &m_currentColumn->m_subColumns[0];
    }
  }
  else
  {
    bool const hasParent = (m_currentColumn->getParent() != nullptr);
    m_currentLayer -= size_t( hasParent );
    m_currentColumn = hasParent ? m_currentColumn->getParent() : nullptr;
  }
  return *this;
}

TableLayout::DeepFirstIterator TableLayout::DeepFirstIterator::operator++( int )
{
  TableLayout::DeepFirstIterator temp = *this;
  ++(*this);
  return temp;
}

TableLayout::DeepFirstIterator TableLayout::beginDeepFirst()
{
  TableLayout::Column * startColumn = &(*m_tableColumnsData.begin());
  size_t idxLayer = 0;
  if( startColumn->hasChild() )
  {
    while( startColumn->hasChild() )
    {
      idxLayer++;
      startColumn = &startColumn->m_subColumns[0];
    }
  }
  return DeepFirstIterator( startColumn, idxLayer );
}

}
