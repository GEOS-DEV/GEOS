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

#include "common/format/StringUtilities.hpp"
#include <numeric>

namespace geos
{

void TableLayout::addColumns( std::vector< string > const & columnNames )
{
  for( auto const & columnName : columnNames )
  {
    addColumn( columnName );
  }
}

void TableLayout::addColumns( std::vector< TableLayout::Column > const & columns )
{
  for( auto const & column : columns )
  {
    addColumn( column );
  }
}

void TableLayout::addColumn( string_view columnName )
{
  TableLayout::Column column = TableLayout::Column().setName( columnName );
  m_tableColumns.emplace_back( column );
}

void TableLayout::addColumn( TableLayout::Column const & column )
{
  m_tableColumns.emplace_back( column );
}

TableLayout & TableLayout::setTitle( string_view title )
{
  m_tableTitleStr = title;
  return *this;
}

TableLayout & TableLayout::enableLineBreak( bool value )
{
  m_lineBreakAtBegin = value;
  return *this;
}

TableLayout & TableLayout::setMargin( MarginValue marginValue )
{
  m_marginValue = marginValue;
  m_borderMargin = marginValue;
  m_columnMargin = integer( marginValue ) * 2 + 1;

  return *this;
}

TableLayout & TableLayout::setMaxColumnWidth( size_t width )
{
  m_maxColumnWidth = width;
  return *this;
}

bool TableLayout::isLineBreakEnabled() const
{ return m_lineBreakAtBegin; }

TableLayout::CellLayout::CellLayout():
  m_cellWidth( 0 ),
  m_cellType( CellType::Header ),
  m_alignment( TableLayout::Alignment::center )
{}

TableLayout::CellLayout::CellLayout( CellType const cellType ):
  m_cellWidth( 0 ),
  m_cellType( cellType ),
  m_alignment( TableLayout::Alignment::center )
{}

TableLayout::CellLayout::CellLayout( CellType type, TableLayout::Alignment alignment ):
  m_cellWidth( 0 ),
  m_cellType( type ),
  m_alignment( alignment )
{}

TableLayout::Column::Column():
  m_headerStr(),
  m_headerLayout( CellType::Header )
{}

TableLayout::Column::Column( string_view name, TableLayout::ColumnAlignement alignment ):
  m_headerStr( name ),
  m_headerLayout( CellType::Header, alignment.headerAlignment ),
  m_alignment( alignment )
{}

TableLayout::Column & TableLayout::Column::setName( string_view name )
{
  m_headerStr = name;
  return *this;
}

TableLayout::Column & TableLayout::Column::setVisibility( bool visible )
{
  m_headerLayout.m_cellType = visible ? CellType::Header : CellType::Hidden;
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
    subColumns.emplace_back( TableLayout::Column( name, alignment ) );
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

TableLayout::Column & TableLayout::Column::addSubColumn( string_view subColName )
{
  m_subColumns.emplace_back( Column( subColName, TableLayout::ColumnAlignement() ) );
  return *this;
}

TableLayout::Column & TableLayout::Column::addSubColumn( TableLayout::Column const & subCol )
{
  m_subColumns.push_back( subCol );
  return *this;
}

TableLayout::Column & TableLayout::Column::setHeaderAlignment( Alignment headerAlignment )
{
  m_alignment.headerAlignment = headerAlignment;
  m_headerLayout.m_alignment = headerAlignment;

  if( !m_subColumns.empty() )
  {
    for( auto & subColumn : m_subColumns )
    {
      subColumn.setHeaderAlignment( headerAlignment );
    }
  }
  return *this;
}

TableLayout::Column & TableLayout::Column::setValuesAlignment( Alignment valueAlignment )
{
  m_alignment.valueAlignment = valueAlignment;

  if( !m_subColumns.empty() )
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
      m_currentColumn = &m_currentColumn->m_subColumns.front();
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

TableLayout::DeepFirstIterator TableLayout::beginDeepFirst() const
{
  TableLayout::Column const * startColumn = &(*m_tableColumns.begin());
  size_t idxLayer = 0;
  if( startColumn->hasChild() )
  {
    while( startColumn->hasChild() )
    {
      idxLayer++;
      startColumn = &startColumn->m_subColumns.front();
    }
  }
  return DeepFirstIterator( startColumn, idxLayer );
}

PreparedTableLayout::PreparedTableLayout(  ):
  TableLayout(),
  m_columnLayersCount( 0 ),
  m_visibleLowermostColumnCount( 0 )
{}

PreparedTableLayout::PreparedTableLayout( TableLayout const & other ):
  TableLayout( other ),
  m_columnLayersCount( 0 ),
  m_totalLowermostColumnCount( 0 ),
  m_visibleLowermostColumnCount( 0 )
{
  prepareLayoutRecusive( m_tableColumns, 0 );

  m_tableTitleLayout.prepareLayout( m_tableTitleStr, noColumnMaxWidth );
}

void PreparedTableLayout::prepareLayoutRecusive( std::vector< TableLayout::Column > & columns,
                                                 size_t level )
{

  for( size_t idxColumn = 0; idxColumn < columns.size(); ++idxColumn )
  {

    if( columns[idxColumn].isVisible())
      m_columnLayersCount = std::max( m_columnLayersCount, level + 1 );

    Column & column = columns[idxColumn];
    CellType cellType = column.m_headerLayout.m_cellType;

    if( !column.hasChild() )
    {
      ++m_totalLowermostColumnCount;
      if( cellType!= CellType::Hidden )
      {
        ++m_visibleLowermostColumnCount;
      }
    }

    column.m_headerLayout.prepareLayout( column.m_headerStr, getMaxColumnWidth() );

    if( idxColumn < columns.size() - 1 )
    {
      column.setNext( &columns[idxColumn + 1] );
    }

    if( !column.m_subColumns.empty())
    {
      for( auto & subCol : column.m_subColumns )
      {
        subCol.setParent( &column );
        subCol.setVisibility( subCol.isVisible() && column.isVisible() );
      }

      prepareLayoutRecusive( column.m_subColumns, level + 1 );

    }
  }
}

void TableLayout::CellLayout::prepareLayout( string_view inputText, size_t maxLineWidth )
{
  m_lines = stringutilities::divideLines< string_view >( m_cellWidth, inputText );

  if( maxLineWidth != noColumnMaxWidth && m_cellWidth > maxLineWidth )
  {
    m_lines = stringutilities::wrapTextToMaxLength( m_lines, maxLineWidth );
    // maxLineWidth has been updated
    m_cellWidth = maxLineWidth;
  }
}

}
