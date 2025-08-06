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

void TableLayout::addColumns( stdVector< string > const & columnNames )
{
  for( auto const & columnName : columnNames )
  {
    addColumn( columnName );
  }
}

void TableLayout::addColumns( stdVector< TableLayout::Column > const & columns )
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
  m_cellType( CellType::Header ),
  m_alignment( TableLayout::Alignment::center ),
  m_cellWidth( 0 )
{}

TableLayout::CellLayout::CellLayout( CellType const cellType ):
  m_cellType( cellType ),
  m_alignment( TableLayout::Alignment::center ),
  m_cellWidth( 0 )
{}

TableLayout::CellLayout::CellLayout( CellType type, TableLayout::Alignment alignment ):
  m_cellType( type ),
  m_alignment( alignment ),
  m_cellWidth( 0 )
{}

TableLayout::Cell::Cell():
  m_layout(),
  m_text()
{}

TableLayout::Cell::Cell( CellType cellType, TableLayout::Alignment alignment ):
  m_layout( cellType, alignment ),
  m_text()
{}

TableLayout::Cell::Cell( CellType cellType, TableLayout::Alignment alignment, string_view value ):
  m_layout( cellType, alignment ),
  m_text( value )
{}

TableLayout::Cell::Cell( TableLayout::Cell const & other ):
  Cell()
{ *this = other; }

TableLayout::Cell::Cell( TableLayout::Cell && other ):
  Cell()
{ *this = other; }

TableLayout::Cell & TableLayout::Cell::operator=( TableLayout::Cell const & other )
{
  if( this != &other )
  {
    if( !m_layout.getLines().empty() || !other.m_layout.getLines().empty() )
      throw std::runtime_error( "Cannot copy from a Cell after its layout has been prepared." );

    m_layout.m_cellType = other.m_layout.m_cellType;
    m_layout.m_alignment = other.m_layout.m_alignment;
    m_text = other.m_text;
  }
  return *this;
}

TableLayout::Cell & TableLayout::Cell::operator=( TableLayout::Cell && other )
{
  if( this != &other )
  {
    if( !m_layout.getLines().empty() || !other.m_layout.getLines().empty() )
      throw std::runtime_error( "Cannot move from a Cell after its layout has been prepared." );

    m_layout.m_cellType = other.m_layout.m_cellType;
    m_layout.m_alignment = other.m_layout.m_alignment;
    m_text = std::move( other.m_text );
  }
  return *this;
}

void TableLayout::Cell::setText( string_view text )
{
  if( !m_layout.getLines().empty() )
    throw std::runtime_error( "Cannot reassign Cell text after its layout has been prepared." );

  m_text = text;
}

TableLayout::Column::Column():
  m_header( CellType::Header, defaultHeaderAlignment )
{}

TableLayout::Column::Column( string_view name, TableLayout::ColumnAlignement alignment ):
  m_header( CellType::Header, alignment.headerAlignment, name ),
  m_alignment( alignment )
{}

TableLayout::Column & TableLayout::Column::setName( string_view name )
{
  m_header.setText( name );
  return *this;
}

TableLayout::Column & TableLayout::Column::setVisibility( bool visible )
{
  m_header.m_layout.m_cellType = visible ? CellType::Header : CellType::Hidden;
  return *this;
}

TableLayout::Column & TableLayout::Column::addSubColumns( std::initializer_list< string > subColNames )
{
  m_subColumns.reserve( m_subColumns.size() + subColNames.size() );
  for( auto const & name : subColNames )
  {
    m_subColumns.emplace_back( TableLayout::Column( name, m_alignment ) );
  }
  return *this;
}

TableLayout::Column & TableLayout::Column::addSubColumns( stdVector< string > const & subColNames )
{
  m_subColumns.reserve( m_subColumns.size() + subColNames.size() );
  for( auto const & name : subColNames )
  {
    m_subColumns.emplace_back( TableLayout::Column( name, m_alignment ) );
  }
  return *this;
}

TableLayout::Column & TableLayout::Column::addSubColumns( std::initializer_list< TableLayout::Column > newSubColumns )
{
  m_subColumns.insert( m_subColumns.end(),
                       std::make_move_iterator( newSubColumns.begin() ),
                       std::make_move_iterator( newSubColumns.end() ) );
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
  m_header.m_layout.m_alignment = headerAlignment;

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

TableLayout::DeepFirstIterator TableLayout::beginDeepFirst() const
{
  TableLayout::Column const * startColumn = &(*m_tableColumns.begin());
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

PreparedTableLayout::PreparedTableLayout(  ):
  TableLayout(),
  m_columnLayersCount( 0 ),
  m_lowermostColumnCount( 0 )
{}

PreparedTableLayout::PreparedTableLayout( TableLayout const & other ):
  TableLayout( other ),
  m_columnLayersCount( 0 ),
  m_lowermostColumnCount( 0 )
{
  prepareLayoutRecusive( m_tableColumns, 0 );

  m_tableTitleLayout.prepareLayout( m_tableTitleStr, noColumnMaxWidth );
}

void PreparedTableLayout::prepareLayoutRecusive( stdVector< TableLayout::Column > & columns,
                                                 size_t level )
{
  for( size_t idxColumn = 0; idxColumn < columns.size(); ++idxColumn )
  {
    Column & column = columns[idxColumn];

    if( column.isVisible() )
    {
      m_columnLayersCount = std::max( m_columnLayersCount, level + 1 );

      if( !column.hasChild() )
      {
        ++m_lowermostColumnCount;
      }
    }

    column.m_header.prepareLayout( getMaxColumnWidth() );

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

void TableLayout::Cell::prepareLayout( size_t maxLineWidth )
{
  m_layout.prepareLayout( m_text, maxLineWidth );
}

}
