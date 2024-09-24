/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
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

namespace geos
{

void TableLayout::addToColumns( const std::vector< string > & columnNames )
{
  for( const auto & columnName : columnNames )
  {
    addToColumns( string( columnName ) );
  }
}

void TableLayout::addToColumns( string_view columnName )
{
  m_columns.push_back( TableLayout::ColumnParam{ columnName, getDefaultAlignment() } );
}

void TableLayout::addToColumns( ColumnParam const & columnParam )
{
  if( !columnParam.subColumns.empty())
  {
    std::vector< TableLayout::Column > subColumns;
    for( const auto & subColumnsName : columnParam.subColumns )
    {
      subColumns.push_back( TableLayout::Column{ TableLayout::ColumnParam{subColumnsName, columnParam.alignment}  } );
    }
    m_columns.push_back( TableLayout::Column{ columnParam, subColumns } );
  }
  else
  {
    m_columns.push_back( TableLayout::Column{ columnParam } );
  }
}

TableLayout & TableLayout::setTitle( string_view title )
{
  m_tableTitle = title;
  return *this;
}

TableLayout & TableLayout::disableLineWrap()
{
  m_wrapLine = false;
  return *this;
}

TableLayout & TableLayout::setMargin( MarginValue marginValue )
{
  m_marginValue = marginValue;
  m_borderMargin = marginValue + 1; // margin + border character
  m_columnMargin = integer( marginValue ) * 2 + 1;

  return *this;
}

TableLayout & TableLayout::setAlignment( TableLayout::Alignment alignment )
{
  m_defaultAlignment = alignment;
  return *this;
}

bool TableLayout::isLineWrapEnabled() const
{
  return m_wrapLine;
}

void TableLayout::removeSubColumn()
{
  for( auto & column : m_columns )
  {
    if( !column.subColumn.empty() )
    {
      column.subColumn.clear();
    }
  }
}

std::vector< TableLayout::Column > const & TableLayout::getColumns() const
{
  return m_columns;
}

string_view TableLayout::getTitle() const
{
  return m_tableTitle;
}

integer const & TableLayout::getBorderMargin() const
{
  return m_borderMargin;
}

integer const & TableLayout::getColumnMargin() const
{
  return m_columnMargin;
}

integer const & TableLayout::getMarginValue() const
{
  return m_marginValue;
}

integer const & TableLayout::getMarginTitle() const
{
  return m_titleMargin;
}

TableLayout::Alignment TableLayout::getDefaultAlignment() const
{
  return m_defaultAlignment;
}

}
