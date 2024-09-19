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

void TableLayout::addColumns( std::vector< TableLayout::ColumnParam > & columnsParam )
{
  for( auto && columnParam : columnsParam )
  {
    addToColumns( std::move( columnParam ) );
  }
}

void TableLayout::addToColumns( const std::vector< std::string > & columnNames )
{
  for( const auto & columnName : columnNames )
  {
    addToColumns( std::string( columnName ) );
  }
}

void TableLayout::addToColumns( std::string && columnName )
{
  m_columns.push_back( TableLayout::ColumnParam{ {columnName} } );
}

void TableLayout::addToColumns( ColumnParam && columnParam )
{
  if( !columnParam.subColumns.empty())
  {
    std::vector< TableLayout::Column > subColumns;
    for( const auto & subColumnsName : columnParam.subColumns )
    {
      subColumns.push_back( TableLayout::Column{ subColumnsName } );
    }
    m_columns.push_back( TableLayout::Column{ columnParam, subColumns } );
  }
  else
  {
    m_columns.push_back( TableLayout::Column{ columnParam } );
  }
}

TableLayout & TableLayout::setTitle( std::string const & title )
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
  m_borderMargin = marginValue + 1;
  m_columnMargin = integer( marginValue ) * 2 + 1;

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

integer const & TableLayout::getMarginTitle() const
{
  return m_titleMargin;
}

}
