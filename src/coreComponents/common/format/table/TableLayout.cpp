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
 * @file TableData.hpp
 */
#include "TableLayout.hpp"

namespace geos
{

void addColumns( std::vector< TableLayout::ColumnParam > const & columnsParam, std::vector< TableLayout::Column > & columns )
{
  for( const auto & columnParam : columnsParam )
  {
    if( !columnParam.subColumns.empty() )
    {
      std::vector< TableLayout::Column > subColumns;
      for( const auto & subColumnsName: columnParam.subColumns )
      {
        subColumns.push_back( TableLayout::Column{subColumnsName} );
      }
      columns.push_back( TableLayout::Column{columnParam, subColumns} );
    }
    else
    {
      columns.push_back( TableLayout::Column{columnParam} );
    }
  }
}

TableLayout::TableLayout( std::vector< string > const & columnNames )
{
  setMargin( MarginValue::medium );

  std::vector< TableLayout::ColumnParam > columnParams;
  m_columns.reserve( columnNames.size() );
  columnParams.reserve( columnNames.size() );

  for( auto const & name : columnNames )
  {
    columnParams.push_back( TableLayout::ColumnParam{{name}} );
  }

  addColumns( columnParams, m_columns );
}

void TableLayout::setTitle( std::string const & title )
{
  m_tableTitle = title;
}

void TableLayout::noWrapLine()
{
  m_wrapLine = false;
}

bool TableLayout::isLineWrapped() const
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

void TableLayout::setMargin( MarginValue marginValue )
{
  m_borderMargin = marginValue + 1;
  m_columnMargin = integer( marginValue ) * 2 + 1;
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
