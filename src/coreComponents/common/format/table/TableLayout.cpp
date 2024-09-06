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

void addColumn( std::vector< TableLayout::ColumnParam > const & columnsParam, std::vector< TableLayout::Column > & columns )
{
  for( const auto & columnParam : columnsParam )
  {
    std::cout << "columnParam -- " << columnParam.columnName << std::endl;
    bool subColum = false;
    std::vector< TableLayout::Column > subColumnsToAdd;
    TableLayout::Column columnToAdd = {columnParam, {}, "", {}};
    for( const auto & subColumnsName: columnParam.subColumns )
    {
      std::cout << "subColumnsParam -- " << subColumnsName <<   std::endl;
      subColum = true;
      TableLayout::Column subColumn= {subColumnsName, {}, "", {}};
      subColumnsToAdd.push_back( subColumn );
      columnToAdd = { columnParam, {}, "", subColumnsToAdd };
    }
    if( subColum )
    {
      columns.push_back( columnToAdd );
    }
    else
    {
      columns.push_back( { columnParam, {}, "", {}} );
    }
  }

  for( auto const & columnt : columns )
  {
    std::cout <<  "regarde après le addColumn  " << columnt.m_parameter.columnName << std::endl ;
  }

}

TableLayout::TableLayout( std::vector< string > const & columnNames, string const & title ):
  m_tableTitle( title )
{
  setMargin( MarginValue::medium );
  m_columns.reserve( columnNames.size() );
  for( auto const & name : columnNames )
  {
    m_columns.push_back( {TableLayout::ColumnParam{{name}, Alignment::right, true}, {}, "", {}} );
  }
}

TableLayout::TableLayout( std::vector< ColumnParam > const & columnParams, string const & title ):
  m_tableTitle( title )
{
  setMargin( MarginValue::medium );
  //m_columns.reserve( columnParams.size() );
  std::cout << "TEST ADD COLUMN " << std::endl;
  addColumn( columnParams, m_columns );

}


void TableLayout::setMargin( MarginValue marginValue )
{
  m_borderMargin = marginValue;
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
