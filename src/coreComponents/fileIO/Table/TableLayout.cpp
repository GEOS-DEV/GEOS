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
 * @file TableData.hpp
 */

#include "TableLayout.hpp"

namespace geos
{

TableLayout::TableLayout( std::vector< string > const & headers, string const & title ):
  m_tableTitle( title )
{
  setMargin( MarginValue::medium );
  for( size_t idx = 0; idx< headers.size(); ++idx )
  {
    m_columns.push_back( {TableLayout::ColumnParam{{headers[idx]}, Alignment::right, true}, {}, ""} );
  }
}

TableLayout::TableLayout( std::vector< ColumnParam > const & columnParameter, string const & title ):
  m_tableTitle( title )
{
  setMargin( MarginValue::medium );
  for( size_t idx = 0; idx< columnParameter.size(); ++idx )
  {
      m_columns.push_back( {columnParameter[idx], {}, ""} );
  }
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
  return m_marginTitle;
}

}
