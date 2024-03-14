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

#include "codingUtilities/TableLayout.hpp"

namespace geos
{

TableLayout::TableLayout( std::vector< string > const & headers )
{
  setMargin( MarginValue::medium );
  for( size_t idx = 0; idx< headers.size(); idx++ )
  {
    m_columns.push_back( {TableLayout::ColumnParam{{headers[idx]}, Alignment::center, true}, {}, ""} );
  }
}

TableLayout::TableLayout( std::vector< ColumnParam > const & columnParameter )
{
  setMargin( MarginValue::medium );

  for( size_t idx = 0; idx< columnParameter.size(); idx++ )
  {
    if( columnParameter[idx].enabled )
    {
      m_columns.push_back( {columnParameter[idx], {}, ""} );
    }

  }
}

void TableLayout::setTitle( string_view title )
{
  tableTitle = title;
}

void TableLayout::setMargin( MarginValue marginType )
{
  borderMargin = marginType;
  columnMargin = integer( marginType ) * 2 + 1;
}

string_view TableLayout::getTitle() const
{
  return tableTitle;
}


integer const & TableLayout::getBorderMargin() const
{
  return borderMargin;
}

integer const & TableLayout::getColumnMargin() const
{
  return columnMargin;
}

integer const & TableLayout::getMarginTitle() const
{
  return marginTitle;
}

}
