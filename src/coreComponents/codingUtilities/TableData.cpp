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
 * @file TableData.cpp
 */

#include "codingUtilities/TableData.hpp"

namespace geos
{

void TableData::addRow( std::vector< string > const & row )
{
  m_rows.push_back( row );
}

void TableData::clear()
{
  m_rows.clear();
}

std::vector< std::vector< string > > const & TableData::getTableDataRows() const
{
  return m_rows;
}

std::set< real64 > const & TableData2D::getColumns() const
{
  return m_columns;
}
std::set< real64 > const & TableData2D::getRows() const
{
  return m_rows;
}

TableData TableData2D::buildTableData() const
{
  TableData tableDataToBeBuilt;
  for( real64 const & rowValue : m_rows )
  {
    std::vector< string > values;
    values.push_back( GEOS_FMT( "{}", rowValue ) );
    for( real64 const & columnValue : m_columns )
    {
      std::pair< real64, real64 > id = std::pair< real64, real64 >( rowValue, columnValue );
      auto const dataIt = m_data.find( id );

      if( dataIt != m_data.end())
      {
        values.push_back( GEOS_FMT( "{}", dataIt->second ));
      }

    }
    tableDataToBeBuilt.addRow( values );
  }
  return tableDataToBeBuilt;
}

}
