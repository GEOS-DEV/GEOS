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

std::vector< std::vector< string > > & TableData::getTableDataRows()
{
  return m_rows;
}

void TableData2D::buildRows( std::vector< std::vector< string > > & tableRows )
{
  for( real64 const & rowValue : rows )
  {
    std::vector< string > values;
    values.push_back( GEOS_FMT( "{}", rowValue ) );
    for( real64 const & columnValue : columns )
    {
      std::pair< real64, real64 > id = std::pair< real64, real64 >( rowValue, columnValue );
      auto const dataIt = data.find( id );

      if( dataIt != data.end())
      {
        values.push_back( GEOS_FMT( "{}", dataIt->second ));
      }
    }
    tableRows.push_back( values );
  }
}

}
