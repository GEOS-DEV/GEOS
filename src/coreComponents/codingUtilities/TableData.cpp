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

// void TableData2D::printCount()
// {
//   int nb=0, nbx=0, nby=0;

//   for( real64 x : columns )
//   {
//     nbx++;
//   }

//   for( real64 y : row )
//   {
//     nby++;
//     for( real64 x : columns )
//     {
//       auto const dataIt = data.find( std::pair< real64, real64 >( x, y ));
//       if( dataIt != data.end())
//       {
//         nb++;
//       }
//     }
//   }

//   std::cout<<"total = "<<nb<<"  ;  "<<"totalx = "<<nbx<<"  ;  "<<"totaly = "<<nby<<std::endl;
// }

void TableData2D::buildRows()
{
  for( real64 y : row )
  {
    std::vector< string > values;
    values.push_back( GEOS_FMT( "{}", y ) );

    for( real64 x : columns )
    {
      std::pair< real64, real64 > id = std::pair< real64, real64 >( x, y );
      auto const dataIt = data.find( id );
      if( dataIt != data.end())
      {
        values.push_back( GEOS_FMT( "{}", dataIt->second ));
      }
    }

    m_rows.push_back( values );
  }
}
}
