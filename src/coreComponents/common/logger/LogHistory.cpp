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
 * @file LogHistory.cpp
 */

#include "LogHistory.hpp"
#include "common/format/EnumStrings.hpp"
#include "common/format/table/TableData.hpp"
#include "common/format/table/TableLayout.hpp"
#include "common/format/table/TableTypes.hpp"
#include <utility>

namespace geos
{

void LogHistory::notifyMsg( LogPart::Type logPartName, DiagnosticMsg const & msgType, integer threadCount )
{
  // TODO reduction before notify and set count as parameter
  NumMsg numMsg = {msgType.m_file, msgType.m_line, threadCount};
  m_messageCounts
    .get_inserted( logPartName )
    .get_inserted( msgType.m_type )
    .emplace_back( numMsg );
}

template<>
string TableTextFormatter::toString< LogHistory >( LogHistory const & messageCounts ) const
{
  TableLayout tableLayoutPerSection;
  tableLayoutPerSection.addColumn( "Types" );

  for( size_t msgTypeIdx = (size_t) MsgType::Error; msgTypeIdx != (size_t)MsgType::Undefined; msgTypeIdx++ )
  {
    tableLayoutPerSection.addColumn( EnumStrings< MsgType >::toString( (MsgType) msgTypeIdx ) );
  }

  TableData data;
  for( auto const & [ logPartName, msgTypesStatistics ] : messageCounts.get() )
  {
    stdVector< TableData::CellData > row;
    row.emplace_back( TableData::CellData{CellType::Value,
                                          EnumStrings< LogPart::Type >::toString( logPartName )} );
    for( size_t msgTypeIdx = (size_t) MsgType::Error; msgTypeIdx != (size_t)MsgType::Undefined; msgTypeIdx++ )
    {
      MsgType const currentType = (MsgType) msgTypeIdx;
      int count = 0;
      auto itStatistics =  msgTypesStatistics.find( currentType );
      if( itStatistics != msgTypesStatistics.end())
      {
        for( auto const & stats : itStatistics->second )
        {
          count += stats.count;
        }
      }
      row.emplace_back( TableData::CellData{CellType::Value, std::to_string( count )} );
    }
    data.addRow( row );
  }

  TableTextFormatter textFormatter( tableLayoutPerSection );
  return textFormatter.toString( data ) + "\n";;
}


};
