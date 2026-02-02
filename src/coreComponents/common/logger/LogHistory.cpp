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
#include "common/format/LogPart.hpp"
#include "common/format/table/TableData.hpp"
#include "common/format/table/TableLayout.hpp"
#include "common/format/table/TableTypes.hpp"
#include "common/logger/MsgType.hpp"
#include <functional>
#include <string>
#include <utility>

namespace geos
{

std::string extractAfterLastOccurrence( string const & str, char delimiter )
{
  size_t pos = str.find_last_of( delimiter );

  if( pos == std::string::npos )
  {
    return "";
  }

  return str.substr( pos + 1 );
}

void LogHistory::notifyMsg( LogPart::Type logPartName, DiagnosticMsg const & msgType )
{
  MsgStatistics::LocationKey locationKey;
  std::strcpy( locationKey.first.data(), extractAfterLastOccurrence( msgType.m_file, '/' ).c_str() );
  locationKey.second = msgType.m_line;

  auto & stats = m_messageCounts.get_inserted( std::make_tuple( logPartName, msgType.m_type, locationKey ));
  stats.locationKey = locationKey;
  stats.count++;
}

void LogHistory::insertBlanckReport( LogPart::Type logPartName, MsgType msgType, MsgStatistics::LocationKey locationKey )
{
  m_messageCounts.get_inserted( std::make_tuple( logPartName, msgType, locationKey ));
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
  stdMap< LogPart::Type, stdMap< MsgType, int > > countsByLogPart;
  for( auto const & [ key, stats ] : messageCounts.getMessageCounts() )
  {

    auto [logPartType, msgType, locationKey] = key;
    countsByLogPart.get_inserted( logPartType ).get_inserted( msgType )++;
  }

  for( const auto & [logPart, msgTypes] : countsByLogPart )
  {
    stdVector< TableData::CellData > row;
    row.emplace_back( TableData::CellData{CellType::Value,
                                          EnumStrings< LogPart::Type >::toString( logPart )} );
    for( const auto & [msgType, count] : msgTypes )
    {
      row.emplace_back( TableData::CellData{CellType::Value, std::to_string( count )} );
    }
    data.addRow( row );
  }

  TableTextFormatter textFormatter( tableLayoutPerSection );
  return textFormatter.toString( data ) + "\n";;
}


};
