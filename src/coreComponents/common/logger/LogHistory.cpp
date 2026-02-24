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
#include "common/StdContainerWrappers.hpp"
#include "common/format/EnumStrings.hpp"
#include "common/format/LogPart.hpp"
#include "common/format/table/TableData.hpp"
#include "common/format/table/TableLayout.hpp"
#include "common/format/table/TableTypes.hpp"
#include "common/logger/MsgType.hpp"
#include "common/MpiWrapper.hpp"
#include "dataRepository/BufferOps.hpp"
#include "dataRepository/BufferOps_inline.hpp"
#include <algorithm>
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

  auto & stats = m_errorHistory.get_inserted( std::make_tuple( logPartName, msgType.m_type, locationKey ));
  stats.locationKey = locationKey;
  stats.count++;
}

template< typename T >
stdVector< T >  gatherAndUnPackArg( buffer_unit_type * localBuffer, localIndex packedSize )
{
  integer const numRanks = MpiWrapper::commSize();
  integer const numValues = packedSize;

  // Allows to know how much data each rank will send
  stdVector< integer > recvCounts;
  // Displacments vector for global alloc
  stdVector< integer > displs;
  stdVector< buffer_unit_type > globalAllocations;

  if( MpiWrapper::commRank() == 0 )
  {
    recvCounts.resize( numRanks );
    displs.resize( numRanks );
  }


  MpiWrapper::gather( &numValues, 1, recvCounts.data(), 1, 0 );

  if( MpiWrapper::commRank() == 0 )
  {
    integer totalSize = 0;
    for( integer i = 0; i < numRanks; ++i )
    {
      displs[i] = totalSize;
      totalSize += recvCounts[i];
    }
    globalAllocations.resize( totalSize );
  }

  MpiWrapper::gatherv( localBuffer,
                       numValues,
                       globalAllocations.data(),
                       recvCounts.data(),
                       displs.data(),
                       0 );
  if( MpiWrapper::commRank() == 0 )
  {
    stdVector< T > m_receiveBufferPtr;
    //TODO CHANGE IT
    m_receiveBufferPtr.resize( 4 );
    buffer_unit_type const * localDestBuffer = globalAllocations.data();

    int const totalIntegers = globalAllocations.size() / sizeof( T );
    m_receiveBufferPtr.resize( totalIntegers );

    for( int i = 0; i < totalIntegers; ++i )
    {
      bufferOps::Unpack( localDestBuffer, m_receiveBufferPtr[i] );
    }
    return m_receiveBufferPtr;
  }
  else
  {
    return {};
  }

}

void LogHistory::errorStatsReport()
{

  std::shared_ptr< buffer_unit_type > localBufferLogPart = std::make_unique< buffer_unit_type >();
  signed char *  localBufferLogPartIter = localBufferLogPart.get();
  localIndex packedSizeLogPart = 0;

  std::shared_ptr< buffer_unit_type > localBufferMsgType = std::make_unique< buffer_unit_type >();
  signed char * localBufferMsgTypeIter = localBufferMsgType.get();
  localIndex packedSizeMsgType = 0;

  std::shared_ptr< buffer_unit_type > localBufferPair = std::make_unique< buffer_unit_type >();
  signed char * localBufferPairIter = localBufferPair.get();
  localIndex packedSizePair = 0;
  for( auto const & [key, stats] : getErrorHistory() )//todo renommer
  {
    auto const & [logPartType, msgType, locationKey] = key;

    packedSizeLogPart += bufferOps::Pack< true >( localBufferLogPartIter, logPartType );
    packedSizeMsgType += bufferOps::Pack< true >( localBufferMsgTypeIter, msgType );
    packedSizePair += bufferOps::Pack< true >( localBufferPairIter, locationKey );

  }

  stdVector< LogPart::Type > testLogPart = gatherAndUnPackArg< LogPart::Type >( localBufferLogPart.get(), packedSizeLogPart );
  stdVector< MsgType > testMsgType = gatherAndUnPackArg< MsgType >( localBufferMsgType.get(), packedSizeMsgType );
  stdVector< MsgStatistics::LocationKey > testPair = gatherAndUnPackArg< MsgStatistics::LocationKey >( localBufferPair.get(), packedSizePair );
  if( MpiWrapper::commRank() == 0 )
  {
    LogHistory & history = ErrorLogger::global().getLoggerReportData();
    for( size_t i = 0; i<testPair.size(); i++ )
    {
      history.insertBlanckReport( testLogPart[i], testMsgType[i], testPair[i] );
    }

    TableTextFormatter tableReportFormatter;
    GEOS_LOG( tableReportFormatter.toString< LogHistory >( GEOS_GLOBAL_LOGGER.getLoggerReportData()));
  }
}

void LogHistory::insertBlanckReport( LogPart::Type logPartName, MsgType msgType, MsgStatistics::LocationKey locationKey )
{
  m_errorHistory.get_inserted( std::make_tuple( logPartName, msgType, locationKey ));
}

template<>
string TableTextFormatter::toString< LogHistory >( LogHistory const & messageCounts ) const
{
  TableLayout tableLayout;
  tableLayout.addColumn( "Types" );

  for( size_t msgTypeIdx = (size_t) MsgType::Error; msgTypeIdx != (size_t)MsgType::Undefined; msgTypeIdx++ )
  {
    tableLayout.addColumn( EnumStrings< MsgType >::toString( (MsgType) msgTypeIdx ) );
  }

  stdMap< std::pair< LogPart::Type, MsgType >, int > countPerPartAndType;
  using CellRow  = stdArray< TableData::CellData, (size_t) MsgType::Undefined >;
  CellRow emptyCellRow;
  emptyCellRow.fill( TableData::CellData{CellType::Value, "0"} );
  stdMap< LogPart::Type, CellRow > rowByPart;


  for( const auto & [tupleKey, msgTypes] : messageCounts.getErrorHistory())
  {
    auto part = std::get< 0 >( tupleKey );
    auto type = std::get< 1 >( tupleKey );
    
    countPerPartAndType.get_inserted( std::make_pair( part, type ))++;

    if( rowByPart.find( part ) == rowByPart.end())
      rowByPart.get_inserted( part ) = emptyCellRow;
  }

  for( auto & [keyPair, count] : countPerPartAndType )
  {
    auto part = std::get< 0 >( keyPair );
    auto type = std::get< 1 >( keyPair );
    rowByPart.get_inserted( part ).at((size_t)type ).value =  std::to_string( count );
  }

  TableData data;
  for( auto const &  [key, cells] : rowByPart )
  {
    stdVector< TableData::CellData >row (
      {
        TableData::CellData{CellType::Value,
                            EnumStrings< LogPart::Type >::toString( key )}
      } );

    row.insert( row.end(), cells.begin(), cells.end());
    data.addRow( row );
  }

  TableTextFormatter textFormatter( tableLayout );
  return textFormatter.toString( data ) + "\n";
}


};
