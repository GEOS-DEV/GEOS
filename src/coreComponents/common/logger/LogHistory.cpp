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
  string fileName =  extractAfterLastOccurrence( msgType.m_file, '/' );
  integer lineCount = msgType.m_line;

  auto & stats = m_errorHistory.get_inserted( std::make_tuple( logPartName, msgType.m_type, fileName, lineCount ));
  stats.count++;
}

template< typename T >
stdVector< T >  gatherAndUnpackRank0( buffer_unit_type * bufferToSend, localIndex packedSize, size_t maxIteration )
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

  MpiWrapper::gatherv( bufferToSend,
                       numValues,
                       globalAllocations.data(),
                       recvCounts.data(),
                       displs.data(),
                       0 );
  if( MpiWrapper::commRank() == 0 )
  {
    stdVector< T > targettedBuf( maxIteration);
    buffer_unit_type const * globalAllocBuffer = globalAllocations.data();

    for( size_t i = 0; i <  maxIteration; ++i )
    {
      if( recvCounts[i] > 0 )
      {
        globalAllocBuffer =  globalAllocations.data();
        bufferOps::Unpack(globalAllocBuffer, targettedBuf[i] );
      }
    }
    return targettedBuf;
  }
  else
  {
    return {};
  }

}

void LogHistory::diagnosticStatsReport()
{
  LogHistory & history = ErrorLogger::global().getLoggerReportData();
  auto maxIteration = MpiWrapper::max( getDiagnosticHistory().size());
  auto it = getDiagnosticHistory().begin();
  for( size_t i = 0; i< maxIteration; i++ )
  {
    localIndex packedSizeLogPart =  0;
    localIndex packedSizeMsgType =  0;
    localIndex packedSizeFilename = 0;
    localIndex packedSizeLineCount =0;

    stdVector< buffer_unit_type > gLogPart( 0 );
    stdVector< buffer_unit_type > gMsgType( 0 );
    stdVector< buffer_unit_type > gFileName( 0 );
    stdVector< buffer_unit_type > gLineCount( 0 );

    if( it != getDiagnosticHistory().end())
    {
      auto key = it->first;
      auto const & [logPartType, msgType, filename, lineCount] = key;

      buffer_unit_type * ptr = nullptr;
      packedSizeLogPart =   bufferOps::Pack< false >( ptr, logPartType );
      packedSizeMsgType =   bufferOps::Pack< false >( ptr, msgType );
      packedSizeFilename =  bufferOps::Pack< false >( ptr, filename );
      packedSizeLineCount = bufferOps::Pack< false >( ptr, lineCount );

      gLogPart.resize( packedSizeLogPart );
      gMsgType.resize( packedSizeMsgType );
      gFileName.resize( packedSizeFilename );
      gLineCount.resize( packedSizeLineCount );

      buffer_unit_type * itLogPart = gLogPart.data();
      buffer_unit_type * itMsgType = gMsgType.data();
      buffer_unit_type * itFilename = gFileName.data();
      buffer_unit_type * itLineCount = gLineCount.data();

      bufferOps::Pack< true >( itLogPart, logPartType );
      bufferOps::Pack< true >( itMsgType, msgType );
      bufferOps::Pack< true >( itFilename, filename );
      bufferOps::Pack< true >( itLineCount, lineCount );
      it++;

    }

    stdVector< LogPart::Type > logsParts = gatherAndUnpackRank0< LogPart::Type >( gLogPart.data(), packedSizeLogPart,maxIteration );
    stdVector< MsgType > msgsTypes = gatherAndUnpackRank0< MsgType >( gMsgType.data(), packedSizeMsgType,maxIteration );
    stdVector< string > filenames = gatherAndUnpackRank0< string >( gFileName.data(), packedSizeFilename,maxIteration );
    stdVector< integer > linesNumber = gatherAndUnpackRank0< integer >( gLineCount.data(), packedSizeLineCount,maxIteration );

    if( MpiWrapper::commRank() == 0 )
    {
      for( size_t indexIteration = 0; indexIteration<maxIteration; indexIteration++ )
      {
        if( linesNumber[indexIteration] != 0 )
          history.insertBlanckReport( logsParts[indexIteration], msgsTypes[indexIteration],
                                      filenames[indexIteration], linesNumber[indexIteration] );
      }
    }
  }

  if( MpiWrapper::commRank() == 0 )
  {
    TableTextFormatter tableReportFormatter;
    GEOS_LOG( tableReportFormatter.toString< LogHistory >( GEOS_GLOBAL_LOGGER.getLoggerReportData()));
  }
}

void LogHistory::insertBlanckReport( LogPart::Type logPartName, MsgType msgType, string const & fileName, integer lineCount )
{
  m_errorHistory.get_inserted( std::make_tuple( logPartName, msgType, fileName, lineCount ));
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


  for( const auto & [tupleKey, msgTypes] : messageCounts.getDiagnosticHistory())
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
