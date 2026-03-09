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
#include "common/DataTypes.hpp"
#include "common/StdContainerWrappers.hpp"
#include "common/format/EnumStrings.hpp"
#include "common/format/table/TableData.hpp"
#include "common/format/table/TableLayout.hpp"
#include "common/format/table/TableTypes.hpp"
#include "common/MpiWrapper.hpp"
#include <cstring>
#include <string>
#include <utility>

namespace geos
{

string_view extractAfterLastOccurrence( string_view str, char delimiter )
{
  size_t pos = str.find_last_of( delimiter );

  if( pos == std::string::npos )
  {
    return str;
  }

  return str.substr( pos + 1 );
}

LogHistory::LogRecord::LogRecord():
  m_key( {} ),
  m_value( {} )
{}
LogHistory::LogRecord::LogRecord( Key const & key, Values const & values ):
  m_key( key ),
  m_value( values )
{}

void LogHistory::LogRecord::deserialize( stdVector< buffer_unit_type > const & logRecordBytes )
{
  buffer_unit_type const * start = logRecordBytes.data();
  buffer_unit_type const * end = logRecordBytes.end().base();
  while( start < end )
  {
    LogRecord unpackRecord;
    deserializeField( m_key.m_filename, start, end );
    deserializeField( m_key.m_lineId, start, end );
    deserializeField( m_value.m_logPart, start, end );
    deserializeField( m_value.m_msgType, start, end );
    m_value.m_count = 0;
  }
}

void LogHistory::LogRecord::serialize( stdVector< buffer_unit_type > & out ) const
{
  auto filenameSize = m_key.m_filename.size();
  auto logPartSize = m_value.m_logPart.size();

  auto const serializeField = [&]( void const * data, size_t size )
  {
    buffer_unit_type const * d = (buffer_unit_type const *) data;
    out.insert( out.end(), d, d + size );
  };
  serializeField( &filenameSize, sizeof(string::size_type));
  serializeField( m_key.m_filename.data(), m_key.m_filename.size());
  serializeField( &m_key.m_lineId, sizeof(integer));
  serializeField( &logPartSize, sizeof(string::size_type));
  serializeField( m_value.m_logPart.data(), m_value.m_logPart.size());
  serializeField( &m_value.m_msgType, sizeof(MsgType) );
}

void LogHistory::insertDiagnosticReport( LogRecord logRecord )
{
  auto & entry =  m_diagnosticHistory.get_inserted( {logRecord.m_key} );
  entry.m_logPart = logRecord.m_value.m_logPart;
  entry.m_msgType = logRecord.m_value.m_msgType;
  entry.m_count +=  1;
}

void LogHistory::notifyMsg( string_view logPartName, DiagnosticMsg const & msgType )
{
  string_view fileName =  extractAfterLastOccurrence( msgType.m_file, '/' );
  integer lineCount = msgType.m_line;
  insertDiagnosticReport( { {string( fileName ), lineCount},
                            {string( logPartName ), msgType.m_type, 1} } );
}

template< typename T >
std::pair< stdVector< T >, stdVector< integer > >
gatherBufferRank0( stdVector< T > const & bufferToSend )
{
  integer const numRanks = MpiWrapper::commSize();
  integer const numValues = bufferToSend.size();

  // Allows to know how much data each rank will send
  stdVector< integer > recvCounts;
  // Displacments vector for global alloc
  stdVector< integer > displs;
  stdVector< T > globalAllocations;

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

  MpiWrapper::gatherv( bufferToSend.data(),
                       numValues,
                       globalAllocations.data(),
                       recvCounts.data(),
                       displs.data(),
                       0 );

  return {globalAllocations, recvCounts};
}

size_t LogHistory::LogRecord::getSerializedSize() const
{
  return
    sizeOfField( m_key.m_filename ) +
    sizeOfField( m_key.m_lineId ) +
    sizeOfField( m_value.m_logPart ) +
    sizeOfField( m_value.m_msgType );
}

void LogHistory::diagnosticStatsReport()
{
  LogHistory & history = ErrorLogger::global().getLoggerReportData();
  stdVector< buffer_unit_type > localLogRecords( 0 );
  integer totalSize = 0;

  //0 - dry run
  for( auto const & [key, value] : getDiagnosticHistory() )
  {
    LogRecord record( key, value );
    totalSize +=  record.getSerializedSize();
  }
  localLogRecords.reserve( totalSize );

  //1 - Packing
  if( getDiagnosticHistory().size() > 0 )
  {
    for( auto const & [key, value] :  getDiagnosticHistory() )
    {
      LogRecord record( key, value );
      record.serialize( localLogRecords );
    }
  }

  auto [globalLogRecords, recvCounts] = gatherBufferRank0( localLogRecords );

  //2 - Unpacking
  if( MpiWrapper::commRank() == 0 )
  {
    auto startGlobalRecord =  globalLogRecords.begin();
    for( size_t idxRank = 0; idxRank <  (size_t)MpiWrapper::commSize(); ++idxRank )
    {
      auto end =  globalLogRecords.begin();
      LogRecord unpackRecord;
      unpackRecord.deserialize( stdVector< buffer_unit_type >( startGlobalRecord, end ));
      history.insertDiagnosticReport( unpackRecord );
      startGlobalRecord += recvCounts[idxRank];
    }

  }

  //3 - Display
  if( MpiWrapper::commRank() == 0 )
  {
    TableTextFormatter tableReportFormatter;
    GEOS_LOG( tableReportFormatter.toString< LogHistory >( *this ));
  }
}

template<>
string TableTextFormatter::toString< LogHistory >( LogHistory const & messageCounts ) const
{
  TableLayout tableLayout;
  tableLayout.addColumn( "Types" );

  // fill header
  for( size_t msgTypeIdx = (size_t) MsgType::Error; msgTypeIdx != (size_t)MsgType::Undefined; msgTypeIdx++ )
  {
    tableLayout.addColumn( EnumStrings< MsgType >::toString( (MsgType) msgTypeIdx ) );
  }

  stdMap< std::pair< string, MsgType >, integer > countPerPartAndType;
  using CellRow  = stdArray< TableData::CellData, (size_t) MsgType::Undefined >;
  CellRow emptyCellRow;
  emptyCellRow.fill( TableData::CellData{CellType::Value, "0"} );

  stdMap< string, CellRow > rowByPart;

  // TODO
  for( const auto & [key, values] : messageCounts.getDiagnosticHistory())
  {
    auto logPart = values.m_logPart;
    MsgType msgType = values.m_msgType;

    countPerPartAndType.get_inserted( std::make_pair( logPart, msgType ))++;

    if( rowByPart.find( logPart ) == rowByPart.end())
      rowByPart.get_inserted( logPart ) = emptyCellRow;
  }

  for( auto & [keyPair, count] : countPerPartAndType )
  {
    auto logPart = std::get< 0 >( keyPair );
    auto msgType = std::get< 1 >( keyPair );
    rowByPart.get_inserted( logPart ).at((size_t)msgType ).value =  std::to_string( count );
  }

  TableData data;
  for( auto const &  [logPart, cells] : rowByPart )
  {
    stdVector< TableData::CellData >row ( {
        TableData::CellData{ CellType::Value, logPart }
      } );

    row.insert( row.end(), cells.begin(), cells.end());
    data.addRow< TableData::CellData >( row );
  }

  TableTextFormatter textFormatter( tableLayout );
  return textFormatter.toString( data ) + "\n";
}


}
