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

bool LogHistory::LogRecord::Key::operator==( Key const & rhs ) const
{
  return this->m_filename == rhs.m_filename  &&
         this->m_lineId == rhs.m_lineId;
}

bool LogHistory::LocationKeyHash::operator()( LogRecord::Key const & key ) const
{
  size_t h1 = std::hash< string >{} (key.m_filename);
  size_t h2 = std::hash< integer >{} (key.m_lineId);
  return h1 ^ (h2 << 1);
}

void LogHistory::LogRecord::deserialize( buffer_unit_type const * & logRecordBytes, buffer_unit_type const * end )
{
  basicSerialization::deserializeString( m_key.m_filename, logRecordBytes, end );
  basicSerialization::deserializePrimitive( m_key.m_lineId, logRecordBytes, end );
  basicSerialization::deserializeString( m_value.m_logPart, logRecordBytes, end );
  basicSerialization::deserializePrimitive( m_value.m_msgType, logRecordBytes, end );

  m_value.m_count = 0;
}

void LogHistory::LogRecord::serialize( stdVector< buffer_unit_type > & out ) const
{
  basicSerialization::serializeString( m_key.m_filename, out );
  basicSerialization::serializePrimitive( m_key.m_lineId, out );
  basicSerialization::serializeString( m_value.m_logPart, out );
  basicSerialization::serializePrimitive( m_value.m_msgType, out );
}

void LogHistory::recordDiagnostic( DiagnosticMsg const & msgType )
{
  string_view fileName =  extractAfterLastOccurrence( msgType.m_file, '/' );
  integer lineCount = msgType.m_line;
  insertDiagnosticReport( {
      /*.m_key = */ {
        /* .m_filename = */ string( fileName ),
        /* .m_lineId = */ lineCount
      },
      /*.m_value =*/ {
        /* .m_logPart = */ string( msgType.m_logPart ),
        /* .m_msgType = */ msgType.m_type,
        /* .m_count = */ 1
      }
    } );
}

void LogHistory::insertDiagnosticReport( LogRecord const & logRecord )
{
  auto it =  m_diagnosticHistory.find( {logRecord.m_key} );
  if( it == m_diagnosticHistory.end())
  {
    m_diagnosticHistory.emplace( logRecord.m_key, logRecord.m_value );
  }
  else
  {
    it->second.m_count +=  logRecord.m_value.m_count;
  }
}



size_t LogHistory::LogRecord::getSerializedSize() const
{
  return
    basicSerialization::sizeOfString( m_key.m_filename ) +
    basicSerialization::sizeOfPrimitive( m_key.m_lineId ) +
    basicSerialization::sizeOfString( m_value.m_logPart ) +
    basicSerialization::sizeOfPrimitive( m_value.m_msgType );
}

void LogHistory::gatherRecordsRank0()
{
  stdVector< buffer_unit_type > localLogRecords( 0 );
  integer totalSize = 0;

  { // allocation
    for( auto const & [key, value] : getDiagnosticHistory() )
    {
      LogRecord record( key, value );
      totalSize +=  record.getSerializedSize();
    }
    localLogRecords.reserve( totalSize );
  }


  { // Packing
    if( getDiagnosticHistory().size() > 0 )
    {
      for( auto const & [key, value] :  getDiagnosticHistory() )
      {
        LogRecord record( key, value );
        record.serialize( localLogRecords );
      }
    }
  }

  auto [globalLogRecords, counts, offsets] =
    MpiWrapper::gatherBufferRank0< stdVector< buffer_unit_type > >( localLogRecords );

  { // Unpacking
    if( MpiWrapper::commRank() == 0 )
    {
      buffer_unit_type const * startGlobalRecord = globalLogRecords.data();
      for( size_t idxRank = 0; idxRank <  (size_t)MpiWrapper::commSize(); ++idxRank )
      {
        integer byteFromThisRank = counts[idxRank];
        buffer_unit_type const * rankEnd= startGlobalRecord + byteFromThisRank;
        while( startGlobalRecord < rankEnd )
        {
          LogRecord unpackRecord;
          unpackRecord.deserialize( startGlobalRecord, rankEnd );
          insertDiagnosticReport( unpackRecord );
        }
      }
    }
  }
}



template<>
string TableTextFormatter::toString< LogHistory >( LogHistory const & logHistory ) const
{
  using CellRow  = stdArray< TableData::CellData, (size_t) MsgType::Count >;

  TableLayout tableLayout;
  tableLayout.addColumn( "Types" );

  // fill header
  for( size_t msgTypeIdx = 0; msgTypeIdx != (size_t)MsgType::Count; msgTypeIdx++ )
  {
    tableLayout.addColumn( EnumStrings< MsgType >::toString( (MsgType) msgTypeIdx ) );
  }

  // prepare logHistory row
  stdMap< std::pair< string, MsgType >, integer > countPerPartAndType;
  CellRow emptyCellRow;
  emptyCellRow.fill( TableData::CellData{CellType::Value, "0"} );

  stdMap< string, CellRow > rowByPart;
  // retrive unique message
  for( const auto & [key, values] : logHistory.getDiagnosticHistory())
  {
    auto logPart = values.m_logPart;
    MsgType msgType = values.m_msgType;
    countPerPartAndType.get_inserted( std::make_pair( logPart, msgType ))++;

    if( rowByPart.find( logPart ) == rowByPart.end())
      rowByPart.get_inserted( logPart ) = emptyCellRow;
  }

  // update rowByPart values
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
    data.addRow( row );
  }

  TableTextFormatter textFormatter( tableLayout );
  return textFormatter.toString( data ) + "\n";
}


}
