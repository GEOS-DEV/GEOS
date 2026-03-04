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
 * @file LogHistory.hpp
 */

#ifndef GEOS_COMMON_LOGGER_MSG_REPORT_DATA_HPP
#define GEOS_COMMON_LOGGER_MSG_REPORT_DATA_HPP

#include "common/StdContainerWrappers.hpp"
#include "common/format/LogPart.hpp"
#include "common/format/table/TableFormatter.hpp"
#include "DiagnosticMessage.hpp"
#include "common/logger/MsgType.hpp"
#include <string>


namespace geos
{

struct Archive
{
  std::vector< buffer_unit_type > buffer;
  size_t offset = 0;
  buffer_unit_type *  bufferIt = buffer.data();

  enum Mode { PACK, UNPACK };
  Mode mode;

  template< typename T >
  void process( T & value )
  {
    if constexpr (std::is_fundamental_v< T > || std::is_enum_v< T >) {
      if( mode == PACK )
      {
        memcpy( bufferIt, &value, sizeof(T) );
        bufferIt += sizeof(T);
      }
      else
      {
        std::memcpy( &value, &buffer[offset], sizeof(T));
        offset += sizeof(T);
      }
    }
  }

  void process( string & s )
  {
    if( mode == PACK )
    {
      size_t size = s.size();
      process( size );
      buffer.insert( buffer.end(), s.begin(), s.end());
    }
    else
    {
      size_t size;
      process( size );
      s = std::string( bufferIt, bufferIt + size );
      offset += size;
    }
  }
};

/**
 * @brief Data for a diagnostic message
 */
struct DiagnosticData
{
  /// Number of times the same message occured during the simulation
  integer m_count;
};


/**
 * @brief Keep track of all diagnostic message occured during the simulation
 */
class LogHistory
{
public:

  /// POD characterizing an unique diagnostic message
  struct DiagnosticKey
  {
    string logPart;
    MsgType msgType;
    string fileName;
    integer lineId;

    template< typename Archive >
    void serialize( Archive & ar )
    {
      ar.process( logPart );
      ar.process( msgType );
      ar.process( fileName );
      ar.process( lineId );
    }
  };

  /**
   * @brief Report a diagnostic message
   * @param logPartName The logPart where the message occured
   * @param diagMsg The DiagnosticMsg to record
   */
  void notifyMsg( string_view logPartName, DiagnosticMsg const & diagMsg );

  /**
   * @brief Display the diagnostic statistics to the log
   */
  void diagnosticStatsReport();

  /**
   * @return The const historical diagnostic
   */
  auto const & getDiagnosticHistory() const
  { return m_diagnosticHistory; }

private:

  /// @cond DO_NOT_DOCUMENT
  struct LocationKeyHash
  {

    size_t operator()( DiagnosticKey const & key ) const noexcept
    {
      auto const & [logPartType, msgType, filename, lineId] = key;

      size_t h1 = std::hash< string >{} (logPartType);
      size_t h2 = std::hash< MsgType >{} (msgType);
      size_t h3 = std::hash< string >{} (filename);
      size_t h4 = std::hash< integer >{} (lineId);

      return h1 ^ (h2 << 1) ^ (h3 << 2) ^ (h4 << 3);
    }

  };
  /// @endcond

  /// @cond DO_NOT_DOCUMENT
  struct LocationKeyEqual
  {

    bool operator()( DiagnosticKey const & lhs,
                     DiagnosticKey const & rhs ) const
    {
      return lhs.logPart == rhs.logPart &&
             lhs.msgType == rhs.msgType &&
             lhs.fileName == rhs.fileName &&
             lhs.lineId == rhs.lineId;
    }
  };
  /// @endcond

  /**
   * @brief Diagnostic history happened during the simulation
   */
  stdUnorderedMap< DiagnosticKey,
                   DiagnosticData,
                   LocationKeyHash, LocationKeyEqual > m_diagnosticHistory;

  /**
   * @brief Insert an element to the diagnostic history container if an equivalent key doesn't exist.
   * @param logPartName The logPart where the diagnostic occured
   * @param msgType The diagnostic message type
   * @param fileName The filement where the diagnostic occured
   * @param lineId The line where the diagnostic occured
   */
  void insertDiagnosticReport( DiagnosticKey const & diagnosticKey )
  {
    m_diagnosticHistory.get_inserted( {diagnosticKey.logPart, diagnosticKey.msgType,
                                       diagnosticKey.fileName, diagnosticKey.lineId} ).m_count++;
  }

  void serialize( stdVector< buffer_unit_type > & gStruct );

  void deserialize( stdVector< buffer_unit_type > const & globalAllocations,
                    stdVector< integer > const & recvCounts );

};

/**
 * @brief Template specialisation to convert a LogHistory to a table string.
 * @param logHistory The LogHistory object to convert.
 * @return The CSV string representation of the logHistory.
 */
template<>
string TableTextFormatter::toString< LogHistory >( LogHistory const & logHistory ) const;

}

#endif
