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

#include "common/DataTypes.hpp"
#include "common/StdContainerWrappers.hpp"
#include "common/format/LogPart.hpp"
#include "common/format/table/TableFormatter.hpp"
#include "DiagnosticMessage.hpp"
#include "common/logger/MsgType.hpp"
#include <string>
#include <unordered_set>


namespace geos
{

/**
 * @brief Keep track of all diagnostic message occured during the simulation
 */
class LogHistory
{
public:

  struct LogRecord
  {
    /// POD characterizing an unique diagnostic message
    struct Key
    {
      string filename;
      integer lineId;
    } m_key;

    struct Values
    {
      string logPart;
      MsgType msgType;
      integer m_count;
    } m_value;

    size_t getSerializedSize() const;
    LogRecord();
    LogRecord( Key const &, Values const & );
    LogRecord( stdVector< buffer_unit_type > & buffer );

    stdVector< buffer_unit_type > serialize() const;
    void deserialize( buffer_unit_type const * & );

    template< typename T >
    unsigned long sizeOfField( T ) const
    { return sizeof(T); }

    unsigned long sizeOfField( string_view str ) const
    { return sizeof(string::size_type) + str.size(); }

    template< typename T >
    typename std::enable_if_t< std::is_trivially_copyable_v< T > >
    writeInField( T & data, buffer_unit_type const * & ptr )
    {
      static_assert( std::is_trivially_copyable_v< T > );
      memcpy( &data, ptr, sizeof(T) );
      ptr += sizeof(T);
    }

    void writeInField( string & str, buffer_unit_type const * & ptr )
    {
      string::size_type strSize = 0;
      memcpy( &strSize, ptr, sizeof(string::size_type));
      ptr += sizeof(string::size_type);

      str.assign( ptr, ptr + strSize );
      ptr += str.size();
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

    size_t operator()( LogRecord::Key const & key ) const noexcept
    {
      size_t h1 = std::hash< string >{} (key.filename);
      size_t h2 = std::hash< integer >{} (key.lineId);

      return h1 ^ (h2 << 1);
    }

  };
  /// @endcond

  /// @cond DO_NOT_DOCUMENT
  struct LocationKeyEqual
  {

    bool operator()( LogRecord::Key const & lhs,
                     LogRecord::Key const & rhs ) const
    {
      return lhs.filename == rhs.filename  &&
             lhs.lineId == rhs.lineId;
    }
  };
  /// @endcond

  void insertDiagnosticReport( LogRecord log );

  /**
   * @brief Diagnostic history happened during the simulation
   */
  stdUnorderedMap< LogRecord::Key, LogRecord::Values,
                   LocationKeyHash, LocationKeyEqual > m_diagnosticHistory;
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
