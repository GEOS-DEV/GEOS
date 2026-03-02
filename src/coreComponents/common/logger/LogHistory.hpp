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

/**
 * @brief Statistics for a diagnostic message at a specific location
 */
struct MsgStatistics
{
  /// Number of times a message occured during the simulation
  integer count;
};


/**
 * @brief Keep track of all diagnostic message occured during the simulation
 */
class LogHistory
{
public:

  /// Alias for the historical diagnostic unordered_map key
  using HistoricalDiagnosticUnorderedMapKey = std::tuple< string, MsgType, string, integer >;

  /**
   * @brief Report a diagnostic message
   * @param logPartName The logPart where the message occured
   * @param diagMsg The diagnostic message to record
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

  /**
   * @brief Insert an element to the diagnostic history container if an equivalent key doesn't exist.
   * @param logPartName The logPart where the diagnostic occured
   * @param msgType The diagnostic message type
   * @param fileName The filement where the diagnostic occured
   * @param lineCount The line where the diagnostic occured
   */
  void insertDiagnosticReport( string_view logPartName, MsgType msgType,
                               string const & fileName, integer lineCount )
  {
    m_diagnosticHistory.get_inserted( std::make_tuple( string( logPartName ), msgType, fileName, lineCount ));
  }

private:

  /// @cond DO_NOT_DOCUMENT
  struct LocationKeyHash
  {

    std::size_t operator()( HistoricalDiagnosticUnorderedMapKey const & key ) const noexcept
    {
      auto const & [logPartType, msgType, filename, lineCount] = key;

      std::size_t h1 = std::hash< std::string >{} (logPartType);
      std::size_t h2 = std::hash< MsgType >{} (msgType);
      std::string str;
      str.assign( filename );
      std::size_t h3 = std::hash< std::string >{} (str);
      std::size_t h4 = std::hash< int >{} (lineCount);

      return h1 ^ (h2 << 1) ^ (h3 << 2) ^ (h4 << 3);
    }

  };
  /// @endcond

  /// @cond DO_NOT_DOCUMENT
  struct LocationKeyEqual
  {

    bool operator()( HistoricalDiagnosticUnorderedMapKey const & lhs,
                     HistoricalDiagnosticUnorderedMapKey const & rhs ) const
    {
      return std::get< 0 >( lhs ) == std::get< 0 >( rhs ) &&
             std::get< 1 >( lhs ) == std::get< 1 >( rhs ) &&
             std::get< 2 >( lhs ) == std::get< 2 >( rhs ) &&
             std::get< 2 >( lhs ) == std::get< 2 >( rhs );
    }
  };
  /// @endcond

  /**
   * @brief Diagnostic history happened during the simulation
   */
  stdUnorderedMap< HistoricalDiagnosticUnorderedMapKey,
                   MsgStatistics,
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
