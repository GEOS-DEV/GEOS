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
  ///
  integer count;
};


/**
 * @brief Keep track of all diagnostic message occured during the simulation
 */
class LogHistory
{
public:

  /// Alias for the historical error unordered_map key
  using HistoricalErrorUnorderedMapKey = std::tuple< LogPart::Type, MsgType, string, integer >;

  /**
   * @brief Report a diagnostic message
   * @param logPartName The logPart the message occured
   * @param diagMsg The diagnostic message to record
   * @param threadCount
   */
  void notifyMsg( LogPart::Type logPartName, DiagnosticMsg const & diagMsg );

  /**
   * @brief Display the error statistics to the log
   */
  void diagnosticStatsReport();

  /**
   * @return The const historical error
   */
  auto const & getDiagnosticHistory() const
  { return m_errorHistory; }

  /**
   * @brief Insert an element to the error history container if an equivalent key doesn't exist.
   * @param logPartName The logPart where the error occured
   * @param msgType The error message type
   * @param locationKey The key identifying the error source location
   */
  void insertBlanckReport( LogPart::Type logPartName, MsgType msgType, string const & fileName, integer lineCount );

private:

  /// @cond DO_NOT_DOCUMENT
  struct LocationKeyHash
  {

    std::size_t operator()( HistoricalErrorUnorderedMapKey const & key ) const noexcept
    {
      auto const & [logPartType, msgType, filename, lineCount] = key;

      std::size_t h1 = std::hash< LogPart::Type >{} (logPartType);

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

    bool operator()( HistoricalErrorUnorderedMapKey const & lhs,
                     HistoricalErrorUnorderedMapKey const & rhs ) const
    {
      return std::get< 0 >( lhs ) == std::get< 0 >( rhs ) &&
             std::get< 1 >( lhs ) == std::get< 1 >( rhs ) &&
             std::get< 2 >( lhs ) == std::get< 2 >( rhs ) &&
             std::get< 2 >( lhs ) == std::get< 2 >( rhs );
    }
  };
  /// @endcond

  /**
   * @brief Historical error happened during the simulation
   */
  stdUnorderedMap< HistoricalErrorUnorderedMapKey,
                   MsgStatistics,
                   LocationKeyHash, LocationKeyEqual > m_errorHistory;
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
