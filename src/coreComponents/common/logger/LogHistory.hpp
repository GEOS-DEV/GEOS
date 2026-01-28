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


namespace geos
{

/**
 * @brief Statistics for a diagnostic message at a specific location
 */
struct MsgStatistics
{
  /// Key identifying the source location
  using LocationKey = std::pair< string, integer >;
  /// Source code location
  LocationKey locationKey;
  /// Number of occurrences on the current rank
  integer count;
  /// Total number of occurrences across all ranks
  integer totalCount;
};


/**
 * @brief Keep track of all diagnostic message occured during the simulation
 */
class LogHistory
{
public:


  /**
   * @brief Report a diagnostic message
   * @param logPartName The logPart the message occured
   * @param diagMsg The diagnostic message to record
   * @param threadCount
   */
  void notifyMsg( LogPart::Type logPartName, DiagnosticMsg const & diagMsg );

  /**
   * @return The const messageCounts
   */
  auto const & getMessageCounts() const
  { return m_messageCounts; }

  /**
   * @brief Sets the total count across all ranks for a specific message location
   * @param logPartType The log part type
   * @param msgType The message type (error, warning, etc.)
   * @param location The source code location (file, line)
   * @param totalCount The aggregated count across all MPI ranks
   */
  void setTotalCount( LogPart::Type logPartType, MsgType msgType,
                      MsgStatistics::LocationKey locationKey, integer totalCount )
  {
    m_messageCounts.at( logPartType ).at( msgType ).at( locationKey ).totalCount = totalCount;
  }

private:
  /**
   * @brief Hierarchical storage of message statistics
   * Structure: LogPart -> MsgType -> SourceLocation -> Statistics
   */
  stdMap< LogPart::Type, stdMap< MsgType, stdMap< MsgStatistics::LocationKey, MsgStatistics > > > m_messageCounts;

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
