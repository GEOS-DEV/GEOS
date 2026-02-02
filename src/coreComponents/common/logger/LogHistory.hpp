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
  /// Key identifying the source location
  using LocationKey = std::pair< std::array< char, 200 >, integer >;
  /// Source code location
  LocationKey locationKey;
  /// Number of occurrences on the current rank
  integer count;
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

  void insertBlanckReport(LogPart::Type logPartName, MsgType msgType, MsgStatistics::LocationKey locationKey );

private:
  struct LocationKeyHash
  {
    std::size_t operator()( std::tuple< LogPart::Type, MsgType, MsgStatistics::LocationKey > const & key ) const noexcept
    {
      auto const & [logPartType, msgType, locationKey] = key;

      std::size_t h1 = std::hash< LogPart::Type >{} (logPartType);

      std::size_t h2 = std::hash< MsgType >{} (msgType);

      std::string str( std::begin( locationKey.first ), std::end( locationKey.first ));
      std::size_t h3 = std::hash< std::string >{} (str);
      std::size_t h4 = std::hash< int >{} (locationKey.second);

      return h1 ^ (h2 << 1) ^ (h3 << 2) ^ (h4 << 3);
    }
  };
  struct LocationKeyEqual
  {
    bool operator()( const std::tuple< LogPart::Type, MsgType, MsgStatistics::LocationKey > & lhs,
                     const std::tuple< LogPart::Type, MsgType, MsgStatistics::LocationKey > & rhs ) const
    {
      return std::get< 0 >( lhs ) == std::get< 0 >( rhs ) &&
             std::get< 1 >( lhs ) == std::get< 1 >( rhs ) &&
             std::get< 2 >( lhs ).first == std::get< 2 >( rhs ).first &&
             std::get< 2 >( lhs ).second == std::get< 2 >( rhs ).second;
    }
  };

  /**
   * @brief Hierarchical storage of message statistics
   * Structure: LogPart -> MsgType -> SourceLocation -> Statistics
   */
  stdUnorderedMap< std::tuple< LogPart::Type, MsgType, MsgStatistics::LocationKey >,
                   MsgStatistics,
                   LocationKeyHash, LocationKeyEqual >m_messageCounts;
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
