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


namespace geos
{

/**
 * @brief Structure containing the statistics of a diagnostic message
 */
struct NumMsg
{
  /// the source location file
  string fileLocation;
  /// the source location line (default is 0)
  integer codeLocation;
  /// TODO
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
  void notifyMsg( LogPart::Type logPartName, DiagnosticMsg const & diagMsg, integer threadCount );

  /**
   * @return The const messageCounts
   */
  auto const & get() const
  { return m_messageCounts; }

private:
  /// Keep track of all classified diagnostic
  stdMap< LogPart::Type, stdMap< MsgType, stdVector< NumMsg > > > m_messageCounts;

};

/**
 * @brief Template specialisation to convert a LogHistory to a table string.
 * @param tableData The LogHistory object to convert.
 * @return The CSV string representation of the logHistory.
 */
template<>
string TableTextFormatter::toString< LogHistory >( LogHistory const & loghistory) const;

}

#endif
