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

#ifndef GEOS_COMMON_MEMORYINFO_HPP_
#define GEOS_COMMON_MEMORYINFO_HPP_

#include "umpire/util/MemoryResourceTraits.hpp"
#include "common/DataTypes.hpp"
#include "common/logger/Logger.hpp"
#include "format/table/TableFormatter.hpp"
#include <unistd.h>
#include <iostream>
#if defined( GEOS_USE_CUDA )
#include <cuda.h>
#endif

namespace geos
{

/**
 * @class MemoryInfos
 * @brief Class to fetch and store memory information for different resource types.
 */
class MemoryInfos
{
public:

  /**
   * @brief Constructor for MemoryInfos.
   * @param resourceType The type of memory resource.
   */
  MemoryInfos( umpire::MemoryResourceTraits::resource_type resourceType );

  /**
   * @brief Get the total memory available for the resource type.
   * @return Total memory in bytes.
   */
  size_t getTotalMemory() const;

  /**
   * @brief Get the available memory for the resource type.
   * @return Available memory in bytes.
   */
  size_t getAvailableMemory() const;

  /**
   * @brief Check if physical memory is handled.
   * @return True if physical memory is handled, false otherwise.
   */
  bool isPhysicalMemoryHandled() const;
private:

  ///total memory available.
  size_t m_totalMemory;
  ///Available memory.
  size_t m_availableMemory;
  ///Flag indicating if physical memory is handled.
  bool m_physicalMemoryHandled;
};

/**
 * @brief Singleton class keeping the application memory logging settings & features.
 * @todo Manage the activation / deactivation of debug LvArray memory logging.
 */
class MemoryLogging
{
public:

  /**
   * @return the instance reference.
   */
  static MemoryLogging & getInstance();

  /**
   * @return the text log level.
   */
  bool isUmpireStatsLogOutputEnabled() const
  { return m_umpireStatsLogReport; }

  /**
   * @return true if the CSV output is enabled.
   * Refer to X "writeCSV" wrapper documentation for more documentation.
   */
  bool isUmpireStatsCsvOutputEnabled() const
  { return m_umpireStatsCsvReport; }

  /**
   * @return the current cycle for the next report to output.
   */
  integer getCurrentCycle() const
  { return m_currentCycle; }

  /**
   * @param enable enable or disable the umpire statistics text logging.
   * @see memoryStatsReport() for more documentation.
   */
  void enableUmpireStatsLogReport( bool enable )
  { m_umpireStatsLogReport = enable; }

  /**
   * @param enable enable or disable the text log level.
   * @param filename the relative file path to the csv output, including the extention. Can be empty if !enable.
   * @note when enabled, start a new csv immediately from zero, filling its header.
   * @see memoryStatsReport() for more documentation.
   */
  void enableUmpireStatsCsvReport( bool enable, string_view filename = "" );

  /**
   * @brief Set the Current Cycle, to identify each CSV entry.
   * @param currentCycle The current cycle id. Must be initialized before first memoryStatsReport() to valid values.
   *                     The entry after the end of the simulation must be 'last_cycle + 1'.
   */
  void setCurrentCycle( integer currentCycle )
  { m_currentCycle = currentCycle; }

  /**
   * @brief Set the Current Time, to identify each CSV entry.
   * @param currentTime The current simulated time. Must be initialized before first memoryStatsReport() to valid values.
   */
  void setCurrentTime( real64 currentTime )
  { m_currentTime = currentTime; }

  /**
   * @brief Output the umpire statistics according to settings set by enableUmpireStatsLogReport() and
   * enableUmpireStatsCsvReport().
   * The statistics are the Umpire total high water mark across all ranks for each umpire allocator.
   */
  void memoryStatsReport() const;

private:

  /// The table formatter of the stats for log output
  std::unique_ptr< TableTextFormatter > m_memoryStatLogFormatter;

  /// The table formatter of the stats for csv output
  std::unique_ptr< TableCSVFormatter > m_memoryStatCsvFormatter;

  /// Enable the umpire statistics text log report.
  bool m_umpireStatsLogReport;

  /// Enable the umpire statistics CSV report.
  bool m_umpireStatsCsvReport;

  /// the filename for the umpire statistics CSV report.
  string m_umpireStatsCsvReportFilename;

  /// The current cycle id.
  integer m_currentCycle;

  /// The current simulated time.
  real64 m_currentTime;


  /// private constructor as the class is a singleton.
  MemoryLogging();

};

}

#endif
