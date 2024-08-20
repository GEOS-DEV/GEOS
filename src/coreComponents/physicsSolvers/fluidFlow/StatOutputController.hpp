/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#ifndef GEOS_PHYSICSSOLVERS_STATOUTPUTCONTROLLER_HPP_
#define GEOS_PHYSICSSOLVERS_STATOUTPUTCONTROLLER_HPP_
#include "physicsSolvers/SolverBase.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBase.hpp"

#include "dataRepository/BufferOpsDevice.hpp"
#include "dataRepository/HistoryDataSpec.hpp"
#include "events/tasks/TaskBase.hpp"
#include "mesh/DomainPartition.hpp"
#include "fileIO/Outputs/TimeHistoryOutput.hpp"
#include "fileIO/Outputs/OutputManager.hpp"
#include "fileIO/timeHistory/PackCollection.hpp"
#include "events/PeriodicEvent.hpp"

#include <functional>

namespace geos
{

class StatOutputController : public TaskBase
{

public:
/**
 * @brief Construct a new Stat Output Controller object
 *
 * @param name
 * @param parent
 */
  StatOutputController( string const & name,
                        Group * const parent );

  /// Accessor for the catalog name
  static string catalogName() { return "StatOutputController"; }

  /**
   * @copydoc dataRepository::Group::createChild()
   *
   * An event may have an arbitrary number of sub-events defined as children in the input xml.
   * e.g.:
   * @code{.unparsed}
   * <Events>
   *   <PeriodicEStatOutputControllervent name="testStats" ...>
   *     <CompositionalMultiphaseStatistics name="sub_event" .../>
   *   </PeriodicEvent>
   * </Events>
   * @endcode
   */
  virtual Group * createChild( string const & childKey, string const & childName ) override;

  /**
   * @brief Expand any catalogs in the data structure.
   */
  virtual void expandObjectCatalogs() override;

  /**
   * @defgroup Tasks Interface Functions
   *
   * This function implements the interface defined by the abstract TaskBase class
   */
  /**@{*/

  virtual bool execute( real64 const time_n,
                        real64 const dt,
                        integer const cycleNumber,
                        integer const eventCounter,
                        real64 const eventProgress,
                        DomainPartition & domain ) override;

private:

  void postInputInitialization() override;

  void initializePreSubGroups() override;

  TimeHistoryOutput * m_timeHistory;
  std::vector< PeriodicEvent * >  m_periodicEvents;
  std::vector< PackCollection * >  m_packCollections;

};

} /* namespace geos */

#endif
