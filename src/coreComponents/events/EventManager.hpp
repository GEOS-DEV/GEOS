/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */


#ifndef GEOSX_EVENTS_EVENTMANAGER_HPP_
#define GEOSX_EVENTS_EVENTMANAGER_HPP_

#include "dataRepository/Group.hpp"
#include "EventBase.hpp"

namespace geosx
{

class DomainPartition;

namespace dataRepository
{
namespace keys
{
string const Events( "Events" );
}
}

/**
 * @class EventManager
 *
 * A class for managing code events.
 */
class EventManager : public dataRepository::Group
{
public:
  /**
   * @brief Constructor for the EventManager
   * @param[in] name the name of the EventManager
   * @param[in] parent group this EventManager
   */
  EventManager( string const & name,
                Group * const parent );

  /**
   * @brief Default destructor for the EventManager
   */
  virtual ~EventManager() override;

  /**
   * @brief Create a child EventManager
   * @param[in] childKey the key of the Event to be added
   * @param[in] childName the name of the child to be added
   * @return the Event
   */
  virtual Group * createChild( string const & childKey, string const & childName ) override;

  /**
   * @brief This method is used to expand any catalogs in the data structure
   */
  virtual void expandObjectCatalogs() override;

  /**
   * @brief The main execution loop for the code.
   * @param[in] domain the current DomainPartition on which the Event will be ran.
   * @return True iff the simulation exited early, and needs to be run again to completion.
   * @details During each cycle, it will:
   *   - Calculate the event forecast (number of cycles until its expected execution)
   *   - Signal an event to prepare (forecast == 1)
   *   - Execute an event (forecast == 0)
   *   - Determine dt for the next cycle
   *   - Advance time, cycle, etc.
   */
  bool run( DomainPartition & domain );

  /**
   * @name viewKeyStruct/groupKeyStruct
   */
  ///@{
  /// @cond DO_NOT_DOCUMENT
  struct viewKeyStruct
  {
    static constexpr char const * maxTimeString() { return "maxTime"; }
    static constexpr char const * maxCycleString() { return "maxCycle"; }

    static constexpr char const * timeString() { return "time"; }
    static constexpr char const * dtString() { return "dt"; }
    static constexpr char const * cycleString() { return "cycle"; }
    static constexpr char const * currentSubEventString() { return "currentSubEvent"; }

    dataRepository::ViewKey time = { "time" };
    dataRepository::ViewKey dt = { "dt" };
    dataRepository::ViewKey cycle = { "cycle" };
    dataRepository::ViewKey maxTime = { "maxTime" };
    dataRepository::ViewKey maxCycle = { "maxCycle" };
    dataRepository::ViewKey currentSubEvent = { "currentSubEvent" };
  } viewKeys;
  /// @endcond
  ///@}

  /// Alias to access the object catalog for EventBase derived types.
  using CatalogInterface = dataRepository::CatalogInterface< EventBase, string const &, Group * const >;

  /// @copydoc dataRepository::Group::getCatalog()
  static CatalogInterface::CatalogType & getCatalog();

private:
  /// Max time for a simulation
  real64 m_maxTime;

  /// Maximum number of cycles for a simulation
  integer m_maxCycle;

  /// Simulation timestamp at the beginning of the cycle
  real64 m_time;

  /// Current timestep request
  real64 m_dt;

  /// Current cycle
  integer m_cycle;

  /// Current subevent index
  integer m_currentSubEvent;
};


} /* namespace geosx */

#endif /* GEOSX_EVENTS_EVENTMANAGER_HPP_ */
