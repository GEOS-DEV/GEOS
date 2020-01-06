/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */


#ifndef GEOSX_MANAGERS_EVENTMANAGER_HPP_
#define GEOSX_MANAGERS_EVENTMANAGER_HPP_

#include "dataRepository/Group.hpp"
#include "managers/Events/EventBase.hpp"


namespace geosx
{
namespace dataRepository
{
namespace keys
{
string const Events("Events");
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
  /// Main constructor
  EventManager( std::string const & name,
                Group * const parent );

  /// Destructor
  virtual ~EventManager() override;

  /// A method to add child events
  virtual Group * CreateChild( string const & childKey, string const & childName ) override;

  /// This function is used to expand any catalogs in the data structure
  virtual void ExpandObjectCatalogs() override;

  /**
   * The main execution loop for the code.  During each cycle, it will:
   *   - Calculate the event forecast (number of cycles until its expected execution)
   *   - Signal an event to prepare (forecast == 1)
   *   - Execute an event (forecast == 0)
   *   - Determine dt for the next cycle
   *   - Advance time, cycle, etc.
   */
  void Run(dataRepository::Group * domain);

  struct viewKeyStruct
  {
    static constexpr auto maxTimeString = "maxTime";
    static constexpr auto maxCycleString = "maxCycle";

    static constexpr auto timeString = "time";
    static constexpr auto dtString = "dt";
    static constexpr auto cycleString = "cycle";
    static constexpr auto currentSubEventString = "currentSubEvent";

    dataRepository::ViewKey time = { "time" };
    dataRepository::ViewKey dt = { "dt" };
    dataRepository::ViewKey cycle = { "cycle" };
    dataRepository::ViewKey maxTime = { "maxTime" };
    dataRepository::ViewKey maxCycle = { "maxCycle" };
    dataRepository::ViewKey currentSubEvent = { "currentSubEvent" };
  } viewKeys;

  /// Catalog interface
  using CatalogInterface = dataRepository::CatalogInterface< EventBase, std::string const &, Group * const >;
  static CatalogInterface::CatalogType& GetCatalog();

private:

  real64 m_maxTime;
  integer m_maxCycle;

  real64 m_time;
  real64 m_dt;
  integer m_cycle;
  integer m_currentSubEvent;
};


} /* namespace geosx */

#endif /* GEOSX_MANAGERS_EVENTMANAGER_HPP_ */
