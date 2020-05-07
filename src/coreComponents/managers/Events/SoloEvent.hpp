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

/**
 * @file SoloEvent.hpp
 */

#ifndef GEOSX_MANAGERS_EVENTS_SOLOEVENT_HPP_
#define GEOSX_MANAGERS_EVENTS_SOLOEVENT_HPP_

#include "managers/Events/EventBase.hpp"

namespace geosx
{

/**
 * @class SoloEvent
 *
 * An event type for events that occur only once.
 */
class SoloEvent : public EventBase
{
public:
  /// Main constructor
  SoloEvent( const std::string & name,
             Group * const parent );

  /// Destructor
  virtual ~SoloEvent() override;

  /// Catalog name interface
  static string CatalogName() { return "SoloEvent"; }

  /**
   * Estimate the expected number of cycles until an event is expected to trigger.
   */
  virtual void EstimateEventTiming( real64 const time,
                                    real64 const dt,
                                    integer const cycle,
                                    dataRepository::Group * domain ) override;

  /**
   * Grab the next time-step.  If requested, then limit the requested
   * dt to exactly match the application time
   */
  virtual real64 GetEventTypeDtRequest( real64 const time ) override;


  struct viewKeyStruct
  {
    static constexpr auto targetTimeString = "targetTime";
    static constexpr auto targetCycleString = "targetCycle";
    static constexpr auto targetExactTimestepString = "targetExactTimestep";

    dataRepository::ViewKey targetTime = { "targetTime" };
    dataRepository::ViewKey targetCycle = { "targetCycle" };
    dataRepository::ViewKey targetExactTimestep = { "targetExactTimestep" };
  } SoloEventViewKeys;

  real64 m_targetTime;
  integer m_targetCycle;
  integer m_targetExactTimestep;

};

} /* namespace geosx */

#endif /* GEOSX_MANAGERS_EVENTS_SOLOEVENT_HPP_ */
