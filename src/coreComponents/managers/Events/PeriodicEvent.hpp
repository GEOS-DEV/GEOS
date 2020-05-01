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
 * @file PeriodicEvent.hpp
 */

#ifndef GEOSX_MANAGERS_EVENTS_PERIODICEVENT_HPP_
#define GEOSX_MANAGERS_EVENTS_PERIODICEVENT_HPP_

#include "managers/Events/EventBase.hpp"

namespace geosx
{

/**
 * @class PeriodicEvent
 *
 * An event type for periodic events (using either time or cycle as a basis).
 */
class PeriodicEvent : public EventBase
{
public:
  /// Main constructor
  PeriodicEvent( const std::string & name,
                 Group * const parent );

  /// Destructor
  virtual ~PeriodicEvent() override;

  /// Catalog name interface
  static string CatalogName() { return "PeriodicEvent"; }

  /**
   * Estimate the expected number of cycles until an event is expected to trigger.
   * The event frequency can be specified in terms of:
   *   - time (timeFrequency > 0, units = seconds)
   *   - or cycle (cycleFrequency >= 0, units = cycles)
   *
   * In addition, there is an optional function input that will be called if the
   * the nominal forecast (based on timing) is zero.
   */
  virtual void EstimateEventTiming( real64 const time,
                                    real64 const dt,
                                    integer const cycle,
                                    dataRepository::Group * domain ) override;

  /**
   * If the event forecast is zero, and an optional function (f) is specified, then
   * this method will be called to see if the event should be triggered or ignored.
   * For example, this could be used to periodically check the condition of the mesh,
   * and trigger a cleanup if necessary.
   *
   * If functionInputObject is not specified:
   *   - The argument to the function will be the current time
   *   - The event will be executed if f(t) >= eventThreshold
   *
   * If functionInputObject is specified:
   *   - The function will be called on the object, with the arguments given by functionInputSetname
   *   - The function manager will return a set of statistics
   *   - functionStatOption selects the statistic to compare against the eventThreshold (0 = min, 1 = average, 2 = max)
   *   - The event will be executed if f(object, arguments)[stat] >= eventThreshold
   */
  void CheckOptionalFunctionThreshold( real64 const time,
                                       real64 const dt,
                                       integer const cycle,
                                       dataRepository::Group * domain );

  /**
   * Grab the next time-step.  If requested, then limit the requested
   * dt to exactly match the time frequency
   */
  virtual real64 GetEventTypeDtRequest( real64 const time ) override;


  /*
   * This method is called as the code exits the main run loop
   */
  virtual void Cleanup( real64 const time_n,
                        integer const cycleNumber,
                        integer const eventCounter,
                        real64 const eventProgress,
                        dataRepository::Group * domain ) override;


  /// A pointer to an optional function
  dataRepository::Group * m_functionTarget;

  struct viewKeyStruct
  {
    static constexpr auto timeFrequencyString = "timeFrequency";
    static constexpr auto cycleFrequencyString = "cycleFrequency";
    static constexpr auto targetExactTimestepString = "targetExactTimestep";
    static constexpr auto functionNameString = "function";
    static constexpr auto functionInputObjectString = "object";
    static constexpr auto functionInputSetnameString = "set";
    static constexpr auto functionSetNamesString = "setNames";
    static constexpr auto functionStatOptionString = "stat";
    static constexpr auto eventThresholdString = "threshold";


    dataRepository::ViewKey timeFrequency = { "timeFrequency" };
    dataRepository::ViewKey cycleFrequency = { "cycleFrequency" };
    dataRepository::ViewKey targetExactTimestep = { "targetExactTimestep" };
    dataRepository::ViewKey functionName = { "function" };
    dataRepository::ViewKey functionInputObject = { "object" };
    dataRepository::ViewKey functionInputSetname = { "set" };
    dataRepository::ViewKey functionSetNames = { "setNames" };
    dataRepository::ViewKey functionStatOption = { "stat" };
    dataRepository::ViewKey eventThreshold = { "threshold" };
  } periodicEventViewKeys;

  real64 m_timeFrequency;
  integer m_cycleFrequency;
  integer m_targetExactTimestep;
  string m_functionName;
  string m_functionInputObject;
  string m_functionInputSetname;
  integer m_functionStatOption;
  real64 m_eventThreshold;



};

} /* namespace geosx */

#endif /* GEOSX_MANAGERS_EVENTS_PERIODICEVENT_HPP_ */
