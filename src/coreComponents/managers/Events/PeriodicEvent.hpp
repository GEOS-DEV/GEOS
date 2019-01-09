/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * @file PeriodicEvent.hpp
 */

#ifndef SRC_COMPONENTS_CORE_SRC_MANAGERS_EVENTS_PERIODICEVENT_HPP_
#define SRC_COMPONENTS_CORE_SRC_MANAGERS_EVENTS_PERIODICEVENT_HPP_

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
  PeriodicEvent(const std::string& name,
                ManagedGroup * const parent);
  
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
  virtual void EstimateEventTiming(real64 const time,
                                   real64 const dt, 
                                   integer const cycle,
                                   dataRepository::ManagedGroup * domain) override;

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
  void CheckOptionalFunctionThreshold(real64 const time,
                                      real64 const dt, 
                                      integer const cycle,
                                      dataRepository::ManagedGroup * domain);

  /**
   * Grab the next time-step.  If requested, then limit the requested
   * dt to exactly match the time frequency
   */
  virtual real64 GetTimestepRequest(real64 const time) override;

  /// A pointer to an optional function
  dataRepository::ManagedGroup * m_functionTarget;

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

#endif /* SRC_COMPONENTS_CORE_SRC_MANAGERS_EVENTS_PERIODICEVENT_HPP_ */
