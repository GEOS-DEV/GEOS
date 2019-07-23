/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
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
 * @file EventBase.hpp
 */

#ifndef SRC_COMPONENTS_CORE_SRC_MANAGERS_EVENTS_EVENTSBASE_HPP_
#define SRC_COMPONENTS_CORE_SRC_MANAGERS_EVENTS_EVENTSBASE_HPP_

#include "dataRepository/ManagedGroup.hpp"
#include "dataRepository/ExecutableGroup.hpp"
#include "fileIO/schema/SchemaUtilities.hpp"


namespace geosx
{

/**
 * @class EventBase
 * A base class for managing code event targets (solver applications, etc.)
 */
class EventBase : public ExecutableGroup
{
public:
  /// Main constructor
  explicit EventBase( std::string const & name,
                       ManagedGroup * const parent );

  /// Destructor
  virtual ~EventBase() override;

  // Catalog name interface
  static string CatalogName() { return "EventBase"; }

  /**
   * If the event forecast is equal to 1, then signal the targets to prepare for execution
   * during the next cycle.
   */
  virtual void SignalToPrepareForExecution(real64 const time,
                                           real64 const dt,  
                                           integer const cycle,
                                           dataRepository::ManagedGroup * domain) override;
  /**
   * If the event forecast is equal to 0, then call the step function on its target and/or children.
   */
  virtual void Execute( real64 const time_n,
                        real64 const dt,
                        integer const cycleNumber,
                        integer const,
                        real64 const,
                        dataRepository::ManagedGroup * domain ) override;

  /**
   * This method will call the execute method on the target
   * and/or children if present.
   */
  void Step(real64 const time,
            real64 const dt,  
            integer const cycle,
            dataRepository::ManagedGroup * domain );

  /*
   * This method is called as the code exits the main run loop
   */
  virtual void Cleanup( real64 const time_n,
                        integer const cycleNumber,
                        integer const eventCounter,
                        real64 const eventProgress,
                        dataRepository::ManagedGroup * domain ) override;

  /**
   * An event may have an arbitrary number of sub-events defined as children in the input xml.
   * e.g.: <Events>
   *         <PeriodicEvent name="base_event" ...>
   *           <PeriodicEvent name="sub_event" .../>
   *           ...
   *         </PeriodicEvent>
   *       </Events>
   */
  virtual ManagedGroup * CreateChild( string const & childKey, string const & childName ) override;


  /// This function is used to expand any catalogs in the data structure
  virtual void ExpandObjectCatalogs() override;

  /**
   * The target object for an event may be specified via the keyword "target" in the input xml.
   * This string is empty by default and uses GetGroupByPath() method in ManagedGroup, which returns
   * a pointer to the target using a unix-style path as an input (both absolute and relative paths work).
   * This involves a lot of string parsing, so we do it once during initialization.
   */
  void GetTargetReferences();

  /**
   * Events are triggered based upon their forecast values, which are defined
   * as the expected number of code cycles before they are executed.  This method
   * will call EstimateEventTiming (defined in each subclass) on this event and
   * its children.
   */
  virtual void CheckEvents(real64 const time, 
                              real64 const dt,
                              integer const cycle,
                              dataRepository::ManagedGroup * domain);

  /// Method to estimate the timing of the event
  virtual void EstimateEventTiming(real64 const time, 
                                      real64 const dt,
                                      integer const cycle,
                                      dataRepository::ManagedGroup * domain) = 0;

  /**
   * This method will collect time-step size requests from its
   * targets and/or children.
   */
  virtual real64 GetTimestepRequest(real64 const time) override;


  /**
   * This method is used to get event-specifit dt requests
   */
  virtual real64 GetEventTypeDtRequest(real64 const time){ return std::numeric_limits<real64>::max(); }


  /// This method is used to count the number of events/sub-events
  void GetExecutionOrder(array1d<integer> & eventCounters);

  /**
   * This method is used to determine how to handle the timestamp for an event
   * If the event occurs after anything targeting a SolverBase object, then
   * set the m_isPostSolverEvent flag.  If set, then the time passed to the target
   * will be time + dt.
   */
  void SetProgressIndicator(array1d<integer> & eventCounters);



  struct viewKeyStruct
  {
    static constexpr auto eventTargetString = "target";
    static constexpr auto beginTimeString = "beginTime";
    static constexpr auto endTimeString = "endTime";
    static constexpr auto forceDtString = "forceDt";
    static constexpr auto maxEventDtString = "maxEventDt";
    static constexpr auto lastTimeString = "lastTime";
    static constexpr auto lastCycleString = "lastCycle";
    static constexpr auto targetExactStartStopString = "targetExactStartStop";
    static constexpr auto currentSubEventString = "currentSubEvent";
    static constexpr auto isTargetExecutingString = "isTargetExecuting";
    static constexpr auto verbosityString = "verbosity";

    dataRepository::ViewKey eventTarget = { "target" };
    dataRepository::ViewKey beginTime = { "beginTime" };
    dataRepository::ViewKey endTime = { "endTime" };
    dataRepository::ViewKey forceDt = { "forceDt" };
    dataRepository::ViewKey maxEventDt = { "maxEventDt" };
    dataRepository::ViewKey lastTime = { "lastTime" };
    dataRepository::ViewKey lastCycle = { "lastCycle" };
    dataRepository::ViewKey targetExactStartStop = { "targetExactStartStop" };
    dataRepository::ViewKey currentSubEvent = { "currentSubEvent" };
    dataRepository::ViewKey isTargetExecuting = { "isTargetExecuting" };
    dataRepository::ViewKey verbosity = { "verbosity" };
    } viewKeys;

  ///Catalog interface
  using CatalogInterface = cxx_utilities::CatalogInterface< EventBase, std::string const &, ManagedGroup * const >;
  static CatalogInterface::CatalogType& GetCatalog();

  /// Access functions
  integer GetForecast(){ return m_eventForecast; }
  void SetForecast(integer forecast){ m_eventForecast = forecast; }

  integer GetExitFlag();
  void SetExitFlag(integer flag){ m_exitFlag = flag; }

  integer GetEventCount() const { return m_eventCount; }
  real64  GetEventProgress() const { return m_eventProgress; }

  real64  GetCurrentEventDtRequest() const { return m_currentEventDtRequest; }


protected:
  real64 m_lastTime;
  integer m_lastCycle;


private:
  string m_eventTarget;
  real64 m_beginTime;
  real64 m_endTime;
  real64 m_forceDt;
  real64 m_maxEventDt;
  integer m_targetExactStartStop;
  integer m_currentSubEvent;
  integer m_targetExecFlag;
  integer m_eventForecast;
  integer m_exitFlag;
  integer m_eventCount;
  integer m_timeStepEventCount;
  real64 m_eventProgress;
  integer m_verbosity;
  real64 m_currentEventDtRequest;

  /// A pointer to the optional event target
  ExecutableGroup * m_target;

};

} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_MANAGERS_EVENTS_EVENTSBASE_HPP_ */
