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

  EventBase() = default;
  EventBase( EventBase const & ) = default;
  EventBase( EventBase &&) = default;
  EventBase& operator=( EventBase const & ) = default;
  EventBase& operator=( EventBase&& ) = default;

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
   * There are three types of time-steps that are allowed:
   *   - Regular steps (default).  This will call execute the solver with the dt specified by this event's parent.
   *   - Superstep (allowSuperstep = 1).  The dt for the step will be set to (dt + time_n - lastTime)
   *   - Substep (allowSubstep = 1, substepFactor >= 1).  This will repeatedly step with timestep=dt/substepFactor
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


  /**
   * This function is used to inform the schema generator of any
   * deviations between the xml and GEOS data structures.
   */
  virtual void SetSchemaDeviations(xmlWrapper::xmlNode schemaRoot,
                                   xmlWrapper::xmlNode schemaParent) override;


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
    static constexpr auto lastTimeString = "lastTime";
    static constexpr auto lastCycleString = "lastCycle";

    static constexpr auto allowSuperstepString = "allowSuperstep";
    static constexpr auto allowSubstepString = "allowSubstep";
    static constexpr auto substepFactorString = "substepFactor";
    static constexpr auto targetExactStartStopString = "targetExactStartStop";

    static constexpr auto currentSubEventString = "currentSubEvent";
    static constexpr auto isTargetExecutingString = "isTargetExecuting";


    dataRepository::ViewKey eventTarget = { "target" };
    dataRepository::ViewKey beginTime = { "beginTime" };
    dataRepository::ViewKey endTime = { "endTime" };
    dataRepository::ViewKey forceDt = { "forceDt" };
    dataRepository::ViewKey lastTime = { "lastTime" };
    dataRepository::ViewKey lastCycle = { "lastCycle" };

    dataRepository::ViewKey allowSuperstep = { "allowSuperstep" };
    dataRepository::ViewKey allowSubstep = { "allowSubstep" };
    dataRepository::ViewKey substepFactor = { "substepFactor" };
    dataRepository::ViewKey targetExactStartStop = { "targetExactStartStop" };

    dataRepository::ViewKey currentSubEvent = { "currentSubEvent" };
    dataRepository::ViewKey isTargetExecuting = { "isTargetExecuting" };
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

private:
  string m_eventTarget;
  real64 m_beginTime;
  real64 m_endTime;
  real64 m_forceDt;
  integer m_allowSuperstep;
  integer m_allowSubstep;
  integer m_substepFactor;
  integer m_targetExactStartStop;

  integer m_currentSubEvent;
  integer m_isTargetExecuting;
  integer m_eventForecast = 0;
  integer m_exitFlag = 0;
  integer m_eventCount = 0;
  integer m_timeStepEventCount = 0;
  real64 m_eventProgress = 0;
  real64 m_lastTime;
  integer m_lastCycle;

  /// A pointer to the optional event target
  ExecutableGroup * m_target;

};

} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_MANAGERS_EVENTS_EVENTSBASE_HPP_ */
