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
 * @file EventBase.hpp
 */

#ifndef GEOSX_MANAGERS_EVENTS_EVENTSBASE_HPP_
#define GEOSX_MANAGERS_EVENTS_EVENTSBASE_HPP_

#include "dataRepository/Group.hpp"
#include "dataRepository/ExecutableGroup.hpp"
#include "fileIO/schema/schemaUtilities.hpp"


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
                      Group * const parent );

  /// Destructor
  virtual ~EventBase() override;

  // Catalog name interface
  static string CatalogName() { return "EventBase"; }

  /**
   * If the event forecast is equal to 1, then signal the targets to prepare for execution
   * during the next cycle.
   */
  virtual void SignalToPrepareForExecution( real64 const time,
                                            real64 const dt,
                                            integer const cycle,
                                            dataRepository::Group * domain ) override;
  /**
   * If the event forecast is equal to 0, then call the step function on its target and/or children.
   */
  virtual void Execute( real64 const time_n,
                        real64 const dt,
                        integer const cycleNumber,
                        integer const,
                        real64 const,
                        dataRepository::Group * domain ) override;

  /**
   * This method will call the execute method on the target
   * and/or children if present.
   */
  void Step( real64 const time,
             real64 const dt,
             integer const cycle,
             dataRepository::Group * domain );

  /**
   * An event may have an arbitrary number of sub-events defined as children in the input xml.
   * e.g.: <Events>
   *         <PeriodicEvent name="base_event" ...>
   *           <PeriodicEvent name="sub_event" .../>
   *           ...
   *         </PeriodicEvent>
   *       </Events>
   */
  virtual Group * CreateChild( string const & childKey, string const & childName ) override;


  /// This function is used to expand any catalogs in the data structure
  virtual void ExpandObjectCatalogs() override;

  /**
   * The target object for an event may be specified via the keyword "target" in the input xml.
   * This string is empty by default and uses GetGroupByPath() method in Group, which returns
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
  virtual void CheckEvents( real64 const time,
                            real64 const dt,
                            integer const cycle,
                            dataRepository::Group * domain );

  /// Method to estimate the timing of the event
  virtual void EstimateEventTiming( real64 const time,
                                    real64 const dt,
                                    integer const cycle,
                                    dataRepository::Group * domain ) = 0;

  /**
   * This method will collect time-step size requests from its
   * targets and/or children.
   */
  virtual real64 GetTimestepRequest( real64 const time ) override;


  /**
   * This method is used to get event-specifit dt requests
   */
  virtual real64 GetEventTypeDtRequest( real64 const GEOSX_UNUSED_PARAM( time ) ) { return std::numeric_limits< real64 >::max(); }


  /// This method is used to count the number of events/sub-events
  void GetExecutionOrder( array1d< integer > & eventCounters );

  /**
   * This method is used to determine how to handle the timestamp for an event
   * If the event occurs after anything targeting a SolverBase object, then
   * set the m_isPostSolverEvent flag.  If set, then the time passed to the target
   * will be time + dt.
   */
  void SetProgressIndicator( array1d< integer > & eventCounters );



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
    static constexpr auto finalDtStretchString = "finalDtStretch";

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
  } viewKeys;

  ///Catalog interface
  using CatalogInterface = dataRepository::CatalogInterface< EventBase, std::string const &, Group * const >;
  static CatalogInterface::CatalogType & GetCatalog();

  /// Access functions
  integer GetForecast(){ return m_eventForecast; }
  void SetForecast( integer forecast ){ m_eventForecast = forecast; }

  integer GetExitFlag();
  void SetExitFlag( integer flag ){ m_exitFlag = flag; }

  integer GetEventCount() const { return m_eventCount; }
  real64  GetEventProgress() const { return m_eventProgress; }

  real64  GetCurrentEventDtRequest() const { return m_currentEventDtRequest; }

  real64  GetBeginTime() const { return m_beginTime; }
  real64  GetEndTime() const { return m_endTime; }
  ExecutableGroup * GetEventTarget() const { return m_target; }


protected:
  real64 m_lastTime;
  integer m_lastCycle;


private:
  string m_eventTarget;
  real64 m_beginTime;
  real64 m_endTime;
  real64 m_forceDt;
  real64 m_maxEventDt;
  real64 m_finalDtStretch;
  integer m_targetExactStartStop;
  integer m_currentSubEvent;
  integer m_targetExecFlag;
  integer m_eventForecast;
  integer m_exitFlag;
  integer m_eventCount;
  integer m_timeStepEventCount;
  real64 m_eventProgress;
  real64 m_currentEventDtRequest;


  /// A pointer to the optional event target
  ExecutableGroup * m_target;

};

} /* namespace geosx */

#endif /* GEOSX_MANAGERS_EVENTS_EVENTSBASE_HPP_ */
