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
  /**
   * @brief Main constructor.
   * @param name The name of the object in the data repository.
   * @param parent The parent of this object in the data repository.
   **/
  explicit EventBase( std::string const & name,
                      Group * const parent );

  /// Destructor
  virtual ~EventBase() override;

  /**
   * @brief Catalog name interface.
   * @return This type's catalog name.
   **/
  static string CatalogName() { return "EventBase"; }

  /**
   * @brief If the event forecast is equal to 1, then signal the targets to prepare for execution
   *        during the next cycle.
   * @param time The current simulation time.
   * @param dt The current time increment.
   * @param cycle The current cycle.
   * @param domain The DomainPartition the event is occuring on up-casted to a Group.
   */
  virtual void SignalToPrepareForExecution( real64 const time,
                                            real64 const dt,
                                            integer const cycle,
                                            dataRepository::Group * domain ) override;
  /**
   * @brief If the event forecast is equal to 0, then call the step function on its target and/or children.
   * @copydoc ExecutableGroup::Execute()
   */
  virtual void Execute( real64 const time_n,
                        real64 const dt,
                        integer const cycleNumber,
                        integer const eventCounter,
                        real64 const eventProgress,
                        dataRepository::Group * domain ) override;

  /**
   * @brief Call the execute method on the target and/or children if present.
   * @param time The current simulation time.
   * @param dt The current simulation time increment.
   * @param cycle The current simulation cycle.
   * @param domain The DomainPartition up-casted to a Group.
   */
  void Step( real64 const time,
             real64 const dt,
             integer const cycle,
             dataRepository::Group * domain );

  /**
   * @copydoc dataRepository::Group::CreateChild()
   *
   * An event may have an arbitrary number of sub-events defined as children in the input xml.
   * e.g.:
   * @code{.unparsed}
   * <Events>
   *   <PeriodicEvent name="base_event" ...>
   *     <PeriodicEvent name="sub_event" .../>
   *     ...
   *   </PeriodicEvent>
   * </Events>
   * @endcode
   */
  virtual Group * CreateChild( string const & childKey, string const & childName ) override;

  /**
   * @brief Expand any catalogs in the data structure.
   */
  virtual void ExpandObjectCatalogs() override;

  /**
   * @brief Process input data to retrieve targeted objects internally.
   * The target object for an event may be specified via the keyword "target" in the input xml.
   * This string is empty by default and uses GetGroupByPath() method in Group, which returns
   * a pointer to the target using a unix-style path as an input (both absolute and relative paths work).
   * This involves a lot of string parsing, so we do it once during initialization.
   */
  void GetTargetReferences();

  /**
   * @brief Events are triggered based upon their forecast values, which are defined
   *        as the expected number of code cycles before they are executed.  This method
   *        will call EstimateEventTiming (defined in each subclass) on this event and
   *        its children.
   * @param time The current simulation time.
   * @param dt The current simulation time increment.
   * @param cycle the current simulation cycle.
   * @param domain The problem domain up-cast to a Group.
   */
  virtual void CheckEvents( real64 const time,
                            real64 const dt,
                            integer const cycle,
                            dataRepository::Group * domain );

  /**
   * @brief Perform the calculations to estimate the timing of the event.
   * @param time The current simulation time.
   * @param dt The current simulation time increment.
   * @param cycle the current simulation cycle.
   * @param domain The problem domain up-cast to a Group.
   */
  virtual void EstimateEventTiming( real64 const time,
                                    real64 const dt,
                                    integer const cycle,
                                    dataRepository::Group * domain ) = 0;

  /**
   * @brief Collect time-step size requests from targets and/or children.
   * @param time The current simulation time.
   * @return The requested time step.
   */
  virtual real64 GetTimestepRequest( real64 const time ) override;

  /**
   * @brief Get event-specifit dt requests.
   * @param time The current simulation time.
   * @return The requested time step.
   */
  virtual real64 GetEventTypeDtRequest( real64 const time )
  {
    GEOSX_UNUSED_VAR( time );
    return std::numeric_limits< real64 >::max();
  }


  /**
   * @brief Count the number of events/sub-events
   * @param[out] eventCounters The event count for each event/sub-event.
   */
  void GetExecutionOrder( array1d< integer > & eventCounters );

  /**
   * @brief Update the event progress for the event/sub-events.
   * @note This method is used to determine how to handle the timestamp for an event
   * @note If the event occurs after anything targeting a SolverBase object, then
   *       set the m_isPostSolverEvent flag.  If set, then the time passed to the target
   *       will be time + dt.
   * @param eventCounters The event count for each event/sub-event.
   */
  void SetProgressIndicator( array1d< integer > & eventCounters );

  /// @cond DO_NOT_DOCUMENT
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
  /// @endcond

  /// Catalog interface
  using CatalogInterface = dataRepository::CatalogInterface< EventBase, std::string const &, Group * const >;
  /// @copydoc dataRepository::Group::GetCatalog()
  static CatalogInterface::CatalogType & GetCatalog();

  /**
   * @brief Get the sum of the exit flags for the event/sub-events from the last execution.
   * @return The sum of the exit flags for the event/sub-events.
   */
  integer GetExitFlag();

  /**
   * @brief Set this event objects exit flag.
   * @param flag The exit flag value.
   */
  void SetExitFlag( integer flag ){ m_exitFlag = flag; }

  /**
   * @brief Get the current time increment request for this event.
   * @return The current time increment request.
   */
  real64  GetCurrentEventDtRequest() const { return m_currentEventDtRequest; }

  /**
   * @brief Forecast statuses.
   *
   * This enum is kept public only for debugging purpose.
   */
  enum class ForeCast : integer
  {
    READY_FOR_EXEC,
    PREPARE_FOR_EXEC,
    IDLE
  };

  /**
   * @brief Forecasts accessors
   * @return The forecast.
   *
   * The `GetForecast` getter only exists for debugging purpose.
   * Prefer the predicate functions below.
   */
  ///@{
  ForeCast getForecast() const
  { return m_eventForecast; }

  bool isReadyForExec() const
  { return m_eventForecast == ForeCast::READY_FOR_EXEC; }

  bool isPreparingForExec() const
  { return m_eventForecast == ForeCast::PREPARE_FOR_EXEC; }

  bool isIdle() const
  { return m_eventForecast == ForeCast::IDLE; }

protected:

  void setForecast( ForeCast forecast )
  { m_eventForecast = forecast; }

  /**
   * @brief Sets the forecast
   * @param forecast The forecast provided as an integer.
   *
   * If the forecast is 0 or below, the event is considered being READY_FOR_EXEC.
   * If it equals 1, it is in PREPARE_FOR_EXEC status. Above, the event is IDLE.
   */
  void setForecast( integer forecast );
  ///@}

  /**
   * @brief Is the event active?
   * @param time The time at which we want to check if the event is active.
   * @return True if acrive, false otherwise.
   */
  bool isActive( real64 const time ) const
  { return ( time >= m_beginTime ) && ( time < m_endTime ); }

  /**
   * @brief Get the target of this event.
   * @return The target of this event.
   */
  ExecutableGroup * GetEventTarget() const
  { return m_target; }

  /// The last time the event occurred.
  real64 m_lastTime;
  /// The last cycle the event occurred.
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
  ForeCast m_eventForecast;
  integer m_exitFlag;
  integer m_eventCount;
  integer m_timeStepEventCount;
  real64 m_eventProgress;
  real64 m_currentEventDtRequest;

  /// A pointer to the optional event target
  ExecutableGroup * m_target;
};

/**
 * @brief Prints a ForeCast to standard streams.
 * @param os The out stream
 * @param obj The forcast
 * @return The stream
 */
std::ostream & operator<<( std::ostream & os,
                           const EventBase::ForeCast & obj );

} /* namespace geosx */

#endif /* GEOSX_MANAGERS_EVENTS_EVENTSBASE_HPP_ */
