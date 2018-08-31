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
 * GEOSX is a free software; you can redistrubute it and/or modify it under
 * the terms of the GNU Lesser General Public Liscense (as published by the
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
  virtual void Execute( real64 const & time_n,
                        real64 const & dt,
                        int const cycleNumber,
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
  virtual void Cleanup( real64 const & time_n,
                        int const cycleNumber,
                        dataRepository::ManagedGroup * domain ) override;


  /// Documentation assignment
  virtual void FillDocumentationNode() override;

  /**
   * An event may have an arbitrary number of sub-events defined as children in the input xml.
   * e.g.: <Events>
   *         <PeriodicEvent name="base_event" ...>
   *           <PeriodicEvent name="sub_event" .../>
   *           ...
   *         </PeriodicEvent>
   *       </Events>
   */
  virtual void CreateChild( string const & childKey, string const & childName ) override;

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

  /// A pointer to the optional event target
  ExecutableGroup * m_target;

  struct viewKeyStruct
  {
    dataRepository::ViewKey eventTarget = { "target" };
    dataRepository::ViewKey beginTime = { "beginTime" };
    dataRepository::ViewKey endTime = { "endTime" };
    dataRepository::ViewKey forceDt = { "forceDt" };
    dataRepository::ViewKey lastTime = { "lastTime" };
    dataRepository::ViewKey lastCycle = { "lastCycle" };

    dataRepository::ViewKey allowSuperstep = { "allowSuperstep" };
    dataRepository::ViewKey allowSubstep = { "allowSubstep" };
    dataRepository::ViewKey substepFactor = { "substepFactor" };
  } viewKeys;

  ///Catalog interface
  using CatalogInterface = cxx_utilities::CatalogInterface< EventBase, std::string const &, ManagedGroup * const >;
  static CatalogInterface::CatalogType& GetCatalog();

  /// Access functions
  integer GetForecast(){ return m_eventForecast; }
  void SetForecast(integer forecast){ m_eventForecast = forecast; }

  integer GetExitFlag();
  void SetExitFlag(integer flag){ m_exitFlag = flag; }

private:
  integer m_eventForecast = 0;
  integer m_exitFlag = 0;
};

} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_MANAGERS_EVENTS_EVENTSBASE_HPP_ */
