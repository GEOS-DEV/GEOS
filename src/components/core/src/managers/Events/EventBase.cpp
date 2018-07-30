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
 * @file EventBase.cpp
 */

#include "EventBase.hpp"
#include <cstring>

namespace geosx
{

using namespace dataRepository;

/*
 * Constructor.
 */
EventBase::EventBase( const std::string& name,
                      ManagedGroup * const parent ):
  ExecutableGroup(name, parent),
  m_target(nullptr)
{}


/**
 * Destructor.
 */
EventBase::~EventBase()
{}


/**
 * Catalog interface
 */
EventBase::CatalogInterface::CatalogType& EventBase::GetCatalog()
{
  static EventBase::CatalogInterface::CatalogType catalog;
  return catalog;
}


/**
 * Common documentation.
 */
void EventBase::FillDocumentationNode()
{
  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();

  docNode->setName("EventBase");
  docNode->setSchemaType("Node");
  docNode->setShortDescription("Describes the timing of the solver application");

  docNode->AllocateChildNode( viewKeys.eventTarget.Key(),
                              viewKeys.eventTarget.Key(),
                              -1,
                              "string",
                              "string",
                              "event target",
                              "event target",
                              "",
                              "",
                              0,
                              1,
                              0 );

  docNode->AllocateChildNode( viewKeys.beginTime.Key(),
                              viewKeys.beginTime.Key(),
                              -1,
                              "real64",
                              "real64",
                              "start time",
                              "start time",
                              "0.0",
                              "",
                              0,
                              1,
                              0 );

  docNode->AllocateChildNode( viewKeys.endTime.Key(),
                              viewKeys.endTime.Key(),
                              -1,
                              "real64",
                              "real64",
                              "end time",
                              "end time",
                              "1.0e9",
                              "",
                              0,
                              1,
                              0 );

  docNode->AllocateChildNode( viewKeys.forceDt.Key(),
                              viewKeys.forceDt.Key(),
                              -1,
                              "real64",
                              "real64",
                              "forced application dt",
                              "forced application dt",
                              "-1.0",
                              "",
                              0,
                              1,
                              0 );

  docNode->AllocateChildNode( viewKeys.lastTime.Key(),
                              viewKeys.lastTime.Key(),
                              -1,
                              "real64",
                              "real64",
                              "last event occurance (time)",
                              "last event occurance (time)",
                              "0",
                              "",
                              0,
                              1,
                              0 );

  docNode->AllocateChildNode( viewKeys.lastCycle.Key(),
                              viewKeys.lastCycle.Key(),
                              -1,
                              "integer",
                              "integer",
                              "last event occurance (time)",
                              "last event occurance (time)",
                              "0",
                              "",
                              0,
                              1,
                              0 );

  docNode->AllocateChildNode( viewKeys.allowSuperstep.Key(),
                              viewKeys.allowSuperstep.Key(),
                              -1,
                              "integer",
                              "integer",
                              "allows event super-stepping",
                              "allows event super-stepping (dt_super=dt+t-t_last)",
                              "0",
                              "",
                              0,
                              1,
                              0 );

  docNode->AllocateChildNode( viewKeys.allowSubstep.Key(),
                              viewKeys.allowSubstep.Key(),
                              -1,
                              "integer",
                              "integer",
                              "allows event sub-stepping",
                              "allows event sub-stepping",
                              "0",
                              "",
                              0,
                              1,
                              0 );

  docNode->AllocateChildNode( viewKeys.substepFactor.Key(),
                              viewKeys.substepFactor.Key(),
                              -1,
                              "integer",
                              "integer",
                              "integer substep factor",
                              "integer substep factor (dt_sub=dt/f)",
                              "1",
                              "",
                              0,
                              1,
                              0 );

}


/**
 * An event may have an arbitrary number of sub-events defined as children in the input xml.
 * e.g.: <Events>
 *         <PeriodicEvent name="base_event" ...>
 *           <PeriodicEvent name="sub_event" .../>
 *           ...
 *         </PeriodicEvent>
 *       </Events>
 */
void EventBase::CreateChild( string const & childKey, string const & childName )
{
  std::cout << "Adding Event: " << childKey << ", " << childName << std::endl;
  std::unique_ptr<EventBase> event = EventBase::CatalogInterface::Factory( childKey, childName, this );
  this->RegisterGroup<EventBase>( childName, std::move(event) );
}


/**
 * The target object for an event may be specified via the keyword "target" in the input xml.
 * This string is empty by default and uses GetGroupByPath() method in ManagedGroup, which returns
 * a pointer to the target using a unix-style path as an input (both absolute and relative paths work).
 * This involves a lot of string parsing, so we do it once during initialization.
 */
void EventBase::GetTargetReferences()
{
  string eventTarget = this->getReference<string>(viewKeys.eventTarget);
  if (!eventTarget.empty())
  {
    ManagedGroup * tmp = this->GetGroupByPath(eventTarget);
    std::cout << "Type of target = " << cxx_utilities::demangle(tmp->get_typeid().name()) << std::endl;
    if (dynamic_cast<ExecutableGroup *>(tmp) != nullptr)
    {
      m_target = ManagedGroup::group_cast<ExecutableGroup*>(tmp);
    }
    else
    {
      throw std::invalid_argument("The target of an event must be executable!");
    }    
  }

  this->forSubGroups<EventBase>([]( EventBase * subEvent ) -> void
  {
    subEvent->GetTargetReferences();
  });
}

/**
 * Events are triggered based upon their forecast values, which are defined
 * as the expected number of code cycles before they are executed.  This method
 * will call EstimateEventTiming (defined in each subclass) on this event and
 * its children.
 */
void EventBase::CheckEvents(real64 const time,
                            real64 const dt, 
                            integer const cycle,
                            ManagedGroup * domain)
{
  real64 const beginTime = this->getReference<real64>(viewKeys.beginTime);
  real64 const endTime = this->getReference<real64>(viewKeys.endTime);
  
  // Check event status
  if (time < beginTime)
  {
    if (dt <= 0)
    {
      m_eventForecast = std::numeric_limits<integer>::max();
    }
    else
    {
      m_eventForecast = int((beginTime - time) / dt);
    }
  }
  else if (time >= endTime)
  {
    m_eventForecast = std::numeric_limits<integer>::max();
  }
  else
  {
    this->EstimateEventTiming(time, dt, cycle, domain);
    
    // Check sub-events
    this->forSubGroups<EventBase>([&]( EventBase * subEvent ) -> void
    {
      subEvent->CheckEvents(time, dt, cycle, domain);
    });
  }
}


/**
 * If the event forecast is equal to 1, then signal the targets to prepare for execution
 * during the next cycle.
 */
void EventBase::SignalToPrepareForExecution(real64 const time,
                                        real64 const dt, 
                                        integer const cycle,
                                        dataRepository::ManagedGroup * domain)
{
  if (m_target != nullptr)
  {
    // Do something
  }

  this->forSubGroups<EventBase>([&]( EventBase * subEvent ) -> void
  {
    if (subEvent->GetForecast() == 1)
    {
      subEvent->SignalToPrepareForExecution(time, dt, cycle, domain);
    }
  });
}



/**
 * If the event forecast is equal to 0, then call the step function on its target and/or children.
 * There are three types of time-steps that are allowed:
 *   - Regular steps (default).  This will call execute the solver with the dt specified by this event's parent.
 *   - Superstep (allowSuperstep = 1).  The dt for the step will be set to (dt + time_n - lastTime)
 *   - Substep (allowSubstep = 1, substepFactor >= 1).  This will repeatedly step with timestep=dt/substepFactor
 */
void EventBase::Execute(real64 const& time_n,
                        real64 const& dt,
                        const int cycleNumber,
                        ManagedGroup * domain)
{
  real64& lastTime = *(this->getData<real64>(viewKeys.lastTime));
  integer& lastCycle = *(this->getData<integer>(viewKeys.lastCycle));
  integer const allowSuperstep = this->getReference<integer>(viewKeys.allowSuperstep);
  integer const allowSubstep = this->getReference<integer>(viewKeys.allowSubstep);
  integer const substepFactor = this->getReference<integer>(viewKeys.substepFactor);

  if (allowSuperstep > 0)
  {
    real64 actualDt = dt + time_n - lastTime;
    Step(lastTime, actualDt, cycleNumber, domain);     // Should we use the lastTime or time here?
  }
  else if (allowSubstep > 0)
  {
    real64 actualDt = dt / substepFactor;
    real64 actualTime = time_n;

    for (integer ii=0; ii<substepFactor; ii++)
    {
      Step(actualTime, actualDt, cycleNumber, domain);
      actualTime += actualDt;
    }
  }
  else
  {
    Step(time_n, dt, cycleNumber, domain);
  }

  lastTime = time_n;
  lastCycle = cycleNumber;
}


/**
 * This method will call the execute method on the target and/or children if present.
 */
void EventBase::Step(real64 const time,
                     real64 const dt,  
                     integer const cycle,
                     dataRepository::ManagedGroup * domain )
{
  if (m_target != nullptr)
  {
    m_target->Execute(time, dt, cycle, domain);
  }

  this->forSubGroups<EventBase>([&]( EventBase * subEvent ) -> void
  {
    if (subEvent->GetForecast() <= 0)
    {
      subEvent->Execute(time, dt, cycle, domain);
    }
  });
}



/**
 * This method will collect time-step size requests from its targets and/or children.
 */
real64 EventBase::GetTimestepRequest(real64 const time)
{
  real64 nextDt = std::numeric_limits<integer>::max();
  real64 requestedDt = std::numeric_limits<integer>::max();
  real64 const forceDt = this->getReference<real64>(viewKeys.forceDt);
  integer const allowSubstep = this->getReference<integer>(viewKeys.allowSubstep);
  integer const substepFactor = this->getReference<integer>(viewKeys.substepFactor);

  if (forceDt > 0)
  {
    nextDt = forceDt;
  }
  else
  {
    if (m_target != nullptr)
    {
      requestedDt = m_target->GetTimestepRequest(time);
      nextDt = std::min(requestedDt, nextDt);
    }

    this->forSubGroups<EventBase>([&]( EventBase * subEvent ) -> void
    {
      if (subEvent->GetForecast() <= 1)
      {
        requestedDt = subEvent->GetTimestepRequest(time);
        nextDt = std::min(requestedDt, nextDt);
      }
    });

    if (allowSubstep > 0)
    {
      nextDt *= substepFactor;
    }
  }

  return nextDt;
}


} /* namespace geosx */
