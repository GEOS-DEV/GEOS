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

#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"

namespace geosx
{

using namespace dataRepository;


EventBase::EventBase( const std::string& name,
                      ManagedGroup * const parent ):
  ExecutableGroup(name, parent),
  m_target(nullptr)
{}


EventBase::~EventBase()
{}


EventBase::CatalogInterface::CatalogType& EventBase::GetCatalog()
{
  static EventBase::CatalogInterface::CatalogType catalog;
  return catalog;
}


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
                              "1.0e100",
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
                              "-1",
                              "",
                              0,
                              0,
                              0 );

  docNode->AllocateChildNode( viewKeys.lastCycle.Key(),
                              viewKeys.lastCycle.Key(),
                              -1,
                              "integer",
                              "integer",
                              "last event occurance (time)",
                              "last event occurance (time)",
                              "-1",
                              "",
                              0,
                              0,
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

  docNode->AllocateChildNode( viewKeys.targetExactStartStop.Key(),
                              viewKeys.targetExactStartStop.Key(),
                              -1,
                              "integer",
                              "integer",
                              "allows timesteps to be truncated to match the start/stop times exactly",
                              "allows timesteps to be truncated to match the start/stop times exactly",
                              "0",
                              "",
                              0,
                              1,
                              0 );

  docNode->AllocateChildNode( viewKeys.currentSubEvent.Key(),
                              viewKeys.currentSubEvent.Key(),
                              -1,
                              "integer",
                              "integer",
                              "index of the current subevent",
                              "index of the current subevent",
                              "0",
                              "",
                              0,
                              0,
                              0 );

  docNode->AllocateChildNode( viewKeys.isTargetExecuting.Key(),
                              viewKeys.isTargetExecuting.Key(),
                              -1,
                              "integer",
                              "integer",
                              "Flag to indicate whether the event target is executing",
                              "Flag to indicate whether the event target is executing.  This helps to avoid double-executions.",
                              "0",
                              "",
                              0,
                              0,
                              0 );

}


void EventBase::CreateChild( string const & childKey, string const & childName )
{
  GEOS_LOG_RANK_0("Adding Event: " << childKey << ", " << childName);
  std::unique_ptr<EventBase> event = EventBase::CatalogInterface::Factory( childKey, childName, this );
  this->RegisterGroup<EventBase>( childName, std::move(event) );
}


void EventBase::InitializePreSubGroups( ManagedGroup * const group )
{
  real64& lastTime = this->getReference<real64>(viewKeys.lastTime);
  integer& lastCycle = this->getReference<integer>(viewKeys.lastCycle);

  lastTime = std::numeric_limits<real64>::min();
  lastCycle = std::numeric_limits<integer>::min();
}


void EventBase::GetTargetReferences()
{
  string eventTarget = this->getReference<string>(viewKeys.eventTarget);
  if (!eventTarget.empty())
  {
    ManagedGroup * tmp = this->GetGroupByPath(eventTarget);
    m_target = ManagedGroup::group_cast<ExecutableGroup*>(tmp);
    GEOS_ERROR_IF(m_target == nullptr, "The target of an event must be executable! " << m_target);
  }

  this->forSubGroups<EventBase>([]( EventBase * subEvent ) -> void
  {
    subEvent->GetTargetReferences();
  });
}


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


void EventBase::SignalToPrepareForExecution(real64 const time,
                                        real64 const dt, 
                                        integer const cycle,
                                        dataRepository::ManagedGroup * domain)
{
  if (m_target != nullptr)
  {
    m_target->SignalToPrepareForExecution(time, dt, cycle, domain);
  }

  this->forSubGroups<EventBase>([&]( EventBase * subEvent ) -> void
  {
    if (subEvent->GetForecast() == 1)
    {
      subEvent->SignalToPrepareForExecution(time, dt, cycle, domain);
    }
  });
}


void EventBase::Execute(real64 const& time_n,
                        real64 const& dt,
                        const integer cycleNumber,
                        integer const ,
                        real64 const & ,
                        ManagedGroup * domain)
{
  GEOSX_MARK_FUNCTION;
  real64& lastTime = this->getReference<real64>(viewKeys.lastTime);
  integer& lastCycle = this->getReference<integer>(viewKeys.lastCycle);
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


void EventBase::Step(real64 const time,
                     real64 const dt,  
                     integer const cycle,
                     dataRepository::ManagedGroup * domain )
{
  GEOSX_MARK_FUNCTION;
  // currentSubEvent indicates which child event was active when the restart was written
  // isTargetExecuting blocks double-execution of the target during restarts, and is useful debug information in outputs
  integer& currentSubEvent = *(this->getData<integer>(viewKeys.currentSubEvent));
  integer& isTargetExecuting = *(this->getData<integer>(viewKeys.isTargetExecuting));

  if ((m_target != nullptr) && (isTargetExecuting == 0))
  {
    isTargetExecuting = 1;
    m_target->Execute(time, dt, cycle, m_eventCount, m_eventProgress, domain);
  }
  isTargetExecuting = 0;

  // Iterage using the managed integer currentSubEvent, which will
  // allow restart runs to pick up where they left off.
  for ( ; currentSubEvent<this->numSubGroups(); ++currentSubEvent)
  {
    EventBase * subEvent = static_cast<EventBase *>( this->GetSubGroups()[currentSubEvent] );

    if (subEvent->GetForecast() <= 0)
    {
      subEvent->Execute(time, dt, cycle, m_eventCount, m_eventProgress, domain);
    }
  }

  currentSubEvent = 0;
}


real64 EventBase::GetTimestepRequest(real64 const time)
{
  real64 nextDt = std::numeric_limits<integer>::max();
  real64 requestedDt = std::numeric_limits<integer>::max();
  real64 const forceDt = this->getReference<real64>(viewKeys.forceDt);
  integer const allowSubstep = this->getReference<integer>(viewKeys.allowSubstep);
  integer const substepFactor = this->getReference<integer>(viewKeys.substepFactor);
  integer const targetExactStartStop = this->getReference<integer>(viewKeys.targetExactStartStop);

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

  if (targetExactStartStop == 1)
  {
    real64 const beginTime = this->getReference<real64>(viewKeys.beginTime);
    real64 const endTime = this->getReference<real64>(viewKeys.endTime);

    if (time < beginTime)
    {
      nextDt = std::min(beginTime - time, nextDt);
    }
    else if (time < endTime)
    {
      nextDt = std::min(endTime - time, nextDt);
    }
  }

  return nextDt;
}


void EventBase::Cleanup(real64 const& time_n,
                        integer const cycleNumber,
                        integer const eventCounter,
                        real64 const & eventProgress,
                        ManagedGroup * domain)
{
  if (m_target != nullptr)
  {
    m_target->Cleanup(time_n, cycleNumber, m_eventCount, m_eventProgress, domain);
  }

  this->forSubGroups<EventBase>([&]( EventBase * subEvent ) -> void
  {
    subEvent->Cleanup(time_n, cycleNumber, m_eventCount, m_eventProgress, domain);
  });
}



integer EventBase::GetExitFlag()
{
  this->forSubGroups<EventBase>([&]( EventBase * subEvent ) -> void
  {
    m_exitFlag += subEvent->GetExitFlag();
  });

  return m_exitFlag;
}



void EventBase::GetExecutionOrder(array1d<integer> & eventCounters)
{
  // The first entry counts all events, the second tracks solver events
  m_eventCount = eventCounters[0];
  m_timeStepEventCount = eventCounters[1];

  // Increment counters
  ++eventCounters[0];
  if (m_target != nullptr)
  {
    if (m_target->GetTimestepBehavior() > 0)
    {
      ++eventCounters[1];
    }
  }

  this->forSubGroups<EventBase>([&]( EventBase * subEvent ) -> void
  {
    subEvent->GetExecutionOrder(eventCounters);
  });
}


void EventBase::SetProgressIndicator(array1d<integer> & eventCounters)
{
  // Calculate the event progress indicator
  // This is defined as the percent completion through the executaion loop
  // with respect to the beginning of the event.
  m_eventProgress = static_cast<real64>(m_timeStepEventCount) / static_cast<real64>(eventCounters[1]);
  
  // Do this for child events
  this->forSubGroups<EventBase>([&]( EventBase * subEvent ) -> void
  {
    subEvent->SetProgressIndicator(eventCounters);
  });
}


} /* namespace geosx */
