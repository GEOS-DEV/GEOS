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
  m_eventTarget(""),
  m_beginTime(0.0),
  m_endTime(1e100),
  m_forceDt(-1.0),
  m_allowSuperstep(0),
  m_allowSubstep(0),
  m_substepFactor(1),
  m_targetExactStartStop(0),
  m_currentSubEvent(0),
  m_isTargetExecuting(0),
  m_eventForecast(0),
  m_exitFlag(0),
  m_eventCount(0),
  m_timeStepEventCount(0),
  m_eventProgress(0),
  m_lastTime(1e100),
  m_lastCycle(0),
  m_target(nullptr)
{
  RegisterViewWrapper(viewKeyStruct::eventTargetString, &m_eventTarget, false )->
      setInputFlag(InputFlags::REQUIRED)->
      setDescription("event target");

  RegisterViewWrapper(viewKeyStruct::beginTimeString, &m_beginTime, false )->
      setDefaultValue(0.0)->
      setInputFlag(InputFlags::OPTIONAL)->
      setDescription("Start time of this event");

  RegisterViewWrapper(viewKeyStruct::endTimeString, &m_endTime, false )->
      setDefaultValue(1e100)->
      setInputFlag(InputFlags::OPTIONAL)->
      setDescription("End time of this event");

  RegisterViewWrapper(viewKeyStruct::forceDtString, &m_forceDt, false )->
      setDefaultValue(-1.0)->setToDefaultValue()->
      setInputFlag(InputFlags::OPTIONAL)->
      setDescription("Forced timestep for this event");

  RegisterViewWrapper(viewKeyStruct::allowSuperstepString, &m_allowSuperstep, false )->
      setDefaultValue(0)->
      setInputFlag(InputFlags::OPTIONAL)->
      setDescription("allows event super-stepping (dt_super=dt+t-t_last)");

  RegisterViewWrapper(viewKeyStruct::allowSubstepString, &m_allowSubstep, false )->
      setDefaultValue(0)->
      setInputFlag(InputFlags::OPTIONAL)->
      setDescription("allows event sub-stepping");

  RegisterViewWrapper(viewKeyStruct::substepFactorString, &m_substepFactor, false )->
      setDefaultValue(1)->
      setInputFlag(InputFlags::OPTIONAL)->
      setDescription("integer substep factor (dt_sub=dt/f)");

  RegisterViewWrapper(viewKeyStruct::targetExactStartStopString, &m_targetExactStartStop, false )->
      setDefaultValue(0)->
      setInputFlag(InputFlags::OPTIONAL)->
      setDescription("allows timesteps to be truncated to match the start/stop times exactly");


  RegisterViewWrapper(viewKeyStruct::lastTimeString, &m_lastTime, false )->
      setDefaultValue(-1.0e100)->setToDefaultValue()->
      setDescription("last event occurrence (time)");

  RegisterViewWrapper(viewKeyStruct::lastCycleString, &m_lastCycle, false )->
      setDefaultValue(-1.0e9)->setToDefaultValue()->
      setDescription("last event occurrence (cycle)");

  RegisterViewWrapper(viewKeyStruct::currentSubEventString, &m_currentSubEvent, false )->
      setDescription("index of the current subevent");

  RegisterViewWrapper(viewKeyStruct::isTargetExecutingString, &m_isTargetExecuting, false )->
      setDescription("index of the current subevent");


}


EventBase::~EventBase()
{}


EventBase::CatalogInterface::CatalogType& EventBase::GetCatalog()
{
  static EventBase::CatalogInterface::CatalogType catalog;
  return catalog;
}

ManagedGroup * EventBase::CreateChild( string const & childKey, string const & childName )
{
  GEOS_LOG_RANK_0("Adding Event: " << childKey << ", " << childName);
  std::unique_ptr<EventBase> event = EventBase::CatalogInterface::Factory( childKey, childName, this );
  return this->RegisterGroup<EventBase>( childName, std::move(event) );
}


/*
void EventBase::InitializePreSubGroups( ManagedGroup * const group )
{
  real64& lastTime = this->getReference<real64>(viewKeys.lastTime);
  integer& lastCycle = this->getReference<integer>(viewKeys.lastCycle);

  lastTime = std::numeric_limits<real64>::min();
  lastCycle = std::numeric_limits<integer>::min();
}
*/


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
  integer& currentSubEvent = this->getReference<integer>(viewKeys.currentSubEvent);
  integer& isTargetExecuting = this->getReference<integer>(viewKeys.isTargetExecuting);

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
  real64 const beginTime = this->getReference<real64>(viewKeys.beginTime);
  real64 const endTime = this->getReference<real64>(viewKeys.endTime);
  real64 const forceDt = this->getReference<real64>(viewKeys.forceDt);
  integer const allowSubstep = this->getReference<integer>(viewKeys.allowSubstep);
  integer const substepFactor = this->getReference<integer>(viewKeys.substepFactor);
  integer const targetExactStartStop = this->getReference<integer>(viewKeys.targetExactStartStop);

  if (time >= endTime)
  {
    // This is the final timestep for this event, don't include it in the calculation
  }
  else if (forceDt > 0)
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
