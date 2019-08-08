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
  m_lastTime(-1.0e100),
  m_lastCycle(-1.0e9),
  m_eventTarget(""),
  m_beginTime(0.0),
  m_endTime(1e100),
  m_forceDt(-1.0),
  m_maxEventDt(-1.0),
  m_targetExactStartStop(0),
  m_currentSubEvent(0),
  m_targetExecFlag(0),
  m_eventForecast(0),
  m_exitFlag(0),
  m_eventCount(0),
  m_timeStepEventCount(0),
  m_eventProgress(0),
  m_verbosity(0),
  m_currentEventDtRequest(0.0),
  m_target(nullptr)
{
  setInputFlags(InputFlags::OPTIONAL_NONUNIQUE);
  
  RegisterViewWrapper(viewKeyStruct::eventTargetString, &m_eventTarget, false )->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Name of the object to be executed when the event criteria are met.");

  RegisterViewWrapper(viewKeyStruct::beginTimeString, &m_beginTime, false )->
    setApplyDefaultValue(0.0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Start time of this event.");

  RegisterViewWrapper(viewKeyStruct::endTimeString, &m_endTime, false )->
    setApplyDefaultValue(1e100)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("End time of this event.");

  RegisterViewWrapper(viewKeyStruct::forceDtString, &m_forceDt, false )->
    setApplyDefaultValue(-1.0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("While active, this event will request this timestep value (ignoring any children/targets requests).");

  RegisterViewWrapper(viewKeyStruct::maxEventDtString, &m_maxEventDt, false )->
    setApplyDefaultValue(-1.0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("While active, this event will request a timestep <= this value (depending upon any child/target requests).");

  RegisterViewWrapper(viewKeyStruct::targetExactStartStopString, &m_targetExactStartStop, false )->
    setApplyDefaultValue(1)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("If this option is set, the event will reduce its timestep requests to match any specified beginTime/endTimes exactly.");
 
 RegisterViewWrapper(viewKeyStruct::lastTimeString, &m_lastTime, false )->
    setApplyDefaultValue(-1.0e100)->
    setDescription("Last event occurrence (time)");

  RegisterViewWrapper(viewKeyStruct::lastCycleString, &m_lastCycle, false )->
    setApplyDefaultValue(-1.0e9)->
    setDescription("Last event occurrence (cycle)");

  RegisterViewWrapper(viewKeyStruct::currentSubEventString, &m_currentSubEvent, false )->
    setDescription("Index of the current subevent");

  RegisterViewWrapper(viewKeyStruct::isTargetExecutingString, &m_targetExecFlag, false )->
    setDescription("Index of the current subevent");

  RegisterViewWrapper(viewKeyStruct::verbosityString, &m_verbosity, false )->
    setApplyDefaultValue(0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Verbosity level");
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


void EventBase::ExpandObjectCatalogs()
{
  // Only add children if the parent is of type EventManager
  // otherwise, this would fall into a loop
  if (strcmp(this->getParent()->getName().c_str(), "Events") == 0)
  {
    for (auto& catalogIter: EventBase::GetCatalog())
    {
      CreateChild( catalogIter.first, catalogIter.first );
    }
  }
}


void EventBase::GetTargetReferences()
{
  if (!m_eventTarget.empty())
  {
    ManagedGroup * tmp = this->GetGroupByPath(m_eventTarget);
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
  // Check event status
  if (time < m_beginTime)
  {
    if (dt <= 0)
    {
      m_eventForecast = std::numeric_limits<integer>::max();
    }
    else
    {
      m_eventForecast = int((m_beginTime - time) / dt);
    }
  }
  else if (time >= m_endTime)
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


void EventBase::Execute(real64 const time_n,
                        real64 const dt,
                        const integer cycleNumber,
                        integer const,
                        real64 const,
                        ManagedGroup * domain)
{
  GEOSX_MARK_FUNCTION;
  
  GEOSX_MARK_BEGIN("EventBase::Execute() 1");
  // If m_targetExecFlag is set, then the code has resumed at a point
  // after the target has executed. 
  if ((m_target != nullptr) && (m_targetExecFlag == 0))
  {
    m_targetExecFlag = 1;
    m_target->Execute(time_n, dt, cycleNumber, m_eventCount, m_eventProgress, domain);
  }
  GEOSX_MARK_END("EventBase::Execute() 1");
  

  GEOSX_MARK_BEGIN("EventBase::Execute() 2");
  // Iterate through the sub-event list using the managed integer m_currentSubEvent
  // This allows for  restart runs to pick up where they left off.
  for ( ; m_currentSubEvent < this->numSubGroups(); ++m_currentSubEvent)
  {
    EventBase * subEvent = static_cast<EventBase *>( this->GetSubGroups()[m_currentSubEvent] );
    integer subEventForecast = subEvent->GetForecast();

    if (m_verbosity > 0)
    {
      GEOS_LOG_RANK_0("          SubEvent: " << m_currentSubEvent << " (" << subEvent->getName() << "), dt_request=" << subEvent->GetCurrentEventDtRequest() << ", forecast=" << subEventForecast);
    }

    if (subEventForecast <= 0)
    {
      subEvent->Execute(time_n, dt, cycleNumber, m_eventCount, m_eventProgress, domain);
    }
  }
  GEOSX_MARK_END("EventBase::Execute() 2");

  // Update the event status
  m_targetExecFlag = 0;
  m_currentSubEvent = 0;
  m_lastTime = time_n;
  m_lastCycle = cycleNumber;
}


real64 EventBase::GetTimestepRequest(real64 const time)
{
  m_currentEventDtRequest = std::numeric_limits<real64>::max();

  // Events and their targets may request a max dt when active
  if ((time >= m_beginTime) && (time < m_endTime))
  {
    if (m_forceDt > 0)
    {
      // Override the event dt request
      m_currentEventDtRequest = m_forceDt;
    }
    else
    {
      if (m_maxEventDt > 0)
      {
        // Limit the event dt request
        m_currentEventDtRequest = std::min(m_maxEventDt, m_currentEventDtRequest);
      }

      // Get the event-specific dt request
      m_currentEventDtRequest = std::min(m_currentEventDtRequest, GetEventTypeDtRequest(time));

      // Get the target's dt request if the event has the potential to execute this cycle
      if ((m_eventForecast <= 1) && (m_target != nullptr))
      {
        m_currentEventDtRequest = std::min(m_currentEventDtRequest, m_target->GetTimestepRequest(time));
      }

      // Get the sub-event dt requests
      this->forSubGroups<EventBase>([&]( EventBase * subEvent ) -> void
      {
        m_currentEventDtRequest = std::min(m_currentEventDtRequest, subEvent->GetTimestepRequest(time));
      });
    }
  }

  // Try to respect the start/stop times of the event window
  if (m_targetExactStartStop == 1)
  {
    // Using this instead of the raw time will prevent falling into dt = 0 loops:
    real64 tmp_t = std::nextafter(time, time + 1.0);

    if (tmp_t < m_beginTime)
    {
      m_currentEventDtRequest = std::min(m_beginTime - time, m_currentEventDtRequest);
    }
    else if (tmp_t < m_endTime)
    {
      m_currentEventDtRequest = std::min(m_endTime - time, m_currentEventDtRequest);
    }
  }

  return m_currentEventDtRequest;
}


void EventBase::Cleanup(real64 const time_n,
                        integer const cycleNumber,
                        integer const eventCounter,
                        real64 const eventProgress,
                        ManagedGroup * domain)
{
  if (m_target != nullptr)
  {
    // Cleanup the target
    m_target->Cleanup(time_n, cycleNumber, m_eventCount, m_eventProgress, domain);
  }

  // Cleanup any sub-events
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
  if (eventCounters[1] > 0)
  {
    m_eventProgress = static_cast<real64>(m_timeStepEventCount) / static_cast<real64>(eventCounters[1]);
  }
  
  // Do this for child events
  this->forSubGroups<EventBase>([&]( EventBase * subEvent ) -> void
  {
    subEvent->SetProgressIndicator(eventCounters);
  });
}


} /* namespace geosx */
