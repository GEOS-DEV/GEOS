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

/*
 * EventBase.cpp
 *
 *  Created on: Jan 26, 2018
 *      Author: sherman
 */

#include "EventBase.hpp"

namespace geosx
{

using namespace dataRepository;

EventBase::EventBase( const std::string& name,
                      ManagedGroup * const parent ):
  ManagedGroup(name, parent),
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


void EventBase::CreateChild( string const & childKey, string const & childName )
{
  std::cout << "Adding Event: " << childKey << ", " << childName << std::endl;
  std::unique_ptr<EventBase> event = EventBase::CatalogInterface::Factory( childKey, childName, this );
  this->RegisterGroup<EventBase>( childName, std::move(event) );
}


void EventBase::GetTargetReferences()
{
  string eventTarget = this->getReference<string>(viewKeys.eventTarget);
  if (!eventTarget.empty())
  {
    m_target = this->GetGroupByPath(eventTarget);
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
    m_eventForecast = (beginTime - time) / dt;
  }
  else if (time >= endTime)
  {
    m_eventForecast = 1e6;
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


void EventBase::Execute(real64 const time,
                        real64 const dt, 
                        integer const cycle,
                        ManagedGroup * domain)
{
  real64& lastTime = *(this->getData<real64>(viewKeys.lastTime));
  integer& lastCycle = *(this->getData<integer>(viewKeys.lastCycle));
  integer const allowSuperstep = this->getReference<integer>(viewKeys.allowSuperstep);
  integer const allowSubstep = this->getReference<integer>(viewKeys.allowSubstep);
  integer const substepFactor = this->getReference<integer>(viewKeys.substepFactor);

  if (allowSuperstep > 0)
  {
    real64 actualDt = dt + time - lastTime;
    Step(lastTime, actualDt, cycle, domain);     // Should we use the lastTime or time here?
  }
  else if (allowSubstep > 0)
  {
    real64 actualDt = dt / substepFactor;
    real64 actualTime = time;

    for (integer ii=0; ii<substepFactor; ii++)
    {
      Step(actualTime, actualDt, cycle, domain);
      actualTime += actualDt;
    }
  }
  else
  {
    Step(time, dt, cycle, domain);
  }

  lastTime = time;
  lastCycle = cycle;
}


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


real64 EventBase::GetTimestepRequest(real64 const time)
{
  real64 nextDt = 1e6;
  real64 requestedDt = 1e6;
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
