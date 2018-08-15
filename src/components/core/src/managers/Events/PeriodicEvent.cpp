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
  * @file PeriodicEvent.cpp
  */

#include "PeriodicEvent.hpp"
#include "DocumentationNode.hpp"
#include "managers/Functions/NewFunctionManager.hpp"

namespace geosx
{

using namespace dataRepository;


PeriodicEvent::PeriodicEvent( const std::string& name,
                              ManagedGroup * const parent ):
  EventBase(name,parent),
  m_functionTarget(nullptr)
{}


PeriodicEvent::~PeriodicEvent()
{}


void PeriodicEvent::FillDocumentationNode()
{
  EventBase::FillDocumentationNode();
  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();

  docNode->setName("PeriodicEvent");
  docNode->setSchemaType("Node");
  docNode->setShortDescription("Describes the timing of the solver application");

  docNode->AllocateChildNode( viewKeys.timeFrequency.Key(),
                              viewKeys.timeFrequency.Key(),
                              -1,
                              "real64",
                              "real64",
                              "event frequency (time)",
                              "event frequency (time)",
                              "-1.0",
                              "",
                              0,
                              1,
                              0 );

  docNode->AllocateChildNode( viewKeys.cycleFrequency.Key(),
                              viewKeys.cycleFrequency.Key(),
                              -1,
                              "integer",
                              "integer",
                              "event frequency (cycle, Default)",
                              "event frequency (cycle, Default)",
                              "1",
                              "",
                              0,
                              1,
                              0 );

  docNode->AllocateChildNode( viewKeys.targetExactTimestep.Key(),
                              viewKeys.targetExactTimestep.Key(),
                              -1,
                              "integer",
                              "integer",
                              "allows timesteps to be truncated to match time frequency perfectly",
                              "allows timesteps to be truncated to match time frequency perfectly",
                              "-1",
                              "",
                              0,
                              1,
                              0 );

  // Function oriented options
  docNode->AllocateChildNode( viewKeys.functionName.Key(),
                              viewKeys.functionName.Key(),
                              -1,
                              "string",
                              "string",
                              "Name of the symbolic math function",
                              "Name of the symbolic math function",
                              "",
                              "",
                              1,
                              1,
                              0 );

  docNode->AllocateChildNode( viewKeys.functionInputObject.Key(),
                              viewKeys.functionInputObject.Key(),
                              -1,
                              "string",
                              "string",
                              "Path of the function input object",
                              "Path of the function input object (directory format)",
                              "",
                              "",
                              1,
                              1,
                              0 );

  docNode->AllocateChildNode( viewKeys.functionInputSetname.Key(),
                              viewKeys.functionInputSetname.Key(),
                              -1,
                              "string",
                              "string",
                              "Setname of the input object",
                              "Setname of the input object (if empty, default to everything)",
                              "",
                              "",
                              1,
                              1,
                              0 );

  docNode->AllocateChildNode( viewKeys.functionStatOption.Key(),
                              viewKeys.functionStatOption.Key(),
                              -1,
                              "integer",
                              "integer",
                              "Selection of the min/avg/max for functions that target vectors",
                              "Selection of the min/avg/max for functions that target vectors",
                              "",
                              "",
                              1,
                              1,
                              0 );

  docNode->AllocateChildNode( viewKeys.eventThreshold.Key(),
                              viewKeys.eventThreshold.Key(),
                              -1,
                              "real64",
                              "real64",
                              "event threshold",
                              "event threshold",
                              "1e10",
                              "",
                              0,
                              1,
                              0 );
}



void PeriodicEvent::EstimateEventTiming(real64 const time,
                                        real64 const dt, 
                                        integer const cycle,
                                        ManagedGroup * domain)
{
  real64 const timeFrequency = this->getReference<real64>(viewKeys.timeFrequency);
  integer const cycleFrequency = this->getReference<integer>(viewKeys.cycleFrequency);
  real64 const lastTime = this->getReference<real64>(EventBase::viewKeys.lastTime);
  integer const lastCycle = this->getReference<integer>(EventBase::viewKeys.lastCycle);
  string const functionName = this->getReference<string>(viewKeys.functionName);
  
  // Check event status
  if (cycle == 0)
  {
    SetForecast(0);
  } 
  else if (timeFrequency >= 0.0)
  {
    if (dt <= 0)
    {
      SetForecast(std::numeric_limits<integer>::max());
    } 
    else
    {
      real64 forecast = (timeFrequency - (time - lastTime)) / dt;
      SetForecast(static_cast<integer>(std::min(forecast, 1e9)));
    }
  }
  else
  {
    SetForecast(cycleFrequency - (cycle - lastCycle));
  }

  if ((this->GetForecast() <= 0) && (functionName.empty() == 0))
  {
    CheckOptionalFunctionThreshold(time, dt, cycle, domain);
  }
}


void PeriodicEvent::CheckOptionalFunctionThreshold(real64 const time,
                                                   real64 const dt, 
                                                   integer const cycle,
                                                   ManagedGroup * domain)
{
  string const functionName = this->getReference<string>(viewKeys.functionName);
  string const functionInputObject = this->getReference<string>(viewKeys.functionInputObject);
  string const functionInputSetname = this->getReference<string>(viewKeys.functionInputSetname);
  integer const functionStatOption = this->getReference<integer>(viewKeys.functionStatOption);
  real64 const eventThreshold = this->getReference<real64>(viewKeys.eventThreshold);

  // Grab the function
  NewFunctionManager * functionManager = NewFunctionManager::Instance();
  FunctionBase * function = functionManager->GetGroup<FunctionBase>(functionName);

  real64 result = 0.0;
  if (functionInputObject.empty())
  {
    // This is a time-only function
    result = function->Evaluate(&time);
  }
  else
  {
    // Link the target object
    if (m_functionTarget == nullptr)
    {
      m_functionTarget = this->GetGroupByPath(functionInputObject);
    }

    // Get the set
    set<localIndex> mySet;
    if (functionInputSetname.empty())
    {
      for(localIndex ii=0; ii<m_functionTarget->size(); ++ii)
      {
        mySet.insert(ii);
      }
    }
    else
    {
      dataRepository::ManagedGroup const * sets = m_functionTarget->GetGroup(viewKeys.functionSetNames);
      mySet = sets->getReference< set<localIndex> >(functionInputSetname);
    }

    // Find the function (min, average, max)
    real64_array stats = function->EvaluateStats(m_functionTarget, time, mySet);
    result = stats[functionStatOption];
  }
  
  // Forcast event
  if (result > eventThreshold)
  {
    SetForecast(0);
  }
  else
  {
    SetForecast(std::numeric_limits<integer>::max());
  }
}


real64 PeriodicEvent::GetTimestepRequest(real64 const time)
{
  real64 const lastTime = this->getReference<real64>(EventBase::viewKeys.lastTime);
  real64 const timeFrequency = this->getReference<real64>(viewKeys.timeFrequency);
  integer const targetExactTimestep = this->getReference<integer>(viewKeys.targetExactTimestep);
  
  real64 requestedDt = EventBase::GetTimestepRequest(time);

  if ((timeFrequency > 0) && (targetExactTimestep > 0))
  {
    real64 nextTargetTimestep = lastTime + timeFrequency - time;
    requestedDt = (nextTargetTimestep < requestedDt) ? nextTargetTimestep : requestedDt;
  }

  return requestedDt;
}



REGISTER_CATALOG_ENTRY( EventBase, PeriodicEvent, std::string const &, ManagedGroup * const )
} /* namespace geosx */
