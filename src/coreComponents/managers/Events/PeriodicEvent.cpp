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
  m_functionTarget(nullptr),
  m_timeFrequency(),
  m_cycleFrequency(),
  m_targetExactTimestep(),
  m_functionName(),
  m_functionInputObject(),
  m_functionInputSetname(),
  m_functionStatOption(),
  m_eventThreshold()
{
  RegisterViewWrapper(viewKeyStruct::timeFrequencyString, &m_timeFrequency, false )->
      setApplyDefaultValue(-1.0)->
      setInputFlag(InputFlags::OPTIONAL)->
      setDescription("event frequency (time)");

  RegisterViewWrapper(viewKeyStruct::cycleFrequencyString, &m_cycleFrequency, false )->
      setApplyDefaultValue(1)->
      setInputFlag(InputFlags::OPTIONAL)->
      setDescription("event frequency (cycle, Default)");

  RegisterViewWrapper(viewKeyStruct::targetExactTimestepString, &m_targetExactTimestep, false )->
      setApplyDefaultValue(-1)->
      setInputFlag(InputFlags::OPTIONAL)->
      setDescription("allows timesteps to be truncated to match time frequency perfectly");

  RegisterViewWrapper(viewKeyStruct::functionNameString, &m_functionName, false )->
      setInputFlag(InputFlags::OPTIONAL)->
      setDescription("Name of the symbolic math function");

  RegisterViewWrapper(viewKeyStruct::functionInputObjectString, &m_functionInputObject, false )->
      setInputFlag(InputFlags::OPTIONAL)->
      setDescription("Path of the function input object (directory format)");

  RegisterViewWrapper(viewKeyStruct::functionInputSetnameString, &m_functionInputSetname, false )->
      setInputFlag(InputFlags::OPTIONAL)->
      setDescription("Setname of the input object (if empty, default to everything)");

  RegisterViewWrapper(viewKeyStruct::functionStatOptionString, &m_functionStatOption, false )->
      setInputFlag(InputFlags::OPTIONAL)->
      setDescription("Selection of the min/avg/max for functions that target vectors");

  RegisterViewWrapper(viewKeyStruct::eventThresholdString, &m_eventThreshold, false )->
      setApplyDefaultValue(10000000000.0)->
      setInputFlag(InputFlags::OPTIONAL)->
      setDescription("event threshold");

}


PeriodicEvent::~PeriodicEvent()
{}






void PeriodicEvent::EstimateEventTiming(real64 const time,
                                        real64 const dt, 
                                        integer const cycle,
                                        ManagedGroup * domain)
{
  real64 const timeFrequency = this->getReference<real64>(periodicEventViewKeys.timeFrequency);
  integer const cycleFrequency = this->getReference<integer>(periodicEventViewKeys.cycleFrequency);
  real64 const lastTime = this->getReference<real64>(EventBase::viewKeys.lastTime);
  integer const lastCycle = this->getReference<integer>(EventBase::viewKeys.lastCycle);
  string const functionName = this->getReference<string>(periodicEventViewKeys.functionName);
  
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
  string const functionName = this->getReference<string>(periodicEventViewKeys.functionName);
  string const functionInputObject = this->getReference<string>(periodicEventViewKeys.functionInputObject);
  string const functionInputSetname = this->getReference<string>(periodicEventViewKeys.functionInputSetname);
  integer const functionStatOption = this->getReference<integer>(periodicEventViewKeys.functionStatOption);
  real64 const eventThreshold = this->getReference<real64>(periodicEventViewKeys.eventThreshold);

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
      dataRepository::ManagedGroup const * sets = m_functionTarget->GetGroup(periodicEventViewKeys.functionSetNames);
      mySet = sets->getReference< set<localIndex> >(functionInputSetname);
    }

    // Find the function (min, average, max)
    real64_array stats = function->EvaluateStats(m_functionTarget, time, mySet);
    result = stats[functionStatOption];

    // Because the function applied to an object may differ by rank, synchronize
    // (Note: this shouldn't occur very often, since it is only called if the base forecast <= 0)
    #if USE_MPI
      real64 result_global;
      MPI_Allreduce(&result, &result_global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      result = result_global;
    #endif
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
  real64 const timeFrequency = this->getReference<real64>(periodicEventViewKeys.timeFrequency);
  integer const targetExactTimestep = this->getReference<integer>(periodicEventViewKeys.targetExactTimestep);
  
  real64 requestedDt = EventBase::GetTimestepRequest(time);

  if ((timeFrequency > 0) && (targetExactTimestep > 0))
  {
    real64 nextTargetTime = lastTime + timeFrequency;
    real64 maxDt = nextTargetTime - time;

    // If the event is executing next cycle, then maxDt is the timeFrequency
    maxDt = (time >= nextTargetTime) ? timeFrequency : maxDt;
    requestedDt = (maxDt < requestedDt) ? maxDt : requestedDt;
  }

  return requestedDt;
}



REGISTER_CATALOG_ENTRY( EventBase, PeriodicEvent, std::string const &, ManagedGroup * const )
} /* namespace geosx */
