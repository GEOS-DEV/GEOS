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
  * @file SoloEvent.cpp
  */

#include "SoloEvent.hpp"

namespace geosx
{

using namespace dataRepository;


SoloEvent::SoloEvent( const std::string& name,
                              ManagedGroup * const parent ):
  EventBase(name,parent),
  m_targetTime(-1.0),
  m_targetCycle(-1),
  m_targetExactTimestep(0)
{
  RegisterViewWrapper(viewKeyStruct::targetTimeString, &m_targetTime, false )->
    setApplyDefaultValue(-1.0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Event time");

  RegisterViewWrapper(viewKeyStruct::targetCycleString, &m_targetCycle, false )->
    setApplyDefaultValue(-1)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Event cycle");

  RegisterViewWrapper(viewKeyStruct::targetExactTimestepString, &m_targetExactTimestep, false )->
    setApplyDefaultValue(1)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Allows timesteps to be truncated to match time frequency perfectly");
}


SoloEvent::~SoloEvent()
{}


void SoloEvent::EstimateEventTiming(real64 const time,
                                    real64 const dt, 
                                    integer const cycle,
                                    ManagedGroup * domain)
{
  // Check event status
  if (m_lastCycle < 0)
  {
    if (m_targetTime >= 0.0)
    {
      if (dt <= 0)
      {
        SetForecast(std::numeric_limits<integer>::max());
      }
      else
      {
        // Note: add a small value to this forecast to account for floating point errors
        real64 forecast = ((m_targetTime - time) / dt) + 1e-10;
        SetForecast(static_cast<integer>(std::min(forecast, 1e9)));
      }
    }
    else
    {
      SetForecast(m_targetCycle - cycle);
    }
  }
  else
  {
    SetForecast(std::numeric_limits<integer>::max());
  }
}


real64 SoloEvent::GetEventTypeDtRequest(real64 const time)
{
  real64 requestedDt = std::numeric_limits<real64>::max();

  // Note: if m_lastCycle is set, then the event has already executed
  if ((m_lastCycle < 0) && (m_targetTime > 0) && (m_targetExactTimestep > 0))
  {
    // This extra step is necessary to prevent the event manager from
    // falling into a dt=0 loop
    real64 tmp_t = std::nextafter(time, time + 1.0);
    if (tmp_t < m_targetTime)
    {
      requestedDt = std::min(requestedDt, m_targetTime - time);
    }
  }

  return requestedDt;
}



REGISTER_CATALOG_ENTRY( EventBase, SoloEvent, std::string const &, ManagedGroup * const )
} /* namespace geosx */
