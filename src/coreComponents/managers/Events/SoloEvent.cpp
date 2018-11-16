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
  * @file SoloEvent.cpp
  */

#include "SoloEvent.hpp"
#include "DocumentationNode.hpp"

namespace geosx
{

using namespace dataRepository;


SoloEvent::SoloEvent( const std::string& name,
                              ManagedGroup * const parent ):
  EventBase(name,parent)
{}


SoloEvent::~SoloEvent()
{}


void SoloEvent::FillDocumentationNode()
{
  EventBase::FillDocumentationNode();
  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();

  docNode->setName("SoloEvent");
  docNode->setSchemaType("Node");
  docNode->setShortDescription("Describes the timing of the solver application");

  docNode->AllocateChildNode( SoloEventViewKeys.targetTime.Key(),
                              SoloEventViewKeys.targetTime.Key(),
                              -1,
                              "real64",
                              "real64",
                              "event time",
                              "event time",
                              "-1",
                              "",
                              0,
                              1,
                              0 );

  docNode->AllocateChildNode( SoloEventViewKeys.targetCycle.Key(),
                              SoloEventViewKeys.targetCycle.Key(),
                              -1,
                              "integer",
                              "integer",
                              "event cycle",
                              "event cycle",
                              "-1",
                              "",
                              0,
                              1,
                              0 );

  docNode->AllocateChildNode( SoloEventViewKeys.targetExactTimestep.Key(),
                              SoloEventViewKeys.targetExactTimestep.Key(),
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
}



void SoloEvent::EstimateEventTiming(real64 const time,
                                    real64 const dt, 
                                    integer const cycle,
                                    ManagedGroup * domain)
{
  real64 const targetTime = this->getReference<real64>(SoloEventViewKeys.targetTime);
  integer const targetCycle = this->getReference<integer>(SoloEventViewKeys.targetCycle);
  integer const lastCycle = this->getReference<integer>(EventBase::viewKeys.lastCycle);
  
  // Check event status
  if (lastCycle < 0)
  {
    if (targetTime >= 0.0)
    {
      if (dt <= 0)
      {
        SetForecast(std::numeric_limits<integer>::max());
      }
      else
      {
        real64 forecast = (targetTime - time) / dt;
        SetForecast(static_cast<integer>(std::min(forecast, 1e9)));
      }
    }
    else
    {
      SetForecast(targetCycle - cycle);
    }
  }
  else
  {
    SetForecast(std::numeric_limits<integer>::max());
  }
}


real64 SoloEvent::GetTimestepRequest(real64 const time)
{
  integer const targetExactTimestep = this->getReference<integer>(SoloEventViewKeys.targetExactTimestep);
  real64 const targetTime = this->getReference<real64>(SoloEventViewKeys.targetTime);
  
  real64 requestedDt = EventBase::GetTimestepRequest(time);

  if ((targetTime > 0) && (targetExactTimestep > 0))
  {
    requestedDt = std::min(requestedDt, targetTime - time);
  }

  return requestedDt;
}



REGISTER_CATALOG_ENTRY( EventBase, SoloEvent, std::string const &, ManagedGroup * const )
} /* namespace geosx */
