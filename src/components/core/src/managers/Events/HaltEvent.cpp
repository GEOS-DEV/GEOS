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
 * HaltEvent.cpp
 *
 *  Created on: Jan 26, 2018
 *      Author: sherman
 */

#include "HaltEvent.hpp"
#include "DocumentationNode.hpp"
#include <sys/time.h>

namespace geosx
{

using namespace dataRepository;

HaltEvent::HaltEvent( const std::string& name,
                        ManagedGroup * const parent ):
  EventBase(name,parent)
{
  timeval tim;
  gettimeofday(&tim, NULL);
  m_startTime = tim.tv_sec + (tim.tv_usec / 1000000.0);
  m_lastTime = m_startTime;  
}

HaltEvent::~HaltEvent()
{}

void HaltEvent::FillDocumentationNode()
{
  EventBase::FillDocumentationNode();
  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();

  docNode->setName("HaltEvent");
  docNode->setSchemaType("Node");
  docNode->setShortDescription("Describes the timing of the solver application");

  docNode->AllocateChildNode( viewKeys.maxRuntime.Key(),
                              viewKeys.maxRuntime.Key(),
                              -1,
                              "real64",
                              "real64",
                              "max runtime",
                              "max runtime",
                              "cycle",
                              "",
                              0,
                              1,
                              0 );

}



void HaltEvent::EstimateEventTiming(real64 const time,
                                     real64 const dt, 
                                     integer const cycle,
                                     ManagedGroup * domain)
{
  real64 const maxRuntime = this->getReference<real64>(viewKeys.maxRuntime);
  
  // Check run time
  timeval tim;
  gettimeofday(&tim, NULL);
  real64 currentTime = tim.tv_sec + (tim.tv_usec / 1000000.0);

  // Update values
  m_realDt = currentTime - m_lastTime;
  m_lastTime = currentTime;
  SetForecast((maxRuntime - (currentTime - m_startTime)) / m_realDt);

  if (this->GetForecast() <= 0)
  {
    this->SetExitFlag(1);
  }
}


REGISTER_CATALOG_ENTRY( EventBase, HaltEvent, std::string const &, ManagedGroup * const )
} /* namespace geosx */
