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
 * @file HaltEvent.hpp
 * An event type that is designed to look at the external clock.
 * This is useful for managing wall time limitations.
 */


#ifndef SRC_COMPONENTS_CORE_SRC_MANAGERS_EVENTS_HALTEVENT_HPP_
#define SRC_COMPONENTS_CORE_SRC_MANAGERS_EVENTS_HALTEVENT_HPP_

#include "managers/Events/EventBase.hpp"

namespace geosx
{

class HaltEvent : public EventBase
{
public:
  HaltEvent(const std::string& name,
                ManagedGroup * const parent);
  
  virtual ~HaltEvent() override;

  static string CatalogName() { return "HaltEvent"; }

  virtual void FillDocumentationNode() override;
  
  virtual void EstimateEventTiming(real64 const time,
                                   real64 const dt, 
                                   integer const cycle,
                                   dataRepository::ManagedGroup * domain) override;

  real64 m_startTime = 0.0;
  real64 m_lastTime = 0.0;
  real64 m_realDt = 1.0;

  struct viewKeyStruct
  {
    dataRepository::ViewKey maxRuntime = { "maxRuntime" };
  } viewKeys;

};

} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_MANAGERS_EVENTS_HALTEVENT_HPP_ */
