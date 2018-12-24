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
 * @file SoloEvent.hpp
 */

#ifndef SRC_COMPONENTS_CORE_SRC_MANAGERS_EVENTS_SOLOEVENT_HPP_
#define SRC_COMPONENTS_CORE_SRC_MANAGERS_EVENTS_SOLOEVENT_HPP_

#include "managers/Events/EventBase.hpp"

namespace geosx
{

/**
 * @class SoloEvent
 *
 * An event type for events that occur only once.
 */
class SoloEvent : public EventBase
{
public:
  /// Main constructor
  SoloEvent(const std::string& name,
                ManagedGroup * const parent);
  
  /// Destructor
  virtual ~SoloEvent() override;

  /// Catalog name interface
  static string CatalogName() { return "SoloEvent"; }

  /**
   * Estimate the expected number of cycles until an event is expected to trigger.
   */
  virtual void EstimateEventTiming(real64 const time,
                                   real64 const dt, 
                                   integer const cycle,
                                   dataRepository::ManagedGroup * domain) override;

  /**
   * Grab the next time-step.  If requested, then limit the requested
   * dt to exactly match the time frequency
   */
  virtual real64 GetTimestepRequest(real64 const time) override;


  struct viewKeyStruct
  {
    dataRepository::ViewKey targetTime = { "targetTime" };
    dataRepository::ViewKey targetCycle = { "targetCycle" };
    dataRepository::ViewKey targetExactTimestep = { "targetExactTimestep" };
  } SoloEventViewKeys;

};

} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_MANAGERS_EVENTS_SOLOEVENT_HPP_ */
