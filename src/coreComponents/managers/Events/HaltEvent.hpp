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
 * @file HaltEvent.hpp
 */


#ifndef SRC_COMPONENTS_CORE_SRC_MANAGERS_EVENTS_HALTEVENT_HPP_
#define SRC_COMPONENTS_CORE_SRC_MANAGERS_EVENTS_HALTEVENT_HPP_

#include "managers/Events/EventBase.hpp"

namespace geosx
{

/**
 * @class HaltEvent
 * An event type that is designed to look at the external clock.
 * This is useful for managing wall time limitations.
 */
class HaltEvent : public EventBase
{
public:
  /// Constructor
  HaltEvent(const std::string& name,
                Group * const parent);
  
  /// Destructor
  virtual ~HaltEvent() override;

  // Catalog name interface
  static string CatalogName() { return "HaltEvent"; }
  
  /**
   * This event is designed to look at the external clock. Currently,
   * if the event is triggered it will set a flag, which will
   * instruct the code to exit.  This is useful for managing walltime
   */
  virtual void EstimateEventTiming(real64 const time,
                                   real64 const dt, 
                                   integer const cycle,
                                   dataRepository::Group * domain) override;

  /// Timing values
  real64 m_externalStartTime;
  real64 m_externalLastTime;
  real64 m_externalDt;
  real64 m_maxRuntime;

  struct viewKeyStruct
  {
    static constexpr auto maxRuntimeString = "maxRuntime";
    
    dataRepository::ViewKey maxRuntime = { "maxRuntime" };
  } haltEventViewKeys;

};

} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_MANAGERS_EVENTS_HALTEVENT_HPP_ */
