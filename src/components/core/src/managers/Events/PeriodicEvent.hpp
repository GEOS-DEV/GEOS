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
 * PeriodicEvent.hpp
 *
 *  Created on: Jan 26, 2018
 *      Author: sherman
 */

#ifndef SRC_COMPONENTS_CORE_SRC_MANAGERS_EVENTS_PERIODICEVENT_HPP_
#define SRC_COMPONENTS_CORE_SRC_MANAGERS_EVENTS_PERIODICEVENT_HPP_

#include "managers/Events/EventBase.hpp"

namespace geosx
{

class PeriodicEvent : public EventBase
{
public:
  PeriodicEvent(const std::string& name,
                ManagedGroup * const parent);
  
  virtual ~PeriodicEvent();

  static string CatalogName() { return "PeriodicEvent"; }

  virtual void FillDocumentationNode() override;
  
  virtual void EstimateEventTiming(real64 const time,
                                   real64 const dt, 
                                   integer const cycle,
                                   dataRepository::ManagedGroup * domain) override;

  void CheckOptionalFunctionThreshold(real64 const time,
                                      real64 const dt, 
                                      integer const cycle,
                                      dataRepository::ManagedGroup * domain);

  virtual real64 GetTimestepRequest(real64 const time) override;

  dataRepository::ManagedGroup * m_functionTarget;

  struct viewKeyStruct
  {
    dataRepository::ViewKey timeFrequency = { "timeFrequency" };
    dataRepository::ViewKey cycleFrequency = { "cycleFrequency" };
    dataRepository::ViewKey targetExactTimestep = { "targetExactTimestep" };
    dataRepository::ViewKey functionName = { "function" };
    dataRepository::ViewKey functionInputObject = { "object" };
    dataRepository::ViewKey functionInputSetname = { "set" };
    dataRepository::ViewKey functionSetNames = { "setNames" };
    dataRepository::ViewKey functionStatOption = { "stat" };
    dataRepository::ViewKey eventThreshold = { "threshold" };
  } viewKeys;

};

} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_MANAGERS_EVENTS_PERIODICEVENT_HPP_ */
