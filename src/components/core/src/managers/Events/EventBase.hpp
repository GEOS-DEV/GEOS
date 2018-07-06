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
 * EventBase.hpp
 *
 *  Created on: Jul 22, 2017
 *      Author: sherman
 */

#ifndef SRC_COMPONENTS_CORE_SRC_MANAGERS_EVENTS_EVENTSBASE_HPP_
#define SRC_COMPONENTS_CORE_SRC_MANAGERS_EVENTS_EVENTSBASE_HPP_

#include "dataRepository/ManagedGroup.hpp"


namespace geosx
{

class EventBase : public dataRepository::ManagedGroup
{
public:
  explicit EventBase( std::string const & name,
                       ManagedGroup * const parent );

  virtual ~EventBase();

  static string CatalogName() { return "EventBase"; }

  EventBase() = default;
  EventBase( EventBase const & ) = default;
  EventBase( EventBase &&) = default;
  EventBase& operator=( EventBase const & ) = default;
  EventBase& operator=( EventBase&& ) = default;


  virtual void SignalToPrepareForExecution(real64 const time,
                                           real64 const dt,  
                                           integer const cycle,
                                           dataRepository::ManagedGroup * domain);

  virtual void Execute( real64 const & time_n,
                        real64 const & dt,
                        int const cycleNumber,
                        dataRepository::ManagedGroup * domain ) override;

  void Step(real64 const time,
            real64 const dt,  
            integer const cycle,
            dataRepository::ManagedGroup * domain );

  virtual void FillDocumentationNode() override;

  virtual void CreateChild( string const & childKey, string const & childName ) override;

  void GetTargetReferences();

  virtual void CheckEvents(real64 const time, 
                              real64 const dt,
                              integer const cycle,
                              dataRepository::ManagedGroup * domain);

  virtual void EstimateEventTiming(real64 const time, 
                                      real64 const dt,
                                      integer const cycle,
                                      dataRepository::ManagedGroup * domain) = 0;

  void CheckSubEvents(real64 const time, 
                      real64 const dt,
                      integer const cycle,
                      dataRepository::ManagedGroup * domain);

  virtual real64 GetTimestepRequest(real64 const time) override;

  dataRepository::ManagedGroup * m_target;

  struct viewKeyStruct
  {
    dataRepository::ViewKey eventTarget = { "target" };
    dataRepository::ViewKey beginTime = { "beginTime" };
    dataRepository::ViewKey endTime = { "endTime" };
    dataRepository::ViewKey forceDt = { "forceDt" };
    dataRepository::ViewKey lastTime = { "lastTime" };
    dataRepository::ViewKey lastCycle = { "lastCycle" };

    dataRepository::ViewKey allowSuperstep = { "allowSuperstep" };
    dataRepository::ViewKey allowSubstep = { "allowSubstep" };
    dataRepository::ViewKey substepFactor = { "substepFactor" };
  } viewKeys;

  using CatalogInterface = cxx_utilities::CatalogInterface< EventBase, std::string const &, ManagedGroup * const >;
  static CatalogInterface::CatalogType& GetCatalog();

  integer GetForecast(){ return m_eventForecast; }
  void SetForecast(integer forecast){ m_eventForecast = forecast; }

  integer GetExitFlag(){ return m_exitFlag; }
  void SetExitFlag(integer flag){ m_exitFlag = flag; }

private:
  integer m_eventForecast = 0;
  integer m_exitFlag = 0;

};

} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_MANAGERS_EVENTS_EVENTSBASE_HPP_ */
