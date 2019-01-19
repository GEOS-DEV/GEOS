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
 * @file EventManager.cpp
 */

#include "EventManager.hpp"
#include "managers/Events/EventBase.hpp"
#include "common/TimingMacros.hpp"
#include "MPI_Communications/CommunicationTools.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace cxx_utilities;


EventManager::EventManager( std::string const & name,
                            ManagedGroup * const parent ):
  ManagedGroup( name, parent),
  m_maxTime(),
  m_maxCycle(),
  m_verbosity(),
  m_time(),
  m_dt(),
  m_cycle(),
  m_currentSubEvent(),
  m_currentMaxDt()
{
  RegisterViewWrapper(viewKeyStruct::maxTimeString, &m_maxTime, false )->
    setApplyDefaultValue(-1.0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Maximum simulation time.");

  RegisterViewWrapper(viewKeyStruct::maxCycleString, &m_maxCycle, false )->
    setApplyDefaultValue(-1.0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Maximum simulation cycle.");

  RegisterViewWrapper(viewKeyStruct::verbosityString, &m_verbosity, false )->
    setApplyDefaultValue(0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Maximum simulation time.");


  RegisterViewWrapper(viewKeyStruct::timeString, &m_time, false )->
    setRestartFlags(RestartFlags::WRITE_AND_READ)->
    setDescription("Current simulation time.");

  RegisterViewWrapper(viewKeyStruct::dtString, &m_dt, false )->
    setRestartFlags(RestartFlags::WRITE_AND_READ)->
    setDescription("Current simulation timestep.");

  RegisterViewWrapper(viewKeyStruct::cycleString, &m_cycle, false )->
    setRestartFlags(RestartFlags::WRITE_AND_READ)->
    setDescription("Current simulation cycle number.");

  RegisterViewWrapper(viewKeyStruct::currentSubEventString, &m_currentSubEvent, false )->
    setRestartFlags(RestartFlags::WRITE_AND_READ)->
    setDescription("index of the current subevent.");

  RegisterViewWrapper(viewKeyStruct::currentMaxDtString, &m_currentMaxDt, false )->
    setRestartFlags(RestartFlags::WRITE_AND_READ)->
    setDescription("Maximum dt request for event loop.");


}


EventManager::~EventManager()
{}




void EventManager::PostProcessInput()
{
  real64 & maxTime = this->getReference<real64>(viewKeys.maxTime);
  integer & maxCycle = this->getReference<integer>(viewKeys.maxCycle);

  // If maxTime, maxCycle are default, set them to their max values
  if (maxTime < 0)
  {
    maxTime = std::numeric_limits<real64>::max();
  }
  if (maxCycle < 0)
  {
    maxCycle = std::numeric_limits<integer>::max();
  }
}


ManagedGroup * EventManager::CreateChild( string const & childKey, string const & childName )
{
  GEOS_LOG_RANK_0("Adding Event: " << childKey << ", " << childName);
  std::unique_ptr<EventBase> event = EventBase::CatalogInterface::Factory( childKey, childName, this );
  return this->RegisterGroup<EventBase>( childName, std::move(event) );
}


void EventManager::Run(dataRepository::ManagedGroup * domain)
{
  GEOSX_MARK_FUNCTION;
  real64& time = this->getReference<real64>(viewKeys.time);
  real64& dt = this->getReference<real64>(viewKeys.dt);
  integer& cycle = this->getReference<integer>(viewKeys.cycle);
  integer& currentSubEvent = this->getReference<integer>(viewKeys.currentSubEvent);
  real64& currentMaxDt = this->getReference<real64>(viewKeys.currentMaxDt);
  real64 const maxTime = this->getReference<real64>(viewKeys.maxTime);
  integer const maxCycle = this->getReference<integer>(viewKeys.maxCycle);
  integer const verbosity = this->getReference<integer>(viewKeys.verbosity);
  integer exitFlag = 0;

  // Setup MPI communication
  integer const rank = CommunicationTools::MPI_Rank(MPI_COMM_GEOSX );
  integer const comm_size = CommunicationTools::MPI_Size(MPI_COMM_GEOSX );;

  real64 send_buffer[2];
  array1d<real64> receive_buffer(2 * comm_size);

  // Setup event targets, sequence indicators
  array1d<integer> eventCounters(2);
  this->forSubGroups<EventBase>([&]( EventBase * subEvent ) -> void
  {
    subEvent->GetTargetReferences();
    subEvent->GetExecutionOrder(eventCounters);
  });

  // Set the progress indicators
  this->forSubGroups<EventBase>([&]( EventBase * subEvent ) -> void
  {
    subEvent->SetProgressIndicator(eventCounters);
  });

  // Inform user if it appears this is a mid-loop restart
  if ((currentSubEvent > 0))
  {
    GEOS_LOG_RANK_0("The restart-file was written during step " << currentSubEvent << " of the event loop.  Resuming from that point.");
  }

  // Run problem
  while((time < maxTime) && (cycle < maxCycle) && (exitFlag == 0))
  {
    GEOS_LOG_RANK_0("Time: " << time << "s, dt:" << dt << "s, Cycle: " << cycle);

    // Iterage using the managed integer currentSubEvent and real64 currentMaxDt,
    // which will allow restart runs to pick up where they left off.
    if (currentSubEvent == 0)
    {
      currentMaxDt = std::numeric_limits<real64>::max();
    }

    for ( ; currentSubEvent<this->numSubGroups(); ++currentSubEvent)
    {
      EventBase * subEvent = static_cast<EventBase *>( this->GetSubGroups()[currentSubEvent] );

      // Calculate the event and sub-event forecasts
      // Note: because events can be nested, the mpi reduce for event
      // forecasts need to happen in EventBase.
      subEvent->CheckEvents(time, dt, cycle, domain);
      integer eventForecast = subEvent->GetForecast();

      // Execute, signal events
      if (eventForecast == 1)
      {
        subEvent->SignalToPrepareForExecution(time, dt, cycle, domain);
      }

      if (eventForecast <= 0)
      {
        subEvent->Execute(time, dt, cycle, 0, 0, domain);
      }

      // Estimate the time-step for the next cycle
      if (eventForecast <= 1)
      {
        real64 requestedDt = subEvent->GetTimestepRequest(time + dt);
        currentMaxDt = std::min(requestedDt, currentMaxDt);
      }

      // Check the exit flag
      exitFlag += subEvent->GetExitFlag();
    
      // Debug information
      if (verbosity > 0)
      {
        GEOS_LOG_RANK_0("     Event: " << currentSubEvent << " (" << subEvent->getName() << "), f=" << eventForecast);
      }
    }

    // Increment the time, cycle
    time += dt;
    ++cycle;
    dt = currentMaxDt;
    dt = (time + dt > maxTime) ? (maxTime - time) : dt;

    // Reset the subevent counter
    currentSubEvent = 0;

#ifdef GEOSX_USE_MPI
//    MPI_Barrier(MPI_COMM_GEOSX);
    GEOSX_MARK_BEGIN("EventManager::MPI calls");

    send_buffer[0] = dt;
    send_buffer[1] = static_cast<real64>(exitFlag);
    MPI_Gather(send_buffer, 2, MPI_DOUBLE, receive_buffer.data(), 2, MPI_DOUBLE, 0, MPI_COMM_GEOSX);

    if (rank == 0)
    {
      for (integer ii=0; ii<comm_size; ii++)
      {
        send_buffer[0] = std::min(send_buffer[0], receive_buffer[2*ii]);
        send_buffer[1] += receive_buffer[2*ii + 1];
      }
    }

    MPI_Bcast(send_buffer, 2, MPI_DOUBLE, 0, MPI_COMM_GEOSX);
    dt = send_buffer[0];
    if (send_buffer[1] > 0.5)
    {
      exitFlag = 1;
    }
    GEOSX_MARK_END("EventManager::MPI calls");

#endif
  }

  // Cleanup
  GEOS_LOG_RANK_0("Cleaning up events");
  
  this->forSubGroups<EventBase>([&]( EventBase * subEvent ) -> void
  {
    subEvent->Cleanup(time, cycle, 0, 0, domain);
  });
}

} /* namespace geosx */
