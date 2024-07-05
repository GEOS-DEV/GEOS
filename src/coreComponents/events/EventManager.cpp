/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file EventManager.cpp
 */

#include "EventManager.hpp"

#include "common/TimingMacros.hpp"
#include "events/EventBase.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "common/Units.hpp"

namespace geos
{

using namespace dataRepository;


EventManager::EventManager( string const & name,
                            Group * const parent ):
  Group( name, parent ),
  m_minTime(),
  m_maxTime(),
  m_maxCycle(),
  m_time(),
  m_dt(),
  m_cycle(),
  m_currentSubEvent(),
  // TODO: default to TimeOutputFormat::full?
  m_timeOutputFormat( TimeOutputFormat::seconds )
{
  setInputFlags( InputFlags::REQUIRED );

  // This enables logLevel filtering
  enableLogLevelInput();

  registerWrapper( viewKeyStruct::minTimeString(), &m_minTime ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Start simulation time for the global event loop." );

  registerWrapper( viewKeyStruct::maxTimeString(), &m_maxTime ).
    setApplyDefaultValue( std::numeric_limits< real64 >::max() ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Maximum simulation time for the global event loop. Disabled by default." );

  registerWrapper( viewKeyStruct::maxCycleString(), &m_maxCycle ).
    setApplyDefaultValue( std::numeric_limits< integer >::max() ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Maximum simulation cycle for the global event loop. Disabled by default." );

  registerWrapper( viewKeyStruct::timeString(), &m_time ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Current simulation time." );

  registerWrapper( viewKeyStruct::dtString(), &m_dt ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Current simulation timestep." );

  registerWrapper( viewKeyStruct::cycleString(), &m_cycle ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Current simulation cycle number." );

  registerWrapper( viewKeyStruct::currentSubEventString(), &m_currentSubEvent ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Index of the current subevent." );

  registerWrapper( viewKeyStruct::timeOutputFormat(), &m_timeOutputFormat ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Format of the time in the GEOS log." );
}


EventManager::~EventManager()
{}



Group * EventManager::createChild( string const & childKey, string const & childName )
{
  GEOS_LOG_RANK_0( "Adding Event: " << childKey << ", " << childName );
  std::unique_ptr< EventBase > event = EventBase::CatalogInterface::factory( childKey, childName, this );
  return &this->registerGroup< EventBase >( childName, std::move( event ) );
}


void EventManager::expandObjectCatalogs()
{
  // During schema generation, register one of each type derived from EventBase here
  for( auto & catalogIter: EventBase::getCatalog())
  {
    createChild( catalogIter.first, catalogIter.first );
  }
}


bool EventManager::run( DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  integer exitFlag = 0;

  // Setup event targets, sequence indicators
  array1d< integer > eventCounters( 2 );
  this->forSubGroups< EventBase >( [&]( EventBase & subEvent )
  {
    subEvent.getTargetReferences();
    subEvent.getExecutionOrder( eventCounters );
    subEvent.validate();
  } );

  // Set the progress indicators
  this->forSubGroups< EventBase >( [&]( EventBase & subEvent )
  {
    subEvent.setProgressIndicator( eventCounters );
  } );

  // Inform user if it appears this is a mid-loop restart
  if( m_currentSubEvent > 0 )
  {
    GEOS_LOG_RANK_0( "Resuming from step " << m_currentSubEvent << " of the event loop." );
  }
  else if( !isZero( m_minTime ) )
  {
    // the user has requested a "non-standard" min time (negative time possible for initialization events for instance)
    // since we are not doing a mid-loop restart, we can set the current time to the min time
    // note: commenting out "if( !isZero( m_minTime ) )" will not break the code, but will break the contactMechanics integrated tests
    // restart
    //       because it is done in an unusual fashion
    m_time = m_minTime;
  }

  // Run problem
  // Note: if currentSubEvent > 0, then we are resuming from a restart file
  while((m_time < m_maxTime) && (m_cycle < m_maxCycle) && (exitFlag == 0))
  {
    // Determine the cycle timestep
    if( m_currentSubEvent == 0 )
    {
      // The max dt request
      m_dt = m_maxTime - m_time;

      // Determine the dt requests for each event
      for(; m_currentSubEvent<this->numSubGroups(); ++m_currentSubEvent )
      {
        EventBase * subEvent = static_cast< EventBase * >( this->getSubGroups()[m_currentSubEvent] );
        m_dt = std::min( subEvent->getTimestepRequest( m_time ), m_dt );
      }
      m_currentSubEvent = 0;

#ifdef GEOSX_USE_MPI
      // Find the min dt across processes
      real64 dt_global;
      MPI_Allreduce( &m_dt, &dt_global, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_GEOSX );
      m_dt = dt_global;
#endif
    }
    Section section;
    outputTime( section );
    section.begin();

    // Execute
    for(; m_currentSubEvent<this->numSubGroups(); ++m_currentSubEvent )
    {
      EventBase * subEvent = static_cast< EventBase * >( this->getSubGroups()[m_currentSubEvent] );

      // Calculate the event and sub-event forecasts
      subEvent->checkEvents( m_time, m_dt, m_cycle, domain );

      // Print debug information for logLevel >= 1
      GEOS_LOG_LEVEL_RANK_0( 1, GEOS_FMT( "Event: {} ({}), dt_request={}, forecast={}",
                                          m_currentSubEvent, subEvent->getName(), subEvent->getCurrentEventDtRequest(), subEvent->getForecast() ) );

      // Execute, signal events
      bool earlyReturn = false;
      if( subEvent->hasToPrepareForExec() )
      {
        subEvent->signalToPrepareForExecution( m_time, m_dt, m_cycle, domain );
      }
      else if( subEvent->isReadyForExec() )
      {
        earlyReturn = subEvent->execute( m_time, m_dt, m_cycle, 0, 0, domain );
      }

      // Check the exit flag
      // Note: Currently, this is only being used by the HaltEvent
      //       If it starts being used elsewhere it may need to be synchronized
      exitFlag += subEvent->getExitFlag();

      if( earlyReturn )
      {
        ++m_currentSubEvent;
        return true;
      }
    }

    section.end();
    // Increment time/cycle, reset the subevent counter
    m_time += m_dt;
    ++m_cycle;
    m_currentSubEvent = 0;
  }

  // Cleanup
  GEOS_LOG_RANK_0( "Cleaning up events" );

  this->forSubGroups< EventBase >( [&]( EventBase & subEvent )
  {
    subEvent.cleanup( m_time, m_cycle, 0, 0, domain );
  } );

  return false;
}

void EventManager::outputTime( Section & section ) const
{
  const bool isTimeLimited = m_maxTime < std::numeric_limits< real64 >::max();
  const bool isCycleLimited = m_maxCycle < std::numeric_limits< integer >::max();
  units::TimeFormatInfo const timeInfo = units::TimeFormatInfo::fromSeconds( m_time );
  units::TimeFormatInfo const maxTimeInfo = units::TimeFormatInfo::fromSeconds( m_maxTime );

  const auto timeCompletionUnfoldedString = [&]() -> std::string {
    return GEOS_FMT( " out of {} ({:.0f}% completed)",
                     maxTimeInfo.toUnfoldedString(),
                     100.0 * (m_time - m_minTime) / ( m_maxTime - m_minTime ) );
  };
  const auto timeCompletionSecondsString = [&]() -> std::string {
    return GEOS_FMT( " / {}", maxTimeInfo.toSecondsString() );
  };
  const auto cycleCompletionString = [&]() -> std::string {
    return GEOS_FMT( " out of {} ({:.0f}% completed)",
                     m_maxCycle, ( 100.0 * m_cycle ) / m_maxCycle );
  };

  string const timeCompletionUnfolded = isTimeLimited ? timeCompletionUnfoldedString() : "";
  string const timeCompletionSecond = isTimeLimited ? timeCompletionSecondsString() : "";
  string const cycleLimited = isCycleLimited ? cycleCompletionString() : "";

  string const timeInfosUnfolded = timeInfo.toUnfoldedString() + timeCompletionUnfolded;
  string const timeCompletionSeconds = timeInfo.toSecondsString() + timeCompletionSecond;

  section.setName( "TIMESTEP START" );
  section.addDescription( "Time", timeInfosUnfolded, timeCompletionSeconds );
  section.addDescription( "Delta Time", units::TimeFormatInfo::fromSeconds( m_dt ).toString() );
  section.addDescription( "Cycle", m_cycle, cycleLimited );
  section.setMinWidth(70);

  // GEOS_LOG_RANK_0( "\n------------------------- TIMESTEP START -------------------------" );
  // GEOS_LOG_RANK_0( GEOS_FMT( "    - Time:       {}{}",
  //                            timeInfo.toUnfoldedString(),
  //                            isTimeLimited ? timeCompletionUnfoldedString() : "" ) );
  // GEOS_LOG_RANK_0( GEOS_FMT( "                  ({}{})",
  //                            timeInfo.toSecondsString(),
  //                            isTimeLimited ? timeCompletionSecondsString() : "" ) );
  // GEOS_LOG_RANK_0( GEOS_FMT( "    - Delta Time: {}", units::TimeFormatInfo::fromSeconds( m_dt ) ) );
  // GEOS_LOG_RANK_0( GEOS_FMT( "    - Cycle:      {}{}",
  //                            m_cycle,
  //                            isCycleLimited ? cycleCompletionString() : "" ) );
  // GEOS_LOG_RANK_0( "--------------------------------------------------------------------\n" );

  // We are keeping the old outputs to keep compatibility with current log reading scripts.
  if( m_timeOutputFormat==TimeOutputFormat::full )
  {
    GEOS_LOG_RANK_0( GEOS_FMT( "Time: {} years, {} days, {} hrs, {} min, {} s, dt: {} s, Cycle: {}\n",
                               timeInfo.m_years, timeInfo.m_days, timeInfo.m_hours, timeInfo.m_minutes, timeInfo.m_seconds, m_dt, m_cycle ) );
  }
  else if( m_timeOutputFormat==TimeOutputFormat::years )
  {
    real64 const yearsOut = m_time / units::YearSeconds;
    GEOS_LOG_RANK_0( GEOS_FMT( "Time: {:.2f} years, dt: {} s, Cycle: {}\n", yearsOut, m_dt, m_cycle ) );
  }
  else if( m_timeOutputFormat==TimeOutputFormat::days )
  {
    real64 const daysOut = m_time / units::DaySeconds;
    GEOS_LOG_RANK_0( GEOS_FMT( "Time: {:.2f} days, dt: {} s, Cycle: {}\n", daysOut, m_dt, m_cycle ) );
  }
  else if( m_timeOutputFormat==TimeOutputFormat::hours )
  {
    real64 const hoursOut = m_time / units::HourSeconds;
    GEOS_LOG_RANK_0( GEOS_FMT( "Time: {:.2f} hrs, dt: {} s, Cycle: {}\n", hoursOut, m_dt, m_cycle ) );
  }
  else if( m_timeOutputFormat==TimeOutputFormat::minutes )
  {
    real64 const minutesOut = m_time / units::MinuteSeconds;
    GEOS_LOG_RANK_0( GEOS_FMT( "Time: {:.2f} min, dt: {} s, Cycle: {}\n", minutesOut, m_dt, m_cycle ) );
  }
  else if( m_timeOutputFormat == TimeOutputFormat::seconds )
  {
    GEOS_LOG_RANK_0( GEOS_FMT( "Time: {:4.2e} s, dt: {} s, Cycle: {}\n", m_time, m_dt, m_cycle ) );
  }
  else
  {
    GEOS_ERROR( "Unknown time output format requested." );
  }
}

} /* namespace geos */
