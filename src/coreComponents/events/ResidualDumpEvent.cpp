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

#include "ResidualDumpEvent.hpp"

namespace geos
{

using namespace dataRepository;

ResidualDumpEvent::ResidualDumpEvent( const string & name, Group * const parent ):
  PeriodicEvent( name, parent ),
  m_secondaryEventTarget( "" ),
  m_secondaryTarget( nullptr ),
  m_flagSecondaryTrigger()
{
  registerWrapper( viewKeyStruct::secondaryTargetString(), &m_secondaryEventTarget ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Name of the object to be executed when the secondary event criteria are met." );

}

void ResidualDumpEvent::getTargetReferences()
{
  EventBase::getTargetReferences();
  if( !m_secondaryEventTarget.empty() )
  {
    try
    {
      m_secondaryTarget = &this->getGroupByPath< ExecutableGroup >( m_secondaryEventTarget );
    }
    catch( std::exception const & e )
    {
      throw InputError( e, GEOS_FMT( "Error while reading {}:\n",
                                     getWrapperDataContext( viewKeyStruct::secondaryTargetString() ) ) );
    }
  }

  this->forSubGroups< EventBase >( []( EventBase & subEvent ){
    subEvent.getTargetReferences();
  } );
}

void ResidualDumpEvent::validate() const
{

  if( m_secondaryTarget == nullptr )
    return;

  PeriodicEvent::validate();

}

bool ResidualDumpEvent::execute( const geos::real64 time_n, const geos::real64 dt, const geos::integer cycleNumber,
                                 const geos::integer, const geos::real64, geos::DomainPartition & domain )
{

  bool earlyReturn = false;

  ExecutableGroup * target = getEventTarget();
  // If m_targetExecFlag is set, then the code has resumed at a point
  // after the target has executed.
  if(( target != nullptr) && (m_targetExecFlag == 0))
  {
    m_targetExecFlag = 1;
    earlyReturn = earlyReturn ||
                  target->execute( time_n, dt, cycleNumber, m_eventCount, m_eventProgress, domain );

    m_flagSecondaryTrigger = dynamicCast< SolverBase * >( target )->getRootFlag();
    if( m_flagSecondaryTrigger->isAnySet( SolverBase::SolverGroupFlags::StruggleCvg ))
      m_secondaryTarget->execute( time_n, dt, cycleNumber, m_eventCount, m_eventProgress, domain );
  }

  // Iterate through the sub-event list using the managed integer m_currentSubEvent
  // This allows for  restart runs to pick up where they left off.
  for(; m_currentSubEvent < this->numSubGroups(); ++m_currentSubEvent )
  {
    EventBase * subEvent = static_cast< EventBase * >( this->getSubGroups()[m_currentSubEvent] );

    // Print debug information for logLevel >= 1
    GEOS_LOG_LEVEL_RANK_0( 1,
                           "          SubEvent: " << m_currentSubEvent << " (" << subEvent->getName() << "), dt_request=" << subEvent->getCurrentEventDtRequest() << ", forecast=" <<
                           subEvent->getForecast() );

    if( subEvent->isReadyForExec() )
    {
      earlyReturn = earlyReturn ||
                    subEvent->execute( time_n, dt, cycleNumber, m_eventCount, m_eventProgress, domain );

      m_flagSecondaryTrigger = dynamicCast< SolverBase * >( target )->getRootFlag();
      if( m_flagSecondaryTrigger->isAnySet( SolverBase::SolverGroupFlags::StruggleCvg ))
        m_secondaryTarget->execute( time_n, dt, cycleNumber, m_eventCount, m_eventProgress, domain );
    }
  }

  // Update the event status
  m_targetExecFlag =  0;
  m_currentSubEvent = 0;
  m_lastTime = time_n;
  m_lastCycle = cycleNumber;

  return earlyReturn;

}


REGISTER_CATALOG_ENTRY( EventBase, ResidualDumpEvent, string const &, Group * const )
} // geos
