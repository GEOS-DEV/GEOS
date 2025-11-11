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
 * @file TemperatureProfileMPMEvent.cpp
 */

#include "TemperatureProfileMPMEvent.hpp"

namespace geos
{

using namespace dataRepository;

TemperatureProfileMPMEvent::TemperatureProfileMPMEvent( const string & name,
                                                        Group * const parent ):
  MPMEventBase( name, parent )
{}

TemperatureProfileMPMEvent::~TemperatureProfileMPMEvent()
{}

void TemperatureProfileMPMEvent::postInputInitialization()
{
  MPMEventBase::postInputInitialization();

  GEOS_LOG_RANK_0( "TemperatureProfileEvent: " <<
                   "Start time=" << m_startTime << ", " <<
                   "Time interval=" << getTimeInterval() );
}

REGISTER_CATALOG_ENTRY( MPMEventBase, TemperatureProfileMPMEvent, string const &, Group * const )

} /* namespace geos */
