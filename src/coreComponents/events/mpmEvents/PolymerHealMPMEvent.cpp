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
 * @file PolymerHealMPMEvent.cpp
 */

#include "PolymerHealMPMEvent.hpp"

namespace geos
{

using namespace dataRepository;

PolymerHealMPMEvent::PolymerHealMPMEvent( const string & name,
                                          Group * const parent ):
  MPMEventBase( name, parent ),
  m_targetRegion( "mat1" )
{
  registerWrapper( viewKeyStruct::targetRegionString(), &m_targetRegion ).
    setInputFlag( InputFlags::REQUIRED ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Particle region to perform polymer heal on" );
}

PolymerHealMPMEvent::~PolymerHealMPMEvent()
{}

void PolymerHealMPMEvent::postInputInitialization()
{
  MPMEventBase::postInputInitialization();
}

REGISTER_CATALOG_ENTRY( MPMEventBase, PolymerHealMPMEvent, string const &, Group * const )

} /* namespace geos */
