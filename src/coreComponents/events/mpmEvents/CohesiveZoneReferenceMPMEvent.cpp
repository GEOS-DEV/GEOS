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
 * @file CohesiveZoneReferenceMPMEvent.cpp
 */

#include "CohesiveZoneReferenceMPMEvent.hpp"

namespace geos
{

  using namespace dataRepository;
  
  CohesiveZoneReferenceMPMEvent::CohesiveZoneReferenceMPMEvent( const string & name,
                                                                Group * const parent ) :
                                                                MPMEventBase(  name, parent ),
                                                                m_czVolumeNormalization( 1 )
  {  
    registerWrapper( "czVolumeNormalization", &m_czVolumeNormalization ).
      setInputFlag( InputFlags::OPTIONAL ).
      setApplyDefaultValue( m_czVolumeNormalization ).
      setRestartFlags( RestartFlags::WRITE_AND_READ ).
      setDescription( "Flag to perform volume normalization of cohesive nodal area for partially filled grid cells" );
  }

  CohesiveZoneReferenceMPMEvent::~CohesiveZoneReferenceMPMEvent() 
  {}

  void CohesiveZoneReferenceMPMEvent::postInputInitialization()
  {
    GEOS_LOG_RANK_0( "CohesiveZoneReferenceEvent: " << 
                     "Start time=" << m_startTime << ", " << 
                     "Time interval=" << getTimeInterval() << ", " << 
                     "czVolumeNormalization=" << m_czVolumeNormalization );
  }

  REGISTER_CATALOG_ENTRY( MPMEventBase, CohesiveZoneReferenceMPMEvent, string const &, Group * const )

} /* namespace geos */
