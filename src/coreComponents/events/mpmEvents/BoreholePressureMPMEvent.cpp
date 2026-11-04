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
 * @file BoreholePressureMPMEvent.cpp
 */

#include "BoreholePressureMPMEvent.hpp"
#include <cmath>

namespace geos
{

using namespace dataRepository;

BoreholePressureMPMEvent::BoreholePressureMPMEvent( const string & name,
                                                    Group * const parent ):
  MPMEventBase( name, parent ),
  m_boreholeRadius( 0.0 ),
  m_startPressure( 0.0 ),
  m_endPressure( 0.0 ),
  m_interpType( 1 )
{
  registerWrapper( viewKeyStruct::boreholeRadiusString(), &m_boreholeRadius ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Positive pressure-region cutoff radius: sqrt(x*x+y*y) in Cartesian, r in Axisymmetric (independent of axial y)" );

  registerWrapper( viewKeyStruct::startPressureString(), &m_startPressure ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Starting pressure to ramp from" );

  registerWrapper( viewKeyStruct::endPressureString(), &m_endPressure ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "End pressure to ramp to" );

  registerWrapper( viewKeyStruct::interpTypeString(), &m_interpType ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_interpType ).
    setDescription( "Interpolation scheme: 0 (Linear), 1 (Cosine), 2 (Smooth-step)" );     // CC: double check this!
}

BoreholePressureMPMEvent::~BoreholePressureMPMEvent()
{}

void BoreholePressureMPMEvent::postInputInitialization()
{
  MPMEventBase::postInputInitialization();

  GEOS_ERROR_IF( !std::isfinite(m_boreholeRadius) || m_boreholeRadius <= 0.0,
                 "Borehole radius must be finite and positive" );
  GEOS_ERROR_IF( !std::isfinite(m_startPressure) || !std::isfinite(m_endPressure),
                 "Borehole pressure endpoints must be finite" );
  GEOS_ERROR_IF( m_interpType < 0 || m_interpType > 2, "Pressure interpType must be 0, 1 or 2" );
}

REGISTER_CATALOG_ENTRY( MPMEventBase, BoreholePressureMPMEvent, string const &, Group * const )

} /* namespace geos */
