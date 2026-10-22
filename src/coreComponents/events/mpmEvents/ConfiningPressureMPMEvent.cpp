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
 * @file ConfiningPressureMPMEvent.cpp
 */

#include "ConfiningPressureMPMEvent.hpp"

namespace geos
{

using namespace dataRepository;

ConfiningPressureMPMEvent::ConfiningPressureMPMEvent( const string & name,
                                                      Group * const parent ):
  MPMEventBase( name, parent ),
  m_confiningPressureBoxMin(  ),
  m_confiningPressureBoxMax(  ),
  m_confiningPressureSphereRadius( 0.0 ),
  m_uniformPressure( 0 ),
  m_startPressure( 0.0 ),
  m_endPressure( 0.0 ),
  m_interpType( 1 )
{
  registerWrapper( viewKeyStruct::confiningPressureBoxMinString(), &m_confiningPressureBoxMin ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Min corner of box defining confining pressure region. Required unless confiningPressureSphereRadius is positive." );

  registerWrapper( viewKeyStruct::confiningPressureBoxMaxString(), &m_confiningPressureBoxMax ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Max corner of box defining confining pressure region. Required unless confiningPressureSphereRadius is positive." );

  registerWrapper( viewKeyStruct::confiningPressureSphereRadiusString(), &m_confiningPressureSphereRadius ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_confiningPressureSphereRadius ).
    setDescription( "Optional origin-centered spherical confining pressure region radius. If positive, nodes outside this radius receive the confining background stress instead of using the box test." );

  // Uniform pressure needs no surface reconstruction or geometric exclusion.
  // With one velocity field, the background-stress weak form loads the external
  // boundary; contributions on conforming internal material interfaces cancel.
  registerWrapper( viewKeyStruct::uniformPressureString(), &m_uniformPressure ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_uniformPressure ).
    setDescription( "Apply pressure everywhere (0/1), without box or sphere exclusion. Required for Axisymmetric." );

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

ConfiningPressureMPMEvent::~ConfiningPressureMPMEvent()
{}

void ConfiningPressureMPMEvent::postInputInitialization()
{
  MPMEventBase::postInputInitialization();

  GEOS_ERROR_IF( m_confiningPressureSphereRadius < 0.0,
                 "confiningPressureSphereRadius must be non-negative." );

  GEOS_ERROR_IF( m_uniformPressure < 0 || m_uniformPressure > 1,
                 "uniformPressure must be 0 or 1." );
  GEOS_ERROR_IF( m_uniformPressure && (m_confiningPressureBoxMin.size() > 0 ||
                  m_confiningPressureBoxMax.size() > 0 || m_confiningPressureSphereRadius > 0.0),
                 "uniformPressure=1 must not specify a confining-pressure box or sphere." );

  bool const useSphere = m_confiningPressureSphereRadius > 0.0;
  bool const hasBoxMin = m_confiningPressureBoxMin.size() > 0;
  bool const hasBoxMax = m_confiningPressureBoxMax.size() > 0;
  GEOS_ERROR_IF( !m_uniformPressure && !useSphere && !( hasBoxMin && hasBoxMax ),
                 "ConfiningPressure requires uniformPressure=1, a positive sphere radius, or both box corners." );
  GEOS_ERROR_IF( hasBoxMin && m_confiningPressureBoxMin.size() != 3,
                 "confiningPressureBoxMin must be of length 3." );
  GEOS_ERROR_IF( hasBoxMax && m_confiningPressureBoxMax.size() != 3,
                 "confiningPressureBoxMax must be of length 3." );
  GEOS_ERROR_IF( !m_uniformPressure && !useSphere && m_confiningPressureBoxMin.size() != 3,
                 "confiningPressureBoxMin must be of length 3 when confiningPressureSphereRadius is not positive." );
  GEOS_ERROR_IF( !m_uniformPressure && !useSphere && m_confiningPressureBoxMax.size() != 3,
                 "confiningPressureBoxMax must be of length 3 when confiningPressureSphereRadius is not positive." );

  GEOS_LOG_RANK_0( "ConfiningPressureEvent: " <<
                   "Start time=" << m_startTime << ", " <<
                   "Time interval=" << getTimeInterval() << ", " <<
                   "confiningPressureBoxMin=" << m_confiningPressureBoxMin << ", " <<
                   "confiningPressureBoxMax=" << m_confiningPressureBoxMax << ", " <<
                   "confiningPressureSphereRadius=" << m_confiningPressureSphereRadius << ", " <<
                   "uniformPressure=" << m_uniformPressure << ", " <<
                   "startPressure=" << m_startPressure << ", " <<
                   "endPressure=" << m_endPressure << ", " <<
                   "interpType=" << m_interpType );
}

REGISTER_CATALOG_ENTRY( MPMEventBase, ConfiningPressureMPMEvent, string const &, Group * const )

} /* namespace geos */
