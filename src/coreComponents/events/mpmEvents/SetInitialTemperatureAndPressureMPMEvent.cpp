/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file SetInitialTemperatureAndPressureMPMEvent.cpp
 */

#include "SetInitialTemperatureAndPressureMPMEvent.hpp"

namespace geos
{

using namespace dataRepository;

SetInitialTemperatureAndPressureMPMEvent::SetInitialTemperatureAndPressureMPMEvent( const string & name,
                                                                                    Group * const parent ):
  MPMEventBase( name, parent ),
  m_pressure( 0.0 ),
  m_temperature( 300.0 ),
  m_targetRegions(),
  m_initialPhase( "auto" ),
  m_initializeStress( 1 )
{
  registerWrapper( viewKeyStruct::pressureString(), &m_pressure ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Target positive-compressive pressure for EOS-consistent initialization." );

  registerWrapper( viewKeyStruct::temperatureString(), &m_temperature ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Target temperature for EOS-consistent initialization." );

  registerWrapper( viewKeyStruct::targetRegionsString(), &m_targetRegions ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Particle material regions to initialize. Use all to initialize every particle region." );

  registerWrapper( viewKeyStruct::initialPhaseString(), &m_initialPhase ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( string( "auto" ) ).
    setDescription( "Optional material-specific phase selector passed to the EOS helper, e.g. auto, solid, or liquid." );

  registerWrapper( viewKeyStruct::initializeStressString(), &m_initializeStress ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 1 ).
    setDescription( "If nonzero, initialize particle and constitutive stress to the returned hydrostatic pressure." );
}

SetInitialTemperatureAndPressureMPMEvent::~SetInitialTemperatureAndPressureMPMEvent()
{}

void SetInitialTemperatureAndPressureMPMEvent::postInputInitialization()
{
  MPMEventBase::postInputInitialization();

  GEOS_ERROR_IF( m_targetRegions.empty(),
                 "SetInitialTemperatureAndPressure requires at least one target region." );
  GEOS_ERROR_IF( m_temperature < 0.0,
                 "SetInitialTemperatureAndPressure requires nonnegative temperature." );
  GEOS_ERROR_IF( !( m_initializeStress == 0 || m_initializeStress == 1 ),
                 "SetInitialTemperatureAndPressure initializeStress must be 0 or 1." );

  GEOS_LOG_RANK_0( "SetInitialTemperatureAndPressureEvent: "
                   << "Start time=" << m_startTime << ", "
                   << "Time interval=" << getTimeInterval() << ", "
                   << "pressure=" << m_pressure << ", "
                   << "temperature=" << m_temperature << ", "
                   << "initialPhase=" << m_initialPhase );
}

REGISTER_CATALOG_ENTRY( MPMEventBase, SetInitialTemperatureAndPressureMPMEvent, string const &, Group * const )

} /* namespace geos */
