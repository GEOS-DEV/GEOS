/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file VanGenuchtenCapillaryPressure.cpp
 */

#include "VanGenuchtenCapillaryPressure.hpp"

#include <cmath>

namespace geosx
{

using namespace dataRepository;

namespace constitutive
{

VanGenuchtenCapillaryPressure::VanGenuchtenCapillaryPressure( string const & name,
                                                              Group * const parent )
  : CapillaryPressureBase( name, parent )
{
  registerWrapper( viewKeyStruct::phaseMinVolumeFractionString(), &m_phaseMinVolumeFraction ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Minimum volume fraction value for each phase" );

  registerWrapper( viewKeyStruct::phaseCapPressureExponentInvString(), &m_phaseCapPressureExponentInv ).
    setApplyDefaultValue( 0.5 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Inverse of capillary power law exponent for each phase" );

  registerWrapper( viewKeyStruct::phaseCapPressureMultiplierString(), &m_phaseCapPressureMultiplier ).
    setApplyDefaultValue( 1.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Entry pressure value for each phase" );

  registerWrapper( viewKeyStruct::capPressureEpsilonString(), &m_capPressureEpsilon ).
    setApplyDefaultValue( 1e-6 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription(
    "Saturation at which the extremum capillary pressure is attained; used to avoid infinite capillary pressure values for saturations close to 0 and 1" );

  registerWrapper( viewKeyStruct::volFracScaleString(), &m_volFracScale ).
    setApplyDefaultValue( 1.0 ).
    setDescription( "Factor used to scale the phase capillary pressure, defined as: one minus the sum of the phase minimum volume fractions." );

}

VanGenuchtenCapillaryPressure::~VanGenuchtenCapillaryPressure()
{}

namespace
{

template< typename ARRAY >
void checkInputSize( ARRAY const & array, localIndex const expected, string const & attr )
{
  GEOSX_THROW_IF_NE_MSG( array.size(), expected,
                         "VanGenuchtenCapillaryPressure: invalid number of entries in " << attr << " attribute",
                         InputError );

}

}

void VanGenuchtenCapillaryPressure::postProcessInput()
{
  CapillaryPressureBase::postProcessInput();

  localIndex const NP = numFluidPhases();

  checkInputSize( m_phaseMinVolumeFraction, NP, viewKeyStruct::phaseMinVolumeFractionString() );
  checkInputSize( m_phaseCapPressureExponentInv, NP, viewKeyStruct::phaseCapPressureExponentInvString() );
  checkInputSize( m_phaseCapPressureMultiplier, NP, viewKeyStruct::phaseCapPressureMultiplierString() );

  m_volFracScale = 1.0;
  for( localIndex ip = 0; ip < NP; ++ip )
  {
    GEOSX_THROW_IF( (m_phaseMinVolumeFraction[ip] < 0.0) || (m_phaseMinVolumeFraction[ip] > 1.0),
                    "VanGenuchtenCapillaryPressure: invalid min volume fraction value: " << m_phaseMinVolumeFraction[ip],
                    InputError );
    m_volFracScale -= m_phaseMinVolumeFraction[ip];

    GEOSX_THROW_IF( ((m_phaseCapPressureExponentInv[ip] < 0.0) || (m_phaseCapPressureExponentInv[ip] > 1.0))
                    && (m_phaseTypes[ip] != CapillaryPressureBase::REFERENCE_PHASE),
                    "VanGenuchtenCapillaryPressure: invalid exponent inverse value: " << m_phaseCapPressureExponentInv[ip],
                    InputError );

    GEOSX_THROW_IF( (m_phaseCapPressureMultiplier[ip] < 0.0)
                    && (m_phaseTypes[ip] != CapillaryPressureBase::REFERENCE_PHASE),
                    "VanGenuchtenCapillaryPressure: invalid entry pressure: " << m_phaseCapPressureMultiplier[ip],
                    InputError );

    GEOSX_THROW_IF( (m_capPressureEpsilon< 0.0 || m_capPressureEpsilon > 0.2)
                    && (m_phaseTypes[ip] != CapillaryPressureBase::REFERENCE_PHASE),
                    "VanGenuchtenCapillaryPressure: invalid epsilon: " << m_capPressureEpsilon,
                    InputError );

  }

  GEOSX_THROW_IF( m_volFracScale < 0.0, "VanGenuchtenCapillaryPressure: sum of min volume fractions exceeds 1.0", InputError );
}

VanGenuchtenCapillaryPressure::KernelWrapper VanGenuchtenCapillaryPressure::createKernelWrapper()
{
  return KernelWrapper( m_phaseMinVolumeFraction,
                        m_phaseCapPressureExponentInv,
                        m_phaseCapPressureMultiplier,
                        m_capPressureEpsilon,
                        m_volFracScale,
                        m_phaseTypes,
                        m_phaseOrder,
                        m_phaseCapPressure,
                        m_dPhaseCapPressure_dPhaseVolFrac );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, VanGenuchtenCapillaryPressure, string const &, Group * const )
} // namespace constitutive

} // namespace geosx
