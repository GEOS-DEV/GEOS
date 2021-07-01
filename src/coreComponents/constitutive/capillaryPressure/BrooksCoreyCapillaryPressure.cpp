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
 * @file BrooksCoreyCapillaryPressure.cpp
 */

#include "BrooksCoreyCapillaryPressure.hpp"

namespace geosx
{

using namespace dataRepository;

namespace constitutive
{


BrooksCoreyCapillaryPressure::BrooksCoreyCapillaryPressure( string const & name,
                                                            Group * const parent )
  : CapillaryPressureBase( name, parent )
{
  registerWrapper( viewKeyStruct::phaseMinVolumeFractionString(), &m_phaseMinVolumeFraction ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Minimum volume fraction value for each phase" );

  registerWrapper( viewKeyStruct::phaseCapPressureExponentInvString(), &m_phaseCapPressureExponentInv ).
    setApplyDefaultValue( 2.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Inverse of capillary power law exponent for each phase" );

  registerWrapper( viewKeyStruct::phaseEntryPressureString(), &m_phaseEntryPressure ).
    setApplyDefaultValue( 1.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Entry pressure value for each phase" );

  registerWrapper( viewKeyStruct::capPressureEpsilonString(), &m_capPressureEpsilon ).
    setApplyDefaultValue( 1e-6 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription(
    "Wetting-phase saturation at which the max cap. pressure is attained; used to avoid infinite cap. pressure values for saturations close to zero" );

  registerWrapper( viewKeyStruct::volFracScaleString(), &m_volFracScale ).
    setApplyDefaultValue( 1.0 ).
    setDescription( "Factor used to scale the phase capillary pressure, defined as: one minus the sum of the phase minimum volume fractions." );
}

BrooksCoreyCapillaryPressure::~BrooksCoreyCapillaryPressure()
{}

namespace
{

template< typename ARRAY >
void checkInputSize( ARRAY const & array, localIndex const expected, string const & attr )
{
  GEOSX_THROW_IF_NE_MSG( array.size(), expected,
                         "BrooksCoreyCapillaryPressure: invalid number of entries in " << attr << " attribute",
                         InputError );

}

}

void BrooksCoreyCapillaryPressure::postProcessInput()
{
  CapillaryPressureBase::postProcessInput();

  localIndex const NP = numFluidPhases();

  checkInputSize( m_phaseMinVolumeFraction, NP, viewKeyStruct::phaseMinVolumeFractionString() );
  checkInputSize( m_phaseCapPressureExponentInv, NP, viewKeyStruct::phaseCapPressureExponentInvString() );
  checkInputSize( m_phaseEntryPressure, NP, viewKeyStruct::phaseEntryPressureString() );

  m_volFracScale = 1.0;
  for( localIndex ip = 0; ip < NP; ++ip )
  {
    GEOSX_THROW_IF( (m_phaseMinVolumeFraction[ip] < 0.0 || m_phaseMinVolumeFraction[ip] > 1.0),
                    "BrooksCoreyCapillaryPressure: invalid min volume fraction value: " << m_phaseMinVolumeFraction[ip],
                    InputError );
    m_volFracScale -= m_phaseMinVolumeFraction[ip];

    GEOSX_THROW_IF( (m_phaseCapPressureExponentInv[ip] < 1.0)
                    && (m_phaseTypes[ip] != CapillaryPressureBase::REFERENCE_PHASE),
                    "BrooksCoreyCapillaryPressure: invalid exponent inverse value: " << m_phaseCapPressureExponentInv[ip],
                    InputError );

    GEOSX_THROW_IF( (m_phaseEntryPressure[ip] < 0.0)
                    && (m_phaseTypes[ip] != CapillaryPressureBase::REFERENCE_PHASE),
                    "BrooksCoreyCapillaryPressure: invalid entry pressure: " << m_phaseEntryPressure[ip],
                    InputError );

    GEOSX_THROW_IF( (m_capPressureEpsilon< 0.0 || m_capPressureEpsilon > 0.2)
                    && (m_phaseTypes[ip] != CapillaryPressureBase::REFERENCE_PHASE),
                    "BrooksCoreyCapillaryPressure: invalid epsilon: " << m_capPressureEpsilon,
                    InputError );

  }

  GEOSX_THROW_IF( m_volFracScale < 0.0, "BrooksCoreyCapillaryPressure: sum of min volume fractions exceeds 1.0", InputError );
}

BrooksCoreyCapillaryPressure::KernelWrapper BrooksCoreyCapillaryPressure::createKernelWrapper()
{
  return KernelWrapper( m_phaseMinVolumeFraction,
                        m_phaseCapPressureExponentInv,
                        m_phaseEntryPressure,
                        m_capPressureEpsilon,
                        m_volFracScale,
                        m_phaseTypes,
                        m_phaseOrder,
                        m_phaseCapPressure,
                        m_dPhaseCapPressure_dPhaseVolFrac );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, BrooksCoreyCapillaryPressure, string const &, Group * const )
} // namespace constitutive

} // namespace geosx
