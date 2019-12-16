/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
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
using namespace cxx_utilities;

namespace constitutive
{


VanGenuchtenCapillaryPressure::VanGenuchtenCapillaryPressure( std::string const & name,
                                                            Group * const parent )
  : CapillaryPressureBase( name, parent )
{
  registerWrapper( viewKeyStruct::phaseMinVolumeFractionString, &m_phaseMinVolumeFraction, false )->
    setApplyDefaultValue(0.0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Minimum volume fraction value for each phase");

  registerWrapper( viewKeyStruct::phaseCapPressureExponentInvString,   &m_phaseCapPressureExponentInv,   false )->
    setApplyDefaultValue(0.5)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Inverse of capillary power law exponent for each phase");

  registerWrapper( viewKeyStruct::phaseCapPressureMultiplierString,   &m_phaseCapPressureMultiplier,   false )->
    setApplyDefaultValue(1.0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Entry pressure value for each phase");

  registerWrapper( viewKeyStruct::capPressureEpsilonString,   &m_capPressureEpsilon,   false )->
    setApplyDefaultValue(1e-6)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Saturation at which the extremum capillary pressure is attained; used to avoid infinite capillary pressure values for saturations close to 0 and 1");
}

VanGenuchtenCapillaryPressure::~VanGenuchtenCapillaryPressure()
{

}

void
VanGenuchtenCapillaryPressure::DeliverClone( string const & name,
                                             Group * const parent,
                                             std::unique_ptr<ConstitutiveBase> & clone ) const
{
  std::unique_ptr< VanGenuchtenCapillaryPressure > newModel = std::make_unique<VanGenuchtenCapillaryPressure>( name, parent );

  newModel->m_phaseNames = this->m_phaseNames;
  newModel->m_phaseTypes = this->m_phaseTypes;
  newModel->m_phaseOrder = this->m_phaseOrder;
  
  newModel->m_phaseMinVolumeFraction      = this->m_phaseMinVolumeFraction;
  newModel->m_phaseCapPressureExponentInv = this->m_phaseCapPressureExponentInv;
  newModel->m_phaseCapPressureMultiplier  = this->m_phaseCapPressureMultiplier;

  newModel->m_capPressureEpsilon = this->m_capPressureEpsilon;
  newModel->m_volFracScale       = this->m_volFracScale;

  clone = std::move( newModel );
}


void VanGenuchtenCapillaryPressure::PostProcessInput()
{
  CapillaryPressureBase::PostProcessInput();

  localIndex const NP = numFluidPhases();

#define COREY_CHECK_INPUT_LENGTH( data, expected, attr ) \
  if (integer_conversion<localIndex>((data).size()) != integer_conversion<localIndex>(expected)) \
  { \
    GEOSX_ERROR( "VanGenuchtenCapillaryPressure: invalid number of entries in " \
                << (attr) << " attribute (" \
                << (data).size() << "given, " \
                << (expected) << " expected)"); \
  }

  COREY_CHECK_INPUT_LENGTH( m_phaseMinVolumeFraction, NP, viewKeyStruct::phaseMinVolumeFractionString )
  COREY_CHECK_INPUT_LENGTH( m_phaseCapPressureExponentInv, NP, viewKeyStruct::phaseCapPressureExponentInvString )
  COREY_CHECK_INPUT_LENGTH( m_phaseCapPressureMultiplier, NP, viewKeyStruct::phaseCapPressureMultiplierString )

#undef COREY_CHECK_INPUT_LENGTH

  m_volFracScale = 1.0;
  for (localIndex ip = 0; ip < NP; ++ip)
  {
    GEOSX_ERROR_IF( m_phaseMinVolumeFraction[ip] < 0.0 || m_phaseMinVolumeFraction[ip] > 1.0,
                   "VanGenuchtenCapillaryPressure: invalid min volume fraction value: " << m_phaseMinVolumeFraction[ip] );
    m_volFracScale -= m_phaseMinVolumeFraction[ip];

    GEOSX_ERROR_IF(    (m_phaseCapPressureExponentInv[ip] < 0 || m_phaseCapPressureExponentInv[ip] > 1.0)
                   && (m_phaseTypes[ip] != CapillaryPressureBase::REFERENCE_PHASE),
                   "VanGenuchtenCapillaryPressure: invalid exponent inverse value: " << m_phaseCapPressureExponentInv[ip] );

    GEOSX_ERROR_IF(    (m_phaseCapPressureMultiplier[ip] < 0.0)
                   && (m_phaseTypes[ip] != CapillaryPressureBase::REFERENCE_PHASE),
                   "VanGenuchtenCapillaryPressure: invalid entry pressure: " << m_phaseCapPressureMultiplier[ip] );

    GEOSX_ERROR_IF(    (m_capPressureEpsilon < 0.0 || m_capPressureEpsilon > 0.2)
                   && (m_phaseTypes[ip] != CapillaryPressureBase::REFERENCE_PHASE),
                   "VanGenuchtenCapillaryPressure: invalid epsilon: " << m_capPressureEpsilon );

  }

  GEOSX_ERROR_IF( m_volFracScale < 0.0, "VanGenuchtenCapillaryPressure: sum of min volume fractions exceeds 1.0" );
}


void VanGenuchtenCapillaryPressure::BatchUpdate( arrayView2d<real64 const> const & phaseVolumeFraction )
{

  arrayView1d<real64 const> const & phaseMinVolumeFraction      = m_phaseMinVolumeFraction;
  arrayView1d<real64 const> const & phaseCapPressureExponentInv = m_phaseCapPressureExponentInv;
  arrayView1d<real64 const> const & phaseCapPressureMultiplier  = m_phaseCapPressureMultiplier;
  real64 const & capPressureEpsilon = m_capPressureEpsilon;

  CapillaryPressureBase::BatchUpdateKernel<VanGenuchtenCapillaryPressure>( phaseVolumeFraction,
                                                                           phaseMinVolumeFraction,
                                                                           phaseCapPressureExponentInv,
                                                                           phaseCapPressureMultiplier,
                                                                           capPressureEpsilon,
                                                                           m_volFracScale );
}


void VanGenuchtenCapillaryPressure::PointUpdate( arraySlice1d<real64 const> const & phaseVolFraction,
                                                 localIndex const k,
                                                 localIndex const q )
{
  arraySlice1d<real64> const capPressure           = m_phaseCapPressure[k][q];
  arraySlice2d<real64> const dCapPressure_dVolFrac = m_dPhaseCapPressure_dPhaseVolFrac[k][q];

  localIndex const NP = numFluidPhases();

  Compute( NP,
           phaseVolFraction,
           capPressure,
           dCapPressure_dVolFrac,
           m_phaseOrder,
           m_phaseMinVolumeFraction,
           m_phaseCapPressureExponentInv,
           m_phaseCapPressureMultiplier,
           m_capPressureEpsilon,
           m_volFracScale );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, VanGenuchtenCapillaryPressure, std::string const &, Group * const )
} // namespace constitutive

} // namespace geosx
