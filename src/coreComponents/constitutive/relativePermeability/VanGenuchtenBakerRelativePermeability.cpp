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
 * @file VanGenuchtenBakerRelativePermeability.cpp
 */

#include "VanGenuchtenBakerRelativePermeability.hpp"

#include <cmath>

namespace geosx
{

using namespace dataRepository;

namespace constitutive
{


VanGenuchtenBakerRelativePermeability::VanGenuchtenBakerRelativePermeability( string const & name,
                                                                              Group * const parent )
  : RelativePermeabilityBase( name, parent )
{
  registerWrapper( viewKeyStruct::phaseMinVolumeFractionString(), &m_phaseMinVolumeFraction ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Minimum volume fraction value for each phase" );


  registerWrapper( viewKeyStruct::waterOilRelPermExponentInvString(), &m_waterOilRelPermExponentInv ).
    setApplyDefaultValue( 0.5 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Rel perm power law exponent inverse for the pair (water phase, oil phase) at residual gas saturation\n"
                    "The expected format is \"{ waterExp, oilExp }\", in that order" );

  registerWrapper( viewKeyStruct::waterOilRelPermMaxValueString(), &m_waterOilRelPermMaxValue ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Maximum rel perm value for the pair (water phase, oil phase) at residual gas saturation\n"
                    "The expected format is \"{ waterMax, oilMax }\", in that order" );


  registerWrapper( viewKeyStruct::gasOilRelPermExponentInvString(), &m_gasOilRelPermExponentInv ).
    setApplyDefaultValue( 0.5 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Rel perm power law exponent inverse for the pair (gas phase, oil phase) at residual water saturation\n"
                    "The expected format is \"{ gasExp, oilExp }\", in that order" );

  registerWrapper( viewKeyStruct::gasOilRelPermMaxValueString(), &m_gasOilRelPermMaxValue ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Maximum rel perm value for the pair (gas phase, oil phase) at residual water saturation\n"
                    "The expected format is \"{ gasMax, oilMax }\", in that order" );

  registerWrapper( viewKeyStruct::volFracScaleString(), &m_volFracScale ).
    setApplyDefaultValue( 1.0 ).
    setDescription( "Factor used to scale the phase capillary pressure, defined as: one minus the sum of the phase minimum volume fractions." );

}

VanGenuchtenBakerRelativePermeability::~VanGenuchtenBakerRelativePermeability()
{}

namespace
{

template< typename ARRAY >
void checkInputSize( ARRAY const & array, localIndex const expected, string const & attr )
{
  GEOSX_THROW_IF_NE_MSG( array.size(), expected,
                         "VanGenuchtenBakerRelativePermeability: invalid number of entries in " << attr << " attribute",
                         InputError );

}

}

void VanGenuchtenBakerRelativePermeability::postProcessInput()
{
  RelativePermeabilityBase::postProcessInput();

  localIndex const numPhases = numFluidPhases();

  GEOSX_THROW_IF( m_phaseOrder[PhaseType::OIL] < 0,
                  "VanGenuchtenBakerRelativePermeability: reference oil phase has not been defined and must be included in model",
                  InputError );

  checkInputSize( m_phaseMinVolumeFraction, numPhases, viewKeyStruct::phaseMinVolumeFractionString() );

  if( m_phaseOrder[PhaseType::WATER] >= 0 )
  {
    checkInputSize( m_waterOilRelPermExponentInv, 2, viewKeyStruct::waterOilRelPermExponentInvString() );
    checkInputSize( m_waterOilRelPermMaxValue, 2, viewKeyStruct::waterOilRelPermMaxValueString() );
  }

  if( m_phaseOrder[PhaseType::GAS] >= 0 )
  {
    checkInputSize( m_gasOilRelPermExponentInv, 2, viewKeyStruct::gasOilRelPermExponentInvString() );
    checkInputSize( m_gasOilRelPermMaxValue, 2, viewKeyStruct::gasOilRelPermMaxValueString() );
  }

  m_volFracScale = 1.0;
  for( localIndex ip = 0; ip < numPhases; ++ip )
  {
    GEOSX_THROW_IF( m_phaseMinVolumeFraction[ip] < 0.0 || m_phaseMinVolumeFraction[ip] > 1.0,
                    "VanGenuchtenBakerRelativePermeability: invalid phase min volume fraction value: " << m_phaseMinVolumeFraction[ip],
                    InputError );
    m_volFracScale -= m_phaseMinVolumeFraction[ip];
  }
  GEOSX_THROW_IF( m_volFracScale < 0.0, "VanGenuchtenBakerRelativePermeability: sum of min volume fractions exceeds 1.0", InputError );


  for( localIndex ip = 0; ip < 2; ++ip )
  {
    if( m_phaseOrder[PhaseType::WATER] >= 0 )
    {
      GEOSX_THROW_IF( m_waterOilRelPermExponentInv[ip] < 0.0,
                      "VanGenuchtenBakerRelativePermeability: invalid water-oil exponent inverse value: " << m_waterOilRelPermExponentInv[ip],
                      InputError );
      GEOSX_THROW_IF( m_waterOilRelPermMaxValue[ip] < 0.0 || m_waterOilRelPermMaxValue[ip] > 1.0,
                      "VanGenuchtenBakerRelativePermeability: invalid maximum value: " << m_waterOilRelPermMaxValue[ip],
                      InputError );
    }

    if( m_phaseOrder[PhaseType::GAS] >= 0 )
    {
      GEOSX_THROW_IF( m_gasOilRelPermExponentInv[ip] < 0.0,
                      "VanGenuchtenBakerRelativePermeability: invalid gas-oil exponent inverse value: " << m_gasOilRelPermExponentInv[ip],
                      InputError );
      GEOSX_THROW_IF( m_gasOilRelPermMaxValue[ip] < 0.0 || m_gasOilRelPermMaxValue[ip] > 1.0,
                      "VanGenuchtenBakerRelativePermeability: invalid maximum value: " << m_gasOilRelPermMaxValue[ip],
                      InputError );
    }
  }

  if( m_phaseOrder[PhaseType::WATER] >= 0 && m_phaseOrder[PhaseType::GAS] >= 0 )
  {
    real64 const mean = 0.5 * ( m_gasOilRelPermMaxValue[GasOilPairPhaseType::OIL]
                                + m_waterOilRelPermMaxValue[WaterOilPairPhaseType::OIL] );
    m_gasOilRelPermMaxValue[GasOilPairPhaseType::OIL]     = mean;
    m_waterOilRelPermMaxValue[WaterOilPairPhaseType::OIL] = mean;
  }

}

VanGenuchtenBakerRelativePermeability::KernelWrapper VanGenuchtenBakerRelativePermeability::createKernelWrapper()
{
  return KernelWrapper( m_phaseMinVolumeFraction,
                        m_waterOilRelPermExponentInv,
                        m_waterOilRelPermMaxValue,
                        m_gasOilRelPermExponentInv,
                        m_gasOilRelPermMaxValue,
                        m_volFracScale,
                        m_phaseTypes,
                        m_phaseOrder,
                        m_phaseRelPerm,
                        m_dPhaseRelPerm_dPhaseVolFrac );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, VanGenuchtenBakerRelativePermeability, string const &, Group * const )
} // namespace constitutive

} // namespace geosx
