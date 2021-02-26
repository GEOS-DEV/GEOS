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
 * @file BrooksCoreyBakerRelativePermeability.cpp
 */

#include "BrooksCoreyBakerRelativePermeability.hpp"

#include <cmath>

namespace geosx
{

using namespace dataRepository;

namespace constitutive
{


BrooksCoreyBakerRelativePermeability::BrooksCoreyBakerRelativePermeability( string const & name,
                                                                            Group * const parent )
  : RelativePermeabilityBase( name, parent )
{
  registerWrapper( viewKeyStruct::phaseMinVolumeFractionString(), &m_phaseMinVolumeFraction ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Minimum volume fraction value for each phase" );


  registerWrapper( viewKeyStruct::waterOilRelPermExponentString(), &m_waterOilRelPermExponent ).
    setApplyDefaultValue( 1.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Rel perm power law exponent for the pair (water phase, oil phase) at residual gas saturation\n"
                    "The expected format is \"{ waterExp, oilExp }\", in that order" );

  registerWrapper( viewKeyStruct::waterOilRelPermMaxValueString(), &m_waterOilRelPermMaxValue ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Maximum rel perm value for the pair (water phase, oil phase) at residual gas saturation\n"
                    "The expected format is \"{ waterMax, oilMax }\", in that order" );


  registerWrapper( viewKeyStruct::gasOilRelPermExponentString(), &m_gasOilRelPermExponent ).
    setApplyDefaultValue( 1.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Rel perm power law exponent for the pair (gas phase, oil phase) at residual water saturation\n"
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

BrooksCoreyBakerRelativePermeability::~BrooksCoreyBakerRelativePermeability()
{}

void BrooksCoreyBakerRelativePermeability::postProcessInput()
{
  RelativePermeabilityBase::postProcessInput();

  localIndex const numPhases = numFluidPhases();

  GEOSX_ERROR_IF( m_phaseOrder[PhaseType::OIL] < 0,
                  "BrooksCoreyBakerRelativePermeability: reference oil phase has not been defined and must be included in model" );

  #define COREY_CHECK_INPUT_LENGTH( data, expected, attr ) \
    if( LvArray::integerConversion< localIndex >((data).size()) != LvArray::integerConversion< localIndex >( expected )) \
    { \
      GEOSX_ERROR( "BrooksCoreyBakerRelativePermeability: invalid number of entries in " \
                   << (attr) << " attribute (" \
                   << (data).size() << " given, " \
                   << (expected) << " expected)" ); \
    }

  COREY_CHECK_INPUT_LENGTH( m_phaseMinVolumeFraction, numPhases, viewKeyStruct::phaseMinVolumeFractionString() )

  if( m_phaseOrder[PhaseType::WATER] >= 0 )
  {
    COREY_CHECK_INPUT_LENGTH( m_waterOilRelPermExponent, 2, viewKeyStruct::waterOilRelPermExponentString() )
    COREY_CHECK_INPUT_LENGTH( m_waterOilRelPermMaxValue, 2, viewKeyStruct::waterOilRelPermMaxValueString() )
  }

  if( m_phaseOrder[PhaseType::GAS] >=0 )
  {
    COREY_CHECK_INPUT_LENGTH( m_gasOilRelPermExponent, 2, viewKeyStruct::gasOilRelPermExponentString() )
    COREY_CHECK_INPUT_LENGTH( m_gasOilRelPermMaxValue, 2, viewKeyStruct::gasOilRelPermMaxValueString() )
  }

#undef COREY_CHECK_INPUT_LENGTH

  m_volFracScale = 1.0;
  for( localIndex ip = 0; ip < numPhases; ++ip )
  {
    GEOSX_ERROR_IF( m_phaseMinVolumeFraction[ip] < 0.0 || m_phaseMinVolumeFraction[ip] > 1.0,
                    "BrooksCoreyBakerRelativePermeability: invalid phase min volume fraction value: " << m_phaseMinVolumeFraction[ip] );
    m_volFracScale -= m_phaseMinVolumeFraction[ip];
  }
  GEOSX_ERROR_IF( m_volFracScale < 0.0, "BrooksCoreyBakerRelativePermeability: sum of min volume fractions exceeds 1.0" );


  for( localIndex ip = 0; ip < 2; ++ip )
  {
    if( m_phaseOrder[PhaseType::WATER] >= 0 )
    {
      GEOSX_ERROR_IF( m_waterOilRelPermExponent[ip] < 0.0,
                      "BrooksCoreyBakerRelativePermeability: invalid water-oil exponent value: " << m_waterOilRelPermExponent[ip] );
      GEOSX_ERROR_IF( m_waterOilRelPermMaxValue[ip] < 0.0 || m_waterOilRelPermMaxValue[ip] > 1.0,
                      "BrooksCoreyBakerRelativePermeability: invalid maximum value: " << m_waterOilRelPermMaxValue[ip] );
    }

    if( m_phaseOrder[PhaseType::GAS] >= 0 )
    {
      GEOSX_ERROR_IF( m_gasOilRelPermExponent[ip] < 0.0,
                      "BrooksCoreyBakerRelativePermeability: invalid gas-oil exponent value: " << m_gasOilRelPermExponent[ip] );
      GEOSX_ERROR_IF( m_gasOilRelPermMaxValue[ip] < 0.0 || m_gasOilRelPermMaxValue[ip] > 1.0,
                      "BrooksCoreyBakerRelativePermeability: invalid maximum value: " << m_gasOilRelPermMaxValue[ip] );
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

BrooksCoreyBakerRelativePermeability::KernelWrapper BrooksCoreyBakerRelativePermeability::createKernelWrapper()
{
  return KernelWrapper( m_phaseMinVolumeFraction,
                        m_waterOilRelPermExponent,
                        m_waterOilRelPermMaxValue,
                        m_gasOilRelPermExponent,
                        m_gasOilRelPermMaxValue,
                        m_volFracScale,
                        m_phaseTypes,
                        m_phaseOrder,
                        m_phaseRelPerm,
                        m_dPhaseRelPerm_dPhaseVolFrac );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, BrooksCoreyBakerRelativePermeability, string const &, Group * const )
} // namespace constitutive

} // namespace geosx
