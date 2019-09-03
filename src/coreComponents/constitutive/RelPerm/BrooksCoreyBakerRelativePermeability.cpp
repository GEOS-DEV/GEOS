/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
  * @file BrooksCoreyBakerRelativePermeability.cpp
  */

#include "BrooksCoreyBakerRelativePermeability.hpp"

#include <cmath>

namespace geosx
{

using namespace dataRepository;
using namespace cxx_utilities;

namespace constitutive
{


BrooksCoreyBakerRelativePermeability::BrooksCoreyBakerRelativePermeability( std::string const & name,
                                                                            Group * const parent )
  : RelativePermeabilityBase( name, parent )
{
  registerWrapper( viewKeyStruct::phaseMinVolumeFractionString, &m_phaseMinVolumeFraction, false )->
    setApplyDefaultValue(0.0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Minimum volume fraction value for each phase");


  registerWrapper( viewKeyStruct::waterOilRelPermExponentString,   &m_waterOilRelPermExponent,   false )->
    setApplyDefaultValue(1.0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Rel perm power law exponent for the pair (water phase, oil phase) at residual gas saturation");

  registerWrapper( viewKeyStruct::waterOilRelPermMaxValueString,   &m_waterOilRelPermMaxValue,   false )->
    setApplyDefaultValue(0.0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Maximum rel perm value for the pair (water phase, oil phase) at residual gas saturation");


  registerWrapper( viewKeyStruct::gasOilRelPermExponentString,   &m_gasOilRelPermExponent,   false )->
    setApplyDefaultValue(1.0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Rel perm power law exponent for the pair (gas phase, oil phase) at residual water saturation");

  registerWrapper( viewKeyStruct::gasOilRelPermMaxValueString,   &m_gasOilRelPermMaxValue,   false )->
    setApplyDefaultValue(0.0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Maximum rel perm value for the pair (gas phase, oil phase) at residual water saturation");

}

BrooksCoreyBakerRelativePermeability::~BrooksCoreyBakerRelativePermeability()
{

}

void
BrooksCoreyBakerRelativePermeability::DeliverClone( string const & name, 
                                                    Group * const parent, 
                                                    std::unique_ptr<ConstitutiveBase> & clone ) const
{
  std::unique_ptr< BrooksCoreyBakerRelativePermeability > newModel = std::make_unique<BrooksCoreyBakerRelativePermeability>( name, parent );

  newModel->m_phaseNames = this->m_phaseNames;
  newModel->m_phaseTypes = this->m_phaseTypes;
  newModel->m_phaseOrder = this->m_phaseOrder;

  newModel->m_phaseMinVolumeFraction = this->m_phaseMinVolumeFraction;

  newModel->m_waterOilRelPermExponent = this->m_waterOilRelPermExponent;
  newModel->m_waterOilRelPermMaxValue = this->m_waterOilRelPermMaxValue;

  newModel->m_gasOilRelPermExponent   = this->m_gasOilRelPermExponent;
  newModel->m_gasOilRelPermMaxValue   = this->m_gasOilRelPermMaxValue;

  newModel->m_volFracScale = this->m_volFracScale;

  clone = std::move( newModel );
}


void BrooksCoreyBakerRelativePermeability::PostProcessInput()
{
  RelativePermeabilityBase::PostProcessInput();

  localIndex const NP = numFluidPhases();

  GEOS_ERROR_IF( m_phaseOrder[PhaseType::OIL] < 0 , "BrooksCoreyBakerRelativePermeability: reference oil phase has not been defined and must be included in model" );

#define COREY_CHECK_INPUT_LENGTH( data, expected, attr ) \
  if (integer_conversion<localIndex>((data).size()) != integer_conversion<localIndex>(expected)) \
  { \
    GEOS_ERROR( "BrooksCoreyBakerRelativePermeability: invalid number of entries in " \
                << (attr) << " attribute (" \
                << (data).size() << " given, " \
                << (expected) << " expected)"); \
  }

  COREY_CHECK_INPUT_LENGTH( m_phaseMinVolumeFraction,  NP,   viewKeyStruct::phaseMinVolumeFractionString )

  if (m_phaseOrder[PhaseType::WATER] >= 0)
    {
      COREY_CHECK_INPUT_LENGTH( m_waterOilRelPermExponent, 2, viewKeyStruct::waterOilRelPermExponentString )
      COREY_CHECK_INPUT_LENGTH( m_waterOilRelPermMaxValue, 2, viewKeyStruct::waterOilRelPermMaxValueString )
    }

  if (m_phaseOrder[PhaseType::GAS] >=0)
    {
      COREY_CHECK_INPUT_LENGTH( m_gasOilRelPermExponent,   2, viewKeyStruct::gasOilRelPermExponentString )
      COREY_CHECK_INPUT_LENGTH( m_gasOilRelPermMaxValue,   2, viewKeyStruct::gasOilRelPermMaxValueString )
    }

#undef COREY_CHECK_INPUT_LENGTH

  m_volFracScale = 1.0;
  for (localIndex ip = 0; ip < NP; ++ip)
  {
    GEOS_ERROR_IF( m_phaseMinVolumeFraction[ip] < 0.0 || m_phaseMinVolumeFraction[ip] > 1.0,
                   "BrooksCoreyBakerRelativePermeability: invalid phase min volume fraction value: " << m_phaseMinVolumeFraction[ip] );
    m_volFracScale -= m_phaseMinVolumeFraction[ip];
  }
  GEOS_ERROR_IF( m_volFracScale < 0.0, "BrooksCoreyBakerRelativePermeability: sum of min volume fractions exceeds 1.0" );


  for (localIndex ip = 0; ip < 2; ++ip)
  {
    if (m_phaseOrder[PhaseType::WATER] >= 0)
      {
        GEOS_ERROR_IF( m_waterOilRelPermExponent[ip] < 0.0,
                       "BrooksCoreyBakerRelativePermeability: invalid water-oil exponent value: " << m_waterOilRelPermExponent[ip] );
        GEOS_ERROR_IF( m_waterOilRelPermMaxValue[ip] < 0.0 || m_waterOilRelPermMaxValue[ip] > 1.0,
                       "BrooksCoreyBakerRelativePermeability: invalid maximum value: " << m_waterOilRelPermMaxValue[ip] );
      }

    if (m_phaseOrder[PhaseType::GAS] >= 0)
      {
        GEOS_ERROR_IF( m_gasOilRelPermExponent[ip] < 0.0,
                       "BrooksCoreyBakerRelativePermeability: invalid gas-oil exponent value: " << m_gasOilRelPermExponent[ip] );
        GEOS_ERROR_IF( m_gasOilRelPermMaxValue[ip] < 0.0 || m_gasOilRelPermMaxValue[ip] > 1.0,
                       "BrooksCoreyBakerRelativePermeability: invalid maximum value: " << m_gasOilRelPermMaxValue[ip] );
      }
  }

  if (m_phaseOrder[PhaseType::WATER] >= 0 && m_phaseOrder[PhaseType::GAS] >= 0)
  {
    real64 const mean = 0.5 * ( m_gasOilRelPermMaxValue[GasOilPairPhaseType::OIL]
                              + m_waterOilRelPermMaxValue[WaterOilPairPhaseType::OIL] );
    m_gasOilRelPermMaxValue[GasOilPairPhaseType::OIL]     = mean;
    m_waterOilRelPermMaxValue[WaterOilPairPhaseType::OIL] = mean;
  }
}


void BrooksCoreyBakerRelativePermeability::BatchUpdate( arrayView2d<real64 const> const & phaseVolumeFraction )
{

  arrayView1d<real64 const> const & phaseMinVolumeFraction = m_phaseMinVolumeFraction;

  arrayView1d<real64 const> const & waterOilRelPermExponent = m_waterOilRelPermExponent;
  arrayView1d<real64 const> const & waterOilRelPermMaxValue = m_waterOilRelPermMaxValue;

  arrayView1d<real64 const> const & gasOilRelPermExponent   = m_gasOilRelPermExponent;
  arrayView1d<real64 const> const & gasOilRelPermMaxValue   = m_gasOilRelPermMaxValue;

  RelativePermeabilityBase::BatchUpdateKernel<BrooksCoreyBakerRelativePermeability>( phaseVolumeFraction,
                                                                                     m_phaseOrder,
                                                                                     phaseMinVolumeFraction,
                                                                                     waterOilRelPermExponent,
                                                                                     waterOilRelPermMaxValue,
                                                                                     gasOilRelPermExponent,
                                                                                     gasOilRelPermMaxValue,
                                                                                     m_volFracScale );
}


void BrooksCoreyBakerRelativePermeability::PointUpdate( arraySlice1d<real64 const> const & phaseVolFraction,
                                                        localIndex const k,
                                                        localIndex const q )
{
  arraySlice1d<real64> const relPerm           = m_phaseRelPerm[k][q];
  arraySlice2d<real64> const dRelPerm_dVolFrac = m_dPhaseRelPerm_dPhaseVolFrac[k][q];

  localIndex const NP = numFluidPhases();

  Compute( NP,
           phaseVolFraction,
           relPerm,
           dRelPerm_dVolFrac,
           m_phaseOrder,
           m_phaseMinVolumeFraction,
           m_waterOilRelPermExponent,
           m_waterOilRelPermMaxValue,
           m_gasOilRelPermExponent,
           m_gasOilRelPermMaxValue,
           m_volFracScale );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, BrooksCoreyBakerRelativePermeability, std::string const &, Group * const )
} // namespace constitutive

} // namespace geosx
