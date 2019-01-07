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
  * @file BrooksCoreyRelativePermeability.cpp
  */

#include "BrooksCoreyRelativePermeability.hpp"

#include <cmath>

namespace geosx
{

using namespace dataRepository;
using namespace cxx_utilities;

namespace constitutive
{


BrooksCoreyRelativePermeability::BrooksCoreyRelativePermeability( std::string const & name,
                                                                  ManagedGroup * const parent )
  : RelativePermeabilityBase( name, parent )
{
  RegisterViewWrapper( viewKeyStruct::phaseMinVolumeFractionString, &m_phaseMinVolumeFraction, false )->
    setApplyDefaultValue(0.0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Minimum volume fraction value for each phase");

  RegisterViewWrapper( viewKeyStruct::phaseRelPermExponentString,   &m_phaseRelPermExponent,   false )->
    setApplyDefaultValue(1.0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("MinimumRel perm power law exponent for each phase");


  RegisterViewWrapper( viewKeyStruct::phaseRelPermMaxValueString,   &m_phaseRelPermMaxValue,   false )->
    setApplyDefaultValue(0.0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Maximum rel perm value for each phase");

}

BrooksCoreyRelativePermeability::~BrooksCoreyRelativePermeability()
{

}

std::unique_ptr<ConstitutiveBase>
BrooksCoreyRelativePermeability::DeliverClone(string const & name, ManagedGroup * const parent) const
{
  std::unique_ptr< BrooksCoreyRelativePermeability > clone = std::make_unique<BrooksCoreyRelativePermeability>( name, parent );

  clone->m_phaseNames = this->m_phaseNames;
  clone->m_phaseTypes = this->m_phaseTypes;
  clone->m_phaseOrder = this->m_phaseOrder;

  clone->m_phaseMinVolumeFraction = this->m_phaseMinVolumeFraction;
  clone->m_phaseRelPermExponent   = this->m_phaseRelPermExponent;
  clone->m_phaseRelPermMaxValue   = this->m_phaseRelPermMaxValue;

  clone->m_satScale = this->m_satScale;

  std::unique_ptr<ConstitutiveBase> rval = std::move( clone );
  return rval;
}


void BrooksCoreyRelativePermeability::PostProcessInput()
{
  RelativePermeabilityBase::PostProcessInput();

  localIndex const NP = numFluidPhases();

#define COREY_CHECK_INPUT_LENGTH( data, expected, attr ) \
  if (integer_conversion<localIndex>((data).size()) != integer_conversion<localIndex>(expected)) \
  { \
    GEOS_ERROR( "BrooksCoreyRelativePermeability: invalid number of entries in " \
                << (attr) << " attribute (" \
                << (data).size() << "given, " \
                << (expected) << " expected)"); \
  }

  COREY_CHECK_INPUT_LENGTH( m_phaseMinVolumeFraction, NP, viewKeyStruct::phaseMinVolumeFractionString );
  COREY_CHECK_INPUT_LENGTH( m_phaseRelPermExponent, NP, viewKeyStruct::phaseRelPermExponentString );
  COREY_CHECK_INPUT_LENGTH( m_phaseRelPermMaxValue, NP, viewKeyStruct::phaseRelPermMaxValueString );

#undef COREY_CHECK_INPUT_LENGTH

  m_satScale = 1.0;
  for (localIndex ip = 0; ip < NP; ++ip)
  {
    GEOS_ERROR_IF( m_phaseMinVolumeFraction[ip] < 0.0 || m_phaseMinVolumeFraction[ip] > 1.0,
                   "BrooksCoreyRelativePermeability: invalid min volume fraction value: " << m_phaseMinVolumeFraction[ip] );
    m_satScale -= m_phaseMinVolumeFraction[ip];

    GEOS_ERROR_IF( m_phaseRelPermExponent[ip] < 0.0,
                   "BrooksCoreyRelativePermeability: invalid exponent value: " << m_phaseRelPermExponent[ip] );

    GEOS_ERROR_IF( m_phaseRelPermMaxValue[ip] < 0.0 || m_phaseRelPermMaxValue[ip] > 1.0,
                   "BrooksCoreyRelativePermeability: invalid maximum value: " << m_phaseRelPermMaxValue[ip] );
  }

  GEOS_ERROR_IF( m_satScale < 0.0, "BrooksCoreyRelativePermeability: sum of min volume fractions exceeds 1.0" );
}


void BrooksCoreyRelativePermeability::BatchUpdate( arrayView2d<real64 const> const & phaseVolumeFraction )
{

  arrayView1d<real64 const> const &  phaseMinVolumeFraction = m_phaseMinVolumeFraction;
  arrayView1d<real64 const> const & phaseRelPermExponent = m_phaseRelPermExponent;
  arrayView1d<real64 const> const & phaseRelPermMaxValue = m_phaseRelPermMaxValue;


  RelativePermeabilityBase::BatchUpdateKernel<BrooksCoreyRelativePermeability>( phaseVolumeFraction,
                                                                                phaseMinVolumeFraction,
                                                                                phaseRelPermExponent,
                                                                                phaseRelPermMaxValue,
                                                                                m_satScale );
}


void BrooksCoreyRelativePermeability::StateUpdatePointRelPerm( arraySlice1d<real64 const> const & phaseVolFraction,
                                                               localIndex const k,
                                                               localIndex const q )
{
  arraySlice1d<real64> const relPerm           = m_phaseRelPerm[k][q];
  arraySlice2d<real64> const dRelPerm_dVolFrac = m_dPhaseRelPerm_dPhaseVolFrac[k][q];

  localIndex const NP = numFluidPhases();

  StateUpdatePointRelPerm( NP,
                           phaseVolFraction,
                           relPerm,
                           dRelPerm_dVolFrac,
                           m_phaseMinVolumeFraction,
                           m_phaseRelPermExponent,
                           m_phaseRelPermMaxValue,
                           m_satScale );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, BrooksCoreyRelativePermeability, std::string const &, ManagedGroup * const )
} // namespace constitutive

} // namespace geosx
