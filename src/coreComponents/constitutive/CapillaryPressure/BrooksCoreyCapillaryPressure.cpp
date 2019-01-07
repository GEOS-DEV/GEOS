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
  * @file BrooksCoreyCapillaryPressure.cpp
  */

#include "BrooksCoreyCapillaryPressure.hpp"

#include <cmath>

namespace geosx
{

using namespace dataRepository;
using namespace cxx_utilities;

namespace constitutive
{


BrooksCoreyCapillaryPressure::BrooksCoreyCapillaryPressure( std::string const & name,
                                                            ManagedGroup * const parent )
  : CapillaryPressureBase( name, parent )
{
  RegisterViewWrapper( viewKeyStruct::phaseMinVolumeFractionString, &m_phaseMinVolumeFraction, false )->
    setApplyDefaultValue(0.0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Minimum volume fraction value for each phase");

  RegisterViewWrapper( viewKeyStruct::phaseCapPressureExponentInvString,   &m_phaseCapPressureExponentInv,   false )->
    setApplyDefaultValue(2.0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Inverse of capillary power law exponent for each phase");

  RegisterViewWrapper( viewKeyStruct::phaseEntryPressureString,   &m_phaseEntryPressure,   false )->
    setApplyDefaultValue(1.0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Entry pressure value for each phase");

  RegisterViewWrapper( viewKeyStruct::capPressureEpsilonString,   &m_capPressureEpsilon,   false )->
    setApplyDefaultValue(1e-7)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Wetting-phase saturation at which the max cap. pressure is attained; used to avoid infinite cap. pressure values for saturations close to zero");

}

BrooksCoreyCapillaryPressure::~BrooksCoreyCapillaryPressure()
{

}

std::unique_ptr<ConstitutiveBase>
BrooksCoreyCapillaryPressure::DeliverClone(string const & name,
					   ManagedGroup * const parent) const
{
  std::unique_ptr< BrooksCoreyCapillaryPressure > clone = std::make_unique<BrooksCoreyCapillaryPressure>( name, parent );

  clone->m_phaseNames = this->m_phaseNames;
  clone->m_phaseTypes = this->m_phaseTypes;
  clone->m_phaseOrder = this->m_phaseOrder;
  
  clone->m_phaseMinVolumeFraction      = this->m_phaseMinVolumeFraction;
  clone->m_phaseCapPressureExponentInv = this->m_phaseCapPressureExponentInv;
  clone->m_phaseEntryPressure          = this->m_phaseEntryPressure;

  clone->m_capPressureEpsilon = this->m_capPressureEpsilon;
  clone->m_satScale           = this->m_satScale;

  std::unique_ptr<ConstitutiveBase> rval = std::move( clone );
  return rval;
}


void BrooksCoreyCapillaryPressure::ProcessInputFile_PostProcess()
{
  CapillaryPressureBase::ProcessInputFile_PostProcess();

  localIndex const NP = numFluidPhases();

#define COREY_CHECK_INPUT_LENGTH( data, expected, attr ) \
  if (integer_conversion<localIndex>((data).size()) != integer_conversion<localIndex>(expected)) \
  { \
    GEOS_ERROR( "BrooksCoreyCapillaryPressure: invalid number of entries in " \
                << (attr) << " attribute (" \
                << (data).size() << "given, " \
                << (expected) << " expected)"); \
  }

  COREY_CHECK_INPUT_LENGTH( m_phaseMinVolumeFraction, NP, viewKeyStruct::phaseMinVolumeFractionString );
  COREY_CHECK_INPUT_LENGTH( m_phaseCapPressureExponentInv, NP, viewKeyStruct::phaseCapPressureExponentInvString );
  COREY_CHECK_INPUT_LENGTH( m_phaseEntryPressure, NP, viewKeyStruct::phaseEntryPressureString );

#undef COREY_CHECK_INPUT_LENGTH

  m_satScale = 1.0;
  for (localIndex ip = 0; ip < NP; ++ip)
  {
    GEOS_ERROR_IF( m_phaseMinVolumeFraction[ip] < 0.0 || m_phaseMinVolumeFraction[ip] > 1.0,
                   "BrooksCoreyCapillaryPressure: invalid min volume fraction value: " << m_phaseMinVolumeFraction[ip] );
    m_satScale -= m_phaseMinVolumeFraction[ip];

    GEOS_ERROR_IF( m_phaseCapPressureExponentInv[ip] < 1.0,
                   "BrooksCoreyCapillaryPressure: invalid exponent inverse value: " << m_phaseCapPressureExponentInv[ip] );

    GEOS_ERROR_IF( m_phaseEntryPressure[ip] < 0.0,
                   "BrooksCoreyCapillaryPressure: invalid entry pressure: " << m_phaseEntryPressure[ip] );

    GEOS_ERROR_IF( m_capPressureEpsilon < 0.0 || m_capPressureEpsilon > 0.2,
                   "BrooksCoreyCapillaryPressure: invalid epsilon: " << m_capPressureEpsilon );

  }

  GEOS_ERROR_IF( m_satScale < 0.0, "BrooksCoreyCapillaryPressure: sum of min volume fractions exceeds 1.0" );
}


void BrooksCoreyCapillaryPressure::BatchUpdate( arrayView2d<real64 const> const & phaseVolumeFraction )
{

  arrayView1d<real64 const> const & phaseMinVolumeFraction      = m_phaseMinVolumeFraction;
  arrayView1d<real64 const> const & phaseCapPressureExponentInv = m_phaseCapPressureExponentInv;
  arrayView1d<real64 const> const & phaseEntryPressure          = m_phaseEntryPressure;
  real64 const & capPressureEpsilon = m_capPressureEpsilon;

  CapillaryPressureBase::BatchUpdateKernel<BrooksCoreyCapillaryPressure>( phaseVolumeFraction,
                                                                          phaseMinVolumeFraction,
                                                                          phaseCapPressureExponentInv,
                                                                          phaseEntryPressure,
									  capPressureEpsilon,
                                                                          m_satScale );
}


void BrooksCoreyCapillaryPressure::PointUpdate( arraySlice1d<real64 const> const & phaseVolFraction,
                                                localIndex const k,
                                                localIndex const q )
{
  arraySlice1d<real64> const capPressure           = m_phaseCapPressure[k][q];
  arraySlice2d<real64> const dCapPressure_dVolFrac = m_dPhaseCapPressure_dPhaseVolFrac[k][q];

  localIndex const NP = numFluidPhases();

  Compute( NP,
	   m_phaseTypes,
           phaseVolFraction,
           capPressure,
           dCapPressure_dVolFrac,
           m_phaseMinVolumeFraction,
           m_phaseCapPressureExponentInv,
           m_phaseEntryPressure,
	   m_capPressureEpsilon,
           m_satScale );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, BrooksCoreyCapillaryPressure, std::string const &, ManagedGroup * const )
} // namespace constitutive

} // namespace geosx
