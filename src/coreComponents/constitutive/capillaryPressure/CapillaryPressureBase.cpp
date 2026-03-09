/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file CapillaryPressureBase.cpp
 */

#include "CapillaryPressureBase.hpp"
#include "CapillaryPressureFields.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

CapillaryPressureBase::CapillaryPressureBase( string const & name,
                                              Group * const parent )
  : ConstitutiveBase( name, parent )
{
  registerWrapper( viewKeyStruct::phaseNamesString(), &m_phaseNames ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRefArray ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "List of fluid phases" );

  registerWrapper( viewKeyStruct::phaseTypesString(), &m_phaseTypes ).
    setSizedFromParent( 0 );

  registerWrapper( viewKeyStruct::phaseOrderString(), &m_phaseOrder ).
    setSizedFromParent( 0 );

  registerWrapper( viewKeyStruct::phaseMinVolumeFractionString(), &m_phaseMinVolumeFraction ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setSizedFromParent( 0 ).
    setDescription( "Minimum volume fraction value for each phase" );

    registerField< fields::cappres::phaseCapPressure >( &m_phaseCapPressure );
    registerField< fields::cappres::dPhaseCapPressure_dPhaseVolFraction >( &m_dPhaseCapPressure_dPhaseVolFrac );
    registerField< fields::cappres::phaseTrappedVolFraction >( &m_phaseTrappedVolFrac );
}

void CapillaryPressureBase::postInputInitialization()
{
  ConstitutiveBase::postInputInitialization();

  integer const numPhases = numFluidPhases();
  GEOS_THROW_IF_LT_MSG( numPhases, 2,
                        "invalid number of phases",
                        InputError, getDataContext() );
  GEOS_THROW_IF_GT_MSG( numPhases, MAX_NUM_PHASES,
                        "invalid number of phases",
                        InputError, getDataContext() );

  m_phaseTypes.resize( numPhases );
  m_phaseOrder.resizeDefault( MAX_NUM_PHASES, -1 );

  auto const toPhaseType = [&]( string const & lookup )
  {
    static unordered_map< string, integer > const phaseDict =
    {
      { "gas", PhaseType::GAS },
      { "oil", PhaseType::OIL },
      { "water", PhaseType::WATER }
    };
    return findOption( phaseDict, lookup, viewKeyStruct::phaseNamesString(), getFullName() );
  };

  for( integer ip = 0; ip < numPhases; ++ip )
  {
    m_phaseTypes[ip] = toPhaseType( m_phaseNames[ip] );
    m_phaseOrder[m_phaseTypes[ip]] = ip;
  }

  // set labels on array wrappers for plottable fields
  setLabels();
}

void CapillaryPressureBase::allocateConstitutiveData( Group & parent, localIndex const numPts )
{
  integer const NP = numFluidPhases();
  //phase trapped for stats
  m_phaseTrappedVolFrac.resize( 0, numPts, NP );
  m_phaseTrappedVolFrac.zero();
  m_phaseCapPressure.resize( 0, numPts, NP );
  m_dPhaseCapPressure_dPhaseVolFrac.resize( 0, numPts, NP, NP );
  
  ConstitutiveBase::allocateConstitutiveData( parent, numPts );
}

void CapillaryPressureBase::setLabels()
{
  getField< fields::cappres::phaseTrappedVolFraction >().
    setDimLabels( 2, m_phaseNames );
  getField< fields::cappres::phaseCapPressure >().
    setDimLabels( 2, m_phaseNames );
}

void CapillaryPressureBase::resizeFields( localIndex const size, localIndex const numPts )
{
  integer const NP = numFluidPhases();
  m_phaseTrappedVolFrac.resize( size, numPts, NP );
  m_phaseCapPressure.resize( size, numPts, NP );
  m_dPhaseCapPressure_dPhaseVolFrac.resize( size, numPts, NP, NP );
}

} // namespace constitutive

} // namespace geos
