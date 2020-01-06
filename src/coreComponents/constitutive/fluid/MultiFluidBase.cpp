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
  * @file MultiFluidBase.cpp
  */

#include "MultiFluidBase.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace cxx_utilities;

namespace constitutive
{

MultiFluidBase::MultiFluidBase( std::string const & name, Group * const parent )
  : ConstitutiveBase( name, parent ),
    m_useMass( false )
{
  // We make base inputs optional here, since derived classes may want to predefine/hardcode
  // components/phases. Models that do need these inputs should change input flags accordingly.

  registerWrapper( viewKeyStruct::componentNamesString, &m_componentNames, false )->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("List of component names");

  registerWrapper( viewKeyStruct::componentMolarWeightString, &m_componentMolarWeight, false )->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Component molar weights");

  registerWrapper( viewKeyStruct::phaseNamesString, &m_phaseNames, false )->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("List of fluid phases");

  registerWrapper( viewKeyStruct::phaseFractionString, &m_phaseFraction, false )->setPlotLevel(PlotLevel::LEVEL_0);
  registerWrapper( viewKeyStruct::dPhaseFraction_dPressureString, &m_dPhaseFraction_dPressure, false );
  registerWrapper( viewKeyStruct::dPhaseFraction_dTemperatureString, &m_dPhaseFraction_dTemperature, false );
  registerWrapper( viewKeyStruct::dPhaseFraction_dGlobalCompFractionString, &m_dPhaseFraction_dGlobalCompFraction, false );

  registerWrapper( viewKeyStruct::phaseDensityString, &m_phaseDensity, false )->setPlotLevel(PlotLevel::LEVEL_0);
  registerWrapper( viewKeyStruct::dPhaseDensity_dPressureString, &m_dPhaseDensity_dPressure, false );
  registerWrapper( viewKeyStruct::dPhaseDensity_dTemperatureString, &m_dPhaseDensity_dTemperature, false );
  registerWrapper( viewKeyStruct::dPhaseDensity_dGlobalCompFractionString, &m_dPhaseDensity_dGlobalCompFraction, false );

  registerWrapper( viewKeyStruct::phaseViscosityString, &m_phaseViscosity, false )->setPlotLevel(PlotLevel::LEVEL_0);
  registerWrapper( viewKeyStruct::dPhaseViscosity_dPressureString, &m_dPhaseViscosity_dPressure, false );
  registerWrapper( viewKeyStruct::dPhaseViscosity_dTemperatureString, &m_dPhaseViscosity_dTemperature, false );
  registerWrapper( viewKeyStruct::dPhaseViscosity_dGlobalCompFractionString, &m_dPhaseViscosity_dGlobalCompFraction, false );

  registerWrapper( viewKeyStruct::phaseCompFractionString, &m_phaseCompFraction, false )->setPlotLevel(PlotLevel::LEVEL_0);
  registerWrapper( viewKeyStruct::dPhaseCompFraction_dPressureString, &m_dPhaseCompFraction_dPressure, false );
  registerWrapper( viewKeyStruct::dPhaseCompFraction_dTemperatureString, &m_dPhaseCompFraction_dTemperature, false );
  registerWrapper( viewKeyStruct::dPhaseCompFraction_dGlobalCompFractionString, &m_dPhaseCompFraction_dGlobalCompFraction, false );

  registerWrapper( viewKeyStruct::totalDensityString, &m_totalDensity, false )->setPlotLevel( PlotLevel::LEVEL_0 );
  registerWrapper( viewKeyStruct::dTotalDensity_dPressureString, &m_dTotalDensity_dPressure, false );
  registerWrapper( viewKeyStruct::dTotalDensity_dTemperatureString, &m_dTotalDensity_dTemperature, false );
  registerWrapper( viewKeyStruct::dTotalDensity_dGlobalCompFractionString, &m_dTotalDensity_dGlobalCompFraction, false );
}

void MultiFluidBase::ResizeFields( localIndex const size, localIndex const numPts )
{
  localIndex const NP = numFluidPhases();
  localIndex const NC = numFluidComponents();

  m_phaseFraction.resize( size, numPts, NP );
  m_dPhaseFraction_dPressure.resize( size, numPts, NP );
  m_dPhaseFraction_dTemperature.resize( size, numPts, NP );
  m_dPhaseFraction_dGlobalCompFraction.resize( size, numPts, NP, NC );

  m_phaseDensity.resize( size, numPts, NP );
  m_dPhaseDensity_dPressure.resize( size, numPts, NP );
  m_dPhaseDensity_dTemperature.resize( size, numPts, NP );
  m_dPhaseDensity_dGlobalCompFraction.resize( size, numPts, NP, NC );

  m_phaseViscosity.resize( size, numPts, NP );
  m_dPhaseViscosity_dPressure.resize( size, numPts, NP );
  m_dPhaseViscosity_dTemperature.resize( size, numPts, NP );
  m_dPhaseViscosity_dGlobalCompFraction.resize( size, numPts, NP, NC );

  m_phaseCompFraction.resize( size, numPts, NP, NC );
  m_dPhaseCompFraction_dPressure.resize( size, numPts, NP, NC );
  m_dPhaseCompFraction_dTemperature.resize( size, numPts, NP, NC );
  m_dPhaseCompFraction_dGlobalCompFraction.resize( size, numPts, NP, NC, NC );

  m_totalDensity.resize( size, numPts );
  m_dTotalDensity_dPressure.resize( size, numPts );
  m_dTotalDensity_dTemperature.resize( size, numPts );
  m_dTotalDensity_dGlobalCompFraction.resize( size, numPts, NC );
}

void MultiFluidBase::AllocateConstitutiveData( dataRepository::Group * const parent,
                                               localIndex const numConstitutivePointsPerParentIndex )
{
  ConstitutiveBase::AllocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );
  ResizeFields( parent->size(), numConstitutivePointsPerParentIndex );
}

MultiFluidBase::~MultiFluidBase()
{

}

void MultiFluidBase::PostProcessInput()
{
  ConstitutiveBase::PostProcessInput();

  localIndex const NC = numFluidComponents();
  localIndex const NP = numFluidPhases();

  GEOSX_ERROR_IF( NC == 0, "MultiFluidBase: No fluid components specified" );

  GEOSX_ERROR_IF( NC > MAX_NUM_COMPONENTS,
                 "MultiFluidBase: Number of fluid components exceeds the maximum of " << MAX_NUM_COMPONENTS );

  GEOSX_ERROR_IF( NP == 0, "MultiFluidBase: No fluid phases specified" );

  GEOSX_ERROR_IF( NP > MAX_NUM_PHASES,
                 "MultiFluidBase: Number of fluid phases exceeds the maximum of " << MAX_NUM_PHASES );

#define MULTIFLUID_CHECK_INPUT_LENGTH( data, expected, attr ) \
  if (integer_conversion<localIndex>((data).size()) != integer_conversion<localIndex>(expected)) \
  { \
    GEOSX_ERROR( "MultiFluidBase: invalid number of entries in " \
                << (attr) << " attribute (" \
                << (data).size() << "given, " \
                << (expected) << " expected)" ); \
  }

  MULTIFLUID_CHECK_INPUT_LENGTH( m_componentMolarWeight, NC,
                                 viewKeyStruct::componentMolarWeightString )

  // call to correctly set member array tertiary sizes on the 'main' material object
  ResizeFields( 0, 0 );
}

localIndex MultiFluidBase::numFluidComponents() const
{
  return integer_conversion<localIndex>(m_componentNames.size());
}

string const & MultiFluidBase::componentName(localIndex ic) const
{
  GEOSX_ERROR_IF( ic >= numFluidComponents(), "Index " << ic << " exceeds number of fluid components" );
  return m_componentNames[ic];
}

localIndex MultiFluidBase::numFluidPhases() const
{
  return integer_conversion<localIndex>(m_phaseNames.size());
}

string const & MultiFluidBase::phaseName(localIndex ip) const
{
  GEOSX_ERROR_IF( ip >= numFluidPhases(), "Index " << ip << " exceeds number of fluid phases" );
  return m_phaseNames[ip];
}

bool MultiFluidBase::getMassFlag() const
{
  return m_useMass;
}

void MultiFluidBase::setMassFlag(bool flag)
{
  m_useMass = flag;
}

} //namespace constitutive

} //namespace geosx
