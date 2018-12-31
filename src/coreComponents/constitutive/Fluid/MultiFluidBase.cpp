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
  * @file MultiFluidBase.cpp
  */

#include "MultiFluidBase.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace cxx_utilities;

namespace constitutive
{

MultiFluidBase::MultiFluidBase( std::string const & name, ManagedGroup * const parent )
  : ConstitutiveBase( name, parent ),
    m_useMass( false )
{
  RegisterViewWrapper( viewKeyStruct::componentNamesString, &m_componentNames, false )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("List of component names");

  RegisterViewWrapper( viewKeyStruct::componentMolarWeightString, &m_componentMolarWeight, false )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Component molar weights");


  RegisterViewWrapper( viewKeyStruct::phaseNamesString, &m_phaseNames, false )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("List of fluid phases");

  RegisterViewWrapper( viewKeyStruct::phaseFractionString, &m_phaseFraction, false )->setPlotLevel(PlotLevel::LEVEL_0);
  RegisterViewWrapper( viewKeyStruct::dPhaseFraction_dPressureString, &m_dPhaseFraction_dPressure, false );
  RegisterViewWrapper( viewKeyStruct::dPhaseFraction_dTemperatureString, &m_dPhaseFraction_dTemperature, false );
  RegisterViewWrapper( viewKeyStruct::dPhaseFraction_dGlobalCompFractionString, &m_dPhaseFraction_dGlobalCompFraction, false );

  RegisterViewWrapper( viewKeyStruct::phaseDensityString, &m_phaseDensity, false )->setPlotLevel(PlotLevel::LEVEL_0);
  RegisterViewWrapper( viewKeyStruct::dPhaseDensity_dPressureString, &m_dPhaseDensity_dPressure, false );
  RegisterViewWrapper( viewKeyStruct::dPhaseDensity_dTemperatureString, &m_dPhaseDensity_dTemperature, false );
  RegisterViewWrapper( viewKeyStruct::dPhaseDensity_dGlobalCompFractionString, &m_dPhaseDensity_dGlobalCompFraction, false );

  RegisterViewWrapper( viewKeyStruct::phaseViscosityString, &m_phaseViscosity, false )->setPlotLevel(PlotLevel::LEVEL_0);
  RegisterViewWrapper( viewKeyStruct::dPhaseViscosity_dPressureString, &m_dPhaseViscosity_dPressure, false );
  RegisterViewWrapper( viewKeyStruct::dPhaseViscosity_dTemperatureString, &m_dPhaseViscosity_dTemperature, false );
  RegisterViewWrapper( viewKeyStruct::dPhaseViscosity_dGlobalCompFractionString, &m_dPhaseViscosity_dGlobalCompFraction, false );

  RegisterViewWrapper( viewKeyStruct::phaseCompFractionString, &m_phaseCompFraction, false )->setPlotLevel(PlotLevel::LEVEL_0);
  RegisterViewWrapper( viewKeyStruct::dPhaseCompFraction_dPressureString, &m_dPhaseCompFraction_dPressure, false );
  RegisterViewWrapper( viewKeyStruct::dPhaseCompFraction_dTemperatureString, &m_dPhaseCompFraction_dTemperature, false );
  RegisterViewWrapper( viewKeyStruct::dPhaseCompFraction_dGlobalCompFractionString, &m_dPhaseCompFraction_dGlobalCompFraction, false );

  RegisterViewWrapper( viewKeyStruct::totalDensityString, &m_totalDensity, false )->setPlotLevel( PlotLevel::LEVEL_0 );
  RegisterViewWrapper( viewKeyStruct::dTotalDensity_dPressureString, &m_dTotalDensity_dPressure, false );
  RegisterViewWrapper( viewKeyStruct::dTotalDensity_dTemperatureString, &m_dTotalDensity_dTemperature, false );
  RegisterViewWrapper( viewKeyStruct::dTotalDensity_dGlobalCompFractionString, &m_dTotalDensity_dGlobalCompFraction, false );
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

void MultiFluidBase::AllocateConstitutiveData( dataRepository::ManagedGroup * const parent,
                                               localIndex const numConstitutivePointsPerParentIndex )
{
  ConstitutiveBase::AllocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );
  ResizeFields( parent->size(), numConstitutivePointsPerParentIndex );
}

MultiFluidBase::~MultiFluidBase()
{

}

void MultiFluidBase::ProcessInputFile_PostProcess()
{
  ConstitutiveBase::ProcessInputFile_PostProcess();

  localIndex const NC = numFluidComponents();
  localIndex const NP = numFluidPhases();

  GEOS_ERROR_IF( NC == 0, "MultiFluidBase: No fluid components specified" );

  GEOS_ERROR_IF( NC > MAX_NUM_COMPONENTS,
                 "MultiFluidBase: Number of fluid components exceeds the maximum of " << MAX_NUM_COMPONENTS );

  GEOS_ERROR_IF( NP == 0, "MultiFluidBase: No fluid phases specified" );

  GEOS_ERROR_IF( NP > MAX_NUM_PHASES,
                 "MultiFluidBase: Number of fluid phases exceeds the maximum of " << MAX_NUM_PHASES );

#define MULTIFLUID_CHECK_INPUT_LENGTH( data, expected, attr ) \
  if (integer_conversion<localIndex>((data).size()) != integer_conversion<localIndex>(expected)) \
  { \
    GEOS_ERROR( "MultiFluidBase: invalid number of entries in " \
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
  GEOS_ERROR_IF( ic >= numFluidComponents(), "Index " << ic << " exceeds number of fluid components" );
  return m_componentNames[ic];
}

localIndex MultiFluidBase::numFluidPhases() const
{
  return integer_conversion<localIndex>(m_phaseNames.size());
}

string const & MultiFluidBase::phaseName(localIndex ip) const
{
  GEOS_ERROR_IF( ip >= numFluidPhases(), "Index " << ip << " exceeds number of fluid phases" );
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
