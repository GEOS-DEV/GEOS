/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
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

namespace constitutive
{

MultiFluidBase::MultiFluidBase( string const & name, Group * const parent )
  : ConstitutiveBase( name, parent ),
  m_useMass( false )
{
  // We make base inputs optional here, since derived classes may want to predefine/hardcode
  // components/phases. Models that do need these inputs should change input flags accordingly.

  registerWrapper( viewKeyStruct::componentNamesString(), &m_componentNames ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "List of component names" );

  registerWrapper( viewKeyStruct::componentMolarWeightString(), &m_componentMolarWeight ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Component molar weights" );

  registerWrapper( viewKeyStruct::phaseNamesString(), &m_phaseNames ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "List of fluid phases" );

  registerWrapper( viewKeyStruct::phaseFractionString(), &m_phaseFraction.value ).
    setPlotLevel( PlotLevel::LEVEL_0 );
  registerWrapper( viewKeyStruct::dPhaseFraction_dPressureString(), &m_phaseFraction.dPres ).
    setRestartFlags( RestartFlags::NO_WRITE );
  registerWrapper( viewKeyStruct::dPhaseFraction_dTemperatureString(), &m_phaseFraction.dTemp ).
    setRestartFlags( RestartFlags::NO_WRITE );
  registerWrapper( viewKeyStruct::dPhaseFraction_dGlobalCompFractionString(), &m_phaseFraction.dComp ).
    setRestartFlags( RestartFlags::NO_WRITE );

  registerWrapper( viewKeyStruct::phaseDensityString(), &m_phaseDensity.value ).
    setPlotLevel( PlotLevel::LEVEL_0 );
  registerWrapper( viewKeyStruct::dPhaseDensity_dPressureString(), &m_phaseDensity.dPres ).
    setRestartFlags( RestartFlags::NO_WRITE );
  registerWrapper( viewKeyStruct::dPhaseDensity_dTemperatureString(), &m_phaseDensity.dTemp ).
    setRestartFlags( RestartFlags::NO_WRITE );
  registerWrapper( viewKeyStruct::dPhaseDensity_dGlobalCompFractionString(), &m_phaseDensity.dComp ).
    setRestartFlags( RestartFlags::NO_WRITE );

  registerWrapper( viewKeyStruct::phaseMassDensityString(), &m_phaseMassDensity.value ).
    setPlotLevel( PlotLevel::LEVEL_0 );
  registerWrapper( viewKeyStruct::dPhaseMassDensity_dPressureString(), &m_phaseMassDensity.dPres ).
    setRestartFlags( RestartFlags::NO_WRITE );
  registerWrapper( viewKeyStruct::dPhaseMassDensity_dTemperatureString(), &m_phaseMassDensity.dTemp ).
    setRestartFlags( RestartFlags::NO_WRITE );
  registerWrapper( viewKeyStruct::dPhaseMassDensity_dGlobalCompFractionString(), &m_phaseMassDensity.dComp ).
    setRestartFlags( RestartFlags::NO_WRITE );

  registerWrapper( viewKeyStruct::phaseViscosityString(), &m_phaseViscosity.value ).
    setPlotLevel( PlotLevel::LEVEL_0 );
  registerWrapper( viewKeyStruct::dPhaseViscosity_dPressureString(), &m_phaseViscosity.dPres ).
    setRestartFlags( RestartFlags::NO_WRITE );
  registerWrapper( viewKeyStruct::dPhaseViscosity_dTemperatureString(), &m_phaseViscosity.dTemp ).
    setRestartFlags( RestartFlags::NO_WRITE );
  registerWrapper( viewKeyStruct::dPhaseViscosity_dGlobalCompFractionString(), &m_phaseViscosity.dComp ).
    setRestartFlags( RestartFlags::NO_WRITE );

  registerWrapper( viewKeyStruct::phaseCompFractionString(), &m_phaseCompFraction.value ).
    setPlotLevel( PlotLevel::LEVEL_0 );
  registerWrapper( viewKeyStruct::dPhaseCompFraction_dPressureString(), &m_phaseCompFraction.dPres ).
    setRestartFlags( RestartFlags::NO_WRITE );
  registerWrapper( viewKeyStruct::dPhaseCompFraction_dTemperatureString(), &m_phaseCompFraction.dTemp ).
    setRestartFlags( RestartFlags::NO_WRITE );
  registerWrapper( viewKeyStruct::dPhaseCompFraction_dGlobalCompFractionString(), &m_phaseCompFraction.dComp ).
    setRestartFlags( RestartFlags::NO_WRITE );

  registerWrapper( viewKeyStruct::totalDensityString(), &m_totalDensity.value )
    .setPlotLevel( PlotLevel::LEVEL_0 );
  registerWrapper( viewKeyStruct::dTotalDensity_dPressureString(), &m_totalDensity.dPres ).
    setRestartFlags( RestartFlags::NO_WRITE );
  registerWrapper( viewKeyStruct::dTotalDensity_dTemperatureString(), &m_totalDensity.dTemp ).
    setRestartFlags( RestartFlags::NO_WRITE );
  registerWrapper( viewKeyStruct::dTotalDensity_dGlobalCompFractionString(), &m_totalDensity.dComp ).
    setRestartFlags( RestartFlags::NO_WRITE );

  registerWrapper( viewKeyStruct::initialTotalMassDensityString(), &m_initialTotalMassDensity );

  registerWrapper( viewKeyStruct::useMassString(), &m_useMass ).
    setRestartFlags( RestartFlags::NO_WRITE );

}

void MultiFluidBase::resizeFields( localIndex const size, localIndex const numPts )
{
  integer const numPhase = numFluidPhases();
  integer const numComp = numFluidComponents();

  m_phaseFraction.value.resize( size, numPts, numPhase );
  m_phaseFraction.dPres.resize( size, numPts, numPhase );
  m_phaseFraction.dTemp.resize( size, numPts, numPhase );
  m_phaseFraction.dComp.resize( size, numPts, numPhase, numComp );

  m_phaseDensity.value.resize( size, numPts, numPhase );
  m_phaseDensity.dPres.resize( size, numPts, numPhase );
  m_phaseDensity.dTemp.resize( size, numPts, numPhase );
  m_phaseDensity.dComp.resize( size, numPts, numPhase, numComp );

  m_phaseMassDensity.value.resize( size, numPts, numPhase );
  m_phaseMassDensity.dPres.resize( size, numPts, numPhase );
  m_phaseMassDensity.dTemp.resize( size, numPts, numPhase );
  m_phaseMassDensity.dComp.resize( size, numPts, numPhase, numComp );

  m_phaseViscosity.value.resize( size, numPts, numPhase );
  m_phaseViscosity.dPres.resize( size, numPts, numPhase );
  m_phaseViscosity.dTemp.resize( size, numPts, numPhase );
  m_phaseViscosity.dComp.resize( size, numPts, numPhase, numComp );

  m_phaseCompFraction.value.resize( size, numPts, numPhase, numComp );
  m_phaseCompFraction.dPres.resize( size, numPts, numPhase, numComp );
  m_phaseCompFraction.dTemp.resize( size, numPts, numPhase, numComp );
  m_phaseCompFraction.dComp.resize( size, numPts, numPhase, numComp, numComp );

  m_totalDensity.value.resize( size, numPts );
  m_totalDensity.dPres.resize( size, numPts );
  m_totalDensity.dTemp.resize( size, numPts );
  m_totalDensity.dComp.resize( size, numPts, numComp );

  m_initialTotalMassDensity.resize( size, numPts );
}

void MultiFluidBase::setLabels()
{
  getWrapper< array3d< real64, multifluid::LAYOUT_PHASE > >( viewKeyStruct::phaseFractionString() ).
    setDimLabels( 2, m_phaseNames );

  getWrapper< array3d< real64, multifluid::LAYOUT_PHASE > >( viewKeyStruct::phaseDensityString() ).
    setDimLabels( 2, m_phaseNames );

  getWrapper< array3d< real64, multifluid::LAYOUT_PHASE > >( viewKeyStruct::phaseMassDensityString() ).
    setDimLabels( 2, m_phaseNames );

  getWrapper< array3d< real64, multifluid::LAYOUT_PHASE > >( viewKeyStruct::phaseViscosityString() ).
    setDimLabels( 2, m_phaseNames );

  getWrapper< array4d< real64, multifluid::LAYOUT_PHASE_COMP > >( viewKeyStruct::phaseCompFractionString() ).
    setDimLabels( 2, m_phaseNames ).
    setDimLabels( 3, m_componentNames );
}

void MultiFluidBase::allocateConstitutiveData( dataRepository::Group & parent,
                                               localIndex const numConstitutivePointsPerParentIndex )
{
  ConstitutiveBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );
  resizeFields( parent.size(), numConstitutivePointsPerParentIndex );
}

void MultiFluidBase::postProcessInput()
{
  ConstitutiveBase::postProcessInput();

  integer const numComp = numFluidComponents();
  integer const numPhase = numFluidPhases();

  GEOSX_THROW_IF_LT_MSG( numComp, 1,
                         GEOSX_FMT( "{}: invalid number of components", getFullName() ),
                         InputError );
  GEOSX_THROW_IF_GT_MSG( numComp, MAX_NUM_COMPONENTS,
                         GEOSX_FMT( "{}: invalid number of components", getFullName() ),
                         InputError );
  GEOSX_THROW_IF_LT_MSG( numPhase, 1,
                         GEOSX_FMT( "{}: invalid number of phases", getFullName() ),
                         InputError );
  GEOSX_THROW_IF_GT_MSG( numPhase, MAX_NUM_PHASES,
                         GEOSX_FMT( "{}: invalid number of phases", getFullName() ),
                         InputError );
  GEOSX_THROW_IF_NE_MSG( m_componentMolarWeight.size(), numComp,
                         GEOSX_FMT( "{}: invalid number of values in attribute '{}'", getFullName(), viewKeyStruct::componentMolarWeightString() ),
                         InputError );

  // call to correctly set member array tertiary sizes on the 'main' material object
  resizeFields( 0, 0 );

  // set labels on array wrappers for plottable fields
  setLabels();
}

} //namespace constitutive

} //namespace geosx
