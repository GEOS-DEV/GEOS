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
#include "MultiFluidExtrinsicData.hpp"

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

  registerWrapper( viewKeyStruct::useMassString(), &m_useMass ).
    setRestartFlags( RestartFlags::NO_WRITE );

  registerExtrinsicData( extrinsicMeshData::multifluid::phaseFraction{}, &m_phaseFraction.value );
  registerExtrinsicData( extrinsicMeshData::multifluid::dPhaseFraction_dPressure{}, &m_phaseFraction.dPres );
  registerExtrinsicData( extrinsicMeshData::multifluid::dPhaseFraction_dTemperature{}, &m_phaseFraction.dTemp );
  registerExtrinsicData( extrinsicMeshData::multifluid::dPhaseFraction_dGlobalCompFraction{}, &m_phaseFraction.dComp );

  registerExtrinsicData( extrinsicMeshData::multifluid::phaseDensity{}, &m_phaseDensity.value );
  registerExtrinsicData( extrinsicMeshData::multifluid::dPhaseDensity_dPressure{}, &m_phaseDensity.dPres );
  registerExtrinsicData( extrinsicMeshData::multifluid::dPhaseDensity_dTemperature{}, &m_phaseDensity.dTemp );
  registerExtrinsicData( extrinsicMeshData::multifluid::dPhaseDensity_dGlobalCompFraction{}, &m_phaseDensity.dComp );

  registerExtrinsicData( extrinsicMeshData::multifluid::phaseMassDensity{}, &m_phaseMassDensity.value );
  registerExtrinsicData( extrinsicMeshData::multifluid::dPhaseMassDensity_dPressure{}, &m_phaseMassDensity.dPres );
  registerExtrinsicData( extrinsicMeshData::multifluid::dPhaseMassDensity_dTemperature{}, &m_phaseMassDensity.dTemp );
  registerExtrinsicData( extrinsicMeshData::multifluid::dPhaseMassDensity_dGlobalCompFraction{}, &m_phaseMassDensity.dComp );

  registerExtrinsicData( extrinsicMeshData::multifluid::phaseViscosity{}, &m_phaseViscosity.value );
  registerExtrinsicData( extrinsicMeshData::multifluid::dPhaseViscosity_dPressure{}, &m_phaseViscosity.dPres );
  registerExtrinsicData( extrinsicMeshData::multifluid::dPhaseViscosity_dTemperature{}, &m_phaseViscosity.dTemp );
  registerExtrinsicData( extrinsicMeshData::multifluid::dPhaseViscosity_dGlobalCompFraction{}, &m_phaseViscosity.dComp );

  registerExtrinsicData( extrinsicMeshData::multifluid::phaseCompFraction{}, &m_phaseCompFraction.value );
  registerExtrinsicData( extrinsicMeshData::multifluid::dPhaseCompFraction_dPressure{}, &m_phaseCompFraction.dPres );
  registerExtrinsicData( extrinsicMeshData::multifluid::dPhaseCompFraction_dTemperature{}, &m_phaseCompFraction.dTemp );
  registerExtrinsicData( extrinsicMeshData::multifluid::dPhaseCompFraction_dGlobalCompFraction{}, &m_phaseCompFraction.dComp );

  registerExtrinsicData( extrinsicMeshData::multifluid::totalDensity{}, &m_totalDensity.value );
  registerExtrinsicData( extrinsicMeshData::multifluid::dTotalDensity_dPressure{}, &m_totalDensity.dPres );
  registerExtrinsicData( extrinsicMeshData::multifluid::dTotalDensity_dTemperature{}, &m_totalDensity.dTemp );
  registerExtrinsicData( extrinsicMeshData::multifluid::dTotalDensity_dGlobalCompFraction{}, &m_totalDensity.dComp );

  registerExtrinsicData( extrinsicMeshData::multifluid::initialTotalMassDensity{}, &m_initialTotalMassDensity );

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
  getExtrinsicData< extrinsicMeshData::multifluid::phaseFraction >().
    setDimLabels( 2, m_phaseNames );

  getExtrinsicData< extrinsicMeshData::multifluid::phaseDensity >().
    setDimLabels( 2, m_phaseNames );

  getExtrinsicData< extrinsicMeshData::multifluid::phaseMassDensity >().
    setDimLabels( 2, m_phaseNames );

  getExtrinsicData< extrinsicMeshData::multifluid::phaseViscosity >().
    setDimLabels( 2, m_phaseNames );

  getExtrinsicData< extrinsicMeshData::multifluid::phaseCompFraction >().
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
