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
  registerExtrinsicData( extrinsicMeshData::multifluid::dPhaseFraction{}, &m_phaseFraction.derivs );

  registerExtrinsicData( extrinsicMeshData::multifluid::phaseDensity{}, &m_phaseDensity.value );
  registerExtrinsicData( extrinsicMeshData::multifluid::dPhaseDensity{}, &m_phaseDensity.derivs );

  registerExtrinsicData( extrinsicMeshData::multifluid::phaseMassDensity{}, &m_phaseMassDensity.value );
  registerExtrinsicData( extrinsicMeshData::multifluid::dPhaseMassDensity{}, &m_phaseMassDensity.derivs );

  registerExtrinsicData( extrinsicMeshData::multifluid::phaseViscosity{}, &m_phaseViscosity.value );
  registerExtrinsicData( extrinsicMeshData::multifluid::dPhaseViscosity{}, &m_phaseViscosity.derivs );

  registerExtrinsicData( extrinsicMeshData::multifluid::phaseEnthalpy{}, &m_phaseEnthalpy.value );
  registerExtrinsicData( extrinsicMeshData::multifluid::dPhaseEnthalpy{}, &m_phaseEnthalpy.derivs );

  registerExtrinsicData( extrinsicMeshData::multifluid::phaseInternalEnergy{}, &m_phaseInternalEnergy.value );
  registerExtrinsicData( extrinsicMeshData::multifluid::dPhaseInternalEnergy{}, &m_phaseInternalEnergy.derivs );

  registerExtrinsicData( extrinsicMeshData::multifluid::phaseCompFraction{}, &m_phaseCompFraction.value );
  registerExtrinsicData( extrinsicMeshData::multifluid::dPhaseCompFraction{}, &m_phaseCompFraction.derivs );

  registerExtrinsicData( extrinsicMeshData::multifluid::totalDensity{}, &m_totalDensity.value );
  registerExtrinsicData( extrinsicMeshData::multifluid::dTotalDensity{}, &m_totalDensity.derivs );

  registerExtrinsicData( extrinsicMeshData::multifluid::initialTotalMassDensity{}, &m_initialTotalMassDensity );

}

void MultiFluidBase::resizeFields( localIndex const size, localIndex const numPts )
{
  integer const numPhase = numFluidPhases();
  integer const numComp = numFluidComponents();
  integer const numDof = numComp + 2;

  m_phaseFraction.value.resize( size, numPts, numPhase );
  m_phaseFraction.derivs.resize( size, numPts, numPhase, numDof );

  m_phaseDensity.value.resize( size, numPts, numPhase );
  m_phaseDensity.derivs.resize( size, numPts, numPhase, numDof );

  m_phaseMassDensity.value.resize( size, numPts, numPhase );
  m_phaseMassDensity.derivs.resize( size, numPts, numPhase, numDof );

  m_phaseViscosity.value.resize( size, numPts, numPhase );
  m_phaseViscosity.derivs.resize( size, numPts, numPhase, numDof );

  m_phaseEnthalpy.value.resize( size, numPts, numPhase );
  m_phaseEnthalpy.derivs.resize( size, numPts, numPhase, numDof );

  m_phaseInternalEnergy.value.resize( size, numPts, numPhase );
  m_phaseInternalEnergy.derivs.resize( size, numPts, numPhase, numDof );

  m_phaseCompFraction.value.resize( size, numPts, numPhase, numComp );
  m_phaseCompFraction.derivs.resize( size, numPts, numPhase, numComp, numDof );

  m_totalDensity.value.resize( size, numPts );
  m_totalDensity.derivs.resize( size, numPts, numDof );

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

  getExtrinsicData< extrinsicMeshData::multifluid::phaseEnthalpy >().
    setDimLabels( 2, m_phaseNames );

  getExtrinsicData< extrinsicMeshData::multifluid::phaseInternalEnergy >().
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

} // namespace constitutive

} // namespace geosx
