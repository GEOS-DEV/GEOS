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
 * @file DiffusionBase.cpp
 */

#include "constitutive/diffusion/DiffusionBase.hpp"
#include "constitutive/diffusion/DiffusionFields.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

DiffusionBase::DiffusionBase( string const & name, Group * const parent )
  : ConstitutiveBase( name, parent )
{
  registerWrapper( viewKeyStruct::phaseNamesString(), &m_phaseNames ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "List of fluid phases" );

  registerWrapper( viewKeyStruct::defaultPhaseDiffusivityMultiplierString(), &m_defaultPhaseDiffusivityMultiplier ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 1.0 ).
    setDescription( "List of phase diffusivity multipliers" );

  registerField< fields::diffusion::diffusivity >( &m_diffusivity );
  registerField< fields::diffusion::dDiffusivity_dTemperature >( &m_dDiffusivity_dTemperature );
  registerField< fields::diffusion::phaseDiffusivityMultiplier >( &m_phaseDiffusivityMultiplier );
}

void DiffusionBase::postInputInitialization()
{
  ConstitutiveBase::postInputInitialization();

  integer const numPhases = numFluidPhases();
  GEOS_THROW_IF_LT_MSG( numPhases, 2,
                        "invalid number of phases",
                        InputError, getDataContext() );
  GEOS_THROW_IF_GT_MSG( numPhases, MAX_NUM_PHASES,
                        "invalid number of phases",
                        InputError, getDataContext() );

  GEOS_THROW_IF( numPhases != m_defaultPhaseDiffusivityMultiplier.size(),
                 GEOS_FMT( "the arrays in `{}` and `{}` must have the same size",
                           viewKeyStruct::phaseNamesString(), viewKeyStruct::defaultPhaseDiffusivityMultiplierString() ),
                 InputError,
                 getWrapperDataContext( viewKeyStruct::phaseNamesString()),
                 getWrapperDataContext( viewKeyStruct::defaultPhaseDiffusivityMultiplierString()),
                 getDataContext() );
}

void DiffusionBase::allocateConstitutiveData( Group & parent, localIndex const numPts )
{
  // NOTE: enforcing 1 quadrature point
  m_diffusivity.resize( 0, 1, 3 );
  m_dDiffusivity_dTemperature.resize( 0, 1, 3 );
  m_phaseDiffusivityMultiplier.resize( 0, 1, 3 );

  ConstitutiveBase::allocateConstitutiveData( parent, numPts );

  for( localIndex ei = 0; ei < parent.size(); ++ei ) // TODO move into initializeState?
  {
    // NOTE: enforcing 1 quadrature point
    for( localIndex q = 0; q < 1; ++q )
    {
      for( integer ip = 0; ip < numFluidPhases(); ++ip )
      {
        m_phaseDiffusivityMultiplier[ei][q][ip] = m_defaultPhaseDiffusivityMultiplier[ip];
      }
    }
  }
}

void DiffusionBase::initializeTemperatureState( arrayView1d< real64 const > const & initialTemperature ) const
{
  GEOS_UNUSED_VAR( initialTemperature );

  localIndex const numE = m_phaseDiffusivityMultiplier.size( 0 );
  integer constexpr numQuad = 1; // NOTE: enforcing 1 quadrature point
  integer const numPhases = numFluidPhases();

  auto phaseDiffusivityMultiplierView = m_phaseDiffusivityMultiplier.toView();
  auto const defaultPhaseDiffusivityMultiplierView = m_defaultPhaseDiffusivityMultiplier.toViewConst();

  forAll< parallelDevicePolicy<> >( numE, [=] GEOS_HOST_DEVICE ( localIndex const ei )
  {
    for( localIndex q = 0; q < numQuad; ++q )
    {
      for( integer ip = 0; ip < numPhases; ++ip )
      {
        phaseDiffusivityMultiplierView[ei][q][ip] = defaultPhaseDiffusivityMultiplierView[ip];
      }
    }
  } );
}

} // namespace constitutive

} // namespace geos
