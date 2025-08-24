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
#include "mesh/ElementSubRegionBase.hpp"

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
}

void DiffusionBase::postInputInitialization()
{
  ConstitutiveBase::postInputInitialization();

  integer const numPhases = numFluidPhases();
  GEOS_THROW_IF_LT_MSG( numPhases, 2,
                        GEOS_FMT( "{}: invalid number of phases", getFullName() ),
                        InputError );
  GEOS_THROW_IF_GT_MSG( numPhases, MAX_NUM_PHASES,
                        GEOS_FMT( "{}: invalid number of phases", getFullName() ),
                        InputError );

  GEOS_THROW_IF( numPhases != m_defaultPhaseDiffusivityMultiplier.size(),
                 GEOS_FMT( "{}: the arrays in `{}` and `{}` must have the same size",
                           getFullName(), viewKeyStruct::phaseNamesString(), viewKeyStruct::defaultPhaseDiffusivityMultiplierString() ),
                 InputError );
}

void DiffusionBase::allocateConstitutiveData( Group & parent, localIndex const numPts )
{
  auto subregion = dynamic_cast< ElementSubRegionBase * >(&parent); // TODO remove

  // NOTE: enforcing 1 quadrature point
  subregion->registerField< fields::diffusion::diffusivity >( getName(), &m_diffusivity ).
    reference().resizeDimension< 1, 2 >( 1, 3 );
  subregion->registerField< fields::diffusion::dDiffusivity_dTemperature >( getName(), &m_dDiffusivity_dTemperature ).
    reference().resizeDimension< 1, 2 >( 1, 3 );
  subregion->registerField< fields::diffusion::phaseDiffusivityMultiplier >( getName(), &m_phaseDiffusivityMultiplier ).
    reference().resizeDimension< 1, 2 >( 1, numFluidPhases());

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

} // namespace constitutive

} // namespace geos
