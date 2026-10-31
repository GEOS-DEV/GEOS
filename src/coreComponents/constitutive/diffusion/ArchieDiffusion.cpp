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
 * @file ArchieDiffusion.cpp
 */

#include "ArchieDiffusion.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

ArchieDiffusion::ArchieDiffusion( string const & name, Group * const parent ):
  DiffusionBase( name, parent ),
  m_fluidDiffusivity( 0.0 ),
  m_matrixCementationExponent( 0.0 ),
  m_fractureCementationExponent( 1.0 )
{
  registerWrapper( viewKeyStruct::fluidDiffusivityString(), &m_fluidDiffusivity ).
    setInputFlag( InputFlags::REQUIRED ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Diffusivity D_w of the pore fluid [m^2/s]; the porosity-weighted effective diffusivity is D_w * phi^m(d)" );

  registerWrapper( viewKeyStruct::matrixCementationExponentString(), &m_matrixCementationExponent ).
    setInputFlag( InputFlags::REQUIRED ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Archie cementation exponent m0 of the intact matrix" );

  registerWrapper( viewKeyStruct::fractureCementationExponentString(), &m_fractureCementationExponent ).
    setApplyDefaultValue( 1.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Archie cementation exponent mc of the fully damaged material (open fracture); "
                    "the exponent is m(d) = m0 - ( m0 - mc ) * damage" );
}

void ArchieDiffusion::postInputInitialization()
{
  // Single-phase model, as DamageDiffusion: skip DiffusionBase::postInputInitialization(), which requires numPhases >= 2.
  ConstitutiveBase::postInputInitialization();

  integer const numPhases = numFluidPhases();
  GEOS_THROW_IF_LT_MSG( numPhases, 1,
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

  GEOS_THROW_IF( m_fractureCementationExponent <= 0.0 || m_fractureCementationExponent > m_matrixCementationExponent,
                 GEOS_FMT( "`{}` must be in ( 0, `{}` ]",
                           viewKeyStruct::fractureCementationExponentString(), viewKeyStruct::matrixCementationExponentString() ),
                 InputError, getDataContext() );
}

void ArchieDiffusion::allocateConstitutiveData( Group & parent, localIndex const numPts )
{
  DiffusionBase::allocateConstitutiveData( parent, numPts );

  // Placeholder until the coupled solid sets D_w * phi^( m0 - 1 ) from the initial porosity (initializeState)
  m_diffusivity.setValues< serialPolicy >( m_fluidDiffusivity );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, ArchieDiffusion, string const &, Group * const )

} // namespace constitutive
} // namespace geos
