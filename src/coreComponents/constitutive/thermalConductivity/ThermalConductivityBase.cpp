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
 * @file ThermalConductivityBase.cpp
 */

#include "ThermalConductivityBase.hpp"
#include "ThermalConductivityExtrinsicData.hpp"

namespace geosx
{

using namespace dataRepository;

namespace constitutive
{

ThermalConductivityBase::ThermalConductivityBase( string const & name, Group * const parent )
  : ConstitutiveBase( name, parent )
{
  registerWrapper( viewKeyStruct::phaseNamesString(), &m_phaseNames ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "List of fluid phases" );

  registerExtrinsicData( extrinsicMeshData::thermalconductivity::effectiveConductivity{}, &m_effectiveConductivity );
  registerExtrinsicData( extrinsicMeshData::thermalconductivity::dEffectiveConductivity_dPhaseVolFraction{}, &m_dEffectiveConductivity_dPhaseVolFrac );
}

void ThermalConductivityBase::postProcessInput()
{
  ConstitutiveBase::postProcessInput();

  integer const numPhases = numFluidPhases();
  GEOSX_THROW_IF_LT_MSG( numPhases, 2,
                         GEOSX_FMT( "{}: invalid number of phases", getFullName() ),
                         InputError );
  GEOSX_THROW_IF_GT_MSG( numPhases, MAX_NUM_PHASES,
                         GEOSX_FMT( "{}: invalid number of phases", getFullName() ),
                         InputError );

  m_effectiveConductivity.resize( 0, 0, 3 );
  m_dEffectiveConductivity_dPhaseVolFrac.resize( 0, 0, 3, numPhases );
}

void ThermalConductivityBase::allocateConstitutiveData( dataRepository::Group & parent,
                                                        localIndex const numConstitutivePointsPerParentIndex )
{
  // NOTE: enforcing 1 quadrature point
  integer const numPhases = numFluidPhases();
  m_effectiveConductivity.resize( 0, 1, 3 );
  m_dEffectiveConductivity_dPhaseVolFrac.resize( 0, 1, 3, numPhases );

  ConstitutiveBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );
}

} // namespace constitutive

} // namespace geosx
