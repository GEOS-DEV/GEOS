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
 * @file VolumeWeightedThermalConductivity.cpp
 */

#include "VolumeWeightedThermalConductivity.hpp"

#include "ThermalConductivityExtrinsicData.hpp"

namespace geosx
{

using namespace dataRepository;

namespace constitutive
{

VolumeWeightedThermalConductivity::VolumeWeightedThermalConductivity( string const & name, Group * const parent ):
  ThermalConductivityBase( name, parent )
{
  registerWrapper( viewKeyStruct::defaultRockThermalConductivityString(), &m_defaultRockThermalConductivity ).
    setInputFlag( InputFlags::REQUIRED ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Default rock thermal conductivity [W/(m.K)]" );

  registerWrapper( viewKeyStruct::phaseThermalConductivityString(), &m_phaseThermalConductivity ).
    setInputFlag( InputFlags::REQUIRED ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Phase thermal conductivity [W/(m.K)]" );

  registerExtrinsicData( extrinsicMeshData::thermalconductivity::rockThermalConductivity{}, &m_rockThermalConductivity );
}

std::unique_ptr< ConstitutiveBase >
VolumeWeightedThermalConductivity::deliverClone( string const & name,
                                                 Group * const parent ) const
{
  return ThermalConductivityBase::deliverClone( name, parent );
}

void VolumeWeightedThermalConductivity::allocateConstitutiveData( dataRepository::Group & parent,
                                                                  localIndex const numConstitutivePointsPerParentIndex )
{
  // NOTE: enforcing 1 quadrature point
  m_rockThermalConductivity.resize( 0, 1 );

  ThermalConductivityBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  for( localIndex ei = 0; ei < parent.size(); ++ei )
  {
    for( localIndex q = 0; q < 1; ++q )
    {
      m_rockThermalConductivity[ei][q] = m_defaultRockThermalConductivity;
    }
  }
}

void VolumeWeightedThermalConductivity::postProcessInput()
{
  GEOSX_THROW_IF( m_defaultRockThermalConductivity <= 0,
                  GEOSX_FMT( "{}: the rock thermal conductivity must be strictly positive",
                             getFullName() ),
                  InputError );

  for( integer ip = 0; ip < numFluidPhases(); ++ip )
  {
    GEOSX_THROW_IF( m_phaseThermalConductivity[ip] <= 0,
                    GEOSX_FMT( "{}: the phase thermal conductivity for phase {} must be strictly positive",
                               getFullName(), ip ),
                    InputError );
  }


}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, VolumeWeightedThermalConductivity, string const &, Group * const )

} // namespace constitutive

} // namespace geosx
