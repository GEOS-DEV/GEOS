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
 * @file ConstantDiffusion.cpp
 */

#include "ConstantDiffusion.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

ConstantDiffusion::ConstantDiffusion( string const & name, Group * const parent ):
  DiffusionBase( name, parent )
{
  registerWrapper( viewKeyStruct::diffusivityComponentsString(), &m_diffusivityComponents ).
    setInputFlag( InputFlags::REQUIRED ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "xx, yy, and zz components of a diffusivity tensor [m^2/s]" );
}

std::unique_ptr< ConstitutiveBase >
ConstantDiffusion::deliverClone( string const & name,
                                 Group * const parent ) const
{
  return DiffusionBase::deliverClone( name, parent );
}

void ConstantDiffusion::allocateConstitutiveData( dataRepository::Group & parent,
                                                  localIndex const numConstitutivePointsPerParentIndex )
{
  DiffusionBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  for( localIndex ei = 0; ei < parent.size(); ++ei )
  {
    for( localIndex ip = 0; ip < numFluidPhases(); ++ip )
    {

      // NOTE: enforcing 1 quadrature point
      for( localIndex q = 0; q < 1; ++q )
      {
        //TODO change input to have diag tensor by phase
        m_diffusivity[ei][q][ip][0] = m_diffusivityComponents[ip][0];
        m_diffusivity[ei][q][ip][1] = m_diffusivityComponents[ip][1];
        m_diffusivity[ei][q][ip][2] = m_diffusivityComponents[ip][2];
      }
    }
  }
}

void ConstantDiffusion::postProcessInput()
{
  GEOS_THROW_IF( m_diffusivityComponents.size( 0 ) != m_phaseNames.size(),
                 GEOS_FMT( "{}: the size of the diffusivity must be equal to the number of phases ({})",
                           getFullName(), m_phaseNames.size() ),
                 InputError );

  GEOS_THROW_IF( m_diffusivityComponents.size( 1 ) != 3,
                 GEOS_FMT( "{}: the size of the diffusivity must be equal to the number of direction (3)",
                           getFullName()),
                 InputError );

  for( localIndex ip = 0; ip < numFluidPhases(); ++ip )
  {

    GEOS_THROW_IF( m_diffusivityComponents[ip][0] < 0 ||
                   m_diffusivityComponents[ip][1] < 0 ||
                   m_diffusivityComponents[ip][2] < 0,
                   GEOS_FMT( "{}: the components of the diffusivity tensor must be non-negative for phase {}",
                             getFullName(), m_phaseNames[ip] ),
                   InputError );
  }
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, ConstantDiffusion, string const &, Group * const )

} // namespace constitutive

} // namespace geos
