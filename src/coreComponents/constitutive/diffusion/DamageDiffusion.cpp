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
 * @file DamageDiffusion.cpp
 */

#include "DamageDiffusion.hpp"
#include "DiffusionFields.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

DamageDiffusion::DamageDiffusion( string const & name, Group * const parent ):
  DiffusionBase( name, parent ),
  m_defaultBulkDiffusivity( 0.0 )
{
  registerWrapper( viewKeyStruct::defaultBulkDiffusivityString(), &m_defaultBulkDiffusivity ).
    setInputFlag( InputFlags::REQUIRED ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Default diffusivity of the intact bulk material [m^2/s]" );

  registerWrapper( viewKeyStruct::damageDependenceConstantString(), &m_damageDependenceConstant ).
    setInputFlag( InputFlags::REQUIRED ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "xx, yy and zz exponents controlling how diffusivity scales with damage in each direction: "
                    "D_dim = bulkDiffusivity * exp( alpha_dim * damage )" );

  registerField< fields::diffusion::bulkDiffusivity >( &m_bulkDiffusivity );
}

void DamageDiffusion::allocateConstitutiveData( Group & parent, localIndex const numPts )
{
  DiffusionBase::allocateConstitutiveData( parent, numPts );

  // Initialize per-element bulk diffusivity from the XML default, and the
  // diffusivity tensor to the corresponding undamaged (bulk) value.
  for( localIndex ei = 0; ei < parent.size(); ++ei )
  {
    m_bulkDiffusivity[ei] = m_defaultBulkDiffusivity;

    for( localIndex q = 0; q < 1; ++q )
    {
      m_diffusivity[ei][q][0] = m_defaultBulkDiffusivity;
      m_diffusivity[ei][q][1] = m_defaultBulkDiffusivity;
      m_diffusivity[ei][q][2] = m_defaultBulkDiffusivity;
    }
  }
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, DamageDiffusion, string const &, Group * const )

} // namespace constitutive
} // namespace geos
