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
 * @file StrainDependentPermeability.cpp
 */

#include "StrainDependentPermeability.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{


StrainDependentPermeability::StrainDependentPermeability( string const & name, Group * const parent ):
  PermeabilityBase( name, parent )
{
  registerWrapper( viewKeyStruct::empiricalConstantString(), &m_empiricalConstant ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "an empirical constant." );

  registerField< fields::permeability::dPerm_dVolStrain >( &m_dPerm_dVolStrain );
}

std::unique_ptr< ConstitutiveBase >
StrainDependentPermeability::deliverClone( string const & name,
                                           Group * const parent ) const
{
  return PermeabilityBase::deliverClone( name, parent );
}

void StrainDependentPermeability::allocateConstitutiveData( Group & parent,
                                                            localIndex const numPts )
{
  // NOTE: enforcing 1 quadrature point
  m_dPerm_dVolStrain.resize( 0, 1, 3 );

  PermeabilityBase::allocateConstitutiveData( parent, numPts );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, StrainDependentPermeability, string const &, Group * const )

}
} /* namespace geos */
