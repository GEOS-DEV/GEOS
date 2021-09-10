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
 * @file StrainDependentPermeability.cpp
 */

#include "StrainDependentPermeability.hpp"

namespace geosx
{

using namespace dataRepository;

namespace constitutive
{


StrainDependentPermeability::StrainDependentPermeability( string const & name, Group * const parent ):
  PermeabilityBase( name, parent )
{
  registerWrapper( viewKeyStruct::strainThresholdString(), &m_strainThreshold ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Threshold of shear strain." );

  registerWrapper( viewKeyStruct::maxPermMultiplierString(), &m_maxPermMultiplier ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Maximum permeability multiplier." );

  registerWrapper( viewKeyStruct::iniPermeabilityString(), &m_permeabilityComponents ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Initial permeability tensor of the fault." );

  registerWrapper( viewKeyStruct::dPerm_dStrainString(), &m_dPerm_dStrain );
}

std::unique_ptr< ConstitutiveBase >
StrainDependentPermeability::deliverClone( string const & name,
                                           Group * const parent ) const
{
  return ConstitutiveBase::deliverClone( name, parent );
}

void StrainDependentPermeability::allocateConstitutiveData( dataRepository::Group & parent,
                                                            localIndex const numConstitutivePointsPerParentIndex )
{
// NOTE: enforcing 1 quadrature point
  m_dPerm_dStrain.resize( 0, 1, 3 );

  PermeabilityBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, StrainDependentPermeability, string const &, Group * const )

}
} /* namespace geosx */
