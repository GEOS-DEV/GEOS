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
 * @file SlipDependentPermeability.cpp
 */

#include "SlipDependentPermeability.hpp"
#include "LvArray/src/tensorOps.hpp"
#include "constitutive/permeability/PermeabilityExtrinsicData.hpp"

namespace geosx
{

using namespace dataRepository;

namespace constitutive
{


SlipDependentPermeability::SlipDependentPermeability( string const & name, Group * const parent ):
  PermeabilityBase( name, parent )
{
  registerWrapper( viewKeyStruct::shearDispThresholdString(), &m_shearDispThreshold ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Threshold of shear displacement." );

  registerWrapper( viewKeyStruct::maxPermMultiplierString(), &m_maxPermMultiplier ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Maximum permeability multiplier." );

  registerWrapper( viewKeyStruct::initialPermeabilityString(), &m_initialPermeability ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( " initial permeability of the fracture." );

  registerField( fields::permeability::dPerm_dDispJump{}, &m_dPerm_dDispJump );
}

std::unique_ptr< ConstitutiveBase >
SlipDependentPermeability::deliverClone( string const & name,
                                         Group * const parent ) const
{
  return ConstitutiveBase::deliverClone( name, parent );
}

void SlipDependentPermeability::allocateConstitutiveData( dataRepository::Group & parent,
                                                          localIndex const numConstitutivePointsPerParentIndex )
{
// NOTE: enforcing 1 quadrature point
  m_dPerm_dDispJump.resize( 0, 1, 3, 3 );

  PermeabilityBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, SlipDependentPermeability, string const &, Group * const )

} /* namespace constitutive */
} /* namespace geosx */
