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
 * @file WillisRichardsPermeability.cpp
 */

#include "WillisRichardsPermeability.hpp"
#include "constitutive/permeability/PermeabilityFields.hpp"
#include "mesh/ElementSubRegionBase.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{


WillisRichardsPermeability::WillisRichardsPermeability( string const & name, Group * const parent ):
  PermeabilityBase( name, parent )
{
  registerWrapper( viewKeyStruct::maxFracApertureString(), &m_maxFracAperture ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Maximum fracture aperture at zero contact stress." );

  registerWrapper( viewKeyStruct::dilationCoefficientString(), &m_dilationCoefficient ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Dilation coefficient (tan of dilation angle)." );

  registerWrapper( viewKeyStruct::refClosureStressString(), &m_refClosureStress ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Effective normal stress causes 90% reduction in aperture." );
}

void WillisRichardsPermeability::allocateConstitutiveData( Group & parent,
                                                           localIndex const numConstitutivePointsPerParentIndex )
{
  PermeabilityBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  auto subregion = dynamic_cast< ElementSubRegionBase * >( &parent ); // TODO remove

  // NOTE: enforcing 1 quadrature point
  subregion->registerField< fields::permeability::dPerm_dDispJump >( getName(), &m_dPerm_dDispJump ).
    reference().resizeDimension< 1, 2, 3 >( 1, 3, 3 );
  subregion->registerField< fields::permeability::dPerm_dTraction >( getName(), &m_dPerm_dTraction ).
    reference().resizeDimension< 1, 2, 3 >( 1, 3, 3 );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, WillisRichardsPermeability, string const &, Group * const )

} /* namespace constitutive */
} /* namespace geos */
