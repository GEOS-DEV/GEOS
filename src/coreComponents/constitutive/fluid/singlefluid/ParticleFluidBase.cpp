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
 * @file ParticleFluidBase.cpp
 */

#include "ParticleFluidBase.hpp"

#include "ParticleFluidFields.hpp"
#include "mesh/ElementSubRegionBase.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

ParticleFluidBase::ParticleFluidBase( string const & name, Group * const parent )
  : ConstitutiveBase( name, parent )
{
  registerWrapper( viewKeyStruct::maxProppantConcentrationString(), &m_maxProppantConcentration ).
    setApplyDefaultValue( 0.6 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Max proppant concentration" );

  registerWrapper( viewKeyStruct::isCollisionalSlipString(), &m_isCollisionalSlip ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Whether the collisional component of the slip velocity is considered" );
}

void ParticleFluidBase::allocateConstitutiveData( dataRepository::Group & parent,
                                                  localIndex const numConstitutivePointsPerParentIndex )
{
  ConstitutiveBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  auto subregion = dynamic_cast< ElementSubRegionBase * >(&parent); // TODO remove

  subregion->registerField< fields::particlefluid::settlingFactor >( getName(), &m_settlingFactor );
  subregion->registerField< fields::particlefluid::dSettlingFactor_dPressure >( getName(), &m_dSettlingFactor_dPressure );
  subregion->registerField< fields::particlefluid::dSettlingFactor_dProppantConcentration >( getName(), &m_dSettlingFactor_dProppantConcentration );
  subregion->registerField< fields::particlefluid::dSettlingFactor_dComponentConcentration >( getName(), &m_dSettlingFactor_dComponentConcentration ).
    reference().resizeDimension< 1 >( MAX_NUM_COMPONENTS );

  subregion->registerField< fields::particlefluid::collisionFactor >( getName(), &m_collisionFactor );
  subregion->registerField< fields::particlefluid::dCollisionFactor_dProppantConcentration >( getName(), &m_dCollisionFactor_dProppantConcentration );

  subregion->registerField< fields::particlefluid::proppantPackPermeability >( getName(), &m_proppantPackPermeability );
}

} //namespace constitutive

} //namespace geos
