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
 * @file ParticleFluidBase.cpp
 */

#include "ParticleFluidBase.hpp"

namespace geosx
{

using namespace dataRepository;

namespace constitutive
{

ParticleFluidBase::ParticleFluidBase( string const & name, Group * const parent )
  : ConstitutiveBase( name, parent )
{

  registerWrapper( viewKeyStruct::settlingFactorString(), &m_settlingFactor );
  registerWrapper( viewKeyStruct::dSettlingFactor_dPressureString(), &m_dSettlingFactor_dPressure );
  registerWrapper( viewKeyStruct::dSettlingFactor_dProppantConcentrationString(), &m_dSettlingFactor_dProppantConcentration );
  registerWrapper( viewKeyStruct::dSettlingFactor_dComponentConcentrationString(), &m_dSettlingFactor_dComponentConcentration );

  registerWrapper( viewKeyStruct::collisionFactorString(), &m_collisionFactor );
  registerWrapper( viewKeyStruct::dCollisionFactor_dProppantConcentrationString(), &m_dCollisionFactor_dProppantConcentration );

  registerWrapper( viewKeyStruct::maxProppantConcentrationString(), &m_maxProppantConcentration ).
    setApplyDefaultValue( 0.6 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Max proppant concentration" );

  registerWrapper( viewKeyStruct::isCollisionalSlipString(), &m_isCollisionalSlip ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Whether the collisional component of the slip velocity is considered" );

  registerWrapper( viewKeyStruct::proppantPackPermeabilityString(), &m_proppantPackPermeability );

}

ParticleFluidBase::~ParticleFluidBase() = default;

void ParticleFluidBase::postProcessInput()
{
  ConstitutiveBase::postProcessInput();
}

void ParticleFluidBase::allocateConstitutiveData( Group & parent,
                                                  localIndex const numConstitutivePointsPerParentIndex )
{
  ConstitutiveBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  this->resize( parent.size() );
  m_dSettlingFactor_dComponentConcentration.resize( parent.size(), MAX_NUM_COMPONENTS );
}

} //namespace constitutive

} //namespace geosx
