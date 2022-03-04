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

#include "ParticleFluidExtrinsicData.hpp"

namespace geosx
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

  registerExtrinsicData( extrinsicMeshData::particlefluid::settlingFactor{}, &m_settlingFactor );
  registerExtrinsicData( extrinsicMeshData::particlefluid::dSettlingFactor_dPressure{}, &m_dSettlingFactor_dPressure );
  registerExtrinsicData( extrinsicMeshData::particlefluid::dSettlingFactor_dProppantConcentration{}, &m_dSettlingFactor_dProppantConcentration );
  registerExtrinsicData( extrinsicMeshData::particlefluid::dSettlingFactor_dComponentConcentration{}, &m_dSettlingFactor_dComponentConcentration );

  registerExtrinsicData( extrinsicMeshData::particlefluid::collisionFactor{}, &m_collisionFactor );
  registerExtrinsicData( extrinsicMeshData::particlefluid::dCollisionFactor_dProppantConcentration{}, &m_dCollisionFactor_dProppantConcentration );

  registerExtrinsicData( extrinsicMeshData::particlefluid::proppantPackPermeability{}, &m_proppantPackPermeability );
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
