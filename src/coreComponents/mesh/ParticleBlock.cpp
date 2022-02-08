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
 * @file ParticleBlock.cpp
 *
 */

#include "ParticleBlock.hpp"
#include "MeshLevel.hpp"

#include "common/GEOS_RAJA_Interface.hpp"

namespace geosx
{
using namespace dataRepository;

ParticleBlock::ParticleBlock( string const & name, Group * const parent ):
  ParticleSubRegionBase( name, parent ),
  m_externalPropertyNames()
{

}

ParticleBlock::~ParticleBlock()
{}

void ParticleBlock::setParticleType( ParticleType const particleType )
{
  m_particleType = particleType;

  switch( m_particleType )
  {
    case ParticleType::SinglePoint:
    {
      // Single Point
      this->setHasRVectors( false );
      break;
    }
    case ParticleType::CPDI:
    {
      // CPDI
      this->setHasRVectors( true );
      break;
    }
    case ParticleType::CPTI:
    {
      // CPTI
      this->setHasRVectors( true );
      break;
    }
    case ParticleType::CPDI2:
    {
      // CPDI2
      this->setHasRVectors( false );
      break;
    }
    default:
    {
      GEOSX_ERROR( "Invalid particle type: " << m_particleType );
    }
  }
}

REGISTER_CATALOG_ENTRY( ObjectManagerBase, ParticleBlock, string const &, Group * const )

}
