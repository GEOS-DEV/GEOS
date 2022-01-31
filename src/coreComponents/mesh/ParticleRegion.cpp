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
 * @file ParticleRegion.cpp
 */

#include "ParticleRegion.hpp"
#include "common/TimingMacros.hpp"
#include "LvArray/src/SparsityPattern.hpp"

namespace geosx
{
using namespace dataRepository;

ParticleRegion::ParticleRegion( string const & name, Group * const parent ):
  ElementRegionBase( name, parent )
{
  registerWrapper( viewKeyStruct::sourceParticleBlockNamesString(), &m_particleBlockNames ).
    setInputFlag( InputFlags::OPTIONAL );
}

ParticleRegion::~ParticleRegion()
{}

void ParticleRegion::generateMesh( Group & particleBlocks )
{
  Group & particleSubRegions = this->getGroup( viewKeyStruct::particleSubRegions() );

  for( string const & particleBlockName : this->m_particleBlockNames )
  {
    ParticleSubRegion & subRegion = particleSubRegions.registerGroup< ParticleSubRegion >( particleBlockName );
    ParticleBlock & source = particleBlocks.getGroup< ParticleBlock >( subRegion.getName() );
    subRegion.copyFromParticleBlock( source );
  }
}

REGISTER_CATALOG_ENTRY( ObjectManagerBase, ParticleRegion, string const &, Group * const )

} /* namespace geosx */
