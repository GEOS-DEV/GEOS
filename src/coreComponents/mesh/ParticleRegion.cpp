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
  ParticleRegionBase( name, parent )
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
    ParticleBlockABC & source = particleBlocks.getGroup< ParticleBlockABC >( subRegion.getName() );
    subRegion.copyFromParticleBlock( source );
  }
}

// TODO This isn't great because really this calculation should only happen once, classic speed/memory trade-off. It just seems silly having multiple versions of the particle coordinates owned by the manager, regions and subregions.
array2d< real64 > ParticleRegion::getParticleCoordinates() const
{
  int size = this->size();
  array2d< real64 > coords(size, 3);
  int index = 0;
  this->forParticleSubRegions( [&]( auto & subRegion )
  {
    arrayView2d< real64 const > particleCoords = subRegion.getParticleCenter(); // I don't understand why the various consts are needed here, but it works so... TODO
    for(int i=0; i<subRegion.size(); i++)
    {
      for(int j=0; j<3; j++)
      {
        coords[index][j] = particleCoords[i][j];
      }
      index++;
    }
  } );
  return coords;
}

REGISTER_CATALOG_ENTRY( ObjectManagerBase, ParticleRegion, string const &, Group * const )

} /* namespace geosx */
