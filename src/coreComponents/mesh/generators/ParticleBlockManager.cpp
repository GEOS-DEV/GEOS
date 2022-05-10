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

#include "ParticleBlockManager.hpp"

#include "mesh/generators/CellBlockUtilities.hpp"
#include "mesh/utilities/MeshMapUtilities.hpp"

#include <algorithm>

namespace geosx
{
using namespace dataRepository;

ParticleBlockManager::ParticleBlockManager( string const & name, Group * const parent ):
  ParticleBlockManagerABC( name, parent )
{
  this->registerGroup< Group >( viewKeyStruct::particleBlocks() );
}

void ParticleBlockManager::resize( integer_array const & numParticles,
                                   string_array const & regionNames )
{
  localIndex const numRegions = LvArray::integerConversion< localIndex >( regionNames.size());
  for( localIndex reg=0; reg<numRegions; ++reg )
  {
    this->getParticleBlock( regionNames[reg] ).resize( numParticles[reg] );
  }
}

Group * ParticleBlockManager::createChild( string const & GEOSX_UNUSED_PARAM( childKey ), string const & GEOSX_UNUSED_PARAM( childName ) )
{
  return nullptr;
}

ParticleBlock & ParticleBlockManager::registerParticleBlock( string name )
{
  return this->getParticleBlocks().registerGroup< ParticleBlock >( name );
}

const Group & ParticleBlockManager::getParticleBlocks() const
{
  return this->getGroup( viewKeyStruct::particleBlocks() );
}

Group & ParticleBlockManager::getParticleBlocks()
{
  return this->getGroup( viewKeyStruct::particleBlocks() );
}

const ParticleBlock & ParticleBlockManager::getParticleBlock( localIndex const blockIndex ) const
{
  return this->getParticleBlocks().getGroup< ParticleBlock >( blockIndex );
}

}
