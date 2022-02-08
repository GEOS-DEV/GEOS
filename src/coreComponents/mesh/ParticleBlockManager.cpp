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
 * @file ParticleBlockManager.cpp
 */

#include "ParticleBlockManager.hpp"

#include <map>
#include <vector>

namespace geosx
{
using namespace dataRepository;

ParticleBlockManager::ParticleBlockManager( string const & name, Group * const parent ):
  ObjectManagerBase( name, parent )
{
  this->registerGroup< Group >( keys::particleBlocks );
}

ParticleBlockManager::~ParticleBlockManager()
{
  // TODO Auto-generated destructor stub
}

void ParticleBlockManager::resize( integer_array const & numParticles,
                               string_array const & regionNames )
{
  localIndex const numRegions = LvArray::integerConversion< localIndex >( regionNames.size());
  for( localIndex reg=0; reg<numRegions; ++reg )
  {
    this->getRegion( regionNames[reg] ).resize( numParticles[reg] );
  }
}

Group * ParticleBlockManager::createChild( string const & GEOSX_UNUSED_PARAM( childKey ), string const & GEOSX_UNUSED_PARAM( childName ) )
{
  return nullptr;
}



REGISTER_CATALOG_ENTRY( ObjectManagerBase, ParticleBlockManager, string const &, Group * const )
}
