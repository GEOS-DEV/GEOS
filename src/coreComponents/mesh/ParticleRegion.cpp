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
#include "common/MpiWrapper.hpp"

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
    int const mpiRank = MpiWrapper::commRank( MPI_COMM_GEOSX );
    subRegion.setGhostRank( mpiRank, subRegion.size() );
  }
}

// TODO This isn't great because really this calculation should only happen once, classic speed/memory trade-off. It just seems silly having multiple versions of the particle coordinates owned by the manager, regions and subregions.
// TODO This should ABSOLUTELY be changed to call a ParticleSubRegion::getParticleCoordinates (and/or getParticleCorners) on each subregion such that we can access those functions directly if needed
array2d< real64 > ParticleRegion::getParticleCoordinates() const
{
  int size = 8*(this->size()); // number of particle corners in this region
  array2d< real64 > coords(size, 3);
  int index = 0;
  int signs[8][3] = { { 1,  1,  1},
                      { 1,  1, -1},
                      { 1, -1,  1},
                      { 1, -1, -1},
                      {-1,  1,  1},
                      {-1,  1, -1},
                      {-1, -1,  1},
                      {-1, -1, -1} };
  this->forParticleSubRegions( [&]( auto & subRegion )
  {
    arrayView2d< real64 const > particleCenter = subRegion.getParticleCenter();
    if(!subRegion.getHasRVectors()) // if no r-vectors, return coordinates of a cube centered at the particle with side length = volume^(1/3)
    {
      arrayView1d< real64 const > particleVolume = subRegion.getParticleVolume();
      for(int p=0; p<subRegion.size(); p++)
      {
        real64 a = 0.5*std::pow(particleVolume[p],1.0/3.0); // cube half-side-length
        for(int corner=0; corner<8; corner++)
        {
          for(int i=0; i<3; i++)
          {
            coords[index][i] = particleCenter[p][i] + signs[corner][i]*a;
          }
          index++;
        }
      }
    }
    else
    {
      arrayView3d< real64 const > particleRVectors = subRegion.getParticleRVectors();
      for(int p=0; p<subRegion.size(); p++)
      {
        for(int corner=0; corner<8; corner++)
        {
          for(int i=0; i<3; i++)
          {
            coords[index][i] = particleCenter[p][i] + signs[corner][0]*particleRVectors[p][0][i] + signs[corner][1]*particleRVectors[p][1][i] + signs[corner][2]*particleRVectors[p][2][i];
          }
          index++;
        }
      }
    }
  } );
  return coords;
}

REGISTER_CATALOG_ENTRY( ObjectManagerBase, ParticleRegion, string const &, Group * const )

} /* namespace geosx */
