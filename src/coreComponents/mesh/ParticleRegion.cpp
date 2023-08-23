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

namespace geos
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

    // Copy data from particle blocks into subregions. This also sets particle type.
    subRegion.copyFromParticleBlock( source );

    // Set the rank of particles on each subregion
    int const mpiRank = MpiWrapper::commRank( MPI_COMM_GEOSX );
    subRegion.setParticleRank( mpiRank );

    // Set the number of vertices of the particles on each subregion
    subRegion.setNumberOfVerticesPerParticle( subRegion.getParticleType() );
  }
}

// TODO This should be changed to call a ParticleSubRegion::getParticleCoordinates (and/or getParticleCorners) on each subregion such that
// we can access those functions directly if needed
array2d< real64 > ParticleRegion::getParticleCorners() const
{
  GEOS_LOG_RANK_0("Particle Region Size " << this->getName() << ", size " << this->size());
  
  // CC: Not sure why size of particle region is used instead of number of particles
  // size seems related to group parent object and number of registered wrappers but someone holds size of particles in region
  // TODO investigate this
  // int const size = 8 * ( this->size() ); // number of particle corners in this region for a parallelpiped CPDI particle type (TODO add
                                         // support for general types)
  int const size = 8 * ( this->getNumberOfParticles() );
  array2d< real64 > coords( size, 3 );
  int index = 0;
  int const signs[8][3] = { { 1, 1, 1},   // These have to be in the order VTK wants
    { 1, 1, -1},
    { 1, -1, 1},
    { 1, -1, -1},
    {-1, 1, 1},
    {-1, 1, -1},
    {-1, -1, 1},
    {-1, -1, -1} };
  this->forParticleSubRegions( [&]( auto & subRegion )
  {
    GEOS_LOG_RANK_0("SubRegion Name: " << subRegion.getName() << ", size " << subRegion.size());
    arrayView2d< real64 const > const particleCenter = subRegion.getParticleCenter();
    GEOS_LOG_RANK_0("Got particle centers");
    if( !subRegion.hasRVectors() ) // if no r-vectors, return coordinates of a cube centered at the particle with side length = volume^(1/3)
    {
      GEOS_LOG_RANK_0("Failed getting volume");
      arrayView1d< real64 const > const particleVolume = subRegion.getParticleVolume();
      forAll< serialPolicy >( subRegion.size(), [=, &coords, &index] GEOS_HOST ( localIndex const p )
        {
          real64 a = 0.5 * pow( particleVolume[p], 1.0/3.0 ); // cube half-side-length
          for( int corner=0; corner<8; corner++ )
          {
            for( int i=0; i<3; i++ )
            {
              coords[index][i] = particleCenter[p][i] + signs[corner][i] * a;
            }
            index++;
          }
        } );
    }
    else
    {
      GEOS_LOG_RANK_0("Before Rvectors");
      arrayView3d< real64 const > const particleRVectors = subRegion.getParticleRVectors();
      GEOS_LOG_RANK_0("After rvectors centers");
      GEOS_LOG_RANK_0(particleRVectors.size() << ", " << particleRVectors.size(0) << ", " << particleRVectors.size(1) << ", " << particleRVectors.size(2));
      forAll< serialPolicy >( subRegion.size(), [=, &coords, &index] GEOS_HOST ( localIndex const p )
        {
          for( int corner=0; corner<8; corner++ )
          {
            for( int i=0; i<3; i++ )
            {
              coords[index][i] = particleCenter[p][i] + signs[corner][0] * particleRVectors[p][0][i] + signs[corner][1] * particleRVectors[p][1][i] + signs[corner][2] * particleRVectors[p][2][i];
            }
            index++;
          }
        } );
    }
  } );
  return coords;
}

REGISTER_CATALOG_ENTRY( ObjectManagerBase, ParticleRegion, string const &, Group * const )

} /* namespace geos */
