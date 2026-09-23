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
 * @file ParticleRegion.cpp
 */

#include "ParticleRegion.hpp"
#include "common/TimingMacros.hpp"
#include "LvArray/src/SparsityPattern.hpp"
#include "common/MpiWrapper.hpp"
#include "physicsSolvers/solidMechanics/MPMSolverFields.hpp"
#include <cmath>

namespace geos
{
using namespace dataRepository;

ParticleRegion::ParticleRegion( string const & name, Group * const parent ):
  ParticleRegionBase( name, parent )
{
  registerWrapper( viewKeyStruct::sourceParticleBlockNamesString(), &m_particleBlockNames ).
    setInputFlag( InputFlags::OPTIONAL );
  registerWrapper( "particleDomainOutput", &m_particleDomainOutput ).
    setInputFlag( InputFlags::OPTIONAL ).setApplyDefaultValue( "native" ).
    setDescription( "native: use optional simplex reference geometry; equivalentBox: legacy particle domains" );
}

ParticleRegion::~ParticleRegion()
{}

void ParticleRegion::generateMesh( Group & particleBlocks )
{
  GEOS_ERROR_IF( m_particleDomainOutput != "native" && m_particleDomainOutput != "equivalentBox",
                 "particleDomainOutput must be native or equivalentBox" );
  Group & particleSubRegions = this->getGroup( viewKeyStruct::particleSubRegions() );

  for( string const & particleBlockName : this->m_particleBlockNames )
  {
    ParticleSubRegion & subRegion = particleSubRegions.registerGroup< ParticleSubRegion >( particleBlockName );
    ParticleBlockABC & source = particleBlocks.getGroup< ParticleBlockABC >( subRegion.getName() );

    // Copy data from particle blocks into subregions. This also sets particle type.
    subRegion.copyFromParticleBlock( source );

    // Set the rank of particles on each subregion
    int const mpiRank = MpiWrapper::commRank( MPI_COMM_GEOS );
    subRegion.setParticleRank( mpiRank );

    // Set the number of vertices of the particles on each subregion
    subRegion.setNumberOfVerticesPerParticle( subRegion.getParticleType() );
  }
}

// TODO This should be changed to call a ParticleSubRegion::getParticleCoordinates (and/or getParticleCorners) on each subregion such that
// we can access those functions directly if needed
array2d< real64 > ParticleRegion::getParticleCorners() const
{
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
    arrayView2d< real64 const > const particleCenter = subRegion.getParticleCenter();
    if( !subRegion.hasRVectors() ) // if no r-vectors, return coordinates of a cube centered at the particle with side length = volume^(1/3)
    {
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
      arrayView3d< real64 const > const particleRVectors = subRegion.getParticleRVectors();
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

bool ParticleRegion::usesNativeParticleGeometry() const
{
  if( m_particleDomainOutput == "equivalentBox" ) return false;
  bool enabled = false;
  forParticleSubRegions( [&]( auto const & subRegion )
  { enabled = enabled || subRegion.hasParticleReferenceGeometry(); } );
  return enabled;
}

ParticleRegion::ParticleDomainMesh ParticleRegion::getParticleDomainMesh() const
{
  ParticleDomainMesh mesh;
  mesh.offsets.push_back( 0 );
  forParticleSubRegions( [&]( auto const & subRegion )
  {
    // Output is host-side. This is synchronization, not a mechanical mutation.
    // WrapperBase::move updates the allocation's memory residency and timestamp.
    auto & hostSubRegion = const_cast< ParticleSubRegion & >( subRegion );
    auto move = [&]( char const * key )
    { hostSubRegion.getWrapperBase( key ).move( LvArray::MemorySpace::host, false ); };
    move( ParticleSubRegion::viewKeyStruct::particleCenterString() );
    move( ParticleSubRegion::viewKeyStruct::particleRankString() );
    move( ParticleSubRegion::viewKeyStruct::particleVolumeString() );
    move( ParticleSubRegion::viewKeyStruct::particleRVectorsString() );
    bool const enabled = subRegion.hasParticleReferenceGeometry();
    bool const hasF = subRegion.hasWrapper( fields::mpm::particleDeformationGradient::key() );
    if( enabled ) { move( particleGeometry::shapeKey() ); move( particleGeometry::vectorsKey() ); }
    if( hasF ) move( fields::mpm::particleDeformationGradient::key() );
    auto const center = subRegion.getParticleCenter();
    auto const rank = subRegion.getParticleRank();
    auto const volume = subRegion.getParticleVolume();
    auto const r = subRegion.getParticleRVectors();
    for( localIndex p=0; p<subRegion.size(); ++p )
    {
      integer shape = 0;
      if( enabled )
      {
        auto const referenceShape = subRegion.template getReference< array1d< integer > >( particleGeometry::shapeKey() );
        shape = referenceShape[p];
      }
      GEOS_ERROR_IF( particleGeometry::vertexCount( shape ) == 0, "Invalid particle reference geometry in output" );
      real64 g[3][3] = {};
      if( shape == 0 )
      {
        // Box particles continue to expose the actual legacy numerical domain.
        if( subRegion.hasRVectors() )
          for( int a=0; a<3; ++a ) for( int j=0; j<3; ++j ) g[a][j] = r[p][a][j];
        else for( int a=0; a<3; ++a ) g[a][a] = 0.5*std::cbrt( volume[p] );
      }
      else
      {
        auto const & reference = subRegion.template getReference< array3d< real64 > >( particleGeometry::vectorsKey() );
        if( hasF )
        {
          auto const deformationGradient = subRegion.template getField< fields::mpm::particleDeformationGradient >();
          particleGeometry::deform( reference[p], deformationGradient[p], g );
        }
        else for( int a=0; a<3; ++a ) for( int j=0; j<3; ++j ) g[a][j] = reference[p][a][j];
      }
      real64 v[8][3] = {};
      particleGeometry::vertices( shape, center[p], g, v );
      for( int k=0; k<particleGeometry::vertexCount( shape ); ++k )
        mesh.points.push_back( {{ v[k][0], v[k][1], v[k][2] }} );
      mesh.shape.push_back( shape );
      mesh.ownerRank.push_back( rank[p] );
      mesh.offsets.push_back( static_cast< localIndex >( mesh.points.size() ) );
    }
  } );
  return mesh;
}

REGISTER_CATALOG_ENTRY( ObjectManagerBase, ParticleRegion, string const &, Group * const )

} /* namespace geos */
