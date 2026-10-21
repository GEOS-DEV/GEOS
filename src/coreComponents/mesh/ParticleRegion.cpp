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
    int const mpiRank = MpiWrapper::commRank( MPI_COMM_GEOS );
    subRegion.setParticleRank( mpiRank );

  }
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
    // Output needs current material geometry only; avoid staging G0 or F.
    move( particleGeometry::shapeKey() );
    move( particleGeometry::currentKey() );
    auto const shapes = subRegion.getParticleGeometryShape();
    auto const generators = subRegion.getParticleGeometryCurrent();
    auto const center = subRegion.getParticleCenter();
    auto const rank = subRegion.getParticleRank();
    for( localIndex p=0; p<subRegion.size(); ++p )
    {
      auto const domain = particleGeometry::affine::Domain::fromGenerators(
        static_cast<particleGeometry::Shape>(shapes[p]), generators[p] );
      integer const shape = static_cast< integer >( domain.shape );
      GEOS_ERROR_IF( particleGeometry::vertexCount( shape ) == 0, "Invalid particle reference geometry in output" );
      real64 v[8][3] = {};
      domain.vertices( center[p], v );
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
