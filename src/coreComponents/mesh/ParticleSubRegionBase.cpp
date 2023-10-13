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
 * @file ParticleSubRegionBase.cpp
 */

#include "ParticleSubRegionBase.hpp"
#include "constitutive/ConstitutiveManager.hpp"

namespace geos
{

using namespace dataRepository;
using namespace constitutive;

ParticleSubRegionBase::ParticleSubRegionBase( string const & name, Group * const parent ):
  ObjectManagerBase( name, parent ),
  m_constitutiveModels( groupKeyStruct::constitutiveModelsString(), this ),
  m_hasRVectors(),
  m_particleRank(),
  m_particleID(),
  m_particleGroup(),
  m_particleSurfaceFlag(),
  m_particleDamage(),
  m_particleStrengthScale(),
  m_particleCenter(),
  m_particleVelocity(),
  m_particleMaterialDirection(),
  m_particleVolume(),
  m_particleType(),
  m_particleRVectors()
{
  registerGroup( groupKeyStruct::constitutiveModelsString(), &m_constitutiveModels ).
    setSizedFromParent( 1 );

  registerWrapper( viewKeyStruct::particleRankString(), &m_particleRank ).
    setPlotLevel( PlotLevel::LEVEL_1 );

  registerWrapper( viewKeyStruct::particleIDString(), &m_particleID ).
    setPlotLevel( PlotLevel::LEVEL_1 );

  registerWrapper( viewKeyStruct::particleGroupString(), &m_particleGroup ).
    setPlotLevel( PlotLevel::LEVEL_1 );

  registerWrapper( viewKeyStruct::particleSurfaceFlagString(), &m_particleSurfaceFlag ).
    setPlotLevel( PlotLevel::LEVEL_1 );

  registerWrapper( viewKeyStruct::particleDamageString(), &m_particleDamage ).
    setPlotLevel( PlotLevel::LEVEL_1 );

  registerWrapper( viewKeyStruct::particleStrengthScaleString(), &m_particleStrengthScale ).
    setPlotLevel( PlotLevel::LEVEL_1 );

  registerWrapper( viewKeyStruct::particleCenterString(), &m_particleCenter ).
    setPlotLevel( PlotLevel::NOPLOT ).
    reference().resizeDimension< 1 >( 3 );

  registerWrapper( viewKeyStruct::particleVelocityString(), &m_particleVelocity ).
    setPlotLevel( PlotLevel::LEVEL_1 ).
    reference().resizeDimension< 1 >( 3 );

  registerWrapper( viewKeyStruct::particleMaterialDirectionString(), &m_particleMaterialDirection ).
    setPlotLevel( PlotLevel::LEVEL_1 ).
    reference().resizeDimension< 1 >( 3 );

  registerWrapper( viewKeyStruct::particleVolumeString(), &m_particleVolume ).
    setPlotLevel( PlotLevel::LEVEL_1 );

  registerWrapper( viewKeyStruct::particleRVectorsString(), &m_particleRVectors ).
    setPlotLevel( PlotLevel::NOPLOT ).
    reference().resizeDimension< 1, 2 >( 3, 3 );
}

ParticleSubRegionBase::~ParticleSubRegionBase()
{}

unsigned int ParticleSubRegionBase::particlePack( buffer_type & buffer,
                                                  arrayView1d< localIndex > const & localIndices,
                                                  bool doPack ) const
{
  // Declarations
  parallelDeviceEvents events; // I have no idea what this thing is
  unsigned int packedSize = 0;

  // Pack particle fields
  if( !doPack ) // doPack == false, so we're just getting the size
  {
    packedSize += this->packSize( localIndices, 0, false, events );
  }
  else // doPack == true, perform the pack
  {
    buffer_unit_type * bufferPtr = buffer.data();
    packedSize += this->pack( bufferPtr, localIndices, 0, false, events );
  }

  return packedSize;
}

void ParticleSubRegionBase::particleUnpack( buffer_type & buffer,
                                            int const & startingIndex,
                                            int const & numberOfIncomingParticles )
{
  // Declarations
  parallelDeviceEvents events; // I have no idea what this thing is
  const buffer_unit_type * receiveBufferPtr = buffer.data(); // needed for const cast

  // Get the indices we're overwriting during unpack.
  // Before unpacking, those indices should contain junk/default data created when subRegion.resize() was called.
  array1d< localIndex > indices( numberOfIncomingParticles );
  for( int i=0; i<numberOfIncomingParticles; i++ )
  {
    indices[i] = startingIndex + i;
  }

  // Unpack
  this->unpack( receiveBufferPtr, indices, 0, false, events );
}

void ParticleSubRegionBase::erase( std::set< localIndex > const & indicesToErase )
{
  if( indicesToErase.size() > 0 )
  {
    // The new subregion size
    int newSize = (this->size()) - indicesToErase.size();

    // Call ObjectManagerBase::eraseObject
    this->eraseObject( indicesToErase );

    // Decrease the size of this subregion
    this->resize( newSize );

    // Reconstruct the list of non-ghost indices
    this->setActiveParticleIndices();
  }
}

void ParticleSubRegionBase::setActiveParticleIndices()
{
  m_activeParticleIndices.move( LvArray::MemorySpace::host ); // TODO: Is this needed?
  m_activeParticleIndices.clear();
  arrayView1d< int const > const particleRank = m_particleRank.toViewConst();
  forAll< serialPolicy >( this->size(), [&, particleRank] GEOS_HOST ( localIndex const p ) // This must be on host since we're dealing with
                                                                                           // a sorted array. Parallelize with atomics?
  {
    if( particleRank[p] == MpiWrapper::commRank( MPI_COMM_GEOSX ) )
    {
      m_activeParticleIndices.insert( p );
    }
  } );
}

void ParticleSubRegionBase::updateMaps()
{
  arrayView1d< globalIndex > const localToGlobalMap = m_localToGlobalMap;
  arrayView1d< globalIndex const > const particleID = m_particleID;
  forAll< serialPolicy >( this->size(), [=] GEOS_HOST ( localIndex const p ) // TODO: must be on host because constructGlobalToLocalMap has
                                                                             // to be on host?
    {
      localToGlobalMap[p] = particleID[p];
    } );
  this->constructGlobalToLocalMap();
}

} /* namespace geos */
