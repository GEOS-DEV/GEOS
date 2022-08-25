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

namespace geosx
{

using namespace dataRepository;

ParticleSubRegionBase::ParticleSubRegionBase( string const & name, Group * const parent ):
  ObjectManagerBase( name, parent ),
  m_constitutiveModels( groupKeyStruct::constitutiveModelsString(), this ),
  m_hasRVectors(),
  m_particleGhostRank(),
  m_particleID(),
  m_particleGroup(),
  m_particleCenter(),
  m_particleVelocity(),
  m_particleVolume(),
  m_particleInitialVolume(),
  m_particleMass(),
  m_particleDeformationGradient(),
  m_particleType(),
  m_particleRVectors(),
  m_particleInitialRVectors()
{
  registerGroup( groupKeyStruct::constitutiveModelsString(), &m_constitutiveModels ).
    setSizedFromParent( 1 );

  registerWrapper( viewKeyStruct::particleGhostRankString(), &m_particleGhostRank ).
    setPlotLevel( PlotLevel::LEVEL_1 );

  registerWrapper( viewKeyStruct::particleIDString(), &m_particleID ).
    setPlotLevel( PlotLevel::LEVEL_1 );

  registerWrapper( viewKeyStruct::particleGroupString(), &m_particleGroup ).
    setPlotLevel( PlotLevel::LEVEL_1 );

  registerWrapper( viewKeyStruct::particleCenterString(), &m_particleCenter ).
    setPlotLevel( PlotLevel::LEVEL_1 ).
    reference().resizeDimension< 1 >( 3 );

  registerWrapper( viewKeyStruct::particleVelocityString(), &m_particleVelocity ).
    setPlotLevel( PlotLevel::LEVEL_1 ).
    reference().resizeDimension< 1 >( 3 );

  registerWrapper( viewKeyStruct::particleVolumeString(), &m_particleVolume ).
    setPlotLevel( PlotLevel::LEVEL_1 );

  registerWrapper( viewKeyStruct::particleInitialVolumeString(), &m_particleInitialVolume ).
    setPlotLevel( PlotLevel::LEVEL_1 );

  registerWrapper( viewKeyStruct::particleMassString(), &m_particleMass ).
    setPlotLevel( PlotLevel::LEVEL_1 );

  registerWrapper( viewKeyStruct::particleRVectorsString(), &m_particleRVectors ).
    setPlotLevel( PlotLevel::NOPLOT ).
    reference().resizeDimension< 1, 2 >( 3, 3 );

  registerWrapper( viewKeyStruct::particleInitialRVectorsString(), &m_particleInitialRVectors ).
    setPlotLevel( PlotLevel::NOPLOT ).
    reference().resizeDimension< 1, 2 >( 3, 3 );

  // The only things that should be registered here are those that are read in from the input files. So e.g. particle mass shouldn't be here since it's not specified in any input file, whereas particle volume is.
  // A solver(?) should then on the first cycle register particle mass and calculate it from volume and density. Same idea for other similarly post-processed and/or solver-specific fields.

  registerWrapper( viewKeyStruct::particleDeformationGradientString(), &m_particleDeformationGradient ).
    setPlotLevel( PlotLevel::LEVEL_1 ).
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
  if(!doPack) // doPack == false, so we're just getting the size
  {
    packedSize += this->packSize( localIndices, 0, false, events );
  }
  else // doPack == true, perform the pack
  {
    buffer_unit_type* bufferPtr = buffer.data();
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
  const buffer_unit_type* receiveBufferPtr = buffer.data(); // needed for const cast

  // Get the indices we're overwriting during unpack.
  // Before unpacking, those indices should contain junk/default data created when subRegion.resize() was called.
  array1d< localIndex > indices(numberOfIncomingParticles);
  for(int i=0; i<numberOfIncomingParticles; i++)
  {
    indices[i] = startingIndex + i;
  }

  // Unpack
  this->unpack( receiveBufferPtr, indices, 0, false, events );
}

void ParticleSubRegionBase::erase(localIndex pp)
{
  int newSize = this->size()-1;

  // Scalar fields:
  m_particleGhostRank.erase(pp); // TODO: Can we automatically loop over all registered wrappers and erase that way?
  m_particleID.erase(pp);
  m_particleVolume.erase(pp);
  m_particleInitialVolume.erase(pp);
  m_particleMass.erase(pp);

  // Vector fields:
  this->eraseVector( m_particleCenter, pp );
  this->eraseVector( m_particleVelocity, pp );

  // Matrix fields:
  this->eraseTensor( m_particleDeformationGradient, pp );
  this->eraseTensor( m_particleRVectors, pp );
  this->eraseTensor( m_particleInitialRVectors, pp );

  // Decrement the size of this subregion
  this->resize(newSize);
}

void ParticleSubRegionBase::eraseVector(array2d< real64 > & vector, localIndex index)
{
  int oldSize = this->size();
  int newSize = this->size()-1;

  array1d< real64 > temp(3*oldSize);
  for(int i=0; i<3*oldSize; i++)
  {
    temp[i] = vector[i/3][i%3];
  }
  temp.erase(3*index+2);
  temp.erase(3*index+1);
  temp.erase(3*index+0);
  vector.resize(newSize,3);
  for(int i=0; i<3*newSize; i++) // TODO: This can maybe be optimized to start from 'index' since everything before 'index' should be unchanged.
  {                              //       Depends on whether 'resize' preserves existing entries.
    vector[i/3][i%3] = temp[i];
  }
}

void ParticleSubRegionBase::eraseTensor(array3d< real64 > & tensor, localIndex index)
{
  int oldSize = this->size();
  int newSize = this->size()-1;

  array1d< real64 > temp(9*oldSize);
  for(int i=0; i<9*oldSize; i++)
  {
    temp[i] = tensor[i/9][i/3%3][i%3];
  }
  temp.erase(9*index+8);
  temp.erase(9*index+7);
  temp.erase(9*index+6);
  temp.erase(9*index+5);
  temp.erase(9*index+4);
  temp.erase(9*index+3);
  temp.erase(9*index+2);
  temp.erase(9*index+1);
  temp.erase(9*index+0);
  tensor.resize(newSize,3,3);
  for(int i=0; i<9*newSize; i++)
  {
    tensor[i/9][i/3%3][i%3] = temp[i];
  }
}

} /* namespace geosx */
