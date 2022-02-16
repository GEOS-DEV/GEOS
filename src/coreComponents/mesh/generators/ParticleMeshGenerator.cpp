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
 * @file ParticleMeshGenerator.cpp
 */

#include "ParticleMeshGenerator.hpp"

#include "mesh/DomainPartition.hpp"
#include "mesh/mpiCommunications/PartitionBase.hpp"
#include "mesh/mpiCommunications/SpatialPartition.hpp"
#include "mesh/MeshBody.hpp"
#include "ParticleBlockManager.hpp"

#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"

#include <cmath>

namespace geosx
{
using namespace dataRepository;

ParticleMeshGenerator::ParticleMeshGenerator( string const & name, Group * const parent ):
  MeshGeneratorBase( name, parent ),
  m_dim( 3 ),
  m_min(),
  m_max(),
  m_pCoord{}
{
  registerWrapper( viewKeyStruct::filePathString(), &m_filePath ).
    setInputFlag( InputFlags::REQUIRED ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "path to the particle file" );

  registerWrapper( viewKeyStruct::particleBlockNamesString(), &m_regionNames ).
    setInputFlag( InputFlags::REQUIRED ).
    setSizedFromParent( 0 ).
    setDescription( "Names of each particle block" );

  registerWrapper( viewKeyStruct::particleTypesString(), &m_particleType ).
    setInputFlag( InputFlags::REQUIRED ).
    setSizedFromParent( 0 ).
    setDescription( "Particle types of each mesh block" );
}

Group * ParticleMeshGenerator::createChild( string const & GEOSX_UNUSED_PARAM( childKey ),
                                            string const & GEOSX_UNUSED_PARAM( childName ) )
{
  return nullptr;
}

/**
 * @param partition
 * @param domain
 */
void ParticleMeshGenerator::generateMesh( DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  MeshBody & meshBody = domain.getMeshBody( this->getName() );
  GEOSX_LOG_RANK( "Generating mesh for " + meshBody.getName() );
  MeshLevel & meshLevel0 = meshBody.getMeshLevel( 0 );
  ParticleManager & particleManager = meshLevel0.getParticleManager();

  ParticleBlockManager & particleBlockManager = meshBody.getGroup< ParticleBlockManager >( keys::particleManager );

  SpatialPartition & partition = dynamic_cast< SpatialPartition & >(domain.getReference< PartitionBase >( keys::partitionManager ) );

  // This should probably handled elsewhere:
  int aa = 0;
  for( auto & particleBlockName : m_regionNames )
  {
    ParticleBlock & particleBlock = particleBlockManager.registerParticleBlock( particleBlockName );
    particleBlock.setParticleType( EnumStrings< ParticleType >::fromString( m_particleType[aa++] ) );
  }

  GEOSX_LOG_RANK_0( "MPM particle file path: " << m_filePath );

  int numParticles = 0;

  // Get MPI rank
  int const mpiRank = MpiWrapper::commRank( MPI_COMM_GEOSX );


  std::vector<std::vector<double> > particleData;
  if(mpiRank==0) // Only rank 0 should read the particle file (for now).
  {
    std::ifstream file(m_filePath);
    std::string line;
    while(std::getline(file, line))
    {
      std::vector<double> lineData;
      std::stringstream lineStream(line);

      double value;
      while(lineStream >> value)
      {
        //std::cout << value << "\t"; // debug
        lineData.push_back(value);
      }
      //std::cout << "\n"; //debug
      particleData.push_back(lineData);
      numParticles++;
    }
  }

  // Rank 0 will have to broadcast the particle information to all ranks at some point, probably here. Let's skip this for now and get serial MPM working.
  // All ranks should, for now, read in the entire particle file and throw out what they don't need. Ultimately, probably better to have separate files for each rank.

  particleManager.resize( numParticles );
  arrayView2d< real64, particles::REFERENCE_POSITION_USD > const & X = particleManager.referencePosition();

  // Assign the read-in particle coordinates to the particle reference positions.
  for(localIndex i=0; i<numParticles; i++)
  {
    for(localIndex j=0; j<3; j++)
    {
      X[i][j] = particleData[i][j];
      std::cout << X[i][j] << "\t";
    }
    std::cout << "\n"; //debug
  }

  GEOSX_LOG_RANK_0( "Total number of particles:" << particleManager.size() );

}

void ParticleMeshGenerator::postProcessInput()
{
  GEOSX_LOG_RANK_0("Someone called ParticleMeshGenerator::postProcessInput!");
}

REGISTER_CATALOG_ENTRY( MeshGeneratorBase, ParticleMeshGenerator, string const &, Group * const )
} /* namespace geosx */
