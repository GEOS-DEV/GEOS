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

#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/MeshBody.hpp"
#include "mesh/mpiCommunications/PartitionBase.hpp"
#include "mesh/mpiCommunications/SpatialPartition.hpp"


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

  GEOSX_LOG_RANK_0( "MPM particle file name: " << m_filePath );

  int numParticles = 0;

  // Get MPI rank
  int const mpiRank = MpiWrapper::commRank( MPI_COMM_GEOSX );

  if(mpiRank==0) // Only rank 0 should read the particle file (for now)
  {
    std::vector<std::vector<double> > data;
    std::ifstream file(m_filePath);
    std::string line;
    while(std::getline(file, line))
    {
      std::vector<double> lineData;
      std::stringstream lineStream(line);

      double value;
      while(lineStream >> value)
      {
        std::cout << value << "\t"; // debug
        lineData.push_back(value);
      }
      std::cout << "\n"; //debug
      data.push_back(lineData);
      numParticles++;
    }
  }

  // Rank 0 will have to broadcast the particle information to all ranks at some point, probably here. Let's skip this for now and get serial MPM working.

  particleManager.resize( numParticles );
  // Do we want a new namespace for particles so we can do particles::REFERENCE_POSITION_USD and particles::STRESS_PERMUTATION?
  // Probably too redundant, but calling things called "nodes" when operating on particles may be confusing.
  // Good documentation may help this.
  arrayView2d< real64, nodes::REFERENCE_POSITION_USD > const & X = particleManager.referencePosition();


  GEOSX_LOG_RANK_0( "Total number of particles:" << particleManager.size() );

}

REGISTER_CATALOG_ENTRY( MeshGeneratorBase, ParticleMeshGenerator, string const &, Group * const )
} /* namespace geosx */
