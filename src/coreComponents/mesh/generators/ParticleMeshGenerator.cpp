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
  registerWrapper( viewKeyStruct::particleFilePathString(), &m_particleFilePath ).
    setInputFlag( InputFlags::REQUIRED ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "path to the particle file" );

  registerWrapper( viewKeyStruct::headerFilePathString(), &m_headerFilePath ).
    setInputFlag( InputFlags::REQUIRED ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "path to the header file" );

  registerWrapper( viewKeyStruct::particleBlockNamesString(), &m_blockNames ).
    setInputFlag( InputFlags::REQUIRED ).
    setSizedFromParent( 0 ).
    setDescription( "Names of each particle block" );

  registerWrapper( viewKeyStruct::particleTypesString(), &m_particleType ).
    setInputFlag( InputFlags::REQUIRED ).
    setSizedFromParent( 0 ).
    setDescription( "Particle types of each particle block" );
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

  ParticleBlockManager & particleBlockManager = meshBody.registerGroup< ParticleBlockManager >( keys::particleManager );

  SpatialPartition & partition = dynamic_cast< SpatialPartition & >(domain.getReference< PartitionBase >( keys::partitionManager ) );

  // This should probably handled elsewhere:
  int aa = 0;
  for( auto & particleBlockName : m_blockNames )
  {
    ParticleBlock & particleBlock = particleBlockManager.registerParticleBlock( particleBlockName );
    particleBlock.setParticleType( EnumStrings< ParticleType >::fromString( m_particleType[aa++] ) );
  }

  GEOSX_LOG_RANK_0( "MPM particle file path: " << m_particleFilePath );
  GEOSX_LOG_RANK_0( "MPM header file path: " << m_headerFilePath );

  int numParticles = 0;
  int numMaterials, numParticleTypes;
  map < std::string, std::vector<std::vector<double>> > particleData;
  map <std::string, int> particleTypeMap;
  std::vector< std::string > particleTypes; // This is needed because the input file format is such that data associated with each particle type is in the same order as the preceding type listing.
                                            // Looping over the particleTypeMap with an iterator does not respect this ordering since a map automatically sorts by its keys.
  map <std::string, int> materialMap;

  // Get MPI rank
  int const mpiRank = MpiWrapper::commRank( MPI_COMM_GEOSX );

  if(mpiRank==0) // Only rank 0 should read the particle and header file (for now).
  {
    // Get and process header and particle files
    std::ifstream headerFile(m_headerFilePath);
    std::ifstream particleFile(m_particleFilePath);
    std::string line; // initialize line variable

    // Read in number of materials and particle types
    std::getline(headerFile, line); // get a line
    std::istringstream iss1(line); // turn the line into a stream
    iss1 >> numMaterials >> numParticleTypes;
    particleTypes.resize(numParticleTypes);

    std::cout << "Number of particle materials: " << numMaterials << std::endl;
    std::cout << "Number of particle types: " << numParticleTypes << std::endl;

    // Read in material key
    for(int i=0; i<numMaterials; i++)
    {
      std::getline(headerFile, line);
      std::istringstream iss2(line);
      std::string key; // Material name
      int value; // Material ID
      iss2 >> key >> value;
      materialMap[key] = value;
      std::cout << "Material name/ID: " + key + "/" << value << std::endl;
    }

    // Read in particle type key
    for(int i=0; i<numParticleTypes; i++)
    {
      std::getline(headerFile, line);
      std::istringstream iss2(line);
      std::string particleType; // Particle type
      int np; // Number of particles of that type
      iss2 >> particleType >> np;
      numParticles += np;
      particleTypeMap[particleType] = np;
      particleTypes[i] = particleType;
    }

    // Read in particle data
    for(size_t i=0; i<particleTypes.size(); i++)
    {
      for(int j=0; j<particleTypeMap[particleTypes[i]]; j++)
      {
        std::getline(particleFile, line);
        std::vector<double> lineData;
        std::istringstream lineStream(line);

        double value;
        int positionComponent = 0;
        bool inPartition = true;
        while(lineStream >> value)
        {
          lineData.push_back(value);
          if(positionComponent<3) // first three values are the particle position, check for partition membership
          {
            inPartition = inPartition && partition.isCoordInPartition(value, positionComponent);
            if(!inPartition) // if the current particle is outside this partition, we can ignore the rest of its data and go to the next line
            {
              break;
            }
          }
          positionComponent++;
        }
        if(inPartition)
        {
          particleData[particleTypes[i]].push_back(lineData);
        }
      }
    }

    // Print out particle data
//    for(size_t i=0; i<particleTypes.size(); i++)
//    {
//      std::string particleType = particleTypes[i];
//      std::cout << "Printing out data for all particles of type: " << particleType << std::endl;
//
//      std::vector<std::vector<double>> & temp = particleData[particleType];
//      for(size_t j=0; j<temp.size(); j++)
//      {
//        for(size_t k=0; k<temp[j].size(); k++)
//        {
//          std::cout << temp[j][k] << "\t";
//        }
//        std::cout << std::endl;
//      }
//    }
  }

  // Rank 0 will have to broadcast the particle information to all ranks at some point, probably here. Let's skip this for now and get serial MPM working.
  // All ranks should, for now, read in the entire particle file and throw out what they don't need. Ultimately, probably better to have separate files for each rank.
  // Each rank will have its own version of particleData containing only the particles in the partition associated with that rank.

  particleManager.resize( numParticles ); // All this does is change m_size for the particleManager, gives a convenient way to get the total number of particles

  // Get map from particle blocks to particle regions (more specifically the regions' associated materials)
  map < std::string, int > blockMaterialMap;
  particleManager.forParticleRegions< ParticleRegion >( [&]( auto & particleRegion )
  {
    string_array particleBlockNames = particleRegion.getParticleBlockNames();
    std::string material = particleRegion.getMaterialList()[0]; // We will assume that the material list for a region contains only one material since MPM will only be doing single phase mechanics for now
    for(auto i=0; i<particleBlockNames.size(); i++)
    {
      blockMaterialMap[particleBlockNames[i]] = materialMap[material];
    }
  } );

  // Distribute particle information to particle blocks
  map < std::string, std::vector<int> > indexMap; // This will keep track of the indices of particleData associated with each particle block. It's built in the loop that checks for which particles belong to a block.
  map < std::string, int > sizeMap; // This keeps track of the size of each particle block so we can resize the ParticleRegions later
  for( auto & particleBlockName : m_blockNames )
  {
    ParticleBlock & particleBlock = particleBlockManager.getParticleBlock( particleBlockName );
    std::string s = particleBlock.hasRVectors() ? " does " : " doesn't ";
    std::string particleType = EnumStrings< ParticleType >::toString( particleBlock.getParticleType() );
    GEOSX_LOG_RANK( "Particle block " << particleBlock.getName() << " is of type " << particleType << " and hence" << s << "have r-vectors!");

    int numThisType = particleData[particleType].size(); // I hope this returns zero when particleType isn't in the map...
    int materialID;
    int np = 0; // Number of particles in this particle block
    for(localIndex i=0; i<numThisType; i++) // Find out which particles belong to the current particle block
    {
      materialID = particleData[particleType][i][6]; // The particle file is configured such that the 7th column has the material ID
      if(materialID == blockMaterialMap[particleBlock.getName()])
      {
        np++;
        indexMap[particleBlockName].push_back(i);
      }
    }
    sizeMap[particleBlock.getName()] = np;
    particleBlock.resize(np);
    //std::cout << "Particle block " << particleBlock.getName() << " contains " << particleBlock.size() << " particles." << std::endl;

    array2d< real64 > particleCenter(np,3);
    array2d< real64 > particleVelocity(np,3);
    array1d< real64 > particleVolume(np);
    array3d< real64 > particleRVectors(np,3,3);

    // Assign particle data to the appropriate block.
    std::vector<int> & indices = indexMap[particleBlockName];
    int index = 0;
    for(int i : indices)
    {
      // Position
      particleCenter[index][0] = particleData[particleType][i][0];
      particleCenter[index][1] = particleData[particleType][i][1];
      particleCenter[index][2] = particleData[particleType][i][2];

       // Velocity
      particleVelocity[index][0] = particleData[particleType][i][3];
      particleVelocity[index][1] = particleData[particleType][i][4];
      particleVelocity[index][2] = particleData[particleType][i][5];

      // Volume and R-Vectors
      if(particleType == "SinglePoint") // I'm sure there's a better way to handle this case switching
      {
        particleVolume[index] = particleData[particleType][i][7];
        double a = std::pow(particleVolume[index],1.0/3.0);
        particleRVectors[index][0][0] = a;
        particleRVectors[index][0][1] = 0.0;
        particleRVectors[index][0][2] = 0.0;
        particleRVectors[index][1][0] = 0.0;
        particleRVectors[index][1][1] = a;
        particleRVectors[index][1][2] = 0.0;
        particleRVectors[index][2][0] = 0.0;
        particleRVectors[index][2][1] = 0.0;
        particleRVectors[index][2][2] = a;
      }
      else if(particleType == "CPDI")
      {
        double x1, y1, z1, x2, y2, z2, x3, y3, z3;
        x1 = particleData[particleType][i][7];
        y1 = particleData[particleType][i][8];
        z1 = particleData[particleType][i][9];
        x2 = particleData[particleType][i][10];
        y2 = particleData[particleType][i][11];
        z2 = particleData[particleType][i][12];
        x3 = particleData[particleType][i][13];
        y3 = particleData[particleType][i][14];
        z3 = particleData[particleType][i][15];
        particleRVectors[index][0][0] = x1;
        particleRVectors[index][0][1] = y1;
        particleRVectors[index][0][2] = z1;
        particleRVectors[index][1][0] = x2;
        particleRVectors[index][1][1] = y2;
        particleRVectors[index][1][2] = z2;
        particleRVectors[index][2][0] = x3;
        particleRVectors[index][2][1] = y3;
        particleRVectors[index][2][2] = z3;
        particleVolume[index] = std::fabs(-(x3*y2*z1) + x2*y3*z1 + x3*y1*z2 - x1*y3*z2 - x2*y1*z3 + x1*y2*z3);
      }
      else
      {
        GEOSX_ERROR("Invalid particle type specification! Cannot determine particle volume, aborting.");
      }

      // Increment index
      index++;
    }
    particleBlock.setParticleCenter(particleCenter);
    particleBlock.setParticleVelocity(particleVelocity);
    particleBlock.setParticleVolume(particleVolume);
    particleBlock.setParticleRVectors(particleRVectors);
  } // loop over particle blocks

  // Resize particle regions
  int numParticlesCheck = 0;
  particleManager.forParticleRegions< ParticleRegion >( [&]( auto & particleRegion )
  {
    string_array particleBlockNames = particleRegion.getParticleBlockNames();
    std::string material = particleRegion.getMaterialList()[0]; // We will assume that the material list for a region contains only one material since MPM will only be doing single phase mechanics for now
    int size = 0;
    for(auto i=0; i<particleBlockNames.size(); i++)
    {
      size += sizeMap[particleBlockNames[i]];
    }
    numParticlesCheck += size;
    particleRegion.resize(size);
    GEOSX_LOG_RANK("Particle region " << particleRegion.getName() << " contains " << size << " particles on this rank.");
  } );


  //GEOSX_ERROR_IF(numParticlesCheck != numParticles, "Inconsistency detected between MPM particle file and GEOSX input file! Check if you've correctly allocated all particle blocks.");

  GEOSX_LOG_RANK( "Total number of particles on this rank: " << particleManager.size() );
}

void ParticleMeshGenerator::postProcessInput()
{
  GEOSX_LOG_RANK_0("Someone called ParticleMeshGenerator::postProcessInput!");
}

REGISTER_CATALOG_ENTRY( MeshGeneratorBase, ParticleMeshGenerator, string const &, Group * const )
} /* namespace geosx */
