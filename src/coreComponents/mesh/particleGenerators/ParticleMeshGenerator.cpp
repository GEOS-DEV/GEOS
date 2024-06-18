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

#include "ParticleBlockManager.hpp"
#include "mesh/mpiCommunications/SpatialPartition.hpp"

#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"

#include <cmath>

namespace geos
{
using namespace dataRepository;

ParticleMeshGenerator::ParticleMeshGenerator( string const & name, Group * const parent ):
  ParticleMeshGeneratorBase( name, parent )
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

  registerWrapper( viewKeyStruct::particleTypesString(), &m_particleTypes ).
    setInputFlag( InputFlags::REQUIRED ).
    setSizedFromParent( 0 ).
    setDescription( "Particle types of each particle block" );
}

Group * ParticleMeshGenerator::createChild( string const & GEOS_UNUSED_PARAM( childKey ),
                                            string const & GEOS_UNUSED_PARAM( childName ) )
{
  return nullptr;
}


void ParticleMeshGenerator::fillParticleBlockManager( ParticleBlockManager & particleBlockManager, ParticleManager & particleManager, SpatialPartition const & partition )
{
  GEOS_MARK_FUNCTION;

  // This should probably handled elsewhere:
  int aa = 0;
  for( auto & particleBlockName : m_blockNames )
  {
    ParticleBlock & particleBlock = particleBlockManager.registerParticleBlock( particleBlockName );
    particleBlock.setParticleType( EnumStrings< ParticleType >::fromString( m_particleTypes[aa++] ) );
  }

  int numMaterials, numParticleTypes;
  map< std::string, std::vector< std::vector< double > > > particleData;
  map< std::string, int > particleTypeMap;
  std::vector< std::string > particleTypes; // This is needed because the input file format is such that data associated with each particle
                                            // type is in the same order as the preceding type listing.
                                            // Looping over the particleTypeMap with an iterator does not respect this ordering since a map
                                            // automatically sorts by its keys.
  map< std::string, int > materialMap;

  // Get and process header and particle files
  std::ifstream headerFile( m_headerFilePath );
  std::ifstream particleFile( m_particleFilePath );
  std::string line; // initialize line variable

  // Read in number of materials and particle types
  // TODO CC remove the header file, all of this should be read from the particle mesh generator xml block
  std::getline( headerFile, line ); // get a line
  std::istringstream iss1( line ); // turn the line into a stream
  iss1 >> numMaterials >> numParticleTypes;
  particleTypes.resize( numParticleTypes );

  // numParticleTypes = m_particleTypes.size();
  // numMaterials = m_materialNames.size(); // As read directly from XML block

  GEOS_LOG_RANK_0( "Number of particle materials: " << numMaterials );
  GEOS_LOG_RANK_0( "Number of particle types: " << numParticleTypes );

  // Read in the material key
  for( int i=0; i<numMaterials; i++ )
  {
    std::getline( headerFile, line );
    std::istringstream iss2( line );
    std::string key; // Material name
    int value; // Material ID
    iss2 >> key >> value;
    materialMap[key] = value;
    // GEOS_LOG_RANK_0( "Material name/ID: " + key + "/" << value );
  }

  // Read in the particle type key
  for( int i=0; i<numParticleTypes; i++ )
  {
    std::getline( headerFile, line );
    std::istringstream iss2( line );
    std::string particleType; // Particle type
    int np; // Number of particles of that type
    iss2 >> particleType >> np;
    particleTypeMap[particleType] = np;
    particleTypes[i] = particleType;
  }

  // Read column headers for particle data
  int numColumnHeaders = 0;
  std::map< int, int > columnHeaderMap;
  {
    std::getline( particleFile, line );
    std::string token;
  
    std::istringstream lineStream( line );
    while( std::getline( lineStream, token, '\t') )
    {
      columnHeaderMap.insert( std::pair< int, int >( static_cast< int >( EnumStrings< ParticleColumnHeaders >::fromString( token ) ), numColumnHeaders ) );
      numColumnHeaders++;
    }
  }

  // Read in particle data
  int lineNumber = 1; // Since colum header takes one line
  for( size_t i=0; i < particleTypes.size(); i++ )
  {
    for( int j=0; j < particleTypeMap[particleTypes[i]]; j++ )
    {
      std::getline( particleFile, line );
      std::istringstream lineStream( line );

      std::vector< double > lineData; // TODO: Not great because we cast all input as doubles, but it all gets re-cast later so maybe it's fine.

      // Read line from particle file and parse columns
      double value;
      int numColumns = 0;
      while( lineStream >> value )
      {
        lineData.push_back( value );
        numColumns++;
      }

      GEOS_ERROR_IF( numColumns != numColumnHeaders, "Particle file line " << lineNumber << " has a different number of terms than the column headers! Was " << numColumns << " but should be " << numColumnHeaders );

      lineNumber++;

      // If particle is inside partition add to particleData otherwise ignore and continue parsing file
      bool inPartition = partition.isCoordInPartition( lineData[ columnHeaderMap[ static_cast< int >( ParticleColumnHeaders::PositionX ) ] ], 0 ) && 
                         partition.isCoordInPartition( lineData[ columnHeaderMap[ static_cast< int >( ParticleColumnHeaders::PositionY ) ] ], 1 ) && 
                         partition.isCoordInPartition( lineData[ columnHeaderMap[ static_cast< int >( ParticleColumnHeaders::PositionZ ) ] ], 2 );
      if( !inPartition )
      {
        continue;
      }

      // Reformat particle data and apply defaults to fields not specified
      std::vector< double > lineDataInside;
      // CC: TODO: Can you get the number of options from the enum directly?
      for(int c = 0; c < EnumSize<ParticleColumnHeaders>; c++)
      {
        if( columnHeaderMap.find( c ) != columnHeaderMap.end() )
        {
          lineDataInside.push_back( lineData[ columnHeaderMap[ c ] ] );
          continue;
        }

        // Apply default value
        double defaultValue;
        switch( static_cast< ParticleColumnHeaders >( c ) )
        {
          case ParticleColumnHeaders::StrengthScale:
          case ParticleColumnHeaders::MaterialDirectionX:
          case ParticleColumnHeaders::SurfaceNormalX:
            defaultValue = 1.0;
            break;
          case ParticleColumnHeaders::Temperature:
            defaultValue = 300.0;
            break;
          case ParticleColumnHeaders::MaterialType:
          case ParticleColumnHeaders::ContactGroup:
          case ParticleColumnHeaders::Damage:
          case ParticleColumnHeaders::Porosity:
          case ParticleColumnHeaders::VelocityX:
          case ParticleColumnHeaders::VelocityY:
          case ParticleColumnHeaders::VelocityZ:
          case ParticleColumnHeaders::MaterialDirectionY:
          case ParticleColumnHeaders::MaterialDirectionZ:
          case ParticleColumnHeaders::SurfaceNormalY:
          case ParticleColumnHeaders::SurfaceNormalZ:
          case ParticleColumnHeaders::SurfacePositionX:
          case ParticleColumnHeaders::SurfacePositionY:
          case ParticleColumnHeaders::SurfacePositionZ:
          case ParticleColumnHeaders::SurfaceTractionX:
          case ParticleColumnHeaders::SurfaceTractionY:
          case ParticleColumnHeaders::SurfaceTractionZ:
            defaultValue = 0.0;
            break;
          default:
            GEOS_ERROR( EnumStrings< ParticleColumnHeaders >::toString( static_cast< ParticleColumnHeaders >( c ) ) << " must be specified in particle file!" );
            break;
        }

        lineDataInside.push_back( defaultValue );
      }

      particleData[particleTypes[i]].push_back( lineDataInside );
    }
  }

  // Construct the map from particle blocks to particle regions (more specifically, the regions' associated materials)
  map< std::string, int > blockMaterialMap;
  particleManager.forParticleRegions< ParticleRegion >( [&]( auto & particleRegion )
  {
    string_array particleBlockNames = particleRegion.getParticleBlockNames();
    std::string material = particleRegion.getMaterialList()[0]; // We will assume that the material list for a region contains only one
                                                                // material since MPM will only be doing single phase mechanics for now
    for( auto i=0; i < particleBlockNames.size(); i++ )
    {
      blockMaterialMap[particleBlockNames[i]] = materialMap[material];
    }
  } );

  // Distribute particle information to particle blocks
  map< std::string, std::vector< int > > indexMap; // This will keep track of the indices of particleData associated with each particle
                                                   // block. It's populated in the loop that checks for which particles belong to a block.
  map< std::string, int > sizeMap;  // This keeps track of the size of each particle block so we can resize the ParticleRegions later
  for( auto & particleBlockName : m_blockNames )
  {
    ParticleBlock & particleBlock = particleBlockManager.getParticleBlock( particleBlockName );
    std::string particleType = EnumStrings< ParticleType >::toString( particleBlock.getParticleType() );

    int numThisType = particleData[particleType].size(); // I hope this returns zero when particleType isn't in the map...
    int materialID;
    int npInBlock = 0; // Number of particles in this particle block
    for( localIndex i=0; i<numThisType; i++ ) // Find out which particles belong to the current particle block
    {
      materialID = particleData[particleType][i][static_cast< int >( ParticleColumnHeaders::MaterialType )];
      if( materialID == blockMaterialMap[particleBlock.getName()] )
      {
        npInBlock++;
        indexMap[particleBlockName].push_back( i );
      }
    }
    sizeMap[particleBlock.getName()] = npInBlock;
    particleBlock.resize( npInBlock );

    array1d< globalIndex > particleID( npInBlock );
    array2d< real64 > particleCenter( npInBlock, 3 );
    array2d< real64 > particleVelocity( npInBlock, 3 );
    array2d< real64 > particleMaterialDirection( npInBlock, 3 );
    array1d< int > particleGroup( npInBlock );
    array1d< int > particleSurfaceFlag( npInBlock );
    array1d< real64 > particleDamage( npInBlock );
    array1d< real64 > particlePorosity( npInBlock );
    array1d< real64 > particleTemperature( npInBlock );
    array1d< real64 > particleVolume( npInBlock );
    array1d< real64 > particleStrengthScale( npInBlock );
    array3d< real64 > particleRVectors( npInBlock, 3, 3 ); // TODO: Flatten the r-vector array into a 1x9 for each particle
    array2d< real64 > particleSurfaceNormal( npInBlock, 3); // TODO:: read from file eventually
    array2d< real64 > particleSurfacePosition( npInBlock, 3 );
    array2d< real64 > particleSurfaceTraction( npInBlock, 3 );

    // Assign particle data to the appropriate block.
    std::vector< int > & indices = indexMap[particleBlockName];
    int index = 0;
    for( int i : indices )
    {
      // Global ID
      particleID[index] = particleData[particleType][i][static_cast< int >( ParticleColumnHeaders::ID )];

      // Position
      particleCenter[index][0] = particleData[particleType][i][static_cast< int >( ParticleColumnHeaders::PositionX )];
      particleCenter[index][1] = particleData[particleType][i][static_cast< int >( ParticleColumnHeaders::PositionY )];
      particleCenter[index][2] = particleData[particleType][i][static_cast< int >( ParticleColumnHeaders::PositionZ )];

      // Velocity
      particleVelocity[index][0] = particleData[particleType][i][static_cast< int >( ParticleColumnHeaders::VelocityX )];
      particleVelocity[index][1] = particleData[particleType][i][static_cast< int >( ParticleColumnHeaders::VelocityY )];
      particleVelocity[index][2] = particleData[particleType][i][static_cast< int >( ParticleColumnHeaders::VelocityZ )];
  
      // Material (set above) is [10]

      // Group
      particleGroup[index] = particleData[particleType][i][static_cast< int >( ParticleColumnHeaders::ContactGroup )];

      // surfaceFlag
      particleSurfaceFlag[index] = particleData[particleType][i][static_cast< int >( ParticleColumnHeaders::SurfaceFlag )];

      // Damage
      particleDamage[index] = particleData[particleType][i][static_cast< int >( ParticleColumnHeaders::Damage )];

      // Porosity
      particlePorosity[index] = particleData[particleType][i][static_cast< int >( ParticleColumnHeaders::Porosity )];

      // Temperature
      particleTemperature[index] = particleData[particleType][i][static_cast< int >( ParticleColumnHeaders::Temperature )];

      // strengthScale
      particleStrengthScale[index] = particleData[particleType][i][static_cast< int >( ParticleColumnHeaders::StrengthScale )];

      // Volume and R-Vectors
      if( particleType == "SinglePoint" )
      {
        particleVolume[index] = particleData[particleType][i][static_cast< int >( ParticleColumnHeaders::RVectorXX )];

        double a = std::pow( particleVolume[index], 1.0/3.0 );
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
      else if( particleType == "CPDI" )
      {
        double x1, y1, z1, x2, y2, z2, x3, y3, z3;
        x1 = particleData[particleType][i][static_cast< int >( ParticleColumnHeaders::RVectorXX )];
        y1 = particleData[particleType][i][static_cast< int >( ParticleColumnHeaders::RVectorXY )];
        z1 = particleData[particleType][i][static_cast< int >( ParticleColumnHeaders::RVectorXZ )];
        x2 = particleData[particleType][i][static_cast< int >( ParticleColumnHeaders::RVectorYX )];
        y2 = particleData[particleType][i][static_cast< int >( ParticleColumnHeaders::RVectorYY )];
        z2 = particleData[particleType][i][static_cast< int >( ParticleColumnHeaders::RVectorYZ )];
        x3 = particleData[particleType][i][static_cast< int >( ParticleColumnHeaders::RVectorZX )];
        y3 = particleData[particleType][i][static_cast< int >( ParticleColumnHeaders::RVectorZY )];
        z3 = particleData[particleType][i][static_cast< int >( ParticleColumnHeaders::RVectorZZ )];

        particleRVectors[index][0][0] = x1;
        particleRVectors[index][0][1] = y1;
        particleRVectors[index][0][2] = z1;
        particleRVectors[index][1][0] = x2;
        particleRVectors[index][1][1] = y2;
        particleRVectors[index][1][2] = z2;
        particleRVectors[index][2][0] = x3;
        particleRVectors[index][2][1] = y3;
        particleRVectors[index][2][2] = z3;
        particleVolume[index] = 8.0*std::fabs( -(x3*y2*z1) + x2*y3*z1 + x3*y1*z2 - x1*y3*z2 - x2*y1*z3 + x1*y2*z3 );
      }
      else
      {
        GEOS_ERROR( "Invalid particle type specification! Cannot determine particle volume, aborting." );
      }

      // Material Direction
      particleMaterialDirection[index][0] = particleData[particleType][i][static_cast< int >( ParticleColumnHeaders::MaterialDirectionX )];
      particleMaterialDirection[index][1] = particleData[particleType][i][static_cast< int >( ParticleColumnHeaders::MaterialDirectionY )];
      particleMaterialDirection[index][2] = particleData[particleType][i][static_cast< int >( ParticleColumnHeaders::MaterialDirectionZ )];

      particleSurfaceNormal[index][0] = particleData[particleType][i][static_cast< int >( ParticleColumnHeaders::SurfaceNormalX )];
      particleSurfaceNormal[index][1] = particleData[particleType][i][static_cast< int >( ParticleColumnHeaders::SurfaceNormalY )];
      particleSurfaceNormal[index][2] = particleData[particleType][i][static_cast< int >( ParticleColumnHeaders::SurfaceNormalZ )];

      particleSurfacePosition[index][0] = particleData[particleType][i][static_cast< int >( ParticleColumnHeaders::SurfacePositionX )];
      particleSurfacePosition[index][1] = particleData[particleType][i][static_cast< int >( ParticleColumnHeaders::SurfacePositionY )];
      particleSurfacePosition[index][2] = particleData[particleType][i][static_cast< int >( ParticleColumnHeaders::SurfacePositionZ )];

      particleSurfaceTraction[index][0] = particleData[particleType][i][static_cast< int >( ParticleColumnHeaders::SurfaceTractionX )];
      particleSurfaceTraction[index][1] = particleData[particleType][i][static_cast< int >( ParticleColumnHeaders::SurfaceTractionY )];
      particleSurfaceTraction[index][2] = particleData[particleType][i][static_cast< int >( ParticleColumnHeaders::SurfaceTractionZ )];

      // Increment index
      index++;
    }
    particleBlock.setParticleID( particleID );
    particleBlock.setParticleCenter( particleCenter );
    particleBlock.setParticleVelocity( particleVelocity );
    particleBlock.setParticleMaterialDirection( particleMaterialDirection );
    particleBlock.setParticleGroup( particleGroup );
    particleBlock.setParticleSurfaceFlag( particleSurfaceFlag );
    particleBlock.setParticleDamage( particleDamage );
    particleBlock.setParticlePorosity( particlePorosity );
    particleBlock.setParticleTemperature( particleTemperature );
    particleBlock.setParticleStrengthScale( particleStrengthScale );
    particleBlock.setParticleVolume( particleVolume );
    particleBlock.setParticleRVectors( particleRVectors );
    particleBlock.setParticleSurfaceNormal( particleSurfaceNormal );
    particleBlock.setParticleSurfacePosition( particleSurfacePosition );
    particleBlock.setParticleSurfaceTraction( particleSurfaceTraction );
  } // loop over particle blocks

  // Resize particle regions
  int numParticles = 0;
  particleManager.forParticleRegions< ParticleRegion >( [&]( auto & particleRegion )
  {
    string_array particleBlockNames = particleRegion.getParticleBlockNames();
    std::string material = particleRegion.getMaterialList()[0]; // We will assume that the material list for a region contains only one
                                                                // material since MPM will only be doing single phase mechanics for now
    int size = 0;
    for( auto i=0; i<particleBlockNames.size(); i++ )
    {
      size += sizeMap[particleBlockNames[i]];
    }
    numParticles += size;
    particleRegion.resize( size );
  } );

  particleManager.resize( numParticles ); // All this does is change m_size for the particleManager, gives a convenient way to get the total
                                          // number of particles
  GEOS_LOG_RANK( "Total number of particles on this rank: " << particleManager.size() );
}

void ParticleMeshGenerator::postProcessInput()
{
  //GEOS_LOG_RANK_0( "Someone called ParticleMeshGenerator::postProcessInput!" );
}

REGISTER_CATALOG_ENTRY( ParticleMeshGeneratorBase, ParticleMeshGenerator, string const &, Group * const )

} /* namespace geos */
