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

#include "mainInterface/GeosxState.hpp"
#include "common/initializeEnvironment.hpp"

#include "ParticleBlockManager.hpp"
#include "mesh/mpiCommunications/SpatialPartition.hpp"

#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"

#include "physicsSolvers/solidMechanics/MPMSolverEnums.hpp"

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

  // registerWrapper( viewKeyStruct::headerFilePathString(), &m_headerFilePath ).
  //   setInputFlag( InputFlags::REQUIRED ).
  //   setRestartFlags( RestartFlags::NO_WRITE ).
  //   setDescription( "path to the header file" );

  registerWrapper( viewKeyStruct::particleBlockNamesString(), &m_blockNames ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRefArray ).
    setInputFlag( InputFlags::REQUIRED ).
    setSizedFromParent( 0 ).
    setDescription( "Names of each particle block" );

  registerWrapper( viewKeyStruct::particleTypesString(), &m_particleTypes ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRefArray ).
    setInputFlag( InputFlags::REQUIRED ).
    setSizedFromParent( 0 ).
    setDescription( "Particle types of each particle block" );
}

Group * ParticleMeshGenerator::createChild( string const & GEOS_UNUSED_PARAM( childKey ),
                                            string const & GEOS_UNUSED_PARAM( childName ) )
{
  return nullptr;
}

void ParticleMeshGenerator::postInputInitialization()
{
  ParticleMeshGeneratorBase::postInputInitialization();

  GEOS_ERROR_IF( m_blockNames.size() == 0, "No particle blocks were specified! Must specify at least one particle block." );
  GEOS_ERROR_IF( m_blockNames.size() != m_particleTypes.size(), "The particle block and type lists must have the same size." );

  for( int i = 0; i < static_cast< int >( m_particleTypes.size() ); ++i )
  {
    for( int j = 0; j < EnumSize< ParticleType >; ++j )
    {
      if( m_particleTypes[i] == EnumStrings< ParticleType >::toString( static_cast< ParticleType >( j ) ) )
      {
        break;
      }

      if( j == EnumSize< ParticleType >-1 )
      {
        GEOS_ERROR( "No particle type of " << m_particleTypes[i] << " in particleTypes input. Available options are SinglePoint (0), SinglePointBSpline (1), CPDI (2), CPTI (3), CPDI2 (4)." ); // TODO
                                                                                                                                                                                                // Make
                                                                                                                                                                                                // available
                                                                                                                                                                                                // options
                                                                                                                                                                                                // dynamic
      }
    }
  }
}

void ParticleMeshGenerator::fillParticleBlockManager( ParticleBlockManager & particleBlockManager, ParticleManager & particleManager, SpatialPartition const & partition )
{
  GEOS_MARK_FUNCTION;

  // TODO: modify this function for particle restarts (e.g. no need to read read from particle file, just initialize particle blocks)
  // CommandLineOptions const & opts = getGlobalState().getCommandLineOptions();
  // if( opts.beginFromRestart )
  // {
  //   return;
  // }

  int numBlocks = m_blockNames.size();

  // This should probably be handled elsewhere:
  std::vector< int > blockMaterialIndices( numBlocks );
  for( int b = 0; b < numBlocks; b++ )
  {
    ParticleBlock & particleBlock = particleBlockManager.registerParticleBlock( m_blockNames[b] );
    particleBlock.setParticleType( EnumStrings< ParticleType >::fromString( m_particleTypes[b] ) );
  }

  // Collect material names for each particle block
  // Assumes particleRegions loop follows material name ordering of input file
  int regionIndex = 0;
  particleManager.forParticleRegions< ParticleRegion >( [&]( auto & particleRegion )
  {
    auto & subRegions = particleRegion.getSubRegions();
    for( int r=0; r < subRegions.size(); ++r )
    {
      ParticleSubRegion & subRegion = dynamicCast< ParticleSubRegion & >( *subRegions[r] );
      GEOS_LOG_RANK( particleRegion.getName() << " | " << subRegion.getName() );
    }

    string_array particleBlockNames = particleRegion.getParticleBlockNames();
    for( auto & particleBlockName : particleBlockNames )
    {
      for( int b = 0; b < numBlocks; b++ )
      {
        if( particleBlockName == m_blockNames[b] )
        {
          blockMaterialIndices[b] = regionIndex;
          break;
        }
      }
    }
    regionIndex++;
  } );

  std::ifstream particleFile( m_particleFilePath );
  
  GEOS_ERROR_IF( !particleFile.is_open(), "Could not open MPM particle file: " << m_particleFilePath );

  // Read column headers for particle data
  int numColumnHeaders = 0;
  std::string line; // initialize line variable
  std::vector< std::vector< std::vector< double > > > particleData( numBlocks );
  std::map< int, int > columnHeaderMap;
  {
    std::getline( particleFile, line );
    std::string token;

    std::istringstream lineStream( line );
    while( std::getline( lineStream, token, '\t' ) )
    {
      columnHeaderMap.insert( std::pair< int, int >( static_cast< int >( EnumStrings< ParticleColumnHeaders >::fromString( token ) ), numColumnHeaders ) );
      numColumnHeaders++;
    }
  }

  auto requireColumn = [&]( ParticleColumnHeaders h )
  {
    GEOS_ERROR_IF( columnHeaderMap.find( static_cast< int >( h ) ) == columnHeaderMap.end(),
                  "Missing required MPM particle-file column: "
                  << EnumStrings< ParticleColumnHeaders >::toString( h ) );
  };

  requireColumn( ParticleColumnHeaders::ID );
  requireColumn( ParticleColumnHeaders::PositionX );
  requireColumn( ParticleColumnHeaders::PositionY );
  requireColumn( ParticleColumnHeaders::PositionZ );
  requireColumn( ParticleColumnHeaders::RVectorXX ); // For single point we currently use this to input particle volume 
  // I don't like that we should optionally check if particleVolume field is defined and use that regardless of rvectors if instead, and warn if the user tries specifying both

  // Read in particle data
  int lineNumber = 1; // Since column header takes one line
  while( particleFile.peek() != EOF )
  {
    std::getline( particleFile, line );
    std::istringstream lineStream( line );

    std::vector< double > lineData; // TODO: Not great because we cast all input as doubles, but it all gets re-cast later so maybe it's
                                    // fine.

    // Read line from particle file and parse columns
    double value;
    int numLineColumns = 0;
    while( lineStream >> value )
    {
      lineData.push_back( value );
      numLineColumns++;
    }

    GEOS_ERROR_IF( numLineColumns != numColumnHeaders,
                   "Particle file line " << lineNumber << " has a different number of terms than the column headers! Was " << numLineColumns << " but should be " << numColumnHeaders );
    lineNumber++;

    // If particle is inside partition add to particleData otherwise ignore and continue parsing file
    bool inPartition = partition.isCoordInPartition( lineData[ columnHeaderMap.at( static_cast< int >( ParticleColumnHeaders::PositionX ) ) ], 0 ) &&
                       partition.isCoordInPartition( lineData[ columnHeaderMap.at( static_cast< int >( ParticleColumnHeaders::PositionY ) ) ], 1 ) &&
                       partition.isCoordInPartition( lineData[ columnHeaderMap.at( static_cast< int >( ParticleColumnHeaders::PositionZ ) ) ], 2 );
    if( !inPartition )
    {
      continue;
    }

    // Reformat particle data and apply defaults to particle fields not specified in particle file
    std::vector< double > lineDataInside;
    // Presently get the number of enums using a sentinel (COUNT), but ideally this should be done more robustly
    for( int c = 0; c < EnumSize< ParticleColumnHeaders >; c++ )
    {
      if( columnHeaderMap.find( c ) != columnHeaderMap.end() )
      {
        lineDataInside.push_back( lineData[ columnHeaderMap.at( c ) ] );
        continue;
      }

      // Apply default value
      double defaultValue;
      switch( static_cast< ParticleColumnHeaders >( c ) )
      {
        case ParticleColumnHeaders::StrengthScale:
        case ParticleColumnHeaders::MaterialDirectionXX:
        case ParticleColumnHeaders::MaterialDirectionYY:
        case ParticleColumnHeaders::MaterialDirectionZZ:
          defaultValue = 1.0;
          break;
        case ParticleColumnHeaders::Temperature:
          defaultValue = 300.0;
          break;
        case ParticleColumnHeaders::ParticleType:
          defaultValue = 2.0;
          break;
        case ParticleColumnHeaders::CZTag:
        case ParticleColumnHeaders::SurfaceFlag:
        case ParticleColumnHeaders::MaterialType:
        case ParticleColumnHeaders::ContactGroup:
        case ParticleColumnHeaders::Damage:
        case ParticleColumnHeaders::Porosity:
        case ParticleColumnHeaders::VelocityX:
        case ParticleColumnHeaders::VelocityY:
        case ParticleColumnHeaders::VelocityZ:
        case ParticleColumnHeaders::MaterialDirectionXY:
        case ParticleColumnHeaders::MaterialDirectionXZ:
        case ParticleColumnHeaders::MaterialDirectionYX:
        case ParticleColumnHeaders::MaterialDirectionYZ:
        case ParticleColumnHeaders::MaterialDirectionZX:
        case ParticleColumnHeaders::MaterialDirectionZY:
        case ParticleColumnHeaders::SurfaceNormalX:
        case ParticleColumnHeaders::SurfaceNormalY:
        case ParticleColumnHeaders::SurfaceNormalZ:
        case ParticleColumnHeaders::SurfacePositionX:
        case ParticleColumnHeaders::SurfacePositionY:
        case ParticleColumnHeaders::SurfacePositionZ:
        case ParticleColumnHeaders::SurfaceTractionX:
        case ParticleColumnHeaders::SurfaceTractionY:
        case ParticleColumnHeaders::SurfaceTractionZ:
        case ParticleColumnHeaders::TemperatureRate:
          defaultValue = 0.0;
          break;
        default:
          GEOS_ERROR( EnumStrings< ParticleColumnHeaders >::toString( static_cast< ParticleColumnHeaders >( c ) ) << " must be specified in particle file!" );
          break;
      }

      lineDataInside.push_back( defaultValue );
    }

    int materialTypeIndex = lineDataInside[static_cast< int >( ParticleColumnHeaders::MaterialType )];
    int particleTypeIndex = lineDataInside[static_cast< int >( ParticleColumnHeaders::ParticleType )];

    // Match particle to block using material and particle types
    int blockIndex = -1;
    for( int b = 0; b < numBlocks; b++ )
    {
      if( materialTypeIndex == blockMaterialIndices[b] &&
          particleTypeIndex == static_cast< int >( EnumStrings< ParticleType >::fromString( m_particleTypes[b] ) ) )
      {
        blockIndex = b;
      }
    }

    GEOS_ERROR_IF( blockIndex < 0,
                   "Particle at line  " << lineNumber <<  " with type and material indices of " << particleTypeIndex << " and " << materialTypeIndex <<
      " respectively, does not match any particle block!" );

    particleData[blockIndex].push_back( lineDataInside );
  }

  // Loop over blocks and assign particle attributes
  for( int b = 0; b < numBlocks; b++ )
  {
    std::string blockName = m_blockNames[b];
    ParticleBlock & particleBlock = particleBlockManager.getParticleBlock( blockName );

    int npInBlock = particleData[b].size();
    particleBlock.resize( npInBlock );
    array1d< globalIndex > particleID( npInBlock );
    array2d< real64 > particleCenter( npInBlock, 3 );
    array2d< real64 > particleVelocity( npInBlock, 3 );
    array3d< real64 > particleMaterialDirection( npInBlock, 3, 3 );
    array1d< int > particleGroup( npInBlock );
    array1d< integer > particleSurfaceFlag( npInBlock );
    array1d< real64 > particleDamage( npInBlock );
    array1d< real64 > particlePorosity( npInBlock );
    array1d< real64 > particleTemperature( npInBlock );
    array1d< real64 > particleTemperatureRate( npInBlock );
    array1d< real64 > particleVolume( npInBlock );
    array1d< real64 > particleStrengthScale( npInBlock );
    array3d< real64 > particleRVectors( npInBlock, 3, 3 ); // TODO: Flatten the r-vector array into a 1x9 for each particle
    array2d< real64 > particleSurfaceNormal( npInBlock, 3 );
    array2d< real64 > particleSurfacePosition( npInBlock, 3 );
    array2d< real64 > particleSurfaceTraction( npInBlock, 3 );
    array1d< int > particleCZTag( npInBlock );

    // Populate particle fields with data
    for( localIndex i = 0; i < npInBlock; i++ )
    {
      // Global ID
      particleID[i] = particleData[b][i][static_cast< int >( ParticleColumnHeaders::ID )];

      // Position
      particleCenter[i][0] = particleData[b][i][static_cast< int >( ParticleColumnHeaders::PositionX )];
      particleCenter[i][1] = particleData[b][i][static_cast< int >( ParticleColumnHeaders::PositionY )];
      particleCenter[i][2] = particleData[b][i][static_cast< int >( ParticleColumnHeaders::PositionZ )];

      // Velocity
      particleVelocity[i][0] = particleData[b][i][static_cast< int >( ParticleColumnHeaders::VelocityX )];
      particleVelocity[i][1] = particleData[b][i][static_cast< int >( ParticleColumnHeaders::VelocityY )];
      particleVelocity[i][2] = particleData[b][i][static_cast< int >( ParticleColumnHeaders::VelocityZ )];

      // Material (set above) is [10]

      // Group
      particleGroup[i] = particleData[b][i][static_cast< int >( ParticleColumnHeaders::ContactGroup )];

      // surfaceFlag
      particleSurfaceFlag[i] = particleData[b][i][static_cast< int >( ParticleColumnHeaders::SurfaceFlag )];

      // Damage
      particleDamage[i] = particleData[b][i][static_cast< int >( ParticleColumnHeaders::Damage )];

      // Porosity
      particlePorosity[i] = particleData[b][i][static_cast< int >( ParticleColumnHeaders::Porosity )];

      // Temperature
      particleTemperature[i] = particleData[b][i][static_cast< int >( ParticleColumnHeaders::Temperature )];

      // Temperature Rate
      particleTemperatureRate[i] = particleData[b][i][static_cast< int >( ParticleColumnHeaders::TemperatureRate )];

      // Strength Scale
      particleStrengthScale[i] = particleData[b][i][static_cast< int >( ParticleColumnHeaders::StrengthScale )];

      // Cohesive Zone ID
      particleCZTag[i] = particleData[b][i][static_cast< int >( ParticleColumnHeaders::CZTag )];

      // Volume and R-Vectors
      switch( particleBlock.getParticleType() )
      {
        case ParticleType::SinglePoint:
        {
          particleVolume[i] = particleData[b][i][static_cast< int >( ParticleColumnHeaders::RVectorXX )];

          GEOS_ERROR_IF( !std::isfinite( particleVolume[i] ) || particleVolume[i] <= 0.0,
               "MPM particle volume must be finite and positive." );

          real64 a = LvArray::math::pow( particleVolume[i], 1.0/3.0 );
          particleRVectors[i][0][0] = a;
          particleRVectors[i][0][1] = 0.0;
          particleRVectors[i][0][2] = 0.0;
          particleRVectors[i][1][0] = 0.0;
          particleRVectors[i][1][1] = a;
          particleRVectors[i][1][2] = 0.0;
          particleRVectors[i][2][0] = 0.0;
          particleRVectors[i][2][1] = 0.0;
          particleRVectors[i][2][2] = a;
          break;
        }
        case ParticleType::SinglePointBSpline:
        {
          // SinglePointBSpline uses particle centers for MPM interpolation.
          // CPDI-style r-vectors are retained for particle-domain output and
          // the existing particle-volume initialization convention.
          double x1, y1, z1, x2, y2, z2, x3, y3, z3;
          x1 = particleData[b][i][static_cast< int >( ParticleColumnHeaders::RVectorXX )];
          y1 = particleData[b][i][static_cast< int >( ParticleColumnHeaders::RVectorXY )];
          z1 = particleData[b][i][static_cast< int >( ParticleColumnHeaders::RVectorXZ )];
          x2 = particleData[b][i][static_cast< int >( ParticleColumnHeaders::RVectorYX )];
          y2 = particleData[b][i][static_cast< int >( ParticleColumnHeaders::RVectorYY )];
          z2 = particleData[b][i][static_cast< int >( ParticleColumnHeaders::RVectorYZ )];
          x3 = particleData[b][i][static_cast< int >( ParticleColumnHeaders::RVectorZX )];
          y3 = particleData[b][i][static_cast< int >( ParticleColumnHeaders::RVectorZY )];
          z3 = particleData[b][i][static_cast< int >( ParticleColumnHeaders::RVectorZZ )];
          particleRVectors[i][0][0] = x1;
          particleRVectors[i][0][1] = y1;
          particleRVectors[i][0][2] = z1;
          particleRVectors[i][1][0] = x2;
          particleRVectors[i][1][1] = y2;
          particleRVectors[i][1][2] = z2;
          particleRVectors[i][2][0] = x3;
          particleRVectors[i][2][1] = y3;
          particleRVectors[i][2][2] = z3;
          particleVolume[i] = 8.0*LvArray::math::abs( -(x3*y2*z1) + x2*y3*z1 + x3*y1*z2 - x1*y3*z2 - x2*y1*z3 + x1*y2*z3 );
          break;
        }
        case ParticleType::CPDI:
        {
          double x1, y1, z1, x2, y2, z2, x3, y3, z3;
          x1 = particleData[b][i][static_cast< int >( ParticleColumnHeaders::RVectorXX )];
          y1 = particleData[b][i][static_cast< int >( ParticleColumnHeaders::RVectorXY )];
          z1 = particleData[b][i][static_cast< int >( ParticleColumnHeaders::RVectorXZ )];
          x2 = particleData[b][i][static_cast< int >( ParticleColumnHeaders::RVectorYX )];
          y2 = particleData[b][i][static_cast< int >( ParticleColumnHeaders::RVectorYY )];
          z2 = particleData[b][i][static_cast< int >( ParticleColumnHeaders::RVectorYZ )];
          x3 = particleData[b][i][static_cast< int >( ParticleColumnHeaders::RVectorZX )];
          y3 = particleData[b][i][static_cast< int >( ParticleColumnHeaders::RVectorZY )];
          z3 = particleData[b][i][static_cast< int >( ParticleColumnHeaders::RVectorZZ )];
          particleRVectors[i][0][0] = x1;
          particleRVectors[i][0][1] = y1;
          particleRVectors[i][0][2] = z1;
          particleRVectors[i][1][0] = x2;
          particleRVectors[i][1][1] = y2;
          particleRVectors[i][1][2] = z2;
          particleRVectors[i][2][0] = x3;
          particleRVectors[i][2][1] = y3;
          particleRVectors[i][2][2] = z3;
          particleVolume[i] = 8.0*LvArray::math::abs( -(x3*y2*z1) + x2*y3*z1 + x3*y1*z2 - x1*y3*z2 - x2*y1*z3 + x1*y2*z3 );
          break;
        }
        case ParticleType::CPTI:
        {
          GEOS_ERROR( "CPTI particle type is not implemented!" );
          break;
        }
        case ParticleType::CPDI2:
        {
          GEOS_ERROR( "CPDI2 particle type is not implemented!" );
          break;
        }
        default:
        {
          GEOS_ERROR( "Invalid particle type: " << particleBlock.getParticleType() );
          break;
        }
      }

      // Material Direction
      // For convenient indexing we store the transpose of the particle material directions so the X axis is easy accessed from [0]
      particleMaterialDirection[i][0][0] = particleData[b][i][static_cast< int >( ParticleColumnHeaders::MaterialDirectionXX )];
      particleMaterialDirection[i][0][1] = particleData[b][i][static_cast< int >( ParticleColumnHeaders::MaterialDirectionXY )];
      particleMaterialDirection[i][0][2] = particleData[b][i][static_cast< int >( ParticleColumnHeaders::MaterialDirectionXZ )];
      particleMaterialDirection[i][1][0] = particleData[b][i][static_cast< int >( ParticleColumnHeaders::MaterialDirectionYX )];
      particleMaterialDirection[i][1][1] = particleData[b][i][static_cast< int >( ParticleColumnHeaders::MaterialDirectionYY )];
      particleMaterialDirection[i][1][2] = particleData[b][i][static_cast< int >( ParticleColumnHeaders::MaterialDirectionYZ )];
      particleMaterialDirection[i][2][0] = particleData[b][i][static_cast< int >( ParticleColumnHeaders::MaterialDirectionZX )];
      particleMaterialDirection[i][2][1] = particleData[b][i][static_cast< int >( ParticleColumnHeaders::MaterialDirectionZY )];
      particleMaterialDirection[i][2][2] = particleData[b][i][static_cast< int >( ParticleColumnHeaders::MaterialDirectionZZ )];

      particleSurfaceNormal[i][0] = particleData[b][i][static_cast< int >( ParticleColumnHeaders::SurfaceNormalX )];
      particleSurfaceNormal[i][1] = particleData[b][i][static_cast< int >( ParticleColumnHeaders::SurfaceNormalY )];
      particleSurfaceNormal[i][2] = particleData[b][i][static_cast< int >( ParticleColumnHeaders::SurfaceNormalZ )];

      particleSurfacePosition[i][0] = particleData[b][i][static_cast< int >( ParticleColumnHeaders::SurfacePositionX )];
      particleSurfacePosition[i][1] = particleData[b][i][static_cast< int >( ParticleColumnHeaders::SurfacePositionY )];
      particleSurfacePosition[i][2] = particleData[b][i][static_cast< int >( ParticleColumnHeaders::SurfacePositionZ )];

      particleSurfaceTraction[i][0] = particleData[b][i][static_cast< int >( ParticleColumnHeaders::SurfaceTractionX )];
      particleSurfaceTraction[i][1] = particleData[b][i][static_cast< int >( ParticleColumnHeaders::SurfaceTractionY )];
      particleSurfaceTraction[i][2] = particleData[b][i][static_cast< int >( ParticleColumnHeaders::SurfaceTractionZ )];
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
    particleBlock.setParticleTemperatureRate( particleTemperatureRate );
    particleBlock.setParticleStrengthScale( particleStrengthScale );
    particleBlock.setParticleVolume( particleVolume );
    particleBlock.setParticleRVectors( particleRVectors );
    particleBlock.setParticleSurfaceNormal( particleSurfaceNormal );
    particleBlock.setParticleSurfacePosition( particleSurfacePosition );
    particleBlock.setParticleSurfaceTraction( particleSurfaceTraction );
    particleBlock.setParticleCZTag( particleCZTag );
  } // loop over particle blocks

  // Resize particle regions
  int numParticles = 0;
  particleManager.forParticleRegions< ParticleRegion >( [&]( auto & particleRegion )
  {
    string_array particleBlockNames = particleRegion.getParticleBlockNames();

    int size = 0;
    for( auto & particleBlockName : particleBlockNames )
    {
      for( int b = 0; b < numBlocks; b++ )
      {
        if( particleBlockName == m_blockNames[b] )
        {
          size += particleData[b].size();
          break;
        }
      }
    }
    numParticles += size;
    particleRegion.resize( size );
  } );

  particleManager.resize( numParticles ); // All this does is change m_size for the particleManager, gives a convenient way to get the total
                                          // number of particles
  GEOS_LOG_RANK( "Total number of particles on this rank: " << particleManager.size() );
}

REGISTER_CATALOG_ENTRY( ParticleMeshGeneratorBase, ParticleMeshGenerator, string const &, Group * const )

} /* namespace geos */
