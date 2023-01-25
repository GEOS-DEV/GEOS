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


#include "ParticleSubRegion.hpp"

#include "common/TypeDispatch.hpp"
#include "mesh/MeshLevel.hpp"

namespace geosx
{
using namespace dataRepository;
using namespace constitutive;

ParticleSubRegion::ParticleSubRegion( string const & name, Group * const parent ):
  ParticleSubRegionBase( name, parent )
{
  registerWrapper( viewKeyStruct::constitutiveGroupingString(), &m_constitutiveGrouping ).
    setSizedFromParent( 0 );
}

ParticleSubRegion::~ParticleSubRegion()
{
  // Left blank
}

void ParticleSubRegion::updateMaps()
{
  for( localIndex p=0; p<this->size(); p++ )
  {
    m_localToGlobalMap[p] = m_particleID[p];
  }
  this->constructGlobalToLocalMap();
}

void ParticleSubRegion::setParticleRank(int rank, int np)
{
  for(int i=0; i<np; i++)
  {
    m_particleRank[i] = rank;
  }
}

int ParticleSubRegion::numNodesMappedTo()
{
  switch( m_particleType ) // TODO: Make this a ParticleSubRegion method
  {
    case ParticleType::SinglePoint:
    {
      return 8;
      break;
    }
    case ParticleType::CPDI:
    {
      return 64;
      break;
    }
    default:
    {
      GEOSX_ERROR( "Particle type \"" << m_particleType << "\" is not yet supported." );
      return -1; // Needed to compile
    }
  }
}

void ParticleSubRegion::copyFromParticleBlock( ParticleBlockABC & particleBlock )
{
  // Defines the (unique) particle type of this cell particle region,
  // and its associated number of nodes, edges, faces.
  m_particleRank.resize(particleBlock.size());
  m_particleType = particleBlock.getParticleType();
  m_particleID = particleBlock.getParticleID();
  m_particleGroup = particleBlock.getParticleGroup();
  m_particleDamage = particleBlock.getParticleDamage();
  m_particleCenter = particleBlock.getParticleCenter();
  m_particleVelocity = particleBlock.getParticleVelocity();
  m_particleVolume = particleBlock.getParticleVolume();
  m_particleRVectors = particleBlock.getParticleRVectors();
  m_hasRVectors = particleBlock.hasRVectors();


  // We call the `resize` member function of the cell to (nodes, edges, faces) relations,
  // before calling the `ParticleSubRegion::resize` in order to keep the first dimension.
  // Be careful when refactoring.
  this->resize( particleBlock.numParticles() );

  this->m_localToGlobalMap = particleBlock.localToGlobalMap();

  this->constructGlobalToLocalMap();
  particleBlock.forExternalProperties( [&]( WrapperBase & wrapper )
  {
    types::dispatch( types::StandardArrays{}, wrapper.getTypeId(), true, [&]( auto array )
    {
      using ArrayType = decltype( array );
      Wrapper< ArrayType > & wrapperT = Wrapper< ArrayType >::cast( wrapper );
      this->registerWrapper( wrapper.getName(), &wrapperT.reference() );
    } );
  } );
}

void ParticleSubRegion::flagOutOfRangeParticles( std::array< real64, 3 > const & xGlobalMin,
                                                 std::array< real64, 3 > const & xGlobalMax,
                                                 std::array< real64, 3 > const & hEl,
                                                 arrayView1d< int > const isBad )
{
  // For CPDI
  int signs[8][3] = { { 1,  1,  1},
                      { 1,  1, -1},
                      { 1, -1,  1},
                      { 1, -1, -1},
                      {-1,  1,  1},
                      {-1,  1, -1},
                      {-1, -1,  1},
                      {-1, -1, -1} };
  
  std::array< real64, 3 > globalMin; // including buffer cells
  std::array< real64, 3 > globalMax; // including buffer cells
  for( int i=0; i<3; i++)
  {
    globalMin[i] = xGlobalMin[i] - hEl[i];
    globalMax[i] = xGlobalMax[i] + hEl[i];
  }

  for( localIndex const p : this->nonGhostIndices() )
  {
    arraySlice1d< real64 > const & p_x = m_particleCenter[p];

    switch( m_particleType )
    {
      case ParticleType::SinglePoint:
      {
        for( int i=0; i<3; i++ )
        {
          if( p_x[i] < globalMin[i] || globalMax[i] < p_x[i] )
          {
            isBad[p] = 1;
            break;
          }
        }
        break;
      }
      case ParticleType::CPDI:
      {
        for(int cornerIndex=0; cornerIndex<8; cornerIndex++)
        {
          for(int i=0; i<3; i++)
          {
            real64 cornerPositionComponent = p_x[i] + signs[cornerIndex][0] * m_particleRVectors[p][0][i] + signs[cornerIndex][1] * m_particleRVectors[p][1][i] + signs[cornerIndex][2] * m_particleRVectors[p][2][i];
            if( cornerPositionComponent < globalMin[i] || globalMax[i] < cornerPositionComponent )
            {
              isBad[p] = 1;
              break;
            }
          }
          if( isBad[p] == 1 )
          {
            break;
          }
        }
        break;
      }
      default:
      {
        GEOSX_ERROR( "Particle type \"" << m_particleType << "\" is not yet supported." );
      }
    }
  }
}

void ParticleSubRegion::computeRVectors( int const p,
                                         arraySlice2d< real64 > const F,
                                         arraySlice2d< real64 > const initialRVectors )
{
  if(m_hasRVectors)
  {
    for(int i=0; i<3; i++)
    {
      for(int j=0; j<3; j++)
      {
        m_particleRVectors[p][i][j] = F[j][0]*initialRVectors[i][0] + F[j][1]*initialRVectors[i][1] + F[j][2]*initialRVectors[i][2];
      }
    }
  }
}

void ParticleSubRegion::cpdiDomainScaling( real64 lCrit,
                                           int m_planeStrain )
{
 for( int p=0; p<this->size(); p++ )
 {
    arraySlice1d< real64 > r1 = m_particleRVectors[p][0];
    arraySlice1d< real64 > r2 = m_particleRVectors[p][1];
    arraySlice1d< real64 > r3 = m_particleRVectors[p][2];

    if( m_planeStrain ) // 2D cpdi domain scaling
    {
      // Initialize l-vectors.  Eq. 8a-d in the CPDI domain scaling paper.
      real64 l[2][3];
      for(int i=0; i<3; i++)
      {
        l[0][i] = r1[i] + r2[i]; // la
        l[1][i] = r1[i] - r2[i]; // lb
      }

      // scale l-vectors if needed.  Eq. 9 in the CPDI domain scaling paper.
      bool scale = false;
      for( int i = 0 ; i < 2 ; i++ )
      {
        real64 lLength = sqrt(l[i][0] * l[i][0] + l[i][1] * l[i][1] + l[i][2] * l[i][2]);
        if( lLength > lCrit )
        {
          l[i][0] *= lCrit / lLength;
          l[i][1] *= lCrit / lLength;
          l[i][2] *= lCrit / lLength;
          scale = true;
        }
      }

      // reconstruct r-vectors.  eq. 11 in the CPDI domain scaling paper.
      if( scale )
      {
        for(int i=0; i<3; i++)
        {
          r1[i] = 0.5 * (l[0][i] + l[1][i]);
          r2[i] = 0.5 * (l[0][i] - l[1][i]);
        }
      }
    }
    else // 3D cpdi domain scaling
    {
      // Initialize l-vectors.  Eq. 8a-d in the CPDI domain scaling paper.
      real64 l[4][3];
      for(int i=0; i<3; i++)
      {
        l[0][i] = r1[i] + r2[i] + r3[i]; // la
        l[1][i] = r1[i] - r2[i] + r3[i]; // lb
        l[2][i] = r2[i] - r1[i] + r3[i]; // lc
        l[3][i] = r3[i] - r1[i] - r2[i]; // ld
      }

      // scale l vectors if needed.  Eq. 9 in the CPDI domain scaling paper.
      bool scale = false;
      for( int i = 0 ; i < 4 ; i++ )
      {
        real64 lLength = sqrt(l[i][0] * l[i][0] + l[i][1] * l[i][1] + l[i][2] * l[i][2]);
        if( lLength > lCrit )
        {
          l[i][0] *= lCrit / lLength;
          l[i][1] *= lCrit / lLength;
          l[i][2] *= lCrit / lLength;
          scale = true;
        }
      }

      // reconstruct r vectors.  eq. 11 in the CPDI domain scaling paper.
      if( scale )
      {
        for(int i=0; i<3; i++)
        {
          r1[i] = 0.25 * ( l[0][i] + l[1][i] - l[2][i] - l[3][i] );
          r2[i] = 0.25 * ( l[0][i] - l[1][i] + l[2][i] - l[3][i] );
          r3[i] = 0.25 * ( l[0][i] + l[1][i] + l[2][i] + l[3][i] );
        }
      }
    }
 }
}

void ParticleSubRegion::getMappedNodes( int const p,
                                        std::array< real64, 3 > const & xMin,
                                        std::array< real64, 3 > const & hx,
                                        array3d< int > const & ijkMap,
                                        arrayView1d< localIndex > const nodeIDs )
{
  arraySlice1d< real64 > const & p_x = m_particleCenter[p];

  switch( m_particleType )
  {
    case ParticleType::SinglePoint:
    {
      // get IJK associated with particle center
      std::vector<int> centerIJK;
      centerIJK.resize(3);
      for(int i=0; i<3; i++)
      {
        centerIJK[i] = std::floor( ( p_x[i] - xMin[i] ) / hx[i] );
      }

      int node = 0;
      for(int i=0; i<2; i++)
      {
        for(int j=0; j<2; j++)
        {
          for(int k=0; k<2; k++)
          {
            nodeIDs[node] = ijkMap[centerIJK[0]+i][centerIJK[1]+j][centerIJK[2]+k] ;
            node++;
          }
        }
      }

      break;
    }
    case ParticleType::CPDI:
    {
      // Precalculated things
      int signs[8][3] = { { -1, -1, -1 },
                          {  1, -1, -1 },
                          {  1,  1, -1 },
                          { -1,  1, -1 },
                          { -1, -1,  1 },
                          {  1, -1,  1 },
                          {  1,  1,  1 },
                          { -1,  1,  1 } };
      real64 p_r1[3], p_r2[3], p_r3[3]; // allowing 1-indexed r-vectors to persist to torture future postdocs >:)

      for(int i=0; i<3; i++)
      {
        p_r1[i] = m_particleRVectors[p][0][i];
        p_r2[i] = m_particleRVectors[p][1][i];
        p_r3[i] = m_particleRVectors[p][2][i];
      }

      // get IJK associated with each corner
      int cornerIJK[8][3]; // CPDI can map to up to 8 cells
      for(int corner=0; corner<8; corner++)
      {
        for(int i=0; i<3; i++)
        {
          real64 cornerPositionComponent = p_x[i] + signs[corner][0] * m_particleRVectors[p][0][i] + signs[corner][1] * m_particleRVectors[p][1][i] + signs[corner][2] * m_particleRVectors[p][2][i];
          cornerIJK[corner][i] = std::floor( ( cornerPositionComponent - xMin[i] ) / hx[i] ); // TODO: Temporarily store the CPDI corners since they're re-used below?
        }
      }

      // get node IDs associated with each corner from IJK map
      int node = 0;
      for(int corner=0; corner<8; corner++)
      {
        for(int i=0; i<2; i++)
        {
          for(int j=0; j<2; j++)
          {
            for(int k=0; k<2; k++)
            {
              nodeIDs[node] = ijkMap[cornerIJK[corner][0]+i][cornerIJK[corner][1]+j][cornerIJK[corner][2]+k];
              node++;
            }
          }
        }
      }

      break;
    }

    default:
    {
      GEOSX_ERROR( "Particle type \"" << m_particleType << "\" is not yet supported." );
    }
  }
}

void ParticleSubRegion::getAllWeights( int const p,
                                       std::array< real64, 3 > const & xMin,
                                       std::array< real64, 3 > const & hx,
                                       array3d< int > const & ijkMap,
                                       arrayView2d< real64, nodes::REFERENCE_POSITION_USD > const g_X,
                                       arrayView1d< localIndex > const nodeIDs,
                                       arrayView1d< real64 > const weights,
                                       arrayView2d< real64 > const gradWeights )
{
  arraySlice1d< real64 > const & p_x = m_particleCenter[p];

  switch( m_particleType )
  {
    case ParticleType::SinglePoint:
    {
      // get IJK associated with particle center
      std::vector<int> centerIJK;
      centerIJK.resize(3);
      for(int i=0; i<3; i++)
      {
        centerIJK[i] = std::floor( ( p_x[i] - xMin[i] ) / hx[i] );
      }

      // get node IDs, weights and grad weights
      int node = 0;
      int corner = ijkMap[centerIJK[0]][centerIJK[1]][centerIJK[2]];
      auto corner_x = g_X[corner];

      real64 xRel = (p_x[0] - corner_x[0]) / hx[0];
      real64 yRel = (p_x[1] - corner_x[1]) / hx[1];
      real64 zRel = (p_x[2] - corner_x[2]) / hx[2];

      for(int i=0; i<2; i++)
      {
        real64 xWeight = i*xRel + (1-i)*(1.0-xRel);
        real64 dxWeight = i/hx[0] - (1-i)/hx[0];
        for(int j=0; j<2; j++)
        {
          real64 yWeight = j*yRel + (1-j)*(1.0-yRel);
          real64 dyWeight = j/hx[1] - (1-j)/hx[1];
          for(int k=0; k<2; k++)
          {
            real64 zWeight = k*zRel + (1-k)*(1.0-zRel);
            real64 dzWeight = k/hx[2] - (1-k)/hx[2];
            nodeIDs[node] = ijkMap[centerIJK[0]+i][centerIJK[1]+j][centerIJK[2]+k] ;
            weights[node] = xWeight * yWeight * zWeight;
            gradWeights[node][0] = dxWeight * yWeight * zWeight;
            gradWeights[node][1] = xWeight * dyWeight * zWeight;
            gradWeights[node][2] = xWeight * yWeight * dzWeight;
            node++;
          }
        }
      }

      break;
    }

    case ParticleType::CPDI:
    {
      // Precalculated things
      int signs[8][3] = { { -1, -1, -1 },
                          {  1, -1, -1 },
                          {  1,  1, -1 },
                          { -1,  1, -1 },
                          { -1, -1,  1 },
                          {  1, -1,  1 },
                          {  1,  1,  1 },
                          { -1,  1,  1 } };

      real64 alpha[8][8];
      real64 cpdiVolume, oneOverV;
      real64 p_r1[3], p_r2[3], p_r3[3]; // allowing 1-indexed r-vectors to persist to torture future postdocs >:)

      for(int i=0; i<3; i++)
      {
        p_r1[i] = m_particleRVectors[p][0][i];
        p_r2[i] = m_particleRVectors[p][1][i];
        p_r3[i] = m_particleRVectors[p][2][i];
      }

      // We need this because of CPDI domain scaling
      cpdiVolume = 8.0*(-(p_r1[2]*p_r2[1]*p_r3[0]) + p_r1[1]*p_r2[2]*p_r3[0] + p_r1[2]*p_r2[0]*p_r3[1] - p_r1[0]*p_r2[2]*p_r3[1] - p_r1[1]*p_r2[0]*p_r3[2] + p_r1[0]*p_r2[1]*p_r3[2]);
      oneOverV = 1.0 / cpdiVolume;

      alpha[0][0]=oneOverV*(p_r1[2]*p_r2[1] - p_r1[1]*p_r2[2] - p_r1[2]*p_r3[1] + p_r2[2]*p_r3[1] + p_r1[1]*p_r3[2] - p_r2[1]*p_r3[2]);
      alpha[0][1]=oneOverV*(-(p_r1[2]*p_r2[0]) + p_r1[0]*p_r2[2] + p_r1[2]*p_r3[0] - p_r2[2]*p_r3[0] - p_r1[0]*p_r3[2] + p_r2[0]*p_r3[2]);
      alpha[0][2]=oneOverV*(p_r1[1]*p_r2[0] - p_r1[0]*p_r2[1] - p_r1[1]*p_r3[0] + p_r2[1]*p_r3[0] + p_r1[0]*p_r3[1] - p_r2[0]*p_r3[1]);
      alpha[1][0]=oneOverV*(p_r1[2]*p_r2[1] - p_r1[1]*p_r2[2] - p_r1[2]*p_r3[1] - p_r2[2]*p_r3[1] + p_r1[1]*p_r3[2] + p_r2[1]*p_r3[2]);
      alpha[1][1]=oneOverV*(-(p_r1[2]*p_r2[0]) + p_r1[0]*p_r2[2] + p_r1[2]*p_r3[0] + p_r2[2]*p_r3[0] - p_r1[0]*p_r3[2] - p_r2[0]*p_r3[2]);
      alpha[1][2]=oneOverV*(p_r1[1]*p_r2[0] - p_r1[0]*p_r2[1] - p_r1[1]*p_r3[0] - p_r2[1]*p_r3[0] + p_r1[0]*p_r3[1] + p_r2[0]*p_r3[1]);
      alpha[2][0]=oneOverV*(p_r1[2]*p_r2[1] - p_r1[1]*p_r2[2] + p_r1[2]*p_r3[1] - p_r2[2]*p_r3[1] - p_r1[1]*p_r3[2] + p_r2[1]*p_r3[2]);
      alpha[2][1]=oneOverV*(-(p_r1[2]*p_r2[0]) + p_r1[0]*p_r2[2] - p_r1[2]*p_r3[0] + p_r2[2]*p_r3[0] + p_r1[0]*p_r3[2] - p_r2[0]*p_r3[2]);
      alpha[2][2]=oneOverV*(p_r1[1]*p_r2[0] - p_r1[0]*p_r2[1] + p_r1[1]*p_r3[0] - p_r2[1]*p_r3[0] - p_r1[0]*p_r3[1] + p_r2[0]*p_r3[1]);
      alpha[3][0]=oneOverV*(p_r1[2]*p_r2[1] - p_r1[1]*p_r2[2] + p_r1[2]*p_r3[1] + p_r2[2]*p_r3[1] - p_r1[1]*p_r3[2] - p_r2[1]*p_r3[2]);
      alpha[3][1]=oneOverV*(-(p_r1[2]*p_r2[0]) + p_r1[0]*p_r2[2] - p_r1[2]*p_r3[0] - p_r2[2]*p_r3[0] + p_r1[0]*p_r3[2] + p_r2[0]*p_r3[2]);
      alpha[3][2]=oneOverV*(p_r1[1]*p_r2[0] - p_r1[0]*p_r2[1] + p_r1[1]*p_r3[0] + p_r2[1]*p_r3[0] - p_r1[0]*p_r3[1] - p_r2[0]*p_r3[1]);
      alpha[4][0]=oneOverV*(-(p_r1[2]*p_r2[1]) + p_r1[1]*p_r2[2] - p_r1[2]*p_r3[1] + p_r2[2]*p_r3[1] + p_r1[1]*p_r3[2] - p_r2[1]*p_r3[2]);
      alpha[4][1]=oneOverV*(p_r1[2]*p_r2[0] - p_r1[0]*p_r2[2] + p_r1[2]*p_r3[0] - p_r2[2]*p_r3[0] - p_r1[0]*p_r3[2] + p_r2[0]*p_r3[2]);
      alpha[4][2]=oneOverV*(-(p_r1[1]*p_r2[0]) + p_r1[0]*p_r2[1] - p_r1[1]*p_r3[0] + p_r2[1]*p_r3[0] + p_r1[0]*p_r3[1] - p_r2[0]*p_r3[1]);
      alpha[5][0]=oneOverV*(-(p_r1[2]*p_r2[1]) + p_r1[1]*p_r2[2] - p_r1[2]*p_r3[1] - p_r2[2]*p_r3[1] + p_r1[1]*p_r3[2] + p_r2[1]*p_r3[2]);
      alpha[5][1]=oneOverV*(p_r1[2]*p_r2[0] - p_r1[0]*p_r2[2] + p_r1[2]*p_r3[0] + p_r2[2]*p_r3[0] - p_r1[0]*p_r3[2] - p_r2[0]*p_r3[2]);
      alpha[5][2]=oneOverV*(-(p_r1[1]*p_r2[0]) + p_r1[0]*p_r2[1] - p_r1[1]*p_r3[0] - p_r2[1]*p_r3[0] + p_r1[0]*p_r3[1] + p_r2[0]*p_r3[1]);
      alpha[6][0]=oneOverV*(-(p_r1[2]*p_r2[1]) + p_r1[1]*p_r2[2] + p_r1[2]*p_r3[1] - p_r2[2]*p_r3[1] - p_r1[1]*p_r3[2] + p_r2[1]*p_r3[2]);
      alpha[6][1]=oneOverV*(p_r1[2]*p_r2[0] - p_r1[0]*p_r2[2] - p_r1[2]*p_r3[0] + p_r2[2]*p_r3[0] + p_r1[0]*p_r3[2] - p_r2[0]*p_r3[2]);
      alpha[6][2]=oneOverV*(-(p_r1[1]*p_r2[0]) + p_r1[0]*p_r2[1] + p_r1[1]*p_r3[0] - p_r2[1]*p_r3[0] - p_r1[0]*p_r3[1] + p_r2[0]*p_r3[1]);
      alpha[7][0]=oneOverV*(-(p_r1[2]*p_r2[1]) + p_r1[1]*p_r2[2] + p_r1[2]*p_r3[1] + p_r2[2]*p_r3[1] - p_r1[1]*p_r3[2] - p_r2[1]*p_r3[2]);
      alpha[7][1]=oneOverV*(p_r1[2]*p_r2[0] - p_r1[0]*p_r2[2] - p_r1[2]*p_r3[0] - p_r2[2]*p_r3[0] + p_r1[0]*p_r3[2] + p_r2[0]*p_r3[2]);
      alpha[7][2]=oneOverV*(-(p_r1[1]*p_r2[0]) + p_r1[0]*p_r2[1] + p_r1[1]*p_r3[0] + p_r2[1]*p_r3[0] - p_r1[0]*p_r3[1] - p_r2[0]*p_r3[1]);

      // get IJK associated with each corner
      int cornerIJK[8][3]; // CPDI can map to up to 8 cells
      for(int corner=0; corner<8; corner++)
      {
        for(int i=0; i<3; i++)
        {
          real64 cornerPositionComponent = p_x[i] + signs[corner][0] * m_particleRVectors[p][0][i] + signs[corner][1] * m_particleRVectors[p][1][i] + signs[corner][2] * m_particleRVectors[p][2][i];
          cornerIJK[corner][i] = std::floor( ( cornerPositionComponent - xMin[i] ) / hx[i] ); // TODO: Temporarily store the CPDI corners since they're re-used below?
        }
      }

      // get node IDs associated with each corner from IJK map, along with weights and grad weights
      // *** The order in which we access the IJK map must match the order we evaluate the shape functions! ***
      int node = 0;
      for(int corner=0; corner<8; corner++)
      {
        int cornerNode = ijkMap[cornerIJK[corner][0]][cornerIJK[corner][1]][cornerIJK[corner][2]];
        auto cornerNodePosition = g_X[cornerNode];

        real64 x, y, z;
        x = p_x[0] + signs[corner][0] * m_particleRVectors[p][0][0] + signs[corner][1] * m_particleRVectors[p][1][0] + signs[corner][2] * m_particleRVectors[p][2][0];
        y = p_x[1] + signs[corner][0] * m_particleRVectors[p][0][1] + signs[corner][1] * m_particleRVectors[p][1][1] + signs[corner][2] * m_particleRVectors[p][2][1];
        z = p_x[2] + signs[corner][0] * m_particleRVectors[p][0][2] + signs[corner][1] * m_particleRVectors[p][1][2] + signs[corner][2] * m_particleRVectors[p][2][2];

        real64 xRel = (x - cornerNodePosition[0]) / hx[0];
        real64 yRel = (y - cornerNodePosition[1]) / hx[1];
        real64 zRel = (z - cornerNodePosition[2]) / hx[2];

        for(int i=0; i<2; i++)
        {
          real64 xWeight = i * xRel + (1 - i) * (1.0 - xRel);
          for(int j=0; j<2; j++)
          {
            real64 yWeight = j * yRel + (1 - j) * (1.0 - yRel);
            for(int k=0; k<2; k++)
            {
              real64 zWeight = k * zRel + (1 - k) * (1.0 - zRel);
              real64 weight = xWeight * yWeight * zWeight;
              nodeIDs[node] = ijkMap[cornerIJK[corner][0]+i][cornerIJK[corner][1]+j][cornerIJK[corner][2]+k];
              weights[node] =  0.125 * weight;
              gradWeights[node][0] = alpha[corner][0] * weight;
              gradWeights[node][1] = alpha[corner][1] * weight;
              gradWeights[node][2] = alpha[corner][2] * weight;
              node++;
            }
          }
        }
      }

      break;
    }

    default:
    {
      GEOSX_ERROR( "Particle type \"" << m_particleType << "\" is not yet supported." );
    }
  }

}

REGISTER_CATALOG_ENTRY( ObjectManagerBase, ParticleSubRegion, string const &, Group * const )

} /* namespace geosx */
