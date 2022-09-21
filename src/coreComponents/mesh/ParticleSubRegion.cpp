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

void ParticleSubRegion::setParticleRank(int rank, int np)
{
  for(int i=0; i<np; i++)
  {
    m_particleRank[i] = rank;
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
  m_particleCenter = particleBlock.getParticleCenter();
  m_particleVelocity = particleBlock.getParticleVelocity();
  m_particleInitialVolume = particleBlock.getParticleInitialVolume();
  m_particleVolume = m_particleInitialVolume;
  //m_particleDeformationGradient.resize(particleBlock.size(),3,3); // handled by ParticleSubRegionBase constructor?
  m_particleInitialRVectors = particleBlock.getParticleInitialRVectors();
  m_particleRVectors = m_particleInitialRVectors;
  //m_particleMass.resize(particleBlock.size()); // handled by ParticleSubRegionBase constructor?
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

void ParticleSubRegion::deleteOutOfRangeParticles( std::array< real64, 3 > const & xGlobalMin,
                                                   std::array< real64, 3 > const & xGlobalMax,
                                                   std::array< real64, 3 > const & hEl )
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
  
  std::array< real64, 3 > globalMin; // including buffer nodes
  std::array< real64, 3 > globalMax; // including buffer nodes
  for( int i=0; i<3; i++)
  {
    globalMin[i] = xGlobalMin[i] - hEl[i];
    globalMax[i] = xGlobalMax[i] + hEl[i];
  }

  for( int p=this->size()-1; p>=0; p-- )
  {
    arraySlice1d< real64 > const & p_x = m_particleCenter[p];

    bool deleted = false;

    switch( m_particleType )
    {
      case ParticleType::SinglePoint:
      {
        for( int i=0; i<3; i++ )
        {
          if( p_x[i] < globalMin[i] || globalMax[i] < p_x[i] )
          {
            this->erase(p);
            deleted = true;
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
              this->erase(p);
              deleted = true;
              break;
            }
          }
          if( deleted )
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

void ParticleSubRegion::applyFtoRVectors( int const p,
                                          LvArray::ArraySlice<double, 2, 1, int> const & p_F )
{
  if(m_hasRVectors)
  {
    for(int i=0; i<3; i++)
    {
      for(int j=0; j<3; j++)
      {
        m_particleRVectors[p][i][j] = p_F[j][0]*m_particleInitialRVectors[p][i][0] + p_F[j][1]*m_particleInitialRVectors[p][i][1] + p_F[j][2]*m_particleInitialRVectors[p][i][2];
      }
    }
  }
}

void ParticleSubRegion::cpdiDomainScaling( int p,
                                           real64 lCrit,
                                           int m_planeStrain )
{
//  for( int p=0; p<this->size(); p++ )
//  {
    arraySlice1d< real64 > r1 = m_particleRVectors[p][0];
    arraySlice1d< real64 > r2 = m_particleRVectors[p][1];
    arraySlice1d< real64 > r3 = m_particleRVectors[p][2];

    if( m_planeStrain ) // 2D cpdi domain scaling
    {
      // Initialize l-vectors.  Eq. 8a-d in the CPDI domain scaling paper.
      array2d< real64 > l;
      l.resize( 2, 3 );

      // l[0] = r1 + r2; // la
      LvArray::tensorOps::copy< 3 >( l[0], r1 );
      LvArray::tensorOps::add< 3 >( l[0], r2 );
      // l[1] = r1 - r2; // lb
      LvArray::tensorOps::copy< 3 >( l[1], r1 );
      LvArray::tensorOps::subtract< 3 >( l[1], r2 );

      // scale l-vectors if needed.  Eq. 9 in the CPDI domain scaling paper.
      bool scale = false;
      for( int i = 0 ; i < 2 ; i++ )
      {
        real64 lLength = LvArray::tensorOps::l2Norm< 3 >( l[i] );
        if( lLength > lCrit )
        {
          LvArray::tensorOps::scale< 3 >( l[i], lCrit / lLength );
          scale = true;
        }
      }

      // reconstruct r-vectors.  eq. 11 in the CPDI domain scaling paper.
      if( scale )
      {
        // r1 = 0.5 * ( l[0] + l[1] );
        LvArray::tensorOps::scaledCopy< 3 >( r1, l[0], 0.5 );
        LvArray::tensorOps::scaledAdd< 3 >( r1, l[1], 0.5 );
        // r2 = 0.5 * ( l[0] - l[1] );
        LvArray::tensorOps::scaledCopy< 3 >( r2, l[0], 0.5 );
        LvArray::tensorOps::scaledAdd< 3 >( r2, l[1], -0.5 );
      }
    }
    else // 3D cpdi domain scaling
    {
      // Initialize l-vectors.  Eq. 8a-d in the CPDI domain scaling paper.
      array2d< real64 > l;
      l.resize( 4, 3 );

      // l[0] = r1 + r2 + r3; // la
      LvArray::tensorOps::copy< 3 >( l[0], r1 );
      LvArray::tensorOps::add< 3 >( l[0], r2 );
      LvArray::tensorOps::add< 3 >( l[0], r3 );
      // l[1] = r1 - r2 + r3; // lb
      LvArray::tensorOps::copy< 3 >( l[1], r1 );
      LvArray::tensorOps::subtract< 3 >( l[1], r2 );
      LvArray::tensorOps::add< 3 >( l[1], r3 );
      // l[2] = r2 - r1 + r3; // lc
      LvArray::tensorOps::copy< 3 >( l[2], r2 );
      LvArray::tensorOps::subtract< 3 >( l[2], r1 );
      LvArray::tensorOps::add< 3 >( l[2], r3 );
      // l[3] = r3 - r1 - r2; // ld
      LvArray::tensorOps::copy< 3 >( l[3], r3 );
      LvArray::tensorOps::subtract< 3 >( l[3], r1 );
      LvArray::tensorOps::subtract< 3 >( l[3], r2 );

      // scale l vectors if needed.  Eq. 9 in the CPDI domain scaling paper.
      bool scale = false;
      for( int i = 0 ; i < 4 ; i++ )
      {
        real64 lLength = LvArray::tensorOps::l2Norm< 3 >( l[i] );
        if( lLength > lCrit )
        {
          LvArray::tensorOps::scale< 3 >( l[i], lCrit / lLength );
          scale = true;
        }
      }

      // reconstruct r vectors.  eq. 11 in the CPDI domain scaling paper.
      if( scale )
      {
        // r1 = 0.25 * ( l[0] + l[1] - l[2] - l[3] );
        LvArray::tensorOps::scaledCopy< 3 >( r1, l[0], 0.25 );
        LvArray::tensorOps::scaledAdd< 3 >( r1, l[1], 0.25 );
        LvArray::tensorOps::scaledAdd< 3 >( r1, l[2], -0.25 );
        LvArray::tensorOps::scaledAdd< 3 >( r1, l[3], -0.25 );
        // r2 = 0.25 * ( l[0] - l[1] + l[2] - l[3] );
        LvArray::tensorOps::scaledCopy< 3 >( r2, l[0], 0.25 );
        LvArray::tensorOps::scaledAdd< 3 >( r2, l[1], -0.25 );
        LvArray::tensorOps::scaledAdd< 3 >( r2, l[2], 0.25 );
        LvArray::tensorOps::scaledAdd< 3 >( r2, l[3], -0.25 );
        // r3 = 0.25 * ( l[0] + l[1] + l[2] + l[3] );
        LvArray::tensorOps::scaledCopy< 3 >( r3, l[0], 0.25 );
        LvArray::tensorOps::scaledAdd< 3 >( r3, l[1], 0.25 );
        LvArray::tensorOps::scaledAdd< 3 >( r3, l[2], 0.25 );
        LvArray::tensorOps::scaledAdd< 3 >( r3, l[3], 0.25 );
      }
    }
//  }
}

void ParticleSubRegion::getAllWeights( int const p,
                                       std::array< real64, 3 > const & xMin,
                                       std::array< real64, 3 > const & hx,
                                       array3d< int > const & ijkMap,
                                       arrayView2d< real64, nodes::REFERENCE_POSITION_USD > const & g_X,
                                       std::vector< int > & nodeIDs,
                                       std::vector< real64 > & weights,
                                       std::vector< std::vector< real64 > > & gradWeights )
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
        centerIJK[i] = std::floor( (p_x[i] - xMin[i])/hx[i] );
      }

      // get node IDs
      for(int i=0; i<2; i++)
      {
        for(int j=0; j<2; j++)
        {
          for(int k=0; k<2; k++)
          {
            nodeIDs.push_back( ijkMap[centerIJK[0]+i][centerIJK[1]+j][centerIJK[2]+k] );
          }
        }
      }

      // get weights and grad weights
      int corner = ijkMap[centerIJK[0]][centerIJK[1]][centerIJK[2]];
      auto corner_x = g_X[corner];

      real64 xRel = (p_x[0] - corner_x[0])/hx[0];
      real64 yRel = (p_x[1] - corner_x[1])/hx[1];
      real64 zRel = (p_x[2] - corner_x[2])/hx[2];

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
            weights.push_back(xWeight*yWeight*zWeight);
            gradWeights[0].push_back(dxWeight*yWeight*zWeight);
            gradWeights[1].push_back(xWeight*dyWeight*zWeight);
            gradWeights[2].push_back(xWeight*yWeight*dzWeight);
          }
        }
      }

      break;
    }

    case ParticleType::CPDI:
    {
      // Precalculated things
      int signs[8][3] = { { -1, -1, -1},
                          {  1, -1, -1},
                          {  1,  1, -1},
                          { -1,  1, -1},
                          { -1, -1,  1},
                          {  1, -1,  1},
                          {  1,  1,  1},
                          { -1,  1,  1} };

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
      std::vector<std::vector<int>> cornerIJK;
      cornerIJK.resize(8); // CPDI can map to up to 8 cells
      for(int corner=0; corner<8; corner++)
      {
        cornerIJK[corner].resize(3);
        for(int i=0; i<3; i++)
        {
          real64 cornerPositionComponent = p_x[i] + signs[corner][0] * m_particleRVectors[p][0][i] + signs[corner][1] * m_particleRVectors[p][1][i] + signs[corner][2] * m_particleRVectors[p][2][i];
          cornerIJK[corner][i] = std::floor( ( cornerPositionComponent - xMin[i] ) / hx[i] ); // TODO: Temporarily store the CPDI corners since they're re-used below?
        }
      }

      // get node IDs associated with each corner from IJK map
      // *** The order in which we access the IJK map must match the order we evaluate the shape functions! ***
      for(int corner=0; corner<8; corner++)
      {
        for(int i=0; i<2; i++)
        {
          for(int j=0; j<2; j++)
          {
            for(int k=0; k<2; k++)
            {
              nodeIDs.push_back( ijkMap[cornerIJK[corner][0]+i][cornerIJK[corner][1]+j][cornerIJK[corner][2]+k] );
            }
          }
        }
      }

      // get weights and grad weights
      for(int corner=0; corner<8; corner++)
      {
        int cornerNode = ijkMap[cornerIJK[corner][0]][cornerIJK[corner][1]][cornerIJK[corner][2]];
        auto cornerNodePosition = g_X[cornerNode];

        real64 x, y, z;
        x = p_x[0] + signs[corner][0] * m_particleRVectors[p][0][0] + signs[corner][1] * m_particleRVectors[p][1][0] + signs[corner][2] * m_particleRVectors[p][2][0];
        y = p_x[1] + signs[corner][0] * m_particleRVectors[p][0][1] + signs[corner][1] * m_particleRVectors[p][1][1] + signs[corner][2] * m_particleRVectors[p][2][1];
        z = p_x[2] + signs[corner][0] * m_particleRVectors[p][0][2] + signs[corner][1] * m_particleRVectors[p][1][2] + signs[corner][2] * m_particleRVectors[p][2][2];

        real64 xRel = ( x - cornerNodePosition[0] ) / hx[0];
        real64 yRel = ( y - cornerNodePosition[1] ) / hx[1];
        real64 zRel = ( z - cornerNodePosition[2] ) / hx[2];

        for(int i=0; i<2; i++)
        {
          real64 xWeight = i * xRel + ( 1 - i ) * ( 1.0 - xRel );
          for(int j=0; j<2; j++)
          {
            real64 yWeight = j * yRel + ( 1 - j ) * ( 1.0 - yRel );
            for(int k=0; k<2; k++)
            {
              real64 zWeight = k * zRel + ( 1 - k ) * ( 1.0 - zRel );
              real64 weight = xWeight * yWeight * zWeight;
              weights.push_back( 0.125 * weight );
              gradWeights[0].push_back( alpha[corner][0] * weight );
              gradWeights[1].push_back( alpha[corner][1] * weight );
              gradWeights[2].push_back( alpha[corner][2] * weight );
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
