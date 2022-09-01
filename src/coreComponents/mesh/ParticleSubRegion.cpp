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

void ParticleSubRegion::updateRVectors(int const p,
                                       LvArray::ArraySlice<double, 2, 1, int> const & p_F)
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
      // get cell IDs
      std::vector<int> cellID;
      cellID.resize(3);
      for(int i=0; i<3; i++)
      {
        cellID[i] = std::floor( (p_x[i] - xMin[i])/hx[i] );
      }

      // get node IDs
      for(int i=0; i<2; i++)
      {
        for(int j=0; j<2; j++)
        {
          for(int k=0; k<2; k++)
          {
            nodeIDs.push_back( ijkMap[cellID[0]+i][cellID[1]+j][cellID[2]+k] );
          }
        }
      }

      // get weights and grad weights
      int corner = ijkMap[cellID[0]][cellID[1]][cellID[2]];
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
      int signs[8][3] = { { 1,  1,  1},
                          { 1,  1, -1},
                          { 1, -1,  1},
                          { 1, -1, -1},
                          {-1,  1,  1},
                          {-1,  1, -1},
                          {-1, -1,  1},
                          {-1, -1, -1} };

      real64 alpha[8][8];
      real64 one_over_V = 1.0 / m_particleVolume[p];
      real64 p_r1[3], p_r2[3], p_r3[3]; // allowing 1-indexed r-vectors to persist to torture future postdocs >:)

      for(int i=0; i<3; i++)
      {
        p_r1[i] = m_particleRVectors[p][0][i];
        p_r2[i] = m_particleRVectors[p][1][i];
        p_r3[i] = m_particleRVectors[p][2][i];
      }

      alpha[0][0] = one_over_V * ( p_r1[2] * p_r2[1] - p_r1[1] * p_r2[2] - p_r1[2] * p_r3[1] + p_r2[2] * p_r3[1] + p_r1[1] * p_r3[2] - p_r2[1] * p_r3[2] );
      alpha[0][1] = one_over_V * ( -( p_r1[2] * p_r2[0] ) + p_r1[0] * p_r2[2] + p_r1[2] * p_r3[0] - p_r2[2] * p_r3[0] - p_r1[0] * p_r3[2] + p_r2[0] * p_r3[2] );
      alpha[0][2] = one_over_V * ( p_r1[1] * p_r2[0] - p_r1[0] * p_r2[1] - p_r1[1] * p_r3[0] + p_r2[1] * p_r3[0] + p_r1[0] * p_r3[1] - p_r2[0] * p_r3[1] );
      alpha[1][0] = one_over_V * ( p_r1[2] * p_r2[1] - p_r1[1] * p_r2[2] - p_r1[2] * p_r3[1] - p_r2[2] * p_r3[1] + p_r1[1] * p_r3[2] + p_r2[1] * p_r3[2] );
      alpha[1][1] = one_over_V * ( -( p_r1[2] * p_r2[0] ) + p_r1[0] * p_r2[2] + p_r1[2] * p_r3[0] + p_r2[2] * p_r3[0] - p_r1[0] * p_r3[2] - p_r2[0] * p_r3[2] );
      alpha[1][2] = one_over_V * ( p_r1[1] * p_r2[0] - p_r1[0] * p_r2[1] - p_r1[1] * p_r3[0] - p_r2[1] * p_r3[0] + p_r1[0] * p_r3[1] + p_r2[0] * p_r3[1] );
      alpha[2][0] = one_over_V * ( p_r1[2] * p_r2[1] - p_r1[1] * p_r2[2] + p_r1[2] * p_r3[1] - p_r2[2] * p_r3[1] - p_r1[1] * p_r3[2] + p_r2[1] * p_r3[2] );
      alpha[2][1] = one_over_V * ( -( p_r1[2] * p_r2[0] ) + p_r1[0] * p_r2[2] - p_r1[2] * p_r3[0] + p_r2[2] * p_r3[0] + p_r1[0] * p_r3[2] - p_r2[0] * p_r3[2] );
      alpha[2][2] = one_over_V * ( p_r1[1] * p_r2[0] - p_r1[0] * p_r2[1] + p_r1[1] * p_r3[0] - p_r2[1] * p_r3[0] - p_r1[0] * p_r3[1] + p_r2[0] * p_r3[1] );
      alpha[3][0] = one_over_V * ( p_r1[2] * p_r2[1] - p_r1[1] * p_r2[2] + p_r1[2] * p_r3[1] + p_r2[2] * p_r3[1] - p_r1[1] * p_r3[2] - p_r2[1] * p_r3[2] );
      alpha[3][1] = one_over_V * ( -( p_r1[2] * p_r2[0] ) + p_r1[0] * p_r2[2] - p_r1[2] * p_r3[0] - p_r2[2] * p_r3[0] + p_r1[0] * p_r3[2] + p_r2[0] * p_r3[2] );
      alpha[3][2] = one_over_V * ( p_r1[1] * p_r2[0] - p_r1[0] * p_r2[1] + p_r1[1] * p_r3[0] + p_r2[1] * p_r3[0] - p_r1[0] * p_r3[1] - p_r2[0] * p_r3[1] );
      alpha[4][0] = -alpha[2][0];
      alpha[4][1] = -alpha[2][1];
      alpha[4][2] = -alpha[2][2];
      alpha[5][0] = -alpha[3][0];
      alpha[5][1] = -alpha[3][1];
      alpha[5][2] = -alpha[3][2];
      alpha[6][0] = -alpha[0][0];
      alpha[6][1] = -alpha[0][1];
      alpha[6][2] = -alpha[0][2];
      alpha[7][0] = -alpha[1][0];
      alpha[7][1] = -alpha[1][1];
      alpha[7][2] = -alpha[1][2];

      // GEOS-to-GEOSX corner mapping, because I'm lazy
      int cornerMap[8] = {0, 4, 3, 7, 1, 5, 2, 6};

      // get cell IDs
      std::vector<std::vector<int>> cellID;
      cellID.resize(8); // CPDI can map to up to 8 cells
      for(int cell=0; cell<8; cell++)
      {
        cellID[cell].resize(3);
        for(int i=0; i<3; i++)
        {
          real64 CPDIcorner = p_x[i] + signs[cell][0]*m_particleRVectors[p][0][i] + signs[cell][1]*m_particleRVectors[p][1][i] + signs[cell][2]*m_particleRVectors[p][2][i];
          cellID[cell][i] = std::floor((CPDIcorner - xMin[i])/hx[i]); // TODO: Temporarily store the CPDI corners since they're re-used below?
        }
      }

      // get node IDs
      for(size_t cell=0; cell<cellID.size(); cell++)
      {
        for(int i=0; i<2; i++)
        {
          for(int j=0; j<2; j++)
          {
            for(int k=0; k<2; k++)
            {
              nodeIDs.push_back(ijkMap[cellID[cell][0]+i][cellID[cell][1]+j][cellID[cell][2]+k]);
            }
          }
        }
      }

      // get weights and grad weights
      for(size_t cell=0; cell<cellID.size(); cell++)
      {
        int corner = ijkMap[cellID[cell][0]][cellID[cell][1]][cellID[cell][2]];
        auto corner_x = g_X[corner];

        real64 x, y, z;
        x = p_x[0] + signs[cell][0]*m_particleRVectors[p][0][0] + signs[cell][1]*m_particleRVectors[p][1][0] + signs[cell][2]*m_particleRVectors[p][2][0];
        y = p_x[1] + signs[cell][0]*m_particleRVectors[p][0][1] + signs[cell][1]*m_particleRVectors[p][1][1] + signs[cell][2]*m_particleRVectors[p][2][1];
        z = p_x[2] + signs[cell][0]*m_particleRVectors[p][0][2] + signs[cell][1]*m_particleRVectors[p][1][2] + signs[cell][2]*m_particleRVectors[p][2][2];

        real64 xRel = (x - corner_x[0])/hx[0];
        real64 yRel = (y - corner_x[1])/hx[1];
        real64 zRel = (z - corner_x[2])/hx[2];

        for(int i=0; i<2; i++)
        {
          real64 xWeight = i*xRel + (1-i)*(1.0-xRel);
          for(int j=0; j<2; j++)
          {
            real64 yWeight = j*yRel + (1-j)*(1.0-yRel);
            for(int k=0; k<2; k++)
            {
              real64 zWeight = k*zRel + (1-k)*(1.0-zRel);
              real64 weight = xWeight*yWeight*zWeight;
              weights.push_back(0.125*weight); // note the built-in factor of 1/8 so we don't need it in the solver
              gradWeights[0].push_back(-alpha[cornerMap[cell]][0]*weight);
              gradWeights[1].push_back(-alpha[cornerMap[cell]][1]*weight);
              gradWeights[2].push_back(-alpha[cornerMap[cell]][2]*weight);
            }
          }
        }
      }

      break;
    }

    default:
    {
      GEOSX_ERROR( "Invalid particle type: " << m_particleType );
    }
  }

}

REGISTER_CATALOG_ENTRY( ObjectManagerBase, ParticleSubRegion, string const &, Group * const )

} /* namespace geosx */
