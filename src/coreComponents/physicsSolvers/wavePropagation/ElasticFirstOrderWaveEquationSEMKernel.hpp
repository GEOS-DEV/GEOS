/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOS Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file ElasticFirstOrderWaveEquationSEMKernel.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ELASTICFIRSTORDERWAVEEQUATIONSEMKERNEL_HPP_
#define GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ELASTICFIRSTORDERWAVEEQUATIONSEMKERNEL_HPP_

#include "finiteElement/kernelInterface/KernelBase.hpp"
#include "WaveSolverUtils.hpp"


namespace geos
{

/// Namespace to contain the elastic wave kernels.
namespace elasticFirstOrderWaveEquationSEMKernels
{

template< typename FE_TYPE >
struct ComputeTimeStep
{

  ComputeTimeStep( FE_TYPE const & finiteElement )
    : m_finiteElement( finiteElement )
  {}

  /**
   * @brief Compute timestep using power iteration method
   */
  template< typename EXEC_POLICY, typename ATOMIC_POLICY >
  real64
  launch( localIndex const sizeElem,
          localIndex const sizeNode,
          arrayView2d< WaveSolverBase::wsCoordType const, nodes::REFERENCE_POSITION_USD > const nodeCoords,
          arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes,
          arrayView1d< real32 const > const pWavespeed,
          arrayView1d< real32 const > const sWavespeed,
          arrayView1d< real32 const > const density,
          arrayView1d< real32 > const lambda,
          arrayView1d< real32 > const mu,
          arrayView1d< real32 > const mass )
  {

    constexpr localIndex numNodesPerElem = FE_TYPE::numNodes;
    constexpr localIndex numQuadraturePointsPerElem = FE_TYPE::numQuadraturePoints;

    real64 const epsilon = 0.00001;
    localIndex const nIterMax = 10000;
    localIndex numberIter = 0;
    localIndex counter = 0;
    real64 lambdaNew = 0.0;

    array1d< real32 > const ux( sizeNode );
    array1d< real32 > const uxAux( sizeNode );
    array1d< real32 > const uy( sizeNode );
    array1d< real32 > const uyAux( sizeNode );
    array1d< real32 > const uz( sizeNode );
    array1d< real32 > const uzAux( sizeNode );
    array2d< real32 > const stressxx( sizeElem, numNodesPerElem );
    array2d< real32 > const stressyy( sizeElem, numNodesPerElem );
    array2d< real32 > const stresszz( sizeElem, numNodesPerElem );
    array2d< real32 > const stressxy( sizeElem, numNodesPerElem );
    array2d< real32 > const stressxz( sizeElem, numNodesPerElem );
    array2d< real32 > const stressyz( sizeElem, numNodesPerElem );

    arrayView1d< real32 > const uxView = ux;
    arrayView1d< real32 > const uyView = uy;
    arrayView1d< real32 > const uzView = uz;
    arrayView1d< real32 > const uxAuxView = uxAux;
    arrayView1d< real32 > const uyAuxView = uyAux;
    arrayView1d< real32 > const uzAuxView = uzAux;
    arrayView2d< real32 > const stressxxView = stressxx;
    arrayView2d< real32 > const stressyyView = stressyy;
    arrayView2d< real32 > const stresszzView = stresszz;
    arrayView2d< real32 > const stressxyView = stressxy;
    arrayView2d< real32 > const stressxzView = stressxz;
    arrayView2d< real32 > const stressyzView = stressyz;


    //Randomize u values
    srand( time( NULL ));
    for( localIndex a = 0; a < sizeNode; ++a )
    {
      uxView[a] = (real64)rand()/(real64) RAND_MAX;
      uyView[a] = (real64)rand()/(real64) RAND_MAX;
      uzView[a] = (real64)rand()/(real64) RAND_MAX;
    }

    //Step 1: Normalize randomized pressure
    real64 normUx= 0.0;
    real64 normUy= 0.0;
    real64 normUz= 0.0;
    WaveSolverUtils::dotProduct( sizeNode, uxView, uxView, normUx );
    WaveSolverUtils::dotProduct( sizeNode, uyView, uyView, normUy );
    WaveSolverUtils::dotProduct( sizeNode, uzView, uzView, normUz );
    real64 normUtot = normUx+normUy+normUz;


    forAll< EXEC_POLICY >( sizeNode, [=] GEOS_HOST_DEVICE ( localIndex const a )
    {
      uxView[a]/= sqrt( normUtot );
      uyView[a]/= sqrt( normUtot );
      uzView[a]/= sqrt( normUtot );
    } );

    //Step 2: Initial iteration of (M^{-1}K)p
    forAll< EXEC_POLICY >( sizeElem, [=] GEOS_HOST_DEVICE ( localIndex const k )
    {
      real64 xLocal[8][3];
      for( localIndex a=0; a< 8; ++a )
      {
        localIndex const nodeIndex = elemsToNodes( k, FE_TYPE::meshIndexToLinearIndex3D( a ) );
        for( localIndex i=0; i<3; ++i )
        {
          xLocal[a][i] = nodeCoords[ nodeIndex ][ i ];
        }
      }

      mu[k] = density[k] * sWavespeed[k] * sWavespeed[k];
      lambda[k] = density[k] * pWavespeed[k] * pWavespeed[k] - 2.0*mu[k];

      real32 uelemxx[numNodesPerElem] = {0.0};
      real32 uelemyy[numNodesPerElem] = {0.0};
      real32 uelemzz[numNodesPerElem] = {0.0};
      real32 uelemxy[numNodesPerElem] = {0.0};
      real32 uelemxz[numNodesPerElem] = {0.0};
      real32 uelemyz[numNodesPerElem]= {0.0};
      real32 auxx[numNodesPerElem] = {0.0};
      real32 auyy[numNodesPerElem] = {0.0};
      real32 auzz[numNodesPerElem] = {0.0};
      real32 auxy[numNodesPerElem] = {0.0};
      real32 auxz[numNodesPerElem] = {0.0};
      real32 auyz[numNodesPerElem] = {0.0};

      for( localIndex q=0; q<numQuadraturePointsPerElem; ++q )
      {

        //Volume integral
        m_finiteElement.computeFirstOrderStiffnessTermX( q, xLocal, [&] ( int i, int j, real32 dfx1, real32 dfx2, real32 dfx3 )
        {
          auxx[j]+= dfx1*uxView[elemsToNodes[k][i]];
          auyy[j]+= dfx2*uyView[elemsToNodes[k][i]];
          auzz[j]+= dfx3*uzView[elemsToNodes[k][i]];
          auxy[j]+= dfx1*uyView[elemsToNodes[k][i]]+dfx2*uxView[elemsToNodes[k][i]];
          auxz[j]+= dfx1*uzView[elemsToNodes[k][i]]+dfx3*uxView[elemsToNodes[k][i]];
          auyz[j]+= dfx2*uzView[elemsToNodes[k][i]]+dfx3*uyView[elemsToNodes[k][i]];

        } );

        m_finiteElement.computeFirstOrderStiffnessTermY( q, xLocal, [&] ( int i, int j, real32 dfy1, real32 dfy2, real32 dfy3 )
        {
          auxx[j]+= dfy1*uxView[elemsToNodes[k][i]];
          auyy[j]+= dfy2*uyView[elemsToNodes[k][i]];
          auzz[j]+= dfy3*uzView[elemsToNodes[k][i]];
          auxy[j]+= dfy1*uyView[elemsToNodes[k][i]]+dfy2*uxView[elemsToNodes[k][i]];
          auxz[j]+= dfy1*uzView[elemsToNodes[k][i]]+dfy3*uxView[elemsToNodes[k][i]];
          auyz[j]+= dfy2*uzView[elemsToNodes[k][i]]+dfy3*uyView[elemsToNodes[k][i]];

        } );

        m_finiteElement.computeFirstOrderStiffnessTermZ( q, xLocal, [&] ( int i, int j, real32 dfz1, real32 dfz2, real32 dfz3 )
        {
          auxx[j]+= dfz1*uxView[elemsToNodes[k][i]];
          auyy[j]+= dfz2*uyView[elemsToNodes[k][i]];
          auzz[j]+= dfz3*uzView[elemsToNodes[k][i]];
          auxy[j]+= dfz1*uyView[elemsToNodes[k][i]]+dfz2*uxView[elemsToNodes[k][i]];
          auxz[j]+= dfz1*uzView[elemsToNodes[k][i]]+dfz3*uxView[elemsToNodes[k][i]];
          auyz[j]+= dfz2*uzView[elemsToNodes[k][i]]+dfz3*uyView[elemsToNodes[k][i]];

        } );

      }
      //Time integration
      for( localIndex i = 0; i < numNodesPerElem; ++i )
      {
        real32 diag = lambda[k]*(auxx[i]+auyy[i]+auzz[i]);
        uelemxx[i]+= (diag+2*mu[k]*auxx[i]);
        uelemyy[i]+= (diag+2*mu[k]*auyy[i]);
        uelemzz[i]+= (diag+2*mu[k]*auzz[i]);
        uelemxy[i]+= mu[k]*auxy[i];
        uelemxz[i]+= mu[k]*auxz[i];
        uelemyz[i]+= mu[k]*auyz[i];
      }

      // Multiplication by inverse mass matrix
      for( localIndex i = 0; i < numNodesPerElem; ++i )
      {
        real32 massLoc = m_finiteElement.computeMassTerm( i, xLocal );
        stressxxView[k][i] = uelemxx[i]/massLoc;
        stressyyView[k][i] = uelemyy[i]/massLoc;
        stresszzView[k][i] = uelemzz[i]/massLoc;
        stressxyView[k][i] = uelemxy[i]/massLoc;
        stressxzView[k][i] = uelemxz[i]/massLoc;
        stressyzView[k][i] = uelemyz[i]/massLoc;
      }

    } );

    uxAuxView.zero();
    uyAuxView.zero();
    uzAuxView.zero();

    forAll< EXEC_POLICY >( sizeElem, [=] GEOS_HOST_DEVICE ( localIndex const k )
    {

      real64 xLocal[8][3];
      for( localIndex a=0; a< 8; ++a )
      {
        localIndex const nodeIndex = elemsToNodes( k, FE_TYPE::meshIndexToLinearIndex3D( a ) );
        for( localIndex i=0; i<3; ++i )
        {
          xLocal[a][i] = nodeCoords[ nodeIndex ][ i ];
        }
      }

      real32 flowx[numNodesPerElem] = {0.0};
      real32 flowy[numNodesPerElem] = {0.0};
      real32 flowz[numNodesPerElem] = {0.0};

      for( localIndex q=0; q<numQuadraturePointsPerElem; ++q )
      {


        // Stiffness part
        m_finiteElement.computeFirstOrderStiffnessTermX( q, xLocal, [&] ( int i, int j, real32 dfx1, real32 dfx2, real32 dfx3 )
        {
          flowx[i] -= stressxxView[k][j]*dfx1 + stressxyView[k][j]*dfx2 + stressxzView[k][j]*dfx3;
          flowy[i] -= stressxyView[k][j]*dfx1 + stressyyView[k][j]*dfx2 + stressyzView[k][j]*dfx3;
          flowz[i] -= stressxzView[k][j]*dfx1 + stressyzView[k][j]*dfx2 + stresszzView[k][j]*dfx3;

        } );

        m_finiteElement.computeFirstOrderStiffnessTermY( q, xLocal, [&] ( int i, int j, real32 dfy1, real32 dfy2, real32 dfy3 )
        {
          flowx[i] -= stressxxView[k][j]*dfy1 + stressxyView[k][j]*dfy2 + stressxzView[k][j]*dfy3;
          flowy[i] -= stressxyView[k][j]*dfy1 + stressyyView[k][j]*dfy2 + stressyzView[k][j]*dfy3;
          flowz[i] -= stressxzView[k][j]*dfy1 + stressyzView[k][j]*dfy2 + stresszzView[k][j]*dfy3;
        } );

        m_finiteElement.computeFirstOrderStiffnessTermZ( q, xLocal, [&] ( int i, int j, real32 dfz1, real32 dfz2, real32 dfz3 )
        {
          flowx[i] -= stressxxView[k][j]*dfz1 + stressxyView[k][j]*dfz2 + stressxzView[k][j]*dfz3;
          flowy[i] -= stressxyView[k][j]*dfz1 + stressyyView[k][j]*dfz2 + stressyzView[k][j]*dfz3;
          flowz[i] -= stressxzView[k][j]*dfz1 + stressyzView[k][j]*dfz2 + stresszzView[k][j]*dfz3;
        } );

      }

      // Mult by inverse mass matrix + damping matrix
      for( localIndex i = 0; i < numNodesPerElem; ++i )
      {
        real32 localIncrement1 = flowx[i]/mass[elemsToNodes[k][i]];
        real32 localIncrement2 = flowy[i]/mass[elemsToNodes[k][i]];
        real32 localIncrement3 = flowz[i]/mass[elemsToNodes[k][i]];
        RAJA::atomicAdd< ATOMIC_POLICY >( &uxAuxView[elemsToNodes[k][i]], localIncrement1 );
        RAJA::atomicAdd< ATOMIC_POLICY >( &uyAuxView[elemsToNodes[k][i]], localIncrement2 );
        RAJA::atomicAdd< ATOMIC_POLICY >( &uzAuxView[elemsToNodes[k][i]], localIncrement3 );
      }

    } );

    real64 lambdaOld = lambdaNew;

    //Compute lambdaNew using two dotProducts
    real64 dotProductUxUxaux = 0.0;
    real64 dotProductUyUyaux = 0.0;
    real64 dotProductUzUzaux = 0.0;

    WaveSolverUtils::dotProduct( sizeNode, uxView, uxAuxView, dotProductUxUxaux );
    WaveSolverUtils::dotProduct( sizeNode, uyView, uyAuxView, dotProductUyUyaux );
    WaveSolverUtils::dotProduct( sizeNode, uzView, uzAuxView, dotProductUzUzaux );
    real64 dotProductUtotUtotAux = dotProductUxUxaux+dotProductUyUyaux+dotProductUzUzaux;

    normUx = 0.0;
    normUy = 0.0;
    normUz = 0.0;

    WaveSolverUtils::dotProduct( sizeNode, uxView, uxView, normUx );
    WaveSolverUtils::dotProduct( sizeNode, uyView, uyView, normUy );
    WaveSolverUtils::dotProduct( sizeNode, uzView, uzView, normUz );
    normUtot = normUx+normUy+normUz;


    lambdaNew = dotProductUtotUtotAux/normUtot;

    real64 normUxaux = 0.0;
    real64 normUyaux = 0.0;
    real64 normUzaux = 0.0;
    WaveSolverUtils::dotProduct( sizeNode, uxAuxView, uxAuxView, normUxaux );
    WaveSolverUtils::dotProduct( sizeNode, uyAuxView, uyAuxView, normUyaux );
    WaveSolverUtils::dotProduct( sizeNode, uzAuxView, uzAuxView, normUzaux );

    real64 normUtotAux = normUxaux+normUyaux+normUzaux;


    forAll< EXEC_POLICY >( sizeNode, [=] GEOS_HOST_DEVICE ( localIndex const a )
    {
      uxView[a]= uxAuxView[a]/sqrt( normUtotAux );
      uyView[a]= uyAuxView[a]/sqrt( normUtotAux );
      uzView[a]= uzAuxView[a]/sqrt( normUtotAux );
    } );

    //Step 3: Do previous algorithm until we found the max eigenvalues
    do
    {

      forAll< EXEC_POLICY >( sizeElem, [=] GEOS_HOST_DEVICE ( localIndex const k )
      {
        real64 xLocal[8][3];
        for( localIndex a=0; a< 8; ++a )
        {
          localIndex const nodeIndex = elemsToNodes( k, FE_TYPE::meshIndexToLinearIndex3D( a ) );
          for( localIndex i=0; i<3; ++i )
          {
            xLocal[a][i] = nodeCoords[ nodeIndex ][ i ];
          }
        }

        mu[k] = density[k] * sWavespeed[k] * sWavespeed[k];
        lambda[k] = density[k] * pWavespeed[k] * pWavespeed[k] - 2.0*mu[k];

        real32 uelemxx[numNodesPerElem] = {0.0};
        real32 uelemyy[numNodesPerElem] = {0.0};
        real32 uelemzz[numNodesPerElem] = {0.0};
        real32 uelemxy[numNodesPerElem] = {0.0};
        real32 uelemxz[numNodesPerElem] = {0.0};
        real32 uelemyz[numNodesPerElem]= {0.0};
        real32 auxx[numNodesPerElem] = {0.0};
        real32 auyy[numNodesPerElem] = {0.0};
        real32 auzz[numNodesPerElem] = {0.0};
        real32 auxy[numNodesPerElem] = {0.0};
        real32 auxz[numNodesPerElem] = {0.0};
        real32 auyz[numNodesPerElem] = {0.0};

        for( localIndex q=0; q<numQuadraturePointsPerElem; ++q )
        {

          //Volume integral
          m_finiteElement.computeFirstOrderStiffnessTermX( q, xLocal, [&] ( int i, int j, real32 dfx1, real32 dfx2, real32 dfx3 )
          {
            auxx[j]+= dfx1*uxView[elemsToNodes[k][i]];
            auyy[j]+= dfx2*uyView[elemsToNodes[k][i]];
            auzz[j]+= dfx3*uzView[elemsToNodes[k][i]];
            auxy[j]+= dfx1*uyView[elemsToNodes[k][i]]+dfx2*uxView[elemsToNodes[k][i]];
            auxz[j]+= dfx1*uzView[elemsToNodes[k][i]]+dfx3*uxView[elemsToNodes[k][i]];
            auyz[j]+= dfx2*uzView[elemsToNodes[k][i]]+dfx3*uyView[elemsToNodes[k][i]];

          } );

          m_finiteElement.computeFirstOrderStiffnessTermY( q, xLocal, [&] ( int i, int j, real32 dfy1, real32 dfy2, real32 dfy3 )
          {
            auxx[j]+= dfy1*uxView[elemsToNodes[k][i]];
            auyy[j]+= dfy2*uyView[elemsToNodes[k][i]];
            auzz[j]+= dfy3*uzView[elemsToNodes[k][i]];
            auxy[j]+= dfy1*uyView[elemsToNodes[k][i]]+dfy2*uxView[elemsToNodes[k][i]];
            auxz[j]+= dfy1*uzView[elemsToNodes[k][i]]+dfy3*uxView[elemsToNodes[k][i]];
            auyz[j]+= dfy2*uzView[elemsToNodes[k][i]]+dfy3*uyView[elemsToNodes[k][i]];

          } );

          m_finiteElement.computeFirstOrderStiffnessTermZ( q, xLocal, [&] ( int i, int j, real32 dfz1, real32 dfz2, real32 dfz3 )
          {
            auxx[j]+= dfz1*uxView[elemsToNodes[k][i]];
            auyy[j]+= dfz2*uyView[elemsToNodes[k][i]];
            auzz[j]+= dfz3*uzView[elemsToNodes[k][i]];
            auxy[j]+= dfz1*uyView[elemsToNodes[k][i]]+dfz2*uxView[elemsToNodes[k][i]];
            auxz[j]+= dfz1*uzView[elemsToNodes[k][i]]+dfz3*uxView[elemsToNodes[k][i]];
            auyz[j]+= dfz2*uzView[elemsToNodes[k][i]]+dfz3*uyView[elemsToNodes[k][i]];

          } );

        }
        //Time integration
        for( localIndex i = 0; i < numNodesPerElem; ++i )
        {
          real32 diag = lambda[k]*(auxx[i]+auyy[i]+auzz[i]);
          uelemxx[i]+= (diag+2*mu[k]*auxx[i]);
          uelemyy[i]+= (diag+2*mu[k]*auyy[i]);
          uelemzz[i]+= (diag+2*mu[k]*auzz[i]);
          uelemxy[i]+= mu[k]*auxy[i];
          uelemxz[i]+= mu[k]*auxz[i];
          uelemyz[i]+= mu[k]*auyz[i];
        }

        // Multiplication by inverse mass matrix
        for( localIndex i = 0; i < numNodesPerElem; ++i )
        {
          real32 massLoc = m_finiteElement.computeMassTerm( i, xLocal );
          stressxxView[k][i] = uelemxx[i]/massLoc;
          stressyyView[k][i] = uelemyy[i]/massLoc;
          stresszzView[k][i] = uelemzz[i]/massLoc;
          stressxyView[k][i] = uelemxy[i]/massLoc;
          stressxzView[k][i] = uelemxz[i]/massLoc;
          stressyzView[k][i] = uelemyz[i]/massLoc;
        }

      } );

      uxAuxView.zero();
      uyAuxView.zero();
      uzAuxView.zero();

      forAll< EXEC_POLICY >( sizeElem, [=] GEOS_HOST_DEVICE ( localIndex const k )
      {

        real64 xLocal[8][3];
        for( localIndex a=0; a< 8; ++a )
        {
          localIndex const nodeIndex = elemsToNodes( k, FE_TYPE::meshIndexToLinearIndex3D( a ) );
          for( localIndex i=0; i<3; ++i )
          {
            xLocal[a][i] = nodeCoords[ nodeIndex ][ i ];
          }
        }

        real32 flowx[numNodesPerElem] = {0.0};
        real32 flowy[numNodesPerElem] = {0.0};
        real32 flowz[numNodesPerElem] = {0.0};

        for( localIndex q=0; q<numQuadraturePointsPerElem; ++q )
        {


          // Stiffness part
          m_finiteElement.computeFirstOrderStiffnessTermX( q, xLocal, [&] ( int i, int j, real32 dfx1, real32 dfx2, real32 dfx3 )
          {
            flowx[i] -= stressxxView[k][j]*dfx1 + stressxyView[k][j]*dfx2 + stressxzView[k][j]*dfx3;
            flowy[i] -= stressxyView[k][j]*dfx1 + stressyyView[k][j]*dfx2 + stressyzView[k][j]*dfx3;
            flowz[i] -= stressxzView[k][j]*dfx1 + stressyzView[k][j]*dfx2 + stresszzView[k][j]*dfx3;

          } );

          m_finiteElement.computeFirstOrderStiffnessTermY( q, xLocal, [&] ( int i, int j, real32 dfy1, real32 dfy2, real32 dfy3 )
          {
            flowx[i] -= stressxxView[k][j]*dfy1 + stressxyView[k][j]*dfy2 + stressxzView[k][j]*dfy3;
            flowy[i] -= stressxyView[k][j]*dfy1 + stressyyView[k][j]*dfy2 + stressyzView[k][j]*dfy3;
            flowz[i] -= stressxzView[k][j]*dfy1 + stressyzView[k][j]*dfy2 + stresszzView[k][j]*dfy3;
          } );

          m_finiteElement.computeFirstOrderStiffnessTermZ( q, xLocal, [&] ( int i, int j, real32 dfz1, real32 dfz2, real32 dfz3 )
          {
            flowx[i] -= stressxxView[k][j]*dfz1 + stressxyView[k][j]*dfz2 + stressxzView[k][j]*dfz3;
            flowy[i] -= stressxyView[k][j]*dfz1 + stressyyView[k][j]*dfz2 + stressyzView[k][j]*dfz3;
            flowz[i] -= stressxzView[k][j]*dfz1 + stressyzView[k][j]*dfz2 + stresszzView[k][j]*dfz3;
          } );

        }

        // Mult by inverse mass matrix + damping matrix
        for( localIndex i = 0; i < numNodesPerElem; ++i )
        {
          real32 localIncrement1 = flowx[i]/mass[elemsToNodes[k][i]];
          real32 localIncrement2 = flowy[i]/mass[elemsToNodes[k][i]];
          real32 localIncrement3 = flowz[i]/mass[elemsToNodes[k][i]];
          RAJA::atomicAdd< ATOMIC_POLICY >( &uxAuxView[elemsToNodes[k][i]], localIncrement1 );
          RAJA::atomicAdd< ATOMIC_POLICY >( &uyAuxView[elemsToNodes[k][i]], localIncrement2 );
          RAJA::atomicAdd< ATOMIC_POLICY >( &uzAuxView[elemsToNodes[k][i]], localIncrement3 );
        }

      } );


      lambdaOld = lambdaNew;

      //Compute lambdaNew using two dotProducts
      dotProductUxUxaux = 0.0;
      dotProductUyUyaux = 0.0;
      dotProductUzUzaux = 0.0;

      WaveSolverUtils::dotProduct( sizeNode, uxView, uxAuxView, dotProductUxUxaux );
      WaveSolverUtils::dotProduct( sizeNode, uyView, uyAuxView, dotProductUyUyaux );
      WaveSolverUtils::dotProduct( sizeNode, uzView, uzAuxView, dotProductUzUzaux );
      dotProductUtotUtotAux = dotProductUxUxaux+dotProductUyUyaux+dotProductUzUzaux;

      normUx = 0.0;
      normUy = 0.0;
      normUz = 0.0;

      WaveSolverUtils::dotProduct( sizeNode, uxView, uxView, normUx );
      WaveSolverUtils::dotProduct( sizeNode, uyView, uyView, normUy );
      WaveSolverUtils::dotProduct( sizeNode, uzView, uzView, normUz );
      normUtot = normUx+normUy+normUz;


      lambdaNew = dotProductUtotUtotAux/normUtot;

      //lambdaNew = LvArray::tensorOps::AiBi<sizeNode>(p,pAux)/LvArray::tensorOps::AiBi<sizeNode>(pAux,pAux);

      normUxaux = 0.0;
      normUyaux = 0.0;
      normUzaux = 0.0;
      WaveSolverUtils::dotProduct( sizeNode, uxAuxView, uxAuxView, normUxaux );
      WaveSolverUtils::dotProduct( sizeNode, uyAuxView, uyAuxView, normUyaux );
      WaveSolverUtils::dotProduct( sizeNode, uzAuxView, uzAuxView, normUzaux );

      normUtotAux = normUxaux+normUyaux+normUzaux;


      forAll< EXEC_POLICY >( sizeNode, [=] GEOS_HOST_DEVICE ( localIndex const a )
      {
        uxView[a]= uxAuxView[a]/sqrt( normUtotAux );
        uyView[a]= uyAuxView[a]/sqrt( normUtotAux );
        uzView[a]= uzAuxView[a]/sqrt( normUtotAux );
      } );

      if( LvArray::math::abs( lambdaNew-lambdaOld )/LvArray::math::abs( lambdaNew )<= epsilon )
      {
        counter++;
      }
      else
      {
        counter=0;
      }

      numberIter++;


    }
    while (counter < 10 && numberIter < nIterMax);

    GEOS_THROW_IF( numberIter> nIterMax, "Power Iteration algorithm does not converge", std::runtime_error );

    real64 dt = 1.99/sqrt( LvArray::math::abs( lambdaNew ));

    return dt;

  }

  /// The finite element space/discretization object for the element type in the subRegion
  FE_TYPE const & m_finiteElement;
};

template< typename FE_TYPE >
struct StressComputation
{

  StressComputation( FE_TYPE const & finiteElement )
    : m_finiteElement( finiteElement )
  {}

  /**
   * Add comments
   */

  template< typename EXEC_POLICY, typename ATOMIC_POLICY >
  void
  launch( localIndex const size,
          localIndex const regionIndex,
          arrayView2d< WaveSolverBase::wsCoordType const, nodes::REFERENCE_POSITION_USD > const nodeCoords,
          arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodes,
          arrayView1d< real32 const > const ux_np1,
          arrayView1d< real32 const > const uy_np1,
          arrayView1d< real32 const > const uz_np1,
          arrayView1d< real32 const > const density,
          arrayView1d< real32 const > const velocityVp,
          arrayView1d< real32 const > const velocityVs,
          arrayView1d< real32 > const lambda,
          arrayView1d< real32 > const mu,
          arrayView2d< real64 const > const sourceConstants,
          arrayView1d< localIndex const > const sourceIsLocal,
          arrayView1d< localIndex const > const sourceElem,
          arrayView1d< localIndex const > const sourceRegion,
          arrayView2d< real32 const > const sourceValue,
          real64 const dt,
          integer const cycleNumber,
          arrayView2d< real32 > const stressxx,
          arrayView2d< real32 > const stressyy,
          arrayView2d< real32 > const stresszz,
          arrayView2d< real32 > const stressxy,
          arrayView2d< real32 > const stressxz,
          arrayView2d< real32 > const stressyz )

  {
    forAll< EXEC_POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const k )
    {
      // only the eight corners of the mesh cell are needed to compute the Jacobian
      real64 xLocal[8][3];
      for( localIndex a=0; a< 8; ++a )
      {
        localIndex const nodeIndex = elemsToNodes( k, FE_TYPE::meshIndexToLinearIndex3D( a ) );
        for( localIndex i=0; i<3; ++i )
        {
          xLocal[a][i] = nodeCoords( nodeIndex, i );
        }
      }

      mu[k] = density[k] * pow( velocityVs[k], 2 );
      lambda[k] = density[k] * pow( velocityVp[k], 2 ) - 2.0*mu[k];

      constexpr localIndex numNodesPerElem = FE_TYPE::numNodes;
      real32 uelemxx[numNodesPerElem] = {0.0};
      real32 uelemyy[numNodesPerElem] = {0.0};
      real32 uelemzz[numNodesPerElem] = {0.0};
      real32 uelemxy[numNodesPerElem] = {0.0};
      real32 uelemxz[numNodesPerElem] = {0.0};
      real32 uelemyz[numNodesPerElem]= {0.0};
      real32 auxx[numNodesPerElem] = {0.0};
      real32 auyy[numNodesPerElem] = {0.0};
      real32 auzz[numNodesPerElem] = {0.0};
      real32 auxy[numNodesPerElem] = {0.0};
      real32 auxz[numNodesPerElem] = {0.0};
      real32 auyz[numNodesPerElem] = {0.0};


      //Pre-multiplication by mass matrix
      for( localIndex i = 0; i < numNodesPerElem; ++i )
      {
        real32 massLoc = m_finiteElement.computeMassTerm( i, xLocal );
        uelemxx[i] = massLoc*stressxx[k][i];
        uelemyy[i] = massLoc*stressyy[k][i];
        uelemzz[i] = massLoc*stresszz[k][i];
        uelemxy[i] = massLoc*stressxy[k][i];
        uelemxz[i] = massLoc*stressxz[k][i];
        uelemyz[i] = massLoc*stressyz[k][i];
      }

      for( localIndex q = 0; q < numNodesPerElem; q++ )
      {

        //Volume integral
        m_finiteElement.template computeFirstOrderStiffnessTermX( q, xLocal, [&] ( int i, int j, real32 dfx1, real32 dfx2, real32 dfx3 )
        {
          auxx[j]+= dfx1*ux_np1[elemsToNodes[k][i]];
          auyy[j]+= dfx2*uy_np1[elemsToNodes[k][i]];
          auzz[j]+= dfx3*uz_np1[elemsToNodes[k][i]];
          auxy[j]+= dfx1*uy_np1[elemsToNodes[k][i]]+dfx2*ux_np1[elemsToNodes[k][i]];
          auxz[j]+= dfx1*uz_np1[elemsToNodes[k][i]]+dfx3*ux_np1[elemsToNodes[k][i]];
          auyz[j]+= dfx2*uz_np1[elemsToNodes[k][i]]+dfx3*uy_np1[elemsToNodes[k][i]];

        } );

        m_finiteElement.template computeFirstOrderStiffnessTermY( q, xLocal, [&] ( int i, int j, real32 dfy1, real32 dfy2, real32 dfy3 )
        {
          auxx[j]+= dfy1*ux_np1[elemsToNodes[k][i]];
          auyy[j]+= dfy2*uy_np1[elemsToNodes[k][i]];
          auzz[j]+= dfy3*uz_np1[elemsToNodes[k][i]];
          auxy[j]+= dfy1*uy_np1[elemsToNodes[k][i]]+dfy2*ux_np1[elemsToNodes[k][i]];
          auxz[j]+= dfy1*uz_np1[elemsToNodes[k][i]]+dfy3*ux_np1[elemsToNodes[k][i]];
          auyz[j]+= dfy2*uz_np1[elemsToNodes[k][i]]+dfy3*uy_np1[elemsToNodes[k][i]];

        } );

        m_finiteElement.template computeFirstOrderStiffnessTermZ( q, xLocal, [&] ( int i, int j, real32 dfz1, real32 dfz2, real32 dfz3 )
        {
          auxx[j]+= dfz1*ux_np1[elemsToNodes[k][i]];
          auyy[j]+= dfz2*uy_np1[elemsToNodes[k][i]];
          auzz[j]+= dfz3*uz_np1[elemsToNodes[k][i]];
          auxy[j]+= dfz1*uy_np1[elemsToNodes[k][i]]+dfz2*ux_np1[elemsToNodes[k][i]];
          auxz[j]+= dfz1*uz_np1[elemsToNodes[k][i]]+dfz3*ux_np1[elemsToNodes[k][i]];
          auyz[j]+= dfz2*uz_np1[elemsToNodes[k][i]]+dfz3*uy_np1[elemsToNodes[k][i]];

        } );

      }
      //Time integration
      for( localIndex i = 0; i < numNodesPerElem; ++i )
      {
        real32 diag = lambda[k]*(auxx[i]+auyy[i]+auzz[i]);
        uelemxx[i]+= dt*(diag+2*mu[k]*auxx[i]);
        uelemyy[i]+= dt*(diag+2*mu[k]*auyy[i]);
        uelemzz[i]+= dt*(diag+2*mu[k]*auzz[i]);
        uelemxy[i]+= dt*mu[k]*auxy[i];
        uelemxz[i]+= dt*mu[k]*auxz[i];
        uelemyz[i]+= dt*mu[k]*auyz[i];
      }

      // Multiplication by inverse mass matrix
      for( localIndex i = 0; i < numNodesPerElem; ++i )
      {
        real32 massLoc = m_finiteElement.computeMassTerm( i, xLocal );
        stressxx[k][i] = uelemxx[i]/massLoc;
        stressyy[k][i] = uelemyy[i]/massLoc;
        stresszz[k][i] = uelemzz[i]/massLoc;
        stressxy[k][i] = uelemxy[i]/massLoc;
        stressxz[k][i] = uelemxz[i]/massLoc;
        stressyz[k][i] = uelemyz[i]/massLoc;
      }


      //Source injection
      for( localIndex isrc = 0; isrc < sourceConstants.size( 0 ); ++isrc )
      {
        if( sourceIsLocal[isrc] == 1 )
        {
          if( sourceElem[isrc]==k && sourceRegion[isrc] == regionIndex )
          {
            for( localIndex i = 0; i < numNodesPerElem; ++i )
            {
              real32 massLoc = m_finiteElement.computeMassTerm( i, xLocal );
              real32 const localIncrement = dt*(sourceConstants[isrc][i]*sourceValue[cycleNumber][isrc])/massLoc;
              RAJA::atomicAdd< ATOMIC_POLICY >( &stressxx[k][i], localIncrement );
              RAJA::atomicAdd< ATOMIC_POLICY >( &stressyy[k][i], localIncrement );
              RAJA::atomicAdd< ATOMIC_POLICY >( &stresszz[k][i], localIncrement );
            }

          }

        }
      }

    } );
  }

  /// The finite element space/discretization object for the element type in the subRegion

  FE_TYPE const & m_finiteElement;
};

template< typename FE_TYPE >
struct VelocityComputation
{

  VelocityComputation( FE_TYPE const & finiteElement )
    : m_finiteElement( finiteElement )
  {}

  /**
   * add doc
   */

  template< typename EXEC_POLICY, typename ATOMIC_POLICY >
  void
  launch( localIndex const size,
          localIndex const size_node,
          arrayView2d< WaveSolverBase::wsCoordType const, nodes::REFERENCE_POSITION_USD > const nodeCoords,
          arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodes,
          arrayView2d< real32 const > const stressxx,
          arrayView2d< real32 const > const stressyy,
          arrayView2d< real32 const > const stresszz,
          arrayView2d< real32 const > const stressxy,
          arrayView2d< real32 const > const stressxz,
          arrayView2d< real32 const > const stressyz,
          arrayView1d< const real32 > const mass,
          arrayView1d< real32 const > const dampingx,
          arrayView1d< real32 const > const dampingy,
          arrayView1d< real32 const > const dampingz,
          real64 const dt,
          arrayView1d< real32 > const ux_np1,
          arrayView1d< real32 > const uy_np1,
          arrayView1d< real32 > const uz_np1 )
  {

    forAll< EXEC_POLICY >( size_node, [=] GEOS_HOST_DEVICE ( localIndex const a )
    {
      ux_np1[a] *= 1.0-((dt/2)*(dampingx[a]/mass[a]));
      uy_np1[a] *= 1.0-((dt/2)*(dampingy[a]/mass[a]));
      uz_np1[a] *= 1.0-((dt/2)*(dampingz[a]/mass[a]));
    } );

    forAll< EXEC_POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const k )
    {
      // only the eight corners of the mesh cell are needed to compute the Jacobian
      real64 xLocal[8][3];
      for( localIndex a=0; a< 8; ++a )
      {
        localIndex const nodeIndex = elemsToNodes( k, FE_TYPE::meshIndexToLinearIndex3D( a ) );
        for( localIndex i=0; i<3; ++i )
        {
          xLocal[a][i] = nodeCoords( nodeIndex, i );
        }
      }

      constexpr localIndex numNodesPerElem = FE_TYPE::numNodes;
      real32 uelemx[numNodesPerElem] = {0.0};
      real32 uelemy[numNodesPerElem] = {0.0};
      real32 uelemz[numNodesPerElem] = {0.0};
      real32 flowx[numNodesPerElem] = {0.0};
      real32 flowy[numNodesPerElem] = {0.0};
      real32 flowz[numNodesPerElem] = {0.0};

      for( localIndex q = 0; q < numNodesPerElem; q++ )
      {


        // Stiffness part
        m_finiteElement.template computeFirstOrderStiffnessTermX( q, xLocal, [&] ( int i, int j, real32 dfx1, real32 dfx2, real32 dfx3 )
        {
          flowx[i] -= stressxx[k][j]*dfx1 + stressxy[k][j]*dfx2 + stressxz[k][j]*dfx3;
          flowy[i] -= stressxy[k][j]*dfx1 + stressyy[k][j]*dfx2 + stressyz[k][j]*dfx3;
          flowz[i] -= stressxz[k][j]*dfx1 + stressyz[k][j]*dfx2 + stresszz[k][j]*dfx3;

        } );

        m_finiteElement.template computeFirstOrderStiffnessTermY( q, xLocal, [&] ( int i, int j, real32 dfy1, real32 dfy2, real32 dfy3 )
        {
          flowx[i] -= stressxx[k][j]*dfy1 + stressxy[k][j]*dfy2 + stressxz[k][j]*dfy3;
          flowy[i] -= stressxy[k][j]*dfy1 + stressyy[k][j]*dfy2 + stressyz[k][j]*dfy3;
          flowz[i] -= stressxz[k][j]*dfy1 + stressyz[k][j]*dfy2 + stresszz[k][j]*dfy3;
        } );

        m_finiteElement.template computeFirstOrderStiffnessTermZ( q, xLocal, [&] ( int i, int j, real32 dfz1, real32 dfz2, real32 dfz3 )
        {
          flowx[i] -= stressxx[k][j]*dfz1 + stressxy[k][j]*dfz2 + stressxz[k][j]*dfz3;
          flowy[i] -= stressxy[k][j]*dfz1 + stressyy[k][j]*dfz2 + stressyz[k][j]*dfz3;
          flowz[i] -= stressxz[k][j]*dfz1 + stressyz[k][j]*dfz2 + stresszz[k][j]*dfz3;
        } );

      }
      // Time update
      for( localIndex i = 0; i < numNodesPerElem; ++i )
      {
        uelemx[i]+=dt*flowx[i];
        uelemy[i]+=dt*flowy[i];
        uelemz[i]+=dt*flowz[i];
      }

      // Mult by inverse mass matrix + damping matrix
      for( localIndex i = 0; i < numNodesPerElem; ++i )
      {
        real32 localIncrement1 = uelemx[i]/mass[elemsToNodes[k][i]];
        real32 localIncrement2 = uelemy[i]/mass[elemsToNodes[k][i]];
        real32 localIncrement3 = uelemz[i]/mass[elemsToNodes[k][i]];
        RAJA::atomicAdd< ATOMIC_POLICY >( &ux_np1[elemsToNodes[k][i]], localIncrement1 );
        RAJA::atomicAdd< ATOMIC_POLICY >( &uy_np1[elemsToNodes[k][i]], localIncrement2 );
        RAJA::atomicAdd< ATOMIC_POLICY >( &uz_np1[elemsToNodes[k][i]], localIncrement3 );
      }

    } );
    forAll< EXEC_POLICY >( size_node, [=] GEOS_HOST_DEVICE ( localIndex const a )
    {
      ux_np1[a] /= 1.0+((dt/2)*(dampingx[a]/mass[a]));
      uy_np1[a] /= 1.0+((dt/2)*(dampingy[a]/mass[a]));
      uz_np1[a] /= 1.0+((dt/2)*(dampingz[a]/mass[a]));
    } );
  }

  /// The finite element space/discretization object for the element type in the subRegion
  FE_TYPE const & m_finiteElement;

};


} // namespace ElasticFirstOrderWaveEquationSEMKernels

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ELASTICFIRSTORDERWAVEEQUATIONSEMKERNEL_HPP_
