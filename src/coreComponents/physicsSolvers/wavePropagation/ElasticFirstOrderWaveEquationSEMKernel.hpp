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

struct PrecomputeSourceAndReceiverKernel
{

  /**
   * @brief Launches the precomputation of the source and receiver terms
   * @tparam EXEC_POLICY execution policy
   * @tparam FE_TYPE finite element type
   * @param[in] size the number of cells in the subRegion
   * @param[in] numFacesPerElem number of face on an element
   * @param[in] nodeCoords coordinates of the nodes
   * @param[in] elemGhostRank array containing the ghost rank
   * @param[in] elemsToNodes map from element to nodes
   * @param[in] elemsToFaces map from element to faces
   * @param[in] elemCenter coordinates of the element centers
   * @param[in] faceNormal array containing the normal of all faces
   * @param[in] faceCenter array containing the center of all faces
   * @param[in] sourceCoordinates coordinates of the source terms
   * @param[in] receiverCoordinates coordinates of the receiver terms
   * @param[in] dt time-step
   * @param[in] timeSourceFrequency Peak frequency of the source
   * @param[in] timeSourceDelay the time delay of the source
   * @param[in] rickerOrder Order of the Ricker wavelet
   * @param[out] sourceNodeIds indices of the nodes of the element where the source is located
   * @param[out] sourceConstants constant part of the source terms
   * @param[out] receiverIsLocal flag indicating whether the receiver is local or not
   * @param[out] receiverNodeIds indices of the nodes of the element where the receiver is located
   * @param[out] sourceValue array containing the value of the time dependent source (Ricker for e.g)
   */
  template< typename EXEC_POLICY, typename FE_TYPE >
  static void
  launch( localIndex const size,
          localIndex const regionIndex,
          localIndex const numNodesPerElem,
          localIndex const numFacesPerElem,
          arrayView2d< WaveSolverBase::wsCoordType const, nodes::REFERENCE_POSITION_USD > const nodeCoords,
          arrayView1d< integer const > const elemGhostRank,
          arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes,
          arrayView2d< localIndex const > const elemsToFaces,
          arrayView2d< real64 const > const & elemCenter,
          arrayView2d< real64 const > const faceNormal,
          arrayView2d< real64 const > const faceCenter,
          arrayView2d< real64 const > const sourceCoordinates,
          arrayView1d< localIndex > const sourceIsAccessible,
          arrayView1d< localIndex > const sourceElem,
          arrayView2d< localIndex > const sourceNodeIds,
          arrayView2d< real64 > const sourceConstants,
          arrayView1d< localIndex > const sourceRegion,
          arrayView2d< real64 const > const receiverCoordinates,
          arrayView1d< localIndex > const receiverIsLocal,
          arrayView1d< localIndex > const rcvElem,
          arrayView2d< localIndex > const receiverNodeIds,
          arrayView2d< real64 > const receiverConstants,
          arrayView1d< localIndex > const receiverRegion,
          arrayView2d< real32 > const sourceValue,
          real64 const dt,
          real32 const timeSourceFrequency,
          real32 const timeSourceDelay,
          localIndex const rickerOrder )
  {

    forAll< EXEC_POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const k )
    {
      real64 const center[3] = { elemCenter[k][0],
                                 elemCenter[k][1],
                                 elemCenter[k][2] };

      // Step 1: locate the sources, and precompute the source term

      /// loop over all the source that haven't been found yet
      for( localIndex isrc = 0; isrc < sourceCoordinates.size( 0 ); ++isrc )
      {
        if( sourceIsAccessible[isrc] == 0 )
        {
          real64 const coords[3] = { sourceCoordinates[isrc][0],
                                     sourceCoordinates[isrc][1],
                                     sourceCoordinates[isrc][2] };

          bool const sourceFound =
            WaveSolverUtils::locateSourceElement( numFacesPerElem,
                                                  center,
                                                  faceNormal,
                                                  faceCenter,
                                                  elemsToFaces[k],
                                                  coords );

          if( sourceFound )
          {
            real64 coordsOnRefElem[3]{};

            WaveSolverUtils::computeCoordinatesOnReferenceElement< FE_TYPE >( coords,
                                                                              elemsToNodes[k],
                                                                              nodeCoords,
                                                                              coordsOnRefElem );
            sourceIsAccessible[isrc] = 1;
            sourceElem[isrc] = k;
            sourceRegion[isrc] = regionIndex;
            real64 Ntest[FE_TYPE::numNodes];
            FE_TYPE::calcN( coordsOnRefElem, Ntest );

            for( localIndex a = 0; a < numNodesPerElem; ++a )
            {
              sourceNodeIds[isrc][a] = elemsToNodes[k][a];
              sourceConstants[isrc][a] = Ntest[a];
            }

            for( localIndex cycle = 0; cycle < sourceValue.size( 0 ); ++cycle )
            {
              sourceValue[cycle][isrc] = WaveSolverUtils::evaluateRicker( cycle * dt, timeSourceFrequency, timeSourceDelay, rickerOrder );
            }

          }
        }
      } // end loop over all sources


      // Step 2: locate the receivers, and precompute the receiver term

      /// loop over all the receivers that haven't been found yet
      for( localIndex ircv = 0; ircv < receiverCoordinates.size( 0 ); ++ircv )
      {
        if( receiverIsLocal[ircv] == 0 )
        {
          real64 const coords[3] = { receiverCoordinates[ircv][0],
                                     receiverCoordinates[ircv][1],
                                     receiverCoordinates[ircv][2] };

          real64 coordsOnRefElem[3]{};
          bool const receiverFound =
            WaveSolverUtils::locateSourceElement( numFacesPerElem,
                                                  center,
                                                  faceNormal,
                                                  faceCenter,
                                                  elemsToFaces[k],
                                                  coords );
          if( receiverFound && elemGhostRank[k] < 0 )
          {
            WaveSolverUtils::computeCoordinatesOnReferenceElement< FE_TYPE >( coords,
                                                                              elemsToNodes[k],
                                                                              nodeCoords,
                                                                              coordsOnRefElem );
            receiverIsLocal[ircv] = 1;
            rcvElem[ircv] = k;
            receiverRegion[ircv] = regionIndex;

            real64 Ntest[FE_TYPE::numNodes];
            FE_TYPE::calcN( coordsOnRefElem, Ntest );

            for( localIndex a = 0; a < numNodesPerElem; ++a )
            {
              receiverNodeIds[ircv][a] = elemsToNodes[k][a];
              receiverConstants[ircv][a] = Ntest[a];
            }
          }
        }
      } // end loop over receivers
    } );

  }
};

template< typename FE_TYPE >
struct MassMatrixKernel
{

  MassMatrixKernel( FE_TYPE const & finiteElement )
    : m_finiteElement( finiteElement )
  {}

  /**
   * @brief Launches the precomputation of the mass matrices
   * @tparam EXEC_POLICY the execution policy
   * @tparam ATOMIC_POLICY the atomic policy
   * @param[in] size the number of cells in the subRegion
   * @param[in] numFacesPerElem number of faces per element
   * @param[in] nodeCoords coordinates of the nodes
   * @param[in] elemsToNodes map from element to nodes
   * @param[in] velocity cell-wise velocity
   * @param[out] mass diagonal of the mass matrix
   */
  template< typename EXEC_POLICY, typename ATOMIC_POLICY >
  void
  launch( localIndex const size,
          arrayView2d< WaveSolverBase::wsCoordType const, nodes::REFERENCE_POSITION_USD > const nodeCoords,
          arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodes,
          arrayView1d< real32 const > const density,
          arrayView1d< real32 > const mass )

  {
    forAll< EXEC_POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const k )
    {

      constexpr localIndex numNodesPerElem = FE_TYPE::numNodes;
      constexpr localIndex numQuadraturePointsPerElem = FE_TYPE::numQuadraturePoints;

      real64 xLocal[ 8 ][ 3 ];
      for( localIndex a = 0; a < 8; ++a )
      {
        localIndex const nodeIndex = elemsToNodes( k, FE_TYPE::meshIndexToLinearIndex3D( a ) );
        for( localIndex i = 0; i < 3; ++i )
        {
          xLocal[a][i] = nodeCoords( nodeIndex, i );
        }
      }

      for( localIndex q = 0; q < numQuadraturePointsPerElem; ++q )
      {
        real32 const localIncrement = density[k] * m_finiteElement.computeMassTerm( q, xLocal );
        RAJA::atomicAdd< ATOMIC_POLICY >( &mass[elemsToNodes[k][q]], localIncrement );
      }
    } ); // end loop over element
  }

  /// The finite element space/discretization object for the element type in the subRegion
  FE_TYPE const & m_finiteElement;

};

template< typename FE_TYPE >
struct DampingMatrixKernel
{

  DampingMatrixKernel( FE_TYPE const & finiteElement )
    : m_finiteElement( finiteElement )
  {}

  /**
   * @brief Launches the precomputation of the damping matrices
   * @tparam EXEC_POLICY the execution policy
   * @tparam ATOMIC_POLICY the atomic policy
   * @param[in] size the number of cells in the subRegion
   * @param[in] nodeCoords coordinates of the nodes
   * @param[in] elemsToFaces map from elements to faces
   * @param[in] facesToNodes map from face to nodes
   * @param[in] facesDomainBoundaryIndicator flag equal to 1 if the face is on the boundary, and to 0 otherwise
   * @param[in] freeSurfaceFaceIndicator flag equal to 1 if the face is on the free surface, and to 0 otherwise
   * @param[in] velocity cell-wise velocity
   * @param[out] damping diagonal of the damping matrix
   */
  template< typename EXEC_POLICY, typename ATOMIC_POLICY >
  void
  launch( localIndex const size,
          arrayView2d< WaveSolverBase::wsCoordType const, nodes::REFERENCE_POSITION_USD > const nodeCoords,
          arrayView2d< localIndex const > const elemsToFaces,
          ArrayOfArraysView< localIndex const > const facesToNodes,
          arrayView1d< integer const > const facesDomainBoundaryIndicator,
          arrayView1d< localIndex const > const freeSurfaceFaceIndicator,
          arrayView2d< real64 const > const faceNormal,
          arrayView1d< real32 const > const density,
          arrayView1d< real32 const > const velocityVp,
          arrayView1d< real32 const > const velocityVs,
          arrayView1d< real32 > const dampingx,
          arrayView1d< real32 > const dampingy,
          arrayView1d< real32 > const dampingz )
  {
    forAll< EXEC_POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const e )
    {
      for( localIndex i = 0; i < elemsToFaces.size( 1 ); ++i )
      {
        localIndex const f = elemsToFaces( e, i );
        // face on the domain boundary and not on free surface
        if( facesDomainBoundaryIndicator[f] == 1 && freeSurfaceFaceIndicator[f] != 1 )
        {
          constexpr localIndex numNodesPerFace = FE_TYPE::numNodesPerFace;
          real64 xLocal[ 4 ][ 3 ];
          for( localIndex a = 0; a < 4; ++a )
          {
            localIndex const nodeIndex = facesToNodes( f, FE_TYPE::meshIndexToLinearIndex2D( a ) );
            for( localIndex d = 0; d < 3; ++d )
            {
              xLocal[a][d] = nodeCoords( nodeIndex, d );
            }
          }

          real32 const nx = faceNormal( f, 0 ), ny = faceNormal( f, 1 ), nz = faceNormal( f, 2 );
          for( localIndex q = 0; q < numNodesPerFace; ++q )
          {
            real32 const aux = density[e] * m_finiteElement.computeDampingTerm( q, xLocal );
            real32 const localIncrementx = density[e] * (velocityVp[e] * abs( nx ) + velocityVs[e] * sqrt( pow( ny, 2 ) + pow( nz, 2 ) ) ) * aux;
            real32 const localIncrementy = density[e] * (velocityVp[e] * abs( ny ) + velocityVs[e] * sqrt( pow( nx, 2 ) + pow( nz, 2 ) ) ) * aux;
            real32 const localIncrementz = density[e] * (velocityVp[e] * abs( nz ) + velocityVs[e] * sqrt( pow( nx, 2 ) + pow( ny, 2 ) ) ) * aux;

            RAJA::atomicAdd< ATOMIC_POLICY >( &dampingx[facesToNodes( f, q )], localIncrementx );
            RAJA::atomicAdd< ATOMIC_POLICY >( &dampingy[facesToNodes( f, q )], localIncrementy );
            RAJA::atomicAdd< ATOMIC_POLICY >( &dampingz[facesToNodes( f, q )], localIncrementz );
          }
        }
      }
    } );
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
      constexpr localIndex numNodesPerElem = FE_TYPE::numNodes;
      constexpr localIndex numQuadraturePointsPerElem = FE_TYPE::numQuadraturePoints;
      real64 xLocal[8][3];
      for( localIndex a=0; a< 8; ++a )
      {
        localIndex const nodeIndex = elemsToNodes( k, FE_TYPE::meshIndexToLinearIndex3D( a ) );
        for( localIndex i=0; i<3; ++i )
        {
          xLocal[a][i] = nodeCoords( nodeIndex, i );
        }
      }

      mu[k] = density[k] * velocityVs[k] * velocityVs[k];
      lambda[k] = density[k] * velocityVp[k] * velocityVp[k] - 2.0*mu[k];

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

      constexpr localIndex numNodesPerElem = FE_TYPE::numNodes;
      constexpr localIndex numQuadraturePointsPerElem = FE_TYPE::numQuadraturePoints;

      real64 xLocal[8][3];
      for( localIndex a=0; a< 8; ++a )
      {
        localIndex const nodeIndex = elemsToNodes( k, FE_TYPE::meshIndexToLinearIndex3D( a ) );
        for( localIndex i=0; i<3; ++i )
        {
          xLocal[a][i] = nodeCoords( nodeIndex, i );
        }
      }

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

#endif //GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ElasticFirstOrderWaveEquationSEMKERNEL_HPP_
