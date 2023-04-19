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
 * @file AcousticFirstOrderWaveEquationSEMKernel.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICFIRSTTORDERWAVEEQUATIONSEMKERNEL_HPP_
#define GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICFIRSTTORDERWAVEEQUATIONSEMKERNEL_HPP_

#include "finiteElement/kernelInterface/KernelBase.hpp"
#include "WaveSolverUtils.hpp"



namespace geos
{

/// Namespace to contain the first order acoustic wave kernels.
namespace acousticFirstOrderWaveEquationSEMKernels
{

struct PrecomputeSourceAndReceiverKernel
{

  /**
   * @brief Launches the precomputation of the source and receiver terms
   * @tparam EXEC_POLICY execution policy
   * @tparam FE_TYPE finite element type
   * @param[in] size the number of cells in the subRegion
   * @param[in] numNodesPerElem number of nodes per element
   * @param[in] numFacesPerElem number of faces per element
   * @param[in] X coordinates of the nodes
   * @param[in] elemGhostRank rank of the ghost element
   * @param[in] elemsToNodes map from element to nodes
   * @param[in] elemsToFaces map from element to faces
   * @param[in] elemCenter coordinates of the element centers
   * @param[in] faceNormal normal of each faces
   * @param[in] faceCenter coordinates of the center of a face
   * @param[in] sourceCoordinates coordinates of the source terms
   * @param[out] sourceIsAccessible flag indicating whether the source is accessible or not
   * @param[out] sourceElem element where a source is located
   * @param[out] sourceNodeIds indices of the nodes of the element where the source is located
   * @param[out] sourceConstants constant part of the source terms
   * @param[in] receiverCoordinates coordinates of the receiver terms
   * @param[out] receiverIsLocal flag indicating whether the receiver is local or not
   * @param[out] rcvElem element where a receiver is located
   * @param[out] receiverNodeIds indices of the nodes of the element where the receiver is located
   * @param[out] receiverConstants constant part of the receiver term
   * @param[out] sourceValue value of the temporal source (eg. Ricker)
   * @param[in] dt time-step
   * @param[in] timeSourceFrequency the central frequency of the source
   * @param[in] rickerOrder order of the Ricker wavelet
   */
  template< typename EXEC_POLICY, typename FE_TYPE >
  static void
  launch( localIndex const size,
          localIndex const numNodesPerElem,
          localIndex const numFacesPerElem,
          arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X,
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
          arrayView2d< real64 const > const receiverCoordinates,
          arrayView1d< localIndex > const receiverIsLocal,
          arrayView1d< localIndex > const rcvElem,
          arrayView2d< localIndex > const receiverNodeIds,
          arrayView2d< real64 > const receiverConstants,
          arrayView2d< real32 > const sourceValue,
          real64 const dt,
          real32 const timeSourceFrequency,
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
                                                                              X,
                                                                              coordsOnRefElem );

            sourceIsAccessible[isrc] = 1;
            sourceElem[isrc] = k;
            real64 Ntest[FE_TYPE::numNodes];
            FE_TYPE::calcN( coordsOnRefElem, Ntest );

            for( localIndex a = 0; a < numNodesPerElem; ++a )
            {
              sourceNodeIds[isrc][a] = elemsToNodes[k][a];
              sourceConstants[isrc][a] = Ntest[a];
            }

            for( localIndex cycle = 0; cycle < sourceValue.size( 0 ); ++cycle )
            {
              real64 const time = cycle*dt;
              sourceValue[cycle][isrc] = WaveSolverUtils::evaluateRicker( time, timeSourceFrequency, rickerOrder );
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
                                                                              X,
                                                                              coordsOnRefElem );
            receiverIsLocal[ircv] = 1;
            rcvElem[ircv] = k;

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
   * @brief Launches the precomputation of the mass matrix
   * @tparam EXEC_POLICY the execution policy
   * @tparam ATOMIC_POLICY the atomic policy
   * @param[in] size the number of cells in the subRegion
   * @param[in] X coordinates of the nodes
   * @param[in] elemsToNodes map from element to nodes
   * @param[in] velocity cell-wise velocity
   * @param[in] density cell-wise density
   * @param[out] mass diagonal of the mass matrix
   */
  template< typename EXEC_POLICY, typename ATOMIC_POLICY >
  void
  launch( localIndex const size,
          arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X,
          arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodes,
          arrayView1d< real32 const > const velocity,
          arrayView1d< real32 const > const density,
          arrayView1d< real32 > const mass )

  {
    forAll< EXEC_POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const k )
    {

      constexpr localIndex numNodesPerElem = FE_TYPE::numNodes;
      constexpr localIndex numQuadraturePointsPerElem = FE_TYPE::numQuadraturePoints;

      real32 const invC2 = 1.0 / ( density[k] * velocity[k] * velocity[k] );
      real64 xLocal[ numNodesPerElem ][ 3 ];
      for( localIndex a = 0; a < numNodesPerElem; ++a )
      {
        for( localIndex i = 0; i < 3; ++i )
        {
          xLocal[a][i] = X( elemsToNodes( k, a ), i );
        }
      }

      for( localIndex q = 0; q < numQuadraturePointsPerElem; ++q )
      {
        real32 const localIncrement = invC2 * m_finiteElement.computeMassTerm( q, xLocal );
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
   * @param[in] X coordinates of the nodes
   * @param[in] facesToElems map from faces to elements
   * @param[in] facesToNodes map from face to nodes
   * @param[in] facesDomainBoundaryIndicator flag equal to 1 if the face is on the boundary, and to 0 otherwise
   * @param[in] freeSurfaceFaceIndicator flag equal to 1 if the face is on the free surface, and to 0 otherwise
   * @param[in] velocity cell-wise velocity
   * @param[out] damping diagonal of the damping matrix
   */
  template< typename EXEC_POLICY, typename ATOMIC_POLICY >
  void
  launch( localIndex const size,
          arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X,
          arrayView2d< localIndex const > const facesToElems,
          ArrayOfArraysView< localIndex const > const facesToNodes,
          arrayView1d< integer const > const facesDomainBoundaryIndicator,
          arrayView1d< localIndex const > const freeSurfaceFaceIndicator,
          arrayView1d< real32 const > const velocity,
          arrayView1d< real32 > const damping )
  {
    forAll< EXEC_POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const f )
    {
      // face on the domain boundary and not on free surface
      if( facesDomainBoundaryIndicator[f] == 1 && freeSurfaceFaceIndicator[f] != 1 )
      {
        localIndex k = facesToElems( f, 0 );
        if( k < 0 )
        {
          k = facesToElems( f, 1 );
        }

        real32 const alpha = 1.0 / velocity[k];

        constexpr localIndex numNodesPerFace = FE_TYPE::numNodesPerFace;

        real64 xLocal[ numNodesPerFace ][ 3 ];
        for( localIndex a = 0; a < numNodesPerFace; ++a )
        {
          for( localIndex i = 0; i < 3; ++i )
          {
            xLocal[a][i] = X( facesToNodes( k, a ), i );
          }
        }

        for( localIndex q = 0; q < numNodesPerFace; ++q )
        {
          real32 const localIncrement = alpha * m_finiteElement.computeDampingTerm( q, xLocal );
          RAJA::atomicAdd< ATOMIC_POLICY >( &damping[facesToNodes[f][q]], localIncrement );
        }
      }
    } ); // end loop over element
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
   * @brief Launches the computation of the velocity for one iteration
   * @tparam EXEC_POLICY the execution policy
   * @tparam ATOMIC_POLICY the atomic policy
   * @param[in] size the number of cells in the subRegion
   * @param[in] X coordinates of the nodes
   * @param[in] elemsToNodes map from element to nodes
   * @param[in] p_np1 pressure array (only used here)
   * @param[in] dt time-step
   * @param[out] velocity_x velocity array in the x direction (updated here)
   * @param[out] velocity_y velocity array in the y direction (updated here)
   * @param[out] velocity_z velocity array in the z direction (updated here)
   */
  template< typename EXEC_POLICY, typename ATOMIC_POLICY >
  void
  launch( localIndex const size,
          arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X,
          arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodes,
          arrayView1d< real32 const > const p_np1,
          arrayView1d< real32 const > const density,
          real64 const dt,
          arrayView2d< real32 > const velocity_x,
          arrayView2d< real32 > const velocity_y,
          arrayView2d< real32 > const velocity_z )
  {
    forAll< EXEC_POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const k )
    {
      constexpr localIndex numNodesPerElem = FE_TYPE::numNodes;
      constexpr localIndex numQuadraturePointsPerElem = FE_TYPE::numQuadraturePoints;

      real64 xLocal[numNodesPerElem][3];
      for( localIndex a=0; a< numNodesPerElem; ++a )
      {
        for( localIndex i=0; i<3; ++i )
        {
          xLocal[a][i] = X( elemsToNodes( k, a ), i );
        }
      }

      real32 uelemx[numNodesPerElem] = {0.0};
      real32 uelemy[numNodesPerElem] = {0.0};
      real32 uelemz[numNodesPerElem] = {0.0};
      real32 flowx[numNodesPerElem] = {0.0};
      real32 flowy[numNodesPerElem] = {0.0};
      real32 flowz[numNodesPerElem] = {0.0};

      for( localIndex i = 0; i < numNodesPerElem; ++i )
      {
        real32 massLoc = m_finiteElement.computeMassTerm( i, xLocal );
        uelemx[i] = massLoc*velocity_x[k][i];
        uelemy[i] = massLoc*velocity_y[k][i];
        uelemz[i] = massLoc*velocity_z[k][i];
      }

      for( localIndex q=0; q<numQuadraturePointsPerElem; ++q )
      {



        m_finiteElement.template computeFirstOrderStiffnessTermX( q, xLocal, [&] ( int i, int j, real32 dfx1, real32 dfx2, real32 dfx3 )
        {
          flowx[j] += dfx1*p_np1[elemsToNodes[k][i]];
          flowy[j] += dfx2*p_np1[elemsToNodes[k][i]];
          flowz[j] += dfx3*p_np1[elemsToNodes[k][i]];
        } );

        m_finiteElement.template computeFirstOrderStiffnessTermY( q, xLocal, [&] ( int i, int j, real32 dfy1, real32 dfy2, real32 dfy3 )
        {
          flowx[j] += dfy1*p_np1[elemsToNodes[k][i]];
          flowy[j] += dfy2*p_np1[elemsToNodes[k][i]];
          flowz[j] += dfy3*p_np1[elemsToNodes[k][i]];
        } );

        m_finiteElement.template computeFirstOrderStiffnessTermZ( q, xLocal, [&] ( int i, int j, real32 dfz1, real32 dfz2, real32 dfz3 )
        {
          flowx[j] += dfz1*p_np1[elemsToNodes[k][i]];
          flowy[j] += dfz2*p_np1[elemsToNodes[k][i]];
          flowz[j] += dfz3*p_np1[elemsToNodes[k][i]];

        } );

      }
      for( localIndex i = 0; i < numNodesPerElem; ++i )
      {
        real32 massLoc = m_finiteElement.computeMassTerm( i, xLocal );
        uelemx[i]+=dt*flowx[i]/density[k];
        uelemy[i]+=dt*flowy[i]/density[k];
        uelemz[i]+=dt*flowz[i]/density[k];

        velocity_x[k][i] = uelemx[i]/massLoc;
        velocity_y[k][i] = uelemy[i]/massLoc;
        velocity_z[k][i] = uelemz[i]/massLoc;
      }
    } );
  }

  /// The finite element space/discretization object for the element type in the subRegion
  FE_TYPE const & m_finiteElement;
};

template< typename FE_TYPE >
struct PressureComputation
{

  PressureComputation( FE_TYPE const & finiteElement )
    : m_finiteElement( finiteElement )
  {}

  /**
   * @brief Launches the computation of the pressure for one iteration
   * @tparam EXEC_POLICY the execution policy
   * @tparam ATOMIC_POLICY the atomic policy
   * @param[in] size the number of cells in the subRegion
   * @param[in] size_node the number of nodes in the subRegion
   * @param[in] X coordinates of the nodes
   * @param[in] elemsToNodes map from element to nodes
   * @param[out] velocity_x velocity array in the x direction (only used here)
   * @param[out] velocity_y velocity array in the y direction (only used here)
   * @param[out] velocity_z velocity array in the z direction (only used here)
   * @param[in] mass the mass matrix
   * @param[in] damping the damping matrix
   * @param[in] sourceConstants constant part of the source terms
   * @param[in] sourceValue value of the temporal source (eg. Ricker)
   * @param[in] sourceIsAccessible flag indicating whether the source is accessible or not
   * @param[in] sourceElem element where a source is located
   * @param[in] cycleNumber the number of cycle
   * @param[in] dt time-step
   * @param[out] p_np1 pressure array (updated here)
   */

  template< typename EXEC_POLICY, typename ATOMIC_POLICY >
  void
  launch( localIndex const size,
          localIndex const size_node,
          arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X,
          arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodes,
          arrayView2d< real32 const > const velocity_x,
          arrayView2d< real32 const > const velocity_y,
          arrayView2d< real32 const > const velocity_z,
          arrayView1d< real32 const > const mass,
          arrayView1d< real32 const > const damping,
          arrayView2d< real64 const > const sourceConstants,
          arrayView2d< real32 const > const sourceValue,
          arrayView1d< localIndex const > const sourceIsAccessible,
          arrayView1d< localIndex const > const sourceElem,
          real64 const dt,
          integer const cycleNumber,
          arrayView1d< real32 > const p_np1 )

  {

    //Pre-mult by the first factor for damping
    forAll< EXEC_POLICY >( size_node, [=] GEOS_HOST_DEVICE ( localIndex const a )
    {
      p_np1[a] *= 1.0-((dt/2)*(damping[a]/mass[a]));
    } );

    forAll< EXEC_POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const k )
    {
      constexpr localIndex numNodesPerElem = FE_TYPE::numNodes;
      constexpr localIndex numQuadraturePointsPerElem = FE_TYPE::numQuadraturePoints;

      real64 xLocal[numNodesPerElem][3];
      for( localIndex a=0; a< numNodesPerElem; ++a )
      {
        for( localIndex i=0; i<3; ++i )
        {
          xLocal[a][i] = X( elemsToNodes( k, a ), i );
        }
      }

      real32 auxx[numNodesPerElem]  = {0.0};
      real32 auyy[numNodesPerElem]  = {0.0};
      real32 auzz[numNodesPerElem]  = {0.0};
      real32 uelemx[numNodesPerElem] = {0.0};


      // Volume integration
      for( localIndex q=0; q<numQuadraturePointsPerElem; ++q )
      {


        m_finiteElement.template computeFirstOrderStiffnessTermX( q, xLocal, [&] ( int i, int j, real32 dfx1, real32 dfx2, real32 dfx3 )
        {
          auxx[i] -= dfx1*velocity_x[k][j];
          auyy[i] -= dfx2*velocity_y[k][j];
          auzz[i] -= dfx3*velocity_z[k][j];
        } );

        m_finiteElement.template computeFirstOrderStiffnessTermY( q, xLocal, [&] ( int i, int j, real32 dfy1, real32 dfy2, real32 dfy3 )
        {
          auxx[i] -= dfy1*velocity_x[k][j];
          auyy[i] -= dfy2*velocity_y[k][j];
          auzz[i] -= dfy3*velocity_z[k][j];
        } );

        m_finiteElement.template computeFirstOrderStiffnessTermZ( q, xLocal, [&] ( int i, int j, real32 dfz1, real32 dfz2, real32 dfz3 )
        {
          auxx[i] -= dfz1*velocity_x[k][j];
          auyy[i] -= dfz2*velocity_y[k][j];
          auzz[i] -= dfz3*velocity_z[k][j];
        } );



      }

      //Time update + multiplication by inverse of the mass matrix
      for( localIndex i = 0; i < numNodesPerElem; ++i )
      {
        real32 diag=(auxx[i]+auyy[i]+auzz[i]);
        uelemx[i]+=dt*diag;

        real32 const localIncrement = uelemx[i]/mass[elemsToNodes[k][i]];
        RAJA::atomicAdd< ATOMIC_POLICY >( &p_np1[elemsToNodes[k][i]], localIncrement );
      }

      //Source Injection
      for( localIndex isrc = 0; isrc < sourceConstants.size( 0 ); ++isrc )
      {
        if( sourceIsAccessible[isrc] == 1 )
        {
          if( sourceElem[isrc]==k )
          {
            for( localIndex i = 0; i < numNodesPerElem; ++i )
            {
              real32 const localIncrement2 = dt*(sourceConstants[isrc][i]*sourceValue[cycleNumber][isrc])/(mass[elemsToNodes[k][i]]);
              RAJA::atomicAdd< ATOMIC_POLICY >( &p_np1[elemsToNodes[k][i]], localIncrement2 );
            }
          }
        }
      }

    } );

    //Pre-mult by the first factor for damping
    forAll< EXEC_POLICY >( size_node, [=] GEOS_HOST_DEVICE ( localIndex const a )
    {
      p_np1[a] /= 1.0+((dt/2)*(damping[a]/mass[a]));
    } );
  }

  /// The finite element space/discretization object for the element type in the subRegion
  FE_TYPE const & m_finiteElement;

};

} // namespace AcousticFirstOrderWaveEquationSEMKernels

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICFIRSTTORDERWAVEEQUATIONSEMKERNEL_HPP_
