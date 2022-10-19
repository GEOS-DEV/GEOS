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

#ifndef GEOSX_PHYSICSSOLVERS_WAVEPROPAGATION_AcousticFirstOrderWaveEquationSEMKERNEL_HPP_
#define GEOSX_PHYSICSSOLVERS_WAVEPROPAGATION_AcousticFirstOrderWaveEquationSEMKERNEL_HPP_

#include "finiteElement/kernelInterface/KernelBase.hpp"


namespace geosx
{

/// Namespace to contain the first order acoustic wave kernels.
namespace acousticFirstOrderWaveEquationSEMKernels
{

struct PrecomputeSourceAndReceiverKernel
{

  /**
   * @brief Check if the sourtc epoint is inside an element or not
   */

  GEOSX_HOST_DEVICE
  static bool
  locateSourceElement( real64 const numFacesPerElem,
                       real64 const (&elemCenter)[3],
                       arrayView2d< real64 const > const faceNormal,
                       arrayView2d< real64 const > const faceCenter,
                       arraySlice1d< localIndex const > const elemsToFaces,
                       real64 const (&coords)[3] )
  {
    //Loop over the element faces
    real64 tmpVector[3]{};
    for( localIndex kfe = 0; kfe < numFacesPerElem; ++kfe )
    {

      localIndex const iface = elemsToFaces[kfe];
      real64 faceCenterOnFace[3] = {faceCenter[iface][0],
                                    faceCenter[iface][1],
                                    faceCenter[iface][2]};
      real64 faceNormalOnFace[3] = {faceNormal[iface][0],
                                    faceNormal[iface][1],
                                    faceNormal[iface][2]};

      //Test to make sure if the normal is outwardly directed
      LvArray::tensorOps::copy< 3 >( tmpVector, faceCenterOnFace );
      LvArray::tensorOps::subtract< 3 >( tmpVector, elemCenter );
      if( LvArray::tensorOps::AiBi< 3 >( tmpVector, faceNormalOnFace ) < 0.0 )
      {
        LvArray::tensorOps::scale< 3 >( faceNormalOnFace, -1 );
      }

      // compute the vector face center to query point
      LvArray::tensorOps::subtract< 3 >( faceCenterOnFace, coords );
      localIndex const s = computationalGeometry::sign( LvArray::tensorOps::AiBi< 3 >( faceNormalOnFace, faceCenterOnFace ));

      // all dot products should be non-negative (we enforce outward normals)
      if( s < 0 )
      {
        return false;
      }

    }
    return true;
  }

  /**
   * @brief Convert a mesh element point coordinate into a coordinate on the reference element
   * @tparam FE_TYPE finite element type
   * @param[in] coords coordinate of the point
   * @param[in] elemCenter coordinate of the cell center
   * @param[in] elemsToNodes map to obtain global nodes from element index
   * @param[in] elemsToFaces map to obtain global faces from element index
   * @param[in] facesToNodes map to obtain global nodes from global faces
   * @param[in] X array of mesh nodes coordinates
   * @param[out] coordsOnRefElem to contain the coordinate computed in the reference element
   * @return true if coords is inside the element
   */
  template< typename FE_TYPE >
  GEOSX_HOST_DEVICE
  static void
  computeCoordinatesOnReferenceElement( real64 const (&coords)[3],
                                        arraySlice1d< localIndex const, cells::NODE_MAP_USD - 1 > const elemsToNodes,
                                        arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X,
                                        real64 (& coordsOnRefElem)[3] )
  {
    real64 xLocal[FE_TYPE::numNodes][3]{};
    for( localIndex a = 0; a < FE_TYPE::numNodes; ++a )
    {
      LvArray::tensorOps::copy< 3 >( xLocal[a], X[ elemsToNodes[a] ] );
    }
    // coordsOnRefElem = invJ*(coords-coordsNode_0)
    real64 invJ[3][3]{};
    FE_TYPE::invJacobianTransformation( 0, xLocal, invJ );
    for( localIndex i = 0; i < 3; ++i )
    {
      // init at (-1,-1,-1) as the origin of the referential elem
      coordsOnRefElem[i] = -1.0;
      for( localIndex j = 0; j < 3; ++j )
      {
        coordsOnRefElem[i] += invJ[i][j] * (coords[j] - xLocal[0][j]);
      }
    }
  }


  GEOSX_HOST_DEVICE
  static real64
  evaluateRicker( real64 const & time_n,
                  real64 const & f0,
                  localIndex const & order )
  {
    real32 const o_tpeak = 1.0/f0;
    real32 pulse = 0.0;
    if((time_n <= -0.9*o_tpeak) || (time_n >= 2.9*o_tpeak))
    {
      return pulse;
    }

    constexpr real32 pi = M_PI;
    real64 const lam = (f0*pi)*(f0*pi);

    switch( order )
    {
      case 2:
      {
        pulse = 2.0*lam*(2.0*lam*(time_n-o_tpeak)*(time_n-o_tpeak)-1.0)*exp( -lam*(time_n-o_tpeak)*(time_n-o_tpeak));
      }
      break;
      case 1:
      {
        pulse = -2.0*lam*(time_n-o_tpeak)*exp( -lam*(time_n-o_tpeak)*(time_n-o_tpeak));
      }
      break;
      case 0:
      {
        pulse = -(time_n-o_tpeak)*exp( -2*lam*(time_n-o_tpeak)*(time_n-o_tpeak) );
      }
      break;
      default:
        GEOSX_ERROR( "This option is not supported yet, rickerOrder must be 0, 1 or 2" );
    }

    return pulse;
  }


  /**
   * @brief Launches the precomputation of the source and receiver terms
   * @tparam EXEC_POLICY execution policy
   * @tparam FE_TYPE finite element type
   * @param[in] size the number of cells in the subRegion
   * @param[in] numNodesPerElem number of nodes per element
   * @param[in] X coordinates of the nodes
   * @param[in] elemsToNodes map from element to nodes
   * @param[in] elemsToFaces map from element to faces
   * @param[in] facesToNodes map from faces to nodes
   * @param[in] elemCenter coordinates of the element centers
   * @param[in] sourceCoordinates coordinates of the source terms
   * @param[out] sourceIsAccessible flag indicating whether the source is accessible or not
   * @param[out] sourceNodeIds indices of the nodes of the element where the source is located
   * @param[out] sourceNodeConstants constant part of the source terms
   * @param[in] receiverCoordinates coordinates of the receiver terms
   * @param[out] receiverIsLocal flag indicating whether the receiver is local or not
   * @param[out] receiverNodeIds indices of the nodes of the element where the receiver is located
   * @param[out] receiverNodeConstants constant part of the receiver term
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
          arrayView2d< localIndex > const receiverNodeIds,
          arrayView2d< real64 > const receiverConstants,
          arrayView2d< real32 > const sourceValue,
          real64 const dt,
          real32 const timeSourceFrequency,
          localIndex const rickerOrder )
  {

    forAll< EXEC_POLICY >( size, [=] GEOSX_HOST_DEVICE ( localIndex const k )
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
            locateSourceElement( numFacesPerElem,
                                 center,
                                 faceNormal,
                                 faceCenter,
                                 elemsToFaces[k],
                                 coords );

          if( sourceFound )
          {
            real64 coordsOnRefElem[3]{};

            computeCoordinatesOnReferenceElement< FE_TYPE >( coords,
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
              sourceValue[cycle][isrc] = evaluateRicker( time, timeSourceFrequency, rickerOrder );
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
            locateSourceElement( numFacesPerElem,
                                 center,
                                 faceNormal,
                                 faceCenter,
                                 elemsToFaces[k],
                                 coords );

          if( receiverFound && elemGhostRank[k] < 0 )
          {
            computeCoordinatesOnReferenceElement< FE_TYPE >( coords,
                                                             elemsToNodes[k],
                                                             X,
                                                             coordsOnRefElem );
            receiverIsLocal[ircv] = 1;

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
struct MassAndDampingMatrixKernel
{

  MassAndDampingMatrixKernel( FE_TYPE const & finiteElement )
    : m_finiteElement( finiteElement )
  {}

  /**
   * @brief Launches the precomputation of the mass and damping matrices
   * @tparam EXEC_POLICY the execution policy
   * @tparam ATOMIC_POLICY the atomic policy
   * @param[in] size the number of cells in the subRegion
   * @param[in] numFacesPerElem number of faces per element
   * @param[in] numNodesPerFace number of nodes per face
   * @param[in] X coordinates of the nodes
   * @param[in] elemsToNodes map from element to nodes
   * @param[in] elemsToFaces map from element to faces
   * @param[in] facesToNodes map from face to nodes
   * @param[in] facesDomainBoundaryIndicator flag equal to 1 if the face is on the boundary, and to 0 otherwise
   * @param[in] freeSurfaceFaceIndicator flag equal to 1 if the face is on the free surface, and to 0 otherwise
   * @param[in] faceNormal normal vectors at the faces
   * @param[in] velocity cell-wise velocity
   * @param[in] density cell-wise density
   * @param[out] mass diagonal of the mass matrix
   * @param[out] damping diagonal of the damping matrix
   */
  template< typename EXEC_POLICY, typename ATOMIC_POLICY >
  void
  launch( localIndex const size,
          localIndex const numFacesPerElem,
          localIndex const numNodesPerFace,
          arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X,
          arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodes,
          arrayView2d< localIndex const > const elemsToFaces,
          ArrayOfArraysView< localIndex const > const facesToNodes,
          arrayView1d< integer const > const facesDomainBoundaryIndicator,
          arrayView1d< localIndex const > const freeSurfaceFaceIndicator,
          arrayView2d< real64 const > const faceNormal,
          arrayView1d< real32 const > const velocity,
          arrayView1d< real32 const > const density,
          arrayView1d< real32 > const mass,
          arrayView1d< real32 > const damping )
  {
    forAll< EXEC_POLICY >( size, [=] GEOSX_HOST_DEVICE ( localIndex const k )
    {

      constexpr localIndex numNodesPerElem = FE_TYPE::numNodes;
      constexpr localIndex numQuadraturePointsPerElem = FE_TYPE::numQuadraturePoints;

      real32 const invC2 = 1.0 / ( velocity[k] * velocity[k] * density[k] );
      real64 xLocal[ numNodesPerElem ][ 3 ];
      for( localIndex a = 0; a < numNodesPerElem; ++a )
      {
        for( localIndex i = 0; i < 3; ++i )
        {
          xLocal[a][i] = X( elemsToNodes( k, a ), i );
        }
      }

      real64 N[ numNodesPerElem ];
      real64 gradN[ numNodesPerElem ][ 3 ];

      for( localIndex q = 0; q < numQuadraturePointsPerElem; ++q )
      {
        FE_TYPE::calcN( q, N );
        real32 const detJ = m_finiteElement.template getGradN< FE_TYPE >( k, q, xLocal, gradN );

        for( localIndex a = 0; a < numNodesPerElem; ++a )
        {
          real32 const localIncrement = invC2 * detJ * N[a];
          RAJA::atomicAdd< ATOMIC_POLICY >( &mass[elemsToNodes[k][a]], localIncrement );
        }
      }

      real32 const alpha = 1.0 / velocity[k];

      for( localIndex kfe = 0; kfe < numFacesPerElem; ++kfe )
      {
        localIndex const iface = elemsToFaces[k][kfe];

        // face on the domain boundary and not on free surface
        if( facesDomainBoundaryIndicator[iface] == 1 && freeSurfaceFaceIndicator[iface] != 1 )
        {
          for( localIndex q = 0; q < numQuadraturePointsPerElem; ++q )
          {
            FE_TYPE::calcN( q, N );
            real32 const detJ = m_finiteElement.template getGradN< FE_TYPE >( k, q, xLocal, gradN );

            real64 invJ[3][3]{};
            FE_TYPE::invJacobianTransformation( q, xLocal, invJ );

            for( localIndex a = 0; a < numNodesPerFace; ++a )
            {
              // compute ds = || detJ*invJ*normalFace_{kfe} ||

              real32 ds = 0.0;
              for( localIndex i = 0; i < 3; ++i )
              {
                real32 tmp = 0.0;
                for( localIndex j = 0; j < 3; ++j )
                {
                  tmp += invJ[j][i] * faceNormal[iface][j];
                }
                ds += tmp * tmp;
              }
              ds = sqrt( ds );

              real32 const localIncrement = alpha * detJ * ds * N[a];
              RAJA::atomicAdd< ATOMIC_POLICY >( &damping[facesToNodes[iface][a]], localIncrement );
            }
          }
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
   * Add Comments
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
    forAll< EXEC_POLICY >( size, [=] GEOSX_HOST_DEVICE ( localIndex const k )
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

      for( localIndex q=0; q<numQuadraturePointsPerElem; ++q )
      {

        real64 N[numNodesPerElem];
        real64 gradN[ numNodesPerElem ][ 3 ];

        real32 uelemx[numNodesPerElem] = {0.0};
        real32 uelemy[numNodesPerElem] = {0.0};
        real32 uelemz[numNodesPerElem] = {0.0};
        real32 flowx[numNodesPerElem] = {0.0};
        real32 flowy[numNodesPerElem] = {0.0};
        real32 flowz[numNodesPerElem] = {0.0};

        FE_TYPE::calcN( q, N );
        real32 const detJ = m_finiteElement.template getGradN< FE_TYPE >( k, q, xLocal, gradN );

        for( localIndex i = 0; i < numNodesPerElem; ++i )
        {
          uelemx[i] = detJ*velocity_x[k][i];
          uelemy[i] = detJ*velocity_y[k][i];
          uelemz[i] = detJ*velocity_z[k][i];
        }

        for( localIndex j = 0; j < numNodesPerElem; ++j )
        {
          for( localIndex i = 0; i < numNodesPerElem; ++i )
          {
            real32 dfx2 = detJ*gradN[j][0]*N[i];
            real32 dfy2 = detJ*gradN[j][1]*N[i];
            real32 dfz2 = detJ*gradN[j][2]*N[i];

            flowx[i] += dfx2*p_np1[elemsToNodes[k][j]];
            flowy[i] += dfy2*p_np1[elemsToNodes[k][j]];
            flowz[i] += dfz2*p_np1[elemsToNodes[k][j]];
          }
        }

        for( localIndex i = 0; i < numNodesPerElem; ++i )
        {
          uelemx[i]+=dt*flowx[i]/density[i];
          uelemy[i]+=dt*flowy[i]/density[i];
          uelemz[i]+=dt*flowz[i]/density[i];
        }

        for( localIndex i = 0; i < numNodesPerElem; ++i )
        {
          velocity_x[k][i] = uelemx[i]/(detJ);
          velocity_y[k][i] = uelemy[i]/(detJ);
          velocity_z[k][i] = uelemz[i]/(detJ);
        }

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
   * Add doc
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
          arrayView1d< real32 const > const mediumVelocity,
          arrayView1d< real32 const > const density,
          arrayView2d< real64 const > const sourceConstants,
          arrayView2d< real32 const > const sourceValue,
          arrayView1d< localIndex const > const sourceIsAccessible,
          arrayView1d< localIndex const > const sourceElem,
          real64 const dt,
          integer const cycleNumber,
          arrayView1d< real32 > const p_np1 )

  {

    //Pre-mult by the first factor for damping
    forAll< EXEC_POLICY >( size_node, [=] GEOSX_HOST_DEVICE ( localIndex const a )
    {
      p_np1[a] *= 1.0-((dt/2)*(damping[a]/mass[a]));
    } );

    forAll< EXEC_POLICY >( size, [=] GEOSX_HOST_DEVICE ( localIndex const k )
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


      for( localIndex q=0; q<numQuadraturePointsPerElem; ++q )
      {
        real64 N[numNodesPerElem];
        real64 gradN[ numNodesPerElem ][ 3 ];

        real32 auxx[numNodesPerElem]  = {0.0};
        real32 auyy[numNodesPerElem]  = {0.0};
        real32 auzz[numNodesPerElem]  = {0.0};
        real32 uelemx[numNodesPerElem] = {0.0};


        FE_TYPE::calcN( q, N );
        real64 const detJ = m_finiteElement.template getGradN< FE_TYPE >( k, q, xLocal, gradN );

        for( localIndex j = 0; j < numNodesPerElem; ++j )
        {
          for( localIndex i = 0; i < numNodesPerElem; ++i )
          {
            real32 dfx = detJ*gradN[i][0]*N[j];
            real32 dfy = detJ*gradN[i][1]*N[j];
            real32 dfz = detJ*gradN[i][2]*N[j];
            auxx[i] -= dfx*velocity_x[k][j];
            auyy[i] -= dfy*velocity_y[k][j];
            auzz[i] -= dfz*velocity_z[k][j];
          }

        }

        // Time update
        for( localIndex i = 0; i < numNodesPerElem; ++i )
        {
          real64 diag=(auxx[i]+auyy[i]+auzz[i]);
          uelemx[i]+=dt*diag;
        }

        for( localIndex i = 0; i < numNodesPerElem; ++i )
        {
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
                real32 const localIncrement2 = dt*(sourceConstants[isrc][i]*sourceValue[cycleNumber][isrc])/(mass[elemsToNodes[k][i]]*mediumVelocity[k]*mediumVelocity[k]*density[k]);
                RAJA::atomicAdd< ATOMIC_POLICY >( &p_np1[elemsToNodes[k][i]], localIncrement2 );
              }
            }
          }
        }

      }

    } );

    //Pre-mult by the first factor for damping
    forAll< EXEC_POLICY >( size_node, [=] GEOSX_HOST_DEVICE ( localIndex const a )
    {
      p_np1[a] *= 1.0+((dt/2)*(damping[a]/mass[a]));
    } );
  }

  /// The finite element space/discretization object for the element type in the subRegion
  FE_TYPE const & m_finiteElement;

};

} // namespace AcousticFirstOrderWaveEquationSEMKernels

} // namespace geosx

#endif //GEOSX_PHYSICSSOLVERS_WAVEPROPAGATION_AcousticFirstOrderWaveEquationSEMKERNEL_HPP_
