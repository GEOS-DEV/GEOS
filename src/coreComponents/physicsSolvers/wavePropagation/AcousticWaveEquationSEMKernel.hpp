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
 * @file AcousticWaveEquationSEMKernel.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICWAVEEQUATIONSEMKERNEL_HPP_
#define GEOSX_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICWAVEEQUATIONSEMKERNEL_HPP_

#include "finiteElement/kernelInterface/KernelBase.hpp"
#include "WaveSolverBase.hpp"


namespace geosx
{

/// Namespace to contain the acoustic wave kernels.
namespace acousticWaveEquationSEMKernels
{

struct PrecomputeSourceAndReceiverKernel
{

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
  static bool
  computeCoordinatesOnReferenceElement( real64 const (&coords)[3],
                                        real64 const (&elemCenter)[3],
                                        arraySlice1d< localIndex const, cells::NODE_MAP_USD - 1 > const elemsToNodes,
                                        arraySlice1d< localIndex const > const elemsToFaces,
                                        ArrayOfArraysView< localIndex const > const & facesToNodes,
                                        arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X,
                                        real64 (& coordsOnRefElem)[3] )
  {
    bool const isInsidePolyhedron =
      computationalGeometry::isPointInsidePolyhedron( X,
                                                      elemsToFaces,
                                                      facesToNodes,
                                                      elemCenter,
                                                      coords );
    if( isInsidePolyhedron )
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
      return true;
    }
    return false;
  }


  GEOSX_HOST_DEVICE
  static real64
  evaluateRicker( real64 const & time_n,
                  real64 const & f0,
                  localIndex const & order )
  {
    real64 const o_tpeak = 1.0/f0;
    real64 pulse = 0.0;
    if((time_n <= -0.9*o_tpeak) || (time_n >= 2.9*o_tpeak))
    {
      return pulse;
    }

    constexpr real64 pi = M_PI;
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
   * @param[out] sourceIsLocal flag indicating whether the source is local or not
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
          arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X,
          arrayView1d< integer const > const elemGhostRank,
          arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes,
          arrayView2d< localIndex const > const elemsToFaces,
          ArrayOfArraysView< localIndex const > const & facesToNodes,
          arrayView2d< real64 const > const & elemCenter,
          arrayView2d< real64 const > const sourceCoordinates,
          arrayView1d< localIndex > const sourceIsLocal,
          arrayView2d< localIndex > const sourceNodeIds,
          arrayView2d< real64 > const sourceConstants,
          arrayView2d< real64 const > const receiverCoordinates,
          arrayView1d< localIndex > const receiverIsLocal,
          arrayView2d< localIndex > const receiverNodeIds,
          arrayView2d< real64 > const receiverConstants,
          arrayView2d< real64 > const sourceValue,
          real64 const dt,
          real64 const timeSourceFrequency,
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
        if( sourceIsLocal[isrc] == 0 )
        {
          real64 const coords[3] = { sourceCoordinates[isrc][0],
                                     sourceCoordinates[isrc][1],
                                     sourceCoordinates[isrc][2] };

          real64 coordsOnRefElem[3]{};
          bool const sourceFound =
            computeCoordinatesOnReferenceElement< FE_TYPE >( coords,
                                                             center,
                                                             elemsToNodes[k],
                                                             elemsToFaces[k],
                                                             facesToNodes,
                                                             X,
                                                             coordsOnRefElem );
          if( sourceFound && elemGhostRank[k] < 0 )
          {
            sourceIsLocal[isrc] = 1;
            real64 Ntest[8];
            finiteElement::LagrangeBasis1::TensorProduct3D::value( coordsOnRefElem, Ntest );

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
            computeCoordinatesOnReferenceElement< FE_TYPE >( coords,
                                                             center,
                                                             elemsToNodes[k],
                                                             elemsToFaces[k],
                                                             facesToNodes,
                                                             X,
                                                             coordsOnRefElem );

          if( receiverFound && elemGhostRank[k] < 0 )
          {
            receiverIsLocal[ircv] = 1;

            real64 Ntest[8];
            finiteElement::LagrangeBasis1::TensorProduct3D::value( coordsOnRefElem, Ntest );
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
          arrayView1d< real64 const > const velocity,
          arrayView1d< real64 > const mass,
          arrayView1d< real64 > const damping )
  {
    forAll< EXEC_POLICY >( size, [=] GEOSX_HOST_DEVICE ( localIndex const k )
    {

      constexpr localIndex numNodesPerElem = FE_TYPE::numNodes;
      constexpr localIndex numQuadraturePointsPerElem = FE_TYPE::numQuadraturePoints;

      real64 const invC2 = 1.0 / ( velocity[k] * velocity[k] );
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
        real64 const detJ = m_finiteElement.template getGradN< FE_TYPE >( k, q, xLocal, gradN );

        for( localIndex a = 0; a < numNodesPerElem; ++a )
        {
          real64 const localIncrement = invC2 * detJ * N[a];
          RAJA::atomicAdd< ATOMIC_POLICY >( &mass[elemsToNodes[k][a]], localIncrement );
        }
      }

      real64 const alpha = 1.0 / velocity[k];

      for( localIndex kfe = 0; kfe < numFacesPerElem; ++kfe )
      {
        localIndex const iface = elemsToFaces[k][kfe];

        // face on the domain boundary and not on free surface
        if( facesDomainBoundaryIndicator[iface] == 1 && freeSurfaceFaceIndicator[iface] != 1 )
        {
          for( localIndex q = 0; q < numQuadraturePointsPerElem; ++q )
          {
            FE_TYPE::calcN( q, N );
            real64 const detJ = m_finiteElement.template getGradN< FE_TYPE >( k, q, xLocal, gradN );

            real64 invJ[3][3]{};
            FE_TYPE::invJacobianTransformation( q, xLocal, invJ );

            for( localIndex a = 0; a < numNodesPerFace; ++a )
            {
              // compute ds = || detJ*invJ*normalFace_{kfe} ||

              real64 ds = 0.0;
              for( localIndex i = 0; i < 3; ++i )
              {
                real64 tmp = 0.0;
                for( localIndex j = 0; j < 3; ++j )
                {
                  tmp += invJ[j][i] * faceNormal[iface][j];
                }
                ds += tmp * tmp;
              }
              ds = sqrt( ds );

              real64 const localIncrement = alpha * detJ * ds * N[a];
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


struct PMLKernelHelper
{
  /**
   * @brief Compute the damping profile for the Perfectly Matched Layer (PML)
   * @param xLocal a given x-y-z coordinates (3-components array)
   * @param xMin coordinate limits of the inner PML boundaries, left-front-top
   * @param xMax coordinate limits of the inner PML boundaries, right-back-bottom
   * @param dMin PML thickness, left-front-top
   * @param dMax PML thickness, right-back-bottom
   * @param cMin PML wave speed, left-front-top
   * @param cMax PML wave speed, right-back-bottom
   * @param r desired reflectivity of the PML
   * @param sigma 3-components array to hold the damping profile in each direction
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static void computeDampingProfilePML( real64 const (&xLocal)[3],
                                        real64 const (&xMin)[3],
                                        real64 const (&xMax)[3],
                                        real64 const (&dMin)[3],
                                        real64 const (&dMax)[3],
                                        real64 const (&cMin)[3],
                                        real64 const (&cMax)[3],
                                        real64 const r,
                                        real64 (&sigma)[3])
  {

    sigma[0] = 0;
    sigma[1] = 0;
    sigma[2] = 0;
    
    if (xLocal[0] < xMin[0])
    {
      real64 const factor =  -3.0/2.0*cMin[0]*log(r)/(dMin[0]*dMin[0]*dMin[0]);
      sigma[0] = factor*(xLocal[0]-xMin[0])*(xLocal[0]-xMin[0]);
    }
    else if (xLocal[0] > xMax[0])
    {
      real64 const factor =  -3.0/2.0*cMax[0]*log(r)/(dMax[0]*dMax[0]*dMax[0]);
      sigma[0] = factor*(xLocal[0]-xMax[0])*(xLocal[0]-xMax[0]);
    }
    if (xLocal[1] < xMin[1])
    {
      real64 const factor =  -3.0/2.0*cMin[1]*log(r)/(dMin[1]*dMin[1]*dMin[1]);
      sigma[1] = factor*(xLocal[1]-xMin[1])*(xLocal[1]-xMin[1]);
    }
    else if (xLocal[1] > xMax[1])
    {
      real64 const factor =  -3.0/2.0*cMax[1]*log(r)/(dMax[1]*dMax[1]*dMax[1]);
      sigma[1] = factor*(xLocal[1]-xMax[1])*(xLocal[1]-xMax[1]);
    }
    if (xLocal[2] < xMin[2])
    {
      real64 const factor =  -3.0/2.0*cMin[2]*log(r)/(dMin[2]*dMin[2]*dMin[2]);
      sigma[2] = factor*(xLocal[2]-xMin[2])*(xLocal[2]-xMin[2]);
    }
    else if (xLocal[2] > xMax[2])
    {
      real64 const factor =  -3.0/2.0*cMax[2]*log(r)/(dMax[2]*dMax[2]*dMax[2]);
      sigma[2] = factor*(xLocal[2]-xMax[2])*(xLocal[2]-xMax[2]);
    }
  }
};

template< typename FE_TYPE >
struct PMLKernel
{

  PMLKernel( FE_TYPE const & finiteElement )
    : m_finiteElement( finiteElement )
  {}

  /**
   * @brief Launches the computation of field gradients and divergence for PML region
   * @tparam EXEC_POLICY the execution policy
   * @tparam ATOMIC_POLICY the atomic policy
   * @param[in] targetSet list of cells in the target set
   * @param[in] X coordinates of the nodes
   * @param[in] elemToNodesViewConst constant array view of map from element to nodes
   * @param[in] velocity cell-wise velocity
   * @param[in] p_n pressure field at time n
   * @param[in] v_n PML auxiliary field at time n
   * @param[in] u_n PML auxiliary field at time n
   * @param[in] xMin coordinate limits of the inner PML boundaries, left-front-top
   * @param[in] xMax coordinate limits of the inner PML boundaries, right-back-bottom
   * @param[in] dMin PML thickness, left-front-top
   * @param[in] dMax PML thickness, right-back-bottom
   * @param[in] cMin PML wave speed, left-front-top
   * @param[in] cMax PML wave speed, right-back-bottom
   * @param[in] r desired reflectivity of the PML
   * @param[in] flagPML switch to flip between PML formulations
   * @param[out] grad_n array holding the gradients at time n
   * @param[out] divV_n array holding the divergence at time n
   */
  template< typename EXEC_POLICY, typename ATOMIC_POLICY >
  void
  launch( SortedArrayView< localIndex const > const targetSet,
          arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X,
          traits::ViewTypeConst< CellElementSubRegion::NodeMapType > const elemToNodesViewConst,
          arrayView1d< real64 const > const velocity,
          arrayView1d< real64 const > const p_n,
          arrayView2d< real64 const > const v_n,
          arrayView1d< real64 const > const u_n,
          real64 const (&xMin)[3],
          real64 const (&xMax)[3],
          real64 const (&dMin)[3],
          real64 const (&dMax)[3],
          real64 const (&cMin)[3],
          real64 const (&cMax)[3],
          real64 const r,
          int const flagPML,
          arrayView2d< real64 > const grad_n,
          arrayView1d< real64 > const divV_n )
  {
    
    /// Loop over elements in the subregion, 'l' is the element index within the target set
    forAll< EXEC_POLICY >( targetSet.size(), [=] GEOSX_HOST_DEVICE ( localIndex const l )
    {
      constexpr localIndex numNodesPerElem = FE_TYPE::numNodes;

      /// global element index
      localIndex const k = targetSet[l];

      /// wave speed at the element
      real64 const c = velocity[k];

      /// coordinates of the element nodes
      real64 xLocal[ numNodesPerElem ][ 3 ];

      /// local arrays to store the pressure at all nodes and its gradient at a given node
      real64 pressure[ numNodesPerElem ];
      real64 pressureGrad[ 3 ];

      /// local arrays to store the PML vectorial auxiliary variable at all nodes and its gradient at a given node
      real64 auxV[3][ numNodesPerElem ];
      real64 auxVGrad[3][3];

      /// local arrays to store the PML scalar auxiliary variable at all nodes and its gradient at a given node
      real64 auxU[ numNodesPerElem ];
      real64 auxUGrad[3];

      /// local array to store the PML damping profile
      real64 sigma[ 3 ];

      /// copy from global to local arrays
      for( localIndex i=0; i<numNodesPerElem; ++i )
      {
        pressure[i] = p_n[elemToNodesViewConst[k][i]];
        auxU[i] = u_n[elemToNodesViewConst[k][i]];
        for( int j=0; j<3; ++j )
        {
          xLocal[i][j] =  X[elemToNodesViewConst[k][i]][j];
          auxV[j][i] = v_n[elemToNodesViewConst[k][i]][j];
        }
      }

      if (flagPML>1)
      {
        for( localIndex i=0; i<numNodesPerElem; ++i )
        {
          /// compute the PML damping profile
          PMLKernelHelper::computeDampingProfilePML( 
            xLocal[i],
            xMin,
            xMax,
            dMin,
            dMax,
            cMin,
            cMax,
            r,
            sigma);
          for( int j=0; j<3; ++j )
          {
            auxV[j][i] *= sigma[j];
          }
        }
      }

      /// local arrays to store shape functions gradients
      real64 gradN[ numNodesPerElem ][ 3 ];
      using GRADIENT_TYPE = TYPEOFREF( gradN );
    
      /// loop over the nodes i in the element k
      /// the nodes are implicitly assumed the same as quadrature points
      for( localIndex i=0; i<numNodesPerElem; ++i )
      {
        
        /// compute the shape functions gradients
        real64 const detJ = m_finiteElement.template getGradN< FE_TYPE >( k, i, xLocal, gradN );
        GEOSX_UNUSED_VAR (detJ);

        /// compute the gradient of the pressure and the PML auxiliary variables at the node
        m_finiteElement.template gradient< numNodesPerElem, GRADIENT_TYPE >(gradN, pressure, pressureGrad );
        m_finiteElement.template gradient< numNodesPerElem, GRADIENT_TYPE >(gradN, auxU, auxUGrad );
        for( int j=0; j<3; ++j )
        {
          m_finiteElement.template gradient< numNodesPerElem, GRADIENT_TYPE >(gradN, auxV[j], auxVGrad[j] );
        }

        if (flagPML==1)
        {
          /// compute the PML damping profile
          PMLKernelHelper::computeDampingProfilePML( 
            xLocal[i],
            xMin,
            xMax,
            dMin,
            dMax,
            cMin,
            cMax,
            r,
            sigma);

          /// compute B.pressureGrad - C.auxUGrad where B and C are functions of the damping profile
          real64 localIncrementArray[3];
          localIncrementArray[0] = (sigma[0]-sigma[1]-sigma[2])*pressureGrad[0] - (sigma[1]*sigma[2])*auxUGrad[0];
          localIncrementArray[1] = (sigma[1]-sigma[0]-sigma[2])*pressureGrad[1] - (sigma[0]*sigma[2])*auxUGrad[1];
          localIncrementArray[2] = (sigma[2]-sigma[0]-sigma[1])*pressureGrad[2] - (sigma[0]*sigma[1])*auxUGrad[2];
          for (int j=0; j<3; ++j)
          {
            RAJA::atomicAdd< ATOMIC_POLICY >( &grad_n[elemToNodesViewConst[k][i]][j], localIncrementArray[j]/numNodesPerElem );
          }
          /// compute beta.pressure + gamma.u - c^2 * divV where beta and gamma are functions of the damping profile
          real64 const beta = sigma[0]*sigma[1]+sigma[0]*sigma[2]+sigma[1]*sigma[2];
          real64 const gamma = sigma[0]*sigma[1]*sigma[2];
          real64 const localIncrement = beta*p_n[elemToNodesViewConst[k][i]]
                                      + gamma*u_n[elemToNodesViewConst[k][i]]
                                      - c*c*( auxVGrad[0][0] + auxVGrad[1][1] + auxVGrad[2][2] );

          RAJA::atomicAdd< ATOMIC_POLICY >( &divV_n[elemToNodesViewConst[k][i]], localIncrement/numNodesPerElem);
        }
        else
        {
          for (int j=0; j<3; ++j)
          {
            RAJA::atomicAdd< ATOMIC_POLICY >( &grad_n[elemToNodesViewConst[k][i]][j], pressureGrad[j]/numNodesPerElem );
          }
          real64 const localIncrement = -c*c*(auxVGrad[0][0] + auxVGrad[1][1] + auxVGrad[2][2]);
          RAJA::atomicAdd< ATOMIC_POLICY >( &divV_n[elemToNodesViewConst[k][i]], localIncrement/numNodesPerElem);
        }
      }
    } );
  }

  /// The finite element space/discretization object for the element type in the subRegion
  FE_TYPE const & m_finiteElement;
};


template< typename FE_TYPE >
struct waveSpeedPMLKernel
{

  waveSpeedPMLKernel( FE_TYPE const & finiteElement )
    : m_finiteElement( finiteElement )
  {}

  /**
   * @brief Launches the computation of average wave speeds in the PML region
   * @tparam EXEC_POLICY the execution policy
   * @tparam ATOMIC_POLICY the atomic policy
   * @param[in] targetSet list of cells in the target set
   * @param[in] X coordinates of the nodes
   * @param[in] elemToNodesViewConst constant array view of map from element to nodes
   * @param[in] velocity cell-wise velocity
   * @param[in] xMin coordinate limits of the inner PML boundaries, left-front-top
   * @param[in] xMax coordinate limits of the inner PML boundaries, right-back-bottom
   * @param[out] cMin PML wave speed, left-front-top
   * @param[out] cMax PML wave speed, right-back-bottom
   * @param[out] counterMin PML wave speed counter, left-front-top
   * @param[out] counterMax PML wave speed counter, left-front-top
   */
  template< typename EXEC_POLICY, typename ATOMIC_POLICY >
  void
  launch( SortedArrayView< localIndex const > const targetSet,
          arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X,
          traits::ViewTypeConst< CellElementSubRegion::NodeMapType > const elemToNodesViewConst,
          arrayView1d< real64 const > const velocity,
          real64 const (&xMin)[3],
          real64 const (&xMax)[3],
          real64 (&cMin)[3],
          real64 (&cMax)[3],
          int (&counterMin)[3],
          int (&counterMax)[3])
  {
    
    RAJA::ReduceSum< parallelDeviceReduce, real64 > subRegionAvgWaveSpeedLeft( 0.0 );
    RAJA::ReduceSum< parallelDeviceReduce, real64 > subRegionAvgWaveSpeedRight( 0.0 );
    RAJA::ReduceSum< parallelDeviceReduce, real64 > subRegionAvgWaveSpeedFront( 0.0 );
    RAJA::ReduceSum< parallelDeviceReduce, real64 > subRegionAvgWaveSpeedBack( 0.0 );
    RAJA::ReduceSum< parallelDeviceReduce, real64 > subRegionAvgWaveSpeedTop( 0.0 );
    RAJA::ReduceSum< parallelDeviceReduce, real64 > subRegionAvgWaveSpeedBottom( 0.0 );
    RAJA::ReduceSum< parallelDeviceReduce, int > subRegionAvgWaveSpeedCounterLeft( 0 );
    RAJA::ReduceSum< parallelDeviceReduce, int > subRegionAvgWaveSpeedCounterRight( 0 );
    RAJA::ReduceSum< parallelDeviceReduce, int > subRegionAvgWaveSpeedCounterFront( 0 );
    RAJA::ReduceSum< parallelDeviceReduce, int > subRegionAvgWaveSpeedCounterBack( 0 );
    RAJA::ReduceSum< parallelDeviceReduce, int > subRegionAvgWaveSpeedCounterTop( 0 );
    RAJA::ReduceSum< parallelDeviceReduce, int > subRegionAvgWaveSpeedCounterBottom( 0 );

    /// Loop over elements in the subregion, 'l' is the element index within the target set
    forAll< EXEC_POLICY >( targetSet.size(), [=] GEOSX_HOST_DEVICE ( localIndex const l )
    {
      constexpr localIndex numNodesPerElem = FE_TYPE::numNodes;

      /// global element index
      localIndex const k = targetSet[l];

      /// wave speed at the element
      real64 const c = velocity[k];

      /// coordinates of the element center
      real64 xLocal[ 3 ] = {0.0, 0.0, 0.0};

      /// compute the coordinates of the element center
      for( int j=0; j<3; ++j )
      {
        for( localIndex i=0; i<numNodesPerElem; ++i )
        {
          xLocal[j] +=  X[elemToNodesViewConst[k][i]][j];
        }
        xLocal[j] /= numNodesPerElem;
      }

      if ( xLocal[0] < xMin[0] 
           && xLocal[1] >= xMin[1] && xLocal[1] <= xMax[1]
           && xLocal[2] >= xMin[2] && xLocal[2] <= xMax[2] )
      {
        subRegionAvgWaveSpeedLeft += c;
        subRegionAvgWaveSpeedCounterLeft += 1;
      }
      else if ( xLocal[0] > xMax[0]
                && xLocal[1] >= xMin[1] && xLocal[1] <= xMax[1]
                && xLocal[2] >= xMin[2] && xLocal[2] <= xMax[2] )
      {
        subRegionAvgWaveSpeedRight += c;
        subRegionAvgWaveSpeedCounterRight += 1;
      }
      if ( xLocal[1] < xMin[1]
           && xLocal[0] >= xMin[0] && xLocal[0] <= xMax[0]
           && xLocal[2] >= xMin[2] && xLocal[2] <= xMax[2] )
      {
        subRegionAvgWaveSpeedFront += c;
        subRegionAvgWaveSpeedCounterFront += 1;
      }
      else if ( xLocal[1] > xMax[1]
                && xLocal[0] >= xMin[0] && xLocal[0] <= xMax[0]
                && xLocal[2] >= xMin[2] && xLocal[2] <= xMax[2] )
      {
        subRegionAvgWaveSpeedBack += c;
        subRegionAvgWaveSpeedCounterBack += 1;
      }
      if ( xLocal[2] < xMin[2]
           && xLocal[0] >= xMin[0] && xLocal[0] <= xMax[0]
           && xLocal[1] >= xMin[1] && xLocal[1] <= xMax[1] )
      {
        subRegionAvgWaveSpeedTop += c;
        subRegionAvgWaveSpeedCounterTop += 1;
      }
      else if ( xLocal[2] > xMax[2]
                && xLocal[0] >= xMin[0] && xLocal[0] <= xMax[0]
                && xLocal[1] >= xMin[1] && xLocal[1] <= xMax[1] )
      {
        subRegionAvgWaveSpeedBottom += c;
        subRegionAvgWaveSpeedCounterBottom += 1;
      }
    } );

    cMin[0]+=subRegionAvgWaveSpeedLeft.get();
    cMin[1]+=subRegionAvgWaveSpeedFront.get();
    cMin[2]+=subRegionAvgWaveSpeedTop.get();
    cMax[0]+=subRegionAvgWaveSpeedRight.get();
    cMax[1]+=subRegionAvgWaveSpeedBack.get();
    cMax[2]+=subRegionAvgWaveSpeedBottom.get();
    counterMin[0]+=subRegionAvgWaveSpeedCounterLeft.get();
    counterMin[1]+=subRegionAvgWaveSpeedCounterFront.get();
    counterMin[2]+=subRegionAvgWaveSpeedCounterTop.get();
    counterMax[0]+=subRegionAvgWaveSpeedCounterRight.get();
    counterMax[1]+=subRegionAvgWaveSpeedCounterBack.get();
    counterMax[2]+=subRegionAvgWaveSpeedCounterBottom.get();
  }

  /// The finite element space/discretization object for the element type in the subRegion
  FE_TYPE const & m_finiteElement;
};




/**
 * @brief Implements kernels for solving the acoustic wave equations
 *   explicit central FD method and SEM
 * @copydoc geosx::finiteElement::KernelBase
 * @tparam SUBREGION_TYPE The type of subregion that the kernel will act on.
 *
 * ### AcousticWaveEquationSEMKernel Description
 * Implements the KernelBase interface functions required for solving
 * the acoustic wave equations using the
 * "finite element kernel application" functions such as
 * geosx::finiteElement::RegionBasedKernelApplication.
 *
 * The number of degrees of freedom per support point for both
 * the test and trial spaces are specified as `1`.
 */


template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
class ExplicitAcousticSEM : public finiteElement::KernelBase< SUBREGION_TYPE,
                                                              CONSTITUTIVE_TYPE,
                                                              FE_TYPE,
                                                              1,
                                                              1 >
{
public:

  /// Alias for the base class;
  using Base = finiteElement::KernelBase< SUBREGION_TYPE,
                                          CONSTITUTIVE_TYPE,
                                          FE_TYPE,
                                          1,
                                          1 >;

  /// Maximum number of nodes per element, which is equal to the maxNumTestSupportPointPerElem and
  /// maxNumTrialSupportPointPerElem by definition. When the FE_TYPE is not a Virtual Element, this
  /// will be the actual number of nodes per element.
  static constexpr int numNodesPerElem = Base::maxNumTestSupportPointsPerElem;

  using Base::numDofPerTestSupportPoint;
  using Base::numDofPerTrialSupportPoint;
  using Base::m_elemsToNodes;
  using Base::m_elemGhostRank;
  using Base::m_constitutiveUpdate;
  using Base::m_finiteElementSpace;

//*****************************************************************************
  /**
   * @brief Constructor
   * @copydoc geosx::finiteElement::KernelBase::KernelBase
   * @param nodeManager Reference to the NodeManager object.
   * @param edgeManager Reference to the EdgeManager object.
   * @param faceManager Reference to the FaceManager object.
   * @param targetRegionIndex Index of the region the subregion belongs to.
   * @param dt The time interval for the step.
   *   elements to be processed during this kernel launch.
   */
  ExplicitAcousticSEM( NodeManager & nodeManager,
                       EdgeManager const & edgeManager,
                       FaceManager const & faceManager,
                       localIndex const targetRegionIndex,
                       SUBREGION_TYPE const & elementSubRegion,
                       FE_TYPE const & finiteElementSpace,
                       CONSTITUTIVE_TYPE & inputConstitutiveType,
                       real64 const dt ):
    Base( elementSubRegion,
          finiteElementSpace,
          inputConstitutiveType ),
    m_X( nodeManager.referencePosition() ),
    m_p_n( nodeManager.getExtrinsicData< extrinsicMeshData::Pressure_n >() ),
    m_stiffnessVector( nodeManager.getExtrinsicData< extrinsicMeshData::StiffnessVector >() ),
    m_dt( dt )
  {
    GEOSX_UNUSED_VAR( edgeManager );
    GEOSX_UNUSED_VAR( faceManager );
    GEOSX_UNUSED_VAR( targetRegionIndex );
  }



  //*****************************************************************************
  /**
   * @copydoc geosx::finiteElement::KernelBase::StackVariables
   *
   * ### ExplicitAcousticSEM Description
   * Adds a stack arrays for the nodal force, primary displacement variable, etc.
   */
  struct StackVariables : Base::StackVariables
  {
public:
    GEOSX_HOST_DEVICE
    StackVariables():
      xLocal()
    {}

    /// C-array stack storage for element local the nodal positions.
    real64 xLocal[ numNodesPerElem ][ 3 ];
  };
  //***************************************************************************


  /**
   * @copydoc geosx::finiteElement::KernelBase::setup
   *
   * Copies the primary variable, and position into the local stack array.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void setup( localIndex const k,
              StackVariables & stack ) const
  {
    /// numDofPerTrialSupportPoint = 1
    for( localIndex a=0; a< numNodesPerElem; ++a )
    {
      localIndex const nodeIndex = m_elemsToNodes( k, a );
      for( int i=0; i< 3; ++i )
      {
        stack.xLocal[ a ][ i ] = m_X[ nodeIndex ][ i ];
      }
    }
  }

  /**
   * @copydoc geosx::finiteElement::KernelBase::quadraturePointKernel
   *
   * ### ExplicitAcousticSEM Description
   * Calculates stiffness vector
   *
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void quadraturePointKernel( localIndex const k,
                              localIndex const q,
                              StackVariables & stack ) const
  {
    real64 gradN[ numNodesPerElem ][ 3 ];

    real64 const detJ = m_finiteElementSpace.template getGradN< FE_TYPE >( k, q, stack.xLocal, gradN );

    for( localIndex i=0; i<numNodesPerElem; ++i )
    {
      for( localIndex j=0; j<numNodesPerElem; ++j )
      {
        real64 const Rh_ij = detJ * LvArray::tensorOps::AiBi< 3 >( gradN[ i ], gradN[ j ] );
        real64 const localIncrement = Rh_ij*m_p_n[m_elemsToNodes[k][j]];

        RAJA::atomicAdd< parallelDeviceAtomic >( &m_stiffnessVector[m_elemsToNodes[k][i]], localIncrement );
      }
    }
  }


protected:
  /// The array containing the nodal position array.
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const m_X;

  /// The array containing the nodal pressure array.
  arrayView1d< real64 const > const m_p_n;

  /// The array containing the product of the stiffness matrix and the nodal pressure.
  arrayView1d< real64 > const m_stiffnessVector;

  /// The time increment for this time integration step.
  real64 const m_dt;


};



/// The factory used to construct a ExplicitAcousticWaveEquation kernel.
using ExplicitAcousticSEMFactory = finiteElement::KernelFactory< ExplicitAcousticSEM,
                                                                 real64 >;


} // namespace acousticWaveEquationSEMKernels

} // namespace geosx

#endif //GEOSX_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICWAVEEQUATIONSEMKERNEL_HPP_
