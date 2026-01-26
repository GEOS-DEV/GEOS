/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file AcousticMatricesSEMKernel.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICMATRICESSEMKERNEL_HPP_
#define GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICMATRICESSEMKERNEL_HPP_

#include "finiteElement/elementFormulations/Qk_Hexahedron_Lagrange_GaussLobatto.hpp"


namespace geos
{
struct AcousticMatricesSEM
{
  // Debug
  template< typename FE_TYPE >

  struct DofArrays
  {
    DofArrays( FE_TYPE const & finiteElement ): m_finiteElement( finiteElement ) {}

    //==============================================================================
    // Public API
    //==============================================================================
public:

    template< typename EXEC_POLICY, typename ATOMIC_POLICY >
    void
    computeDofArraysVTI( localIndex const size,
                         arrayView2d< WaveSolverBase::wsCoordType const, nodes::REFERENCE_POSITION_USD > const nodeCoords,
                         arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodes,
                         arrayView1d< real32 const > const vtiEpsilon,
                         arrayView1d< real32 const > const vtiDelta,
                         arrayView1d< real32 > const dofEpsilon,
                         arrayView1d< real32 > const dofDelta,
                         arrayView1d< real32 > const dofOrder,
                         arrayView2d< real64 const > const sourceCoordinates,
                         real64 const radiusIsoAroundSource )
    {
      zero_out_array< EXEC_POLICY >( dofEpsilon );
      zero_out_array< EXEC_POLICY >( dofDelta );
      zero_out_array< EXEC_POLICY >( dofOrder );

      computeDofOrder< EXEC_POLICY, ATOMIC_POLICY >( size, elemsToNodes, dofOrder );

      computeCornerNodeValuesVTI< EXEC_POLICY, ATOMIC_POLICY >( size, elemsToNodes, vtiEpsilon, vtiDelta,
                                                                dofEpsilon, dofDelta, dofOrder );

      applySourceTapering< EXEC_POLICY >( nodeCoords, sourceCoordinates, radiusIsoAroundSource,
                                          dofEpsilon, dofDelta );

      if( FE_TYPE::numNodes > 8 ) // A Hex8 element has 8 nodes. Anything more is high-order.
      {
        interpolateInternalGLLNodes< EXEC_POLICY, ATOMIC_POLICY >( size, nodeCoords, elemsToNodes, dofOrder,
                                                                   dofEpsilon, dofDelta );
      }
    }

    template< typename EXEC_POLICY, typename ATOMIC_POLICY >
    void
    computeDofArraysTTI( localIndex const size,
                         arrayView2d< WaveSolverBase::wsCoordType const, nodes::REFERENCE_POSITION_USD > const nodeCoords,
                         arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodes,
                         arrayView1d< real32 const > const tti_dipx,
                         arrayView1d< real32 const > const tti_dipy,
                         arrayView1d< real32 > const dofTilt,
                         arrayView1d< real32 > const dofAzimuth,
                         arrayView1d< real32 > const dofOrder,
                         arrayView2d< real64 const > const sourceCoordinates,
                         real64 const radiusIsoAroundSource )
    {
      zero_out_array< EXEC_POLICY >( dofTilt );
      zero_out_array< EXEC_POLICY >( dofAzimuth );

      computeCornerNodeValuesTTI< EXEC_POLICY, ATOMIC_POLICY >( size, elemsToNodes, tti_dipx, tti_dipy,
                                                                dofTilt, dofAzimuth, dofOrder );

      applySourceTapering< EXEC_POLICY >( nodeCoords, sourceCoordinates, radiusIsoAroundSource,
                                          dofTilt, dofAzimuth );

      if( FE_TYPE::numNodes > 8 ) // A Hex8 element has 8 nodes. Anything more is high-order.
      {
        interpolateInternalGLLNodes< EXEC_POLICY, ATOMIC_POLICY >( size, nodeCoords, elemsToNodes, dofOrder,
                                                                   dofTilt, dofAzimuth );
      }
    }


public:

    template< typename EXEC_POLICY, typename ViewType >
    void zero_out_array( ViewType const & view ) const
    {
      const localIndex size = view.size();
      forAll< EXEC_POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const i )
      {
        view[i] = 0;
      } );
    }

    template< typename EXEC_POLICY, typename ATOMIC_POLICY >
    void
    computeDofOrder( localIndex const size,
                     arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodes,
                     arrayView1d< real32 > const dofOrder ) const
    {
      forAll< EXEC_POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const e )
      {
        constexpr localIndex numNodesPerElem = FE_TYPE::numNodes;
        for( localIndex a = 0; a < numNodesPerElem; ++a )
        {
          localIndex const nodeIndex = elemsToNodes( e, a );
          RAJA::atomicAdd< ATOMIC_POLICY >( &dofOrder[nodeIndex], 1.0f );
        }
      } );
    }

    template< typename EXEC_POLICY, typename ATOMIC_POLICY >
    void
    computeCornerNodeValuesVTI( localIndex const size,
                                arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodes,
                                arrayView1d< real32 const > const vtiEpsilon,
                                arrayView1d< real32 const > const vtiDelta,
                                arrayView1d< real32 > const dofEpsilon,
                                arrayView1d< real32 > const dofDelta,
                                arrayView1d< real32 > const dofOrder ) const
    {
      forAll< EXEC_POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const e )
      {
        real32 localEpsilon = vtiEpsilon[e];
        real32 localDelta = vtiDelta[e];
        if( localEpsilon < 1e-5f ) localEpsilon = 0.f;
        if( localDelta < 1e-5f ) localDelta = 0.f;
        if( localDelta > localEpsilon ) localDelta = localEpsilon;

        constexpr localIndex numCorners = 8;
        for( localIndex a = 0; a < numCorners; ++a )
        {
          localIndex const cornerNodeIndex = FE_TYPE::meshIndexToLinearIndex3D( a );
          localIndex const globalNodeIndex = elemsToNodes( e, cornerNodeIndex );
          real32 const invOrder = 1.f / dofOrder[globalNodeIndex];

          RAJA::atomicAdd< ATOMIC_POLICY >( &dofEpsilon[globalNodeIndex], localEpsilon * invOrder );
          RAJA::atomicAdd< ATOMIC_POLICY >( &dofDelta[globalNodeIndex], localDelta * invOrder );
        }
      } );
    }

    template< typename EXEC_POLICY, typename ATOMIC_POLICY >
    void
    computeCornerNodeValuesTTI( localIndex const size,
                                arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodes,
                                arrayView1d< real32 const > const tti_dipx,
                                arrayView1d< real32 const > const tti_dipy,
                                arrayView1d< real32 > const dofTilt,
                                arrayView1d< real32 > const dofAzimuth,
                                arrayView1d< real32 > const dofOrder ) const
    {
      forAll< EXEC_POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const e )
      {
        real32 const dipx = tti_dipx[e];
        real32 const dipy = tti_dipy[e];
        real32 const deg_to_rad = geos::constants::pi / 180.f;
        real32 const localTilt = std::atan( std::sqrt( dipx * dipx + dipy * dipy ) );

        real32 localAzimuth = 0.f;
        real32 const ftmp = std::atan2( dipy, dipx );
        if( localTilt >= ( 0.001f * deg_to_rad ) )
        {
          localAzimuth = ( ftmp <= 0.f ) ? ( ftmp + geos::constants::pi )
                                         : ( ftmp - geos::constants::pi );
          if( localAzimuth < 0. ) localAzimuth += 2 * geos::constants::pi;
        }

        constexpr localIndex numCorners = 8;
        for( localIndex a = 0; a < numCorners; ++a )
        {
          localIndex const cornerNodeIndex = FE_TYPE::meshIndexToLinearIndex3D( a );
          localIndex const globalNodeIndex = elemsToNodes( e, cornerNodeIndex );
          real32 const invOrder = 1.f / dofOrder[globalNodeIndex];
          RAJA::atomicAdd< ATOMIC_POLICY >( &dofTilt[globalNodeIndex], localTilt * invOrder );
          RAJA::atomicAdd< ATOMIC_POLICY >( &dofAzimuth[globalNodeIndex], localAzimuth * invOrder );
        }
      } );
    }

    template< typename EXEC_POLICY, typename View1, typename View2 >
    void
    applySourceTapering( arrayView2d< WaveSolverBase::wsCoordType const, nodes::REFERENCE_POSITION_USD > const nodeCoords,
                         arrayView2d< real64 const > const sourceCoordinates,
                         real64 const radiusIsoAroundSource,
                         View1 const dofField1,
                         View2 const dofField2 ) const
    {
      const localIndex numNodes = nodeCoords.size( 0 );
      forAll< EXEC_POLICY >( numNodes, [=] GEOS_HOST_DEVICE ( localIndex const nodeIndex )
      {
        real64 minDistSq = -1.0;

        for( localIndex isrc = 0; isrc < sourceCoordinates.size( 0 ); ++isrc )
        {
          real64 const dx = nodeCoords( nodeIndex, 0 ) - sourceCoordinates[isrc][0];
          real64 const dy = nodeCoords( nodeIndex, 1 ) - sourceCoordinates[isrc][1];
          real64 const dz = nodeCoords( nodeIndex, 2 ) - sourceCoordinates[isrc][2];
          real64 const distSq = dx * dx + dy * dy + dz * dz;
          if( ( minDistSq < 0.0 ) || ( distSq < minDistSq ) )
          {
            minDistSq = distSq;
          }
        }

        real64 const radiusSq = radiusIsoAroundSource * radiusIsoAroundSource;
        if( ( minDistSq >= 0.0 ) && ( minDistSq < radiusSq ) )
        {
          real64 const minDist = std::sqrt( minDistSq );
          real32 const taper = static_cast< real32 >( cosine_taper( minDist, 0.0, radiusIsoAroundSource ) );
          auto const taper_squared = taper * taper;

          dofField1[nodeIndex] *= taper_squared * taper_squared;
          dofField2[nodeIndex] *= taper_squared * taper_squared;
        }
      } );
    }

    template< typename EXEC_POLICY, typename ATOMIC_POLICY, typename View1, typename View2 >
    void
    interpolateInternalGLLNodes( localIndex const size,
                                 arrayView2d< WaveSolverBase::wsCoordType const, nodes::REFERENCE_POSITION_USD > const nodeCoords,
                                 arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodes,
                                 arrayView1d< real32 > const dofOrder,
                                 View1 const dofField1,
                                 View2 const dofField2 ) const
    {
      forAll< EXEC_POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const e )
      {
        constexpr localIndex numCorners = 8;

        real64 coordCorners[numCorners][3];
        real32 vals1Corners[numCorners];
        real32 vals2Corners[numCorners];

        for( localIndex a = 0; a < numCorners; ++a )
        {
          localIndex const cornerNodeIndex = FE_TYPE::meshIndexToLinearIndex3D( a );
          localIndex const globalNodeIndex = elemsToNodes( e, cornerNodeIndex );

          coordCorners[a][0] = nodeCoords( globalNodeIndex, 0 );
          coordCorners[a][1] = nodeCoords( globalNodeIndex, 1 );
          coordCorners[a][2] = nodeCoords( globalNodeIndex, 2 );

          vals1Corners[a] = dofField1[globalNodeIndex];
          vals2Corners[a] = dofField2[globalNodeIndex];
        }

        constexpr localIndex numNodesPerElem = FE_TYPE::numNodes;
        for( localIndex q = 0; q < numNodesPerElem; ++q )
        {
          if( !is_corner_node( q ) )
          {
            localIndex const globalNodeIndex = elemsToNodes( e, q );
            real64 coordTarget[3] = { nodeCoords( globalNodeIndex, 0 ), nodeCoords( globalNodeIndex, 1 ), nodeCoords( globalNodeIndex, 2 ) };
            real32 const invOrder = 1.f / dofOrder[globalNodeIndex];

            real32 const interp_val1 = trilinear_interpolate( coordTarget, coordCorners, vals1Corners );
            RAJA::atomicAdd< ATOMIC_POLICY >( &dofField1[globalNodeIndex], interp_val1 * invOrder );

            real32 const interp_val2 = trilinear_interpolate( coordTarget, coordCorners, vals2Corners );
            RAJA::atomicAdd< ATOMIC_POLICY >( &dofField2[globalNodeIndex], interp_val2 * invOrder );
          }
        }
      } );
    }

    //==============================================================================
    // Low-level static helpers
    //==============================================================================

    GEOS_HOST_DEVICE static bool is_corner_node( localIndex const q )
    {
      for( localIndex a = 0; a < 8; ++a )
      {
        if( FE_TYPE::meshIndexToLinearIndex3D( a ) == q )
        {
          return true;
        }
      }
      return false;
    }

    GEOS_HOST_DEVICE static real32
    trilinear_interpolate( real64 const (&coordTarget)[3],
                           real64 const (&coordCorners)[8][3],
                           real32 const (&valCorners)[8] )
    {
      // This assumes the 8 corners are ordered logically like a unit cube.
      real64 const xd = ( coordTarget[0] - coordCorners[0][0] ) / ( coordCorners[7][0] - coordCorners[0][0] );
      real64 const yd = ( coordTarget[1] - coordCorners[0][1] ) / ( coordCorners[7][1] - coordCorners[0][1] );
      real64 const zd = ( coordTarget[2] - coordCorners[0][2] ) / ( coordCorners[7][2] - coordCorners[0][2] );

      real32 const c00 = valCorners[0] * ( 1.f - xd ) + valCorners[1] * xd;
      real32 const c10 = valCorners[2] * ( 1.f - xd ) + valCorners[3] * xd;
      real32 const c01 = valCorners[4] * ( 1.f - xd ) + valCorners[5] * xd;
      real32 const c11 = valCorners[6] * ( 1.f - xd ) + valCorners[7] * xd;

      real32 const c0 = c00 * ( 1.f - yd ) + c10 * yd;
      real32 const c1 = c01 * ( 1.f - yd ) + c11 * yd;

      return c0 * ( 1.f - zd ) + c1 * zd;
    }

    GEOS_HOST_DEVICE double cosine_taper( double r, double rmin, double rmax ) const
    {
      double const expr = ( rmax - rmin > 1e-9 ) ? ( ( r - rmin ) / ( rmax - rmin ) ) : 0.0;
      double const t = (expr < 0.0) ? 0.0 : ( (1.0 < expr) ? 1.0 : expr );
      return 0.5 * ( 1.0 - std::cos( geos::constants::pi * t ) );
    }

    // Member data
    FE_TYPE const & m_finiteElement;
  };

  // End debug

  template< typename FE_TYPE >
  struct MassMatrix
  {

    MassMatrix( FE_TYPE const & finiteElement )
      : m_finiteElement( finiteElement )
    {}
    /**
     * @brief Launches the precomputation of the mass matrices
     * @tparam EXEC_POLICY the execution policy
     * @tparam ATOMIC_POLICY the atomic policy
     * @tparam FE_TYPE the type of discretization
     * @param[in] finiteElement The finite element discretization used
     * @param[in] size the number of cells in the subRegion
     * @param[in] numFacesPerElem number of faces per element
     * @param[in] nodeCoords coordinates of the nodes
     * @param[in] elemsToNodes map from element to nodes
     * @param[in] velocity cell-wise velocity
     * @param[in] density cell-wise density
     * @param[out] mass diagonal of the mass matrix
     */
    template< typename EXEC_POLICY, typename ATOMIC_POLICY >
    void
    computeMassMatrix( localIndex const size,
                       arrayView2d< WaveSolverBase::wsCoordType const, nodes::REFERENCE_POSITION_USD > const nodeCoords,
                       arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodes,
                       arrayView1d< real32 const > const velocity,
                       arrayView1d< real32 const > const density,
                       arrayView1d< real32 > const mass )

    {
      forAll< EXEC_POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const e )
      {

        real32 const invC2 = 1.0 / ( density[e] * pow( velocity[e], 2 ) );
        // only the eight corners of the mesh cell are needed to compute the Jacobian
        real64 xLocal[ 8 ][ 3 ];
        for( localIndex a = 0; a < 8; ++a )
        {
          localIndex const nodeIndex = elemsToNodes( e, FE_TYPE::meshIndexToLinearIndex3D( a ) );
          for( localIndex i = 0; i < 3; ++i )
          {
            xLocal[a][i] = nodeCoords( nodeIndex, i );
          }
        }
        constexpr localIndex numQuadraturePointsPerElem = FE_TYPE::numQuadraturePoints;
        for( localIndex q = 0; q < numQuadraturePointsPerElem; ++q )
        {
          real32 const localIncrement = invC2 * m_finiteElement.computeMassTerm( q, xLocal );
          RAJA::atomicAdd< ATOMIC_POLICY >( &mass[elemsToNodes( e, q )], localIncrement );
        }
      } );    // end loop over element
    }

    FE_TYPE const & m_finiteElement;
  };
  template< typename FE_TYPE >
  struct DampingMatrix
  {

    DampingMatrix( FE_TYPE const & finiteElement )
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
     * @param[in] density cell-wise density
     * @param[out] damping diagonal of the damping matrix
     */
    template< typename EXEC_POLICY, typename ATOMIC_POLICY >
    void
    computeDampingMatrix( localIndex const size,
                          arrayView2d< WaveSolverBase::wsCoordType const, nodes::REFERENCE_POSITION_USD > const nodeCoords,
                          arrayView2d< localIndex const > const elemsToFaces,
                          ArrayOfArraysView< localIndex const > const facesToNodes,
                          arrayView1d< integer const > const facesDomainBoundaryIndicator,
                          arrayView1d< localIndex const > const freeSurfaceFaceIndicator,
                          arrayView1d< real32 const > const velocity,
                          arrayView1d< real32 const > const density,
                          arrayView1d< real32 > const damping )
    {
      forAll< EXEC_POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const e )
      {
        for( localIndex i = 0; i < elemsToFaces.size( 1 ); ++i )
        {
          localIndex const f = elemsToFaces( e, i );
          // face on the domain boundary and not on free surface
          if( facesDomainBoundaryIndicator[f] == 1 && freeSurfaceFaceIndicator[f] != 1 )
          {
            // only the four corners of the mesh face are needed to compute the Jacobian
            real64 xLocal[ 4 ][ 3 ];
            for( localIndex a = 0; a < 4; ++a )
            {
              localIndex const nodeIndex = facesToNodes( f, FE_TYPE::meshIndexToLinearIndex2D( a ) );
              for( localIndex d = 0; d < 3; ++d )
              {
                xLocal[a][d] = nodeCoords( nodeIndex, d );
              }
            }
            real32 const alpha = 1.0 / (density[e] * velocity[e]);
            constexpr localIndex numNodesPerFace = FE_TYPE::numNodesPerFace;
            for( localIndex q = 0; q < numNodesPerFace; ++q )
            {
              real32 const localIncrement = alpha * m_finiteElement.computeDampingTerm( q, xLocal );
              RAJA::atomicAdd< ATOMIC_POLICY >( &damping[facesToNodes( f, q )], localIncrement );
            }
          }
        }
      } );
    }

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
     * @param[in] lateralSurfaceFaceIndicator flag equal to 1 if the face is on the lateral surface, and to 0 otherwise
     * @param[in] bottomSurfaceFaceIndicator flag equal to 1 if the face is on the bottom surface, and to 0 otherwise
     * @param[in] velocity cell-wise velocity
     * @param[in] density cell-wise density
     * @param[in] vtiEpsilon cell-wise epsilon (Thomsen parameter)
     * @param[in] vtiDelta density cell-wise delta (Thomsen parameter)
     * @param[in] vti_sigma sigma cell-wise parameter
     * @param[out] damping_pp Damping matrix D^{pp}
     * @param[out] damping_pq Damping matrix D^{pq}
     * @param[out] damping_qp Damping matrix D^{qp}
     * @param[out] damping_qq Damping matrix D^{qq}
     */
    template< typename EXEC_POLICY, typename ATOMIC_POLICY >
    void
    computeVTIFletcherDampingMatrices( localIndex const size,
                                       arrayView2d< WaveSolverBase::wsCoordType const, nodes::REFERENCE_POSITION_USD > const nodeCoords,
                                       arrayView2d< localIndex const > const elemsToFaces,
                                       ArrayOfArraysView< localIndex const > const facesToNodes,
                                       arrayView1d< integer const > const facesDomainBoundaryIndicator,
                                       arrayView1d< localIndex const > const freeSurfaceFaceIndicator,
                                       arrayView1d< localIndex const > const lateralSurfaceFaceIndicator,
                                       arrayView1d< localIndex const > const bottomSurfaceFaceIndicator,
                                       arrayView1d< real32 const > const velocity,
                                       arrayView1d< real32 const > const density,
                                       arrayView1d< real32 const > const vtiEpsilon,
                                       arrayView1d< real32 const > const vtiDelta,
                                       arrayView1d< real32 const > const vti_sigma,
                                       arrayView1d< real32 > const damping_pp,
                                       arrayView1d< real32 > const damping_pq,
                                       arrayView1d< real32 > const damping_qp,
                                       arrayView1d< real32 > const damping_qq )
    {
      forAll< EXEC_POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const e )
      {
        for( localIndex i = 0; i < elemsToFaces.size( 1 ); ++i )
        {
          localIndex const f = elemsToFaces( e, i );
          // face on the domain boundary and not on free surface
          if( facesDomainBoundaryIndicator[f] == 1 && freeSurfaceFaceIndicator[f] != 1 )
          {
            // only the four corners of the mesh face are needed to compute the Jacobian
            real64 xLocal[ 4 ][ 3 ];
            for( localIndex a = 0; a < 4; ++a )
            {
              localIndex const nodeIndex = facesToNodes( f, FE_TYPE::meshIndexToLinearIndex2D( a ) );
              for( localIndex d = 0; d < 3; ++d )
              {
                xLocal[a][d] = nodeCoords( nodeIndex, d );
              }
            }
            constexpr localIndex numNodesPerFace = FE_TYPE::numNodesPerFace;
            real32 vti_f = 1 - (vtiEpsilon[e] - vtiDelta[e]) / vti_sigma[e];
            if( lateralSurfaceFaceIndicator[f] == 1 )
            {
              // ABC coefficients
              real32 alpha = 1.0 / (velocity[e] * density[e] * sqrt( 1+2*vtiEpsilon[e] ));
              // VTI coefficients
              real32 vti_p_xy = 0;
              real32 vti_q_xy = 0;
              real32 vti_qp_xy= 0;

              vti_p_xy  = (1+2*vtiEpsilon[e]);
              vti_q_xy  = -(vti_f - 1);
              vti_qp_xy = (vti_f+2*vtiDelta[e]);
              for( localIndex q = 0; q < numNodesPerFace; ++q )
              {
                real32 const aux = m_finiteElement.computeDampingTerm( q, xLocal );
                real32 const localIncrement_p = alpha * vti_p_xy  * aux;
                RAJA::atomicAdd< ATOMIC_POLICY >( &damping_pp[facesToNodes( f, q )], localIncrement_p );

                real32 const localIncrement_q = alpha*vti_q_xy * aux;
                RAJA::atomicAdd< ATOMIC_POLICY >( &damping_qq[facesToNodes( f, q )], localIncrement_q );

                real32 const localIncrement_qp = alpha*vti_qp_xy * aux;
                RAJA::atomicAdd< ATOMIC_POLICY >( &damping_qp[facesToNodes( f, q )], localIncrement_qp );
              }
            }
            if( bottomSurfaceFaceIndicator[f] == 1 )
            {
              // ABC coefficients updated to fit horizontal velocity
              real32 alpha = 1.0 / (velocity[e] * density[e]);
              // VTI coefficients
              real32 vti_p_z  = 0;
              real32 vti_pq_z = 0;
              real32 vti_q_z  = 0;
              vti_p_z  = -(vti_f - 1);
              vti_pq_z = vti_f;
              vti_q_z  = 1;
              for( localIndex q = 0; q < numNodesPerFace; ++q )
              {
                real32 const aux = m_finiteElement.computeDampingTerm( q, xLocal );
                real32 const localIncrement_p = alpha * vti_p_z * aux;
                RAJA::atomicAdd< ATOMIC_POLICY >( &damping_pp[facesToNodes( f, q )], localIncrement_p );

                real32 const localIncrement_pq = alpha*vti_pq_z * aux;
                RAJA::atomicAdd< ATOMIC_POLICY >( &damping_pq[facesToNodes( f, q )], localIncrement_pq );

                real32 const localIncrement_q = alpha * vti_q_z * aux;
                RAJA::atomicAdd< ATOMIC_POLICY >( &damping_qq[facesToNodes( f, q )], localIncrement_q );
              }
            }
          }
        }
      } );
    }

#if 0
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
     * @param[in] lateralSurfaceFaceIndicator flag equal to 1 if the face is on the lateral surface, and to 0 otherwise
     * @param[in] bottomSurfaceFaceIndicator flag equal to 1 if the face is on the bottom surface, and to 0 otherwise
     * @param[in] velocity cell-wise velocity
     * @param[in] density cell-wise density
     * @param[in] vtiEpsilon cell-wise epsilon (Thomsen parameter)
     * @param[in] vtiDelta density cell-wise delta (Thomsen parameter)
     * @param[out] damping_pp Damping matrix D^{pp}
     * @param[out] damping_pq Damping matrix D^{pq}
     * @param[out] damping_qp Damping matrix D^{qp}
     * @param[out] damping_qq Damping matrix D^{qq}
     */
    template< typename EXEC_POLICY, typename ATOMIC_POLICY >
    void
    computeVTIZhangDampingMatrices( localIndex const size,
                                    arrayView2d< WaveSolverBase::wsCoordType const, nodes::REFERENCE_POSITION_USD > const nodeCoords,
                                    arrayView2d< localIndex const > const elemsToFaces,
                                    ArrayOfArraysView< localIndex const > const facesToNodes,
                                    arrayView1d< integer const > const facesDomainBoundaryIndicator,
                                    arrayView1d< localIndex const > const freeSurfaceFaceIndicator,
                                    arrayView1d< localIndex const > const lateralSurfaceFaceIndicator,
                                    arrayView1d< localIndex const > const bottomSurfaceFaceIndicator,
                                    arrayView1d< real32 const > const velocity,
                                    arrayView1d< real32 const > const density,
                                    arrayView1d< real32 const > const dofEpsilon,
                                    arrayView1d< real32 const > const dofDelta,
                                    arrayView1d< real32 > const damping_pp,
                                    arrayView1d< real32 > const damping_pq,
                                    arrayView1d< real32 > const damping_qp,
                                    arrayView1d< real32 > const damping_qq )
    {
      forAll< EXEC_POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const e )
      {
        for( localIndex i = 0; i < elemsToFaces.size( 1 ); ++i )
        {
          localIndex const f = elemsToFaces( e, i );
          // face on the domain boundary and not on free surface
          if( facesDomainBoundaryIndicator[f] == 1 && freeSurfaceFaceIndicator[f] != 1 )
          {
            // only the four corners of the mesh face are needed to compute the Jacobian
            real64 xLocal[ 4 ][ 3 ];
            for( localIndex a = 0; a < 4; ++a )
            {
              localIndex const nodeIndex = facesToNodes( f, FE_TYPE::meshIndexToLinearIndex2D( a ) );
              for( localIndex d = 0; d < 3; ++d )
              {
                xLocal[a][d] = nodeCoords( nodeIndex, d );
              }
            }
            constexpr localIndex numNodesPerFace = FE_TYPE::numNodesPerFace;
            // debug

            for( localIndex q = 0; q < numNodesPerFace; ++q )
            {
              //            real32 epsi = std::fabs( vtiEpsilon[e] );
              real32 epsi = std::fabs( dofEpsilon[facesToNodes( f, q )] );
              // end debug

              if( std::fabs( epsi ) < 1e-5 )
                epsi = 0;
              // debug
              //         real32 delt = std::fabs( vtiDelta[e] );
              real32 delt = std::fabs( dofDelta[facesToNodes( f, q )] );
              // end debug
              if( std::fabs( delt ) < 1e-5 )
                delt = 0;
              if( delt > epsi )
                delt = epsi;
              real32 sqrtEpsi = sqrt( 1 + 2 * epsi );
              // debug
              //            real32 sqrtDelta = sqrt( 1 + 2 * vtiDelta[e] );
              real32 sqrtDelta = sqrt( 1 + 2 * dofDelta[facesToNodes( f, q )] );
              // end debug
              if( lateralSurfaceFaceIndicator[f] == 1 )
              {
                // ABC coefficients updated to fit horizontal velocity
                real32 alpha = 1.0 / (velocity[e] * density[e] * sqrtEpsi);
                // VTI coefficients
                real32 vti_p_xy  = 1 + 2 * epsi;
                real32 vti_qp_xy = sqrtDelta;

                for( localIndex qq = 0; qq < numNodesPerFace; ++qq )
                {
                  real32 const aux = m_finiteElement.computeDampingTerm( qq, xLocal );
                  real32 const localIncrement_p = alpha* vti_p_xy  * aux;
                  RAJA::atomicAdd< ATOMIC_POLICY >( &damping_pp[facesToNodes( f, qq )], localIncrement_p );

                  real32 const localIncrement_qp = alpha * vti_qp_xy * aux;
                  RAJA::atomicAdd< ATOMIC_POLICY >( &damping_qp[facesToNodes( f, qq )], localIncrement_qp );
                }
              }
              if( bottomSurfaceFaceIndicator[f] == 1 )
              {
                // ABC coefficients updated to fit horizontal velocity
                real32 alpha = 1.0 / (velocity[e] * density[e]);
                // VTI coefficients
                real32 vti_pq_z = sqrtDelta;
                real32 vti_q_z  = 1;
                for( localIndex qq = 0; qq < numNodesPerFace; ++qq )
                {
                  real32 const aux = m_finiteElement.computeDampingTerm( qq, xLocal );

                  real32 const localIncrement_pq = alpha * vti_pq_z * aux;
                  RAJA::atomicAdd< ATOMIC_POLICY >( &damping_pq[facesToNodes( f, qq )], localIncrement_pq );

                  real32 const localIncrement_q = alpha * vti_q_z * aux;
                  RAJA::atomicAdd< ATOMIC_POLICY >( &damping_qq[facesToNodes( f, qq )], localIncrement_q );
                }
              }
            }// Debug
          }
        }
      } );
    }

#else
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
     * @param[in] lateralSurfaceFaceIndicator flag equal to 1 if the face is on the lateral surface, and to 0 otherwise
     * @param[in] bottomSurfaceFaceIndicator flag equal to 1 if the face is on the bottom surface, and to 0 otherwise
     * @param[in] velocity cell-wise velocity
     * @param[in] density cell-wise density
     * @param[in] vtiEpsilon cell-wise epsilon (Thomsen parameter)
     * @param[in] vtiDelta density cell-wise delta (Thomsen parameter)
     * @param[out] damping_pp Damping matrix D^{pp}
     * @param[out] damping_pq Damping matrix D^{pq}
     * @param[out] damping_qp Damping matrix D^{qp}
     * @param[out] damping_qq Damping matrix D^{qq}
     */
    template< typename EXEC_POLICY, typename ATOMIC_POLICY >
    void
    computeVTIZhangDampingMatrices( localIndex const size,
                                    arrayView2d< WaveSolverBase::wsCoordType const, nodes::REFERENCE_POSITION_USD > const nodeCoords,
                                    arrayView2d< localIndex const > const elemsToFaces,
                                    ArrayOfArraysView< localIndex const > const facesToNodes,
                                    arrayView1d< integer const > const facesDomainBoundaryIndicator,
                                    arrayView1d< localIndex const > const freeSurfaceFaceIndicator,
                                    arrayView1d< localIndex const > const lateralSurfaceFaceIndicator,
                                    arrayView1d< localIndex const > const bottomSurfaceFaceIndicator,
                                    arrayView1d< real32 const > const velocity,
                                    arrayView1d< real32 const > const density,
                                    arrayView1d< real32 const > const vtiEpsilon,
                                    arrayView1d< real32 const > const vtiDelta,
                                    arrayView1d< real32 > const damping_pp,
                                    arrayView1d< real32 > const damping_pq,
                                    arrayView1d< real32 > const damping_qp,
                                    arrayView1d< real32 > const damping_qq )
    {
      forAll< EXEC_POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const e )
      {
        for( localIndex i = 0; i < elemsToFaces.size( 1 ); ++i )
        {
          localIndex const f = elemsToFaces( e, i );
          // face on the domain boundary and not on free surface
          if( facesDomainBoundaryIndicator[f] == 1 && freeSurfaceFaceIndicator[f] != 1 )
          {
            // only the four corners of the mesh face are needed to compute the Jacobian
            real64 xLocal[ 4 ][ 3 ];
            for( localIndex a = 0; a < 4; ++a )
            {
              localIndex const nodeIndex = facesToNodes( f, FE_TYPE::meshIndexToLinearIndex2D( a ) );
              for( localIndex d = 0; d < 3; ++d )
              {
                xLocal[a][d] = nodeCoords( nodeIndex, d );
              }
            }
            constexpr localIndex numNodesPerFace = FE_TYPE::numNodesPerFace;
            // debug


            {
              real32 epsi = std::fabs( vtiEpsilon[e] );

              // end debug

              if( std::fabs( epsi ) < 1e-5 )
                epsi = 0;
              // debug
              real32 delt = std::fabs( vtiDelta[e] );

              // end debug
              if( std::fabs( delt ) < 1e-5 )
                delt = 0;
              if( delt > epsi )
                delt = epsi;
              real32 sqrtEpsi = sqrt( 1 + 2 * epsi );
              // debug
              real32 sqrtDelta = sqrt( 1 + 2 * vtiDelta[e] );

              // end debug
              if( lateralSurfaceFaceIndicator[f] == 1 )
              {
                // ABC coefficients updated to fit horizontal velocity
                real32 alpha = 1.0 / (velocity[e] * density[e] * sqrtEpsi);
                // VTI coefficients
                real32 vti_p_xy  = 1 + 2 * epsi;
                real32 vti_qp_xy = sqrtDelta;

                for( localIndex qq = 0; qq < numNodesPerFace; ++qq )
                {
                  real32 const aux = m_finiteElement.computeDampingTerm( qq, xLocal );
                  real32 const localIncrement_p = alpha* vti_p_xy  * aux;
                  RAJA::atomicAdd< ATOMIC_POLICY >( &damping_pp[facesToNodes( f, qq )], localIncrement_p );

                  real32 const localIncrement_qp = alpha * vti_qp_xy * aux;
                  RAJA::atomicAdd< ATOMIC_POLICY >( &damping_qp[facesToNodes( f, qq )], localIncrement_qp );
                }
              }
              if( bottomSurfaceFaceIndicator[f] == 1 )
              {
                // ABC coefficients updated to fit horizontal velocity
                real32 alpha = 1.0 / (velocity[e] * density[e]);
                // VTI coefficients
                real32 vti_pq_z = sqrtDelta;
                real32 vti_q_z  = 1;
                for( localIndex qq = 0; qq < numNodesPerFace; ++qq )
                {
                  real32 const aux = m_finiteElement.computeDampingTerm( qq, xLocal );

                  real32 const localIncrement_pq = alpha * vti_pq_z * aux;
                  RAJA::atomicAdd< ATOMIC_POLICY >( &damping_pq[facesToNodes( f, qq )], localIncrement_pq );

                  real32 const localIncrement_q = alpha * vti_q_z * aux;
                  RAJA::atomicAdd< ATOMIC_POLICY >( &damping_qq[facesToNodes( f, qq )], localIncrement_q );
                }
              }
            }// Debug
          }
        }
      } );
    }
#endif

    /// The finite element space/discretization object for the element type in the subRegion
    FE_TYPE const & m_finiteElement;
  };

  template< typename FE_TYPE >
  struct GradientKappaBuoyancy
  {

    GradientKappaBuoyancy( FE_TYPE const & finiteElement )
      : m_finiteElement( finiteElement )
    {}

    /**
     * @brief Launch the computation of the 2 gradients relative to the coeff of the wave equation K=1/rho*c2 and b=1/rho
     * @tparam EXEC_POLICY the execution policy
     * @tparam ATOMIC_POLICY the atomic policy
     * @param[in] size the number of cells in the subRegion
     * @param[in] nodeCoords coordinates of the nodes
     * @param[in] elemsToNodes map from element to nodes
     * @param[in] q_dt2 second order derivative in time of backward
     * @param[in] q_n current time step of backward
     * @param[in] p_n current time step of forward
     * @param[out] grad first part of gradient vector with respect to K=1/rho*c2
     * @param[out] grad2 second part of gradient vector with respact to b=1/rho
     */
    template< typename EXEC_POLICY, typename ATOMIC_POLICY >
    void
    computeGradient( localIndex const size,
                     arrayView2d< WaveSolverBase::wsCoordType const, nodes::REFERENCE_POSITION_USD > const nodeCoords,
                     arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodes,
                     arrayView1d< integer const > const elemGhostRank,
                     arrayView1d< real32 const > const q_dt2,
                     arrayView1d< real32 const > const q_n,
                     arrayView1d< real32 const > const p_n,
                     arrayView1d< real32 > const grad,
                     arrayView1d< real32 > const grad2 )

    {
      forAll< EXEC_POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const e )
      {
        if( elemGhostRank[e]<0 )
        {
          // only the eight corners of the mesh cell are needed to compute the Jacobian
          real64 xLocal[ 8 ][ 3 ];
          for( localIndex a = 0; a < 8; ++a )
          {
            localIndex const nodeIndex = elemsToNodes( e, FE_TYPE::meshIndexToLinearIndex3D( a ) );
            for( localIndex i = 0; i < 3; ++i )
            {
              xLocal[a][i] = nodeCoords( nodeIndex, i );
            }
          }
          constexpr localIndex numQuadraturePointsPerElem = FE_TYPE::numQuadraturePoints;
          for( localIndex q = 0; q < numQuadraturePointsPerElem; ++q )
          {
            localIndex nodeIdx = elemsToNodes( e, q );
            grad[e] += q_dt2[nodeIdx] * p_n[nodeIdx] * m_finiteElement.computeMassTerm( q, xLocal );
            m_finiteElement.template computeStiffnessTerm<>( q, xLocal, [&] ( const int i, const int j, const real64 val )
            {
              grad2[e] += val* q_n[elemsToNodes( e, j )] * p_n[elemsToNodes( e, i )];
            } );
          }
        }
      } );    // end loop over element
    }
    /// The finite element space/discretization object for the element type in the subRegion
    FE_TYPE const & m_finiteElement;
  };

  template< typename FE_TYPE >
  struct ImagingCondition
  {

    ImagingCondition( FE_TYPE const & finiteElement )
      : m_finiteElement( finiteElement )
    {}

    /**
     * @brief Launch the computation of the imaging condition for RTM
     * @tparam EXEC_POLICY the execution policy
     * @tparam ATOMIC_POLICY the atomic policy
     * @param[in] size the number of cells in the subRegion
     * @param[in] nodeCoords coordinates of the nodes
     * @param[in] elemsToNodes map from element to nodes
     * @param[in] q_n current time step of backward
     * @param[in] p_n current time step of forward
     * @param[out] imag will contain the imaging condition value
     */
    template< typename EXEC_POLICY, typename ATOMIC_POLICY >
    void
    computeImagingCondition( localIndex const size,
                             arrayView2d< WaveSolverBase::wsCoordType const, nodes::REFERENCE_POSITION_USD > const nodeCoords,
                             arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodes,
                             arrayView1d< integer const > const elemGhostRank,
                             arrayView1d< real32 const > const q_n,
                             arrayView1d< real32 const > const p_n,
                             arrayView1d< real32 > const imag )
    {
      forAll< EXEC_POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const e )
      {
        if( elemGhostRank[e]<0 )
        {
          // only the eight corners of the mesh cell are needed to compute the Jacobian
          real64 xLocal[ 8 ][ 3 ];
          for( localIndex a = 0; a < 8; ++a )
          {
            localIndex const nodeIndex = elemsToNodes( e, FE_TYPE::meshIndexToLinearIndex3D( a ) );
            for( localIndex i = 0; i < 3; ++i )
            {
              xLocal[a][i] = nodeCoords( nodeIndex, i );
            }
          }
          constexpr localIndex numQuadraturePointsPerElem = FE_TYPE::numQuadraturePoints;
          for( localIndex q = 0; q < numQuadraturePointsPerElem; ++q )
          {
            localIndex nodeIdx = elemsToNodes( e, q );
            imag[e] += q_n[nodeIdx] * p_n[nodeIdx] * m_finiteElement.computeMassTerm( q, xLocal );
          }
        }
      } );    // end loop over element
    }
    /// The finite element space/discretization object for the element type in the subRegion
    FE_TYPE const & m_finiteElement;
  };
};

} // namespace geos

#endif // GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICMATRICESSEMKERNEL_HPP_
