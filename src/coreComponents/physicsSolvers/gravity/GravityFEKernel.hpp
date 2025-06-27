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
 * @file GravityFEKernel.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_GRAVITY_GRAVITYFEKERNEL_HPP_
#define GEOS_PHYSICSSOLVERS_GRAVITY_GRAVITYFEKERNEL_HPP_

//#include "finiteElement/kernelInterface/ImplicitKernelBase.hpp"

namespace geos
{

/// Namespace to contain the gravity kernels.
namespace gravityFEKernel
{

template< typename FE_TYPE >
struct DensityVolumeIntegralKernel
{
  explicit DensityVolumeIntegralKernel( FE_TYPE const & finiteElement )
    : m_finiteElement( finiteElement )
  {}

  /**
   * @brief Launches the precomputation of the DensityVolumeIntegral matrix.
   * @tparam EXEC_POLICY the execution policy
   * @tparam ATOMIC_POLICY the atomic policy
   * @param[in] size Number of elements in the subregion
   * @param[in] X Node coordinates
   * @param[in] elemsToNodes Element-to-node connectivity
   * @param[in] density Cell-wise density
   * @param[out] volumeIntegral Output mass matrix diagonal
   */
  template< typename EXEC_POLICY, typename ATOMIC_POLICY >
  void launch( localIndex const size,
               arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X,
               arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodes,
               arrayView1d< real64 const > const density,
               arrayView1d< real64 > const volumeIntegral ) const
  {
    forAll< EXEC_POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const k )
    {
      constexpr localIndex numNodesPerElem = FE_TYPE::numNodes;
      constexpr localIndex numQuadraturePoints = FE_TYPE::numQuadraturePoints;

      auto const rho = density[k];

      // Gather nodal coordinates for this element
      real64 xLocal[numNodesPerElem][3];
      for( localIndex a = 0; a < numNodesPerElem; ++a )
      {
        localIndex const nodeIdx = elemsToNodes( k, a );
        for( localIndex i = 0; i < 3; ++i )
        {
          xLocal[a][i] = X( nodeIdx, i );
        }
      }

      real64 N[numNodesPerElem];
      real64 gradN[numNodesPerElem][3];

      for( localIndex q = 0; q < numQuadraturePoints; ++q )
      {
        FE_TYPE::calcN( q, N );
        real64 const detJ = m_finiteElement.template getGradN< FE_TYPE >( k, q, xLocal, gradN );
        real64 const weightedJacobian = rho * detJ;

        for( localIndex a = 0; a < numNodesPerElem; ++a )
        {
          localIndex const nodeIdx = elemsToNodes( k, a );
          real64 const localIncrement = weightedJacobian * N[a];
          RAJA::atomicAdd< ATOMIC_POLICY >( &volumeIntegral[nodeIdx], localIncrement );
        }
      }
    } );
  }

private:
  FE_TYPE const & m_finiteElement;
};



template< typename FE_TYPE >
struct VolumeIntegralKernel_uni2
{
  explicit VolumeIntegralKernel_uni2( FE_TYPE const & finiteElement )
    : m_finiteElement( finiteElement )
  {}

  /**
   * @brief Launches the precomputation of the volume integral matrix.
   * @tparam EXEC_POLICY the execution policy
   * @tparam ATOMIC_POLICY the atomic policy
   * @param[in] size Number of elements in the subregion
   * @param[in] X Node coordinates
   * @param[in] elemsToNodes Element-to-node connectivity
   * @param[out] volumeIntegral Output volume integral matrix
   */
  template< typename EXEC_POLICY, typename ATOMIC_POLICY >
  void launch( localIndex const size,
               arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X,
               arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodes,
               arrayView2d< real64 > const volumeIntegral ) const
  {
    forAll< EXEC_POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const k )
    {
      constexpr localIndex numNodes = FE_TYPE::numNodes;
      constexpr localIndex numQuadPoints = FE_TYPE::numQuadraturePoints;

      // Gather nodal coordinates for this element
      real64 xLocal[numNodes][3];
      for( localIndex a = 0; a < numNodes; ++a )
      {
        localIndex const nodeIdx = elemsToNodes( k, a );
        for( localIndex i = 0; i < 3; ++i )
        {
          xLocal[a][i] = X( nodeIdx, i );
        }
      }

      real64 N[numNodes];
      real64 gradN[numNodes][3];

      for( localIndex q = 0; q < numQuadPoints; ++q )
      {
        FE_TYPE::calcN( q, N );
        real64 const detJ = m_finiteElement.template getGradN< FE_TYPE >( k, q, xLocal, gradN );

        for( localIndex a = 0; a < numNodes; ++a )
        {
          real64 const localIncrement = detJ * N[a];
          RAJA::atomicAdd< ATOMIC_POLICY >( &volumeIntegral( k, a ), localIncrement );
        }
      }
    } );
  }

private:
  FE_TYPE const & m_finiteElement;
};



} // namespace gravityFEKernel

} // namespace geos

#endif // GEOS_PHYSICSSOLVERS_GRAVITY_GRAVITYFEKERNEL_HPP_
