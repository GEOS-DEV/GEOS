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
 * @file ExplictSmallStrain_impl.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_KERNELS_SMALLSTRAINRESIDUAL_IMPL_HPP_
#define GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_KERNELS_SMALLSTRAINRESIDUAL_IMPL_HPP_

#include "SmallStrainResidual.hpp"


namespace geosx
{

/// Namespace to contain the solid mechanics kernels.
namespace solidMechanicsLagrangianFEMKernels
{


template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
SmallStrainResidual< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::SmallStrainResidual( NodeManager & nodeManager,
                                                                                        EdgeManager const & edgeManager,
                                                                                        FaceManager const & faceManager,
                                                                                        localIndex const targetRegionIndex,
                                                                                        SUBREGION_TYPE const & elementSubRegion,
                                                                                        FE_TYPE const & finiteElementSpace,
                                                                                        CONSTITUTIVE_TYPE & inputConstitutiveType,
                                                                                        arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const inputSrc,
                                                                                        arrayView2d< real64, nodes::TOTAL_DISPLACEMENT_USD > const inputDst,
                                                                                        real64 const dt,
                                                                                        string const elementListName ):
  Base( elementSubRegion,
        finiteElementSpace,
        inputConstitutiveType ),
  m_X( nodeManager.referencePosition()),
  m_input( inputSrc ),
  m_res( inputDst )
{
  GEOSX_UNUSED_VAR( edgeManager );
  GEOSX_UNUSED_VAR( faceManager );
  GEOSX_UNUSED_VAR( targetRegionIndex );
}


template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void SmallStrainResidual< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::setup( localIndex const k,
                                                                               StackVariables & stack ) const
{
  #pragma unroll
  for( localIndex a=0; a< numNodesPerElem; ++a )
  {
    localIndex const nodeIndex = m_elemsToNodes( k, a );
    #pragma unroll
    for( int i=0; i<numDofPerTrialSupportPoint; ++i )
    {
      stack.xLocal[ a ][ i ] = m_X( nodeIndex, i );
      stack.varLocal[ a ][ i ] = m_input( nodeIndex, i );
    }
  }
}

template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void SmallStrainResidual< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::quadraturePointKernel( localIndex const k,
                                                                                               localIndex const q,
                                                                                               StackVariables & stack ) const
{
//#define USE_JACOBIAN
#if !defined( USE_JACOBIAN )
  real64 dNdX[ numNodesPerElem ][ 3 ];
  real64 const detJ = FE_TYPE::calcGradN( q, stack.xLocal, dNdX );

  /// Macro to substitute in the shape function derivatives.
  real64 strain[6] = {0};
  FE_TYPE::symmetricGradient( dNdX, stack.varLocal, strain );

  real64 stressLocal[ 6 ] = {0};
  m_constitutiveUpdate.smallStrainNoStateUpdate_StressOnly( k, q, strain, stressLocal );
  #pragma unroll
  for( localIndex c = 0; c < 6; ++c )
  {
    stressLocal[ c ] *= -detJ;
  }
  FE_TYPE::plusGradNajAij( dNdX, stressLocal, stack.fLocal );
  
#else //defined( USE_JACOBIAN )
  real64 invJ[3][3];
  real64 const detJ = FE_TYPE::invJacobianTransformation( q, stack.xLocal, invJ );

  real64 strain[6] = {0};
  FE_TYPE::symmetricGradient( q, invJ, stack.varLocal, strain );

  real64 stressLocal[ 6 ] = {0};
  m_constitutiveUpdate.smallStrainNoStateUpdate_StressOnly( k, q, strain, stressLocal );
  for( localIndex c = 0; c < 6; ++c )
  {
    stressLocal[ c ] *= detJ;
  }
  FE_TYPE::plusGradNajAij( q, invJ, stressLocal, stack.fLocal );
#endif // !defined( USE_JACOBIAN )
}

template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64 SmallStrainResidual< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::complete( localIndex const k,
                                                                                    StackVariables const & stack ) const
{
  for( localIndex a = 0; a < numNodesPerElem; ++a )
  {
    localIndex const nodeIndex = m_elemsToNodes( k, a );
    for( int i = 0; i < numDofPerTestSupportPoint; ++i )
    {
      RAJA::atomicAdd< parallelDeviceAtomic >( &m_res( nodeIndex, i ), stack.fLocal[ a ][ i ] );
    }
  }
  return 0;
}
















template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void SmallStrainResidual< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::setup( localIndex const k,
                                                                               real64 (&xLocal) [ numNodesPerElem ][ numDofPerTrialSupportPoint ],
                                                                               real64 (&varLocal) [ numNodesPerElem ][ numDofPerTrialSupportPoint ] ) const
{
  #pragma unroll
  for( localIndex a=0; a< numNodesPerElem; ++a )
  {
    localIndex const nodeIndex = m_elemsToNodes( k, a );
    #pragma unroll
    for( int i=0; i<numDofPerTrialSupportPoint; ++i )
    {
      xLocal[ a ][ i ] = m_X( nodeIndex, i );
      varLocal[ a ][ i ] = m_input( nodeIndex, i );
    }
  }
}

template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void SmallStrainResidual< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::quadraturePointKernel( localIndex const k,
                                                                                               localIndex const qa,
                                                                                               localIndex const qb,
                                                                                               localIndex const qc,
                                                                                               real64 const (&xLocal) [ numNodesPerElem ][ numDofPerTrialSupportPoint ],
                                                                                               real64 const (&varLocal) [ numNodesPerElem ][ numDofPerTrialSupportPoint ],
                                                                                               real64 (&fLocal) [ numNodesPerElem ][ numDofPerTrialSupportPoint ] ) const
{
  real64 dNdX[ numNodesPerElem ][ 3 ];
  real64 const detJ = FE_TYPE::calcGradN( qa, qb, qc, xLocal, dNdX );

  /// Macro to substitute in the shape function derivatives.
  real64 strain[6] = {0};
  FE_TYPE::symmetricGradient( dNdX, varLocal, strain );

  real64 stressLocal[ 6 ] = {0};
  m_constitutiveUpdate.smallStrainNoStateUpdate_StressOnly( k, qa+2*qb+4*qc, strain, stressLocal );
  #pragma unroll
  for( localIndex c = 0; c < 6; ++c )
  {
    stressLocal[ c ] *= -detJ;
  }
  FE_TYPE::plusGradNajAij( dNdX, stressLocal, fLocal );
}

template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64 SmallStrainResidual< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::complete( localIndex const k,
                                                                                    real64 const (&fLocal) [ numNodesPerElem ][ numDofPerTrialSupportPoint ] ) const
{
  for( localIndex a = 0; a < numNodesPerElem; ++a )
  {
    localIndex const nodeIndex = m_elemsToNodes( k, a );
    for( int i = 0; i < numDofPerTestSupportPoint; ++i )
    {
      RAJA::atomicAdd< parallelDeviceAtomic >( &m_res( nodeIndex, i ), fLocal[ a ][ i ] );
    }
  }
  return 0;
}









template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
template< typename POLICY,
          typename KERNEL_TYPE >
real64
SmallStrainResidual< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::
kernelLaunch( localIndex const numElems,
              KERNEL_TYPE const & kernelComponent )
{
  GEOSX_MARK_FUNCTION;


#if 1
  forAll< POLICY >( numElems,
                    [=] GEOSX_DEVICE ( localIndex const k )
  {

    typename KERNEL_TYPE::StackVariables stack;

    kernelComponent.setup( k, stack );
//    for( integer q=0; q<KERNEL_TYPE::numQuadraturePointsPerElem; ++q )
      for( integer qc=0; qc<2; ++qc )
      for( integer qb=0; qb<2; ++qb )
      for( integer qa=0; qa<2; ++qa )
    {
      int const q = qa + 2 * qb + 4*qc;
      kernelComponent.quadraturePointKernel( k, q, stack );
    }
    kernelComponent.complete( k, stack );
  } );
#elif 1    
    forAll< POLICY >( numElems,
                      [=] GEOSX_DEVICE ( localIndex const k )
    {
      real64 fLocal[ KERNEL_TYPE::numNodesPerElem ][ 3 ];
      real64 varLocal[ KERNEL_TYPE::numNodesPerElem ][ 3 ];
      real64 xLocal[ KERNEL_TYPE::numNodesPerElem ][ 3 ];

      kernelComponent.setup( k, xLocal, varLocal );
      for( integer qc=0; qc<2; ++qc )
      for( integer qb=0; qb<2; ++qb )
      for( integer qa=0; qa<2; ++qa )
      {
        kernelComponent.quadraturePointKernel( k, qa, qb, qc, xLocal, varLocal, fLocal );
      }
      kernelComponent.complete( k, fLocal );

    } );
#else
   traits::ViewTypeConst< typename SUBREGION_TYPE::NodeMapType::base_type > const elemsToNodes = kernelComponent.m_elemsToNodes;
   arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X = kernelComponent.m_X;
   arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const src = kernelComponent.m_input;
   arrayView2d< real64, nodes::TOTAL_DISPLACEMENT_USD > const dst = kernelComponent.m_res;
   typename CONSTITUTIVE_TYPE::KernelWrapper const & constitutiveUpdate = kernelComponent.m_constitutiveUpdate;

  constexpr int numNodesPerElem = KERNEL_TYPE::numNodesPerElem;
  constexpr int dims = 3;
    forAll< POLICY >( numElems,
                      [=] GEOSX_DEVICE ( localIndex const k )
    {

      real64 xLocal[numNodesPerElem][dims];
      real64 srcLocal[numNodesPerElem][dims];
      real64 dstLocal[numNodesPerElem][dims] = {};
      localIndex nodeIndices[numNodesPerElem];

      for( localIndex a=0; a<numNodesPerElem; ++a )
      {
        nodeIndices[a] = elemsToNodes( k, a );
        for( int i=0; i<dims; ++i )
        {
          xLocal[ a ][ i ] = X[ nodeIndices[a] ][i];
          srcLocal[ a ][ i ] = src[ nodeIndices[a] ][i];
        }
      }


      for( integer q=0; q<KERNEL_TYPE::numQuadraturePointsPerElem; ++q )
      {
        real64 dNdX[ numNodesPerElem ][ dims ];
        real64 const detJ = FE_TYPE::calcGradN( q, xLocal, dNdX );
        /// Macro to substitute in the shape function derivatives.
        real64 strain[6] = {0};
        FE_TYPE::symmetricGradient( dNdX, srcLocal, strain );

        real64 stressLocal[ 6 ] = {0};
        constitutiveUpdate.smallStrainNoStateUpdate_StressOnly( k, q, strain, stressLocal );

        for( localIndex c = 0; c < 6; ++c )
        {
          stressLocal[ c ] *= -detJ;
        }
        FE_TYPE::plusGradNajAij( dNdX, stressLocal, dstLocal );
      }

      
      for( localIndex a = 0; a < numNodesPerElem; ++a )
      {
        for( int b = 0; b < numDofPerTestSupportPoint; ++b )
        {
          RAJA::atomicAdd< parallelDeviceAtomic >( &(dst[ nodeIndices[a] ][b]), dstLocal[ a ][ b ] );
        }
      }

    } );


#endif
  return 0;
}


#undef UPDATE_STRESS

} // namespace solidMechanicsLagrangianFEMKernels

} // namespace geosx

#endif //GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_KERNELS_SMALLSTRAINRESIDUAL_IMPL_HPP_
