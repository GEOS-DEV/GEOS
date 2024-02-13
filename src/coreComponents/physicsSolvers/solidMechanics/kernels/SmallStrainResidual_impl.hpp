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
 * @file SmallStrainResidual_impl.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_KERNELS_SMALLSTRAINRESIDUAL_IMPL_HPP_
#define GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_KERNELS_SMALLSTRAINRESIDUAL_IMPL_HPP_

#include "SmallStrainResidual.hpp"


namespace geos
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
                                                                                        real64 const GEOS_UNUSED_PARAM( dt ),
                                                                                        string const GEOS_UNUSED_PARAM( elementListName ),
                                                                                        int const kernelOptimizationOption ):
  Base( elementSubRegion,
        finiteElementSpace,
        inputConstitutiveType ),
  m_X( nodeManager.referencePosition()),
  m_input( inputSrc ),
  m_res( inputDst ),
  m_kernelOptimizationOption( kernelOptimizationOption )
{
  GEOS_UNUSED_VAR( edgeManager );
  GEOS_UNUSED_VAR( faceManager );
  GEOS_UNUSED_VAR( targetRegionIndex );
}


template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void SmallStrainResidual< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::setup( localIndex const k,
                                                                               StackVariables & stack ) const
{
  RAJA_UNROLL
  for( localIndex a=0; a< numNodesPerElem; ++a )
  {
    localIndex const nodeIndex = m_elemsToNodes( k, a );
    RAJA_UNROLL
    for( int i=0; i<numDofPerTrialSupportPoint; ++i )
    {
#if defined(CALC_FEM_SHAPE_IN_KERNEL)
      stack.xLocal[ a ][ i ] = m_X( nodeIndex, i );
#endif
      stack.varLocal[ a ][ i ] = m_input( nodeIndex, i );
    }
  }
}

template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void SmallStrainResidual< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::quadraturePointKernel( localIndex const k,
                                                                                               localIndex const q,
                                                                                               StackVariables & stack ) const
{
  real64 invJ[3][3] = {{0}};
  real64 const detJ = FE_TYPE::invJacobianTransformation( q, stack.xLocal, invJ );

  real64 strain[6] = {0};
  FE_TYPE::symmetricGradient( q, invJ, stack.varLocal, strain );

  real64 stressLocal[ 6 ] = {0};
  m_constitutiveUpdate.smallStrainNoStateUpdate_StressOnly( k, q, strain, stressLocal );

  for( localIndex c = 0; c < 6; ++c )
  {
    stressLocal[ c ] *= -detJ;
  }

  FE_TYPE::plusGradNajAij( q, invJ, stressLocal, stack.fLocal );
}

template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
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
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void SmallStrainResidual< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::setup( localIndex const k,
                                                                               real64 (& xLocal) [ numNodesPerElem ][ numDofPerTrialSupportPoint ],
                                                                               real64 (& varLocal) [ numNodesPerElem ][ numDofPerTrialSupportPoint ] ) const
{
  RAJA_UNROLL
  for( localIndex a=0; a< numNodesPerElem; ++a )
  {
    localIndex const nodeIndex = m_elemsToNodes( k, a );
    RAJA_UNROLL
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
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void SmallStrainResidual< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::setup( localIndex const k,
                                                                               localIndex ( & elemToNodeMap )[numNodesPerElem],
                                                                               real64 (& xLocal) [ numNodesPerElem ][ numDofPerTrialSupportPoint ],
                                                                               real64 (& varLocal) [ numNodesPerElem ][ numDofPerTrialSupportPoint ] ) const
{
  RAJA_UNROLL
  for( localIndex a=0; a< numNodesPerElem; ++a )
  {
    elemToNodeMap[a] = m_elemsToNodes( k, a );
    RAJA_UNROLL
    for( int i=0; i<numDofPerTrialSupportPoint; ++i )
    {
      xLocal[ a ][ i ] = m_X( elemToNodeMap[a], i );
      varLocal[ a ][ i ] = m_input( elemToNodeMap[a], i );
    }
  }
}


template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void SmallStrainResidual< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::quadraturePointKernel( localIndex const k,
                                                                                               localIndex const qa,
                                                                                               localIndex const qb,
                                                                                               localIndex const qc,
                                                                                               real64 const (&xLocal) [ numNodesPerElem ][ numDofPerTrialSupportPoint ],
                                                                                               real64 const (&varLocal) [ numNodesPerElem ][ numDofPerTrialSupportPoint ],
                                                                                               real64 (& fLocal) [ numNodesPerElem ][ numDofPerTrialSupportPoint ] ) const
{
  #define USE_JACOBIAN 2

#if USE_JACOBIAN==1
  real64 invJ[3][3] = {{0}};
  real64 const detJ = FE_TYPE::invJacobianTransformation( qa, qb, qc, xLocal, invJ );

  real64 strain[6] = {0};
  FE_TYPE::symmetricGradient( qa, qb, qc, invJ, varLocal, strain );

  real64 stressLocal[ 6 ] = {0};
  m_constitutiveUpdate.smallStrainNoStateUpdate_StressOnly( k, qa+2*qb+4*qc, strain, stressLocal );

  for( localIndex c = 0; c < 6; ++c )
  {
    stressLocal[ c ] *= -detJ;
  }

  FE_TYPE::plusGradNajAij( qa, qb, qc, invJ, stressLocal, fLocal );
#else

  real64 invJ[3][3] = {{0}};
  real64 parentGradVar[3][3] = {{0}};

  FE_TYPE::parentGradient2( qa, qb, qc, xLocal, varLocal, invJ, parentGradVar );
  real64 const detJ = LvArray::tensorOps::invert< 3 >( invJ );

  real64 gradVar[3][3] = {{0}};

  RAJA_UNROLL
  for( int i = 0; i < 3; ++i )
  {
    RAJA_UNROLL
    for( int j = 0; j < 3; ++j )
    {
      RAJA_UNROLL
      for( int kk = 0; kk < 3; ++kk )
      {
        gradVar[i][j] = gradVar[i][j] + parentGradVar[i][kk] * invJ[kk][j];
      }
    }
  }
  real64 strain[6] = {0};
  strain[0] = gradVar[0][0];
  strain[1] = gradVar[1][1];
  strain[2] = gradVar[2][2];
  strain[3] = gradVar[2][1] + gradVar[1][2];
  strain[4] = gradVar[2][0] + gradVar[0][2];
  strain[5] = gradVar[1][0] + gradVar[0][1];


  real64 stressLocal[ 6 ] = {0};
  m_constitutiveUpdate.smallStrainNoStateUpdate_StressOnly( k, qa+2*qb+4*qc, strain, stressLocal );

  for( localIndex c = 0; c < 6; ++c )
  {
    stressLocal[ c ] *= -detJ;
  }

  FE_TYPE::plusGradNajAij( qa, qb, qc, invJ, stressLocal, fLocal );
#endif
#undef USE_JACOBIAN

}



template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
template< int qa, int qb, int qc >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void SmallStrainResidual< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::quadraturePointKernel( localIndex const k,
                                                                                               real64 const (&xLocal) [ numNodesPerElem ][ numDofPerTrialSupportPoint ],
                                                                                               real64 const (&varLocal) [ numNodesPerElem ][ numDofPerTrialSupportPoint ],
                                                                                               real64 (& fLocal) [ numNodesPerElem ][ numDofPerTrialSupportPoint ] ) const
{
  real64 invJ[3][3] = {{0}};
  real64 parentGradVar[3][3] = {{0}};

  FE_TYPE::template parentGradient2< qa, qb, qc >( xLocal, varLocal, invJ, parentGradVar );
  real64 const detJ = LvArray::tensorOps::invert< 3 >( invJ );

//  FE_TYPE::template parentGradient< qa, qb, qc >( varLocal, parentGradVar);
  real64 gradVar[3][3] = {{0}};

  RAJA_UNROLL
  for( int i = 0; i < 3; ++i )
  {
    RAJA_UNROLL
    for( int j = 0; j < 3; ++j )
    {
      RAJA_UNROLL
      for( int kk = 0; kk < 3; ++kk )
      {
        gradVar[i][j] = gradVar[i][j] + parentGradVar[i][kk] * invJ[kk][j];
      }
    }
  }
  real64 strain[6] = {0};
  strain[0] = gradVar[0][0];
  strain[1] = gradVar[1][1];
  strain[2] = gradVar[2][2];
  strain[3] = gradVar[2][1] + gradVar[1][2];
  strain[4] = gradVar[2][0] + gradVar[0][2];
  strain[5] = gradVar[1][0] + gradVar[0][1];


  real64 stressLocal[ 6 ] = {0};
  m_constitutiveUpdate.smallStrainNoStateUpdate_StressOnly( k, qa+2*qb+4*qc, strain, stressLocal );

  RAJA_UNROLL
  for( localIndex c = 0; c < 6; ++c )
  {
    stressLocal[ c ] *= -detJ;
  }

  FE_TYPE::template plusGradNajAij< qa, qb, qc >( invJ, stressLocal, fLocal );

}


template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 SmallStrainResidual< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::complete( localIndex const k,
                                                                                    real64 const (&fLocal) [ numNodesPerElem ][ numDofPerTestSupportPoint ] ) const
{
  RAJA_UNROLL
  for( localIndex a = 0; a < numNodesPerElem; ++a )
  {
    localIndex const nodeIndex = m_elemsToNodes( k, a );
    RAJA_UNROLL
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
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 SmallStrainResidual< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::complete( localIndex const ( &elemToNodeMap )[numNodesPerElem],
                                                                                    real64 const (&fLocal) [ numNodesPerElem ][ numDofPerTestSupportPoint ] ) const
{
  RAJA_UNROLL
  for( localIndex a = 0; a < numNodesPerElem; ++a )
  {
    RAJA_UNROLL
    for( int i = 0; i < numDofPerTestSupportPoint; ++i )
    {
      RAJA::atomicAdd< parallelDeviceAtomic >( &m_res( elemToNodeMap[a], i ), fLocal[ a ][ i ] );
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
  GEOS_MARK_FUNCTION;
  int const kernelOptimizationOption = kernelComponent.m_kernelOptimizationOption;

  if( kernelOptimizationOption == 1 )
  {
    forAll< POLICY >( numElems,
                      [=] GEOS_DEVICE ( localIndex const k )
    {
      typename KERNEL_TYPE::StackVariables stack;

      kernelComponent.setup( k, stack );
      //    for( integer q=0; q<KERNEL_TYPE::numQuadraturePointsPerElem; ++q )
      RAJA_UNROLL
      for( integer qa=0; qa<2; ++qa )
        RAJA_UNROLL
        for( integer qb=0; qb<2; ++qb )
          RAJA_UNROLL
          for( integer qc=0; qc<2; ++qc )
          {
            //  int qa, qb, qc;
            //FE_TYPE::LagrangeBasis1::TensorProduct3D::multiIndex( q, qa, qb, qc );

            int const q = qa + 2 * qb + 4*qc;
            kernelComponent.quadraturePointKernel( k, q, stack );
          }
      kernelComponent.complete( k, stack );
    } );
  }
  else if( kernelOptimizationOption == 2 )
  {
    forAll< POLICY >( numElems,
                      [=] GEOS_DEVICE ( localIndex const k )
    {
      real64 fLocal[ KERNEL_TYPE::numNodesPerElem ][ 3 ] = {{0}};
      real64 varLocal[ KERNEL_TYPE::numNodesPerElem ][ 3 ];
      real64 xLocal[ KERNEL_TYPE::numNodesPerElem ][ 3 ];

      kernelComponent.setup( k, xLocal, varLocal );
      for( integer qc=0; qc<2; ++qc )
      {
        for( integer qb=0; qb<2; ++qb )
        {
          for( integer qa=0; qa<2; ++qa )
          {
            kernelComponent.quadraturePointKernel( k, qa, qb, qc, xLocal, varLocal, fLocal );
          }
        }
      }
      kernelComponent.complete( k, fLocal );

    } );
  }
  else if( kernelOptimizationOption == 21 )
  {
    forAll< POLICY >( numElems,
                      [=] GEOS_DEVICE ( localIndex const k )
    {
      real64 fLocal[ KERNEL_TYPE::numNodesPerElem ][ 3 ] = {{0}};
      real64 varLocal[ KERNEL_TYPE::numNodesPerElem ][ 3 ];
      real64 xLocal[ KERNEL_TYPE::numNodesPerElem ][ 3 ];
      localIndex elementToNodeMap[ KERNEL_TYPE::numNodesPerElem ];

      kernelComponent.setup( k, elementToNodeMap, xLocal, varLocal );
      for( integer qc=0; qc<2; ++qc )
      {
        for( integer qb=0; qb<2; ++qb )
        {
          for( integer qa=0; qa<2; ++qa )
          {
            kernelComponent.quadraturePointKernel( k, qa, qb, qc, xLocal, varLocal, fLocal );
          }
        }
      }
      kernelComponent.complete( elementToNodeMap, fLocal );

    } );
  }
  else if( kernelOptimizationOption == 3 )
  {
    forAll< POLICY >( numElems,
                      [=] GEOS_DEVICE ( localIndex const k )
    {
      real64 fLocal[ KERNEL_TYPE::numNodesPerElem ][ 3 ] = {{0}};
      real64 varLocal[ KERNEL_TYPE::numNodesPerElem ][ 3 ];
      real64 xLocal[ KERNEL_TYPE::numNodesPerElem ][ 3 ];
//      localIndex elementToNodeMap[ KERNEL_TYPE::numNodesPerElem ] = {0};

      kernelComponent.setup( k, xLocal, varLocal );
      kernelComponent.template quadraturePointKernel< 0, 0, 0 >( k, xLocal, varLocal, fLocal );
      kernelComponent.template quadraturePointKernel< 0, 0, 1 >( k, xLocal, varLocal, fLocal );
      kernelComponent.template quadraturePointKernel< 0, 1, 0 >( k, xLocal, varLocal, fLocal );
      kernelComponent.template quadraturePointKernel< 0, 1, 1 >( k, xLocal, varLocal, fLocal );
      kernelComponent.template quadraturePointKernel< 1, 0, 0 >( k, xLocal, varLocal, fLocal );
      kernelComponent.template quadraturePointKernel< 1, 0, 1 >( k, xLocal, varLocal, fLocal );
      kernelComponent.template quadraturePointKernel< 1, 1, 0 >( k, xLocal, varLocal, fLocal );
      kernelComponent.template quadraturePointKernel< 1, 1, 1 >( k, xLocal, varLocal, fLocal );
      kernelComponent.complete( k, fLocal );

    } );
  }
  return 0;
}


#undef UPDATE_STRESS

} // namespace solidMechanicsLagrangianFEMKernels

} // namespace geos

#endif //GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_KERNELS_SMALLSTRAINRESIDUAL_IMPL_HPP_