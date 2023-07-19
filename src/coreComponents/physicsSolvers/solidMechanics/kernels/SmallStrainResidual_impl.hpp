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
                                                                                        string const GEOS_UNUSED_PARAM( elementListName ) ):
  Base( elementSubRegion,
        finiteElementSpace,
        inputConstitutiveType ),
  m_X( nodeManager.referencePosition()),
  m_input( inputSrc ),
  m_res( inputDst )
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
  #pragma unroll
  for( localIndex a=0; a< numNodesPerElem; ++a )
  {
    localIndex const nodeIndex = m_elemsToNodes( k, a );
    #pragma unroll
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
  #define USE_JACOBIAN
#if !defined( USE_JACOBIAN )
  real64 dNdX[ numNodesPerElem ][ 3 ];
  real64 const detJ = FE_TYPE::calcGradN( q, stack.xLocal, dNdX );

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
#endif // !defined( USE_JACOBIAN )
}
#undef USE_JACOBIAN

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
#if !defined(USE_JACOBIAN)
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
#elif USE_JACOBIAN==1
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

  #pragma unroll
  for( int i = 0; i < 3; ++i )
  {
    #pragma unroll
    for( int j = 0; j < 3; ++j )
    {
      #pragma unroll
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

//  FE_TYPE::template parentGradient< qa, qb, qc >( xLocal, invJ);
  real64 const detJ = LvArray::tensorOps::invert< 3 >( invJ );

//  FE_TYPE::template parentGradient< qa, qb, qc >( varLocal, parentGradVar);
  real64 gradVar[3][3] = {{0}};

  #pragma unroll
  for( int i = 0; i < 3; ++i )
  {
    #pragma unroll
    for( int j = 0; j < 3; ++j )
    {
      #pragma unroll
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
  GEOS_MARK_FUNCTION;

  #define KERNEL_OPTION 3
#if KERNEL_OPTION == 1
  forAll< POLICY >( numElems,
                    [=] GEOS_DEVICE ( localIndex const k )
  {
    typename KERNEL_TYPE::StackVariables stack;

    kernelComponent.setup( k, stack );
//    for( integer q=0; q<KERNEL_TYPE::numQuadraturePointsPerElem; ++q )
  #pragma unroll
    for( integer qa=0; qa<2; ++qa )
                                     #pragma unroll
      for( integer qb=0; qb<2; ++qb )
                                       #pragma unroll
        for( integer qc=0; qc<2; ++qc )
        {
          //  int qa, qb, qc;
          //FE_TYPE::LagrangeBasis1::TensorProduct3D::multiIndex( q, qa, qb, qc );

          int const q = qa + 2 * qb + 4*qc;
          kernelComponent.quadraturePointKernel( k, q, stack );
        }
    kernelComponent.complete( k, stack );
  } );
#elif KERNEL_OPTION == 2 // no stack

  forAll< POLICY >( numElems,
                    [=] GEOS_DEVICE ( localIndex const k )
  {
    real64 fLocal[ KERNEL_TYPE::numNodesPerElem ][ 3 ] = {{0}};
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
  forAll< POLICY >( numElems,
                    [=] GEOS_DEVICE ( localIndex const k )
  {
    real64 fLocal[ KERNEL_TYPE::numNodesPerElem ][ 3 ] = {{0}};
    real64 varLocal[ KERNEL_TYPE::numNodesPerElem ][ 3 ];
    real64 xLocal[ KERNEL_TYPE::numNodesPerElem ][ 3 ];

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

#endif
  return 0;
}


#undef UPDATE_STRESS

} // namespace solidMechanicsLagrangianFEMKernels

} // namespace geos

#endif //GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_KERNELS_SMALLSTRAINRESIDUAL_IMPL_HPP_
