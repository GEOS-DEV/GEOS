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
 * @file ImplicitSmallStrainQuasiStatic_impl.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_KERNELS_IMPLCITSMALLSTRAINQUASISTATIC_IMPL_HPP_
#define GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_KERNELS_IMPLCITSMALLSTRAINQUASISTATIC_IMPL_HPP_

#include "ImplicitSmallStrainQuasiStatic.hpp"

namespace geosx
{

namespace solidMechanicsLagrangianFEMKernels
{


template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
ImplicitSmallStrainQuasiStatic< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::
ImplicitSmallStrainQuasiStatic( NodeManager const & nodeManager,
                                EdgeManager const & edgeManager,
                                FaceManager const & faceManager,
                                localIndex const targetRegionIndex,
                                SUBREGION_TYPE const & elementSubRegion,
                                FE_TYPE const & finiteElementSpace,
                                CONSTITUTIVE_TYPE & inputConstitutiveType,
                                arrayView1d< globalIndex const > const inputDofNumber,
                                globalIndex const rankOffset,
                                CRSMatrixView< real64, globalIndex const > const inputMatrix,
                                arrayView1d< real64 > const inputRhs,
                                real64 const (&inputGravityVector)[3] ):
  Base( nodeManager,
        edgeManager,
        faceManager,
        targetRegionIndex,
        elementSubRegion,
        finiteElementSpace,
        inputConstitutiveType,
        inputDofNumber,
        rankOffset,
        inputMatrix,
        inputRhs ),
  m_X( nodeManager.referencePosition()),
  m_disp( nodeManager.getExtrinsicData< extrinsicMeshData::solidMechanics::totalDisplacement >() ),
  m_uhat( nodeManager.getExtrinsicData< extrinsicMeshData::solidMechanics::incrementalDisplacement >() ),
  m_gravityVector{ inputGravityVector[0], inputGravityVector[1], inputGravityVector[2] },
  m_density( inputConstitutiveType.getDensity() )
{}


template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void ImplicitSmallStrainQuasiStatic< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::
setup( localIndex const k,
       StackVariables & stack ) const
{
  for( localIndex a=0; a<numNodesPerElem; ++a )
  {
    localIndex const localNodeIndex = m_elemsToNodes( k, a );

    for( int i=0; i<3; ++i )
    {
#if defined(CALC_FEM_SHAPE_IN_KERNEL)
      stack.xLocal[ a ][ i ] = m_X[ localNodeIndex ][ i ];
#endif
      stack.u_local[ a ][i] = m_disp[ localNodeIndex ][i];
      stack.uhat_local[ a ][i] = m_uhat[ localNodeIndex ][i];
      stack.localRowDofIndex[a*3+i] = m_dofNumber[localNodeIndex]+i;
      stack.localColDofIndex[a*3+i] = m_dofNumber[localNodeIndex]+i;
    }
  }
}

template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
template< typename STRESS_MODIFIER >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void ImplicitSmallStrainQuasiStatic< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::quadraturePointKernel( localIndex const k,
                                                                                                          localIndex const q,
                                                                                                          StackVariables & stack,
                                                                                                          STRESS_MODIFIER && stressModifier ) const
{
  real64 dNdX[ numNodesPerElem ][ 3 ];
  real64 const detJ = m_finiteElementSpace.template getGradN< FE_TYPE >( k, q, stack.xLocal, dNdX );

  real64 strainInc[6] = {0};
  real64 stress[6] = {0};

  typename CONSTITUTIVE_TYPE::KernelWrapper::DiscretizationOps stiffness;

  FE_TYPE::symmetricGradient( dNdX, stack.uhat_local, strainInc );

  m_constitutiveUpdate.smallStrainUpdate( k, q, strainInc, stress, stiffness );

  stressModifier( stress );
  for( localIndex i=0; i<6; ++i )
  {
    stress[i] *= -detJ;
  }

  real64 const gravityForce[3] = { m_gravityVector[0] * m_density( k, q )* detJ,
                                   m_gravityVector[1] * m_density( k, q )* detJ,
                                   m_gravityVector[2] * m_density( k, q )* detJ };

  real64 N[numNodesPerElem];
  FE_TYPE::calcN( q, N );
  FE_TYPE::plusGradNajAijPlusNaFi( dNdX,
                                   stress,
                                   N,
                                   gravityForce,
                                   reinterpret_cast< real64 (&)[numNodesPerElem][3] >(stack.localResidual) );
  stiffness.template upperBTDB< numNodesPerElem >( dNdX, -detJ, stack.localJacobian );
}


template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64 ImplicitSmallStrainQuasiStatic< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::complete( localIndex const k,
                                                                                               StackVariables & stack ) const
{
  GEOSX_UNUSED_VAR( k );
  real64 maxForce = 0;

  // TODO: Does this work if BTDB is non-symmetric?
  CONSTITUTIVE_TYPE::KernelWrapper::DiscretizationOps::template fillLowerBTDB< numNodesPerElem >( stack.localJacobian );

  for( int localNode = 0; localNode < numNodesPerElem; ++localNode )
  {
    for( int dim = 0; dim < numDofPerTestSupportPoint; ++dim )
    {
      localIndex const dof =
        LvArray::integerConversion< localIndex >( stack.localRowDofIndex[ numDofPerTestSupportPoint * localNode + dim ] - m_dofRankOffset );
      if( dof < 0 || dof >= m_matrix.numRows() )
        continue;
      m_matrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( dof,
                                                                              stack.localRowDofIndex,
                                                                              stack.localJacobian[ numDofPerTestSupportPoint * localNode + dim ],
                                                                              numNodesPerElem * numDofPerTrialSupportPoint );

      RAJA::atomicAdd< parallelDeviceAtomic >( &m_rhs[ dof ], stack.localResidual[ numDofPerTestSupportPoint * localNode + dim ] );
      maxForce = fmax( maxForce, fabs( stack.localResidual[ numDofPerTestSupportPoint * localNode + dim ] ) );
    }
  }


  return maxForce;
}

template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
template< typename POLICY,
          typename KERNEL_TYPE >
real64
ImplicitSmallStrainQuasiStatic< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::kernelLaunch( localIndex const numElems,
                                                                                            KERNEL_TYPE const & kernelComponent )
{
  return Base::template kernelLaunch< POLICY, KERNEL_TYPE >( numElems, kernelComponent );
}


} // namespace solidMechanicsLagrangianFEMKernels

} // namespace geosx

#endif // GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_KERNELS_IMPLCITSMALLSTRAINQUASISTATIC_IMPL_HPP_
