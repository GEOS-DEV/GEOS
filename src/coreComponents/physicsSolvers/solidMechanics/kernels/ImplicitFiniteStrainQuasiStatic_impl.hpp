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
 * @file ImplicitFinitStrainQuasiStatic_impl.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_SOLIDMECHANICS_KERNELS_IMPLCITFINITSTRAINQUASISTATIC_IMPL_HPP_
#define GEOS_PHYSICSSOLVERS_SOLIDMECHANICS_KERNELS_IMPLCITFINITSTRAINQUASISTATIC_IMPL_HPP_

#include "ImplicitFinitStrainQuasiStatic.hpp"
namespace geos
{

namespace solidMechanicsLagrangianFEMKernels
{

template<typename SUBREGION_TYPE, typename CONSTITUTIVE_TYPE, typename FE_TYPE>
ImplicitFiniteStrainQuasis<SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE> :: ImplicitFinitStrainQuasiStatic(
    FaceManager const & faceManager, localIndex const targetRegionIndex,
    SUBREGION_TYPE const & elementSubRegion, FE_TYPE const & finiteElementSpace,
    CONSTITUTIVE_TYPE const & inputConstitutiveType, arrayView1d<globalIndex const> const inputDofNumber,
    globalIndex const rankOffset, CRSMatrixView<real64, globalIndex const> const inputMatrix,
    arrayView1d<real64> const inputRhs, real64 const inputDt, real64 const (&inputGravityVector)[3])
  : Base(nodeManager, edgeManager, faceManager, targetRegionIndex, elementSubRegion, finiteElementSpace,
        inputConstitutiveType, inputDofNumber, rankOffset, inputMatrix, inputRhs, inputDt),
    m_X(nodeManager.referencePosition()),
    m_disp(nodeManager.getField<fields::solidMechanics::totalDisplacement>()),
    m_uhat(nodeManager.getField<fields::solidMechanics::incrementalDisplacement>()),
    m_gravityVector{inputGravityVector[0], inputGravityVector[1], inputGravityVector[2]},
    m_density(inputConstitutiveType.getDensity())
  {}

template<typename SUBREGION_TYPE, typename CONSTITUTIVE_TYPE, typename FE_TYPE>
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void ImplicitFiniteStrainQuasiStatic<SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE> :: setup(
    localIndex const k, StackVariables& stack) const
{
  m_finiteElementSpace.template setup<FE_TYPE>(k, m_meshData, stack.feStack);

  localIndex const numSupportPoints = m_finiteElementSpace.template numSupportPoints<FE_TYPE>(stack.feStack);

  stack.numRows = 3 * numSupportPoints;
  stack.numRows = stack.numRows;

  for (localIndex a = 0; a < numSupportPoints; ++a)
  {
    localIndex const localNodeIndex = m_elemsToNodes(k, a);

    for (int i = 0; i < numDofPerTestSupportPoint; ++i)
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

  // Hanyu: not adding stabilization to the local jacobian since this is for FEM
}

template<typename SUBREGION_TYPE, typename CONSTITUTIVE_TYPE, typename FE_TYPE>
template<typename STRESS_MODIFIER>
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void ImplicitFiniteStrainQuasiStatic<SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE>::quadraturePointKernel(
    localIndex const k, localIndex const q, StackVariables & stack, STRESS_MODIFIER && stressModifier) const
{
  real64 dNdX[numNodesPerElem][3];
  real64 const detJxW = m_finiteElementSpace.template getGradN<FE_TYPE>(k, q, stack.xLocal, stack.feStack, dNdX);

  real64 stress[6] = {0};

  typename CONSTITUTIVE_TYPE::KernelWrapper::DiscretizationOps stiffness;

  real64 dUhatdX[3][3] = { {0} };
  real64 dUdX[3][3] = { {0} };
  real64 F[3][3] = { {0} };

  FE_TYPE::gradient(dNdX, stack.u_local, dUdX);
  FE_TYPE::gradient(dNdX, stack.uhat_local, dUhatdX);

  // calculate deformation gradient
  LvArray::tensorOps::copy<3, 3>(F, dUdX);
  LvArray::tensorOps::add<3, 3>(F, dUhatdX);
  LvArray::tensorOps::addIdentity< 3 >(F, 1.0);

  m_constitutiveUpdate.finiteStrainUpdate(k, q, m_dt, F, stress, stiffness);
}

} // namespace solidMechanicsLagrangianFEMKernels

} // namespace geos

#endif // GEOS_PHYSICSSOLVERS_SOLIDMECHANICS_KERNELS_IMPLCITFINITSTRAINQUASISTATIC_IMPL_HPP_