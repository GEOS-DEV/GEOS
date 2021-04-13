/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file PhaseFieldDamageKernels.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_SIMPLEPDE_PHASEFIELDDAMAGEKERNELS_HPP_
#define GEOSX_PHYSICSSOLVERS_SIMPLEPDE_PHASEFIELDDAMAGEKERNELS_HPP_

#include "finiteElement/kernelInterface/ImplicitKernelBase.hpp"

namespace geosx
{
//*****************************************************************************
/**
 * @brief Implements kernels for solving the Damage(or phase-field) equation
 * in a phase-field fracture problem.
 * @copydoc geosx::finiteElement::KernelBase
 * @tparam NUM_NODES_PER_ELEM The number of nodes per element for the
 *                            @p SUBREGION_TYPE.
 * @tparam UNUSED An unused parameter since we are assuming that the test and
 *                trial space have the same number of support points.
 *
 * ### PhaseFieldDamageKernel Description
 * Implements the KernelBase interface functions required for solving the
 * Damage(or phase-field) equation in a phase-field fracture problem.
 * It uses the finite element kernel application functions such as
 * geosx::finiteElement::RegionBasedKernelApplication.
 *
 * In this implementation, the template parameter @p NUM_NODES_PER_ELEM is used
 * in place of both @p NUM_TEST_SUPPORT_POINTS_PER_ELEM and
 * @p NUM_TRIAL_SUPPORT_POINTS_PER_ELEM, which are assumed to be equal. This
 * results in the @p UNUSED template parameter as only the NUM_NODES_PER_ELEM
 * is passed to the ImplicitKernelBase template to form the base class.
 *
 * Additionally, the number of degrees of freedom per support point for both
 * the test and trial spaces are specified as `1` when specifying the base
 * class.
 */
template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
class PhaseFieldDamageKernel :
  public finiteElement::ImplicitKernelBase< SUBREGION_TYPE,
                                            CONSTITUTIVE_TYPE,
                                            FE_TYPE,
                                            1,
                                            1 >
{
public:
  /// An alias for the base class.
  using Base = finiteElement::ImplicitKernelBase< SUBREGION_TYPE,
                                                  CONSTITUTIVE_TYPE,
                                                  FE_TYPE,
                                                  1,
                                                  1 >;

  using Base::numDofPerTestSupportPoint;
  using Base::numDofPerTrialSupportPoint;
  using Base::m_dofNumber;
  using Base::m_dofRankOffset;
  using Base::m_matrix;
  using Base::m_rhs;
  using Base::m_elemsToNodes;
  using Base::m_constitutiveUpdate;
  using Base::m_finiteElementSpace;

  /// The number of nodes per element.
  static constexpr int numNodesPerElem = Base::numTestSupportPointsPerElem;

  /**
   * @brief Constructor
   * @copydoc geosx::finiteElement::ImplicitKernelBase::ImplicitKernelBase
   * @param fieldName The name of the primary field
   *                  (i.e. Temperature, Pressure, etc.)
   */
  PhaseFieldDamageKernel( NodeManager const & nodeManager,
                          EdgeManager const & edgeManager,
                          FaceManager const & faceManager,
                          localIndex const targetRegionIndex,
                          SUBREGION_TYPE const & elementSubRegion,
                          FE_TYPE const & finiteElementSpace,
                          CONSTITUTIVE_TYPE & inputConstitutiveType,
                          arrayView1d< globalIndex const > const & inputDofNumber,
                          globalIndex const rankOffset,
                          CRSMatrixView< real64, globalIndex const > const & inputMatrix,
                          arrayView1d< real64 > const & inputRhs,
                          string const & fieldName,
                          int const localDissipationOption ):
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
    m_nodalDamage( nodeManager.template getReference< array1d< real64 > >( fieldName )),
    m_localDissipationOption( localDissipationOption )
  {}

  //***************************************************************************
  /**
   * @class StackVariables
   * @copydoc geosx::finiteElement::ImplicitKernelBase::StackVariables
   *
   * Adds a stack array for the primary field.
   */
  struct StackVariables : Base::StackVariables
  {
public:

    /**
     * @brief Constructor
     */
    GEOSX_HOST_DEVICE
    StackVariables():
      Base::StackVariables(),
            xLocal(),
            nodalDamageLocal{ 0.0 }
    {}

#if !defined(CALC_FEM_SHAPE_IN_KERNEL)
    /// Dummy
    int xLocal;
#else
    /// C-array stack storage for element local the nodal positions.
    real64 xLocal[ numNodesPerElem ][ 3 ];
#endif

    /// C-array storage for the element local primary field variable.
    real64 nodalDamageLocal[numNodesPerElem];
  };


  /**
   * @brief Copy global values from primary field to a local stack array.
   * @copydoc geosx::finiteElement::ImplicitKernelBase::setup
   *
   * For the PhaseFieldDamageKernel implementation, global values from the
   * primaryField, and degree of freedom numbers are placed into element local
   * stack storage.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void setup( localIndex const k,
              StackVariables & stack ) const
  {
    for( localIndex a=0; a<numNodesPerElem; ++a )
    {
      localIndex const localNodeIndex = m_elemsToNodes( k, a );

#if defined(CALC_FEM_SHAPE_IN_KERNEL)
      LvArray::tensorOps::copy< 3 >( stack.xLocal[ a ], m_X[ localNodeIndex ] );
#endif

      stack.nodalDamageLocal[ a ] = m_nodalDamage[ localNodeIndex ];
      stack.localRowDofIndex[a] = m_dofNumber[localNodeIndex];
      stack.localColDofIndex[a] = m_dofNumber[localNodeIndex];
    }
  }

  /**
   * @copydoc geosx::finiteElement::ImplicitKernelBase::quadraturePointJacobianContribution
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void quadraturePointKernel( localIndex const k,
                              localIndex const q,
                              StackVariables & stack ) const
  {

    real64 const strainEnergyDensity = m_constitutiveUpdate.getStrainEnergyDensity( k, q );
    real64 const ell = m_constitutiveUpdate.getRegularizationLength();
    real64 const Gc = m_constitutiveUpdate.getCriticalFractureEnergy();
    real64 const threshold = m_constitutiveUpdate.getEnergyThreshold();

    real64 D = 0;                                                                   //max between threshold and
                                                                                    // Elastic energy
    if( m_localDissipationOption == 1 )
    {
      D = fmax( threshold, strainEnergyDensity );
    }

    //Interpolate d and grad_d
    real64 N[ numNodesPerElem ];
    real64 dNdX[ numNodesPerElem ][ 3 ];
    real64 const detJ = m_finiteElementSpace.template getGradN< FE_TYPE >( k, q, stack.xLocal, dNdX );
    FE_TYPE::calcN( q, N );

    real64 qp_damage = 0.0;
    real64 qp_grad_damage[3] = {0, 0, 0};
    FE_TYPE::valueAndGradient( N, dNdX, stack.nodalDamageLocal, qp_damage, qp_grad_damage );

    for( localIndex a = 0; a < numNodesPerElem; ++a )
    {
      if( m_localDissipationOption == 1 )
      {
        stack.localResidual[ a ] += detJ * ( -3 * N[a] / 16  -
                                             0.375*pow( ell, 2 ) * LvArray::tensorOps::AiBi< 3 >( qp_grad_damage, dNdX[a] ) -
                                             (0.5 * ell * D/Gc) * N[a] * m_constitutiveUpdate.getDegradationDerivative( qp_damage ));
      }
      else
      {
        stack.localResidual[ a ] -= detJ * ( N[a] * qp_damage +
                                             (pow( ell, 2 ) * LvArray::tensorOps::AiBi< 3 >( qp_grad_damage, dNdX[a] ) +
                                              N[a] * m_constitutiveUpdate.getDegradationDerivative( qp_damage ) * (ell*strainEnergyDensity/Gc)) );

      }
      for( localIndex b = 0; b < numNodesPerElem; ++b )
      {
        if( m_localDissipationOption == 1 )
        {
          stack.localJacobian[ a ][ b ] -= detJ *
                                           (0.375*pow( ell, 2 ) * LvArray::tensorOps::AiBi< 3 >( dNdX[a], dNdX[b] ) +
                                            (0.5 * ell * D/Gc) * m_constitutiveUpdate.getDegradationSecondDerivative( qp_damage ) * N[a] * N[b]);

        }
        else
        {
          stack.localJacobian[ a ][ b ] -= detJ *
                                           ( pow( ell, 2 ) * LvArray::tensorOps::AiBi< 3 >( dNdX[a], dNdX[b] ) +
                                             N[a] * N[b] * (1 + m_constitutiveUpdate.getDegradationSecondDerivative( qp_damage ) * ell*strainEnergyDensity/Gc )
                                           );
        }
      }
    }
  }

  /**
   * @copydoc geosx::finiteElement::ImplicitKernelBase::complete
   *
   * Form element residual from the fully formed element Jacobian dotted with
   * the primary field and map the element local Jacobian/Residual to the
   * global matrix/vector.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  real64 complete( localIndex const k,
                   StackVariables & stack ) const
  {
    GEOSX_UNUSED_VAR( k );
    real64 maxForce = 0;

    for( int a = 0; a < numNodesPerElem; ++a )
    {
      localIndex const dof = LvArray::integerConversion< localIndex >( stack.localRowDofIndex[ a ] - m_dofRankOffset );
      if( dof < 0 || dof >= m_matrix.numRows() ) continue;
      m_matrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( dof,
                                                                              stack.localColDofIndex,
                                                                              stack.localJacobian[ a ],
                                                                              numNodesPerElem );

      RAJA::atomicAdd< parallelDeviceAtomic >( &m_rhs[ dof ], stack.localResidual[ a ] );
      maxForce = fmax( maxForce, fabs( stack.localResidual[ a ] ) );
    }

    return maxForce;
  }



protected:
  /// The array containing the nodal position array.
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const m_X;

  /// The global primary field array.
  arrayView1d< real64 const > const m_nodalDamage;

  int const m_localDissipationOption;

};


} // namespace geosx

#include "finiteElement/kernelInterface/SparsityKernelBase.hpp"

#endif // GEOSX_PHYSICSSOLVERS_SIMPLEPDE_PHASEFIELDDAMAGEKERNELS_HPP_
