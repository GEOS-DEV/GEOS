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
 * @file SolidMechanicsSmallStrainQuasiStaticKernel.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSSMALLSTRAINQUASISTATIC_HPP_
#define GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSSMALLSTRAINQUASISTATIC_HPP_

#include "finiteElement/kernelInterface/ImplicitKernelBase.hpp"

namespace geosx
{

namespace solidMechanicsLagrangianFEMKernels
{

/**
 * @brief Implements kernels for solving quasi-static equilibrium.
 * @copydoc geosx::finiteElement::ImplicitKernelBase
 * @tparam NUM_NODES_PER_ELEM The number of nodes per element for the
 *                            @p SUBREGION_TYPE.
 * @tparam UNUSED An unused parameter since we are assuming that the test and
 *                trial space have the same number of support points.
 *
 * ### QuasiStatic Description
 * Implements the KernelBase interface functions required for solving the
 * quasi-static equilibrium equations using one of the
 * "finite element kernel application" functions such as
 * geosx::finiteElement::RegionBasedKernelApplication.
 *
 * In this implementation, the template parameter @p NUM_NODES_PER_ELEM is used
 * in place of both @p NUM_TEST_SUPPORT_POINTS_PER_ELEM and
 * @p NUM_TRIAL_SUPPORT_POINTS_PER_ELEM, which are assumed to be equal. This
 * results in the @p UNUSED template parameter as only the NUM_NODES_PER_ELEM
 * is passed to the ImplicitKernelBase template to form the base class.
 *
 * Additionally, the number of degrees of freedom per support point for both
 * the test and trial spaces are specified as `3` when specifying the base
 * class.
 */
template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
class QuasiStatic :
  public finiteElement::ImplicitKernelBase< SUBREGION_TYPE,
                                            CONSTITUTIVE_TYPE,
                                            FE_TYPE,
                                            3,
                                            3 >
{
public:
  /// Alias for the base class;
  using Base = finiteElement::ImplicitKernelBase< SUBREGION_TYPE,
                                                  CONSTITUTIVE_TYPE,
                                                  FE_TYPE,
                                                  3,
                                                  3 >;

  /// Maximum number of nodes per element, which is equal to the maxNumTestSupportPointPerElem and
  /// maxNumTrialSupportPointPerElem by definition. When the FE_TYPE is not a Virtual Element, this
  /// will be the actual number of nodes per element.
  static constexpr int numNodesPerElem = Base::maxNumTestSupportPointsPerElem;
  using Base::numDofPerTestSupportPoint;
  using Base::numDofPerTrialSupportPoint;
  using Base::m_dofNumber;
  using Base::m_dofRankOffset;
  using Base::m_matrix;
  using Base::m_rhs;
  using Base::m_elemsToNodes;
  using Base::m_constitutiveUpdate;
  using Base::m_finiteElementSpace;
  using Base::m_meshData;


  /**
   * @brief Constructor
   * @copydoc geosx::finiteElement::ImplicitKernelBase::ImplicitKernelBase
   * @param inputGravityVector The gravity vector.
   */
  QuasiStatic( NodeManager const & nodeManager,
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
    m_disp( nodeManager.totalDisplacement()),
    m_uhat( nodeManager.incrementalDisplacement()),
    m_gravityVector{ inputGravityVector[0], inputGravityVector[1], inputGravityVector[2] },
    m_density( inputConstitutiveType.getDensity() )
  {}


  //*****************************************************************************
  /**
   * @class StackVariables
   * @copydoc geosx::finiteElement::ImplicitKernelBase::StackVariables
   *
   * Adds a stack array for the displacement, incremental displacement, and the
   * constitutive stiffness.
   */
  struct StackVariables : public Base::StackVariables
  {
public:

    /// Constructor.
    GEOSX_HOST_DEVICE
    StackVariables():
      Base::StackVariables(),
                                       xLocal(),
                                       u_local(),
                                       uhat_local(),
                                       constitutiveStiffness()
    {}

#if !defined(CALC_FEM_SHAPE_IN_KERNEL)
    /// Dummy
    int xLocal;
#else
    /// C-array stack storage for element local the nodal positions.
    real64 xLocal[ numNodesPerElem ][ 3 ];
#endif

    /// Stack storage for the element local nodal displacement
    real64 u_local[numNodesPerElem][numDofPerTrialSupportPoint];

    /// Stack storage for the element local nodal incremental displacement
    real64 uhat_local[numNodesPerElem][numDofPerTrialSupportPoint];

    /// Stack storage for the constitutive stiffness at a quadrature point.
    real64 constitutiveStiffness[ 6 ][ 6 ];
  };
  //*****************************************************************************

  /**
   * @brief Copy global values from primary field to a local stack array.
   * @copydoc ::geosx::finiteElement::ImplicitKernelBase::setup
   *
   * For the QuasiStatic implementation, global values from the displacement,
   * incremental displacement, and degree of freedom numbers are placed into
   * element local stack storage.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void setup( localIndex const k,
              StackVariables & stack ) const
  {
    m_finiteElementSpace.template setup< FE_TYPE >( k, m_meshData, stack.feStack );
    localIndex const numSupportPoints =
      m_finiteElementSpace.template numSupportPoints< FE_TYPE >( stack.feStack );
    stack.numRows =  3 * numSupportPoints;
    stack.numCols = stack.numRows;
    for( localIndex a = 0; a < numSupportPoints; ++a )
    {
      localIndex const localNodeIndex = m_elemsToNodes( k, a );

      for( int i = 0; i < 3; ++i )
      {
#if defined(CALC_FEM_SHAPE_IN_KERNEL)
        stack.xLocal[ a ][ i ] = m_X[ localNodeIndex ][ i ];
#endif
        stack.u_local[ a ][ i ] = m_disp[ localNodeIndex ][ i ];
        stack.uhat_local[ a ][ i ] = m_uhat[ localNodeIndex ][ i ];
        stack.localRowDofIndex[ a*3+i ] = m_dofNumber[ localNodeIndex ] + i;
        stack.localColDofIndex[ a*3+i ] = m_dofNumber[ localNodeIndex ] + i;
      }
    }
    // Add stabilization to block diagonal parts of the local jacobian
    // (this is a no-operation with FEM classes)
    real64 const stabilizationScaling = computeStabilizationScaling( k );
    m_finiteElementSpace.template addGradGradStabilizationMatrix
    < FE_TYPE, numDofPerTrialSupportPoint, true >( stack.feStack,
                                                   stack.localJacobian,
                                                   -stabilizationScaling );
  }


  /**
   * @brief Internal struct to provide no-op defaults used in the inclusion
   *   of lambda functions into kernel component functions.
   * @struct NoOpFunctors
   */
  struct NoOpFunctors
  {
    /**
     * @brief operator() no-op used for adding an additional dynamics term
     *   inside the jacobian assembly loop.
     * @param a Node index for the row.
     * @param b Node index for the col.
     */
    GEOSX_HOST_DEVICE GEOSX_FORCE_INLINE constexpr
    void operator() ( localIndex const a, localIndex const b )
    {
      GEOSX_UNUSED_VAR( a );
      GEOSX_UNUSED_VAR( b );
    }

    /**
     * @brief operator() no-op used for modifying the stress tensor prior to
     *   integrating the divergence to produce nodal forces.
     * @param stress The stress array.
     */
    GEOSX_HOST_DEVICE GEOSX_FORCE_INLINE constexpr
    void operator() ( real64 (& stress)[6] )
    {
      GEOSX_UNUSED_VAR( stress );
    }
  };

  /**
   * @copydoc geosx::finiteElement::KernelBase::quadraturePointKernel
   * @tparam STRESS_MODIFIER Type of optional functor to allow for the
   * modification of stress prior to integration.
   * @param stressModifier An optional functor to allow for the modification
   *  of stress prior to integration.
   * For solid mechanics kernels, the strain increment is calculated, and the
   * constitutive update is called. In addition, the constitutive stiffness
   * stack variable is filled by the constitutive model.
   */
  template< typename STRESS_MODIFIER = NoOpFunctors >
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void quadraturePointKernel( localIndex const k,
                              localIndex const q,
                              StackVariables & stack,
                              STRESS_MODIFIER && stressModifier = NoOpFunctors{} ) const
  {
    real64 dNdX[ numNodesPerElem ][ 3 ];
    real64 const detJxW = m_finiteElementSpace.template getGradN< FE_TYPE >( k, q, stack.xLocal,
                                                                             stack.feStack, dNdX );

    real64 strainInc[6] = {0};
    real64 stress[6] = {0};

    typename CONSTITUTIVE_TYPE::KernelWrapper::DiscretizationOps stiffness;

    FE_TYPE::symmetricGradient( dNdX, stack.uhat_local, strainInc );

    m_constitutiveUpdate.smallStrainUpdate( k, q, strainInc, stress, stiffness );

    stressModifier( stress );
    for( localIndex i=0; i<6; ++i )
    {
      stress[i] *= -detJxW;
    }

    real64 const gravityForce[3] = { m_gravityVector[0] * m_density( k, q )* detJxW,
                                     m_gravityVector[1] * m_density( k, q )* detJxW,
                                     m_gravityVector[2] * m_density( k, q )* detJxW };

    real64 N[numNodesPerElem];
    FE_TYPE::calcN( q, stack.feStack, N );
    FE_TYPE::plusGradNajAijPlusNaFi( dNdX,
                                     stress,
                                     N,
                                     gravityForce,
                                     reinterpret_cast< real64 (&)[numNodesPerElem][3] >(stack.localResidual) );
    real64 const stabilizationScaling = computeStabilizationScaling( k );
    m_finiteElementSpace.template
    addEvaluatedGradGradStabilizationVector< FE_TYPE,
                                             numDofPerTrialSupportPoint >( stack.feStack,
                                                                           stack.uhat_local,
                                                                           reinterpret_cast< real64 (&)[numNodesPerElem][3] >(stack.localResidual),
                                                                           -stabilizationScaling );
    stiffness.template upperBTDB< numNodesPerElem >( dNdX, -detJxW, stack.localJacobian );
  }

  /**
   * @copydoc geosx::finiteElement::ImplicitKernelBase::complete
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  real64 complete( localIndex const k,
                   StackVariables & stack ) const
  {
    GEOSX_UNUSED_VAR( k );
    real64 maxForce = 0;

    // TODO: Does this work if BTDB is non-symmetric?
    CONSTITUTIVE_TYPE::KernelWrapper::DiscretizationOps::template fillLowerBTDB< numNodesPerElem >( stack.localJacobian );
    localIndex const numSupportPoints =
      m_finiteElementSpace.template numSupportPoints< FE_TYPE >( stack.feStack );
    for( int localNode = 0; localNode < numSupportPoints; ++localNode )
    {
      for( int dim = 0; dim < numDofPerTestSupportPoint; ++dim )
      {
        localIndex const dof =
          LvArray::integerConversion< localIndex >( stack.localRowDofIndex[ numDofPerTestSupportPoint * localNode + dim ] - m_dofRankOffset );
        if( dof < 0 || dof >= m_matrix.numRows() ) continue;
        m_matrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( dof,
                                                                                stack.localRowDofIndex,
                                                                                stack.localJacobian[ numDofPerTestSupportPoint * localNode + dim ],
                                                                                stack.numRows );

        RAJA::atomicAdd< parallelDeviceAtomic >( &m_rhs[ dof ], stack.localResidual[ numDofPerTestSupportPoint * localNode + dim ] );
        maxForce = fmax( maxForce, fabs( stack.localResidual[ numDofPerTestSupportPoint * localNode + dim ] ) );
      }
    }


    return maxForce;
  }



protected:
  /// The array containing the nodal position array.
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const m_X;

  /// The rank-global displacement array.
  arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const m_disp;

  /// The rank-global incremental displacement array.
  arrayView2d< real64 const, nodes::INCR_DISPLACEMENT_USD > const m_uhat;

  /// The gravity vector.
  real64 const m_gravityVector[3];

  /// The rank global density
  arrayView2d< real64 const > const m_density;

  /**
   * @brief Get a parameter representative of the stiffness, used as physical scaling for the
   * stabilization matrix.
   * @param[in] k Element index.
   * @return A parameter representative of the stiffness matrix dstress/dstrain
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  real64 computeStabilizationScaling( localIndex const k ) const
  {
    // TODO: generalize this to other constitutive models (currently we assume linear elasticity).
    return 2.0 * m_constitutiveUpdate.getShearModulus( k );
  }
};

/// The factory used to construct a QuasiStatic kernel.
using QuasiStaticFactory = finiteElement::KernelFactory< QuasiStatic,
                                                         arrayView1d< globalIndex const > const,
                                                         globalIndex,
                                                         CRSMatrixView< real64, globalIndex const > const,
                                                         arrayView1d< real64 > const,
                                                         real64 const (&)[3] >;

} // namespace solidMechanicsLagrangianFEMKernels

} // namespace geosx

#include "finiteElement/kernelInterface/SparsityKernelBase.hpp"

#endif // GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSSMALLSTRAINQUASISTATIC_HPP_
