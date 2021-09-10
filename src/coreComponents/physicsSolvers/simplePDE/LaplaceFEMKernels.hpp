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
 * @file LaplaceFEMKernels.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_SIMPLEPDE_LAPLACEFEMKERNELS_HPP_
#define GEOSX_PHYSICSSOLVERS_SIMPLEPDE_LAPLACEFEMKERNELS_HPP_

#include "finiteElement/kernelInterface/ImplicitKernelBase.hpp"

namespace geosx
{

//*****************************************************************************
/**
 * @brief Implements kernels for solving Laplace's equation.
 * @copydoc geosx::finiteElement::KernelBase
 * @tparam NUM_NODES_PER_ELEM The number of nodes per element for the
 *                            @p SUBREGION_TYPE.
 * @tparam UNUSED An unused parameter since we are assuming that the test and
 *                trial space have the same number of support points.
 *
 * ### LaplaceFEMKernel Description
 * Implements the KernelBase interface functions required for solving Laplace's
 * equation using on of the finite element kernel application functions such as
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
class LaplaceFEMKernel :
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
  using Base::m_finiteElementSpace;

  /// The number of nodes per element.
  static constexpr int numNodesPerElem = Base::numTestSupportPointsPerElem;

  /**
   * @brief Constructor
   * @copydoc geosx::finiteElement::ImplicitKernelBase::ImplicitKernelBase
   * @param fieldName The name of the primary field
   *                  (i.e. Temperature, Pressure, etc.)
   */
  LaplaceFEMKernel( NodeManager const & nodeManager,
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
                    string const & fieldName ):
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
    m_X( nodeManager.referencePosition() ),
    m_primaryField( nodeManager.template getReference< array1d< real64 > >( fieldName ))
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
            primaryField_local{ 0.0 }
    {}

#if !defined(CALC_FEM_SHAPE_IN_KERNEL)
    /// Dummy
    int xLocal;
#else
    /// C-array stack storage for element local the nodal positions.
    real64 xLocal[ numNodesPerElem ][ 3 ];
#endif

    /// C-array storage for the element local primary field variable.
    real64 primaryField_local[numNodesPerElem];
  };


  /**
   * @brief Copy global values from primary field to a local stack array.
   * @copydoc geosx::finiteElement::ImplicitKernelBase::setup
   *
   * For the LaplaceFEMKernel implementation, global values from the
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
      for( int i=0; i<3; ++i )
      {
        stack.xLocal[ a ][ i ] = m_X[ localNodeIndex ][ i ];
      }
#endif

      stack.primaryField_local[ a ] = m_primaryField[ localNodeIndex ];
      stack.localRowDofIndex[a] = m_dofNumber[localNodeIndex];
      stack.localColDofIndex[a] = m_dofNumber[localNodeIndex];
    }
  }

  /**
   * @copydoc geosx::finiteElement::ImplicitKernelBase::quadraturePointKernel
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void quadraturePointKernel( localIndex const k,
                              localIndex const q,
                              StackVariables & stack ) const
  {
    real64 dNdX[ numNodesPerElem ][ 3 ];
    real64 const detJ = m_finiteElementSpace.template getGradN< FE_TYPE >( k, q, stack.xLocal, dNdX );

    for( localIndex a=0; a<numNodesPerElem; ++a )
    {
      for( localIndex b=0; b<numNodesPerElem; ++b )
      {
        stack.localJacobian[ a ][ b ] += LvArray::tensorOps::AiBi< 3 >( dNdX[a], dNdX[b] ) * detJ;
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

    for( localIndex a = 0; a < numNodesPerElem; ++a )
    {
      for( localIndex b = 0; b < numNodesPerElem; ++b )
      {
        stack.localResidual[ a ] += stack.localJacobian[ a ][ b ] * stack.primaryField_local[ b ];
      }
    }

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
  arrayView1d< real64 const > const m_primaryField;

};

/// The factory used to construct a LaplaceFEMKernel.
using LaplaceFEMKernelFactory = finiteElement::KernelFactory< LaplaceFEMKernel,
                                                              arrayView1d< globalIndex const > const &,
                                                              globalIndex, CRSMatrixView< real64, globalIndex const > const &,
                                                              arrayView1d< real64 > const &, string const & >;

} // namespace geosx

#include "finiteElement/kernelInterface/SparsityKernelBase.hpp"

#endif // GEOSX_PHYSICSSOLVERS_SIMPLEPDE_LAPLACEFEMKERNELS_HPP_
