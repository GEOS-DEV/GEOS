/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
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

/**
 * @struct ConstructorParams
 * @brief Temporary object to help with the delivery of the variadic
 *        arguments through all the lambdas prior to the reaching the
 *        constructor for the Kernel. This should not be necessary once
 *        we are able to use the perfect forwarding capabilites of c++20.
 */
struct LaplaceFEMKernelConstructorParams : finiteElement::ImplicitKernelBaseConstructorParams
{
  LaplaceFEMKernelConstructorParams( arrayView1d< globalIndex const > const & inputDofNumber,
                                     ParallelMatrix & inputMatrix,
                                     ParallelVector & inputRhs,
                                     string const & fieldName ):
    finiteElement::ImplicitKernelBaseConstructorParams( inputDofNumber,
                                                        inputMatrix,
                                                        inputRhs ),
    m_fieldName(fieldName)
  {}
  string const & m_fieldName;
};

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
          int NUM_NODES_PER_ELEM,
          int UNUSED >
class LaplaceFEMKernel :
  public finiteElement::ImplicitKernelBase< SUBREGION_TYPE,
                                            CONSTITUTIVE_TYPE,
                                            NUM_NODES_PER_ELEM,
                                            NUM_NODES_PER_ELEM,
                                            1,
                                            1 >
{
public:
  /// An alias for the base class.
  using Base = finiteElement::ImplicitKernelBase< SUBREGION_TYPE,
                                                  CONSTITUTIVE_TYPE,
                                                  NUM_NODES_PER_ELEM,
                                                  NUM_NODES_PER_ELEM,
                                                  1,
                                                  1 >;

  /// Number of nodes per element...which is equal to the
  /// numTestSupportPointPerElem and numTrialSupportPointPerElem by definition.
  static constexpr int numNodesPerElem = NUM_NODES_PER_ELEM;

  /// @copydoc geosx::finiteElement::ImplicitKernelBase::numDofPerTestSupportPoint
  using Base::numDofPerTestSupportPoint;

  /// @copydoc geosx::finiteElement::ImplicitKernelBase::numDofPerTrialSupportPoint
  using Base::numDofPerTrialSupportPoint;

  /// @copydoc geosx::finiteElement::ImplicitKernelBase::m_dofNumber
  using Base::m_dofNumber;

  /// @copydoc geosx::finiteElement::ImplicitKernelBase::m_matrix
  using Base::m_matrix;

  /// @copydoc geosx::finiteElement::ImplicitKernelBase::m_rhs
  using Base::m_rhs;

  /// @copydoc geosx::finiteElement::ImplicitKernelBase::elemsToNodes
  using Base::elemsToNodes;

  /// @copydoc geosx::finiteElement::ImplicitKernelBase::elemGhostRank
  using Base::elemGhostRank;

  /// @copydoc geosx::finiteElement::ImplicitKernelBase::constitutiveUpdate
  using Base::constitutiveUpdate;

  /// @copydoc geosx::finiteElement::ImplicitKernelBase::m_finiteElementSpace
  using Base::m_finiteElementSpace;

  /// @copydoc geosx::finiteElement::ImplicitKernelBase::Launch
  using Base::Launch;

  using ConstructorParams = LaplaceFEMKernelConstructorParams;

  /**
   * @brief Constructor
   * @copydoc geosx::finiteElement::ImplicitKernelBase::ImplicitKernelBase
   * @param fieldName The name of the primary field
   *                  (i.e. Temperature, Pressure, etc.)
   */
  LaplaceFEMKernel( NodeManager const & nodeManager,
                    EdgeManager const & edgeManager,
                    FaceManager const & faceManager,
                    SUBREGION_TYPE const & elementSubRegion,
                    FiniteElementBase const * const finiteElementSpace,
                    CONSTITUTIVE_TYPE * const inputConstitutiveType,
                    arrayView1d< globalIndex const > const & inputDofNumber,
                    ParallelMatrix & inputMatrix,
                    ParallelVector & inputRhs,
                    string const & fieldName ):
    Base( nodeManager,
          edgeManager,
          faceManager,
          elementSubRegion,
          finiteElementSpace,
          inputConstitutiveType,
          inputDofNumber,
          inputMatrix,
          inputRhs ),
    m_primaryField( nodeManager.template getReference< array1d< real64 > >( fieldName )),
    dNdX( elementSubRegion.template getReference< array3d< R1Tensor > >( dataRepository::keys::dNdX )),
    detJ( elementSubRegion.template getReference< array2d< real64 > >( dataRepository::keys::detJ ) )  //,
  {}

  LaplaceFEMKernel( NodeManager const & nodeManager,
                     EdgeManager const & edgeManager,
                     FaceManager const & faceManager,
                     SUBREGION_TYPE const & elementSubRegion,
                     FiniteElementBase const * const finiteElementSpace,
                     CONSTITUTIVE_TYPE * const inputConstitutiveType,
                     ConstructorParams & params ):
    LaplaceFEMKernel( nodeManager,
                      edgeManager,
                      faceManager,
                      elementSubRegion,
                      finiteElementSpace,
                      inputConstitutiveType,
                      params.m_dofNumber,
                      params.m_matrix,
                      params.m_rhs,
                      params.m_fieldName )
  {}


  /**
   * @brief Copy global values from primary field to a local stack array.
   * @copydoc geosx::finiteElement::ImplicitKernelBase::setup.
   *
   * For the LaplaceFEMKernel implementation, global values from the
   * primaryField, and degree of freedom numbers are placed into element local
   * stack storage.
   */
  template< typename STACK_VARIABLE_TYPE >
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void setup( localIndex const k,
              STACK_VARIABLE_TYPE & stack ) const
  {
    for( localIndex a=0; a<NUM_NODES_PER_ELEM; ++a )
    {
      localIndex const localNodeIndex = elemsToNodes( k, a );

      stack.primaryField_local[ a ] = m_primaryField[ localNodeIndex ];
      stack.localRowDofIndex[a] = m_dofNumber[localNodeIndex];
      stack.localColDofIndex[a] = m_dofNumber[localNodeIndex];
    }
  }

  /**
   * @copydoc geosx::finiteElement::ImplicitKernelBase::quadraturePointJacobianContribution.
   */
  template< typename STACK_VARIABLE_TYPE,
            typename DYNAMICS_LAMBDA = std::function< void( localIndex, localIndex) > >
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void quadraturePointJacobianContribution( localIndex const k,
                                            localIndex const q,
                                            STACK_VARIABLE_TYPE & stack ) const
  {
    for( localIndex a=0; a<NUM_NODES_PER_ELEM; ++a )
    {
      for( localIndex b=0; b<NUM_NODES_PER_ELEM; ++b )
      {
        stack.localJacobian[ a ][ b ] += Dot( dNdX( k, q, a ), dNdX( k, q, b ) ) * detJ( k, q );
      }
    }
  }

  /**
   * @copydoc geosx::finiteElement::ImplicitKernelBase::complete.
   *
   * Form element residual from the fully formed element Jacobian dotted with
   * the primary field and map the element local Jacobian/Residual to the
   * global matrix/vector.
   */
  template< typename STACK_VARIABLE_TYPE >
  //GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  real64 complete( localIndex const GEOSX_UNUSED_PARAM( k ),
                   STACK_VARIABLE_TYPE & stack ) const
  {
    for( localIndex a = 0; a < NUM_NODES_PER_ELEM; ++a )
    {
      for( localIndex b = 0; b < NUM_NODES_PER_ELEM; ++b )
      {
        stack.localResidual[ a ] += stack.localJacobian[ a ][ b ] * stack.primaryField_local[ b ];
      }
    }

    m_matrix.add( stack.localRowDofIndex,
                  stack.localColDofIndex,
                  &(stack.localJacobian[0][0]),
                  stack.numRows,
                  stack.numCols );

    m_rhs.add( stack.localRowDofIndex,
               stack.localResidual,
               stack.numRows );

    return 1.0;
  }


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

    /// @copydoc geosx::finiteElement::ImplicitKernelBase::numRows
    using Base::StackVariables::numRows;

    /// @copydoc geosx::finiteElement::ImplicitKernelBase::numCols
    using Base::StackVariables::numCols;

    /**
     * @brief Constructor
     */
    GEOSX_HOST_DEVICE
    StackVariables():
      Base::StackVariables(),
      primaryField_local{ 0.0 }
    {}

    /// C-array storage for the element local primary field variable.
    real64 primaryField_local[numNodesPerElem];
  };


protected:
  /// The global primary field array.
  arrayView1d< real64 const > const m_primaryField;

  /// The global shape function derivatives array.
  arrayView3d< R1Tensor const > const dNdX;

  /// The global determinant of the parent/physical Jacobian.
  arrayView2d< real64 const > const detJ;

};


} // namespace geosx
#endif // GEOSX_PHYSICSSOLVERS_SIMPLEPDE_LAPLACEFEMKERNELS_HPP_
