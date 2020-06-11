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

#include "KernelBase.hpp"

/**
 * @file ImplicitKernelBase.hpp
 */

#ifndef GEOSX_FINITEELEMENT_IMPLICITKERNELBASE_HPP_
#define GEOSX_FINITEELEMENT_IMPLICITKERNELBASE_HPP_



namespace geosx
{

namespace finiteElement
{

//*****************************************************************************
//*****************************************************************************
//*****************************************************************************
/**
 * @class ImplicitKernelBase
 * @brief Define the base interface for implicit finite element kernels.
 * @copydoc geosx::finiteElement::KernelBase
 *
 * ### ImplicitKernelBase Description
 * Provides a common base for kernels that require the assembly of a system of
 * equations. The types required to assemble the system, such as DOF
 * information, the Matrix and Vector object, etc., are declared and set here.
 */
template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          int NUM_TEST_SUPPORT_POINTS_PER_ELEM,
          int NUM_TRIAL_SUPPORT_POINTS_PER_ELEM,
          int NUM_DOF_PER_TEST_SP,
          int NUM_DOF_PER_TRIAL_SP >
class ImplicitKernelBase : public KernelBase< SUBREGION_TYPE,
                                              CONSTITUTIVE_TYPE,
                                              NUM_TEST_SUPPORT_POINTS_PER_ELEM,
                                              NUM_TRIAL_SUPPORT_POINTS_PER_ELEM,
                                              NUM_DOF_PER_TEST_SP,
                                              NUM_DOF_PER_TRIAL_SP >
{
public:
  /// Alias for the base class. (i.e. #geosx::finiteElement::KernelBase)
  using Base = KernelBase< SUBREGION_TYPE,
                           CONSTITUTIVE_TYPE,
                           NUM_TEST_SUPPORT_POINTS_PER_ELEM,
                           NUM_TRIAL_SUPPORT_POINTS_PER_ELEM,
                           NUM_DOF_PER_TEST_SP,
                           NUM_DOF_PER_TRIAL_SP >;

  using Base::numTestSupportPointsPerElem;
  using Base::numTrialSupportPointsPerElem;
  using Base::numDofPerTestSupportPoint;
  using Base::numDofPerTrialSupportPoint;
  using Base::m_elemsToNodes;


  /**
   * @brief Constructor
   * @param nodeManager Reference to the NodeManager object.
   * @param edgeManager Reference to the EdgeManager object.
   * @param faceManager Reference to the FaceManager object.
   * @param inputDofNumber The dof number for the primary field.
   * @param inputMatrix Reference to the Jacobian matrix.
   * @param inputRhs Reference to the RHS vector.
   * @copydoc geosx::finiteElement::KernelBase::KernelBase
   */
  ImplicitKernelBase( NodeManager const & nodeManager,
                      EdgeManager const & edgeManager,
                      FaceManager const & faceManager,
                      SUBREGION_TYPE const & elementSubRegion,
                      FiniteElementBase const * const finiteElementSpace,
                      CONSTITUTIVE_TYPE * const inputConstitutiveType,
                      arrayView1d< globalIndex const > const & inputDofNumber,
                      ParallelMatrix & inputMatrix,
                      ParallelVector & inputRhs ):
    Base( elementSubRegion,
          finiteElementSpace,
          inputConstitutiveType ),
    m_dofNumber( inputDofNumber ),
    m_matrix( inputMatrix ),
    m_rhs( inputRhs )
  {
    GEOSX_UNUSED_VAR( nodeManager );
    GEOSX_UNUSED_VAR( edgeManager );
    GEOSX_UNUSED_VAR( faceManager );
  }


  //***************************************************************************
  /**
   * @copydoc finiteElement::KernelBase::StackVariables
   */
  struct StackVariables : public Base::StackVariables
  {
    /// The number of rows in the element local jacobian matrix.
    static constexpr int numRows = numTestSupportPointsPerElem *numDofPerTestSupportPoint;

    /// The number of columns in the element local jacobian matrix.
    static constexpr int numCols = numTrialSupportPointsPerElem *numDofPerTrialSupportPoint;

    /**
     * Default constructor
     */
    GEOSX_HOST_DEVICE
    StackVariables():
      localRowDofIndex{ 0 },
      localColDofIndex{ 0 },
      localResidual{ 0.0 },
      localJacobian{ {0.0} }
    {}

    /// C-array storage for the element local row degrees of freedom.
    globalIndex localRowDofIndex[numRows];

    /// C-array storage for the element local column degrees of freedom.
    globalIndex localColDofIndex[numCols];

    /// C-array storage for the element local residual vector.
    real64 localResidual[numRows];

    /// C-array storage for the element local Jacobian matrix.
    real64 localJacobian[numRows][numCols];
  };
  //***************************************************************************

  /**
   * @copydoc geosx::finiteElement::KernelBase::setup
   *
   * ### ImplicitKernelBase::setup() Description
   *
   * In this implementation, the element local Row and Column DOF stack arrays
   * are filled for when we fill the global matrix and rhs.
   *
   * @note This seems like a waste of register space. We should do this in
   *       complete() unless we actually need these dof somewhere else in the kernel.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void setup( localIndex const k,
              StackVariables & stack ) const
  {
    for( localIndex a=0; a<numTestSupportPointsPerElem; ++a )
    {
      localIndex const localNodeIndex = m_elemsToNodes[k][a];
      for( int i=0; i<numDofPerTestSupportPoint; ++i )
      {
        stack.localRowDofIndex[a*numDofPerTestSupportPoint+i] = m_dofNumber[localNodeIndex]+i;
      }
    }

    for( localIndex a=0; a<numTrialSupportPointsPerElem; ++a )
    {
      localIndex const localNodeIndex = m_elemsToNodes[k][a];
      for( int i=0; i<numDofPerTrialSupportPoint; ++i )
      {
        stack.localColDofIndex[a*numDofPerTrialSupportPoint+i] = m_dofNumber[localNodeIndex]+i;
      }
    }
  }

  /**
   * @copydoc geosx::finiteElement::KernelBase::complete
   *
   * In this implementation, only the matrix values are inserted, making this
   * implementation appropriate for generating the sparsity pattern.
   */
//    GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  real64 complete( localIndex const k,
                   StackVariables & stack ) const
  {
    GEOSX_UNUSED_VAR( k );
    m_matrix.insert( stack.localRowDofIndex,
                     stack.localColDofIndex,
                     &(stack.localJacobian[0][0]),
                     stack.numRows,
                     stack.numCols );
    return 0;
  }



protected:
  /// The global degree of freedom number
  arrayView1d< globalIndex const > const m_dofNumber;

  /// The global Jacobian matrix.
  ParallelMatrix & m_matrix;

  /// The global residaul vector.
  ParallelVector & m_rhs;

};


//*****************************************************************************
//*****************************************************************************
//*****************************************************************************
/**
 * @brief Helper struct to define a specialization of
 *   #::geosx::finiteElement::ImplicitKernelBase that may be used to generate the sparsity pattern.
 * @tparam KERNEL_TEMPLATE Templated class that defines the physics kernel.
 *                         Most likely derives from ImplicitKernelBase.
 */
template< template< typename,
                    typename,
                    int,
                    int > class KERNEL_TEMPLATE >
struct SparsityHelper
{

  /**
   * Defines an alias for the specialization of
   * #geosx::finiteElement::ImplicitKernelBase from the compile time constants
   * defined in @p KERNEL_TEMPLATE. Specifically, the
   * NUM_TEST_SUPPORT_POINTS_PER_ELEM and NUM_TRIAL_SUPPORT_POINTS_PER_ELEM
   * parameters are specified by the physics solver.
   */
  template< typename SUBREGION_TYPE,
            typename CONSTITUTIVE_TYPE,
            int NUM_TEST_SUPPORT_POINTS_PER_ELEM,
            int NUM_TRIAL_SUPPORT_POINTS_PER_ELEM >
  using Kernel = ImplicitKernelBase< SUBREGION_TYPE,
                                     CONSTITUTIVE_TYPE,
                                     NUM_TEST_SUPPORT_POINTS_PER_ELEM,
                                     NUM_TRIAL_SUPPORT_POINTS_PER_ELEM,
                                     KERNEL_TEMPLATE< SUBREGION_TYPE,
                                                      CONSTITUTIVE_TYPE,
                                                      NUM_TEST_SUPPORT_POINTS_PER_ELEM,
                                                      NUM_TRIAL_SUPPORT_POINTS_PER_ELEM >::numDofPerTestSupportPoint,
                                     KERNEL_TEMPLATE< SUBREGION_TYPE,
                                                      CONSTITUTIVE_TYPE,
                                                      NUM_TEST_SUPPORT_POINTS_PER_ELEM,
                                                      NUM_TRIAL_SUPPORT_POINTS_PER_ELEM >::numDofPerTrialSupportPoint >;
};


//*****************************************************************************
//*****************************************************************************
//*****************************************************************************
/**
 * @brief Fills matrix sparsity.
 * @tparam POLICY The RAJA launch policy to pass to the kernel launch.
 * @tparam REGION_TYPE The type of region to loop over. TODO make this a
 *                     parameter pack?
 * @tparam KERNEL_TEMPLATE The type of template for the physics kernel, which
 *                         conforms to the interface specified by KernelBase.
 * @param mesh The MeshLevel object.
 * @param targetRegions The names of the target regions(of type @p REGION_TYPE)
 *                      to apply the @p KERNEL_TEMPLATE.
 * @param feDiscretization A pointer to the finite element discretization/space
 *                         object.
 * @param inputDofNumber The global degree of freedom numbers.
 * @param inputMatrix The global Jacobian Matrix.
 * @param inputRhs The global residual vector (unused)....perhaps remove?
 * @return 0
 *
 * Fills matrix sparsity using information from physics specific implementation
 * of #geosx::finiteElement::KernelBase interface using the
 * #geosx::finiteElement::KernelBase alias to specialize
 * #geosx::finiteElement::ImplicitKernelBase to conform with the template
 * pattern specified in the physics kernels.
 */
template< typename POLICY,
          typename REGION_TYPE,
          template< typename SUBREGION_TYPE,
                    typename CONSTITUTIVE_TYPE,
                    int NUM_TEST_SUPPORT_POINTS_PER_ELEM,
                    int NUM_TRIAL_SUPPORT_POINTS_PER_ELEM > class KERNEL_TEMPLATE >
static
real64 fillSparsity( MeshLevel & mesh,
                     arrayView1d< string const > const & targetRegions,
                     FiniteElementDiscretization const * const feDiscretization,
                     arrayView1d< globalIndex const > const & inputDofNumber,
                     ParallelMatrix & inputMatrix,
                     ParallelVector & inputRhs )
{
  real64 rval = 0;

  rval = regionBasedKernelApplication< POLICY,
                                       constitutive::NullModel,
                                       REGION_TYPE,
                                       SparsityHelper< KERNEL_TEMPLATE >::template Kernel
                                       >( mesh,
                                          targetRegions,
                                          array1d< string >(),
                                          feDiscretization,
                                          inputDofNumber,
                                          inputMatrix,
                                          inputRhs );

  return rval;
}

}
}



#endif /* GEOSX_FINITEELEMENT_IMPLICITKERNELBASE_HPP_ */
