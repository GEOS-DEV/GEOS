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

#include "ImplicitKernelBase.hpp"

/**
 * @file SparsityKernelBase.hpp
 */

#ifndef GEOSX_FINITEELEMENT_SPARSITYKERNELBASE_HPP_
#define GEOSX_FINITEELEMENT_SPARSITYKERNELBASE_HPP_



namespace geosx
{

namespace finiteElement
{

//*****************************************************************************
//*****************************************************************************
//*****************************************************************************
/**
 * @class SparsityKernelBase
 * @brief Define the base interface for implicit finite element kernels.
 * @copydoc geosx::finiteElement::KernelBase
 *
 * ### SparsityKernelBase Description
 * Provide common kernels for generation of CRS Sparsity patterns.
 */
template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          int NUM_TEST_SUPPORT_POINTS_PER_ELEM,
          int NUM_TRIAL_SUPPORT_POINTS_PER_ELEM,
          int NUM_DOF_PER_TEST_SP,
          int NUM_DOF_PER_TRIAL_SP >
class SparsityKernelBase : public ImplicitKernelBase< SUBREGION_TYPE,
                                                      CONSTITUTIVE_TYPE,
                                                      NUM_TEST_SUPPORT_POINTS_PER_ELEM,
                                                      NUM_TRIAL_SUPPORT_POINTS_PER_ELEM,
                                                      NUM_DOF_PER_TEST_SP,
                                                      NUM_DOF_PER_TRIAL_SP,
                                                      ParallelMatrix >
{
public:
  /// Alias for the base class. (i.e. #geosx::finiteElement::ImplicitKernelBase)
  using Base = ImplicitKernelBase< SUBREGION_TYPE,
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
  using Base::m_matrix;

  using typename Base::StackVariables;
  using Base::setup;


  /**
   * @copydoc geosx::finiteElement::ImplicitKernelBase::ImplicitKernelBase
   */
  SparsityKernelBase( NodeManager const & nodeManager,
                      EdgeManager const & edgeManager,
                      FaceManager const & faceManager,
                      SUBREGION_TYPE const & elementSubRegion,
                      FiniteElementBase const * const finiteElementSpace,
                      CONSTITUTIVE_TYPE * const inputConstitutiveType,
                      arrayView1d< globalIndex const > const & inputDofNumber,
                      ParallelMatrix & inputMatrix,
                      ParallelVector & inputRhs ):
    Base( nodeManager,
          edgeManager,
          faceManager,
          elementSubRegion,
          finiteElementSpace,
          inputConstitutiveType,
          inputDofNumber,
          inputMatrix,
          inputRhs )
  {}


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

};


//*****************************************************************************
//*****************************************************************************
//*****************************************************************************
/**
 * @brief Helper struct to define a specialization of
 *   #::geosx::finiteElement::SparsityKernelBase that may be used to generate the sparsity pattern.
 * @tparam KERNEL_TEMPLATE Templated class that defines the physics kernel.
 *                         Most likely derives from SparsityKernelBase.
 */
template< template< typename,
                    typename,
                    int,
                    int > class KERNEL_TEMPLATE >
struct SparsityHelper
{

  /**
   * Defines an alias for the specialization of
   * #geosx::finiteElement::SparsityKernelBase from the compile time constants
   * defined in @p KERNEL_TEMPLATE. Specifically, the
   * NUM_TEST_SUPPORT_POINTS_PER_ELEM and NUM_TRIAL_SUPPORT_POINTS_PER_ELEM
   * parameters are specified by the physics solver.
   */
  template< typename SUBREGION_TYPE,
            typename CONSTITUTIVE_TYPE,
            int NUM_TEST_SUPPORT_POINTS_PER_ELEM,
            int NUM_TRIAL_SUPPORT_POINTS_PER_ELEM >
  using Kernel = SparsityKernelBase< SUBREGION_TYPE,
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
 * #geosx::finiteElement::ImplicitKernelBase alias to specialize
 * #geosx::finiteElement::SparsityKernelBase to conform with the template
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



#endif /* GEOSX_FINITEELEMENT_SPARSITYKERNELBASE_HPP_ */
