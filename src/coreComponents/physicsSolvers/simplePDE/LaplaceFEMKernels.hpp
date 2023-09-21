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

#ifndef GEOS_PHYSICSSOLVERS_SIMPLEPDE_LAPLACEFEMKERNELSNEW_HPP_
#define GEOS_PHYSICSSOLVERS_SIMPLEPDE_LAPLACEFEMKERNELSNEW_HPP_

#include "finiteElement/kernelInterface/ImplicitKernelBase.hpp"
#include "BasisFunctionUtilities.hpp"
#include "SubRegionMeshUtilities.hpp"
#include "QuadratureUtilities.hpp"
#include "finiteElement/BilinearFormUtilities.hpp"

namespace geos
{

//*****************************************************************************
/**
 * @brief Implements kernels for solving Laplace's equation.
 * @copydoc geos::finiteElement::KernelBase
 * @tparam NUM_NODES_PER_ELEM The number of nodes per element for the
 *                            @p SUBREGION_TYPE.
 * @tparam UNUSED An unused parameter since we are assuming that the test and
 *                trial space have the same number of support points.
 *
 * ### LaplaceFEMKernel Description
 * Implements the KernelBase interface functions required for solving Laplace's
 * equation using on of the finite element kernel application functions such as
 * geos::finiteElement::RegionBasedKernelApplication.
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
  using Base::maxNumTestSupportPointsPerElem;
  using Base::numQuadraturePointsPerElem;
  using Base::m_dofNumber;
  using Base::m_dofRankOffset;
  using Base::m_matrix;
  using Base::m_rhs;
  using Base::m_elemsToNodes;
  using Base::m_finiteElementSpace;
  using Base::m_meshData;

  /**
   * @brief Constructor
   * @copydoc geos::finiteElement::ImplicitKernelBase::ImplicitKernelBase
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
                    arrayView1d< globalIndex const > const inputDofNumber,
                    globalIndex const rankOffset,
                    CRSMatrixView< real64, globalIndex const > const inputMatrix,
                    arrayView1d< real64 > const inputRhs,
                    real64 const inputDt,
                    string const fieldName ):
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
          inputRhs,
          inputDt ),
    // m_X( nodeManager.referencePosition() ),
    m_primaryField( nodeManager.template getReference< array1d< real64 > >( fieldName )),
    // m_subregionMesh( selectIsoparametricMesh< traits::ViewTypeConst< typename SUBREGION_TYPE::NodeMapType::base_type > >(
    //                    elementSubRegion.getElementType(),
    //                    nodeManager.referencePosition(),
    //                    elementSubRegion.nodeList().toViewConst()) )
    m_subregionMesh( nodeManager.referencePosition(),
                     elementSubRegion.nodeList().toViewConst() )
  {}

  //***************************************************************************
  /**
   * @class StackVariables
   * @copydoc geos::finiteElement::ImplicitKernelBase::StackVariables
   *
   * Adds a stack array for the primary field.
   */
  struct StackVariables : Base::StackVariables
  {
public:

    /**
     * @brief Constructor
     */
    GEOS_HOST_DEVICE
    StackVariables():
      Base::StackVariables(),
      // xLocal(),
            primaryField_local{ 0.0 }
    {}

// #if !defined(CALC_FEM_SHAPE_IN_KERNEL)
//     /// Dummy
//     int xLocal;
// #else
//     /// C-array stack storage for element local the nodal positions.
//     real64 xLocal[ maxNumTestSupportPointsPerElem ][ 3 ];
// #endif

    /// C-array storage for the element local primary field variable.
    real64 primaryField_local[ maxNumTestSupportPointsPerElem ];
  };


  /**
   * @brief Copy global values from primary field to a local stack array.
   * @copydoc geos::finiteElement::ImplicitKernelBase::setup
   *
   * For the LaplaceFEMKernel implementation, global values from the
   * primaryField, and degree of freedom numbers are placed into element local
   * stack storage.
   */
  GEOS_HOST_DEVICE
  inline
  void setup( localIndex const k,
              StackVariables & stack ) const
  {
    m_finiteElementSpace.template setup< FE_TYPE >( k, m_meshData, stack.feStack );
    stack.numRows = m_finiteElementSpace.template numSupportPoints< FE_TYPE >( stack.feStack );
    stack.numCols = stack.numRows;
    for( localIndex a = 0; a < stack.numRows; ++a )
    {
      localIndex const localNodeIndex = m_elemsToNodes( k, a );

      stack.primaryField_local[ a ] = m_primaryField[ localNodeIndex ];
      stack.localRowDofIndex[a] = m_dofNumber[localNodeIndex];
      stack.localColDofIndex[a] = m_dofNumber[localNodeIndex];
    }
    m_finiteElementSpace.template
    addGradGradStabilizationMatrix< FE_TYPE, numDofPerTrialSupportPoint >( stack.feStack,
                                                                           stack.localJacobian );
  }

  /**
   * @copydoc geos::finiteElement::ImplicitKernelBase::quadraturePointKernel
   */
  GEOS_HOST_DEVICE
  inline
  void quadraturePointKernel( localIndex const k, //CellIndexType const k,
                              localIndex const q,
                              StackVariables & stack ) const
  {
    using namespace PDEUtilities;

    constexpr PDEUtilities::FunctionSpace TrialSpace = FE_TYPE::template getFunctionSpace< numDofPerTrialSupportPoint >();
    constexpr PDEUtilities::FunctionSpace TestSpace = TrialSpace;

    // ... Get info for cell k
    CellType cell = m_subregionMesh.getCell( k );

    // ... Get coordinates and weight for quadrature point q
    QuadratureUtilities::Data quadratureData = QuadratureUtilities::getData< CellType,
                                                                             QuadratureUtilities::Rule::Gauss,
                                                                             numQuadraturePointsPerElem >( q );

    // ... Evaluate Jacobian determinant and Jacobian inverse
    auto [ detJ, Jinv ] = CellUtilities::getJacobianDeterminantAndJacobianInverse( cell, quadratureData.Xiq );

    // ... Evaluate basis function gradients
    real64 dNdX[maxNumTestSupportPointsPerElem][3]{{}};
    BasisFunctionUtilities::getGradient< CellType,
                                         BasisFunctionUtilities::BasisFunction::Lagrange,
                                         maxNumTestSupportPointsPerElem >( quadratureData.Xiq, Jinv, dNdX );

    // ... Compute local stiffness matrix
    real64 const detJxW = detJ * quadratureData.wq;
    BilinearFormUtilities::compute< TestSpace,
                                    TrialSpace,
                                    DifferentialOperator::Gradient,
                                    DifferentialOperator::Gradient >
    (
      stack.localJacobian,
      dNdX,
      1.0,
      dNdX,
      detJxW );
  }

  /**
   * @copydoc geos::finiteElement::ImplicitKernelBase::complete
   *
   * Form element residual from the fully formed element Jacobian dotted with
   * the primary field and map the element local Jacobian/Residual to the
   * global matrix/vector.
   */
  GEOS_HOST_DEVICE
  inline
  real64 complete( localIndex const k,
                   StackVariables & stack ) const
  {
    GEOS_UNUSED_VAR( k );
    real64 maxForce = 0;

    for( localIndex a = 0; a < stack.numRows; ++a )
    {
      for( localIndex b = 0; b < stack.numCols; ++b )
      {
        stack.localResidual[ a ] += stack.localJacobian[ a ][ b ] * stack.primaryField_local[ b ];
      }
    }

    for( int a = 0; a < stack.numRows; ++a )
    {
      localIndex const dof = LvArray::integerConversion< localIndex >( stack.localRowDofIndex[ a ] - m_dofRankOffset );
      if( dof < 0 || dof >= m_matrix.numRows() ) continue;
      m_matrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( dof,
                                                                              stack.localColDofIndex,
                                                                              stack.localJacobian[ a ],
                                                                              stack.numCols );

      RAJA::atomicAdd< parallelDeviceAtomic >( &m_rhs[ dof ], stack.localResidual[ a ] );
      maxForce = fmax( maxForce, fabs( stack.localResidual[ a ] ) );
    }

    return maxForce;
  }

protected:
  /// The global primary field array.
  arrayView1d< real64 const > const m_primaryField;

  // using subregionMeshType = isoparametricMesh< traits::ViewTypeConst< typename SUBREGION_TYPE::NodeMapType::base_type > >;
  using SubregionMeshType = typename NumVertexToSubregionMesh< traits::ViewTypeConst< typename SUBREGION_TYPE::NodeMapType::base_type >,
                                                               FE_TYPE::numNodes >::type;
  SubregionMeshType m_subregionMesh;

  using CellIndexType = typename SubregionMeshType::CellIndexType;
  using CellType = typename SubregionMeshType::CellType;
};

/// The factory used to construct a LaplaceFEMKernel.
using LaplaceFEMKernelFactory = finiteElement::KernelFactory< LaplaceFEMKernel,
                                                              arrayView1d< globalIndex const > const,
                                                              globalIndex const,
                                                              CRSMatrixView< real64, globalIndex const > const,
                                                              arrayView1d< real64 > const,
                                                              real64 const,
                                                              string const >;

} // namespace geos

#include "finiteElement/kernelInterface/SparsityKernelBase.hpp"

#endif // GEOS_PHYSICSSOLVERS_SIMPLEPDE_LAPLACEFEMKERNELSNEW_HPP_
