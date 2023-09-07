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
 * @file LaplaceFEMKernelsNew.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_SIMPLEPDE_LAPLACEFEMKERNELSNEW_HPP_
#define GEOS_PHYSICSSOLVERS_SIMPLEPDE_LAPLACEFEMKERNELSNEW_HPP_

#define GEOSX_DISPATCH_VEM /// enables VEM in FiniteElementDispatch

#include "finiteElement/kernelInterface/ImplicitKernelBase.hpp"
#include "BasisFunctionUtilities.hpp"
#include "IsoparametricMeshUtilities.hpp"
#include "QuadratureUtilities.hpp"

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
 * ### LaplaceFEMKernelNew Description
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
class LaplaceFEMKernelNew :
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
  LaplaceFEMKernelNew( NodeManager const & nodeManager,
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
          inputRhs ),
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
            xLocal(),
            primaryField_local{ 0.0 }
    {}

#if !defined(CALC_FEM_SHAPE_IN_KERNEL)
    /// Dummy
    int xLocal;
#else
    /// C-array stack storage for element local the nodal positions.
    real64 xLocal[ maxNumTestSupportPointsPerElem ][ 3 ];
#endif

    /// C-array storage for the element local primary field variable.
    real64 primaryField_local[ maxNumTestSupportPointsPerElem ];
  };


  /**
   * @brief Copy global values from primary field to a local stack array.
   * @copydoc geos::finiteElement::ImplicitKernelBase::setup
   *
   * For the LaplaceFEMKernelNew implementation, global values from the
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

// #if defined(CALC_FEM_SHAPE_IN_KERNEL)
//       for( int i=0; i<3; ++i )
//       {
//         stack.xLocal[ a ][ i ] = m_X[ localNodeIndex ][ i ];
//       }
// #endif

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
  void quadraturePointKernel( localIndex const k,
                              localIndex const q,
                              StackVariables & stack ) const
  {
    CellType cell = m_subregionMesh.getCell( k );
#if defined(CALC_FEM_SHAPE_IN_KERNEL)
      cell.getLocalCoordinates( stack.xLocal );
#endif

    // Temporary implementation for isoparametric hexahedron
    
    // ... Get parent element coordinates for quadrature point q 
    //     (quadrature rule responsability)

    // real64 const val = 0.5773502691896257645092;
    // int const a = q & 1;
    // int const b = ( q & 2 ) >> 1;
    // int const c = ( q & 4 ) >> 2;
    // real64 const Xiq[3] = { ( 2 * a - 1 ) * val,
    //                         ( 2 * b - 1 ) * val,
    //                         ( 2 * c - 1 ) * val };
    // real64 const wq = 1.0; // weight

    QuadratureUtilities::Data quadratureData = QuadratureUtilities::getData< CellType, 
                                                                             QuadratureUtilities::Rule::Gauss,
                                                                             numQuadraturePointsPerElem >( q );

    real64 const wq = quadratureData.wq;
    real64 const Xiq[3] = { quadratureData.Xiq[0],
                            quadratureData.Xiq[1],
                            quadratureData.Xiq[2] };

    // real64 & wq = quadratureData.wq;
    // real64 (& Xiq)[3] = quadratureData.Xiq;

    // [Xiq, wq] = QuadratureRule::getQuadraturePointCoordsAndWeight( q )
    // getQuadratureData< CELL_TYPE, INTEGRATION_RULE, INTEGRATION_ORDER >
    // {
    //   return Helper< CELL_TYPE, INTEGRATION_RULE, INTEGRATION_ORDER >::getQuadrat... ;
    // }

    //   Helper get...
    //   {
    //     return std::make_pair( Xiq, wq );
    //   }
    // ... 
    // const auto [ Xiq, wq ] = getQuadratureData< CellType, IntegrationRule::Gauss, 1>( q ) 

    // onst auto [ Xiq, wq ] = getQuadratureData< HexahedronCell, IntegrationRule::Gauss, 1>()
    // {

    // } 
    // const auto [ Xiq, wq ] = getQuadratureData< CellType, IntegrationRule::Gauss, 1>() 

    // ... Evaluate Jacobian
    typename CellType::JacobianType J = cell.getJacobian( quadratureData.Xiq );
    real64 const detJxW = quadratureData.wq * LvArray::tensorOps::invert< 3 >( J.val );
    typename CellType::JacobianType const & Jinv = J;
    
     // ... Evaluate basis function gradients
    real64 dNdX[maxNumTestSupportPointsPerElem][3]{{}};
    for( int i = 0; i < maxNumTestSupportPointsPerElem; ++i )
    {
      // Basis function in parent space
      BasisFunctionUtilities::Gradient dNdXi = BasisFunctionUtilities::getGradient< CellType,
                                                                                    BasisFunctionUtilities::BasisFunction::Lagrange >( i, quadratureData.Xiq );

      // Basis function in physical space
      dNdX[i][0] = dNdXi.val[0] * Jinv.val[0][0] + dNdXi.val[1] * Jinv.val[1][0] + dNdXi.val[2] * Jinv.val[2][0];
      dNdX[i][1] = dNdXi.val[0] * Jinv.val[0][1] + dNdXi.val[1] * Jinv.val[1][1] + dNdXi.val[2] * Jinv.val[2][1];
      dNdX[i][2] = dNdXi.val[0] * Jinv.val[0][2] + dNdXi.val[1] * Jinv.val[1][2] + dNdXi.val[2] * Jinv.val[2][2];
    }

    // dNdX = m_finiteElementSpace.template getGradShapeFunctions< FE_TYPE >( Xiq, Jinv );

    // ... Compute local stiffness matrix
    for( localIndex i = 0; i < stack.numRows; ++i )
    {
      for( localIndex j = 0; j < stack.numCols; ++j )
      {
        stack.localJacobian[ i ][ j ] += LvArray::tensorOps::AiBi< 3 >( dNdX[i], dNdX[j] ) * detJxW;
      }
    }

    // real64 dNdX[ maxNumTestSupportPointsPerElem ][ 3 ];
    // real64 const detJxW = m_finiteElementSpace.template getGradN< FE_TYPE >( k, q, stack.xLocal,
                                                                          //  stack.feStack, dNdX );
    // for( localIndex a = 0; a < stack.numRows; ++a )
    // {
      // for( localIndex b = 0; b < stack.numCols; ++b )
      // {
        // stack.localJacobian[ a ][ b ] += LvArray::tensorOps::AiBi< 3 >( dNdX[a], dNdX[b] ) * detJxW;
      // }
    // }
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
  /// The array containing the nodal position array.
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const m_X;

  /// The global primary field array.
  arrayView1d< real64 const > const m_primaryField;

  // using subregionMeshType = isoparametricMesh< traits::ViewTypeConst< typename SUBREGION_TYPE::NodeMapType::base_type > >;
  using SubregionMeshType = typename NumVertexToSubregionMesh< traits::ViewTypeConst< typename SUBREGION_TYPE::NodeMapType::base_type >,
                                                               FE_TYPE::numNodes  >::type;
  SubregionMeshType m_subregionMesh;
   
  using CellType = typename SubregionMeshType::CellType;
};

/// The factory used to construct a LaplaceFEMKernelNew.
using LaplaceFEMKernelNewFactory = finiteElement::KernelFactory< LaplaceFEMKernelNew,
                                                              arrayView1d< globalIndex const > const,
                                                              globalIndex const,
                                                              CRSMatrixView< real64, globalIndex const > const,
                                                              arrayView1d< real64 > const,
                                                              string const >;

} // namespace geos

#include "finiteElement/kernelInterface/SparsityKernelBase.hpp"

#endif // GEOS_PHYSICSSOLVERS_SIMPLEPDE_LAPLACEFEMKERNELSNEW_HPP_
