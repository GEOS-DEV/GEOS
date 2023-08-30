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

namespace geos
{

//namespace cells TODO?
struct Dense3x3Tensor
{
  real64 val[3][3]{{}};
};

class HexadronCell
{
public:
  constexpr static int numVertex = 8;

  using JacobianType = Dense3x3Tensor;

  GEOS_HOST_DEVICE
  HexadronCell( real64 nodeCoords[numVertex][3] )
  {
    for( int i = 0; i < numVertex; ++i )
    {
      for( int j = 0; j < 3; ++j )
      {
        m_nodeCoords[i][j] = nodeCoords[i][j];
      }
    }
  }

  GEOS_HOST_DEVICE
  JacobianType getJacobian( real64 refPointCoords[3] ) const
  {
    JacobianType J;
    
    // Compute Jacobian
    real64 dPhiLin[2] = { -1.0, 1.0 };
    for( int k = 0; k < 2; ++k )
    {
      for( int j = 0; j < 2; ++j )
      {
        for( int i = 0; i < 2; ++i )
        {
          real64 gradPhi[3]{ 0.125 * dPhiLin[i] * ( 1.0 - refPointCoords[1]) * ( 1.0 - refPointCoords[2]),
                             0.125 * ( 1.0 - refPointCoords[0]) * dPhiLin[j]  * ( 1.0 - refPointCoords[2]),
                             0.125 * ( 1.0 - refPointCoords[0]) * ( 1.0 - refPointCoords[1]) * dPhiLin[k] };

          int vertexInd = 4 * k + 2 * j + i;

          J.val[0][0] += m_nodeCoords[vertexInd][0] * gradPhi[0]; 
          J.val[0][1] += m_nodeCoords[vertexInd][0] * gradPhi[1]; 
          J.val[0][2] += m_nodeCoords[vertexInd][0] * gradPhi[2]; 
          J.val[1][0] += m_nodeCoords[vertexInd][1] * gradPhi[0]; 
          J.val[1][1] += m_nodeCoords[vertexInd][1] * gradPhi[1]; 
          J.val[1][2] += m_nodeCoords[vertexInd][1] * gradPhi[2]; 
          J.val[2][0] += m_nodeCoords[vertexInd][2] * gradPhi[0]; 
          J.val[2][1] += m_nodeCoords[vertexInd][2] * gradPhi[1]; 
          J.val[2][2] += m_nodeCoords[vertexInd][2] * gradPhi[2]; 
        }
      }
    }

    return J;
  }

  // GEOS_HOST_DEVICE real64[3] mapping( real64 refPointCoords[3] ) const
  // {

  // }

private:
  real64 m_nodeCoords[numVertex][3];
};

// class HexadronIJKCell
// {
// public:
  // constexpr static int numVertex = 8;
// 
  // using JacobianType = real64[3][3];
// 
  // GEOS_HOST_DEVICE HexadronIJKCell( real64[3] h )
    // :
    // m_h( h )
  // {
// 
  // }
// 
  // GEOS_HOST_DEVICE JacobianType getJacobian( real64[3] refPointCoords )
  // {
    // return 0.5 * m_h;
  // }
// 
  // GEOS_HOST_DEVICE real64[3] mapping( real64[3] refPointCoords ) const
  // {
// 
  // }
// 
// private:
  // real64[3] m_h;
// };


template< typename ARRAY_VIEW_TYPE  >
class isoparametricMesh
{
public:
  isoparametricMesh( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X,
                     ARRAY_VIEW_TYPE const elementToNodes )
   : m_X( X ), m_elementToNodes( elementToNodes )
  {

  }
  
  localIndex numCells() const
  {
    return m_elementToNodes.size( 0 );
  }

protected:
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const m_X;
  // traits::ViewTypeConst< typename SUBREGION_TYPE::NodeMapType::base_type > const m_elementToNodes;
  ARRAY_VIEW_TYPE const m_elementToNodes;
};

template< typename ARRAY_VIEW_TYPE  >
class isoparametricHexahedronMesh : public isoparametricMesh< ARRAY_VIEW_TYPE >
{
public :
  // using const_iterator = CellIterator< HexadronCell >;
  //...

 using CellType = HexadronCell;
  // ... 
  constexpr static int numCellVertex = HexadronCell::numVertex;

  // ...
  // constructor

  // 
  HexadronCell getCell( localIndex k ) const
  {
    return HexadronCell( this->m_X[k] );
  }

  // localIndex numCells() const;
};

// class ijkMesh
// {
  // ...
  // constexpr static nVertex = 8;
  // ...
  // constructor
// 
  
  // HexadronIJKCell getCell( localIndex k )
  // {
    // return ...;
  // }
// private:
  // ArrayOfArrays< real64 const > const m_h;
// 
// };

// class isoparametricWedgeMesh : public isoparametricMesh
// {
//   // ... 
//   constexpr static nVertex = 6;
//   // ...
//   WedgeCell getCell( localIndex k )
//   {
//     return ...;
//   }
// };

// Factory selectign dynamically the appropriate isopametric mesh
template< typename ARRAY_VIEW_TYPE  >
isoparametricMesh< ARRAY_VIEW_TYPE > selectIsoparametricMesh( ElementType elemType,
                                           arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const nodePositions,
                                           ARRAY_VIEW_TYPE const elementToNodes )
{
  if( elemType == ElementType::Hexahedron)
  {
    return isoparametricHexahedronMesh< ARRAY_VIEW_TYPE >( nodePositions, elementToNodes );
  }
  // else if( elemType == ElementType::Wedge) )
  // {
    // return isoparametricWedgenMesh( nodePositions, elementToNodes );
  // }
  else
  {
    GEOS_ERROR( "finiteElement::dispatchlowOrder3D() is not implemented for input of " << elemType );
  }
}


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
    m_subregionMesh( selectIsoparametricMesh< traits::ViewTypeConst< typename SUBREGION_TYPE::NodeMapType::base_type > >(
                       elementSubRegion.getElementType(),
                       nodeManager.referencePosition(),
                       elementSubRegion.nodeList().toViewConst()) ),
    m_primaryField( nodeManager.template getReference< array1d< real64 > >( fieldName ))
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
//    m_subregionMesh::CellType c = m_subregionMesh.getCell( k );
    m_finiteElementSpace.template setup< FE_TYPE >( k, m_meshData, stack.feStack );
    stack.numRows = m_finiteElementSpace.template numSupportPoints< FE_TYPE >( stack.feStack );
    stack.numCols = stack.numRows;
    for( localIndex a = 0; a < stack.numRows; ++a )
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
    typename subregionMeshType::CellType c = m_subregionMesh.getCell( k );
    // real64[3] pts = quadRule::getCoords< RULE >( q ); // ???
    // m_subregionMesh::JacobianType j = c.getJacobian( pts );
    // real64 const detJ = LvArray::tensorOps::determinant<...>( j );


    real64 dNdX[ maxNumTestSupportPointsPerElem ][ 3 ];
    real64 const detJ = m_finiteElementSpace.template getGradN< FE_TYPE >( k, q, stack.xLocal,
                                                                           stack.feStack, dNdX );
    for( localIndex a = 0; a < stack.numRows; ++a )
    {
      for( localIndex b = 0; b < stack.numCols; ++b )
      {
        stack.localJacobian[ a ][ b ] += LvArray::tensorOps::AiBi< 3 >( dNdX[a], dNdX[b] ) * detJ;
      }
    }
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

  using subregionMeshType = isoparametricMesh< traits::ViewTypeConst< typename SUBREGION_TYPE::NodeMapType::base_type > >;
  subregionMeshType m_subregionMesh;

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
