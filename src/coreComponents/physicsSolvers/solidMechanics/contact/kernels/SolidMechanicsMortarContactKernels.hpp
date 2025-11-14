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
 * @file SolidMechanicsALMKernelsBase.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_CONTACT_KERNELS_SOLIDMECHANICSMORTARCONTACTKERNELS_HPP_
#define GEOS_PHYSICSSOLVERS_CONTACT_KERNELS_SOLIDMECHANICSMORTARCONTACTKERNELS_HPP_

#include "finiteElement/kernelInterface/InterfaceKernelBase.hpp"
#include "codingUtilities/Utilities.hpp"

namespace geos
{

namespace solidMechanicsMortarContactKernels
{

/**
 * @brief Implements kernels for ALM.
 * @copydoc geos::finiteElement::InterfaceKernelBase
 *
 */
template< typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
class MortarContactKernels :
  public finiteElement::InterfaceKernelBase< CONSTITUTIVE_TYPE,
                                             FE_TYPE,
                                             3, 3 >
{
public:
  /// Alias for the base class;
  using Base = finiteElement::InterfaceKernelBase< CONSTITUTIVE_TYPE,
                                                   FE_TYPE,
                                                   3, 3 >;

  /// Alias for the sub triangle element type
  using feTriangleCell = finiteElement::H1_TriangleFace_Lagrange1_Gauss4;

  /// Number of nodes per element...which is equal to the
  /// numTestSupportPointPerElem and numTrialSupportPointPerElem by definition.
  static constexpr int numNodesPerElem = Base::maxNumTestSupportPointsPerElem;

  /// Compile time value for the number of quadrature points per sub triangle.
  static constexpr int numQuadraturePointsPerElem = feTriangleCell::numQuadraturePoints;

  /// The number of displacement dofs per element.
  static constexpr int numUdofs = numNodesPerElem * 3;

  /// The number of lagrange multiplier dofs per element.
  static constexpr int numTdofs = 3;

  using Base::m_dofNumber;
  using Base::m_dofRankOffset;
  using Base::m_finiteElementSpace;
  using Base::m_matrix;
  using Base::m_rhs;

  /**
   * @brief Constructor
   * @copydoc geos::finiteElement::InterfaceKernelBase::InterfaceKernelBase
   */
  MortarContactKernels( NodeManager const & nodeManager,
                        EdgeManager const & edgeManager,
                        FaceManager const & faceManager,
                        localIndex const targetRegionIndex,
                        FaceElementSubRegion & subRegion,
                        FE_TYPE const & finiteElementSpace,
                        CONSTITUTIVE_TYPE & inputConstitutiveType,
                        arrayView1d< globalIndex const > const uDofNumber,
                        globalIndex const rankOffset,
                        CRSMatrixView< real64, globalIndex const > const inputMatrix,
                        arrayView1d< real64 > const inputRhs,
                        real64 const inputDt,
                        FaceElementSubRegion const & multipliersRegion,
                        FaceElementSubRegion const & primaryVarRegion,
                        arrayView1d< localIndex const > const faceElementList1,
                        arrayView1d< localIndex const > const faceElementList2,
                        arrayView1d< real64 const > const subTriangleDeterminants,
                        arrayView3d< real64 const > const gpLocalCoords,
                        string const tractionDofKey ):
    Base( nodeManager,
          edgeManager,
          faceManager,
          targetRegionIndex,
          subRegion,
          finiteElementSpace,
          inputConstitutiveType,
          uDofNumber,
          rankOffset,
          inputMatrix,
          inputRhs,
          inputDt ),
    //m_X( nodeManager.referencePosition()),
    m_faceToNodes( faceManager.nodeList().toViewConst()),
    m_elemsToFaces( primaryVarRegion.faceList().toViewConst()),
    m_faceElementList1( faceElementList1 ),
    m_faceElementList2( faceElementList2 ),
    m_subTriangleDeterminants( subTriangleDeterminants ),
    m_gpLocalCoords( gpLocalCoords ),
    m_traction( multipliersRegion.getField< fields::contact::traction >().toViewConst() ),
    m_tDofNumber( multipliersRegion.getReference< globalIndex_array >( tractionDofKey ).toViewConst() ),
    m_displacement( nodeManager.getField< fields::solidMechanics::totalDisplacement >() ),
    m_rotationMatrix( multipliersRegion.getField< fields::contact::rotationMatrix >().toViewConst())
  {}

  //***************************************************************************
  /**
   * @copydoc finiteElement::InterfaceKernelBase::StackVariables
   */
  struct StackVariables
  {

public:

    GEOS_HOST_DEVICE
    StackVariables():
      dispEqnRowIndices{},
      dispColIndices{},
      tEqnRowIndices{},
      tColIndices{},
      localRu{},
      localRt{},
      localAtu{ {} },
      localAut{ {} },
      uLocal{},
      det{},
      localRotationMatrix{ {} }
    {}

    /// C-array storage for the element local row degrees of freedom.
    localIndex dispEqnRowIndices[numUdofs];

    /// C-array storage for the element local column degrees of freedom.
    globalIndex dispColIndices[numUdofs];

    /// C-array storage for the traction local row degrees of freedom.
    localIndex tEqnRowIndices[numTdofs];

    /// C-array storage for the element local column degrees of freedom.
    globalIndex tColIndices[numTdofs];

    /// C-array storage for the element local Ru residual vector.
    real64 localRu[numUdofs];

    /// C-array storage for the element local Rt residual vector.
    real64 localRt[numTdofs];

    /// C-array storage for the element local Atu matrix.
    real64 localAtu[numTdofs][numUdofs];

    /// C-array storage for the element local Aut matrix.
    real64 localAut[numUdofs][numTdofs];

    /// C-array storage for the element local displacement vector
    real64 uLocal[numUdofs];

    /// C-array storage for the subtriangle cell determinants
    real64 det[numQuadraturePointsPerElem];

    /// C-array storage for rotation matrix
    real64 localRotationMatrix[3][3];


  };

  //***************************************************************************

  /**
   * @copydoc ::geos::finiteElement::InterfaceKernelBase::kernelLaunch
   *
   */
  //START_kernelLauncher
  template< typename POLICY,
            typename KERNEL_TYPE >
  static
  real64
  kernelLaunch( localIndex const numElems,
                KERNEL_TYPE const & kernelComponent )
  {
    GEOS_MARK_FUNCTION;
    GEOS_UNUSED_VAR( numElems );

    // Define a RAJA reduction variable to get the maximum residual contribution.
    RAJA::ReduceMax< ReducePolicy< POLICY >, real64 > maxResidual( 0 );

    // Loop over all existing triangular integration subcells
    forAll< POLICY >( kernelComponent.m_faceElementList1.size(),
                      [=] GEOS_HOST_DEVICE ( localIndex const k )
    {
      typename KERNEL_TYPE::StackVariables stack;

      kernelComponent.setup( k, stack );
      for( integer q=0; q<numQuadraturePointsPerElem; ++q )
      {
        kernelComponent.quadraturePointKernel( k, q, stack );
      }
      maxResidual.max( kernelComponent.complete( k, stack ) );
    } );

    return maxResidual.get();
  }
  //END_kernelLauncher


  /**
   * @brief Copy global values from primary field to a local stack array.
   * @copydoc ::geos::finiteElement::InterfaceKernelBase::setup
   */
  GEOS_HOST_DEVICE
  inline
  void setup( localIndex const k,
              StackVariables & stack ) const
  {

    localIndex const kt = m_faceElementList1[k];
    localIndex const ku = m_faceElementList2[k];

    localIndex const kfu = m_elemsToFaces[ku][0];

    for( localIndex a=0; a<numNodesPerElem; ++a )
    {
      localIndex const knu = m_faceToNodes( kfu, a );

      for( int i=0; i<3; ++i )
      {
        stack.dispEqnRowIndices[a*3+i] = m_dofNumber[knu]+i-m_dofRankOffset;
        stack.dispColIndices[a*3+i] = m_dofNumber[knu]+i;
        stack.uLocal[a*3+i] = m_displacement[knu][i];
      }
    }

    for( int j=0; j<3; ++j )
    {
      for( int i=0; i<3; ++i )
      {
        stack.localRotationMatrix[ i ][ j ] = m_rotationMatrix( kt, i, j );
      }
    }

    for( int i=0; i<3; ++i )
    {
      stack.tEqnRowIndices[i] = m_tDofNumber[kt] + i - m_dofRankOffset;
      stack.tColIndices[i] = m_tDofNumber[kt] + i;
    }

    for( int i=0; i<numQuadraturePointsPerElem; ++i )
    {
      stack.det[i] = m_subTriangleDeterminants[k]*feTriangleCell::getQuadratureWeight( i );
    }

  }


  GEOS_HOST_DEVICE
  inline
  void quadraturePointKernel( localIndex const k,
                              localIndex const q,
                              StackVariables & stack ) const
  {

    int permutation[numNodesPerElem];
    m_finiteElementSpace.getPermutation( permutation );

    //constexpr int nUdof = numNodesPerElem*3;

    real64 gpLocalCoords[2];
    gpLocalCoords[0] = m_gpLocalCoords( k, q, 0 );
    gpLocalCoords[1] = m_gpLocalCoords( k, q, 1 );

    real64 Nu[ numNodesPerElem ];
    m_finiteElementSpace.calcN( gpLocalCoords, Nu );

    // accumulate local stack matrix
    for( int a=0; a < numNodesPerElem; ++a )
    {
      for( int i=0; i < 3; ++i )
      {
        stack.localAtu[i][ a*3 + i ] += Nu[ permutation[ a ] ] * stack.det[q];
      }
    }

  }


  GEOS_HOST_DEVICE
  inline
  real64 complete( localIndex const k,
                   StackVariables & stack ) const
  {
    localIndex ke = m_faceElementList1[k];

    real64 rhsU[numUdofs];
    real64 rhsT[numTdofs];

    real64 matRRtAtu[3][numUdofs];

    // transp(R) * Atu
    LvArray::tensorOps::Rij_eq_AkiBkj< 3, numUdofs, 3 >( matRRtAtu, stack.localRotationMatrix, stack.localAtu );
    LvArray::tensorOps::copy< numTdofs, numUdofs >( stack.localAtu, matRRtAtu );
    LvArray::tensorOps::transpose< numUdofs, numTdofs >( stack.localAut, stack.localAtu );

    // Compute the traction contribute of the local residuals
    LvArray::tensorOps::Ri_eq_AijBj< numUdofs, numTdofs >( rhsU, stack.localAut, m_traction[ke] );
    LvArray::tensorOps::scaledAdd< numUdofs >( stack.localRu, rhsU, 1.0 );

    // Compute the displacemement contribute of the local residual
    LvArray::tensorOps::Ri_eq_AijBj< numTdofs, numUdofs >( rhsT, stack.localAtu, stack.uLocal );
    LvArray::tensorOps::scaledAdd< numTdofs >( stack.localRt, rhsT, 1.0 );

    fillGlobalMatrix( stack );

    return 0.0;
  }


  /**
   * @brief Fill global matrix and residual vector
   *
   * @param stack stack variables
   */
  GEOS_HOST_DEVICE
  void fillGlobalMatrix( StackVariables & stack ) const
  {


    for( localIndex i=0; i < numTdofs; ++i )
    {
      localIndex const dof = LvArray::integerConversion< localIndex >( stack.tEqnRowIndices[ i ] );

      if( dof < 0 || dof >= m_matrix.numRows() ) continue;

      // TODO: May not need to be an atomic operation
      RAJA::atomicAdd< parallelDeviceAtomic >( &m_rhs[dof], stack.localRt[i] );

      // Fill in matrix block Atu
      m_matrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( dof,
                                                                              stack.dispColIndices,
                                                                              stack.localAtu[i],
                                                                              numUdofs );

    }

    for( localIndex i=0; i < numUdofs; ++i )
    {
      localIndex const dof = LvArray::integerConversion< localIndex >( stack.dispEqnRowIndices[ i ] );

      if( dof < 0 || dof >= m_matrix.numRows() ) continue;

      // Is it necessary? Each row should be indepenedent
      RAJA::atomicAdd< parallelDeviceAtomic >( &m_rhs[dof], stack.localRu[i] );

      // Fill in matrix
      m_matrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( dof,
                                                                              stack.tColIndices,
                                                                              stack.localAut[i],
                                                                              numTdofs );

    }

  }


protected:

  /// The array of array containing the face to node map.
  ArrayOfArraysView< localIndex const > const m_faceToNodes;

  /// The array of array containing the element to face map.
  arrayView2d< localIndex const > const m_elemsToFaces;

  /// The array containing the list of face element of the same type for the traction side
  /// for each mortar subtriangle.
  arrayView1d< localIndex const > const m_faceElementList1;

  /// The array containing the list of face element of the same type for the displacement side
  /// for each mortar subtriangle.
  arrayView1d< localIndex const > const m_faceElementList2;

  /// The array containing the determinants of the subtriangles for mortar integration.
  arrayView1d< real64 const > const m_subTriangleDeterminants;

  /// The array containing the local gauss point coordinates projected on the displacement side
  arrayView3d< real64 const > const m_gpLocalCoords;

  /// The array containing the traction field
  arrayView2d< real64 const > const m_traction;

  /// The array holding the traction degrees of freedom
  arrayView1d< globalIndex const > const m_tDofNumber;

  /// The array containing the displacement field
  arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const m_displacement;


  /// The array containing the rotation matrix for each element.
  arrayView3d< real64 const > const m_rotationMatrix;

};

/// The factory used to construct a Mortar kernel.
using MortarFactory = finiteElement::InterfaceKernelFactory< MortarContactKernels,
                                                             arrayView1d< globalIndex const > const,
                                                             globalIndex const,
                                                             CRSMatrixView< real64, globalIndex const > const,
                                                             arrayView1d< real64 > const,
                                                             real64 const,
                                                             FaceElementSubRegion const &,
                                                             FaceElementSubRegion const &,
                                                             arrayView1d< localIndex const > const &,
                                                             arrayView1d< localIndex const > const &,
                                                             arrayView1d< real64 const > const &,
                                                             arrayView3d< real64 const > const &,
                                                             string const >;


} // namespace solidMechanicsMortarContactKernels

} // namespace geos


#endif /* GEOS_PHYSICSSOLVERS_CONTACT_KERNELS_SOLIDMECHANICSMORTARCONTACTKERNELS_HPP_ */
