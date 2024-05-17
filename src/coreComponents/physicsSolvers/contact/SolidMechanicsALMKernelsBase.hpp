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
 * @file SolidMechanicsALMKernelsBase.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_CONTACT_SOLIDMECHANICSALMKERNELSBASE_HPP_
#define GEOS_PHYSICSSOLVERS_CONTACT_SOLIDMECHANICSALMKERNELSBASE_HPP_

#include "finiteElement/kernelInterface/InterfaceKernelBase.hpp"
#include "SolidMechanicsALMKernelsHelper.hpp"

namespace geos
{

namespace solidMechanicsALMKernels
{

template< typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
class ALMKernelsBase :
  public finiteElement::InterfaceKernelBase< CONSTITUTIVE_TYPE,
                                             FE_TYPE,
                                             3, 3 >
{
public:
  /// Alias for the base class;
  using Base = finiteElement::InterfaceKernelBase< CONSTITUTIVE_TYPE,
                                                   FE_TYPE, 
                                                   3, 3 >;

  static constexpr int numNodesPerElem = Base::maxNumTestSupportPointsPerElem;

  static constexpr int numQuadraturePointsPerElem = FE_TYPE::numQuadraturePoints;

  using Base::m_dofNumber;
  using Base::m_dofRankOffset;
  using Base::m_finiteElementSpace;
  using Base::m_matrix;
  using Base::m_rhs;

  ALMKernelsBase( NodeManager const & nodeManager,
                  EdgeManager const & edgeManager,
                  FaceManager const & faceManager,
                  localIndex const targetRegionIndex,
                  FaceElementSubRegion & elementSubRegion,
                  FE_TYPE const & finiteElementSpace,
                  CONSTITUTIVE_TYPE & inputConstitutiveType,
                  arrayView1d< globalIndex const > const uDofNumber,
                  arrayView1d< globalIndex const > const bDofNumber,
                  globalIndex const rankOffset,
                  CRSMatrixView< real64, globalIndex const > const inputMatrix,
                  arrayView1d< real64 > const inputRhs,
                  real64 const inputDt, 
                  arrayView1d< localIndex const > const & faceElementList ):
    Base( nodeManager,
          edgeManager,
          faceManager,
          targetRegionIndex,
          elementSubRegion,
          finiteElementSpace,
          inputConstitutiveType,
          uDofNumber,
          rankOffset,
          inputMatrix,
          inputRhs,
          inputDt ),
    m_X( nodeManager.referencePosition()),
    m_faceToNodes(faceManager.nodeList().toViewConst()),
    m_elemsToFaces(elementSubRegion.faceList().toViewConst()),
    m_faceElementList(faceElementList),
    m_bDofNumber( bDofNumber ),
    m_rotationMatrix(elementSubRegion.getField< fields::contact::rotationMatrix >().toViewConst()),
    m_traction(elementSubRegion.getField< fields::contact::traction >().toViewConst()),
    m_dispJump( elementSubRegion.getField< fields::contact::dispJump >().toView() ),
    m_oldDispJump( elementSubRegion.getField< fields::contact::oldDispJump >().toViewConst() ),
    m_penalty( elementSubRegion.getField< fields::contact::penalty >().toViewConst() )
{}

  struct StackVariables  
  {
    /// The number of displacement dofs per element.
    static constexpr int numUdofs = numNodesPerElem * 3 * 2;

    /// The number of lagrange multiplier dofs per element.
    static constexpr int numTdofs = 3;

    /// The number of bubble dofs per element.
    static constexpr int numBdofs = 3*2;

  public:
  
    GEOS_HOST_DEVICE
    StackVariables():
      dispEqnRowIndices{},
      dispColIndices{},
      bEqnRowIndices{},
      bColIndices{},
      localRu{},
      localRb{},
      localAutAtu{{}},
      localAbtAtb{{}},
      localAbtAtu{{}},
      localAutAtb{{}},
      localAtu{{}},
      localAtb{{}},
      localRotationMatrix{{}},
      localPenalty{{}},
      tLocal{},
      dispJumpLocal{},
      oldDispJumpLocal{},
      X{{}}
    {}

    /// C-array storage for the element local row degrees of freedom.
    globalIndex dispEqnRowIndices[numUdofs];

    /// C-array storage for the element local column degrees of freedom.
    globalIndex dispColIndices[numUdofs];

    /// C-array storage for the element local row degrees of freedom.
    globalIndex bEqnRowIndices[numBdofs];

    /// C-array storage for the element local column degrees of freedom.
    globalIndex bColIndices[numBdofs];

    /// C-array storage for the element local Ru residual vector.
    real64 localRu[numUdofs];

    /// C-array storage for the element local Rb residual vector.
    real64 localRb[numBdofs];

    /// C-array storage for the element local AutAtu matrix.
    real64 localAutAtu[numUdofs][numUdofs];

    /// C-array storage for the element local AbtAtb matrix.
    real64 localAbtAtb[numBdofs][numBdofs];

    /// C-array storage for the element local AbtAtu matrix.
    real64 localAbtAtu[numBdofs][numUdofs];

    /// C-array storage for the element local AbtAtu matrix.
    real64 localAutAtb[numUdofs][numBdofs];

    /// C-array storage for the element local Atu matrix.
    real64 localAtu[numTdofs][numUdofs];

    /// C-array storage for the element local Atb matrix.
    real64 localAtb[numTdofs][numBdofs];

    /// C-array storage for rotation matrix
    real64 localRotationMatrix[3][3];

    /// C-array storage for penalty matrix
    real64 localPenalty[3][3];

    /// Stack storage for the element local lagrange multiplier vector
    real64 tLocal[numTdofs];

    /// Stack storage for the element local displacement jump vector
    real64 dispJumpLocal[numTdofs];

    /// Stack storage for the element local old displacement jump vector
    real64 oldDispJumpLocal[numTdofs];

    /// local nodal coordinates
    real64 X[ numNodesPerElem ][ 3 ];

  };

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

    forAll< POLICY >( kernelComponent.m_faceElementList.size(),
                      [=] GEOS_HOST_DEVICE ( localIndex const i )
    {
      //std::cout << "# QuadPoints: " << numQuadraturePointsPerElem << std::endl;
      //std::cout << "# NodesperElem: " << numNodesPerElem << std::endl;

      localIndex k = kernelComponent.m_faceElementList[i];
      typename KERNEL_TYPE::StackVariables stack;

      //std::cout << i << " " << k << std::endl;

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

  GEOS_HOST_DEVICE
  inline
  void quadraturePointKernel( localIndex const k,
                              localIndex const q,
                              StackVariables & stack ) const
  {
    GEOS_UNUSED_VAR( k );
    real64 const detJ = m_finiteElementSpace.template transformedQuadratureWeight( q, stack.X );

    //std::cout << "detJ: " << detJ << std::endl;

    real64 N[ numNodesPerElem ]; 
    m_finiteElementSpace.template calcN( q, N );

    real64 BubbleN[1]; 
    // It is needed only because I inserted a placeholder for calcBubbleN in some finite elements
    BubbleN[0]=0.0;  //make 0
    constexpr int bperm[1] = {0};
    m_finiteElementSpace.template calcBubbleN( q, BubbleN );

    solidMechanicsALMKernelsHelper::accumulateAtuLocalOperator<stack.numTdofs, 
                                                               stack.numUdofs, 
                                                               numNodesPerElem>(stack.localAtu, 
                                                                                N, 
                                                                                FE_TYPE::permutation,
                                                                                detJ);

    solidMechanicsALMKernelsHelper::accumulateAtuLocalOperator<stack.numTdofs, 
                                                               stack.numBdofs, 
                                                               1>(stack.localAtb, 
                                                                  BubbleN, 
                                                                  bperm,
                                                                  detJ);
  }

protected:

  /// The array containing the nodal position array.
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const m_X;

  ArrayOfArraysView< localIndex const > const m_faceToNodes; 

  ArrayOfArraysView< localIndex const > const m_elemsToFaces;

  arrayView1d< localIndex const > const m_faceElementList;

  /// The global degree of freedom number of bubble
  arrayView1d< globalIndex const > const m_bDofNumber;

  arrayView3d< real64 const > const m_rotationMatrix;

  arrayView2d< real64 const > const m_traction;

  arrayView2d< real64 > const m_dispJump;

  arrayView2d< real64 const > const m_oldDispJump;

  arrayView2d< real64 const > const m_penalty;

};

/*
using ALMFactory = finiteElement::InterfaceKernelFactory< ALMKernelsBase,
                                                          arrayView1d< globalIndex const > const,
                                                          globalIndex const,
                                                          CRSMatrixView< real64, globalIndex const > const,
                                                          arrayView1d< real64 > const,
                                                          real64 const, 
                                                          arrayView1d< localIndex const > const >;


struct ComputeRotationMatricesKernel
{
  template< typename POLICY >
  static void
  launch( localIndex const size,
          arrayView2d< real64 const > const & faceNormal,
          ArrayOfArraysView< localIndex const > const & elemsToFaces,
          arrayView3d< real64 > const & rotationMatrix )
  {

    forAll< POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const k )
    {

      localIndex const & f0 = elemsToFaces[k][0];
      localIndex const & f1 = elemsToFaces[k][1];

      //stackArray1d< real64, 3 > Nbar( 3 );
      real64 Nbar[3];
      Nbar[0] = faceNormal[f0][0] - faceNormal[f1][0];
      Nbar[1] = faceNormal[f0][1] - faceNormal[f1][1];
      Nbar[2] = faceNormal[f0][2] - faceNormal[f1][2];

      LvArray::tensorOps::normalize< 3 >( Nbar );
      computationalGeometry::RotationMatrix_3D( Nbar, rotationMatrix[k] );

    } );
  }

};
*/

} // namespace SolidMechanicsALMKernels

} // namespace geos


#endif /* GEOS_PHYSICSSOLVERS_CONTACT_SOLIDMECHANICSALMKERNELSBASE_HPP_ */