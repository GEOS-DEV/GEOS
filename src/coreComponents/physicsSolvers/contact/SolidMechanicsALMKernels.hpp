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
 * @file SolidMechanicsALMKernel.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_CONTACT_SOLIDMECHANICSALMKERNELS_HPP_
#define GEOS_PHYSICSSOLVERS_CONTACT_SOLIDMECHANICSALMKERNELS_HPP_

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

  ALMKernelsBase( NodeManager const & nodeManager,
                  EdgeManager const & edgeManager,
                  FaceManager const & faceManager,
                  localIndex const targetRegionIndex,
                  FaceElementSubRegion const & elementSubRegion,
                  FE_TYPE const & finiteElementSpace,
                  CONSTITUTIVE_TYPE & inputConstitutiveType,
                  arrayView1d< globalIndex const > const inputDofNumber,
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
          inputDofNumber,
          rankOffset,
          inputMatrix,
          inputRhs,
          inputDt ),
    m_X( nodeManager.referencePosition()),
    m_faceToNodes(faceManager.nodeList().toViewConst()),
    m_elemsToFaces(elementSubRegion.faceList().toViewConst()),
    m_faceElementList(faceElementList),
    m_rotationMatrix(elementSubRegion.getField< fields::contact::rotationMatrix >().toViewConst()),
    m_deltaTraction(elementSubRegion.getField< fields::contact::deltaTraction >().toViewConst())
{}

  struct StackVariables  
  {
    /// The number of displacement dofs per element.
    static constexpr int numUdofs = numNodesPerElem * 3 * 2;

    /// The number of lagrange multiplier dofs per element.
    static constexpr int numTdofs = 3;

  public:
  
    GEOS_HOST_DEVICE
    StackVariables():
      dispEqnRowIndices{},
      dispColIndices{},
      localRu{},
      localAutAtu{{}},
      localAtu{},
      tLocal{},
      localRotationMatrix{{}},
      X{}
    {}

    /// C-array storage for the element local row degrees of freedom.
    globalIndex dispEqnRowIndices[numUdofs];

    /// C-array storage for the element local column degrees of freedom.
    globalIndex dispColIndices[numUdofs];

    /// C-array storage for the element local Ru residual vector.
    real64 localRu[numUdofs];

    /// C-array storage for the element local AutAtu matrix.
    real64 localAutAtu[numUdofs][numUdofs];

    /// C-array storage for the element local Atu matrix.
    real64 localAtu[numTdofs][numUdofs];

    /// Stack storage for the element local lagrange multiplier vector
    real64 tLocal[numTdofs];

    /// C-array storage for rotation matrix
    real64 localRotationMatrix[3][3];

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
  void setup( localIndex const k,
              StackVariables & stack ) const
  {
    constexpr int shift = numNodesPerElem * 3;

    localIndex const kf0 = m_elemsToFaces[k][0];
    localIndex const kf1 = m_elemsToFaces[k][1];
    for( localIndex a=0; a<numNodesPerElem; ++a )
    {
      localIndex const kn0 = m_faceToNodes( kf0, FE_TYPE::permutation[ a ] );
      localIndex const kn1 = m_faceToNodes( kf1, FE_TYPE::permutation[ a ] );

      //std::cout << kn0 << " " << kn1 << std::endl;
      for( int i=0; i<3; ++i )
      {
        stack.dispEqnRowIndices[a*3+i] = m_dofNumber[kn0]+i-m_dofRankOffset;
        stack.dispEqnRowIndices[shift + a*3+i] = m_dofNumber[kn1]+i-m_dofRankOffset;
        stack.dispColIndices[a*3+i] = m_dofNumber[kn0]+i;
        stack.dispColIndices[shift + a*3+i] = m_dofNumber[kn1]+i;
        stack.X[ a ][ i ] = m_X[  kn0 ][ i ];
      }
    }

    for( int j=0; j<3; ++j )
    {
       for( int i=0; i<3; ++i )
       {
         stack.localRotationMatrix[ i ][ j ] = m_rotationMatrix( k, i, j );
       }
    }

    LvArray::tensorOps::fill< stack.numTdofs, stack.numUdofs >( stack.localAtu, 0.0 );  //make 0
    LvArray::tensorOps::fill< stack.numUdofs, stack.numUdofs >( stack.localAutAtu, 0.0 );  //make 0
  }

  GEOS_HOST_DEVICE
  inline
  void quadraturePointKernel( localIndex const k,
                              localIndex const q,
                              StackVariables & stack ) const
  {
    GEOS_UNUSED_VAR( k );
    real64 const detJ = m_finiteElementSpace.template transformedQuadratureWeight( q, stack.X );

    std::cout << "detJ: " << detJ << std::endl;

    real64 N[ numNodesPerElem ]; 
    m_finiteElementSpace.template calcN( q, N );

    solidMechanicsALMKernelsHelper::accumulateAtuLocalOperator<stack.numTdofs, 
                                                               stack.numUdofs, 
                                                               numNodesPerElem>(stack.localAtu, 
                                                                                N, 
                                                                                FE_TYPE::permutation,
                                                                                detJ);
  }

  GEOS_HOST_DEVICE
  inline
  real64 complete( localIndex const k,
                   StackVariables & stack ) const
  {
    GEOS_UNUSED_VAR( k );

    real64 matRtAtu[3][stack.numUdofs], matRRtAtu[3][stack.numUdofs];

    // transp(R) * Atu
    LvArray::tensorOps::Rij_eq_AkiBkj< 3, stack.numUdofs, 3 >( matRtAtu, stack.localRotationMatrix, 
                                                               stack.localAtu );

    //std::cout << "matrixAtu: " << std::endl;
    //for (int i=0; i<3; ++i)
    //{
    //  for (int j=0; j<stack.numUdofs; ++j)
    //  {
    //    std::cout << matRtAtu[ i ] [ j ] << " ";
    //  }
    //  std::cout << std::endl;
    //}
    //std::cout << std::endl;

    // R*RtAtu 
    LvArray::tensorOps::Rij_eq_AikBkj< 3, stack.numUdofs, 3 >( matRRtAtu, stack.localRotationMatrix, 
                                                               matRtAtu );

    // transp(Atu)*RRtAtu
    LvArray::tensorOps::Rij_eq_AkiBkj< stack.numUdofs, stack.numUdofs, 3 >( stack.localAutAtu, stack.localAtu, 
                                                                            matRRtAtu); 
                                                                          
    //for (int i=0; i<stack.numUdofs; ++i)
    //{
    //  for (int j=0; j<stack.numUdofs; ++j)
    //  {
    //    std::cout << stack.localAutAtu[ i ] [ j ] << " ";
    //  }
    //  std::cout << std::endl;
    //}
    //std::cout << std::endl;
    //abort();

    for( localIndex i=0; i < stack.numUdofs; ++i )
    {
      localIndex const dof = LvArray::integerConversion< localIndex >( stack.dispEqnRowIndices[ i ] );

      if( dof < 0 || dof >= m_matrix.numRows() ) continue;

      //std::cout << "Add elements: " << dof << " "; 
      //for( localIndex j=0; j < stack.numUdofs; ++j )
      //  std::cout << stack.dispColIndices[j] << " "; 
      //std::cout << std::endl;
      // fill in matrix
      m_matrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( dof,
                                                                              stack.dispColIndices,
                                                                              stack.localAutAtu[i],
                                                                              stack.numUdofs );
    }
   
    return 0.0;
  }

protected:

  /// The array containing the nodal position array.
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const m_X;

  ArrayOfArraysView< localIndex const > const m_faceToNodes; 

  ArrayOfArraysView< localIndex const > const m_elemsToFaces;

  arrayView1d< localIndex const > const m_faceElementList;

  arrayView3d< real64 const > const m_rotationMatrix;

  arrayView2d< real64 const > const m_deltaTraction;

};

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

    GEOS_UNUSED_VAR(faceNormal, elemsToFaces, rotationMatrix);

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

} // namespace SolidMechanicsALMKernels

} // namespace geos


#endif /* GEOS_PHYSICSSOLVERS_CONTACT_SOLIDMECHANICSALMKERNELS_HPP_ */