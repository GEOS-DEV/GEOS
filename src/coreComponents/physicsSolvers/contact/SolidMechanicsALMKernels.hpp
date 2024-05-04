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
 * @file SolidMechanicsALMKernels.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_CONTACT_SOLIDMECHANICSALMKERNELS_HPP_
#define GEOS_PHYSICSSOLVERS_CONTACT_SOLIDMECHANICSALMKERNELS_HPP_

#include "SolidMechanicsALMKernelsBase.hpp"

namespace geos
{

namespace solidMechanicsALMKernels
{

template< typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
class ALM :
  public ALMKernelsBase< CONSTITUTIVE_TYPE,
                         FE_TYPE >
{
public:
  /// Alias for the base class;
  using Base = ALMKernelsBase< CONSTITUTIVE_TYPE,
                               FE_TYPE >;
  
  
  static constexpr int numNodesPerElem = Base::maxNumTestSupportPointsPerElem;

  static constexpr int numQuadraturePointsPerElem = FE_TYPE::numQuadraturePoints;

  using Base::m_elemsToFaces;
  using Base::m_faceToNodes;
  using Base::m_dofNumber;
  using Base::m_dofRankOffset;
  using Base::m_X;
  using Base::m_rotationMatrix;
  using Base::m_penalty;
  using Base::m_traction;
  using Base::m_dispJump;
  using Base::m_oldDispJump;
  using Base::m_matrix;
  using Base::m_rhs;
/*
  using Base::m_finiteElementSpace;
  */

  ALM( NodeManager const & nodeManager,
       EdgeManager const & edgeManager,
       FaceManager const & faceManager,
       localIndex const targetRegionIndex,
       FaceElementSubRegion & elementSubRegion,
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
          inputDt,
          faceElementList )
{}

  struct StackVariables: public Base::StackVariables
  {
  public:

    GEOS_HOST_DEVICE
    StackVariables():
      Base::StackVariables()
    {}
/*    /// The number of displacement dofs per element.
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
      localAtu{{}},
      localRotationMatrix{{}},
      localPenalty{{}},
      tLocal{},
      dispJumpLocal{},
      oldDispJumpLocal{},
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

*/
  };

  template< typename POLICY,
            typename KERNEL_TYPE >
  static
  real64
  kernelLaunch( localIndex const numElems,
                KERNEL_TYPE const & kernelComponent )
  {
    return Base::template kernelLaunch< POLICY, KERNEL_TYPE >( numElems, kernelComponent );
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
      //localIndex const kn0_p = m_faceToNodes( kf0, FE_TYPE::permutation[ a ] );
      localIndex const kn0 = m_faceToNodes( kf0, a );
      localIndex const kn1 = m_faceToNodes( kf1, a );

      //std::cout << kn0 << " " << kn1 << std::endl;
      for( int i=0; i<3; ++i )
      {
        stack.dispEqnRowIndices[a*3+i] = m_dofNumber[kn0]+i-m_dofRankOffset;
        stack.dispEqnRowIndices[shift + a*3+i] = m_dofNumber[kn1]+i-m_dofRankOffset;
        stack.dispColIndices[a*3+i] = m_dofNumber[kn0]+i;
        stack.dispColIndices[shift + a*3+i] = m_dofNumber[kn1]+i;
        stack.X[ a ][ i ] = m_X[ m_faceToNodes( kf0, FE_TYPE::permutation[ a ] ) ][ i ];
        //stack.X[ a ][ i ] = m_X[ kn0 ][ i ];
      }
    }

    for( int j=0; j<3; ++j )
    {
       for( int i=0; i<3; ++i )
       {
         stack.localRotationMatrix[ i ][ j ] = m_rotationMatrix( k, i, j );
       }
    }

    //LvArray::tensorOps::fill< stack.numTdofs, stack.numUdofs >( stack.localAtu, 0.0 );  //make 0
    //LvArray::tensorOps::fill< stack.numUdofs, stack.numUdofs >( stack.localAutAtu, 0.0 );  //make 0

    stack.localPenalty[0][0] = -m_penalty(k, 0);
    stack.localPenalty[1][1] = -m_penalty(k, 1);
    stack.localPenalty[2][2] = -m_penalty(k, 1);

    for( int i=0; i<stack.numTdofs; ++i )
    {
      stack.tLocal[i] = m_traction(k, i);
      stack.dispJumpLocal[i] = m_dispJump(k, i);
      stack.oldDispJumpLocal[i] = m_oldDispJump(k, i);
    }
  }

  GEOS_HOST_DEVICE
  inline
  real64 complete( localIndex const k,
                   StackVariables & stack ) const
  {
    GEOS_UNUSED_VAR( k );

    //real64 matRtAtu[3][stack.numUdofs];
    real64 matRRtAtu[3][stack.numUdofs], matDRtAtu[3][stack.numUdofs];
    real64 dispJumpR[stack.numUdofs];
    real64 oldDispJumpR[stack.numUdofs];
    real64 tractionR[stack.numUdofs];

    // transp(R) * Atu
    LvArray::tensorOps::Rij_eq_AkiBkj< 3, stack.numUdofs, 3 >( matRRtAtu, stack.localRotationMatrix, 
                                                               stack.localAtu );

    LvArray::tensorOps::Ri_eq_AjiBj< stack.numUdofs, 3 >( tractionR, matRRtAtu, stack.tLocal );

    //std::cout << "matrixAtu: " << std::endl;
    //for (int i=0; i<3; ++i)
    //{
    //  for (int j=0; j<stack.numUdofs; ++j)
    //  {
    //    std::cout << matRRtAtu[ i ] [ j ] << " ";
    //  }
    //  std::cout << std::endl;
    //}
    //std::cout << std::endl;
    //for( localIndex i=0; i < stack.numUdofs; ++i )
    //{
    //  localIndex const dof = LvArray::integerConversion< localIndex >( stack.dispEqnRowIndices[ i ] );
    //  std::cout << dof << " ";
    //}
    //std::cout << std::endl;
    //abort();

    // D*RtAtu 
    LvArray::tensorOps::Rij_eq_AikBkj< 3, stack.numUdofs, 3 >( matDRtAtu, stack.localPenalty, 
                                                               matRRtAtu );

    // R*RtAtu 
    LvArray::tensorOps::Rij_eq_AikBkj< 3, stack.numUdofs, 3 >( matRRtAtu, stack.localRotationMatrix, 
                                                               matDRtAtu );

    // transp(Atu)*RRtAtu
    LvArray::tensorOps::Rij_eq_AkiBkj< stack.numUdofs, stack.numUdofs, 3 >( stack.localAutAtu, stack.localAtu, 
                                                                            matRRtAtu); 

    // Compute the local residuals
    LvArray::tensorOps::Ri_eq_AjiBj< stack.numUdofs, 3 >( dispJumpR, matDRtAtu, stack.dispJumpLocal );
    LvArray::tensorOps::Ri_eq_AjiBj< stack.numUdofs, 3 >( oldDispJumpR, matDRtAtu, stack.oldDispJumpLocal );

    LvArray::tensorOps::scaledAdd< stack.numUdofs >( stack.localRu, tractionR, -1 );
    LvArray::tensorOps::scaledAdd< stack.numUdofs >( stack.localRu, dispJumpR,  1 );
    LvArray::tensorOps::scaledAdd< stack.numUdofs >( stack.localRu, oldDispJumpR, -1 );
                                                                          
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

      // Is it needed?
      RAJA::atomicAdd< parallelDeviceAtomic >( &m_rhs[dof], stack.localRu[i] );

      //if (k==1)
      //  std::cout << "Add elements: " << dof << " "; 
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

};

using ALMFactory = finiteElement::InterfaceKernelFactory< ALM,
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

} // namespace SolidMechanicsALMKernels

} // namespace geos


#endif /* GEOS_PHYSICSSOLVERS_CONTACT_SOLIDMECHANICSALMKERNELS_HPP_ */