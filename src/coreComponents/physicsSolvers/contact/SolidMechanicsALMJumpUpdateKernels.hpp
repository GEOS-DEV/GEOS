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
 * @file SolidMechanicsALMUpdateKernels.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_CONTACT_SOLIDMECHANICSALMUPDATEKERNELS_HPP_
#define GEOS_PHYSICSSOLVERS_CONTACT_SOLIDMECHANICSALMUPDATEKERNELS_HPP_

#include "SolidMechanicsALMKernelsBase.hpp"

namespace geos
{

namespace solidMechanicsALMKernels
{

template< typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
class ALMJumpUpdate :
  public ALMKernelsBase< CONSTITUTIVE_TYPE,
                         FE_TYPE >
{
public:
  /// Alias for the base class;
  using Base = ALMKernelsBase< CONSTITUTIVE_TYPE,
                               FE_TYPE >;

  static constexpr int numNodesPerElem = Base::maxNumTestSupportPointsPerElem;

  using Base::m_X;
  using Base::m_dofNumber;
  using Base::m_dofRankOffset;
  using Base::m_elemsToFaces;
  using Base::m_faceToNodes;
  using Base::m_rotationMatrix;
  using Base::m_dispJump;

  ALMJumpUpdate( NodeManager const & nodeManager,
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
          bDofNumber,
          rankOffset,
          inputMatrix,
          inputRhs,
          inputDt,
          faceElementList ),
    m_displacement(nodeManager.getField< fields::solidMechanics::totalDisplacement >()),
    m_bubbleDisp( faceManager.getField< fields::solidMechanics::totalBubbleDisplacement >() ),
    m_incrDisp( nodeManager.getField< fields::solidMechanics::incrementalDisplacement >() ),
    m_deltaDispJump( elementSubRegion.getField< fields::contact::deltaDispJump >().toView() )
  {}

  //***************************************************************************

  struct StackVariables : public Base::StackVariables
  {

    /// The number of displacement dofs per element.
    static constexpr int numUdofs = numNodesPerElem * 3 * 2;

    static constexpr int numBdofs = 3 * 2;

    static constexpr int numTdofs = 3;

  public:

    GEOS_HOST_DEVICE
    StackVariables():
      Base::StackVariables(),
      uLocal{},
      bLocal{},
      duLocal{},
      dbLocal{},
      deltaDispJumpLocal{}
    {}

    real64 uLocal[numUdofs];

    real64 bLocal[numBdofs];

    real64 duLocal[numUdofs];

    real64 dbLocal[numBdofs];

    real64 deltaDispJumpLocal[numTdofs];

  };
  //***************************************************************************

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
      //localIndex const kn0 = m_faceToNodes( kf0, FE_TYPE::permutation[ a ] );
      //localIndex const kn1 = m_faceToNodes( kf1, FE_TYPE::permutation[ a ] );
      localIndex const kn0 = m_faceToNodes( kf0, a );
      localIndex const kn1 = m_faceToNodes( kf1, a );

      //std::cout << kn0 << " " << kn1 << std::endl;
      for( int i=0; i<3; ++i )
      {
        //stack.dispEqnRowIndices[a*3+i] = m_dofNumber[kn0]+i-m_dofRankOffset;
        //stack.dispEqnRowIndices[shift + a*3+i] = m_dofNumber[kn1]+i-m_dofRankOffset;
        //stack.dispColIndices[a*3+i] = m_dofNumber[kn0]+i;
        //stack.dispColIndices[shift + a*3+i] = m_dofNumber[kn1]+i;
        stack.X[ a ][ i ] = m_X[ m_faceToNodes( kf0, FE_TYPE::permutation[ a ]) ][ i ];
        stack.uLocal[a*3+i] = m_displacement[kn0][i];
        stack.uLocal[shift + a*3+i] = m_displacement[kn1][i];
        stack.duLocal[a*3+i] = m_incrDisp[kn0][i];
        stack.duLocal[shift + a*3+i] = m_incrDisp[kn1][i];
      }
    }

    for( int j=0; j<3; ++j )
    {
       for( int i=0; i<3; ++i )
       {
         stack.localRotationMatrix[ i ][ j ] = m_rotationMatrix( k, i, j );
       }
    }

    for( int i=0; i<3; ++i )
    {
      stack.bLocal[ i ] = m_bubbleDisp[ kf0 ][i];
      stack.bLocal[ 3 + i ] = m_bubbleDisp[ kf1 ][i];
    }

    //for( int i=0; i<stack.numTdofs; ++i )
    //{
    //  stack.dispJumpLocal[i] = m_dispJump(k, i);
    //  stack.oldDispJumpLocal[i] = m_oldDispJump(k, i);
    //}
  }

  GEOS_HOST_DEVICE
  inline
  real64 complete( localIndex const k,
                   StackVariables & stack ) const
  {
    GEOS_UNUSED_VAR( k );

    //real64 matRtAtu[3][stack.numUdofs];
    real64 matRtAtu[3][stack.numUdofs];
    real64 matRtAtb[3][stack.numBdofs];

    // transp(R) * Atu
    LvArray::tensorOps::Rij_eq_AkiBkj< 3, stack.numUdofs, 3 >( matRtAtu, stack.localRotationMatrix, 
                                                               stack.localAtu );
    LvArray::tensorOps::Rij_eq_AkiBkj< 3, stack.numBdofs, 3 >( matRtAtb, stack.localRotationMatrix, 
                                                               stack.localAtb );

    LvArray::tensorOps::Ri_eq_AijBj< 3, stack.numUdofs >( stack.dispJumpLocal, matRtAtu, stack.uLocal );
    LvArray::tensorOps::Ri_eq_AijBj< 3, stack.numUdofs >( stack.deltaDispJumpLocal, matRtAtu, stack.duLocal );

    LvArray::tensorOps::Ri_add_AijBj< 3, stack.numBdofs >( stack.dispJumpLocal, matRtAtb, stack.bLocal );
    LvArray::tensorOps::Ri_add_AijBj< 3, stack.numBdofs >( stack.deltaDispJumpLocal, matRtAtb, stack.dbLocal );

    for( int i=0; i<3; ++i )
    {
      m_dispJump[ k ][ i ] = stack.dispJumpLocal[ i ];
      m_deltaDispJump[ k ][ i ] = stack.deltaDispJumpLocal[ i ];
    }

    return 0.0;
  }

protected:

  /// The rank-global displacement array.
  arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const m_displacement;

  arrayView2d< real64 const > const m_bubbleDisp;

  /// The rank-global incremental displacement array.
  arrayView2d< real64 const, nodes::INCR_DISPLACEMENT_USD > const m_incrDisp;

  arrayView2d< real64 > const m_deltaDispJump;

};

using ALMJumpUpdateFactory = finiteElement::InterfaceKernelFactory< ALMJumpUpdate,
                                                                    arrayView1d< globalIndex const > const,
                                                                    arrayView1d< globalIndex const > const,
                                                                    globalIndex const,
                                                                    CRSMatrixView< real64, globalIndex const > const,
                                                                    arrayView1d< real64 > const,
                                                                    real64 const, 
                                                                    arrayView1d< localIndex const > const >;

} // namespace SolidMechanicsALMKernels

} // namespace geos

#endif /* GEOS_PHYSICSSOLVERS_CONTACT_SOLIDMECHANICSALMUPDATEKERNELS_HPP_ */