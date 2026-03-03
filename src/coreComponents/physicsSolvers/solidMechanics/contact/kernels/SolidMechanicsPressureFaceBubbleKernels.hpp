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
 * @file SolidMechanicsPressureFaceBubbleKernels.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_CONTACT_KERNELS_SOLIDMECHANICSPRESSUREFACEBUBBLEKERNELS_HPP_
#define GEOS_PHYSICSSOLVERS_CONTACT_KERNELS_SOLIDMECHANICSPRESSUREFACEBUBBLEKERNELS_HPP_

#include "finiteElement/kernelInterface/ImplicitKernelBase.hpp"
#include "SolidMechanicsConformingContactKernelsHelper.hpp"

namespace geos
{

namespace solidMechanicsConformingContactKernels
{

/**
 * @brief Implements kernels for computing the pressure contribution
 *        given by the bubble face functions to the balance of momentum equation.
 * @copydoc geos::finiteElement::ImplicitKernelBase
 *
 */
template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE,
          typename MATRIX_VIEW = CRSMatrixView< real64, globalIndex const > >
class PressureFaceBubbleKernels :
  public finiteElement::ImplicitKernelBase< SUBREGION_TYPE,
                                            CONSTITUTIVE_TYPE,
                                            FE_TYPE,
                                            3,
                                            3,
                                            MATRIX_VIEW >
{
public:

  /// Alias for the base class;
  using Base = finiteElement::ImplicitKernelBase< SUBREGION_TYPE,
                                                  CONSTITUTIVE_TYPE,
                                                  FE_TYPE,
                                                  3,
                                                  3,
                                                  MATRIX_VIEW >;

  /// Number of nodes per element, which is equal to the
  /// numTestSupportPointPerElem and numTrialSupportPointPerElem by definition.
  static constexpr int numNodesPerElem = Base::maxNumTestSupportPointsPerElem;

  /// Compile time value for the number of faces per element.
  static constexpr int numFacesPerElem = FE_TYPE::numFaces;

  /// Compile time value for the number of quadrature points per element.
  static constexpr int numQuadraturePointsPerElem = FE_TYPE::numQuadraturePoints;

  using Base::m_elemsToNodes;
  using Base::m_finiteElementSpace;
  using Base::m_constitutiveUpdate;
  using Base::m_dofNumber;
  using Base::m_dofRankOffset;
  using Base::m_matrix;
  using Base::m_rhs;
  using Base::m_dt;

  /**
   * @brief Constructor
   * @copydoc geos::finiteElement::ImplicitKernelBase::ImplicitKernelBase
   */
  PressureFaceBubbleKernels( NodeManager const & nodeManager,
                             EdgeManager const & edgeManager,
                             FaceManager const & faceManager,
                             localIndex const targetRegionIndex,
                             SUBREGION_TYPE const & elementSubRegion,
                             FE_TYPE const & finiteElementSpace,
                             CONSTITUTIVE_TYPE & inputConstitutiveType,
                             arrayView1d< globalIndex const > const uDofNumber,
                             arrayView1d< globalIndex const > const bDofNumber,
                             globalIndex const rankOffset,
                             MATRIX_VIEW const inputMatrix,
                             arrayView1d< real64 > const inputRhs,
                             real64 const inputDt ):
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
    m_incrementalDisp( nodeManager.getField< fields::solidMechanics::incrementalDisplacement >() ),
    m_bDofNumber( bDofNumber ),
    m_bubbleElems( elementSubRegion.bubbleElementsList() ),
    m_elemsToFaces( elementSubRegion.faceElementsList() ),
    m_pressure( elementSubRegion.template getField< fields::flow::pressure >().toViewConst() ),
    m_pressure_n( elementSubRegion.template getField< fields::flow::pressure_n >().toViewConst() ),
    m_incrementalBubbleDisp( faceManager.getField< fields::contact::incrementalBubbleDisplacement >().toViewConst() )
  {}

  //***************************************************************************

  /**
   * @copydoc finiteElement::ImplicitKernelBase::StackVariables
   */
  struct StackVariables
  {
public:

    /// The number of displacement dofs per element.
    static constexpr int numUdofs = numNodesPerElem * 3;

    /// The number of jump dofs per element.
    static constexpr int numBubbleUdofs = numFacesPerElem * 3;

    /// The number of pressure dofs per element.
    static constexpr int numPdofs = 1;

    /**
     * Default constructor
     */
    GEOS_HOST_DEVICE
    StackVariables():
      bEqnRowIndices{},
      localRb{},
      X{ {} },
      pLocal{},
      uIncrLocal{},
      bIncrLocal{}
    {}

    /// C-array storage for the element local row degrees of freedom.
    globalIndex bEqnRowIndices[3];

    /// C-array storage for the element local Rb residual vector.
    real64 localRb[numBubbleUdofs];

    /// local nodal coordinates
    real64 X[ numNodesPerElem ][ 3 ];

    /// local pressure
    real64 pLocal[numPdofs];

    /// local incremental displacement
    real64 uIncrLocal[numUdofs];

    /// incremental bubble displacement
    real64 bIncrLocal[3];

  };

  //***************************************************************************

  /**
   * @copydoc ::geos::finiteElement::KernelBase::kernelLaunch
   *
   * @detail it uses the kernelLaunch interface of KernelBase but it only launches the kernel
   * on the set of elements that have bubble dof within the subregion.
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

    forAll< POLICY >( kernelComponent.m_bubbleElems.size(),
                      [=] GEOS_HOST_DEVICE ( localIndex const i )
    {
      typename KERNEL_TYPE::StackVariables stack;

      kernelComponent.setup( i, stack );
      for( integer q=0; q<numQuadraturePointsPerElem; ++q )
      {
        kernelComponent.quadraturePointKernel( i, q, stack );
      }
      maxResidual.max( kernelComponent.complete( i, stack ) );
    } );

    return maxResidual.get();
  }
  //END_kernelLauncher

  GEOS_HOST_DEVICE
  inline
  void setup( localIndex const kk,
              StackVariables & stack ) const
  {

    localIndex k = m_bubbleElems[kk];

    for( localIndex a=0; a<numNodesPerElem; ++a )
    {
      localIndex const localNodeIndex = m_elemsToNodes( k, a );

      for( int i=0; i<3; ++i )
      {
        stack.X[ a ][ i ] = m_X[ localNodeIndex ][ i ];
        stack.uIncrLocal[ a*3 + i ] = m_incrementalDisp[ localNodeIndex ][i];
      }
    }

    localIndex const localFaceIndex = m_elemsToFaces[kk][0];

    for( int i=0; i<3; ++i )
    {
      // need to grab the index.
      stack.bEqnRowIndices[i] = m_bDofNumber[localFaceIndex] + i - m_dofRankOffset;
      stack.bIncrLocal[ i ] = m_incrementalBubbleDisp[ localFaceIndex ][i];
    }

    stack.pLocal[0] = m_pressure( k );

  }


  GEOS_HOST_DEVICE
  inline
  void quadraturePointKernel( localIndex const kk,
                              localIndex const q,
                              StackVariables & stack ) const
  {

    localIndex k = m_bubbleElems[kk];

    constexpr int nBubbleUdof = numFacesPerElem*3;

    constexpr int nUdof = numNodesPerElem*3;

    real64 dBubbleNdX[ numFacesPerElem ][ 3 ];
    // Next line is needed because I only inserted a placeholder for calcGradFaceBubbleN in some finite elements
    LvArray::tensorOps::fill< numFacesPerElem, 3 >( dBubbleNdX, 0 );

    real64 detJ = m_finiteElementSpace.calcGradFaceBubbleN( q, stack.X, dBubbleNdX );

    real64 dNdX[ numNodesPerElem ][ 3 ];
    detJ = m_finiteElementSpace.calcGradN( q, stack.X, dNdX );

    real64 biotCoefficient;
    m_constitutiveUpdate.getBiotCoefficient( k, biotCoefficient );

    real64 strainBubbleMatrix[6][nBubbleUdof];
    solidMechanicsConformingContactKernelsHelper::assembleStrainOperator< 6, nBubbleUdof, numFacesPerElem >( strainBubbleMatrix, dBubbleNdX );

    real64 strainMatrix[6][nUdof];
    solidMechanicsConformingContactKernelsHelper::assembleStrainOperator< 6, nUdof, numNodesPerElem >( strainMatrix, dNdX );

    real64 biotPressure[6] = {0};
    LvArray::tensorOps::symAddIdentity< 3 >( biotPressure, -biotCoefficient * stack.pLocal[0] );

    // transp(Bb)(Identity biot pressure)
    real64 Rb_gauss[nBubbleUdof];
    LvArray::tensorOps::Ri_eq_AjiBj< nBubbleUdof, 6 >( Rb_gauss, strainBubbleMatrix, biotPressure );
    LvArray::tensorOps::scaledAdd< nBubbleUdof >( stack.localRb, Rb_gauss, -detJ );

    real64 localStrainBubbleMatrix[6][3];
    localIndex const parentFaceIndex = m_elemsToFaces[kk][1];
    for( localIndex i = 0; i < 6; ++i )
    {
      for( localIndex j = 0; j < 3; ++j )
      {
        localStrainBubbleMatrix[i][j] = strainBubbleMatrix[i][parentFaceIndex*3+j];
      }
    }
    real64 strainIncBubble[6] = {0};
    real64 strainInc[6] = {0};

    LvArray::tensorOps::Ri_eq_AijBj< 6, 3 >( strainIncBubble, localStrainBubbleMatrix, stack.bIncrLocal );
    LvArray::tensorOps::Ri_eq_AijBj< 6, nUdof >( strainInc, strainMatrix, stack.uIncrLocal );

    //LvArray::tensorOps::add< 6 >( strainInc, strainIncBubble );

    real64 totalStress[6] = {0};
    typename CONSTITUTIVE_TYPE::KernelWrapper::DiscretizationOps stiffness;

    m_constitutiveUpdate.smallStrainUpdatePoromechanicsFixedStress( k, q,
                                                                    m_dt,
                                                                    m_pressure[k],
                                                                    m_pressure_n[k],
                                                                    0.0, 0.0, 0.0,
                                                                    strainInc,
                                                                    totalStress,
                                                                    stiffness );

  }

  GEOS_HOST_DEVICE
  inline
  real64 complete( localIndex const kk,
                   StackVariables & stack ) const
  {

    localIndex const parentFaceIndex = m_elemsToFaces[kk][1];

    // Extract only the submatrix and known term corresponding to the index of the local face
    // on which the bubble function was applied.
    real64 localRb[3];
    for( localIndex i = 0; i < 3; ++i )
    {
      localRb[i] = stack.localRb[parentFaceIndex*3+i];
    }

    for( localIndex i=0; i < 3; ++i )
    {
      localIndex const dof = LvArray::integerConversion< localIndex >( stack.bEqnRowIndices[ i ] );

      if( dof < 0 || dof >= m_matrix.numRows() ) continue;

      RAJA::atomicAdd< parallelDeviceAtomic >( &m_rhs[dof], localRb[i] );

    }

    return 0.0;
  }

protected:

  /// The array containing the nodal position array.
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const m_X;

  /// The rank-global incremental displacement array.
  arrayView2d< real64 const, nodes::INCR_DISPLACEMENT_USD > const m_incrementalDisp;

  /// The global degree of freedom number of bubble
  arrayView1d< globalIndex const > const m_bDofNumber;

  /// The array containing the list of bubble elements.
  arrayView1d< localIndex const > const m_bubbleElems;

  /// The array containing the element to bubble face map.
  /// The bubble face is the face of the element to which the bubble is applied.
  /// Both the local and global face index are stored in this array.
  arrayView2d< localIndex const > const m_elemsToFaces;

  /// The array containing the pressure of each element.
  arrayView1d< real64 const > const m_pressure;

  /// The array containing the pressure at the previous time step of each element.
  arrayView1d< real64 const > const m_pressure_n;

  /// The incremental bubble displacement field.
  arrayView2d< real64 const > const m_incrementalBubbleDisp;

};

/// The factory used to construct a QuasiStatic kernel.
using PressureFaceBubbleFactory = finiteElement::KernelFactory< PressureFaceBubbleKernels,
                                                                arrayView1d< globalIndex const > const,
                                                                arrayView1d< globalIndex const > const,
                                                                globalIndex const,
                                                                CRSMatrixView< real64, globalIndex const > const,
                                                                arrayView1d< real64 > const,
                                                                real64 const >;

} // namespace SolidMechanicsPressureFaceBubbleKernels

} // namespace geos


#endif /* GEOS_PHYSICSSOLVERS_CONTACT_KERNELS_SOLIDMECHANICSPRESSUREFACEBUBBLEKERNELS_HPP_ */
