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
 * @file SolidMechanicsEFEMKernels.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSEFEMKERNELS_HPP_
#define GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSEFEMKERNELS_HPP_

#include "SolidMechanicsSmallStrainQuasiStaticKernel.hpp"
#include "EFEMKernelsHelper.hpp"

namespace geosx
{

namespace SolidMechanicsEFEMKernels
{

/**
 * @brief Implements kernels for solving quasi-static equilibrium.
 * @copydoc geosx::finiteElement::ImplicitKernelBase
 * @tparam NUM_NODES_PER_ELEM The number of nodes per element for the
 *                            @p SUBREGION_TYPE.
 * @tparam UNUSED An unused parameter since we are assuming that the test and
 *                trial space have the same number of support points.
 *
 * ### QuasiStatic Description
 *
 */
template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
class AssumedEnhancedStrain :
  public SolidMechanicsLagrangianFEMKernels::QuasiStatic< SUBREGION_TYPE,
                                                          CONSTITUTIVE_TYPE,
                                                          FE_TYPE >
{
public:
  /// Alias for the base class;
  using Base = SolidMechanicsLagrangianFEMKernels::QuasiStatic< SUBREGION_TYPE,
                                                                CONSTITUTIVE_TYPE,
                                                                FE_TYPE >;

  /// Number of nodes per element...which is equal to the
  /// numTestSupportPointPerElem and numTrialSupportPointPerElem by definition.
  static constexpr int numNodesPerElem = Base::numTestSupportPointsPerElem;
  using Base::numDofPerTestSupportPoint;
  using Base::numDofPerTrialSupportPoint;
  using Base::m_dofNumber;
  using Base::m_dofRankOffset;
  using Base::m_matrix;
  using Base::m_rhs;
  using Base::m_elemsToNodes;
  using Base::m_constitutiveUpdate;
  using Base::m_finiteElementSpace;


  /**
   * @brief Constructor
   * @copydoc geosx::finiteElement::ImplicitKernelBase::ImplicitKernelBase
   * @param inputGravityVector The gravity vector.
   */
  AssumedEnhancedStrain( NodeManager const & nodeManager,
                         EdgeManager const & edgeManager,
                         FaceManager const & faceManager,
                         SUBREGION_TYPE const & elementSubRegion,
                         FE_TYPE const & finiteElementSpace,
                         CONSTITUTIVE_TYPE * const inputConstitutiveType,
                         EmbeddedSurfaceSubRegion const & embeddedSurfSubRegion,
                         arrayView1d< globalIndex const > const & uDofNumber,
                         arrayView1d< globalIndex const > const & wDofNumber,
                         globalIndex const rankOffset,
                         CRSMatrixView< real64, globalIndex const > const & inputMatrix,
                         arrayView1d< real64 > const & inputRhs,
                         real64 const (&inputGravityVector)[3] ):
    Base( nodeManager,
          edgeManager,
          faceManager,
          elementSubRegion,
          finiteElementSpace,
          inputConstitutiveType,
          uDofNumber,
          rankOffset,
          inputMatrix,
          inputRhs,
          inputGravityVector ),
    m_wDofNumber( wDofNumber ),
    m_nVec( embeddedSurfSubRegion.getNormalVector() ),
    m_tVec1( embeddedSurfSubRegion.getTangentVector1() ),
    m_tVec2( embeddedSurfSubRegion.getTangentVector2() ),
    m_surfaceCenter( embeddedSurfSubRegion.getElementCenter()),
    m_surfaceArea( embeddedSurfSubRegion.getElementArea()),
    m_elementVolume( elementSubRegion.getElementVolume()),
    m_fracturedElems( elementSubRegion.fracturedElementsList() ),
    m_cellsToEmbeddedSurfaces( elementSubRegion.embeddedSurfacesList() )
  {}

  //***************************************************************************
  /**
   * @copydoc finiteElement::KernelBase::StackVariables
   */
  struct StackVariables : public Base::StackVariables
  {
public:
    /// The number of displacement dofs per element.
    static constexpr int numUdofs = numNodesPerElem * 3;

    /// The number of jump dofs per element.
    static constexpr int numWdofs = 3;

    /**
     * Default constructor
     */
    GEOSX_HOST_DEVICE
    StackVariables():
      Base::StackVariables(),
            dispEqnRowIndices{ 0 },
      dispColIndices{ 0 },
      jumpEqnRowIndices{ 0 },
      jumpColIndices{ 0 },
      localRu{ 0.0 },
      localRw{ 0.0 },
      localKww{ { 0.0 } },
      localKwu{ { 0.0 } },
      localKuw{ { 0.0 } },
      wLocal(),
	  uLocal(),
      hInv()
    {}

    /// C-array storage for the element local row degrees of freedom.
    globalIndex dispEqnRowIndices[numUdofs];

    /// C-array storage for the element local column degrees of freedom.
    globalIndex dispColIndices[numUdofs];

    /// C-array storage for the element local row degrees of freedom.
    globalIndex jumpEqnRowIndices[numWdofs];

    /// C-array storage for the element local column degrees of freedom.
    globalIndex jumpColIndices[numWdofs];

    /// C-array storage for the element local Ru residual vector.
    real64 localRu[numUdofs];

    /// C-array storage for the element local Rw residual vector.
    real64 localRw[numWdofs];

    /// C-array storage for the element local Kww matrix.
    real64 localKww[numWdofs][numWdofs];

    /// C-array storage for the element local Kwu matrix.
    real64 localKwu[numWdofs][numUdofs];

    /// C-array storage for the element local Kuw matrix.
    real64 localKuw[numUdofs][numWdofs];

    /// Stack storage for the element local jump vector
    real64 wLocal[3];

    /// Stack storage for the elenta displacement vector.
    real64 uLocal[numUdofs];

    /// Stack storage for Area/Volume
    real64 hInv;
  };
  //***************************************************************************

  /**
   * @copydoc ::geosx::finiteElement::KernelBase::kernelLaunch
   *
   * @detail it uses the kernelLaunch interface of KernelBase but it only launches the kernel
   * on the set of fractured elements within the subregion.
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
    GEOSX_MARK_FUNCTION;

    GEOSX_UNUSED_VAR( numElems );

    // Define a RAJA reduction variable to get the maximum residual contribution.
    RAJA::ReduceMax< ReducePolicy< POLICY >, real64 > maxResidual( 0 );

    forAll< POLICY >( m_fracturedElems.size(),
                      [=] GEOSX_HOST_DEVICE ( localIndex const i )
    {
      localIndex k = m_fracturedElems[i];
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
   * @copydoc ::geosx::finiteElement::ImplicitKernelBase::setup
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void setup( localIndex const k,
              StackVariables & stack ) const
  {

    localIndex const embSurfIndex = m_cellsToEmbSurfaces[k];

    stack.hInv = m_surfaceArea[embSurfIndex] / m_elementVolume[k];
    for( localIndex a=0; a<numNodesPerElem; ++a )
    {
      localIndex const localNodeIndex = m_elemsToNodes( k, a );

      for( int i=0; i<3; ++i )
      {
        stack.dispEqnRowIndices[a*3+i] = m_dofNumber[localNodeIndex]+i-m_dofRankOffset;
        stack.dispColIndices[a*3+i]    = m_dofNumber[localNodeIndex]+i;
        stack.xLocal[ a ][ i ] = m_X[ localNodeIndex ][ i ];
        stack.uLocal[ a*3 + i ] = m_disp[localNodeIndex][i];

      }
    }

    for( int i=0; i<3; ++i )
    {
      // need to grab the index.
      stack.jumpEqnRowIndices = m_wDofNumber[embSurfIndex] + i - m_dofRankOffset;
      stack.jumpColIndices    = m_wDofNumber[embSurfIndex] + i;
      stack.wLocal[ i ] = m_w[ m_wDofNumber[embSurfIndex] + i ];
    }
  }

  /**
   * @copydoc geosx::finiteElement::KernelBase::quadraturePointKernel
   * @tparam STRESS_MODIFIER Type of optional functor to allow for the
   * modification of stress prior to integration.
   * @param stressModifier An optional functor to allow for the modification
   *  of stress prior to integration.
   */
  template< typename STRESS_MODIFIER = NoOpFunctors >
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void quadraturePointKernel( localIndex const k,
                              localIndex const q,
                              StackVariables & stack,
                              STRESS_MODIFIER && stressModifier = NoOpFunctors{} ) const
  {
    real64 dNdX[ numNodesPerElem ][ 3 ];
    real64 const detJ = m_finiteElementSpace.template getGradN< FE_TYPE >( k, q, stack.xLocal, dNdX );

    constexpr int nUdof = numNodesPerElem*3;

    real64 Kww_gauss[3][3], Kwu_gauss[3][nUdof], Kuw_gauss[nUdof][3];
    real64 compMatrix[6][3], strainMatrix[6][nUdof], eqMatrix[3][6];
    real64 matBD[nUdof][6], matED[3][6];

    int Heaviside[ numNodesPerElem ];

    m_constitutiveUpdate.GetStiffness( k, q, stack.constitutiveStiffness );

    EFEMKernelsHelper::ComputeHeavisideFunction< numNodesPerElem >( Heaviside, stack.xLocal, m_nVec[k], m_surfaceCenter[k] );


    EFEMKernelsHelper::AssembleEquilibriumOperator( eqMatrix,
                                                    m_nVec[k],
                                                    m_tVec1[k],
                                                    m_tVec2[k],
                                                    stack.hInv );

    EFEMKernelsHelper::AssembleCompatibilityOperator< numNodesPerElem >( compMatrix,
                                                                         m_nVec[k],
                                                                         m_tVec1[k],
                                                                         m_tVec2[k],
                                                                         Heaviside,
                                                                         dNdX );

    EFEMKernelsHelper::AssembleStrainOperator< 6, nUdof, numNodesPerElem >( strainMatrix, dNdX );

    // transp(B)D
    LvArray::tensorOps::Rij_eq_AkiBkj< nUdof, 6, nUdof >( matBD, strainMatrix, stack.constitutiveStiffness );
    // ED
    LvArray::tensorOps::Rij_eq_AikBkj< 3, 6, 3 >( matED, eqMatrix, stack.constitutiveStiffness );
    // EDC
    LvArray::tensorOps::Rij_eq_AikBkj< 3, 3, 6 >( Kww_gauss, matED, compMatrix );
    // EDB
    LvArray::tensorOps::Rij_eq_AikBkj< 3, nUdof, 6 >( Kwu_gauss, matED, strainMatrix );
    // transp(B)DB
    LvArray::tensorOps::Rij_eq_AikBkj< nUdof, 3, 6 >( Kuw_gauss, matBD, compMatrix );

    // multiply by determinant
    LvArray::tensorOps::scale< 3, 3 >( Kww_gauss, -detJ );
    LvArray::tensorOps::scale< 3, nUdof >( Kwu_gauss, -detJ );
    LvArray::tensorOps::scale< nUdof, 3 >( Kuw_gauss, -detJ );

    // TODO add a scale add for matrices to the tensorOps.
    // Add Gauss point contribution to element matrix
    LvArray::tensorOps::add< 3, 3 >( stack.localKww, Kww_gauss );
    LvArray::tensorOps::add< 3, nUdof >( stack.localKwu, Kwu_gauss );
    LvArray::tensorOps::add< nUdof, 3 >( stack.localKuw, Kuw_gauss );

//     real64 strainInc[6] = {0};
//     // compute strainInc due to the fracture jump;
//     LvArray::tensorOps::Ri_add_AijBj(strainInc, compMatrix, stack.wlocal);
//     m_constitutiveUpdate.SmallStrain( k, q, strainInc );
//
//     typename CONSTITUTIVE_TYPE::KernelWrapper::DiscretizationOps stiffnessHelper;
//     m_constitutiveUpdate.setDiscretizationOps( k, q, stiffnessHelper );
  }

  /**
   * @copydoc geosx::finiteElement::ImplicitKernelBase::complete
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  real64 complete( localIndex const k,
                   StackVariables & stack ) const
  {
    GEOSX_UNUSED_VAR( k );
    real64 maxForce = 0;
    constexpr int nUdof = numNodesPerElem*3;

    // Compute the local residuals
    LvArray::tensorOps::Ri_add_AijBj< 3, 3 >( stack.localRw, stack.localKww, stack.wLocal );
    LvArray::tensorOps::Ri_add_AijBj< 3, nUdof >( stack.localRw, stack.localKwu, stack.uLocal );
    LvArray::tensorOps::Ri_add_AijBj< nUdof, 3 >( stack.localRu, stack.localKuw, stack.wLocal );

    for( localIndex i = 0; i < nUdof; ++i )
    {
      if( dispEqnRowIndices[i] >= 0 && dispEqnRowIndices[i] < m_matrix.numRows() )
      {
        RAJA::atomicAdd< parallelDeviceAtomic >( &localRhs[dispEqnRowIndices[i]], R0[i] );

        m_matrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( stack.dispEqnRowIndices[i],
                                                                       stack.jumpColIndices.data(),
                                                                       stack.localKuw[i],
                                                                       3 );
      }
    }

    for( localIndex i=0; i < 3; ++i )
    {
      if( jumpEqnRowIndices[i] >= 0 && jumpEqnRowIndices[i] < m_matrix.numRows() )
      {
        RAJA::atomicAdd< parallelDeviceAtomic >( &localRhs[jumpEqnRowIndices[i]], R1[i] );

        // fill in matrix
        m_matrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( stack.jumpEqnRowIndices[i],
                                                                       stack.jumpColIndices.data(),
                                                                       stack.localKww[i],
                                                                       3 );

        m_matrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( stack.jumpEqnRowIndices[i],
                                                                       stack.dispColIndices.data(),
                                                                       stack.localKwu[i],
                                                                       numNodesPerElem*3 );
      }
    }

  }

protected:

  arrayView1d< globalIndex const > const m_wDofNumber;

  arrayView2d< real64 const > const m_nVec;

  arrayView2d< real64 const > const m_tVec1;

  arrayView2d< real64 const > const m_tVec2;

  arrayView2d< real64 const > const m_surfaceCenter;

  arrayView1d< real64 const > const m_surfaceArea;

  arrayView1d< real64 const > const m_elementVolume;

  SortedArrayView< localIndex const > const m_fracturedElems;

  ArrayOfArraysView< localIndex const > const m_cellsToEmbeddedSurfaces;
};

} // namespace SolidMechanicsEFEMKernels

} // namespace geosx


#endif /* GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSEFEMKERNELS_HPP_ */
