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
 * @file SinglePhasePoroelasticKernel.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASEPOROELASTICEFEMKERNEL_HPP_
#define GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASEPOROELASTICEFEMKERNEL_HPP_
#include "finiteElement/kernelInterface/ImplicitKernelBase.hpp"
#include "SinglePhasePoroelasticKernel.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBase.hpp"


namespace geosx
{

namespace PoroelasticEFEMKernels
{

/**
 * @brief Implements kernels for solving quasi-static single-phase poromechanics.
 * @copydoc geosx::finiteElement::ImplicitKernelBase
 * @tparam NUM_NODES_PER_ELEM The number of nodes per element for the
 *                            @p SUBREGION_TYPE.
 * @tparam UNUSED An unused parameter since we are assuming that the test and
 *                trial space have the same number of support points.
 *
 * ### SinglePhasePoroelastic Description
 * Implements the KernelBase interface functions required for solving the
 * quasi-static single-phase poromechanics problem using one of the
 * "finite element kernel application" functions such as
 * geosx::finiteElement::RegionBasedKernelApplication.
 *
 */
template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
class SinglePhase :
  public PoroelasticKernels::SinglePhase< SUBREGION_TYPE,
                                          CONSTITUTIVE_TYPE,
                                          FE_TYPE,
                                          3,
                                          3 >
{
public:
  /// Alias for the base class;
  using Base = PoroelasticKernels::SinglePhase< SUBREGION_TYPE,
                                                CONSTITUTIVE_TYPE,
                                                FE_TYPE,
                                                3,
                                                3 >;

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


  SinglePhase( NodeManager const & nodeManager,
               EdgeManager const & edgeManager,
               FaceManager const & faceManager,
               localIndex const targetRegionIndex,
               SUBREGION_TYPE const & elementSubRegion,
               FE_TYPE const & finiteElementSpace,
               CONSTITUTIVE_TYPE & inputConstitutiveType,
               arrayView1d< globalIndex const > const & inputDispDofNumber,
               string const & inputFlowDofKey,
               globalIndex const rankOffset,
               CRSMatrixView< real64, globalIndex const > const & inputMatrix,
               arrayView1d< real64 > const & inputRhs,
               real64 const (&inputGravityVector)[3],
               arrayView1d< string const > const & fluidModelNames ):
    Base( nodeManager,
          edgeManager,
          faceManager,
          targetRegionIndex,
          elementSubRegion,
          finiteElementSpace,
          inputConstitutiveType,
          inputDispDofNumber,
          rankOffset,
          inputMatrix,
          inputRhs,
          inputGravityVector,
          fluidModelNames )
  {}

  //*****************************************************************************
  /**
   * @class StackVariables
   * @copydoc geosx::finiteElement::ImplicitKernelBase::StackVariables
   *
   * Adds a stack array for the displacement, incremental displacement, and the
   * constitutive stiffness.
   */
  struct StackVariables : public Base::StackVariables
  {
public:

    static constexpr int numDispDofPerElem =  Base::StackVariables::numRows;

    // The number of displacement dofs per element.
    static constexpr int numUdofs = numDispDofPerElem * 3;

    /// The number of jump dofs per element.
    static constexpr int numWdofs = 3;

    /// Constructor.
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
      hInv(),
      X()
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

    /// local nodal coordinates
    real64 X[ numNodesPerElem ][ 3 ];

    /// Stack storage for the traction
    real64 tractionVec[3];

    /// Stack storage for the derivative of the traction
    real64 dTractiondw[3][3];
  };
  //*****************************************************************************

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

      forAll< POLICY >( kernelComponent.m_fracturedElems.size(),
                        [=] GEOSX_HOST_DEVICE ( localIndex const i )
      {
        localIndex k = kernelComponent.m_fracturedElems[i];
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
   *
   * For the SinglePhase implementation, global values from the displacement,
   * incremental displacement, and degree of freedom numbers are placed into
   * element local stack storage.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void setup( localIndex const k,
              StackVariables & stack ) const
  {
    for( localIndex a=0; a<numNodesPerElem; ++a )
    {
      localIndex const localNodeIndex = m_elemsToNodes( k, a );

      for( int i=0; i<3; ++i )
      {
#if defined(CALC_FEM_SHAPE_IN_KERNEL)
        stack.xLocal[a][i] = m_X[localNodeIndex][i];
#endif
        stack.u_local[a][i] = m_disp[localNodeIndex][i];
        stack.uhat_local[a][i] = m_uhat[localNodeIndex][i];
        stack.localRowDofIndex[a*3+i] = m_dofNumber[localNodeIndex]+i;
        stack.localColDofIndex[a*3+i] = m_dofNumber[localNodeIndex]+i;
      }
    }

    for( int flowDofIndex=0; flowDofIndex<1; ++flowDofIndex )
    {
      stack.localFlowDofIndex[flowDofIndex] = m_flowDofNumber[k] + flowDofIndex;
    }

  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void quadraturePointKernel( localIndex const k,
                              localIndex const q,
                              StackVariables & stack ) const
  {
    // Get displacement: (i) basis functions (N), (ii) basis function
    // derivatives (dNdX), and (iii) determinant of the Jacobian transformation
    // matrix times the quadrature weight (detJxW)
    real64 N[numNodesPerElem];
    real64 dNdX[numNodesPerElem][3];
    FE_TYPE::calcN( q, N );
    real64 const detJxW = m_finiteElementSpace.template getGradN< FE_TYPE >( k, q, stack.xLocal, dNdX );

    // Evaluate total stress tensor
    real64 strainIncrement[6] = {0};
    real64 totalStress[6];
    real64 porosityNew, dPorosity_dPressure, dPorosity_dVolStrainIncrement;

    // --- Update total stress tensor
    typename CONSTITUTIVE_TYPE::KernelWrapper::DiscretizationOps stiffness;
    FE_TYPE::symmetricGradient( dNdX, stack.uhat_local, strainIncrement );

    m_constitutiveUpdate.smallStrainUpdate_porosity( k, q,
                                                     m_fluidPressure[k],
                                                     m_deltaFluidPressure[k],
                                                     strainIncrement,
                                                     porosityNew,
                                                     dPorosity_dPressure,
                                                     dPorosity_dVolStrainIncrement,
                                                     totalStress,
                                                     stiffness );

    // EFEM part starts here
    constexpr int nUdof = numNodesPerElem*3;

    // Gauss contribution to Kww, Kwu and Kuw blocks
    real64 Kww_gauss[3][3], Kwu_gauss[3][nUdof], Kuw_gauss[nUdof][3];

    //  Compatibility, equilibrium and strain operators. The compatibility operator is constructed as
    //  a 3 x 6 because it is more convenient for construction purposes (reduces number of local var).
    real64 compMatrix[3][6], strainMatrix[6][nUdof], eqMatrix[3][6];
    real64 matBD[nUdof][6], matED[3][6];

    int Heaviside[ numNodesPerElem ];

    // TODO: asking for the stiffness here will only work for elastic models.  most other models
    //       need to know the strain increment to compute the current stiffness value.

    m_constitutiveUpdate.getElasticStiffness( k, stack.constitutiveStiffness );

    SolidMechanicsEFEMKernelsHelper::computeHeavisideFunction< numNodesPerElem >( Heaviside,
                                                                                  stack.X,
                                                                                  m_nVec[embSurfIndex],
                                                                                  m_surfaceCenter[embSurfIndex] );


    SolidMechanicsEFEMKernelsHelper::assembleEquilibriumOperator( eqMatrix,
                                                                  m_nVec[embSurfIndex],
                                                                  m_tVec1[embSurfIndex],
                                                                  m_tVec2[embSurfIndex],
                                                                  stack.hInv );

    SolidMechanicsEFEMKernelsHelper::assembleCompatibilityOperator< numNodesPerElem >( compMatrix,
                                                                                       m_nVec[embSurfIndex],
                                                                                       m_tVec1[embSurfIndex],
                                                                                       m_tVec2[embSurfIndex],
                                                                                       Heaviside,
                                                                                       dNdX );

    SolidMechanicsEFEMKernelsHelper::assembleStrainOperator< 6, nUdof, numNodesPerElem >( strainMatrix, dNdX );

    // transp(B)D
    LvArray::tensorOps::Rij_eq_AkiBkj< nUdof, 6, 6 >( matBD, strainMatrix, stack.constitutiveStiffness );
    // ED
    LvArray::tensorOps::Rij_eq_AikBkj< 3, 6, 6 >( matED, eqMatrix, stack.constitutiveStiffness );
    // EDC
    LvArray::tensorOps::Rij_eq_AikBjk< 3, 3, 6 >( Kww_gauss, matED, compMatrix );
    // EDB
    LvArray::tensorOps::Rij_eq_AikBkj< 3, nUdof, 6 >( Kwu_gauss, matED, strainMatrix );
    // transp(B)DB
    LvArray::tensorOps::Rij_eq_AikBjk< nUdof, 3, 6 >( Kuw_gauss, matBD, compMatrix );

    // multiply by determinant and add to element matrix
    LvArray::tensorOps::scaledAdd< 3, 3 >( stack.localKww, Kww_gauss, -detJ );
    LvArray::tensorOps::scaledAdd< 3, nUdof >( stack.localKwu, Kwu_gauss, -detJ );
    LvArray::tensorOps::scaledAdd< nUdof, 3 >( stack.localKuw, Kuw_gauss, -detJ );


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

    // Add traction contribution tranction
    LvArray::tensorOps::scaledAdd< 3 >( stack.localRw, stack.tractionVec, -1 );
    LvArray::tensorOps::scaledAdd< 3, 3 >( stack.localKww, stack.dTractiondw, -1 );

    for( localIndex i = 0; i < nUdof; ++i )
    {
      localIndex const dof = LvArray::integerConversion< localIndex >( stack.dispEqnRowIndices[ i ] );
      if( dof < 0 || dof >= m_matrix.numRows() ) continue;

      RAJA::atomicAdd< parallelDeviceAtomic >( &m_rhs[dof], stack.localRu[i] );

      m_matrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( dof,
                                                                              stack.jumpColIndices,
                                                                              stack.localKuw[i],
                                                                              3 );

    }

    for( localIndex i=0; i < 3; ++i )
    {
      localIndex const dof = LvArray::integerConversion< localIndex >( stack.jumpEqnRowIndices[ i ] );

      if( dof < 0 || dof >= m_matrix.numRows() ) continue;

      RAJA::atomicAdd< parallelDeviceAtomic >( &m_rhs[dof], stack.localRw[i] );

      // fill in matrix
      m_matrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( dof,
                                                                              stack.jumpColIndices,
                                                                              stack.localKww[i],
                                                                              3 );

      m_matrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( dof,
                                                                              stack.dispColIndices,
                                                                              stack.localKwu[i],
                                                                              numNodesPerElem*3 );
    }


    return maxForce;
  }



protected:
  arrayView1d< globalIndex const > const m_wDofNumber;

  arrayView2d< real64 const > const m_w;

  arrayView2d< real64 const > const m_tractionVec;

  arrayView3d< real64 const > const m_dTraction_dJump;

  arrayView2d< real64 const > const m_nVec;

  arrayView2d< real64 const > const m_tVec1;

  arrayView2d< real64 const > const m_tVec2;

  arrayView2d< real64 const > const m_surfaceCenter;

  arrayView1d< real64 const > const m_surfaceArea;

  arrayView1d< real64 const > const m_elementVolume;

  SortedArrayView< localIndex const > const m_fracturedElems;

  ArrayOfArraysView< localIndex const > const m_cellsToEmbeddedSurfaces;

};


} // namespace PoroelasticKernels

} // namespace geosx

#include "finiteElement/kernelInterface/SparsityKernelBase.hpp"

#endif // GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASEPOROELASTICKERNEL_HPP_
