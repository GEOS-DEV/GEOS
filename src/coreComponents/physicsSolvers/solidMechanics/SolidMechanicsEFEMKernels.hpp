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
                                                          FE_TYPE,
                                                          3,
                                                          3 >
{
public:
  /// Alias for the base class;
  using Base = SolidMechanicsLagrangianFEMKernels::QuasiStatic< SUBREGION_TYPE,
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
    m_wDofNumber( wDofNumber )
    {}

  //***************************************************************************
    /**
     * @copydoc finiteElement::KernelBase::StackVariables
     */
    struct StackVariables
    {
      /// The number of displacement dofs per element.
      static constexpr int numUdofs = numTestSupportPointsPerElem * numDofPerTestSupportPoint;

      /// The number of jump dofs per element.
      static constexpr int numWdofs = 3;

      /**
       * Default constructor
       */
      GEOSX_HOST_DEVICE
      StackVariables():
        localRowDofIndex{ 0 },
        localColDofIndex{ 0 },
        localResidual{ 0.0 },
        localJacobian{ {0.0} }
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
    };
    //***************************************************************************

  /**
   * @brief Copy global values from primary field to a local stack array.
   * @copydoc ::geosx::finiteElement::ImplicitKernelBase::setup
   *
   * For the QuasiStatic implementation, global values from the displacement,
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
        stack.dispEqnRowIndices[a*3+i] = m_dofNumber[localNodeIndex]+i-m_dofRankOffset;
        stack.dispColIndices[a*3+i]    = m_dofNumber[localNodeIndex]+i;
      }
    }

    for (int k=0; k<numFracturedElems; ++k)
    {
    	for (int i=0; i<3; ++i)
    	{
    		// need to grab the index.
    		stack.jumpEqnRowIndices = m_wDofNumber[k] + i - m_dofRankOffset;
    		stack.jumpColIndices    = m_wDofNumber[k] + i;
    	}
    }
  }


  /**
    * @copydoc geosx::finiteElement::KernelBase::quadraturePointKernel
    * @tparam STRESS_MODIFIER Type of optional functor to allow for the
    * modification of stress prior to integration.
    * @param stressModifier An optional functor to allow for the modification
    *  of stress prior to integration.
    * For solid mechanics kernels, the strain increment is calculated, and the
    * constitutive update is called. In addition, the constitutive stiffness
    * stack variable is filled by the constitutive model.
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

     real64 strainInc[6] = {0};
     FE_TYPE::symmetricGradient( dNdX, stack.uhat_local, strainInc );

     m_constitutiveUpdate.SmallStrain( k, q, strainInc );

     typename CONSTITUTIVE_TYPE::KernelWrapper::DiscretizationOps stiffnessHelper;
     m_constitutiveUpdate.setDiscretizationOps( k, q, stiffnessHelper );

     stiffnessHelper.template upperBTDB< numNodesPerElem >( dNdX, -detJ, stack.localJacobian );

     real64 stress[6];

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

     for( localIndex i = 0; i < numNodesPerElem*3; ++i )
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
};

} // namespace SolidMechanicsEFEMKernels

} // namespace geosx


#endif /* GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSEFEMKERNELS_HPP_ */
