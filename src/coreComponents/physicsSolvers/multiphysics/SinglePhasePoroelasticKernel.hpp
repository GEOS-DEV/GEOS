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

#ifndef GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASEPOROELASTICKERNEL_HPP_
#define GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASEPOROELASTICKERNEL_HPP_
#include "finiteElement/kernelInterface/ImplicitKernelBase.hpp"

namespace geosx
{

namespace PoroelasticKernels
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
 * Implements the KernelBase interface functions required for solving the
 * quasi-static equilibrium equations using one of the
 * "finite element kernel application" functions such as
 * geosx::finiteElement::RegionBasedKernelApplication.
 *
 * In this implementation, the template parameter @p NUM_NODES_PER_ELEM is used
 * in place of both @p NUM_TEST_SUPPORT_POINTS_PER_ELEM and
 * @p NUM_TRIAL_SUPPORT_POINTS_PER_ELEM, which are assumed to be equal. This
 * results in the @p UNUSED template parameter as only the NUM_NODES_PER_ELEM
 * is passed to the ImplicitKernelBase template to form the base class.
 *
 * Additionally, the number of degrees of freedom per support point for both
 * the test and trial spaces are specified as `3` when specifying the base
 * class.
 */
template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
class SinglePhase :
  public finiteElement::ImplicitKernelBase< SUBREGION_TYPE,
                                            CONSTITUTIVE_TYPE,
                                            FE_TYPE,
                                            3,
                                            3 >
{
public:
  /// Alias for the base class;
  using Base = finiteElement::ImplicitKernelBase< SUBREGION_TYPE,
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
  SinglePhase( NodeManager const & nodeManager,
               EdgeManager const & edgeManager,
               FaceManager const & faceManager,
               SUBREGION_TYPE const & elementSubRegion,
               FE_TYPE const & finiteElementSpace,
               CONSTITUTIVE_TYPE * const inputConstitutiveType,
               arrayView1d< globalIndex const > const & inputDispDofNumber,
               string const & inputFlowDofKey,
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
          inputDispDofNumber,
          rankOffset,
          inputMatrix,
          inputRhs ),
    m_X( nodeManager.referencePosition()),
    m_disp( nodeManager.totalDisplacement()),
    m_uhat( nodeManager.incrementalDisplacement()),
    m_gravityVector{ inputGravityVector[0], inputGravityVector[1], inputGravityVector[2] },
    m_solidDensity( inputConstitutiveType->getDensity() ),
    m_fluidDensity( inputConstitutiveType->getDensity() ), //*** TODO THIS IS INCORRECT !!!!!!
    m_flowDofNumber(elementSubRegion.template getReference< array1d< globalIndex > >( inputFlowDofKey )),
    m_fluidPressure( elementSubRegion.template getReference< array1d< real64 > >( "pressure" ) ),
    m_deltaFluidPressure( elementSubRegion.template getReference< array1d< real64 > >( "deltaPressure" ) )
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

    /// Constructor.
    GEOSX_HOST_DEVICE
    StackVariables():
      Base::StackVariables(),
                                       xLocal(),
                                       u_local(),
                                       uhat_local()
    {}

#if !defined(CALC_FEM_SHAPE_IN_KERNEL)
    /// Dummy
    int xLocal;
#else
    /// C-array stack storage for element local the nodal positions.
    real64 xLocal[ numNodesPerElem ][ 3 ];
#endif

    /// Stack storage for the element local nodal displacement
    real64 u_local[numNodesPerElem][numDofPerTrialSupportPoint];

    /// Stack storage for the element local nodal incremental displacement
    real64 uhat_local[numNodesPerElem][numDofPerTrialSupportPoint];

    real64 localFlowResidual[1];
    real64 localDispFlowJacobian[numDispDofPerElem][1];
    real64 localFlowDispJacobian[1][numDispDofPerElem];

    /// C-array storage for the element local row degrees of freedom.
    globalIndex localFlowDofIndex[1];

  };
  //*****************************************************************************

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
        stack.xLocal[ a ][ i ] = m_X[ localNodeIndex ][ i ];
#endif
        stack.u_local[ a ][i] = m_disp[ localNodeIndex ][i];
        stack.uhat_local[ a ][i] = m_uhat[ localNodeIndex ][i];
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

    constexpr real64 biotCoefficient = 1.0;
    real64 dNdX[ numNodesPerElem ][ 3 ];
    real64 const detJ = m_finiteElementSpace.template getGradN< FE_TYPE >( k, q, stack.xLocal, dNdX );

    real64 strainInc[6] = {0};
    FE_TYPE::symmetricGradient( dNdX, stack.uhat_local, strainInc );

    m_constitutiveUpdate.smallStrain( k, q, strainInc );

    typename CONSTITUTIVE_TYPE::KernelWrapper::DiscretizationOps stiffnessHelper;
    m_constitutiveUpdate.setDiscretizationOps( k, q, stiffnessHelper );

    stiffnessHelper.template upperBTDB< numNodesPerElem >( dNdX, -detJ, stack.localJacobian );

    real64 stress[6];

    m_constitutiveUpdate.getStress( k, q, stress );

    real64 const pressure = biotCoefficient * ( m_fluidPressure[k] + m_deltaFluidPressure[k] );
    stress[0] -= pressure;
    stress[1] -= pressure;
    stress[2] -= pressure;


    real64 const gravityForce[3] = { m_gravityVector[0] * m_solidDensity( k, q )* detJ,
                                     m_gravityVector[1] * m_solidDensity( k, q )* detJ,
                                     m_gravityVector[2] * m_solidDensity( k, q )* detJ };

    for( localIndex i=0; i<6; ++i )
    {
      stress[i] *= -detJ;
    }

    real64 N[numNodesPerElem];
    FE_TYPE::calcN( q, N );
    FE_TYPE::plusGradNajAijPlusNaFi( dNdX,
                                     stress,
                                     N,
                                     gravityForce,
                                     reinterpret_cast< real64 (&)[numNodesPerElem][3] >(stack.localResidual) );


    for( integer a = 0; a < numNodesPerElem; ++a )
    {
      stack.localDispFlowJacobian[ a * 3 + 0][0] += biotCoefficient * dNdX[a][0] * detJ;
      stack.localDispFlowJacobian[ a * 3 + 1][0] += biotCoefficient * dNdX[a][1] * detJ;
      stack.localDispFlowJacobian[ a * 3 + 2][0] += biotCoefficient * dNdX[a][2] * detJ;
      stack.localFlowDispJacobian[ 0][a * 3 + 0] += m_fluidDensity( k, q ) * biotCoefficient * dNdX[a][0] * detJ;
      stack.localFlowDispJacobian[ 0][a * 3 + 1] += m_fluidDensity( k, q ) * biotCoefficient * dNdX[a][1] * detJ;
      stack.localFlowDispJacobian[ 0][a * 3 + 2] += m_fluidDensity( k, q ) * biotCoefficient * dNdX[a][2] * detJ;


      real64 Rf_tmp =   dNdX[a][0] * stack.uhat_local[a][0]
                      + dNdX[a][1] * stack.uhat_local[a][1]
                      + dNdX[a][2] * stack.uhat_local[a][2];
      Rf_tmp *= m_fluidDensity( k, q ) * biotCoefficient * detJ;
      stack.localFlowResidual[0] += Rf_tmp;
    }

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

    CONSTITUTIVE_TYPE::KernelWrapper::DiscretizationOps::template fillLowerBTDB< numNodesPerElem >( stack.localJacobian );

    constexpr int nPDof = 1;
    constexpr int nUDof = numNodesPerElem * numDofPerTestSupportPoint;

    for( int localNode = 0; localNode < numNodesPerElem; ++localNode )
    {
      for( int dim = 0; dim < numDofPerTestSupportPoint; ++dim )
      {
        localIndex const dof = LvArray::integerConversion< localIndex >( stack.localRowDofIndex[ numDofPerTestSupportPoint * localNode + dim ] - m_dofRankOffset );
        if( dof < 0 || dof >= m_matrix.numRows() ) continue;
        m_matrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( dof,
                                                                                stack.localRowDofIndex,
                                                                                stack.localJacobian[ numDofPerTestSupportPoint * localNode + dim ],
                                                                                numNodesPerElem * numDofPerTrialSupportPoint );

        RAJA::atomicAdd< parallelDeviceAtomic >( &m_rhs[ dof ], stack.localResidual[ numDofPerTestSupportPoint * localNode + dim ] );
        maxForce = fmax( maxForce, fabs( stack.localResidual[ numDofPerTestSupportPoint * localNode + dim ] ) );

        m_matrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( dof,
                                                                                stack.localFlowDofIndex,
                                                                                stack.localDispFlowJacobian[numDofPerTestSupportPoint * localNode + dim],
                                                                                nPDof );

      }
    }

    for( localIndex i = 0; i < nPDof; ++i )
    {
      localIndex const dof = LvArray::integerConversion< localIndex >( stack.localFlowDofIndex[ i ] - m_dofRankOffset );
      if( dof < 0 || dof >= m_matrix.numRows() )
        continue;
      m_matrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( dof,
                                                                              stack.localRowDofIndex,
                                                                              stack.localFlowDispJacobian[i],
                                                                              nUDof );

      RAJA::atomicAdd< parallelDeviceAtomic >( &m_rhs[ dof ], stack.localFlowResidual[i] );
    }



    return maxForce;
  }



protected:
  /// The array containing the nodal position array.
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const m_X;

  /// The rank-global displacement array.
  arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const m_disp;

  /// The rank-global incremental displacement array.
  arrayView2d< real64 const, nodes::INCR_DISPLACEMENT_USD > const m_uhat;

  /// The gravity vector.
  real64 const m_gravityVector[3];

  /// The rank global density
  arrayView2d< real64 const > const m_solidDensity;
  arrayView2d< real64 const > const m_fluidDensity;

  /// The global degree of freedom number
  arrayView1d< globalIndex const > const m_flowDofNumber;

  /// The rank-global fluid pressure array.
  arrayView1d< real64 const > const m_fluidPressure;

  /// The rank-global delta-fluid pressure array.
  arrayView1d< real64 const > const m_deltaFluidPressure;


};


} // namespace SolidMechanicsLagrangianFEMKernels

} // namespace geosx

#include "finiteElement/kernelInterface/SparsityKernelBase.hpp"

#endif // GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASEPOROELASTICKERNEL_HPP_
