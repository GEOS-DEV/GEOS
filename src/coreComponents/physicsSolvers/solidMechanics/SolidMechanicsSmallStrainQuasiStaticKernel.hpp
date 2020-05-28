/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file SolidMechanicsSmallStrainQuasiStaticKernels.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSSMALLSTRAINQUASISTATIC_HPP_
#define GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSSMALLSTRAINQUASISTATIC_HPP_

#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/ConstitutiveBase.hpp"
#include "finiteElement/ElementLibrary/FiniteElementBase.h"
#include "finiteElement/FiniteElementShapeFunctionKernel.hpp"
#include "finiteElement/Kinematics.h"
#include "finiteElement/kernelInterface/ImplicitKernelBase.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"
#include "TimeIntegrationOption.hpp"

namespace geosx
{

namespace SolidMechanicsLagrangianFEMKernels
{


/**
 * @struct QuasiStaticConstructorParams
 * @copydoc geosx::finiteElement::ImplicitKernelBaseConstructorParams
 */
struct QuasiStaticConstructorParams : finiteElement::ImplicitKernelBaseConstructorParams
{
  /**
   * @brief Constructor
   * @param inputDofNumber The dof number for the primary field.
   * @param inputMatrix Reference to the Jacobian matrix.
   * @param inputRhs Reference to the RHS vector.
   * @param inputGravityVector The gravity vector.
   */
  QuasiStaticConstructorParams( arrayView1d< globalIndex const > const & inputDofNumber,
                                ParallelMatrix & inputMatrix,
                                ParallelVector & inputRhs,
                                real64 const inputGravityVector[3] ):
    finiteElement::ImplicitKernelBaseConstructorParams( inputDofNumber,
                             inputMatrix,
                             inputRhs ),
    m_gravityVector{ inputGravityVector[0], inputGravityVector[1], inputGravityVector[2] }
  {}
  real64 const m_gravityVector[3];
};

/**
 * @brief Implements kernels for solving quasi-static equilibrium.
 * @copydoc geosx::finiteElement::KernelBase
 * @tparam NUM_NODES_PER_ELEM The number of nodes per element for the
 *                            @p SUBREGION_TYPE.
 * @tparam UNUSED An unused parameter since we are assuming that the test and
 *                trial space have the same number of support points.
 *
 * ### QuasiStatic Description
 * Implements the KernelBase interface functions required for solving the
 * quasi-static equilibrium equations using on of the
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
          int NUM_NODES_PER_ELEM,
          int >
class QuasiStatic :
  public finiteElement::ImplicitKernelBase< SUBREGION_TYPE,
                                            CONSTITUTIVE_TYPE,
                                            NUM_NODES_PER_ELEM,
                                            NUM_NODES_PER_ELEM,
                                            3,
                                            3 >
{
public:
  using Base = finiteElement::ImplicitKernelBase< SUBREGION_TYPE,
                                                  CONSTITUTIVE_TYPE,
                                                  NUM_NODES_PER_ELEM,
                                                  NUM_NODES_PER_ELEM,
                                                  3,
                                                  3 >;

  /// Number of nodes per element...which is equal to the
  /// numTestSupportPointPerElem and numTrialSupportPointPerElem by definition.
  static constexpr int numNodesPerElem = NUM_NODES_PER_ELEM;

  /// @copydoc geosx::finiteElement::ImplicitKernelBase::numDofPerTestSupportPoint
  using Base::numDofPerTestSupportPoint;

  /// @copydoc geosx::finiteElement::ImplicitKernelBase::numDofPerTrialSupportPoint
  using Base::numDofPerTrialSupportPoint;


  /// @copydoc geosx::finiteElement::ImplicitKernelBase::m_dofNumber
  using Base::m_dofNumber;

  /// @copydoc geosx::finiteElement::ImplicitKernelBase::m_matrix
  using Base::m_matrix;

  /// @copydoc geosx::finiteElement::ImplicitKernelBase::m_rhs
  using Base::m_rhs;

  /// @copydoc geosx::finiteElement::ImplicitKernelBase::elemsToNodes
  using Base::elemsToNodes;

  /// @copydoc geosx::finiteElement::ImplicitKernelBase::elemGhostRank
  using Base::elemGhostRank;

  /// @copydoc geosx::finiteElement::ImplicitKernelBase::constitutiveUpdate
  using Base::constitutiveUpdate;

  /// @copydoc geosx::finiteElement::ImplicitKernelBase::m_finiteElementSpace
  using Base::m_finiteElementSpace;

  /// @copydoc geosx::finiteElement::ImplicitKernelBase::Launch
  using Base::Launch;

  /// Alias for the struct that holds the constructor parameters
  using ConstructorParams = QuasiStaticConstructorParams;


  /**
   * @brief Constructor
   * @copydoc geosx::finiteElement::ImplicitKernelBase::ImplicitKernelBase
   * @param inputGravityVector The gravity vector.
   */
  QuasiStatic( NodeManager const & nodeManager,
               EdgeManager const & edgeManager,
               FaceManager const & faceManager,
               SUBREGION_TYPE const & elementSubRegion,
               FiniteElementBase const * const finiteElementSpace,
               CONSTITUTIVE_TYPE * const inputConstitutiveType,
               arrayView1d< globalIndex const > const & inputDofNumber,
               ParallelMatrix & inputMatrix,
               ParallelVector & inputRhs,
               real64 const inputGravityVector[3] ):
    Base( nodeManager,
          edgeManager,
          faceManager,
          elementSubRegion,
          finiteElementSpace,
          inputConstitutiveType,
          inputDofNumber,
          inputMatrix,
          inputRhs ),
    m_disp( nodeManager.totalDisplacement()),
    m_uhat( nodeManager.incrementalDisplacement()),
    dNdX( elementSubRegion.template getReference< array3d< R1Tensor > >( dataRepository::keys::dNdX )),
    detJ( elementSubRegion.template getReference< array2d< real64 > >( dataRepository::keys::detJ ) ),
    m_gravityVector{ inputGravityVector[0], inputGravityVector[1], inputGravityVector[2] }
  {}

  /**
   * @copydoc QuasiStatic
   * @brief Constructor that holds variadic components in a struct.
   * @param params Hold the variadic components of primary constructor.
   */
  QuasiStatic( NodeManager const & nodeManager,
               EdgeManager const & edgeManager,
               FaceManager const & faceManager,
               SUBREGION_TYPE const & elementSubRegion,
               FiniteElementBase const * const finiteElementSpace,
               CONSTITUTIVE_TYPE * const inputConstitutiveType,
               ConstructorParams & params ):
    QuasiStatic( nodeManager,
                 edgeManager,
                 faceManager,
                 elementSubRegion,
                 finiteElementSpace,
                 inputConstitutiveType,
                 params.m_dofNumber,
                 params.m_matrix,
                 params.m_rhs,
                 params.m_gravityVector )
  {}

  /**
   * @brief Copy global values from primary field to a local stack array.
   * @copydoc geosx::finiteElement::ImplicitKernelBase::setup.
   *
   * For the QuasiStatic implementation, global values from the displacement,
   * incremental displacement, and degree of freedom numbers are placed into
   * element local stack storage.
   */

  template< typename STACK_VARIABLE_TYPE >
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void setup( localIndex const k,
              STACK_VARIABLE_TYPE & stack ) const
  {
    for( localIndex a=0; a<NUM_NODES_PER_ELEM; ++a )
    {
      localIndex const localNodeIndex = elemsToNodes( k, a );

      stack.u_local[ a ] = m_disp[ localNodeIndex ];
      stack.uhat_local[ a ] = m_uhat[ localNodeIndex ];

      for( int i=0; i<3; ++i )
      {
        stack.localRowDofIndex[a*3+i] = m_dofNumber[localNodeIndex]+i;
        stack.localColDofIndex[a*3+i] = m_dofNumber[localNodeIndex]+i;
      }
    }

  }

  /**
   * @copydoc geosx::finiteElement::KernelBase::quadraturePointStateUpdate.
   *
   * For solid mechanics kernels, the strain increment is calculated, and the
   * constitutive update is called. In addition, the constitutive stiffness
   * stack variable is filled by the constitutive model.
   */
  template< typename STACK_VARIABLE_TYPE >
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void quadraturePointStateUpdate( localIndex const k,
                                   localIndex const q,
                                   STACK_VARIABLE_TYPE & stack ) const
  {
    real64 strainInc[6] = {0};
    for( localIndex a = 0; a < NUM_NODES_PER_ELEM; ++a )
    {
      strainInc[0] = strainInc[0] + dNdX( k, q, a )[0] * stack.uhat_local[a][0];
      strainInc[1] = strainInc[1] + dNdX( k, q, a )[1] * stack.uhat_local[a][1];
      strainInc[2] = strainInc[2] + dNdX( k, q, a )[2] * stack.uhat_local[a][2];
      strainInc[3] = strainInc[3] + dNdX( k, q, a )[2] * stack.uhat_local[a][1] +
                     dNdX( k, q, a )[1] * stack.uhat_local[a][2];
      strainInc[4] = strainInc[4] + dNdX( k, q, a )[2] * stack.uhat_local[a][0] +
                     dNdX( k, q, a )[0] * stack.uhat_local[a][2];
      strainInc[5] = strainInc[5] + dNdX( k, q, a )[1] * stack.uhat_local[a][0] +
                     dNdX( k, q, a )[0] * stack.uhat_local[a][1];
    }

    constitutiveUpdate.SmallStrain( k, q, strainInc );

    GEOSX_UNUSED_VAR( q )
    constitutiveUpdate.GetStiffness( k, stack.constitutiveStiffness );
  }


  /**
   * @copydoc geosx::finiteElement::KernelBase::quadraturePointJacobianContribution.
   *
   * For solid mechanics kernels, the derivative of the force residual wrt
   * the incremental displacement is filled into the local element jacobian.
   */
  template< typename STACK_VARIABLE_TYPE,
            typename DYNAMICS_LAMBDA = STD_FUNCTION< void( localIndex, localIndex) > >
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void quadraturePointJacobianContribution( localIndex const k,
                                            localIndex const q,
                                            STACK_VARIABLE_TYPE & stack,
                                            DYNAMICS_LAMBDA && dynamicsTerms = [] GEOSX_HOST_DEVICE ( localIndex,
                                                                                                 localIndex){} ) const
  {
    for( localIndex a=0; a<NUM_NODES_PER_ELEM; ++a )
    {
      for( localIndex b=0; b<NUM_NODES_PER_ELEM; ++b )
      {
        real64 const (&c)[6][6] = stack.constitutiveStiffness;
        stack.localJacobian[ a*3+0 ][ b*3+0 ] -= ( c[0][0]*dNdX( k, q, a )[0]*dNdX( k, q, b )[0] +
                                                   c[5][5]*dNdX( k, q, a )[1]*dNdX( k, q, b )[1] +
                                                   c[4][4]*dNdX( k, q, a )[2]*dNdX( k, q, b )[2] ) * detJ( k, q );

        stack.localJacobian[ a*3+0 ][ b*3+1 ] -= ( c[5][5]*dNdX( k, q, a )[1]*dNdX( k, q, b )[0] +
                                                   c[0][1]*dNdX( k, q, a )[0]*dNdX( k, q, b )[1] ) * detJ( k, q );

        stack.localJacobian[ a*3+0 ][ b*3+2 ] -= ( c[4][4]*dNdX( k, q, a )[2]*dNdX( k, q, b )[0] +
                                                   c[0][2]*dNdX( k, q, a )[0]*dNdX( k, q, b )[2] ) * detJ( k, q );

        stack.localJacobian[ a*3+1 ][ b*3+1 ] -= ( c[5][5]*dNdX( k, q, a )[0]*dNdX( k, q, b )[0] +
                                                   c[1][1]*dNdX( k, q, a )[1]*dNdX( k, q, b )[1] +
                                                   c[3][3]*dNdX( k, q, a )[2]*dNdX( k, q, b )[2] ) * detJ( k, q );

        stack.localJacobian[ a*3+1 ][ b*3+0 ] -= ( c[0][1]*dNdX( k, q, a )[1]*dNdX( k, q, b )[0] +
                                                   c[5][5]*dNdX( k, q, a )[0]*dNdX( k, q, b )[1] ) * detJ( k, q );

        stack.localJacobian[ a*3+1 ][ b*3+2 ] -= ( c[3][3]*dNdX( k, q, a )[2]*dNdX( k, q, b )[1] +
                                                   c[1][2]*dNdX( k, q, a )[1]*dNdX( k, q, b )[2] ) * detJ( k, q );

        stack.localJacobian[ a*3+2 ][ b*3+0 ] -= ( c[0][2]*dNdX( k, q, a )[2]*dNdX( k, q, b )[0] +
                                                   c[4][4]*dNdX( k, q, a )[0]*dNdX( k, q, b )[2] ) * detJ( k, q );

        stack.localJacobian[ a*3+2 ][ b*3+1 ] -= ( c[1][2]*dNdX( k, q, a )[2]*dNdX( k, q, b )[1] +
                                                   c[3][3]*dNdX( k, q, a )[1]*dNdX( k, q, b )[2] ) * detJ( k, q );

        stack.localJacobian[ a*3+2 ][ b*3+2 ] -= ( c[4][4]*dNdX( k, q, a )[0]*dNdX( k, q, b )[0] +
                                                   c[3][3]*dNdX( k, q, a )[1]*dNdX( k, q, b )[1] +
                                                   c[2][2]*dNdX( k, q, a )[2]*dNdX( k, q, b )[2] ) * detJ( k, q );

        dynamicsTerms( a, b );
      }
    }
  }

  /**
   * @copydoc geosx::finiteElement::KernelBase::quadraturePointResidualContribution.
   *
   * The divergence of the stress is integrated over the volume of the element,
   * yielding the nodal force (residual) contributions.
   */
  template< typename STACK_VARIABLE_TYPE,
            typename STRESS_MODIFIER = STD_FUNCTION< void( real64 * ) > >
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void quadraturePointResidualContribution( localIndex const k,
                                            localIndex const q,
                                            STACK_VARIABLE_TYPE & stack,
                                            STRESS_MODIFIER && stressModifier = [] GEOSX_DEVICE ( real64 * ) {} ) const
  {
    real64 stress[6] = { constitutiveUpdate.m_stress( k, q, 0 ),
                         constitutiveUpdate.m_stress( k, q, 1 ),
                         constitutiveUpdate.m_stress( k, q, 2 ),
                         constitutiveUpdate.m_stress( k, q, 3 ),
                         constitutiveUpdate.m_stress( k, q, 4 ),
                         constitutiveUpdate.m_stress( k, q, 5 ) };

    stressModifier( stress );

    for( localIndex a = 0; a < NUM_NODES_PER_ELEM; ++a )
    {
      stack.localResidual[ a * 3 + 0 ] -= ( stress[ 0 ] * dNdX( k, q, a )[ 0 ] +
                                            stress[ 5 ] * dNdX( k, q, a )[ 1 ] +
                                            stress[ 4 ] * dNdX( k, q, a )[ 2 ] -
                                            m_gravityVector[0] ) * detJ( k, q );
      stack.localResidual[ a * 3 + 1 ] -= ( stress[ 5 ] * dNdX( k, q, a )[ 0 ] +
                                            stress[ 1 ] * dNdX( k, q, a )[ 1 ] +
                                            stress[ 3 ] * dNdX( k, q, a )[ 2 ] -
                                            m_gravityVector[1] ) * detJ( k, q );
      stack.localResidual[ a * 3 + 2 ] -= ( stress[ 4 ] * dNdX( k, q, a )[ 0 ] +
                                            stress[ 3 ] * dNdX( k, q, a )[ 1 ] +
                                            stress[ 2 ] * dNdX( k, q, a )[ 2 ] -
                                            m_gravityVector[2] ) * detJ( k, q );
    }
  }

  /**
   * @copydoc geosx::finiteElement::ImplicitKernelBase::complete.
   */
  template< typename STACK_VARIABLE_TYPE >
  //GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  real64 complete( localIndex const GEOSX_UNUSED_PARAM( k ),
                   STACK_VARIABLE_TYPE & stack ) const
  {
    real64 meanForce = 0;
    for( localIndex a=0; a<stack.numRows; ++a )
    {
//        RAJA::atomicMax< RAJA::auto_atomic >( &meanForce, stack.localResidual[a] );
      meanForce = std::max( meanForce, stack.localResidual[a] );
//                meanForce += fabs( stack.localResidual[a] );
    }
//            meanForce /= stack.ndof;

    m_matrix.add( stack.localRowDofIndex,
                  stack.localColDofIndex,
                  &(stack.localJacobian[0][0]),
                  stack.numRows,
                  stack.numCols );

    m_rhs.add( stack.localRowDofIndex,
               stack.localResidual,
               stack.numRows );

    return meanForce;
  }

  //*****************************************************************************
  /**
   * @class StackVariables
   * @copydoc geosx::finiteElement::ImplicitKernelBase::StackVariables
   *
   * Adds a stack array for the displacement, incremental displacement, and the
   * constitutive stiffness.
   */
  struct StackVariables : Base::StackVariables
  {
public:

    /// Constructor.
    GEOSX_HOST_DEVICE
    StackVariables():
      Base::StackVariables(),
            u_local(),
            uhat_local(),
            constitutiveStiffness{ {0.0} }
    {}

    /// Stack storage for the element local nodal displacement
    R1Tensor u_local[numNodesPerElem];

    /// Stack storage for the element local nodal incremental displacement
    R1Tensor uhat_local[numNodesPerElem];

    /// Stack storage for the constitutive stiffness at a quadrature point.
    real64 constitutiveStiffness[6][6];
  };
  //*****************************************************************************


protected:
  /// The rank-global displacement array.
  arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const m_disp;

  /// The rank-global incremental displacement array.
  arrayView2d< real64 const, nodes::INCR_DISPLACEMENT_USD > const m_uhat;

  /// The shape function derivative for each quadrature point.
  arrayView3d< R1Tensor const > const dNdX;

  /// The parent->physical jacobian determinant for each quadrature point.
  arrayView2d< real64 const > const detJ;

  /// The gravity vector.
  real64 const m_gravityVector[3];


};


} // namespace SolidMechanicsLagrangianFEMKernels

} // namespace geosx

#endif // GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSSMALLSTRAINQUASISTATIC_HPP_
