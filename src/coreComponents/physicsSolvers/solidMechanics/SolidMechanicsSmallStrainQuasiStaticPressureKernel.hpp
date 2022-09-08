/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file SolidMechanicsSmallStrainQuasiStaticPressureKernel.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSSMALLSTRAINQUASIPRESSURESTATIC_HPP_
#define GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSSMALLSTRAINQUASIPRESSURESTATIC_HPP_

#include "finiteElement/BilinearFormUtilities.hpp"
#include "finiteElement/LinearFormUtilities.hpp"
#include "finiteElement/kernelInterface/ImplicitKernelBase.hpp"

namespace geosx
{

namespace solidMechanicsLagrangianFEMKernels
{

/**
 * @brief Implements kernels for solving quasi-static equilibrium with background pressure effects, including fractures.
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
 * geosx::finiteElement::RegionBasedKernelApplication. Pressure effects in
 * the rock and its fractures are accounted for in the fracturePressure and 
 * matrixPressure terms. Currently (Sept22), only use case for this kernel is 
 * the multi-resolution phase-field method.
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
class QuasiStaticPressure :
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

  /// Maximum number of nodes per element, which is equal to the maxNumTestSupportPointPerElem and
  /// maxNumTrialSupportPointPerElem by definition. When the FE_TYPE is not a Virtual Element, this
  /// will be the actual number of nodes per element.
  static constexpr int numNodesPerElem = Base::maxNumTestSupportPointsPerElem;
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
  QuasiStaticPressure( NodeManager const & nodeManager,
                       EdgeManager const & edgeManager,
                       FaceManager const & faceManager,
                       localIndex const targetRegionIndex,
                       SUBREGION_TYPE const & elementSubRegion,
                       FE_TYPE const & finiteElementSpace,
                       CONSTITUTIVE_TYPE & inputConstitutiveType,
                       arrayView1d< globalIndex const > const inputDispDofNumber,
                       globalIndex const rankOffset,
                       CRSMatrixView< real64, globalIndex const > const inputMatrix,
                       arrayView1d< real64 > const inputRhs,
                       real64 const (&inputGravityVector)[3] ):
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
          inputRhs ),
    m_X( nodeManager.referencePosition()),
    m_disp( nodeManager.totalDisplacement()),
    m_uhat( nodeManager.incrementalDisplacement()),
    m_gravityVector{ inputGravityVector[0], inputGravityVector[1], inputGravityVector[2] },
    m_gravityAcceleration( LvArray::tensorOps::l2Norm< 3 >( inputGravityVector ) ),
    m_solidDensity( inputConstitutiveType.getDensity() ),
    m_pressureMatrix( elementSubRegion.template getReference< array1d< real64 > >( "hardCodedPMatrixName" ) ),
    m_pressureFracture( elementSubRegion.template getReference< array1d< real64 > >( "hardCodedPFractureName" ) )
    //m_pressureMatrix( elementSubRegion.template getExtrinsicData< extrinsicMeshData::flow::matrixPressure >() ),
    //m_pressureFracture( elementSubRegion.template getExtrinsicData< extrinsicMeshData::flow::fracturePressure >() )
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

    static constexpr int numDispDofPerElem =  Base::StackVariables::maxNumRows;

    /// Constructor.
    GEOSX_HOST_DEVICE
    StackVariables():
      Base::StackVariables(),
            xLocal(),
            u_local(),
            uhat_local(),
            localResidualMomentum( Base::StackVariables::localResidual ),
            dLocalResidualMomentum_dDisplacement( Base::StackVariables::localJacobian )
    {}

#if !defined(CALC_FEM_SHAPE_IN_KERNEL)
    /// Dummy
    int xLocal;
#else
    /// C-array stack storage for element local the nodal positions.
    real64 xLocal[numNodesPerElem][3];
#endif

    /// Stack storage for the element local nodal displacement
    real64 u_local[numNodesPerElem][numDofPerTrialSupportPoint];

    /// Stack storage for the element local nodal incremental displacement
    real64 uhat_local[numNodesPerElem][numDofPerTrialSupportPoint];

    real64 ( &localResidualMomentum )[numDispDofPerElem];
    real64 ( &dLocalResidualMomentum_dDisplacement )[numDispDofPerElem][numDispDofPerElem];
    
  };
  //*****************************************************************************

  /**
   * @brief Copy global values from primary field to a local stack array.
   * @copydoc ::geosx::finiteElement::ImplicitKernelBase::setup
   *
   * Global values from the displacement, incremental displacement, 
   * and degree of freedom numbers are placed into element local stack storage.
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

  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void quadraturePointKernel( localIndex const k,
                              localIndex const q,
                              StackVariables & stack ) const
  {
    // Governing equations (strong form)
    // ---------------------------------
    //
    //   divergence( stress - biotCoefficient*matrixPressure*I ) + fracturePressureTerm + bodyForces = 0   (quasi-static linear momentum balance)
    //  
    //   the pressures are assumed to be constant and the term with the biotCoefficient is simply implemented as a body force.
    //   Damage dependency is ignored here since the damage variable is not fully coupled with the displacement or pressure.
    //
    //   stress               = stress( strainIncrement )
    //   fracturePressureTerm = fracturePressureTerm( fracturePressure )
    //   bodyForce            = bodyForce( strainIncrement, pressure)
    //
    // Using a weak formulation of the governing equation the following terms are assembled in this kernel
    //
    //   Rmom = - \int symmetricGradient( \eta ) : stress +\int biotCoefficient*matrixPressure*grad(\eta) \cdot I + \int \eta \cdot fracturePressureTerm + \int \eta \cdot bodyForce = 0
    //       //
    //   dRmom_dVolStrain = - \int_Omega symmetricGradient( \eta ) : dstress_dVolStrain
    //                      + \int \eta \cdot dBodyForce_dVolStrain
    //
    // with \eta test basis functions for the displacement.
    // A continuous interpolation is used for the displacement, with \eta continuous finite element
    // basis functions. Matrix and Fracture pressures are taken from FV simulation, so, they are constant over cells.

    using namespace PDEUtilities;

    constexpr FunctionSpace displacementTrialSpace = FE_TYPE::template getFunctionSpace< numDofPerTrialSupportPoint >();
    constexpr FunctionSpace displacementTestSpace = displacementTrialSpace;

    real64 strainIncrement[6]{};
    real64 stress[6]{};
    typename CONSTITUTIVE_TYPE::KernelWrapper::DiscretizationOps stiffness; // Could this be called dTotalStress_dStrainIncrement?
    real64 const bodyForce[3] = { m_gravityVector[0] * m_solidDensity( k, q ),
                                  m_gravityVector[1] * m_solidDensity( k, q ),
                                  m_gravityVector[2] * m_solidDensity( k, q ) };

    //compute terms that account for the background pressure effects
    real64 fracturePressureTerm[3]{};
    computeFracturePressureTerm( k, q, fracturePressureTerm );
    real64 matrixPressureTerm = computeMatrixPressureTerm( k,q );
    //Dont need this because secondOrder identity is implemented in PDEUtilities
    //real64 const identityTensor[6] = {1, 1, 1, 0, 0, 0};
    //real64 matrixPressureTerm[6]{};
    //LvArray::tensorOps::scaledCopy< 6 >(  matrixPressureTerm, matrixPressureTermScalar, identityTensor );

    // Displacement finite element basis functions (N), basis function derivatives (dNdX), and
    // determinant of the Jacobian transformation matrix times the quadrature weight (detJxW)
    real64 N[numNodesPerElem];
    real64 dNdX[numNodesPerElem][3];
    FE_TYPE::calcN( q, N );
    real64 const detJxW = m_finiteElementSpace.template getGradN< FE_TYPE >( k, q, stack.xLocal, dNdX );

    // Compute strain increment
    FE_TYPE::symmetricGradient( dNdX, stack.uhat_local, strainIncrement );

    // Evaluate conserved quantities (total stress and fluid mass content) and their derivatives
    m_constitutiveUpdate.smallStrainUpdate( k,
                                            q,
                                            strainIncrement,
                                            stress,
                                            stiffness );

    // Compute local linear momentum balance residual
    LinearFormUtilities::compute< displacementTestSpace,
                                  DifferentialOperator::SymmetricGradient >
    (
      stack.localResidualMomentum,
      dNdX,
      stress,
      -detJxW );

    if( m_gravityAcceleration > 0.0 )
    {
      LinearFormUtilities::compute< displacementTestSpace,
                                    DifferentialOperator::Identity >
      (
        stack.localResidualMomentum,
        N,
        bodyForce,
        detJxW );
    }
    
    LinearFormUtilities::compute< displacementTestSpace,
                                  DifferentialOperator::SymmetricGradient >
    (
      stack.localResidualMomentum,
      dNdX,
      matrixPressureTerm,
      -detJxW );

    LinearFormUtilities::compute< displacementTestSpace,
                                  DifferentialOperator::Identity >
    (
      stack.localResidualMomentum,
      N,
      fracturePressureTerm,
      detJxW );

    // Compute local linear momentum balance residual derivatives with respect to displacement
    BilinearFormUtilities::compute< displacementTestSpace,
                                    displacementTrialSpace,
                                    DifferentialOperator::SymmetricGradient,
                                    DifferentialOperator::SymmetricGradient >
    (
      stack.dLocalResidualMomentum_dDisplacement,
      dNdX,
      stiffness, // fourth-order tensor handled via DiscretizationOps
      dNdX,
      -detJxW );

    //ASSUME THAT DENSITY IS CONSTANT
    // if( m_gravityAcceleration > 0.0 )
    // {
    //   BilinearFormUtilities::compute< displacementTestSpace,
    //                                   displacementTrialSpace,
    //                                   DifferentialOperator::Identity,
    //                                   DifferentialOperator::Divergence >
    //   (
    //     stack.dLocalResidualMomentum_dDisplacement,
    //     N,
    //     dBodyForce_dVolStrainIncrement,
    //     dNdX,
    //     detJxW );
    // }

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

    constexpr int nUDof = numNodesPerElem * numDofPerTestSupportPoint;

    for( int localNode = 0; localNode < numNodesPerElem; ++localNode )
    {
      for( int dim = 0; dim < numDofPerTestSupportPoint; ++dim )
      {
        localIndex const dof = LvArray::integerConversion< localIndex >( stack.localRowDofIndex[numDofPerTestSupportPoint*localNode + dim] - m_dofRankOffset );
        if( dof < 0 || dof >= m_matrix.numRows() ) continue;
        m_matrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( dof,
                                                                                stack.localRowDofIndex,
                                                                                stack.dLocalResidualMomentum_dDisplacement[numDofPerTestSupportPoint * localNode + dim],
                                                                                numNodesPerElem * numDofPerTrialSupportPoint );

        RAJA::atomicAdd< parallelDeviceAtomic >( &m_rhs[dof], stack.localResidualMomentum[numDofPerTestSupportPoint * localNode + dim] );
        maxForce = fmax( maxForce, fabs( stack.localResidualMomentum[numDofPerTestSupportPoint * localNode + dim] ) );

      }
    }

    return maxForce;
  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void computeFracturePressureTerm( localIndex const k,
                                    localIndex const q,
                                    real64 ( & fracturePressureTerm )[3]) const
  {

    real64 damageGrad[3]{};
    real64 pressureDamageGrad[3]{};

    m_constitutiveUpdate.getDamageGrad( k, q, damageGrad );

    real64 const damage = m_constitutiveUpdate.getDamage( k, q );
    real64 const pressureDamageDeriv = m_constitutiveUpdate.pressureDamageFunctionDerivative( damage ); 

    LvArray::tensorOps::scaledCopy< 3 >( pressureDamageGrad, damageGrad, pressureDamageDeriv );
    LvArray::tensorOps::scaledCopy< 3 >( fracturePressureTerm, pressureDamageGrad, m_pressureFracture[k] );

  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  real64 computeMatrixPressureTerm( localIndex const k,
                                    localIndex const q ) const
  {

    real64 const biotCoefficient = m_constitutiveUpdate.getBiotCoefficient( k );
    real64 const damage = m_constitutiveUpdate.getDamage( k, q );
    real64 const pressureDamageFunction = m_constitutiveUpdate.pressureDamageFunction( k, q ); 
    return pressureDamageFunction*biotCoefficient*m_pressureMatrix[k];

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
  real64 const m_gravityAcceleration;

  /// The rank global densities
  arrayView2d< real64 const > const m_solidDensity;

  /// The rank-global background pressure arrays.
  arrayView1d< real64 const > const m_pressureMatrix;
  arrayView1d< real64 const > const m_pressureFracture;

};

/// The factory used to construct a QuasiStaticPressure kernel.
using QuasiStaticPressureFactory = finiteElement::KernelFactory< QuasiStaticPressure,
                                                                 arrayView1d< globalIndex const > const,
                                                                 globalIndex,
                                                                 CRSMatrixView< real64, globalIndex const > const,
                                                                 arrayView1d< real64 > const,
                                                                 real64 const (&)[3] >;


} // namespace solidMechanicsLagrangianFEMKernels

} // namespace geosx

#include "finiteElement/kernelInterface/SparsityKernelBase.hpp"

#endif // GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSSMALLSTRAINQUASISTATICPRESSURE_HPP_
