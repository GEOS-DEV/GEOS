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
 * @file PhaseFieldDamagePressureKernels.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_SIMPLEPDE_PHASEFIELDDAMAGEPRESSUREKERNELS_HPP_
#define GEOSX_PHYSICSSOLVERS_SIMPLEPDE_PHASEFIELDDAMAGEPRESSUREKERNELS_HPP_

#include "finiteElement/BilinearFormUtilities.hpp"
#include "finiteElement/LinearFormUtilities.hpp"
#include "finiteElement/kernelInterface/ImplicitKernelBase.hpp"

namespace geosx
{
//*****************************************************************************
/**
 * @brief Implements kernels for solving the Damage(or phase-field) equation
 * in a phase-field fracture problem, with background pressure effects.
 * @copydoc geosx::finiteElement::KernelBase
 * @tparam NUM_NODES_PER_ELEM The number of nodes per element for the
 *                            @p SUBREGION_TYPE.
 * @tparam UNUSED An unused parameter since we are assuming that the test and
 *                trial space have the same number of support points.
 *
 * ### PhaseFieldDamageKernel Description
 * Implements the KernelBase interface functions required for solving the
 * Damage(or phase-field) equation in a phase-field fracture problem.
 * It uses the finite element kernel application functions such as
 * geosx::finiteElement::RegionBasedKernelApplication.
 *
 * In this implementation, the template parameter @p NUM_NODES_PER_ELEM is used
 * in place of both @p NUM_TEST_SUPPORT_POINTS_PER_ELEM and
 * @p NUM_TRIAL_SUPPORT_POINTS_PER_ELEM, which are assumed to be equal. This
 * results in the @p UNUSED template parameter as only the NUM_NODES_PER_ELEM
 * is passed to the ImplicitKernelBase template to form the base class.
 *
 * Additionally, the number of degrees of freedom per support point for both
 * the test and trial spaces are specified as `1` when specifying the base
 * class.
 */
template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
class PhaseFieldDamagePressureKernel :
  public finiteElement::ImplicitKernelBase< SUBREGION_TYPE,
                                            CONSTITUTIVE_TYPE,
                                            FE_TYPE,
                                            1,
                                            1 >
{
public:
  /// An alias for the base class.
  using Base = finiteElement::ImplicitKernelBase< SUBREGION_TYPE,
                                                  CONSTITUTIVE_TYPE,
                                                  FE_TYPE,
                                                  1,
                                                  1 >;

  using Base::numDofPerTestSupportPoint;
  using Base::numDofPerTrialSupportPoint;
  using Base::m_dofNumber;
  using Base::m_dofRankOffset;
  using Base::m_matrix;
  using Base::m_rhs;
  using Base::m_elemsToNodes;
  using Base::m_constitutiveUpdate;
  using Base::m_finiteElementSpace;

  /// Maximum number of nodes per element, which is equal to the maxNumTestSupportPointPerElem and
  /// maxNumTrialSupportPointPerElem by definition. When the FE_TYPE is not a Virtual Element, this
  /// will be the actual number of nodes per element.
  static constexpr int numNodesPerElem = Base::maxNumTestSupportPointsPerElem;

  /**
   * @brief Constructor
   * @copydoc geosx::finiteElement::ImplicitKernelBase::ImplicitKernelBase
   * @param damageName The name of the damage field, usually just "Damage"
   *                  (i.e. Temperature, Pressure, etc.)
   */
  PhaseFieldDamagePressureKernel( NodeManager const & nodeManager,
                                  EdgeManager const & edgeManager,
                                  FaceManager const & faceManager,
                                  localIndex const targetRegionIndex,
                                  SUBREGION_TYPE const & elementSubRegion,
                                  FE_TYPE const & finiteElementSpace,
                                  CONSTITUTIVE_TYPE & inputConstitutiveType,
                                  arrayView1d< globalIndex const > const inputDofNumber,
                                  globalIndex const rankOffset,
                                  CRSMatrixView< real64, globalIndex const > const inputMatrix,
                                  arrayView1d< real64 > const inputRhs,
                                  string const damageName,
                                  int const localDissipationOption):
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
          inputRhs ),
    m_X( nodeManager.referencePosition()),
    m_disp( nodeManager.totalDisplacement()),
    m_nodalDamage( nodeManager.template getReference< array1d< real64 > >( damageName )),
    m_localDissipationOption( localDissipationOption ),
    //m_pressureMatrix( elementSubRegion.template getExtrinsicData< extrinsicMeshData::flow::matrixPressure >() ),
    //m_pressureFracture( elementSubRegion.template getExtrinsicData< extrinsicMeshData::flow::fracturePressure >() )
    //this should compile
    m_pressureMatrix( elementSubRegion.template getReference< array1d< real64 > >( "hardCodedPMatrixName" ) ),
    m_pressureFracture( elementSubRegion.template getReference< array1d< real64 > >( "hardCodedPFractureName" ) )

  {}

  //***************************************************************************
  /**
   * @class StackVariables
   * @copydoc geosx::finiteElement::ImplicitKernelBase::StackVariables
   *
   * Adds a stack array for the primary field.
   */
  struct StackVariables : Base::StackVariables
  {
public:

    /**
     * @brief Constructor
     */
    GEOSX_HOST_DEVICE
    StackVariables():
      Base::StackVariables(),
            xLocal(),
            nodalDamageLocal{ 0.0 }
    {}

#if !defined(CALC_FEM_SHAPE_IN_KERNEL)
    /// Dummy
    int xLocal;
#else
    /// C-array stack storage for element local the nodal positions.
    real64 xLocal[ numNodesPerElem ][ 3 ];
#endif

    /// C-array storage for the element local primary field variable.
    real64 nodalDamageLocal[ numNodesPerElem ];

    /// Stack storage for the element displacement vector.
    real64 dispLocal[ numNodesPerElem ][ 3 ];
  };


  /**
   * @brief Copy global values from primary field to a local stack array.
   * @copydoc geosx::finiteElement::ImplicitKernelBase::setup
   *
   * For the PhaseFieldDamageKernel implementation, global values from the
   * primaryField, and degree of freedom numbers are placed into element local
   * stack storage.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void setup( localIndex const k,
              StackVariables & stack ) const
  {
    for( localIndex a=0; a<numNodesPerElem; ++a )
    {
      localIndex const localNodeIndex = m_elemsToNodes( k, a );

#if defined(CALC_FEM_SHAPE_IN_KERNEL)
      LvArray::tensorOps::copy< 3 >( stack.xLocal[ a ], m_X[ localNodeIndex ] );
#endif
      
      for( localIndex i=0; i<3; ++i)
      {
        stack.dispLocal[a][i] = m_disp[ localNodeIndex ][ i ];
      }
      stack.nodalDamageLocal[ a ] = m_nodalDamage[ localNodeIndex ];
      stack.localRowDofIndex[a] = m_dofNumber[localNodeIndex];
      stack.localColDofIndex[a] = m_dofNumber[localNodeIndex];
    }
  }

  /**
   * @copydoc geosx::finiteElement::ImplicitKernelBase::quadraturePointJacobianContribution
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void quadraturePointKernel( localIndex const k,
                              localIndex const q,
                              StackVariables & stack ) const
  {

    using namespace PDEUtilities;

    constexpr FunctionSpace damageTrialSpace = FE_TYPE::template getFunctionSpace< numDofPerTrialSupportPoint >();
    constexpr FunctionSpace damageTestSpace = damageTrialSpace;

    real64 const strainEnergyDensity = m_constitutiveUpdate.getStrainEnergyDensity( k, q );
    real64 const ell = m_constitutiveUpdate.getRegularizationLength();
    real64 const Gc = m_constitutiveUpdate.getCriticalFractureEnergy();
    real64 const threshold = m_constitutiveUpdate.getEnergyThreshold();
    //real64 const extDrivingForce = m_constitutiveUpdate.getExtDrivingForce( k, q );
    real64 const volStrain = m_constitutiveUpdate.getVolStrain( k, q );
    real64 const biotCoeff = m_constitutiveUpdate.getBiotCoefficient( k );
    real64 c0 = 2;

    //Interpolate d and grad_d
    real64 N[ numNodesPerElem ];
    real64 dNdX[ numNodesPerElem ][ 3 ];
    real64 const detJxW = m_finiteElementSpace.template getGradN< FE_TYPE >( k, q, stack.xLocal, dNdX );
    FE_TYPE::calcN( q, N );

    real64 qpDamage = 0.0;
    real64 qpGradDamage[3] = {0, 0, 0};
    real64 qpDisp[3] = {0, 0, 0};
    FE_TYPE::valueAndGradient( N, dNdX, stack.nodalDamageLocal, qpDamage, qpGradDamage );
    FE_TYPE::value(N, stack.dispLocal, qpDisp);

    real64 D = 0;                                                                   //max between threshold and
                                                                                    // Elastic energy
    if( m_localDissipationOption == 1 )
    {
      c0 = 8/3;
      D = fmax( threshold, strainEnergyDensity );
    }
    
    if( m_localDissipationOption == 1 ) //AT1 KERNELS
    {

        // Compute local dissipation term AT1
        LinearFormUtilities::compute< damageTestSpace,
                                      DifferentialOperator::Identity >
        (
          stack.localResidual,
          N,
          (Gc/(c0*ell))*1,
          detJxW );

        // Compute driving force term
        LinearFormUtilities::compute< damageTestSpace,
                                      DifferentialOperator::Identity >
        (
          stack.localResidual,
          N,
          m_constitutiveUpdate.getDegradationDerivative( qpDamage )*D,
          detJxW );  
    }
    else //AT2 KERNELS
    {
        // Compute local dissipation term AT2
        LinearFormUtilities::compute< damageTestSpace,
                                      DifferentialOperator::Identity >
        (
          stack.localResidual,
          N,
          (Gc/(c0*ell))*2*qpDamage,
          detJxW );


        // Compute driving force term
        LinearFormUtilities::compute< damageTestSpace,
                                      DifferentialOperator::Identity >
        (
          stack.localResidual,
          N,
          m_constitutiveUpdate.getDegradationDerivative( qpDamage )*strainEnergyDensity,
          detJxW );  

    }

    // Compute non-local term
    real64 nonLocalTermIntegrand[3]; 
    LvArray::tensorOps::scaledCopy< 3 >( nonLocalTermIntegrand, qpGradDamage, (2*Gc*ell/c0) );
    LinearFormUtilities::compute< damageTestSpace,
                                  DifferentialOperator::Gradient >
    (
        stack.localResidual,
        dNdX,
        nonLocalTermIntegrand,
        detJxW );


    // Compute matrix pressure contribution
    LinearFormUtilities::compute< damageTestSpace,
                                  DifferentialOperator::Identity >
    (
        stack.localResidual,
        N,
        biotCoeff*m_pressureMatrix[k]*volStrain,
        -detJxW );

    real64 pressureFractureIntegrand[3]; 
    LvArray::tensorOps::scaledCopy< 3 >( pressureFractureIntegrand, qpDisp, m_pressureFracture[k] );
    // Compute fracture pressure contribuition
    LinearFormUtilities::compute< damageTestSpace,
                                  DifferentialOperator::Gradient >
    (
        stack.localResidual,
        dNdX,
        pressureFractureIntegrand,
        detJxW );

    // // Compute jacobian terms

    if( m_localDissipationOption == 1 )
    {
        //no contribution to Jacobian
        BilinearFormUtilities::compute< damageTestSpace,
                                        damageTrialSpace,
                                        DifferentialOperator::Identity,
                                        DifferentialOperator::Identity >
        (
            stack.localJacobian,
            N,
            D*m_constitutiveUpdate.getDegradationSecondDerivative( qpDamage ), 
            N,
            detJxW );

    }
    else
    {

        BilinearFormUtilities::compute< damageTestSpace,
                                        damageTrialSpace,
                                        DifferentialOperator::Identity,
                                        DifferentialOperator::Identity >
        (
            stack.localJacobian,
            N,
            strainEnergyDensity*m_constitutiveUpdate.getDegradationSecondDerivative( qpDamage ), 
            N,
            detJxW );      

        BilinearFormUtilities::compute< damageTestSpace,
                                        damageTrialSpace,
                                        DifferentialOperator::Identity,
                                        DifferentialOperator::Identity >
        (
            stack.localJacobian,
            N,
            (2*Gc/(c0*ell)), 
            N,
            detJxW );

    }
    BilinearFormUtilities::compute< damageTestSpace,
                                    damageTrialSpace,
                                    DifferentialOperator::Gradient,
                                    DifferentialOperator::Gradient >
    (
        stack.localJacobian,
        dNdX,
        (2*Gc*ell/c0), 
        dNdX,
        detJxW );

  }

  /**
   * @copydoc geosx::finiteElement::ImplicitKernelBase::complete
   *
   * Form element residual from the fully formed element Jacobian dotted with
   * the primary field and map the element local Jacobian/Residual to the
   * global matrix/vector.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  real64 complete( localIndex const k,
                   StackVariables & stack ) const
  {
    GEOSX_UNUSED_VAR( k );
    real64 maxForce = 0;

    for( int a = 0; a < numNodesPerElem; ++a )
    {
      localIndex const dof = LvArray::integerConversion< localIndex >( stack.localRowDofIndex[ a ] - m_dofRankOffset );
      if( dof < 0 || dof >= m_matrix.numRows() ) continue;
      m_matrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( dof,
                                                                              stack.localColDofIndex,
                                                                              stack.localJacobian[ a ],
                                                                              numNodesPerElem );

      RAJA::atomicAdd< parallelDeviceAtomic >( &m_rhs[ dof ], stack.localResidual[ a ] );
      maxForce = fmax( maxForce, fabs( stack.localResidual[ a ] ) );
    }

    return maxForce;
  }



protected:
  /// The array containing the nodal position array.
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const m_X;
  arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const m_disp;


  /// The global primary field array.
  arrayView1d< real64 const > const m_nodalDamage;

  int const m_localDissipationOption;

  arrayView1d< real64 const > const m_pressureMatrix;
  arrayView1d< real64 const > const m_pressureFracture;
};

using PhaseFieldDamagePressureKernelFactory = finiteElement::KernelFactory< PhaseFieldDamagePressureKernel,
                                                                            arrayView1d< globalIndex const > const,
                                                                            globalIndex,
                                                                            CRSMatrixView< real64, globalIndex const > const,
                                                                            arrayView1d< real64 > const,
                                                                            string const,
                                                                            int >;

} // namespace geosx

#include "finiteElement/kernelInterface/SparsityKernelBase.hpp"

#endif // GEOSX_PHYSICSSOLVERS_SIMPLEPDE_PHASEFIELDDAMAGEPRESSUREKERNELS_HPP_
