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
 * @file SinglePhasePoroelasticKernel.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASEPOROMECHANICSKERNEL_HPP_
#define GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASEPOROMECHANICSKERNEL_HPP_

#include "finiteElement/BilinearFormUtilities.hpp"
#include "finiteElement/LinearFormUtilities.hpp"
#include "finiteElement/kernelInterface/ImplicitKernelBase.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseExtrinsicData.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBaseExtrinsicData.hpp"

namespace geosx
{

namespace PoromechanicsKernels
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
               arrayView1d< string const > const fluidModelNames ):
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
    m_fluidDensity( elementSubRegion.template getConstitutiveModel< constitutive::SingleFluidBase >( fluidModelNames[targetRegionIndex] ).density() ),
    m_fluidDensityOld( elementSubRegion.template getExtrinsicData< extrinsicMeshData::flow::densityOld >() ),
    m_initialFluidDensity( elementSubRegion.template getConstitutiveModel< constitutive::SingleFluidBase >( fluidModelNames[targetRegionIndex] ).initialDensity() ),
    m_dFluidDensity_dPressure( elementSubRegion.template getConstitutiveModel< constitutive::SingleFluidBase >( fluidModelNames[targetRegionIndex] ).dDensity_dPressure() ),
    m_flowDofNumber( elementSubRegion.template getReference< array1d< globalIndex > >( inputFlowDofKey )),
    m_initialFluidPressure( elementSubRegion.template getExtrinsicData< extrinsicMeshData::flow::initialPressure >() ),
    m_fluidPressureOld( elementSubRegion.template getExtrinsicData< extrinsicMeshData::flow::pressure >() ),
    m_deltaFluidPressure( elementSubRegion.template getExtrinsicData< extrinsicMeshData::flow::deltaPressure >() )
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
            uhat_local(),
            localResidualMomentum( Base::StackVariables::localResidual ),
      dLocalResidualMomentum_dDisplacement( Base::StackVariables::localJacobian ),
      dLocalResidualMomentum_dPressure{ {0.0} },
      localResidualMass{ 0.0 },
      dLocalResidualMass_dDisplacement{ {0.0} },
      dLocalResidualMass_dPressure{ {0.0} },
      localFlowDofIndex{ 0 }
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
    real64 dLocalResidualMomentum_dPressure[numDispDofPerElem][1];

    real64 localResidualMass[1];
    real64 dLocalResidualMass_dDisplacement[1][numDispDofPerElem];
    real64 dLocalResidualMass_dPressure[1][1];

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
        stack.xLocal[a][i] = m_X[localNodeIndex][i];
#endif
        stack.u_local[a][i] = m_disp[localNodeIndex][i];
        stack.uhat_local[a][i] = m_uhat[localNodeIndex][i];
        stack.localRowDofIndex[a*3+i] = m_dofNumber[localNodeIndex]+i;
        stack.localColDofIndex[a*3+i] = m_dofNumber[localNodeIndex]+i;
      }
    }

    stack.localFlowDofIndex[0] = m_flowDofNumber[k];

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
    //   divergence( totalStress ) + bodyForce = 0                       (quasi-static linear momentum balance)
    //   dFluidMassContent_dTime + divergence( fluidMassFlux ) = source  (fluid phase mass balance)
    //
    // with currently the following dependencies on the strainIncrement tensor and pressure
    //
    //   totalStress      = totalStress( strainIncrement, pressure)
    //   bodyForce        = bodyForce( strainIncrement, pressure)
    //   fluidMassContent = fluidMassContent( strainIncrement, pressure)
    //   fluidMassFlux    = fludiMassFlux( pressure)
    //
    // Note that the fluidMassFlux will depend on the straiIncrement if a stress-dependent constitutive
    // model is assumed. A dependency on pressure can also occur in the source term of the mass
    // balance equation, e.g. if a Peaceman well model is used.
    //
    // In this kernel cell-based contributions to Jacobian matrix and residual
    // vector are computed. The face-based contributions, namely associated with the fluid
    // mass flux term are computed in a different kernel. The source term in the mass balance
    // equation is also treated elsewhere.
    //
    // Integration in time is performed using a backward Euler scheme (in the mass balance
    // equation LHS and RHS are multiplied by the timestep).
    //
    // Using a weak formulation of the governing equation the following terms are assembled in this kernel
    //
    //   Rmom = - \int symmetricGradient( \eta ) : totalStress + \int \eta \cdot bodyForce = 0
    //   Rmas = \int \chi ( fluidMassContent - fluidMassContentOld) = 0
    //
    //   dRmom_dVolStrain = - \int_Omega symmetricGradient( \eta ) : dTotalStress_dVolStrain
    //                      + \int \eta \cdot dBodyForce_dVolStrain
    //   dRmom_dPressure  = - \int_Omega symmetricGradient( \eta ) : dTotalStress_dPressure
    //                      + \int \eta \cdot dBodyForce_dPressure
    //   dRmas_dVolStrain = \int \chi dFluidMassContent_dVolStrain
    //   dRmas_dPressure  = \int \chi dFluidMassContent_dPressure
    //
    // with \eta and \chi test basis functions for the displacement and pressure field, respectively.
    // A continuous interpolation is used for the displacement, with \eta continuous finite element
    // basis functions. A piecewise-constant approximation is used for the pressure.

    using namespace PDEUtilities;

    constexpr FunctionSpace displacementTrialSpace = FE_TYPE::template getFunctionSpace< numDofPerTrialSupportPoint >();
    constexpr FunctionSpace displacementTestSpace = displacementTrialSpace;
    constexpr FunctionSpace pressureTrialSpace = FunctionSpace::P0;
    constexpr FunctionSpace pressureTestSpace = pressureTrialSpace;

    real64 strainIncrement[6]{};
    real64 totalStress[6]{};
    typename CONSTITUTIVE_TYPE::KernelWrapper::DiscretizationOps stiffness; // Could this be called dTotalStress_dStrainIncrement?
    real64 dTotalStress_dPressure[6]{};
    real64 bodyForce[3]{};
    real64 dBodyForce_dVolStrainIncrement[3]{};
    real64 dBodyForce_dPressure[3]{};
    real64 fluidMassContentIncrement;
    real64 dFluidMassContent_dPressure;
    real64 dFluidMassContent_dVolStrainIncrement;

    // Displacement finite element basis functions (N), basis function derivatives (dNdX), and
    // determinant of the Jacobian transformation matrix times the quadrature weight (detJxW)
    real64 N[numNodesPerElem];
    real64 dNdX[numNodesPerElem][3];
    FE_TYPE::calcN( q, N );
    real64 const detJxW = m_finiteElementSpace.template getGradN< FE_TYPE >( k, q, stack.xLocal, dNdX );

    // Compute strain increment
    FE_TYPE::symmetricGradient( dNdX, stack.uhat_local, strainIncrement );

    // Evaluate conserved quantities (total stress and fluid mass content) and their derivatives
    m_constitutiveUpdate.smallStrainUpdateSinglePhase( k,
                                                       q,
                                                       m_initialFluidPressure[k],
                                                       m_fluidPressureOld[k],
                                                       m_deltaFluidPressure[k],
                                                       strainIncrement,
                                                       m_gravityAcceleration,
                                                       m_gravityVector,
                                                       m_solidDensity( k, q ),
                                                       m_initialFluidDensity( k, q ),
                                                       m_fluidDensityOld( k ),
                                                       m_fluidDensity( k, q ),
                                                       m_dFluidDensity_dPressure( k, q ),
                                                       totalStress,
                                                       dTotalStress_dPressure,
                                                       bodyForce,
                                                       dBodyForce_dVolStrainIncrement,
                                                       dBodyForce_dPressure,
													   fluidMassContentIncrement,
                                                       dFluidMassContent_dPressure,
                                                       dFluidMassContent_dVolStrainIncrement,
                                                       stiffness );

    // Compute local linear momentum balance residual
    LinearFormUtilities::compute< displacementTestSpace,
                                  DifferentialOperator::SymmetricGradient >
    (
      stack.localResidualMomentum,
      dNdX,
      totalStress,
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

    // --- dBodyForce_dVoldStrain derivative neglected

    // Compute local linear momentum balance residual derivatives with respect to pressure
    real64 Np[1] = { 1.0 };

    BilinearFormUtilities::compute< displacementTestSpace,
                                    pressureTrialSpace,
                                    DifferentialOperator::SymmetricGradient,
                                    DifferentialOperator::Identity >
    (
      stack.dLocalResidualMomentum_dPressure,
      dNdX,
      dTotalStress_dPressure,
      Np,
      -detJxW );

    if( m_gravityAcceleration > 0.0 )
    {
      BilinearFormUtilities::compute< displacementTestSpace,
                                      pressureTrialSpace,
                                      DifferentialOperator::Identity,
                                      DifferentialOperator::Identity >
      (
        stack.dLocalResidualMomentum_dPressure,
        N,
        dBodyForce_dPressure,
        Np,
        detJxW );
    }


    // Compute local mass balance residual
    LinearFormUtilities::compute< pressureTestSpace,
                                  DifferentialOperator::Identity >
    (
      stack.localResidualMass,
      Np,
	  fluidMassContentIncrement,
      detJxW );

    // Compute local mass balance residual derivatives with respect to displacement
    BilinearFormUtilities::compute< pressureTestSpace,
                                    displacementTrialSpace,
                                    DifferentialOperator::Identity,
                                    DifferentialOperator::Divergence >
    (
      stack.dLocalResidualMass_dDisplacement,
      Np,
      dFluidMassContent_dVolStrainIncrement,
      dNdX,
      detJxW );

    // Compute local mass balance residual derivatives with respect to pressure
    BilinearFormUtilities::compute< pressureTestSpace,
                                    pressureTrialSpace,
                                    DifferentialOperator::Identity,
                                    DifferentialOperator::Identity >
    (
      stack.dLocalResidualMass_dPressure,
      Np,
      dFluidMassContent_dPressure,
      Np,
      detJxW );
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

        m_matrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( dof,
                                                                                stack.localFlowDofIndex,
                                                                                stack.dLocalResidualMomentum_dPressure[numDofPerTestSupportPoint * localNode + dim],
                                                                                1 );

      }
    }


    localIndex const dof = LvArray::integerConversion< localIndex >( stack.localFlowDofIndex[0] - m_dofRankOffset );
    if( 0 <= dof && dof < m_matrix.numRows() )
    {
      m_matrix.template addToRowBinarySearchUnsorted< serialAtomic >( dof,
                                                                      stack.localRowDofIndex,
                                                                      stack.dLocalResidualMass_dDisplacement[0],
                                                                      nUDof );
      m_matrix.template addToRow< serialAtomic >( dof,
                                                  stack.localFlowDofIndex,
                                                  stack.dLocalResidualMass_dPressure[0],
                                                  1 );
      RAJA::atomicAdd< serialAtomic >( &m_rhs[dof], stack.localResidualMass[0] );
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
  real64 const m_gravityAcceleration;

  /// The rank global densities
  arrayView2d< real64 const > const m_solidDensity;
  arrayView2d< real64 const > const m_fluidDensity;
  arrayView1d< real64 const > const m_fluidDensityOld;
  arrayView2d< real64 const > const m_initialFluidDensity;
  arrayView2d< real64 const > const m_dFluidDensity_dPressure;

  /// The global degree of freedom number
  arrayView1d< globalIndex const > const m_flowDofNumber;

  /// The rank-global initial fluid pressure array
  arrayView1d< real64 const > const m_initialFluidPressure;

  /// The rank-global fluid pressure array.
  arrayView1d< real64 const > const m_fluidPressureOld;

  /// The rank-global delta-fluid pressure array.
  arrayView1d< real64 const > const m_deltaFluidPressure;

};

using SinglePhaseKernelFactory = finiteElement::KernelFactory< SinglePhase,
                                                               arrayView1d< globalIndex const > const &,
                                                               string const &,
                                                               globalIndex const,
                                                               CRSMatrixView< real64, globalIndex const > const &,
                                                               arrayView1d< real64 > const &,
                                                               real64 const (&)[3],
                                                               arrayView1d< string const > const >;

} // namespace PoroelasticKernels

} // namespace geosx

#include "finiteElement/kernelInterface/SparsityKernelBase.hpp"

#endif // GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASEPOROMECHANICSKERNEL_HPP_
