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
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBaseFields.hpp"

namespace geosx
{

namespace poromechanicsKernels
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
  SinglePhase( NodeManager const & nodeManager,
               EdgeManager const & edgeManager,
               FaceManager const & faceManager,
               localIndex const targetRegionIndex,
               SUBREGION_TYPE const & elementSubRegion,
               FE_TYPE const & finiteElementSpace,
               CONSTITUTIVE_TYPE & inputConstitutiveType,
               arrayView1d< globalIndex const > const inputDispDofNumber,
               string const inputFlowDofKey,
               globalIndex const rankOffset,
               CRSMatrixView< real64, globalIndex const > const inputMatrix,
               arrayView1d< real64 > const inputRhs,
               real64 const (&inputGravityVector)[3],
               string const fluidModelKey ):
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
    m_X( nodeManager.referencePosition() ),
    m_disp( nodeManager.getField< fields::solidMechanics::totalDisplacement >() ),
    m_uhat( nodeManager.getField< fields::solidMechanics::incrementalDisplacement >() ),
    m_gravityVector{ inputGravityVector[0], inputGravityVector[1], inputGravityVector[2] },
    m_gravityAcceleration( LvArray::tensorOps::l2Norm< 3 >( inputGravityVector ) ),
    m_solidDensity( inputConstitutiveType.getDensity() ),
    m_fluidDensity( elementSubRegion.template getConstitutiveModel< constitutive::SingleFluidBase >( elementSubRegion.template getReference< string >( fluidModelKey ) ).density() ),
    m_fluidDensity_n( elementSubRegion.template getConstitutiveModel< constitutive::SingleFluidBase >( elementSubRegion.template getReference< string >( fluidModelKey ) ).density_n() ),
    m_dFluidDensity_dPressure( elementSubRegion.template getConstitutiveModel< constitutive::SingleFluidBase >( elementSubRegion.template getReference< string >( fluidModelKey ) ).dDensity_dPressure() ),
    m_flowDofNumber( elementSubRegion.template getReference< array1d< globalIndex > >( inputFlowDofKey )),
    m_fluidPressure_n( elementSubRegion.template getField< fields::flow::pressure_n >() ),
    m_fluidPressure( elementSubRegion.template getField< fields::flow::pressure >() )
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

    // Storage for displacements

    /// Stack storage for the element local nodal displacement
    real64 u_local[numNodesPerElem][numDofPerTrialSupportPoint];
    /// Stack storage for the element local nodal incremental displacement
    real64 uhat_local[numNodesPerElem][numDofPerTrialSupportPoint];


    // Storage for helper variables used in the quadrature point kernel

    /// Total stress
    real64 totalStress[6]{};
    /// Derivative of total stress wrt displacement
    typename CONSTITUTIVE_TYPE::KernelWrapper::DiscretizationOps stiffness;
    /// Derivative of total stress wrt pressure
    real64 dTotalStress_dPressure[6]{};
    /// Derivative of total stress wrt temperature
    real64 dTotalStress_dTemperature[6]{};

    /// Body force
    real64 bodyForce[3]{};
    /// Derivative of body force wrt volumetric strain increment
    real64 dBodyForce_dVolStrainIncrement[3]{};
    /// Derivative of body force wrt pressure
    real64 dBodyForce_dPressure[3]{};

    /// Mass accumulation
    real64 fluidMassIncrement;
    /// Derivative of mass accumulation wrt volumetric strain increment
    real64 dFluidMassIncrement_dVolStrainIncrement;
    /// Derivative of mass accumulation wrt pressure
    real64 dFluidMassIncrement_dPressure;

    /// Delta temperature
    real64 deltaTemperature = 0.0;

    // Storage for residual and degrees of freedom

    /// Linear momentum balance residual
    real64 ( &localResidualMomentum )[numDispDofPerElem];
    /// Derivative of linear momentum balance residual wrt displacement
    real64 ( &dLocalResidualMomentum_dDisplacement )[numDispDofPerElem][numDispDofPerElem];
    /// Derivative of linear momentum balance residual wrt pressure
    real64 dLocalResidualMomentum_dPressure[numDispDofPerElem][1];

    /// Mass balance residual
    real64 localResidualMass[1];
    /// Derivative of mass balance residual wrt displacement
    real64 dLocalResidualMass_dDisplacement[1][numDispDofPerElem];
    /// Derivative of mass balance residual wrt pressure
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

  /**
   * @brief Helper function to compute 1) the total stress, 2) the body force term, and 3) the fluidMassIncrement
   * using quantities returned by the PorousSolid constitutive model.
   * This function also computes the derivatives of these three quantities wrt primary variables
   * @param[in] k the element index
   * @param[in] q the quadrature point index
   * @param[in] strainIncrement the strain increment used in total stress and porosity computation
   * @param[inout] stack the stack variables
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void smallStrainUpdate( localIndex const k,
                          localIndex const q,
                          real64 const ( &strainIncrement )[6],
                          StackVariables & stack ) const
  {
    real64 porosity = 0.0;
    real64 porosity_n = 0.0;
    real64 dPorosity_dVolStrain = 0.0;
    real64 dPorosity_dPressure = 0.0;
    real64 dPorosity_dTemperature = 0.0;
    real64 dSolidDensity_dPressure = 0.0;

    // Step 1: call the constitutive model to evaluate the total stress and compute porosity
    m_constitutiveUpdate.smallStrainUpdatePoromechanics( k, q,
                                                         m_fluidPressure_n[k],
                                                         m_fluidPressure[k],
                                                         stack.deltaTemperature,
                                                         strainIncrement,
                                                         stack.totalStress,
                                                         stack.dTotalStress_dPressure,
                                                         stack.dTotalStress_dTemperature,
                                                         stack.stiffness,
                                                         porosity,
                                                         porosity_n,
                                                         dPorosity_dVolStrain,
                                                         dPorosity_dPressure,
                                                         dPorosity_dTemperature,
                                                         dSolidDensity_dPressure );

    // Step 2: compute the body force
    computeBodyForce( k, q,
                      porosity,
                      dPorosity_dVolStrain,
                      dPorosity_dPressure,
                      dPorosity_dTemperature,
                      dSolidDensity_dPressure,
                      stack );

    // Step 3: compute fluid mass increment
    computeFluidIncrement( k, q,
                           porosity,
                           porosity_n,
                           dPorosity_dVolStrain,
                           dPorosity_dPressure,
                           dPorosity_dTemperature,
                           stack );
  }

  /**
   * @brief Helper function to compute the body force term and its derivatives wrt primary variables
   * @param[in] k the element index
   * @param[in] q the quadrature point index
   * @param[in] porosity the element porosity
   * @param[in] dPorosity_dVolStrain the derivative of porosity wrt volumetric strain increment
   * @param[in] dPorosity_dPressure the derivative of porosity wrt pressure
   * @param[in] dPorosity_dTemperature the derivative of porosity wrt temperature
   * @param[in] dSolidDensity_dPressure the derivative of solid density wrt pressure
   * @param[inout] stack the stack variables
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void computeBodyForce( localIndex const k,
                         localIndex const q,
                         real64 const & porosity,
                         real64 const & dPorosity_dVolStrain,
                         real64 const & dPorosity_dPressure,
                         real64 const & dPorosity_dTemperature,
                         real64 const & dSolidDensity_dPressure,
                         StackVariables & stack ) const
  {
    GEOSX_UNUSED_VAR( dPorosity_dTemperature );

    real64 const mixtureDensity = ( 1.0 - porosity ) * m_solidDensity( k, q ) + porosity * m_fluidDensity( k, q );
    real64 const dMixtureDens_dVolStrainIncrement = dPorosity_dVolStrain * ( -m_solidDensity( k, q ) + m_fluidDensity( k, q ) );
    real64 const dMixtureDens_dPressure = dPorosity_dPressure * ( -m_solidDensity( k, q ) + m_fluidDensity( k, q ) )
                                          + ( 1.0 - porosity ) * dSolidDensity_dPressure
                                          + porosity * m_dFluidDensity_dPressure( k, q );

    LvArray::tensorOps::scaledCopy< 3 >( stack.bodyForce, m_gravityVector, mixtureDensity );
    LvArray::tensorOps::scaledCopy< 3 >( stack.dBodyForce_dVolStrainIncrement, m_gravityVector, dMixtureDens_dVolStrainIncrement );
    LvArray::tensorOps::scaledCopy< 3 >( stack.dBodyForce_dPressure, m_gravityVector, dMixtureDens_dPressure );
  }

  /**
   * @brief Helper function to compute the fluid mass increment and its derivatives wrt primary variables
   * @param[in] k the element index
   * @param[in] q the quadrature point index
   * @param[in] porosity the element porosity
   * @param[in] porosity_n the element porosity at the previous converged time step
   * @param[in] dPorosity_dVolStrain the derivative of porosity wrt volumetric strain increment
   * @param[in] dPorosity_dPressure the derivative of porosity wrt pressure
   * @param[in] dPorosity_dTemperature the derivative of porosity wrt temperature
   * @param[inout] stack the stack variables
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void computeFluidIncrement( localIndex const k,
                              localIndex const q,
                              real64 const & porosity,
                              real64 const & porosity_n,
                              real64 const & dPorosity_dVolStrain,
                              real64 const & dPorosity_dPressure,
                              real64 const & dPorosity_dTemperature,
                              StackVariables & stack ) const
  {
    GEOSX_UNUSED_VAR( dPorosity_dTemperature );

    stack.fluidMassIncrement = porosity * m_fluidDensity( k, q ) - porosity_n * m_fluidDensity_n( k, q );
    stack.dFluidMassIncrement_dVolStrainIncrement = dPorosity_dVolStrain * m_fluidDensity( k, q );
    stack.dFluidMassIncrement_dPressure = dPorosity_dPressure * m_fluidDensity( k, q ) + porosity * m_dFluidDensity_dPressure( k, q );
  }

  /**
   * @brief Assemble the local linear momentum balance residual using total stress and body force terms
   * @param[in] N displacement finite element basis functions
   * @param[in] dNdX basis function derivatives
   * @param[in] detJxW determinant of the Jacobian transformation matrix times the quadrature weight
   * @param[inout] stack the stack variables
   * @detail This function assembles the discretized form of the following equation
   *   divergence( totalStress ) + bodyForce = 0
   * with the following dependencies on the strainIncrement tensor and pressure
   *   totalStress = totalStress( strainIncrement, pressure)
   *   bodyForce   = bodyForce( strainIncrement, pressure)
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void assembleMomentumBalanceResidual( real64 const ( &N )[numNodesPerElem],
                                        real64 const ( &dNdX )[numNodesPerElem][3],
                                        real64 const & detJxW,
                                        StackVariables & stack ) const
  {
    using namespace PDEUtilities;

    constexpr FunctionSpace displacementTrialSpace = FE_TYPE::template getFunctionSpace< numDofPerTrialSupportPoint >();
    constexpr FunctionSpace displacementTestSpace = displacementTrialSpace;
    constexpr FunctionSpace pressureTrialSpace = FunctionSpace::P0;

    // Step 1: compute local linear momentum balance residual
    LinearFormUtilities::compute< displacementTestSpace,
                                  DifferentialOperator::SymmetricGradient >
    (
      stack.localResidualMomentum,
      dNdX,
      stack.totalStress,
      -detJxW );

    if( m_gravityAcceleration > 0.0 )
    {
      LinearFormUtilities::compute< displacementTestSpace,
                                    DifferentialOperator::Identity >
      (
        stack.localResidualMomentum,
        N,
        stack.bodyForce,
        detJxW );
    }

    // Step 2: compute local linear momentum balance residual derivatives with respect to displacement
    BilinearFormUtilities::compute< displacementTestSpace,
                                    displacementTrialSpace,
                                    DifferentialOperator::SymmetricGradient,
                                    DifferentialOperator::SymmetricGradient >
    (
      stack.dLocalResidualMomentum_dDisplacement,
      dNdX,
      stack.stiffness, // fourth-order tensor handled via DiscretizationOps
      dNdX,
      -detJxW );

    if( m_gravityAcceleration > 0.0 )
    {
      BilinearFormUtilities::compute< displacementTestSpace,
                                      displacementTrialSpace,
                                      DifferentialOperator::Identity,
                                      DifferentialOperator::Divergence >
      (
        stack.dLocalResidualMomentum_dDisplacement,
        N,
        stack.dBodyForce_dVolStrainIncrement,
        dNdX,
        detJxW );
    }

    // Step 3: compute local linear momentum balance residual derivatives with respect to pressure
    BilinearFormUtilities::compute< displacementTestSpace,
                                    pressureTrialSpace,
                                    DifferentialOperator::SymmetricGradient,
                                    DifferentialOperator::Identity >
    (
      stack.dLocalResidualMomentum_dPressure,
      dNdX,
      stack.dTotalStress_dPressure,
      1.0,
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
        stack.dBodyForce_dPressure,
        1.0,
        detJxW );
    }
  }

  /**
   * @brief Assemble the local mass balance using fluid mass increment
   * @param[in] dNdX basis function derivatives
   * @param[in] detJxW determinant of the Jacobian transformation matrix times the quadrature weight
   * @param[inout] stack the stack variables
   * @detail This function assembles the discretized form of the temporal derivative in the following equation
   *   dFluidMass_dTime + divergence( fluidMassFlux ) = source  (fluid phase mass balance)
   * with the following dependencies on the strainIncrement tensor and pressure
   *   fluidMass = fluidMass( strainIncrement, pressure)
   *   fluidMassFlux = fluidMassFlux( pressure)
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void assembleElementBasedFlowResidual( real64 const ( &dNdX )[numNodesPerElem][3],
                                         real64 const & detJxW,
                                         StackVariables & stack ) const
  {

    using namespace PDEUtilities;

    constexpr FunctionSpace displacementTrialSpace = FE_TYPE::template getFunctionSpace< numDofPerTrialSupportPoint >();
    constexpr FunctionSpace pressureTrialSpace = FunctionSpace::P0;
    constexpr FunctionSpace pressureTestSpace = pressureTrialSpace;

    // Step 1: compute local mass balance residual
    LinearFormUtilities::compute< pressureTestSpace,
                                  DifferentialOperator::Identity >
    (
      stack.localResidualMass,
      1.0,
      stack.fluidMassIncrement,
      detJxW );

    // Step 2: compute local mass balance residual derivatives with respect to displacement
    BilinearFormUtilities::compute< pressureTestSpace,
                                    displacementTrialSpace,
                                    DifferentialOperator::Identity,
                                    DifferentialOperator::Divergence >
    (
      stack.dLocalResidualMass_dDisplacement,
      1.0,
      stack.dFluidMassIncrement_dVolStrainIncrement,
      dNdX,
      detJxW );

    // Step 3: compute local mass balance residual derivatives with respect to pressure
    BilinearFormUtilities::compute< pressureTestSpace,
                                    pressureTrialSpace,
                                    DifferentialOperator::Identity,
                                    DifferentialOperator::Identity >
    (
      stack.dLocalResidualMass_dPressure,
      1.0,
      stack.dFluidMassIncrement_dPressure,
      1.0,
      detJxW );
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
    //   dFluidMass_dTime + divergence( fluidMassFlux ) = source  (fluid phase mass balance)
    //
    // with currently the following dependencies on the strainIncrement tensor and pressure
    //
    //   totalStress      = totalStress( strainIncrement, pressure)
    //   bodyForce        = bodyForce( strainIncrement, pressure)
    //   fluidMass = fluidMass( strainIncrement, pressure)
    //   fluidMassFlux    = fludiMassFlux( pressure)
    //
    // Note that the fluidMassFlux will depend on the strainIncrement if a stress-dependent constitutive
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
    //   Rmas = \int \chi ( fluidMass - fluidMass_n) = 0
    //
    //   dRmom_dVolStrain = - \int_Omega symmetricGradient( \eta ) : dTotalStress_dVolStrain
    //                      + \int \eta \cdot dBodyForce_dVolStrain
    //   dRmom_dPressure  = - \int_Omega symmetricGradient( \eta ) : dTotalStress_dPressure
    //                      + \int \eta \cdot dBodyForce_dPressure
    //   dRmas_dVolStrain = \int \chi dFluidMass_dVolStrain
    //   dRmas_dPressure  = \int \chi dFluidMass_dPressure
    //
    // with \eta and \chi test basis functions for the displacement and pressure field, respectively.
    // A continuous interpolation is used for the displacement, with \eta continuous finite element
    // basis functions. A piecewise-constant approximation is used for the pressure.

    // Step 1: compute displacement finite element basis functions (N), basis function derivatives (dNdX), and
    // determinant of the Jacobian transformation matrix times the quadrature weight (detJxW)
    real64 N[numNodesPerElem];
    real64 dNdX[numNodesPerElem][3];
    FE_TYPE::calcN( q, N );
    real64 const detJxW = m_finiteElementSpace.template getGradN< FE_TYPE >( k, q, stack.xLocal, dNdX );

    // Step 2: compute strain increment
    real64 strainIncrement[6]{};
    FE_TYPE::symmetricGradient( dNdX, stack.uhat_local, strainIncrement );

    // Step 3: compute 1) the total stress, 2) the body force terms, and 3) the fluidMassIncrement
    // using quantities returned by the PorousSolid constitutive model.
    // This function also computes the derivatives of these three quantities wrt primary variables
    smallStrainUpdate( k, q, strainIncrement, stack );

    // Step 4: use the total stress and the body force to increment the local momentum balance residual
    // This function also fills the local Jacobian rows corresponding to the momentum balance.
    assembleMomentumBalanceResidual( N, dNdX, detJxW, stack );

    // Step 5: use the fluid mass increment to increment the local mass balance residual
    // This function also fills the local Jacobian rows corresponding to the mass balance.
    assembleElementBasedFlowResidual( dNdX, detJxW, stack );
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
  arrayView2d< real64 const > const m_fluidDensity_n;
  arrayView2d< real64 const > const m_dFluidDensity_dPressure;

  /// The global degree of freedom number
  arrayView1d< globalIndex const > const m_flowDofNumber;

  /// The rank-global fluid pressure arrays.
  arrayView1d< real64 const > const m_fluidPressure_n;
  arrayView1d< real64 const > const m_fluidPressure;

};

using SinglePhaseKernelFactory = finiteElement::KernelFactory< SinglePhase,
                                                               arrayView1d< globalIndex const > const,
                                                               string const,
                                                               globalIndex const,
                                                               CRSMatrixView< real64, globalIndex const > const,
                                                               arrayView1d< real64 > const,
                                                               real64 const (&)[3],
                                                               string const >;

} // namespace poromechanicsKernels

} // namespace geosx

#include "finiteElement/kernelInterface/SparsityKernelBase.hpp"

#endif // GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASEPOROMECHANICSKERNEL_HPP_
