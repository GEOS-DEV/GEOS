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
 * @file ThermalSinglePhasePoroelasticKernel.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_THERMALSINGLEPHASEPOROMECHANICSKERNEL_HPP_
#define GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_THERMALSINGLEPHASEPOROMECHANICSKERNEL_HPP_

#include "physicsSolvers/multiphysics/SinglePhasePoromechanicsKernel.hpp"

namespace geosx
{

namespace thermalPoromechanicsKernels
{

/**
 * @brief Implements kernels for solving quasi-static thermal single-phase poromechanics.
 * @copydoc geosx::finiteElement::ImplicitKernelBase
 * @tparam NUM_NODES_PER_ELEM The number of nodes per element for the
 *                            @p SUBREGION_TYPE.
 * @tparam UNUSED An unused parameter since we are assuming that the test and
 *                trial space have the same number of support points.
 *
 */
template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
class ThermalSinglePhasePoromechanicsKernel :
  public poromechanicsKernels::SinglePhasePoromechanicsKernel< SUBREGION_TYPE,
                                                               CONSTITUTIVE_TYPE,
                                                               FE_TYPE >
{
public:

  /// Alias for the base class;
  using Base = poromechanicsKernels::SinglePhasePoromechanicsKernel< SUBREGION_TYPE,
                                                                     CONSTITUTIVE_TYPE,
                                                                     FE_TYPE >;

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
  using Base::m_gravityAcceleration;
  using Base::m_gravityVector;
  using Base::m_elemsToNodes;
  using Base::m_constitutiveUpdate;
  using Base::m_finiteElementSpace;
  using Base::m_pressure;
  using Base::m_pressure_n;
  using Base::m_fluidDensity;
  using Base::m_fluidDensity_n;
  using Base::m_dFluidDensity_dPressure;
  using Base::m_solidDensity;
  using Base::m_flowDofNumber;

  /**
   * @brief Constructor
   * @copydoc geosx::finiteElement::ImplicitKernelBase::ImplicitKernelBase
   * @param inputGravityVector The gravity vector.
   */
  ThermalSinglePhasePoromechanicsKernel( NodeManager const & nodeManager,
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
          inputFlowDofKey,
          rankOffset,
          inputMatrix,
          inputRhs,
          inputGravityVector,
          fluidModelKey ),
    m_dFluidDensity_dTemperature( elementSubRegion.template getConstitutiveModel< constitutive::SingleFluidBase >( elementSubRegion.template getReference< string >(
                                                                                                                     fluidModelKey ) ).dDensity_dTemperature() ),
    m_fluidInternalEnergy_n( elementSubRegion.template getConstitutiveModel< constitutive::SingleFluidBase >( elementSubRegion.template getReference< string >( fluidModelKey ) ).internalEnergy_n() ),
    m_fluidInternalEnergy( elementSubRegion.template getConstitutiveModel< constitutive::SingleFluidBase >( elementSubRegion.template getReference< string >( fluidModelKey ) ).internalEnergy() ),
    m_dFluidInternalEnergy_dPressure( elementSubRegion.template getConstitutiveModel< constitutive::SingleFluidBase >( elementSubRegion.template getReference< string >(
                                                                                                                         fluidModelKey ) ).dInternalEnergy_dPressure() ),
    m_dFluidInternalEnergy_dTemperature( elementSubRegion.template getConstitutiveModel< constitutive::SingleFluidBase >( elementSubRegion.template getReference< string >(
                                                                                                                            fluidModelKey ) ).dInternalEnergy_dTemperature() ),
    m_temperature_n( elementSubRegion.template getField< fields::flow::temperature_n >() ),
    m_initialTemperature( elementSubRegion.template getField< fields::flow::initialTemperature >() ),
    m_temperature( elementSubRegion.template getField< fields::flow::temperature >() )
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
            dLocalResidualMomentum_dTemperature{ {0.0} },
      dLocalResidualMass_dTemperature{ {0.0} },
      localTemperatureDofIndex{}
    {}

    // Storage for helper variables used in the quadrature point kernel

    /// Derivative of body force wrt temperature
    real64 dBodyForce_dTemperature[3]{};
    /// Derivative of mass accumulation wrt temperature
    real64 dFluidMassIncrement_dTemperature{};

    /// Energy accumulation
    real64 fluidEnergyIncrement{};
    /// Derivative of energy accumulation wrt volumetric strain increment
    real64 dFluidEnergyIncrement_dVolStrainIncrement{};
    /// Derivative of energy accumulation wrt pressure
    real64 dFluidEnergyIncrement_dPressure{};
    /// Derivative of energy accumulation wrt temperature
    real64 dFluidEnergyIncrement_dTemperature{};

    // Storage for mass/energy residual and degrees of freedom

    /// Derivative of linear momentum balance residual wrt temperature
    real64 dLocalResidualMomentum_dTemperature[numDispDofPerElem][1]{};
    /// Derivative of mass balance residual wrt pressure
    real64 dLocalResidualMass_dTemperature[1][1]{};

    /// Energy balance residual
    real64 localResidualEnergy[1]{};
    /// Derivative of energy balance residual wrt displacement
    real64 dLocalResidualEnergy_dDisplacement[1][numDispDofPerElem]{};
    /// Derivative of energy balance residual wrt pressure
    real64 dLocalResidualEnergy_dPressure[1][1]{};
    /// Derivative of energy balance residual wrt temperature
    real64 dLocalResidualEnergy_dTemperature[1][1]{};

    /// C-array storage for the element local row degrees of freedom.
    globalIndex localTemperatureDofIndex{};

  };
  //*****************************************************************************

  /**
   * @brief Copy global values from primary field to a local stack array.
   * @copydoc ::geosx::finiteElement::ImplicitKernelBase::setup
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void setup( localIndex const k,
              StackVariables & stack ) const
  {
    Base::setup( k, stack );
    stack.localTemperatureDofIndex = m_flowDofNumber[k]+1;
    stack.deltaTemperatureFromInit = m_temperature[k] - m_initialTemperature[k];
    stack.deltaTemperatureFromLastStep = m_temperature[k] - m_temperature_n[k];
  }

  /**
   * @brief Helper function to compute 1) the total stress, 2) the body force term, and 3) the fluidMass/EnergyIncrement
   * using quantities returned by the PorousSolid constitutive model.
   * This function also computes the derivatives of these three quantities wrt primary variables
   * @param[in] k the element index
   * @param[in] q the quadrature point index
   * @param[in] strainIncrement the strain increment used in total stress and porosity computation
   * @param[inout] stack the stack variables
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void
  smallStrainUpdate( localIndex const k,
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
                                                         m_pressure_n[k],
                                                         m_pressure[k],
                                                         stack.deltaTemperatureFromInit,
                                                         stack.deltaTemperatureFromLastStep,
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
    if( m_gravityAcceleration > 0.0 )
    {
      computeBodyForce( k, q,
                        porosity,
                        dPorosity_dVolStrain,
                        dPorosity_dPressure,
                        dPorosity_dTemperature,
                        dSolidDensity_dPressure,
                        stack );
    }

    // Step 3: compute fluid mass increment
    computeFluidIncrement( k, q,
                           porosity,
                           porosity_n,
                           dPorosity_dVolStrain,
                           dPorosity_dPressure,
                           dPorosity_dTemperature,
                           stack );
  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void
  computeBodyForce( localIndex const k,
                    localIndex const q,
                    real64 const & porosity,
                    real64 const & dPorosity_dVolStrain,
                    real64 const & dPorosity_dPressure,
                    real64 const & dPorosity_dTemperature,
                    real64 const & dSolidDensity_dPressure,
                    StackVariables & stack ) const
  {
    Base::computeBodyForce( k, q,
                            porosity,
                            dPorosity_dVolStrain,
                            dPorosity_dPressure,
                            dPorosity_dTemperature,
                            dSolidDensity_dPressure,
                            stack );

    real64 const dMixtureDens_dTemperature =
      dPorosity_dTemperature * ( -m_solidDensity( k, q ) + m_fluidDensity( k, q ) )
      + porosity * m_dFluidDensity_dTemperature( k, q );

    LvArray::tensorOps::scaledCopy< 3 >( stack.dBodyForce_dTemperature, m_gravityVector, dMixtureDens_dTemperature );
  }

  /**
   * @brief Helper function to compute the fluid mass/energy increment and its derivatives wrt primary variables
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
  void
  computeFluidIncrement( localIndex const k,
                         localIndex const q,
                         real64 const & porosity,
                         real64 const & porosity_n,
                         real64 const & dPorosity_dVolStrain,
                         real64 const & dPorosity_dPressure,
                         real64 const & dPorosity_dTemperature,
                         StackVariables & stack ) const
  {
    // Step 1: compute fluid mass increment and its derivatives wrt vol strain and pressure
    Base::computeFluidIncrement( k, q,
                                 porosity,
                                 porosity_n,
                                 dPorosity_dVolStrain,
                                 dPorosity_dPressure,
                                 dPorosity_dTemperature,
                                 stack );

    // Step 2: compute derivative of fluid mass increment wrt temperature
    stack.dFluidMassIncrement_dTemperature = dPorosity_dTemperature * m_fluidDensity( k, q ) + porosity * m_dFluidDensity_dTemperature( k, q );

    // Step 3: compute fluid energy increment and its derivatives wrt vol strain, pressure, and temperature
    real64 const fluidMass = porosity * m_fluidDensity( k, q );
    real64 const fluidEnergy = fluidMass * m_fluidInternalEnergy( k, q );
    real64 const fluidEnergy_n = porosity_n * m_fluidDensity_n( k, q ) * m_fluidInternalEnergy_n( k, q );
    stack.fluidEnergyIncrement = fluidEnergy - fluidEnergy_n;

    stack.dFluidEnergyIncrement_dVolStrainIncrement = stack.dFluidMassIncrement_dVolStrainIncrement * m_fluidInternalEnergy( k, q );
    stack.dFluidEnergyIncrement_dPressure = stack.dFluidMassIncrement_dPressure * m_fluidInternalEnergy( k, q )
                                            + fluidMass * m_dFluidInternalEnergy_dPressure( k, q );
    stack.dFluidEnergyIncrement_dTemperature = stack.dFluidMassIncrement_dTemperature * m_fluidInternalEnergy( k, q )
                                               + fluidMass * m_dFluidInternalEnergy_dTemperature( k, q );
  }

  /**
   * @brief Assemble the local linear momentum balance residual and derivatives using total stress and body force terms
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
  void
  assembleMomentumBalanceTerms( real64 const ( &N )[numNodesPerElem],
                                real64 const ( &dNdX )[numNodesPerElem][3],
                                real64 const & detJxW,
                                StackVariables & stack ) const
  {
    using namespace PDEUtilities;

    Base::assembleMomentumBalanceTerms( N, dNdX, detJxW, stack );

    constexpr FunctionSpace displacementTrialSpace = FE_TYPE::template getFunctionSpace< numDofPerTrialSupportPoint >();
    constexpr FunctionSpace displacementTestSpace = displacementTrialSpace;
    constexpr FunctionSpace pressureTrialSpace = FunctionSpace::P0;

    // compute local linear momentum balance residual derivatives with respect to temperature

    BilinearFormUtilities::compute< displacementTestSpace,
                                    pressureTrialSpace,
                                    DifferentialOperator::SymmetricGradient,
                                    DifferentialOperator::Identity >
    (
      stack.dLocalResidualMomentum_dTemperature,
      dNdX,
      stack.dTotalStress_dTemperature,
      1.0,
      -detJxW );

    if( m_gravityAcceleration > 0.0 )
    {
      BilinearFormUtilities::compute< displacementTestSpace,
                                      pressureTrialSpace,
                                      DifferentialOperator::Identity,
                                      DifferentialOperator::Identity >
      (
        stack.dLocalResidualMomentum_dTemperature,
        N,
        stack.dBodyForce_dTemperature,
        1.0,
        detJxW );
    }
  }

  /**
   * @brief Assemble the local mass/energy balance residual and derivatives using fluid mass/energy increment
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
  void
  assembleElementBasedFlowTerms( real64 const ( &dNdX )[numNodesPerElem][3],
                                 real64 const & detJxW,
                                 StackVariables & stack ) const
  {

    using namespace PDEUtilities;

    Base::assembleElementBasedFlowTerms( dNdX, detJxW, stack );

    constexpr FunctionSpace displacementTrialSpace = FE_TYPE::template getFunctionSpace< numDofPerTrialSupportPoint >();
    constexpr FunctionSpace pressureTrialSpace = FunctionSpace::P0;
    constexpr FunctionSpace pressureTestSpace = pressureTrialSpace;

    // Step 1: compute local mass balance residual derivatives with respect to temperature

    BilinearFormUtilities::compute< pressureTestSpace,
                                    pressureTrialSpace,
                                    DifferentialOperator::Identity,
                                    DifferentialOperator::Identity >
    (
      stack.dLocalResidualMass_dTemperature,
      1.0,
      stack.dFluidMassIncrement_dTemperature,
      1.0,
      detJxW );

    // Step 2: compute local energy balance residual and its derivatives

    LinearFormUtilities::compute< pressureTestSpace,
                                  DifferentialOperator::Identity >
    (
      stack.localResidualEnergy,
      1.0,
      stack.fluidEnergyIncrement,
      detJxW );

    BilinearFormUtilities::compute< pressureTestSpace,
                                    displacementTrialSpace,
                                    DifferentialOperator::Identity,
                                    DifferentialOperator::Divergence >
    (
      stack.dLocalResidualEnergy_dDisplacement,
      1.0,
      stack.dFluidEnergyIncrement_dVolStrainIncrement,
      dNdX,
      detJxW );

    BilinearFormUtilities::compute< pressureTestSpace,
                                    pressureTrialSpace,
                                    DifferentialOperator::Identity,
                                    DifferentialOperator::Identity >
    (
      stack.dLocalResidualEnergy_dPressure,
      1.0,
      stack.dFluidEnergyIncrement_dPressure,
      1.0,
      detJxW );

    BilinearFormUtilities::compute< pressureTestSpace,
                                    pressureTrialSpace,
                                    DifferentialOperator::Identity,
                                    DifferentialOperator::Identity >
    (
      stack.dLocalResidualEnergy_dTemperature,
      1.0,
      stack.dFluidEnergyIncrement_dTemperature,
      1.0,
      detJxW );
  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void quadraturePointKernel( localIndex const k,
                              localIndex const q,
                              StackVariables & stack ) const
  {
    // Step 1: compute displacement finite element basis functions (N), basis function derivatives (dNdX), and
    // determinant of the Jacobian transformation matrix times the quadrature weight (detJxW)
    real64 N[numNodesPerElem]{};
    real64 dNdX[numNodesPerElem][3]{};
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
    assembleMomentumBalanceTerms( N, dNdX, detJxW, stack );

    // Step 5: use the fluid mass increment to increment the local mass balance residual
    // This function also fills the local Jacobian rows corresponding to the mass balance.
    assembleElementBasedFlowTerms( dNdX, detJxW, stack );
  }

  /**
   * @copydoc geosx::finiteElement::ImplicitKernelBase::complete
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  real64 complete( localIndex const k,
                   StackVariables & stack ) const
  {
    real64 const maxForce = Base::complete( k, stack );

    constexpr integer nUDof = numNodesPerElem * numDofPerTestSupportPoint;

    // Step 1: assemble the derivatives of linear momentum balance wrt temperature into the global matrix

    for( integer localNode = 0; localNode < numNodesPerElem; ++localNode )
    {
      for( integer dim = 0; dim < numDofPerTestSupportPoint; ++dim )
      {
        localIndex const dof = LvArray::integerConversion< localIndex >( stack.localRowDofIndex[numDofPerTestSupportPoint*localNode + dim] - m_dofRankOffset );
        if( dof < 0 || dof >= m_matrix.numRows() ) continue;

        m_matrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( dof,
                                                                                &stack.localTemperatureDofIndex,
                                                                                stack.dLocalResidualMomentum_dTemperature[numDofPerTestSupportPoint * localNode + dim],
                                                                                1 );
      }
    }

    // Step 2: assemble the derivatives of mass balance residual wrt temperature into the global matrix

    localIndex const massDof = LvArray::integerConversion< localIndex >( stack.localPressureDofIndex - m_dofRankOffset );
    if( 0 <= massDof && massDof < m_matrix.numRows() )
    {
      m_matrix.template addToRow< serialAtomic >( massDof,
                                                  &stack.localTemperatureDofIndex,
                                                  stack.dLocalResidualMass_dTemperature[0],
                                                  1 );
    }

    // Step 3: assemble the energy balance and its derivatives into the global matrix

    localIndex const energyDof = LvArray::integerConversion< localIndex >( stack.localTemperatureDofIndex - m_dofRankOffset );
    if( 0 <= energyDof && energyDof < m_matrix.numRows() )
    {
      m_matrix.template addToRowBinarySearchUnsorted< serialAtomic >( energyDof,
                                                                      stack.localRowDofIndex,
                                                                      stack.dLocalResidualEnergy_dDisplacement[0],
                                                                      nUDof );
      m_matrix.template addToRow< serialAtomic >( energyDof,
                                                  &stack.localPressureDofIndex,
                                                  stack.dLocalResidualEnergy_dPressure[0],
                                                  1 );
      m_matrix.template addToRow< serialAtomic >( energyDof,
                                                  &stack.localTemperatureDofIndex,
                                                  stack.dLocalResidualEnergy_dTemperature[0],
                                                  1 );

      RAJA::atomicAdd< serialAtomic >( &m_rhs[energyDof], stack.localResidualEnergy[0] );
    }

    return maxForce;
  }

protected:

  /// The rank global density derivative wrt temperature
  arrayView2d< real64 const > const m_dFluidDensity_dTemperature;

  /// The rank global internal energy
  arrayView2d< real64 const > const m_fluidInternalEnergy_n;
  arrayView2d< real64 const > const m_fluidInternalEnergy;
  arrayView2d< real64 const > const m_dFluidInternalEnergy_dPressure;
  arrayView2d< real64 const > const m_dFluidInternalEnergy_dTemperature;

  /// The rank-global fluid temperature arrays.
  arrayView1d< real64 const > const m_temperature_n;
  arrayView1d< real64 const > const m_initialTemperature;
  arrayView1d< real64 const > const m_temperature;

};

using ThermalSinglePhasePoromechanicsKernelFactory =
  finiteElement::KernelFactory< ThermalSinglePhasePoromechanicsKernel,
                                arrayView1d< globalIndex const > const,
                                string const,
                                globalIndex const,
                                CRSMatrixView< real64, globalIndex const > const,
                                arrayView1d< real64 > const,
                                real64 const (&)[3],
                                string const >;

} // namespace thermalPoromechanicsKernels

} // namespace geosx

#include "finiteElement/kernelInterface/SparsityKernelBase.hpp"

#endif // GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_THERMALSINGLEPHASEPOROMECHANICSKERNEL_HPP_
