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
 * @file ThermalMultiphasePoromechanicsKernel.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_THERMALMULTIPHASEPOROMECHANICSKERNEL_HPP_
#define GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_THERMALMULTIPHASEPOROMECHANICSKERNEL_HPP_

#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBaseFields.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseUtilities.hpp"
#include "physicsSolvers/multiphysics/PoromechanicsKernelBase.hpp"

namespace geosx
{

namespace thermalPoromechanicsKernels
{

/**
 * @brief Implements kernels for solving quasi-static thermal multiphase poromechanics.
 * @copydoc geosx::finiteElement::ImplicitKernelBase
 * @tparam NUM_NODES_PER_ELEM The number of nodes per element for the
 *                            @p SUBREGION_TYPE.
 * @tparam UNUSED An unused parameter since we are assuming that the test and
 *                trial space have the same number of support points.
 *
 * ### ThermalMultiphasePoromechanics Description
 * Implements the KernelBase interface functions required for solving the
 * quasi-static multiphase poromechanics problem using one of the
 * "finite element kernel application" functions such as
 * geosx::finiteElement::RegionBasedKernelApplication.
 *
 */
template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
class ThermalMultiphasePoromechanicsKernel :
  public poromechanicsKernels::MultiphasePoromechanicsKernel< SUBREGION_TYPE,
                                                              CONSTITUTIVE_TYPE,
                                                              FE_TYPE >
{
public:

  /// Alias for the base class;
  using Base = poromechanicsKernels::MultiphasePoromechanicsKernel< SUBREGION_TYPE,
                                                                    CONSTITUTIVE_TYPE,
                                                                    FE_TYPE >;


  /// Maximum number of nodes per element, which is equal to the maxNumTestSupportPointPerElem and
  /// maxNumTrialSupportPointPerElem by definition. When the FE_TYPE is not a Virtual Element, this
  /// will be the actual number of nodes per element.
  static constexpr int numNodesPerElem = Base::maxNumTestSupportPointsPerElem;
  static constexpr int maxNumComponents = 3;
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
  using Base::m_numComponents;
  using Base::m_numPhases;
  using Base::m_solidDensity;
  using Base::m_fluidPhaseMassDensity;
  using Base::m_dFluidPhaseMassDensity;
  using Base::m_fluidPhaseDensity;
  using Base::m_dFluidPhaseDensity;
  using Base::m_fluidPhaseVolFrac;
  using Base::m_dFluidPhaseVolFrac;
  using Base::m_fluidPhaseCompFrac;
  using Base::m_dFluidPhaseCompFrac;
  using Base::m_dGlobalCompFraction_dGlobalCompDensity;

  /**
   * @brief Constructor
   * @copydoc geosx::finiteElement::ImplicitKernelBase::ImplicitKernelBase
   * @param inputGravityVector The gravity vector.
   */
  ThermalMultiphasePoromechanicsKernel( NodeManager const & nodeManager,
                                        EdgeManager const & edgeManager,
                                        FaceManager const & faceManager,
                                        localIndex const targetRegionIndex,
                                        SUBREGION_TYPE const & elementSubRegion,
                                        FE_TYPE const & finiteElementSpace,
                                        CONSTITUTIVE_TYPE & inputConstitutiveType,
                                        arrayView1d< globalIndex const > const inputDispDofNumber,
                                        string const inputFlowDofKey,
                                        globalIndex const rankOffset,
                                        localIndex const numComponents,
                                        localIndex const numPhases,
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
          numComponents,
          numPhases,
          inputMatrix,
          inputRhs,
          inputGravityVector,
          fluidModelKey ),
    m_rockInternalEnergy_n( inputConstitutiveType.getInternalEnergy_n() ),
    m_rockInternalEnergy( inputConstitutiveType.getInternalEnergy() ),
    m_dRockInternalEnergy_dTemperature( inputConstitutiveType.getDinternalEnergy_dTemperature() ),
    m_temperature_n( elementSubRegion.template getField< fields::flow::temperature_n >() ),
    m_temperature( elementSubRegion.template getField< fields::flow::temperature >() )
  {
    // extract fluid constitutive data views
    {
      string const fluidModelName = elementSubRegion.template getReference< string >( fluidModelKey );
      constitutive::MultiFluidBase const & fluid =
        elementSubRegion.template getConstitutiveModel< constitutive::MultiFluidBase >( fluidModelName );

      m_fluidPhaseInternalEnergy_n = fluid.phaseInternalEnergy_n();
      m_fluidPhaseInternalEnergy = fluid.phaseInternalEnergy();
      m_dFluidPhaseInternalEnergy = fluid.dPhaseInternalEnergy();
    }
  }

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

    /// Constructor.
    GEOSX_HOST_DEVICE
    StackVariables():
      Base::StackVariables()
    {}

    /// Derivative of mass accumulation wrt temperature
    real64 dCompMassIncrement_dTemperature[maxNumComponents]{};
    /// Derivative of pore volume constraint wrt temperature
    real64 dPoreVolConstraint_dTemperature;
    /// Derivative of body force wrt temperature
    real64 dBodyForce_dTemperature[3]{};

    /// Energy accumulation
    real64 fluidEnergyIncrement{};
    /// Derivative of energy accumulation wrt volumetric strain increment
    real64 dFluidEnergyIncrement_dVolStrainIncrement{};
    /// Derivative of energy accumulation wrt pressure
    real64 dFluidEnergyIncrement_dPressure{};
    /// Derivative of energy accumulation wrt temperature
    real64 dFluidEnergyIncrement_dTemperature{};
    /// Derivative of energy accumulation wrt temperature
    real64 dFluidEnergyIncrement_dComponents[maxNumComponents]{};

    // Storage for mass/energy residual and degrees of freedom

    /// Derivative of linear momentum balance residual wrt temperature
    real64 dLocalResidualMomentum_dTemperature[Base::StackVariables::numDispDofPerElem][1]{};
    /// Derivative of mass balance residual wrt pressure
    real64 dLocalResidualMass_dTemperature[maxNumComponents][1]{};
    /// Derivative of pore volume constraint residual wrt pressure
    real64 dLocalResidualPoreVolConstraint_dTemperature[1][1]{};

    /// Energy balance residual
    real64 localResidualEnergy[1]{};
    /// Derivative of energy balance residual wrt displacement
    real64 dLocalResidualEnergy_dDisplacement[1][Base::StackVariables::numDispDofPerElem]{};
    /// Derivative of energy balance residual wrt pressure
    real64 dLocalResidualEnergy_dPressure[1][1]{};
    /// Derivative of energy balance residual wrt temperature
    real64 dLocalResidualEnergy_dTemperature[1][1]{};
    /// Derivative of energy balance residual wrt components
    real64 dLocalResidualEnergy_dComponents[1][maxNumComponents]{};

    /// C-array storage for the element local row degrees of freedom.
    globalIndex localTemperatureDofIndex[1]{};

  };
  //*****************************************************************************

  /**
   * @brief Copy global values from primary field to a local stack array.
   * @copydoc ::geosx::finiteElement::ImplicitKernelBase::setup
   *
   * For the ThermalMultiphase implementation, global values from the displacement,
   * incremental displacement, and degree of freedom numbers are placed into
   * element local stack storage.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void setup( localIndex const k,
              StackVariables & stack ) const
  {
    // initialize displacement dof and pressure/components dofs
    Base::setup( k, stack );

    stack.localTemperatureDofIndex[0] = stack.localPressureDofIndex[0] + m_numComponents + 1;
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

    // Step 4: compute pore volume constraint
    computePoreVolumeConstraint( k,
                                 porosity,
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
    using Deriv = constitutive::multifluid::DerivativeOffset;

    Base::computeBodyForce( k, q,
                            porosity,
                            dPorosity_dVolStrain,
                            dPorosity_dPressure,
                            dPorosity_dTemperature,
                            dSolidDensity_dPressure,
                            stack, [&]( real64 const & totalMassDensity,
                                        real64 const & mixtureDensity )
    {
      GEOSX_UNUSED_VAR( mixtureDensity );

      arraySlice1d< real64 const, constitutive::multifluid::USD_PHASE - 2 > const phaseMassDensity = m_fluidPhaseMassDensity[k][q];
      arraySlice2d< real64 const, constitutive::multifluid::USD_PHASE_DC - 2 > const dPhaseMassDensity = m_dFluidPhaseMassDensity[k][q];
      arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const phaseVolFrac = m_fluidPhaseVolFrac[k];
      arraySlice2d< real64 const, compflow::USD_PHASE_DC - 1 > const dPhaseVolFrac = m_dFluidPhaseVolFrac[k];

      // Step 1: compute fluid total mass density and its derivatives

      real64 dTotalMassDensity_dTemperature = 0.0;

      for( integer ip = 0; ip < m_numPhases; ++ip )
      {
        dTotalMassDensity_dTemperature = dTotalMassDensity_dTemperature
                                         + dPhaseVolFrac( ip, Deriv::dT ) * phaseMassDensity( ip )
                                         + phaseVolFrac( ip ) * dPhaseMassDensity( ip, Deriv::dT );
      }

      // Step 2: compute mixture density as an average between total mass density and solid density

      real64 const dMixtureDens_dTemperature = dPorosity_dTemperature * ( -m_solidDensity( k, q ) + totalMassDensity )
                                               + porosity * dTotalMassDensity_dTemperature;

      // Step 3: finally, get the body force

      LvArray::tensorOps::scaledCopy< 3 >( stack.dBodyForce_dPressure, m_gravityVector, dMixtureDens_dTemperature );

    } );
  }

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
    using Deriv = constitutive::multifluid::DerivativeOffset;

    Base::computeFluidIncrement( k, q,
                                 porosity,
                                 porosity_n,
                                 dPorosity_dVolStrain,
                                 dPorosity_dPressure,
                                 dPorosity_dTemperature,
                                 stack, [&]( real64 const & ip,
                                             real64 const & phaseAmount,
                                             real64 const & phaseAmount_n,
                                             real64 const & dPhaseAmount_dP,
                                             real64 const (&dPhaseAmount_dC)[maxNumComponents] )
    {

      // We are in the loop over phases, ip provides the current phase index.
      // We have to do two things:
      //   1- Assemble the derivatives of the component mass balance equations with respect to temperature
      //   2- Assemble the phase-dependent part of the accumulation term of the energy equation

      real64 dPhaseInternalEnergy_dC[maxNumComponents]{};

      // construct the slices
      arraySlice1d< real64 const, constitutive::multifluid::USD_PHASE - 2 > const phaseInternalEnergy_n = m_fluidPhaseInternalEnergy_n[k][q];
      arraySlice1d< real64 const, constitutive::multifluid::USD_PHASE - 2 > const phaseInternalEnergy = m_fluidPhaseInternalEnergy[k][q];
      arraySlice2d< real64 const, constitutive::multifluid::USD_PHASE_DC - 2 > const dPhaseInternalEnergy = m_dFluidPhaseInternalEnergy[k][q];
      arraySlice1d< real64 const, constitutive::multifluid::USD_PHASE - 2 > const phaseDensity = m_fluidPhaseDensity[k][q];
      arraySlice2d< real64 const, constitutive::multifluid::USD_PHASE_DC - 2 > const dPhaseDensity = m_dFluidPhaseDensity[k][q];
      arraySlice2d< real64 const, constitutive::multifluid::USD_PHASE_COMP - 2 > const phaseCompFrac = m_fluidPhaseCompFrac[k][q];
      arraySlice3d< real64 const, constitutive::multifluid::USD_PHASE_COMP_DC -2 > const dPhaseCompFrac = m_dFluidPhaseCompFrac[k][q];
      arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const phaseVolFrac = m_fluidPhaseVolFrac[k];
      arraySlice2d< real64 const, compflow::USD_PHASE_DC - 1 > const dPhaseVolFrac = m_dFluidPhaseVolFrac[k];
      arraySlice2d< real64 const, compflow::USD_COMP_DC - 1 > const dGlobalCompFrac_dGlobalCompDensity = m_dGlobalCompFraction_dGlobalCompDensity[k];

      // Step 1: assemble the derivatives of the component mass balance equations with respect to temperature

      real64 const dPhaseAmount_dT = dPorosity_dTemperature * phaseVolFrac( ip ) * phaseDensity( ip )
                                     + porosity * ( dPhaseVolFrac( ip, Deriv::dT ) * phaseDensity( ip ) + phaseVolFrac( ip ) * dPhaseDensity( ip, Deriv::dT ) );

      for( integer ic = 0; ic < m_numComponents; ++ic )
      {
        stack.dCompMassIncrement_dTemperature[ic] += dPhaseAmount_dT * phaseCompFrac( ip, ic )
                                                     + phaseAmount * dPhaseCompFrac( ip, ic, Deriv::dT );
      }

      // Step 2: assemble the phase-dependent part of the accumulation term of the energy equation

      // local accumulation
      stack.fluidEnergyIncrement += phaseAmount * phaseInternalEnergy( ip ) - phaseAmount_n * phaseInternalEnergy_n( ip );

      // derivatives w.r.t. pressure and temperature
      stack.dFluidEnergyIncrement_dPressure += dPhaseAmount_dP * phaseInternalEnergy( ip )
                                               + phaseAmount * dPhaseInternalEnergy( ip, Deriv::dP );
      stack.dFluidEnergyIncrement_dTemperature += dPhaseAmount_dT * phaseInternalEnergy( ip )
                                                  + phaseAmount * dPhaseInternalEnergy( ip, Deriv::dT );

      // derivatives w.r.t. component densities
      applyChainRule( m_numComponents, dGlobalCompFrac_dGlobalCompDensity, dPhaseInternalEnergy[ip], dPhaseInternalEnergy_dC, Deriv::dC );
      for( integer jc = 0; jc < m_numComponents; ++jc )
      {
        stack.dFluidEnergyIncrement_dComponents[jc] += dPhaseAmount_dC[jc] * phaseInternalEnergy( ip )
                                                       + phaseAmount * dPhaseInternalEnergy_dC[jc];
      }
    } );

    // Step 3: assemble the solid part of the accumulation term

    real64 const oneMinusPoro = 1 - porosity;

    // local accumulation and derivatives w.r.t. pressure and temperature
    stack.fluidEnergyIncrement += oneMinusPoro * m_rockInternalEnergy( k, q ) - ( 1 - porosity_n ) * m_rockInternalEnergy_n( k, q );
    stack.dFluidEnergyIncrement_dPressure += -dPorosity_dPressure * m_rockInternalEnergy( k, q );
    stack.dFluidEnergyIncrement_dTemperature += -dPorosity_dTemperature * m_rockInternalEnergy( k, q ) + oneMinusPoro * m_dRockInternalEnergy_dTemperature( k, q );

  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void
  computePoreVolumeConstraint( localIndex const k,
                               real64 const & porosity,
                               real64 const & dPorosity_dVolStrain,
                               real64 const & dPorosity_dPressure,
                               real64 const & dPorosity_dTemperature,
                               StackVariables & stack ) const
  {
    using Deriv = constitutive::multifluid::DerivativeOffset;

    Base::computePoreVolumeConstraint( k,
                                       porosity,
                                       dPorosity_dVolStrain,
                                       dPorosity_dPressure,
                                       dPorosity_dTemperature,
                                       stack );

    arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const phaseVolFrac = m_fluidPhaseVolFrac[k];
    arraySlice2d< real64 const, compflow::USD_PHASE_DC - 1 > const dPhaseVolFrac = m_dFluidPhaseVolFrac[k];

    stack.dPoreVolConstraint_dTemperature = 0.0;
    for( integer ip = 0; ip < m_numPhases; ++ip )
    {
      stack.dPoreVolConstraint_dTemperature = stack.dPoreVolConstraint_dTemperature
                                              - dPhaseVolFrac( ip, Deriv::dT ) * porosity
                                              - phaseVolFrac( ip ) * dPorosity_dTemperature;
    }
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
      stack.dCompMassIncrement_dTemperature,
      1.0,
      detJxW );

    // Step 2: compute local pore volume constraint residual derivatives with respect to temperature

    BilinearFormUtilities::compute< pressureTestSpace,
                                    pressureTrialSpace,
                                    DifferentialOperator::Identity,
                                    DifferentialOperator::Identity >
    (
      stack.dLocalResidualPoreVolConstraint_dTemperature,
      1.0,
      stack.dPoreVolConstraint_dTemperature,
      1.0,
      detJxW );

    // Step 3: compute local energy balance residual and its derivatives

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
    using namespace compositionalMultiphaseUtilities;

    real64 const maxForce = Base::complete( k, stack );

    constexpr integer nUDof = numNodesPerElem * numDofPerTestSupportPoint;

    // Apply equation/variable change transformation(s)
    real64 work[maxNumComponents + 1]{};
    shiftRowsAheadByOneAndReplaceFirstRowWithColumnSum( m_numComponents, 1, stack.dLocalResidualMass_dTemperature, work );

    // Step 1: assemble the derivatives of linear momentum balance wrt temperature into the global matrix

    for( integer localNode = 0; localNode < numNodesPerElem; ++localNode )
    {
      for( integer dim = 0; dim < numDofPerTestSupportPoint; ++dim )
      {
        localIndex const dof = LvArray::integerConversion< localIndex >( stack.localRowDofIndex[numDofPerTestSupportPoint*localNode + dim] - m_dofRankOffset );
        if( dof < 0 || dof >= m_matrix.numRows() ) continue;

        m_matrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( dof,
                                                                                stack.localTemperatureDofIndex,
                                                                                stack.dLocalResidualMomentum_dTemperature[numDofPerTestSupportPoint * localNode + dim],
                                                                                1 );
      }
    }

    // Step 2: assemble the derivatives of mass balance residual wrt temperature into the global matrix

    localIndex const massDof = LvArray::integerConversion< localIndex >( stack.localPressureDofIndex[0] - m_dofRankOffset );
    if( 0 <= massDof && massDof < m_matrix.numRows() )
    {
      for( localIndex i = 0; i < m_numComponents; ++i )
      {
        m_matrix.template addToRow< serialAtomic >( massDof + i,
                                                    stack.localTemperatureDofIndex,
                                                    stack.dLocalResidualMass_dTemperature[i],
                                                    1 );
      }
    }

    // Step 3: assemble the derivatives of pore-volume constraint residual wrt temperature into the global matrix

    m_matrix.template addToRow< serialAtomic >( massDof + m_numComponents,
                                                stack.localTemperatureDofIndex,
                                                stack.dLocalResidualPoreVolConstraint_dTemperature[0],
                                                1 );


    // Step 4: assemble the energy balance and its derivatives into the global matrix

    localIndex const energyDof = LvArray::integerConversion< localIndex >( stack.localTemperatureDofIndex[0] - m_dofRankOffset );
    if( 0 <= energyDof && energyDof < m_matrix.numRows() )
    {
      m_matrix.template addToRowBinarySearchUnsorted< serialAtomic >( energyDof,
                                                                      stack.localRowDofIndex,
                                                                      stack.dLocalResidualEnergy_dDisplacement[0],
                                                                      nUDof );
      m_matrix.template addToRow< serialAtomic >( energyDof,
                                                  stack.localPressureDofIndex,
                                                  stack.dLocalResidualEnergy_dPressure[0],
                                                  1 );
      m_matrix.template addToRow< serialAtomic >( energyDof,
                                                  stack.localTemperatureDofIndex,
                                                  stack.dLocalResidualEnergy_dTemperature[0],
                                                  1 );
      m_matrix.template addToRow< serialAtomic >( energyDof,
                                                  stack.localComponentDofIndices,
                                                  stack.dLocalResidualEnergy_dComponents[0],
                                                  m_numComponents );

      RAJA::atomicAdd< serialAtomic >( &m_rhs[energyDof], stack.localResidualEnergy[0] );
    }

    return maxForce;
  }

protected:

  /// Views on phase internal energy
  arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > m_fluidPhaseInternalEnergy_n;
  arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > m_fluidPhaseInternalEnergy;
  arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > m_dFluidPhaseInternalEnergy;

  /// Views on rock internal energy
  arrayView2d< real64 const > const m_rockInternalEnergy_n;
  arrayView2d< real64 const > const m_rockInternalEnergy;
  arrayView2d< real64 const > const m_dRockInternalEnergy_dTemperature;

  /// Views on temperature
  arrayView1d< real64 const > const m_temperature_n;
  arrayView1d< real64 const > const m_initialTemperature;
  arrayView1d< real64 const > const m_temperature;

};

using ThermalMultiphasePoromechanicsKernelFactory =
  finiteElement::KernelFactory< ThermalMultiphasePoromechanicsKernel,
                                arrayView1d< globalIndex const > const,
                                string const,
                                globalIndex const,
                                localIndex const,
                                localIndex const,
                                CRSMatrixView< real64, globalIndex const > const,
                                arrayView1d< real64 > const,
                                real64 const (&)[3],
                                string const >;

} // namespace thermalporomechanicsKernels

} // namespace geosx

#include "finiteElement/kernelInterface/SparsityKernelBase.hpp"

#endif // GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_THERMALMULTIPHASEPOROMECHANICSKERNEL_HPP_
