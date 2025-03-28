/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file ThermalFluxComputeKernel.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEREACTIVE_THERMALFLUXCOMPUTEKERNEL_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEREACTIVE_THERMALFLUXCOMPUTEKERNEL_HPP

#include "physicsSolvers/fluidFlow/kernels/singlePhase/reactive/FluxComputeKernel.hpp"

#include "constitutive/thermalConductivity/SinglePhaseThermalConductivityBase.hpp"
#include "constitutive/thermalConductivity/ThermalConductivityFields.hpp"

namespace geos
{

namespace thermalSinglePhaseReactiveFVMKernels
{
/******************************** FluxComputeKernel ********************************/

/**
 * @class FluxComputeKernel
 * @tparam NUM_SPECIES number of fluid primary species
 * @tparam NUM_EQN number of equations
 * @tparam NUM_DOF number of degrees of freedom
 * @tparam STENCILWRAPPER the type of the stencil wrapper
 * @brief Define the interface for the assembly kernel in charge of flux terms
 */
template< integer NUM_SPECIES, integer NUM_EQN, integer NUM_DOF, typename STENCILWRAPPER >
class FluxComputeKernel : public singlePhaseReactiveFVMKernels::FluxComputeKernel< NUM_SPECIES, NUM_EQN, NUM_DOF, STENCILWRAPPER >
{
public:

  /**
   * @brief The type for element-based data. Consists entirely of ArrayView's.
   *
   * Can be converted from ElementRegionManager::ElementViewConstAccessor
   * by calling .toView() or .toViewConst() on an accessor instance
   */
  template< typename VIEWTYPE >
  using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;

  using AbstractBase = singlePhaseFVMKernels::FluxComputeKernelBase;
  using DofNumberAccessor = AbstractBase::DofNumberAccessor;
  using SinglePhaseFlowAccessors = AbstractBase::SinglePhaseFlowAccessors;
  using SinglePhaseFluidAccessors = AbstractBase::SinglePhaseFluidAccessors;
  using PermeabilityAccessors = AbstractBase::PermeabilityAccessors;

  using AbstractBase::m_dt;
  using AbstractBase::m_rankOffset;
  using AbstractBase::m_dofNumber;
  using AbstractBase::m_gravCoef;
  using AbstractBase::m_mob;
  using AbstractBase::m_dMob;
  using AbstractBase::m_dens;
  using AbstractBase::m_dDens;

  using Base = singlePhaseReactiveFVMKernels::FluxComputeKernel< NUM_SPECIES, NUM_EQN, NUM_DOF, STENCILWRAPPER >;
  using ReactiveSinglePhaseFlowAccessors = typename Base::ReactiveSinglePhaseFlowAccessors;
  using ReactiveSinglePhaseFluidAccessors = typename Base::ReactiveSinglePhaseFluidAccessors;
  using DiffusionAccessors = typename Base::DiffusionAccessors;
  using PorosityAccessors = typename Base::PorosityAccessors;
  using Base::numSpecies;
  using Base::numFluxSupportPoints;
  using Base::numDof;
  using Base::numEqn;
  using Base::maxNumElems;
  using Base::maxNumConns;
  using Base::maxStencilSize;
  using Base::m_stencilWrapper;
  using Base::m_seri;
  using Base::m_sesri;
  using Base::m_sei;
  using Base::m_primarySpeciesAggregateConc;
  using Base::m_referencePorosity;

  using ThermalSinglePhaseFlowAccessors =
    StencilAccessors< fields::flow::temperature >;

  using ThermalReactiveSinglePhaseFluidAccessors =
    StencilMaterialAccessors< constitutive::ReactiveSingleFluid,
                              fields::singlefluid::enthalpy,
                              fields::singlefluid::dEnthalpy >;

  using ThermalConductivityAccessors =
    StencilMaterialAccessors< constitutive::SinglePhaseThermalConductivityBase,
                              fields::thermalconductivity::effectiveConductivity,
                              fields::thermalconductivity::dEffectiveConductivity_dT >;


  /**
   * @brief Constructor for the kernel interface
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] stencilWrapper reference to the stencil wrapper
   * @param[in] dofNumberAccessor accessor for the dofs numbers
   * @param[in] singlePhaseFlowAccessors accessor for wrappers registered by the solver
   * @param[in] reactiveSinglePhaseFlowAccessors accessor for *reactive* wrappers registered by the solver
   * @param[in] thermalSinglePhaseFlowAccessors accessor for *thermal* wrappers registered by the solver
   * @param[in] singlePhaseFluidAccessors accessor for wrappers registered by the single fluid model
   * @param[in] reactiveSinglePhaseFluidAccessors accessor for *reactive* wrappers registered by the single fluid model
   * @param[in] thermalReactiveSinglePhaseFluidAccessors accessor for *thermal reactive* wrappers registered by the single fluid model
   * @param[in] permeabilityAccessors accessor for wrappers registered by the permeability model
   * @param[in] diffusionAccessors accessor for wrappers registered by the diffusion model
   * @param[in] porosityAccessors accessor for wrappers registered by the porosity model
   * @param[in] thermalConductivityAccessors accessor for wrappers registered by the thermal conductivity model
   * @param[in] hasDiffusion the flag to turn on diffusion calculation
   * @param[in] dt time step size
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  FluxComputeKernel( globalIndex const rankOffset,
                     STENCILWRAPPER const & stencilWrapper,
                     DofNumberAccessor const & dofNumberAccessor,
                     SinglePhaseFlowAccessors const & singlePhaseFlowAccessors,
                     ReactiveSinglePhaseFlowAccessors const & reactiveSinglePhaseFlowAccessors,
                     ThermalSinglePhaseFlowAccessors const & thermalSinglePhaseFlowAccessors,
                     SinglePhaseFluidAccessors const & singlePhaseFluidAccessors,
                     ReactiveSinglePhaseFluidAccessors const & reactiveSinglePhaseFluidAccessors,
                     ThermalReactiveSinglePhaseFluidAccessors const & thermalReactiveSinglePhaseFluidAccessors,
                     PermeabilityAccessors const & permeabilityAccessors,
                     DiffusionAccessors const & diffusionAccessors,
                     PorosityAccessors const & porosityAccessors,
                     ThermalConductivityAccessors const & thermalConductivityAccessors,
                     integer const & hasDiffusion,
                     real64 const & dt,
                     CRSMatrixView< real64, globalIndex const > const & localMatrix,
                     arrayView1d< real64 > const & localRhs )
    : Base( rankOffset,
            stencilWrapper,
            dofNumberAccessor,
            singlePhaseFlowAccessors,
            reactiveSinglePhaseFlowAccessors,
            singlePhaseFluidAccessors,
            reactiveSinglePhaseFluidAccessors,
            permeabilityAccessors,
            diffusionAccessors,
            porosityAccessors,
            hasDiffusion,
            dt,
            localMatrix,
            localRhs ),
    m_temp( thermalSinglePhaseFlowAccessors.get( fields::flow::temperature {} ) ),
    m_enthalpy( thermalReactiveSinglePhaseFluidAccessors.get( fields::singlefluid::enthalpy {} ) ),
    m_dEnthalpy( thermalReactiveSinglePhaseFluidAccessors.get( fields::singlefluid::dEnthalpy {} ) ),
    // m_dPrimarySpeciesAggregateConcentration_dTemp( fluid.dPrimarySpeciesAggregateConcentration_dTemp() ),
    m_thermalConductivity( thermalConductivityAccessors.get( fields::thermalconductivity::effectiveConductivity {} ) ),
    m_dThermalCond_dT( thermalConductivityAccessors.get( fields::thermalconductivity::dEffectiveConductivity_dT {} ) )
  {}

  struct StackVariables : public Base::StackVariables
  {
public:

    GEOS_HOST_DEVICE
    StackVariables( localIndex const size, localIndex numElems )
      : Base::StackVariables( size, numElems ),
      energyFlux( 0.0 ),
      dEnergyFlux_dP( size ),
      dEnergyFlux_dT( size )
    {}

    using Base::StackVariables::stencilSize;
    using Base::StackVariables::numFluxElems;
    using Base::StackVariables::transmissibility;
    using Base::StackVariables::dTrans_dPres;
    using Base::StackVariables::dofColIndices;
    using Base::StackVariables::localFlux;
    using Base::StackVariables::localFluxJacobian;
    using Base::StackVariables::diffusionTransmissibility;
    using Base::StackVariables::dDiffusionTrans_dT;


    // Thermal transmissibility
    real64 thermalTransmissibility[maxNumConns][2]{};

    /// Derivatives of thermal transmissibility with respect to temperature
    real64 dThermalTrans_dT[maxNumConns][2]{};

    // Energy fluxes and derivatives

    /// Energy fluxes
    real64 energyFlux;
    /// Derivatives of energy fluxes wrt pressure
    stackArray1d< real64, maxStencilSize > dEnergyFlux_dP;
    /// Derivatives of energy fluxes wrt temperature
    stackArray1d< real64, maxStencilSize > dEnergyFlux_dT;

  };

  /**
   * @brief Compute the local flux contributions to the residual and Jacobian
   * @param[in] iconn the connection index
   * @param[inout] stack the stack variables
   */
  GEOS_HOST_DEVICE
  void computeFlux( localIndex const iconn,
                    StackVariables & stack ) const
  {
    using DerivOffset = constitutive::singlefluid::DerivativeOffsetC< 1 >;
    // ***********************************************
    // First, we call the base computeFlux to compute:
    //  1) massFlux and speciesFlux and their derivatives (including derivatives wrt temperature),
    //  2) enthalpy part of energyFlux and its derivatives (including derivatives wrt temperature)
    //
    // Computing dFlux_dT and the enthalpy flux requires quantities already computed in the base computeFlux,
    // such as potGrad, fluxVal, and the indices of the upwind cell
    // We use the lambda below (called **inside** the phase loop of the base computeFlux) to access these variables
    Base::computeFlux( iconn, stack, [&] ( localIndex const (&k)[2],
                                           localIndex const (&seri)[2],
                                           localIndex const (&sesri)[2],
                                           localIndex const (&sei)[2],
                                           localIndex const connectionIndex,
                                           real64 const alpha,
                                           real64 const mobility,
                                           real64 const & potGrad,
                                           real64 const & fluxVal,
                                           real64 const (&dFlux_dP)[2],
                                           real64 const fluidDens_up )
    {
      // Step 1: compute the derivatives of the (upwinded) massFlux wrt temperature
      // --------------------------------------------------------------------------
      // Step 1.1: compute the derivatives of the mean density at the interface wrt temperature
      real64 dDensMean_dT[numFluxSupportPoints]{0.0, 0.0};

      real64 const trans[numFluxSupportPoints] = { stack.transmissibility[connectionIndex][0], stack.transmissibility[connectionIndex][1] };

      for( integer ke = 0; ke < numFluxSupportPoints; ++ke )
      {
        dDensMean_dT[ke] = 0.5 * m_dDens[seri[ke]][sesri[ke]][sei[ke]][0][DerivOffset::dT];
      }

      // Step 1.2: compute the derivatives of the potential difference wrt temperature
      real64 dGravHead_dT[numFluxSupportPoints]{0.0, 0.0};

      // compute derivative of gravity potential difference wrt temperature
      for( integer ke = 0; ke < numFluxSupportPoints; ++ke )
      {
        localIndex const er  = seri[ke];
        localIndex const esr = sesri[ke];
        localIndex const ei  = sei[ke];

        real64 const gravD = trans[ke] * m_gravCoef[er][esr][ei];

        for( integer i = 0; i < numFluxSupportPoints; ++i )
        {
          dGravHead_dT[i] += dDensMean_dT[i] * gravD;
        }
      }

      real64 dFlux_dT[numFluxSupportPoints]{0.0, 0.0};

      for( integer ke = 0; ke < numFluxSupportPoints; ++ke )
      {
        dFlux_dT[ke] -= dGravHead_dT[ke];
      }

      for( integer ke = 0; ke < numFluxSupportPoints; ++ke )
      {
        dFlux_dT[ke] *= mobility;
      }

      // compute the derivatives of the mobility wrt temperature
      // *** upwinding ***
      real64 dMob_dT[numFluxSupportPoints]{};

      if( alpha <= 0.0 || alpha >= 1.0 )
      {
        localIndex const k_up = 1 - localIndex( fmax( fmin( alpha, 1.0 ), 0.0 ) );
        dMob_dT[k_up] = m_dMob[seri[k_up]][sesri[k_up]][sei[k_up]][DerivOffset::dT];
      }
      else
      {
        real64 const mobWeights[numFluxSupportPoints] = { alpha, 1.0 - alpha };
        for( integer ke = 0; ke < numFluxSupportPoints; ++ke )
        {
          dMob_dT[ke] = mobWeights[ke] * m_dMob[seri[ke]][sesri[ke]][sei[ke]][DerivOffset::dT];
        }
      }

      // add contribution from upstream cell mobility derivatives
      for( integer ke = 0; ke < numFluxSupportPoints; ++ke )
      {
        dFlux_dT[ke] += dMob_dT[ke] * potGrad;
      }

      // Step 1.3: populate local jacobian
      for( integer ke = 0; ke < numFluxSupportPoints; ++ke )
      {
        localIndex const localDofIndexTemp = k[ke] * numDof + numDof - numSpecies - 1;
        stack.localFluxJacobian[k[0]*numEqn][localDofIndexTemp] += m_dt * dFlux_dT[ke];
        stack.localFluxJacobian[k[1]*numEqn][localDofIndexTemp] -= m_dt * dFlux_dT[ke];
      }

      // Step 2: compute the derivatives of the speciesFlux wrt temperature
      // -------------------------------------------------------------------
      real64 dSpeciesFlux_dT[numFluxSupportPoints][numSpecies]{};

      {
        // Step 2.1: compute the derivatives of the upstream density wrt temperature
        // choose upstream cell
        localIndex const k_up = (potGrad >= 0) ? 0 : 1;

        localIndex const er_up  = seri[k_up];
        localIndex const esr_up = sesri[k_up];
        localIndex const ei_up  = sei[k_up];

        real64 const dDens_dTemp = m_dDens[er_up][esr_up][ei_up][0][DerivOffset::dT];

        // Step 2.2: compute speciesFlux derivative wrt temperature
        for( integer is = 0; is < numSpecies; ++is )
        {
          real64 const aggregateConc_i = m_primarySpeciesAggregateConc[er_up][esr_up][ei_up][is];

          // real64 const dAggregateConc_i_dTemp = m_dPrimarySpeciesAggregateConcentration_dTemp[er_up][esr_up][ei_up][is];
          // dSpeciesFlux_dT[k_up][is] += dAggregateConc_i_dTemp * fluxVal / fluidDens_up;
          dSpeciesFlux_dT[k_up][is] += -aggregateConc_i * fluxVal * dDens_dTemp / (fluidDens_up * fluidDens_up);

          for( integer ke = 0; ke < numFluxSupportPoints; ++ke )
          {
            dSpeciesFlux_dT[ke][is] += aggregateConc_i / fluidDens_up * dFlux_dT[ke];
          }
        }
      }

      // Step 2.3: populate local jacobian
      for( integer is = 0; is < numSpecies; ++is )
      {
        integer const eqIndex0 = k[0] * numEqn + numEqn - numSpecies + is;
        integer const eqIndex1 = k[1] * numEqn + numEqn - numSpecies + is;

        for( integer ke = 0; ke < numFluxSupportPoints; ++ke )
        {
          localIndex const localDofIndexTemp = k[ke] * numDof + numDof - numSpecies - 1;

          stack.localFluxJacobian[eqIndex0][localDofIndexTemp] += m_dt * dSpeciesFlux_dT[ke][is];
          stack.localFluxJacobian[eqIndex1][localDofIndexTemp] -= m_dt * dSpeciesFlux_dT[ke][is];
        }
      }

      // Step 3: compute the enthalpy flux
      // ----------------------------------
      real64 enthalpy = 0.0;
      real64 dEnthalpy_dP[numFluxSupportPoints]{0.0, 0.0};
      real64 dEnthalpy_dT[numFluxSupportPoints]{0.0, 0.0};
      // Todo: to add the enthalpy derivatives wrt speciesConc if needed
      // real64 dEnthalpy_dLogConc[numFluxSupportPoints][numSpecies]{};

      if( alpha <= 0.0 || alpha >= 1.0 )
      {
        localIndex const k_up = 1 - localIndex( fmax( fmin( alpha, 1.0 ), 0.0 ) );

        enthalpy = m_enthalpy[seri[k_up]][sesri[k_up]][sei[k_up]][0];
        dEnthalpy_dP[k_up] = m_dEnthalpy[seri[k_up]][sesri[k_up]][sei[k_up]][0][DerivOffset::dP];
        dEnthalpy_dT[k_up] = m_dEnthalpy[seri[k_up]][sesri[k_up]][sei[k_up]][0][DerivOffset::dT];
      }
      else
      {
        real64 const mobWeights[numFluxSupportPoints] = { alpha, 1.0 - alpha };
        for( integer ke = 0; ke < numFluxSupportPoints; ++ke )
        {
          enthalpy += mobWeights[ke] * m_enthalpy[seri[ke]][sesri[ke]][sei[ke]][0];
          dEnthalpy_dP[ke] = mobWeights[ke] * m_dEnthalpy[seri[ke]][sesri[ke]][sei[ke]][0][DerivOffset::dP];
          dEnthalpy_dT[ke] = mobWeights[ke] * m_dEnthalpy[seri[ke]][sesri[ke]][sei[ke]][0][DerivOffset::dT];
        }
      }

      stack.energyFlux += fluxVal * enthalpy;

      for( integer ke = 0; ke < numFluxSupportPoints; ++ke )
      {
        stack.dEnergyFlux_dP[ke] += dFlux_dP[ke] * enthalpy;
        stack.dEnergyFlux_dT[ke] += dFlux_dT[ke] * enthalpy;
      }

      for( integer ke = 0; ke < numFluxSupportPoints; ++ke )
      {
        stack.dEnergyFlux_dP[ke] += fluxVal * dEnthalpy_dP[ke];
        stack.dEnergyFlux_dT[ke] += fluxVal * dEnthalpy_dT[ke];
      }

    } );

    // *****************************************************
    // Computation of the conduction term in the energy flux
    // Note that the enthalpy term in the energy was computed above
    // Note that this term is computed using an explicit treatment of conductivity for now

    // Step 1: compute the thermal transmissibilities at this face
    // We follow how the thermal compositional multi-phase solver does to update the thermal transmissibility
    m_stencilWrapper.computeWeights( iconn,
                                     m_thermalConductivity,
                                     m_dThermalCond_dT,
                                     stack.thermalTransmissibility,
                                     stack.dThermalTrans_dT );

    localIndex k[numFluxSupportPoints];
    localIndex connectionIndex = 0;

    for( k[0] = 0; k[0] < stack.numFluxElems; ++k[0] )
    {
      for( k[1] = k[0] + 1; k[1] < stack.numFluxElems; ++k[1] )
      {
        real64 const thermalTrans[numFluxSupportPoints] = { stack.thermalTransmissibility[connectionIndex][0], stack.thermalTransmissibility[connectionIndex][1] };
        real64 const dThermalTrans_dT[numFluxSupportPoints] = { stack.dThermalTrans_dT[connectionIndex][0], stack.dThermalTrans_dT[connectionIndex][1] };

        localIndex const seri[numFluxSupportPoints]  = {m_seri( iconn, k[0] ), m_seri( iconn, k[1] )};
        localIndex const sesri[numFluxSupportPoints] = {m_sesri( iconn, k[0] ), m_sesri( iconn, k[1] )};
        localIndex const sei[numFluxSupportPoints]   = {m_sei( iconn, k[0] ), m_sei( iconn, k[1] )};

        // Step 2: compute temperature difference at the interface
        for( integer ke = 0; ke < numFluxSupportPoints; ++ke )
        {
          localIndex const er  = seri[ke];
          localIndex const esr = sesri[ke];
          localIndex const ei  = sei[ke];

          stack.energyFlux += thermalTrans[ke] * m_temp[er][esr][ei];
          stack.dEnergyFlux_dT[ke] += thermalTrans[ke] + dThermalTrans_dT[ke] * m_temp[er][esr][ei];
        }

        integer const eqIndex0 = k[0] * numEqn + numEqn - numSpecies - 1;
        integer const eqIndex1 = k[1] * numEqn + numEqn - numSpecies - 1;

        // add energyFlux and its derivatives to localFlux and localFluxJacobian
        stack.localFlux[eqIndex0] += m_dt * stack.energyFlux;
        stack.localFlux[eqIndex1] -= m_dt * stack.energyFlux;

        for( integer ke = 0; ke < numFluxSupportPoints; ++ke )
        {
          integer const localDofIndexPres = k[ke] * numDof;
          stack.localFluxJacobian[eqIndex0][localDofIndexPres] =  m_dt * stack.dEnergyFlux_dP[ke];
          stack.localFluxJacobian[eqIndex1][localDofIndexPres] = -m_dt * stack.dEnergyFlux_dP[ke];
          integer const localDofIndexTemp = localDofIndexPres + numDof - numSpecies - 1;
          stack.localFluxJacobian[eqIndex0][localDofIndexTemp] =  m_dt * stack.dEnergyFlux_dT[ke];
          stack.localFluxJacobian[eqIndex1][localDofIndexTemp] = -m_dt * stack.dEnergyFlux_dT[ke];
        }

        connectionIndex++;
      }
    }
  }

  /**
   * @brief Compute the local flux contributions to the residual and Jacobian
   * @param[in] iconn the connection index
   * @param[inout] stack the stack variables
   */
  GEOS_HOST_DEVICE
  void computeDiffusion( localIndex const iconn,
                         StackVariables & stack ) const
  {
    Base::computeDiffusion( iconn, stack, [&] ( integer const is,
                                                localIndex const (&k)[2],
                                                localIndex const (&seri)[2],
                                                localIndex const (&sesri)[2],
                                                localIndex const (&sei)[2],
                                                localIndex const connectionIndex,
                                                localIndex const k_up )
    {
      real64 dDiffusionFlux_dT[numFluxSupportPoints]{};
      real64 dSpeciesGrad_dT[numFluxSupportPoints]{};

      // Calculate diffusion derivative wrt temperature
      for( integer ke = 0; ke < numFluxSupportPoints; ++ke )
      {
        localIndex const er  = seri[ke];
        localIndex const esr = sesri[ke];
        localIndex const ei  = sei[ke];

        // dSpeciesGrad_dT[ke] += stack.diffusionTransmissibility[connectionIndex][ke]
        //                        * m_dPrimarySpeciesAggregateConcentration_dTemp[er][esr][ei][is];

        dSpeciesGrad_dT[ke] += stack.dDiffusionTrans_dT[connectionIndex][ke] * m_primarySpeciesAggregateConc[er][esr][ei][is];
      }

      for( integer ke = 0; ke < numFluxSupportPoints; ke++ )
      {
        localIndex const er_up  = seri[k_up];
        localIndex const esr_up = sesri[k_up];
        localIndex const ei_up  = sei[k_up];

        dDiffusionFlux_dT[ke] += m_referencePorosity[er_up][esr_up][ei_up] * dSpeciesGrad_dT[ke];
      }

      // populate local Jacobian
      integer const eqIndex0 = k[0] * numEqn + numEqn - numSpecies + is;
      integer const eqIndex1 = k[1] * numEqn + numEqn - numSpecies + is;

      for( integer ke = 0; ke < numFluxSupportPoints; ++ke )
      {
        localIndex const localDofIndexTemp = k[ke] * numDof + numDof - numSpecies - 1;
        stack.localFluxJacobian[eqIndex0][localDofIndexTemp] += m_dt * dDiffusionFlux_dT[ke];
        stack.localFluxJacobian[eqIndex1][localDofIndexTemp] -= m_dt * dDiffusionFlux_dT[ke];
      }
    } );
  }

  /**
   * @brief Performs the complete phase for the kernel.
   * @param[in] iconn the connection index
   * @param[inout] stack the stack variables
   */
  GEOS_HOST_DEVICE
  void complete( localIndex const iconn,
                 StackVariables & stack ) const
  {
    // Call Case::complete to assemble the mass balance equations
    // In the lambda, add contribution to residual and jacobian into the energy balance equation
    Base::complete( iconn, stack, [&] ( integer const i,
                                        localIndex const localRow )
    {
      // The no. of fluxes is equal to the no. of equations in m_localRhs and m_localMatrix
      // Different from the one in compositional multi-phase flow, which has a volume balance eqn.
      RAJA::atomicAdd( parallelDeviceAtomic{}, &AbstractBase::m_localRhs[localRow + numEqn - numSpecies - 1], stack.localFlux[i * numEqn + numEqn - numSpecies - 1] );

      AbstractBase::m_localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( localRow + numEqn - numSpecies - 1,
                                                                                        stack.dofColIndices.data(),
                                                                                        stack.localFluxJacobian[i * numEqn + numEqn - numSpecies - 1].dataIfContiguous(),
                                                                                        stack.stencilSize * numDof );

    } );
  }

protected:

  /// Views on temperature
  ElementViewConst< arrayView1d< real64 const > > const m_temp;

  /// Views on enthalpies
  ElementViewConst< arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID > > const m_enthalpy;
  ElementViewConst< arrayView3d< real64 const, constitutive::singlefluid::USD_FLUID_DER > > const m_dEnthalpy;

  // /// Views on the derivative of primary species aggregate concentration wrt temperature
  // ElementViewConst< arrayView2d< real64 const, compflow::USD_COMP > > const m_dPrimarySpeciesAggregateConc_dTemp;

  /// View on thermal conductivity
  ElementViewConst< arrayView3d< real64 const > > m_thermalConductivity;

  /// View on derivatives of thermal conductivity w.r.t. temperature
  ElementViewConst< arrayView3d< real64 const > > m_dThermalCond_dT;

};

/**
 * @class FluxComputeKernelFactory
 */
class FluxComputeKernelFactory
{
public:

  /**
   * @brief Create a new kernel and launch
   * @tparam POLICY the policy used in the RAJA kernel
   * @tparam STENCILWRAPPER the type of the stencil wrapper
   * @param[in] numSpecies the number of primary species
   * @param[in] hasDiffusion the flag of adding diffusion term
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofKey string to get the element degrees of freedom numbers
   * @param[in] solverName name of the solver (to name accessors)
   * @param[in] elemManager reference to the element region manager
   * @param[in] stencilWrapper reference to the stencil wrapper
   * @param[in] dt time step size
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  template< typename POLICY, typename STENCILWRAPPER >
  static void
  createAndLaunch( integer const numSpecies,
                   integer const hasDiffusion,
                   globalIndex const rankOffset,
                   string const & dofKey,
                   string const & solverName,
                   ElementRegionManager const & elemManager,
                   STENCILWRAPPER const & stencilWrapper,
                   real64 const & dt,
                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                   arrayView1d< real64 > const & localRhs )
  {
    singlePhaseReactiveBaseKernels::internal::kernelLaunchSelectorCompSwitch( numSpecies, [&]( auto NS )
    {
      integer constexpr NUM_SPECIES = NS();
      integer constexpr NUM_DOF = 2+NS();
      integer constexpr NUM_EQN = 2+NS();

      ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > dofNumberAccessor =
        elemManager.constructArrayViewAccessor< globalIndex, 1 >( dofKey );
      dofNumberAccessor.setName( solverName + "/accessors/" + dofKey );

      using KernelType = FluxComputeKernel< NUM_SPECIES, NUM_EQN, NUM_DOF, STENCILWRAPPER >;
      typename KernelType::SinglePhaseFlowAccessors flowAccessors( elemManager, solverName );
      typename KernelType::ReactiveSinglePhaseFlowAccessors reactiveFlowAccessors( elemManager, solverName );
      typename KernelType::ThermalSinglePhaseFlowAccessors thermalFlowAccessors( elemManager, solverName );
      typename KernelType::SinglePhaseFluidAccessors fluidAccessors( elemManager, solverName );
      typename KernelType::ReactiveSinglePhaseFluidAccessors reactiveFluidAccessors( elemManager, solverName );
      typename KernelType::ThermalReactiveSinglePhaseFluidAccessors thermalFluidAccessors( elemManager, solverName );
      typename KernelType::PermeabilityAccessors permAccessors( elemManager, solverName );
      typename KernelType::DiffusionAccessors diffusionAccessors( elemManager, solverName );
      typename KernelType::PorosityAccessors porosityAccessors( elemManager, solverName );
      typename KernelType::ThermalConductivityAccessors thermalConductivityAccessors( elemManager, solverName );

      KernelType kernel( rankOffset, stencilWrapper, dofNumberAccessor,
                         flowAccessors, reactiveFlowAccessors, thermalFlowAccessors, fluidAccessors, reactiveFluidAccessors, thermalFluidAccessors,
                         permAccessors, diffusionAccessors, porosityAccessors, thermalConductivityAccessors,
                         hasDiffusion, dt, localMatrix, localRhs );
      KernelType::template launch< POLICY >( stencilWrapper.size(), kernel );
    } );
  }
};

} // namespace thermalSinglePhaseReactiveFVMKernels

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEREACTIVE_THERMALFLUXCOMPUTEKERNEL_HPP
