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
 * @file FluxComputeKernel.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_REACTIVE_FLUXCOMPUTEKERNEL_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_REACTIVE_FLUXCOMPUTEKERNEL_HPP

#include "constitutive/diffusion/DiffusionFields.hpp"
#include "constitutive/diffusion/DiffusionBase.hpp"
#include "constitutive/solid/porosity/PorosityBase.hpp"
#include "constitutive/solid/porosity/PorosityFields.hpp"
#include "constitutive/fluid/singlefluid/reactive/ReactiveSingleFluid.hpp"
#include "constitutive/fluid/multifluid/reactive/ReactiveMultiFluidFields.hpp"
#include "physicsSolvers/fluidFlow/kernels/singlePhase/FluxComputeKernel.hpp"
#include "physicsSolvers/fluidFlow/kernels/singlePhase/reactive/KernelLaunchSelectors.hpp"


namespace geos
{

namespace singlePhaseReactiveFVMKernels
{

/**
 * @class FluxComputeKernel
 * @tparam NUM_SPECIES number of fluid primary species
 * @tparam NUM_DOF number of degrees of freedom
 * @tparam STENCILWRAPPER the type of the stencil wrapper
 * @brief Define the interface for the assembly kernel in charge of flux terms
 */
template< integer NUM_SPECIES, integer NUM_DOF, typename STENCILWRAPPER >
class FluxComputeKernel : public singlePhaseFVMKernels::FluxComputeKernel< NUM_SPECIES+1, NUM_DOF, STENCILWRAPPER >
{
public:

  /// Compile time value for the number of primary species
  static constexpr integer numSpecies = NUM_SPECIES;

  /// Number of flux support points (hard-coded for TFPA)
  static constexpr integer numFluxSupportPoints = 2;

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
  using AbstractBase::m_dens;
  using AbstractBase::m_dDens;

  using Base = singlePhaseFVMKernels::FluxComputeKernel< NUM_SPECIES+1, NUM_DOF, STENCILWRAPPER >;
  using Base::numDof;
  using Base::numEqn;
  using Base::maxNumElems;
  using Base::maxNumConns;
  using Base::maxStencilSize;
  using Base::m_stencilWrapper;
  using Base::m_seri;
  using Base::m_sesri;
  using Base::m_sei;

  using ReactiveSinglePhaseFlowAccessors =
    StencilAccessors< fields::flow::logPrimarySpeciesConcentration,
                      fields::flow::dMobility_dLogPrimaryConc >;

  using ReactiveSinglePhaseFluidAccessors =
    StencilMaterialAccessors< constitutive::ReactiveSingleFluid,
                              fields::reactivefluid::primarySpeciesAggregateConcentration,
                              fields::reactivefluid::dPrimarySpeciesAggregateConcentration_dLogPrimaryConc >;

  using DiffusionAccessors =
    StencilMaterialAccessors< constitutive::DiffusionBase,
                              fields::diffusion::diffusivity,
                              fields::diffusion::dDiffusivity_dTemperature >;

  using PorosityAccessors =
    StencilMaterialAccessors< constitutive::PorosityBase,
                              fields::porosity::referencePorosity >;

  /**
   * @brief Constructor for the kernel interface
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] stencilWrapper reference to the stencil wrapper
   * @param[in] dofNumberAccessor
   * @param[in] singlePhaseFlowAccessors
   * @param[in] singlePhaseFluidAccessors
   * @param[in] permeabilityAccessors
   * @param[in] dt time step size
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  FluxComputeKernel( globalIndex const rankOffset,
                     STENCILWRAPPER const & stencilWrapper,
                     DofNumberAccessor const & dofNumberAccessor,
                     SinglePhaseFlowAccessors const & singlePhaseFlowAccessors,
                     ReactiveSinglePhaseFlowAccessors const & reactiveSinglePhaseFlowAccessors,
                     SinglePhaseFluidAccessors const & singlePhaseFluidAccessors,
                     ReactiveSinglePhaseFluidAccessors const & reactiveSinglePhaseFluidAccessors,
                     PermeabilityAccessors const & permeabilityAccessors,
                     DiffusionAccessors const & diffusionAccessors,
                     PorosityAccessors const & porosityAccessors,
                     integer const & hasDiffusion,
                     real64 const & dt,
                     CRSMatrixView< real64, globalIndex const > const & localMatrix,
                     arrayView1d< real64 > const & localRhs )
    : Base( rankOffset,
            stencilWrapper,
            dofNumberAccessor,
            singlePhaseFlowAccessors,
            singlePhaseFluidAccessors,
            permeabilityAccessors,
            dt,
            localMatrix,
            localRhs ),
    m_logPrimarySpeciesConc( reactiveSinglePhaseFlowAccessors.get( fields::flow::logPrimarySpeciesConcentration {} ) ),
    m_dMob_dLogPrimaryConc( reactiveSinglePhaseFlowAccessors.get( fields::flow::dMobility_dLogPrimaryConc {} ) ),
    m_primarySpeciesAggregateConc( reactiveSinglePhaseFluidAccessors.get( fields::reactivefluid::primarySpeciesAggregateConcentration {} ) ),
    m_dPrimarySpeciesAggregateConc_dLogPrimaryConc( reactiveSinglePhaseFluidAccessors.get( fields::reactivefluid::dPrimarySpeciesAggregateConcentration_dLogPrimaryConc {} ) ),
    m_diffusivity( diffusionAccessors.get( fields::diffusion::diffusivity {} ) ),
    m_dDiffusivity_dTemp( diffusionAccessors.get( fields::diffusion::dDiffusivity_dTemperature {} ) ),
    m_referencePorosity( porosityAccessors.get( fields::porosity::referencePorosity {} ) ),
    m_hasDiffusion( hasDiffusion )
  {}

  /**
   * @struct StackVariables
   * @brief Kernel variables (dof numbers, jacobian and residual) located on the stack
   */
  struct StackVariables : public Base::StackVariables
  {
public:

    /**
     * @brief Constructor for the stack variables
     * @param[in] size size of the stencil for this connection
     * @param[in] numElems number of elements for this connection
     */
    GEOS_HOST_DEVICE
    StackVariables( localIndex const size, localIndex numElems )
      : Base::StackVariables( size, numElems )
    {}

    using Base::StackVariables::stencilSize;
    using Base::StackVariables::numFluxElems;
    using Base::StackVariables::transmissibility;
    using Base::StackVariables::dTrans_dPres;
    using Base::StackVariables::dofColIndices;
    using Base::StackVariables::localFlux;
    using Base::StackVariables::localFluxJacobian;

    /// Diffusion transmissibility
    real64 diffusionTransmissibility[maxNumConns][numFluxSupportPoints]{};
    /// Derivatives of diffusion transmissibility with respect to temperature
    real64 dDiffusionTrans_dT[maxNumConns][numFluxSupportPoints]{};

  };

  /**
   * @brief Compute the local flux contributions to the residual and Jacobian
   * @tparam FUNC the type of the function that can be used to customize the computation of the flux
   * @param[in] iconn the connection index
   * @param[inout] stack the stack variables
   * @param[in] NoOpFunc the function used to customize the computation of the flux
   */
  template< typename FUNC = NoOpFunc >
  GEOS_HOST_DEVICE
  void computeFlux( localIndex const iconn,
                    StackVariables & stack,
                    FUNC && kernelOp = NoOpFunc{} ) const
  {
    using DerivOffset = constitutive::singlefluid::DerivativeOffsetC< 1 >;
    // ***********************************************
    // First, we call the base computeFlux to compute:
    //  1) massFlux and its derivatives,
    //  2) speciesFlux and its derivatives
    Base::computeFlux( iconn, stack, [&] ( localIndex const (&k)[2],
                                           localIndex const (&seri)[2],
                                           localIndex const (&sesri)[2],
                                           localIndex const (&sei)[2],
                                           localIndex const connectionIndex,
                                           real64 const alpha,
                                           real64 const mobility,
                                           real64 const & potGrad,
                                           real64 const & fluxVal,
                                           real64 const (&dFlux_dP)[2] )
    {
      GEOS_UNUSED_VAR( connectionIndex, alpha, mobility );
      // Step 1: compute the derivatives of the fluid density, potential difference,
      // and the massFlux wrt log of primary species concentration (to complete)

      // Step 2: compute the speciesFlux
      real64 speciesFlux[numSpecies]{};
      real64 dSpeciesFlux_dP[numFluxSupportPoints][numSpecies]{};
      real64 dSpeciesFlux_dLogConc[numFluxSupportPoints][numSpecies][numSpecies]{};
      // real64 dSpeciesFlux_dTrans[numSpecies]{};

      // choose upstream cell
      localIndex const k_up = (potGrad >= 0) ? 0 : 1;

      localIndex const er_up  = seri[k_up];
      localIndex const esr_up = sesri[k_up];
      localIndex const ei_up  = sei[k_up];

      real64 const fluidDens_up = m_dens[er_up][esr_up][ei_up][0];
      real64 const dDens_dPres = m_dDens[er_up][esr_up][ei_up][0][DerivOffset::dP];

      // compute species fluxes and derivatives using upstream cell composition
      for( integer is = 0; is < numSpecies; ++is )
      {
        real64 const aggregateConc_i = m_primarySpeciesAggregateConc[er_up][esr_up][ei_up][is];
        speciesFlux[is] = aggregateConc_i / fluidDens_up * fluxVal;

        for( integer ke = 0; ke < numFluxSupportPoints; ++ke )
        {
          dSpeciesFlux_dP[ke][is] += aggregateConc_i / fluidDens_up * dFlux_dP[ke];
        }

        dSpeciesFlux_dP[k_up][is] += -aggregateConc_i * fluxVal * dDens_dPres / (fluidDens_up * fluidDens_up);

        for( integer js = 0; js < numSpecies; ++js )
        {
          real64 const dAggregateConc_i_dLogConc_j = m_dPrimarySpeciesAggregateConc_dLogPrimaryConc[er_up][esr_up][ei_up][is][js];
          dSpeciesFlux_dLogConc[k_up][is][js] += dAggregateConc_i_dLogConc_j / fluidDens_up * fluxVal;
        }
      }

      /// populate local flux vector and derivatives
      for( integer is = 0; is < numSpecies; ++is )
      {
        integer const eqIndex0 = k[0] * numEqn + is + 1;
        integer const eqIndex1 = k[1] * numEqn + is + 1;

        stack.localFlux[eqIndex0] +=  m_dt * speciesFlux[is];
        stack.localFlux[eqIndex1] -=  m_dt * speciesFlux[is];

        for( integer ke = 0; ke < numFluxSupportPoints; ++ke )
        {
          localIndex const localDofIndexPres = k[ke] * numDof;
          stack.localFluxJacobian[eqIndex0][localDofIndexPres] += m_dt * dSpeciesFlux_dP[ke][is];
          stack.localFluxJacobian[eqIndex1][localDofIndexPres] -= m_dt * dSpeciesFlux_dP[ke][is];

          for( integer js = 0; js < numSpecies; ++js )
          {
            localIndex const localDofIndexSpecies = localDofIndexPres + js + 1;
            stack.localFluxJacobian[eqIndex0][localDofIndexSpecies] += m_dt * dSpeciesFlux_dLogConc[ke][is][js];
            stack.localFluxJacobian[eqIndex1][localDofIndexSpecies] -= m_dt * dSpeciesFlux_dLogConc[ke][is][js];
          }
        }
      }

      // Customize the kernel with this lambda
      kernelOp( k, seri, sesri, sei, connectionIndex, alpha, mobility, potGrad, fluxVal, dFlux_dP );
    } );

    // *****************************************************
    // Computation of the diffusion term in the species flux

    if( m_hasDiffusion )
    {
      // Step 1: compute the diffusion transmissibilities at this face
      m_stencilWrapper.computeWeights( iconn,
                                       m_diffusivity,
                                       m_dDiffusivity_dTemp,
                                       stack.diffusionTransmissibility,
                                       stack.dDiffusionTrans_dT );

      localIndex k[numFluxSupportPoints]{};
      localIndex connectionIndex = 0;

      for( k[0] = 0; k[0] < stack.numFluxElems; ++k[0] )
      {
        for( k[1] = k[0] + 1; k[1] < stack.numFluxElems; ++k[1] )
        {
          localIndex const seri[numFluxSupportPoints]  = {m_seri( iconn, k[0] ), m_seri( iconn, k[1] )};
          localIndex const sesri[numFluxSupportPoints] = {m_sesri( iconn, k[0] ), m_sesri( iconn, k[1] )};
          localIndex const sei[numFluxSupportPoints]   = {m_sei( iconn, k[0] ), m_sei( iconn, k[1] )};

          // clear working arrays
          real64 diffusionFlux[numSpecies]{};
          // real64 dDiffusionFlux_dP[numFluxSupportPoints][numSpecies]{}; // Turn on if diffusionFlux is pressure-dependent
          real64 dDiffusionFlux_dLogConc[numFluxSupportPoints][numSpecies][numSpecies]{};

          real64 const diffusionTrans[numFluxSupportPoints] = { stack.diffusionTransmissibility[connectionIndex][0],
                                                                stack.diffusionTransmissibility[connectionIndex][1] };

          //***** calculation of flux *****
          // loop over primary species
          for( integer is = 0; is < numSpecies; ++is )
          {
            real64 speciesGrad_i = 0.0;
            // real64 dSpeciesGrad_i_dP[numFluxSupportPoints]{}; // Turn on if speciesGrad is pressure-dependent
            real64 dSpeciesGrad_i_dLogConc[numFluxSupportPoints][numSpecies]{};

            // Step 2: compute species gradient
            for( integer ke = 0; ke < numFluxSupportPoints; ++ke )
            {
              localIndex const er  = seri[ke];
              localIndex const esr = sesri[ke];
              localIndex const ei  = sei[ke];

              real64 const aggregateConc_i = m_primarySpeciesAggregateConc[er][esr][ei][is];

              speciesGrad_i += diffusionTrans[ke] * aggregateConc_i;

              for( integer js = 0; js < numSpecies; ++js )
              {
                real64 const dAggregateConc_i_dLogConc_j = m_dPrimarySpeciesAggregateConc_dLogPrimaryConc[er][esr][ei][is][js];

                dSpeciesGrad_i_dLogConc[ke][js] += diffusionTrans[ke] * dAggregateConc_i_dLogConc_j;
              }
            }

            // choose upstream cell for species upwinding
            localIndex const k_up = (speciesGrad_i >= 0) ? 0 : 1;

            localIndex const er_up  = seri[k_up];
            localIndex const esr_up = sesri[k_up];
            localIndex const ei_up  = sei[k_up];

            // computation of the upwinded species flux
            diffusionFlux[is] += m_referencePorosity[er_up][esr_up][ei_up] * speciesGrad_i;

            // add contributions of the derivatives of component fractions wrt pressure/component fractions
            for( integer ke = 0; ke < numFluxSupportPoints; ke++ )
            {
              for( integer js = 0; js < numSpecies; ++js )
              {
                dDiffusionFlux_dLogConc[ke][is][js] += m_referencePorosity[er_up][esr_up][ei_up] * dSpeciesGrad_i_dLogConc[ke][js];
              }
            }
          } // loop over primary species

          // Add the local diffusion flux contribution to the residual and Jacobian
          // loop over primary species
          for( integer is = 0; is < numSpecies; ++is )
          {
            integer const eqIndex0 = k[0] * numEqn + is + 1;
            integer const eqIndex1 = k[1] * numEqn + is + 1;

            stack.localFlux[eqIndex0] += m_dt * diffusionFlux[is];
            stack.localFlux[eqIndex1] -= m_dt * diffusionFlux[is];

            for( integer ke = 0; ke < numFluxSupportPoints; ++ke )
            {
              localIndex const localDofIndexPres = k[ke] * numDof;
              // stack.localFluxJacobian[eqIndex0][localDofIndexPres] += m_dt * dDiffusionFlux_dP[ke][is];
              // stack.localFluxJacobian[eqIndex1][localDofIndexPres] -= m_dt * dDiffusionFlux_dP[ke][is];

              for( integer js = 0; js < numSpecies; ++js )
              {
                localIndex const localDofIndexComp = localDofIndexPres + js + 1;
                stack.localFluxJacobian[eqIndex0][localDofIndexComp] += m_dt * dDiffusionFlux_dLogConc[ke][is][js];
                stack.localFluxJacobian[eqIndex1][localDofIndexComp] -= m_dt * dDiffusionFlux_dLogConc[ke][is][js];
              }
            }
          }

          connectionIndex++;
        }
      }
    }
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
    // Call Base::complete to assemble the total mass balance equation
    // In the lambda, add contribution to residual and jacobian into the species amount balance equation
    Base::complete( iconn, stack, [&] ( integer const i,
                                        localIndex const localRow )
    {
      // The no. of fluxes is equal to the no. of equations in m_localRhs and m_localMatrix
      for( integer is = 0; is < numSpecies; ++is )
      {
        RAJA::atomicAdd( parallelDeviceAtomic{}, &AbstractBase::m_localRhs[localRow + is + 1],
                         stack.localFlux[i * numEqn + is + 1] );
        AbstractBase::m_localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >
          ( localRow + is + 1,
          stack.dofColIndices.data(),
          stack.localFluxJacobian[i * numEqn + is + 1].dataIfContiguous(),
          stack.stencilSize * numDof );
      }
    } );
  }

protected:

  /// Views on log of primary species concentration
  ElementViewConst< arrayView2d< real64 const, compflow::USD_COMP > > const m_logPrimarySpeciesConc;

  /// Views on derivatives of fluid mobilities
  ElementViewConst< arrayView2d< real64 const, compflow::USD_FLUID_DC > > const m_dMob_dLogPrimaryConc;

  /// Views on primary species aggregate concentration
  ElementViewConst< arrayView2d< real64 const, compflow::USD_COMP > > const m_primarySpeciesAggregateConc;

  /// Views on primary species aggregate concentration
  ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const m_dPrimarySpeciesAggregateConc_dLogPrimaryConc;

  /// Views on diffusivity
  ElementViewConst< arrayView3d< real64 const > > const m_diffusivity;
  ElementViewConst< arrayView3d< real64 const > > const m_dDiffusivity_dTemp;

  /// View on the reference porosity
  ElementViewConst< arrayView1d< real64 const > > const m_referencePorosity;

  /// Flag of adding the diffusion term
  integer const m_hasDiffusion;
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
      integer constexpr NUM_DOF = 1+NS();

      ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > dofNumberAccessor =
        elemManager.constructArrayViewAccessor< globalIndex, 1 >( dofKey );
      dofNumberAccessor.setName( solverName + "/accessors/" + dofKey );

      using KernelType = FluxComputeKernel< NUM_SPECIES, NUM_DOF, STENCILWRAPPER >;
      typename KernelType::SinglePhaseFlowAccessors flowAccessors( elemManager, solverName );
      typename KernelType::ReactiveSinglePhaseFlowAccessors reactiveFlowAccessors( elemManager, solverName );
      typename KernelType::SinglePhaseFluidAccessors fluidAccessors( elemManager, solverName );
      typename KernelType::ReactiveSinglePhaseFluidAccessors reactiveFluidAccessors( elemManager, solverName );
      typename KernelType::PermeabilityAccessors permAccessors( elemManager, solverName );
      typename KernelType::DiffusionAccessors diffusionAccessors( elemManager, solverName );
      typename KernelType::PorosityAccessors porosityAccessors( elemManager, solverName );

      KernelType kernel( rankOffset, stencilWrapper, dofNumberAccessor,
                         flowAccessors, reactiveFlowAccessors, fluidAccessors, reactiveFluidAccessors,
                         permAccessors, diffusionAccessors, porosityAccessors, hasDiffusion,
                         dt, localMatrix, localRhs );
      KernelType::template launch< POLICY >( stencilWrapper.size(), kernel );
    } );
  }
};

} // namespace singlePhaseReactiveFVMKernels

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_REACTIVE_FLUXCOMPUTEKERNEL_HPP
