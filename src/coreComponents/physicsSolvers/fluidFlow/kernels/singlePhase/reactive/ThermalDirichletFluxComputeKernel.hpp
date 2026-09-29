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
 * @file ThermalDirichletFluxComputeKernel.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEREACTIVE_THERMALDIRICHLETFLUXCOMPUTEKERNEL_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEREACTIVE_THERMALDIRICHLETFLUXCOMPUTEKERNEL_HPP

#include "physicsSolvers/fluidFlow/kernels/singlePhase/reactive/DirichletFluxComputeKernel.hpp"

#include "constitutive/thermalConductivity/SinglePhaseThermalConductivityBase.hpp"
#include "constitutive/thermalConductivity/ThermalConductivityFields.hpp"

namespace geos
{

namespace thermalSinglePhaseReactiveFVMKernels
{

/******************************** DirichletFluxComputeKernel ********************************/

/**
 * @class DirichletFluxComputeKernel
 * @tparam NUM_SPECIES number of fluid primary species
 * @tparam NUM_EQN number of equations
 * @tparam NUM_DOF number of degrees of freedom
 * @tparam FLUIDWRAPPER the type of the fluid wrapper
 * @tparam BASE_FLUID_TYPE the type of the base model for the reactive fluid model
 * @brief Define the interface for the assembly kernel in charge of Dirichlet face flux terms
 */
template< integer NUM_SPECIES, integer NUM_EQN, integer NUM_DOF, typename FLUIDWRAPPER, typename BASE_FLUID_TYPE >
class DirichletFluxComputeKernel : public singlePhaseReactiveFVMKernels::DirichletFluxComputeKernel< NUM_SPECIES, NUM_EQN, NUM_DOF, FLUIDWRAPPER, BASE_FLUID_TYPE >
{
public:

/**
 * @brief The type for element-based data. Consists entirely of ArrayView's.
 *
 * Can be converted from ElementRegionManager::ElementViewConstAccessor
 * by calling .toView() or .toViewConst() on an accessor instance
 */

  using DerivOffset = constitutive::singlefluid::DerivativeOffsetC< 1 >;
  template< typename VIEWTYPE >
  using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;

  using AbstractBase = singlePhaseFVMKernels::FluxComputeKernelBase;
  using DofNumberAccessor = AbstractBase::DofNumberAccessor;
  using PermeabilityAccessors = AbstractBase::PermeabilityAccessors;
  using SinglePhaseFlowAccessors = AbstractBase::SinglePhaseFlowAccessors;
  using SinglePhaseFluidAccessors = AbstractBase::SinglePhaseFluidAccessors;

  using AbstractBase::m_dt;
  using AbstractBase::m_rankOffset;
  using AbstractBase::m_dofNumber;
  using AbstractBase::m_ghostRank;
  using AbstractBase::m_gravCoef;
  using AbstractBase::m_mob;
  using AbstractBase::m_pres;
  using AbstractBase::m_permeability;
  using AbstractBase::m_dPerm_dPres;
  using AbstractBase::m_dDens;
  using AbstractBase::m_dMob;

  using Base = singlePhaseReactiveFVMKernels::DirichletFluxComputeKernel< NUM_SPECIES, NUM_EQN, NUM_DOF, FLUIDWRAPPER, BASE_FLUID_TYPE >;
  using Base::numSpecies;
  using Base::numDof;
  using Base::numEqn;
  using Base::m_stencilWrapper;
  using Base::m_seri;
  using Base::m_sesri;
  using Base::m_sei;
  using Base::m_facePres;
  using Base::m_faceGravCoef;

  using ReactiveSinglePhaseFlowAccessors = typename Base::ReactiveSinglePhaseFlowAccessors;
  using ReactiveSinglePhaseFluidAccessors = typename Base::ReactiveSinglePhaseFluidAccessors;

  using ThermalSinglePhaseFlowAccessors =
    StencilAccessors< fields::flow::temperature >;

  using ThermalReactiveSinglePhaseFluidAccessors =
    StencilMaterialAccessors< constitutive::reactivefluid::ReactiveSinglePhaseFluid< BASE_FLUID_TYPE >,
                              fields::singlefluid::enthalpy,
                              fields::singlefluid::dEnthalpy >;

  using ThermalConductivityAccessors =
    StencilMaterialAccessors< constitutive::SinglePhaseThermalConductivityBase,
                              fields::thermalconductivity::effectiveConductivity,
                              fields::thermalconductivity::dEffectiveConductivity_dT >;

  /**
   * @brief Constructor for the kernel interface
   * @param[in] rankOffset the offset of the MPI rank
   * @param[in] faceManager the face manager
   * @param[in] stencilWrapper reference to the stencil wrapper
   * @param[in] fluidWrapper reference to the fluid wrapper
   * @param[in] dofNumberAccessor the degree of freedom number accessor
   * @param[in] singlePhaseFlowAccessors the single phase flow accessor
   * @param[in] reactiveSinglePhaseFlowAccessors accessor for *reactive* wrappers registered by the solver
   * @param[in] thermalSinglePhaseFlowAccessors accessor for *thermal* wrappers registered by the solver
   * @param[in] singlePhaseFluidAccessors the single phase fluid accessor
   * @param[in] reactiveSinglePhaseFluidAccessors accessor for *reactive* wrappers registered by the single fluid model
   * @param[in] thermalReactiveSinglePhaseFluidAccessors accessor for *thermal reactive* wrappers registered by the single fluid model
   * @param[in] permeabilityAccessors the permeability accessor
   * @param[in] thermalConductivityAccessors the thermal conductivity accessor
   * @param[in] mobilePrimarySpeciesFlags the array of flags to indicate mobile primary species
   * @param[in] dt the time step size
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  DirichletFluxComputeKernel( globalIndex const rankOffset,
                              FaceManager const & faceManager,
                              BoundaryStencilWrapper const & stencilWrapper,
                              FLUIDWRAPPER const & fluidWrapper,
                              DofNumberAccessor const & dofNumberAccessor,
                              SinglePhaseFlowAccessors const & singlePhaseFlowAccessors,
                              ReactiveSinglePhaseFlowAccessors const & reactiveSinglePhaseFlowAccessors,
                              ThermalSinglePhaseFlowAccessors const & thermalSinglePhaseFlowAccessors,
                              SinglePhaseFluidAccessors const & singlePhaseFluidAccessors,
                              ReactiveSinglePhaseFluidAccessors const & reactiveSinglePhaseFluidAccessors,
                              ThermalReactiveSinglePhaseFluidAccessors const & thermalReactiveSinglePhaseFluidAccessors,
                              PermeabilityAccessors const & permeabilityAccessors,
                              ThermalConductivityAccessors const & thermalConductivityAccessors,
                              arrayView1d< integer const > const & mobilePrimarySpeciesFlags,
                              real64 const & dt,
                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                              arrayView1d< real64 > const & localRhs )

    : Base( rankOffset,
            faceManager,
            stencilWrapper,
            fluidWrapper,
            dofNumberAccessor,
            singlePhaseFlowAccessors,
            reactiveSinglePhaseFlowAccessors,
            singlePhaseFluidAccessors,
            reactiveSinglePhaseFluidAccessors,
            permeabilityAccessors,
            mobilePrimarySpeciesFlags,
            dt,
            localMatrix,
            localRhs ),
    m_temp( thermalSinglePhaseFlowAccessors.get( fields::flow::temperature {} ) ),
    m_faceTemp( faceManager.getField< fields::flow::faceTemperature >() ),
    m_enthalpy( thermalReactiveSinglePhaseFluidAccessors.get( fields::singlefluid::enthalpy {} ) ),
    m_dEnthalpy( thermalReactiveSinglePhaseFluidAccessors.get( fields::singlefluid::dEnthalpy {} ) ),
    m_thermalConductivity( thermalConductivityAccessors.get( fields::thermalconductivity::effectiveConductivity {} ) ),
    m_dThermalCond_dT( thermalConductivityAccessors.get( fields::thermalconductivity::dEffectiveConductivity_dT {} ) )
  {}


  /**
   * @struct StackVariables
   * @brief Kernel variables (dof numbers, jacobian and residual) located on the stack
   */
  struct StackVariables : Base::StackVariables
  {
public:

    /**
     * @brief Constructor for the stack variables
     * @param[in] size size of the stencil for this connection
     * @param[in] numElems number of elements for this connection
     */
    GEOS_HOST_DEVICE
    StackVariables( localIndex const size,
                    localIndex numElems ):
      Base::StackVariables( size,
                            numElems )
    {}

    using Base::StackVariables::localFlux;
    using Base::StackVariables::localFluxJacobian;
    using Base::StackVariables::dofColIndices;
    using Base::StackVariables::transmissibility;

    /// Energy fluxes and derivatives wrt pressure and temperature
    real64 energyFlux = 0.0;
    real64 dEnergyFlux_dP = 0.0;
    real64 dEnergyFlux_dT = 0.0;
  };

  /**
   * @brief Compute the local Dirichlet face flux contributions to the residual and Jacobian
   * @tparam FUNC the type of the function that can be used to customize the computation of the phase fluxes
   * @param[in] iconn the connection index
   * @param[inout] stack the stack variables
   * @param[in] compFluxKernelOp the function used to customize the computation of the component fluxes
   */
  GEOS_HOST_DEVICE
  void computeFlux( localIndex const iconn,
                    StackVariables & stack ) const
  {
    Base::computeFlux( iconn, stack, [&] ( localIndex const er,
                                           localIndex const esr,
                                           localIndex const ei,
                                           localIndex const kf,
                                           real64 const & f,
                                           real64 const & dF_dP,
                                           real64 const & mobility_up,
                                           real64 const & dMobility_dP_up )
    {
      // Compute the derivatives of the density wrt temperature

      real64 const dDens_dT = 0.5 * m_dDens[er][esr][ei][0][DerivOffset::dT];
      // Compute the derivatives of the phase potential difference wrt temperature

      real64 const dF_dT = -stack.transmissibility * dDens_dT * ( m_gravCoef[er][esr][ei] - m_faceGravCoef[kf] );

      // Compute the (upwinded) energy flux

      real64 const flux = mobility_up * f;
      real64 const enthalpy = m_enthalpy[er][esr][ei][0];
      stack.energyFlux += flux * enthalpy;

      // Compute the derivatives of the (upwinded) energy flux wrt pressure and temperature

      if( f >= 0 ) // the element is upstream
      {
        real64 const dFlux_dP = mobility_up * dF_dP + dMobility_dP_up * f;
        real64 const dFlux_dT = mobility_up * dF_dT + m_dMob[er][esr][ei][DerivOffset::dT] * f;

        stack.dEnergyFlux_dP += dFlux_dP * enthalpy + flux * m_dEnthalpy[er][esr][ei][0][DerivOffset::dP];
        stack.dEnergyFlux_dT += dFlux_dT * enthalpy + flux * m_dEnthalpy[er][esr][ei][0][DerivOffset::dT];
      }
      else
      {
        real64 const dFlux_dP = mobility_up * dF_dP;
        real64 const dFlux_dT = mobility_up * dF_dT;

        stack.dEnergyFlux_dP += dFlux_dP * enthalpy;
        stack.dEnergyFlux_dT += dFlux_dT * enthalpy;
      }

      // Contribution of energy conduction through the solid phase
      real64 thermalTrans = 0.0;
      real64 dThermalTrans_dThermalCond[3]{};
      m_stencilWrapper.computeWeights( iconn,
                                       m_thermalConductivity,
                                       thermalTrans,
                                       dThermalTrans_dThermalCond );

      real64 const dThermalTrans_dT = LvArray::tensorOps::AiBi< 3 >( dThermalTrans_dThermalCond, m_dThermalCond_dT[er][esr][ei][0] );

      real64 const deltaT = m_temp[er][esr][ei] - m_faceTemp[kf];
      stack.energyFlux += thermalTrans * deltaT;
      stack.dEnergyFlux_dT += thermalTrans + dThermalTrans_dT * deltaT;

      // Add energyFlux and its derivatives to localFlux and localFluxJacobian
      integer const localRowIndexEnergy = numEqn - numSpecies - 1;
      stack.localFlux[localRowIndexEnergy] =  m_dt * stack.energyFlux;

      stack.localFluxJacobian[localRowIndexEnergy][0] =  m_dt * stack.dEnergyFlux_dP;
      stack.localFluxJacobian[localRowIndexEnergy][numDof-numSpecies-1] =  m_dt * stack.dEnergyFlux_dT;
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
    Base::complete( iconn, stack, [&] ( localIndex const localRow )
    {
      RAJA::atomicAdd( parallelDeviceAtomic{}, &AbstractBase::m_localRhs[localRow + numEqn - numSpecies - 1],
                       stack.localFlux[numEqn - numSpecies - 1] );

      AbstractBase::m_localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >
        ( localRow + numEqn - numSpecies - 1,
        stack.dofColIndices,
        stack.localFluxJacobian[numEqn - numSpecies - 1],
        numDof );
    } );
  }

protected:

  /// Views on temperature
  ElementViewConst< arrayView1d< real64 const > > const m_temp;

  /// Views on face temperature
  arrayView1d< real64 const > const m_faceTemp;

  /// Views on enthalpies
  ElementViewConst< arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID > > const m_enthalpy;
  ElementViewConst< arrayView3d< real64 const, constitutive::singlefluid::USD_FLUID_DER > > const m_dEnthalpy;

  /// View on thermal conductivity
  ElementViewConst< arrayView3d< real64 const > > m_thermalConductivity;

  /// View on derivatives of thermal conductivity w.r.t. temperature
  ElementViewConst< arrayView3d< real64 const > > m_dThermalCond_dT;

};

/**
 * @class DirichletFluxComputeKernelFactory
 */
class DirichletFluxComputeKernelFactory
{
public:

  /**
   * @brief Create a new kernel and launch
   * @tparam POLICY the policy used in the RAJA kernel
   * @param[in] numSpecies the number of primary species
   * @param[in] mobilePrimarySpeciesFlags the array of flags to indicate mobile primary species
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofKey string to get the element degrees of freedom numbers
   * @param[in] solverName name of the solver (to name accessors)
   * @param[in] faceManager reference to the face manager
   * @param[in] elemManager reference to the element region manager
   * @param[in] stencilWrapper reference to the boundary stencil wrapper
   * @param[in] reactiveFluid the single phase reactive fluid constitutive model
   * @param[in] dt time step size
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  template< typename POLICY >
  static void
  createAndLaunch( integer const numSpecies,
                   arrayView1d< integer const > const mobilePrimarySpeciesFlags,
                   globalIndex const rankOffset,
                   string const & dofKey,
                   string const & solverName,
                   FaceManager const & faceManager,
                   ElementRegionManager const & elemManager,
                   BoundaryStencilWrapper const & stencilWrapper,
                   constitutive::reactivefluid::ReactiveThermalCompressibleSinglePhaseFluid & reactiveFluid,
                   real64 const & dt,
                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                   arrayView1d< real64 > const & localRhs )
  {
    constitutiveUpdatePassThru( reactiveFluid, [&]( auto & fluid )
    {
      using FluidType = TYPEOFREF( fluid );
      typename FluidType::KernelWrapper fluidWrapper = fluid.createKernelWrapper();

      singlePhaseReactiveBaseKernels::internal::kernelLaunchSelectorCompSwitch( numSpecies, [&]( auto NS )
      {
        integer constexpr NUM_SPECIES = NS();
        integer constexpr NUM_DOF = 2+NS();
        integer constexpr NUM_EQN = 2+NS();

        using kernelType = DirichletFluxComputeKernel< NUM_SPECIES, NUM_EQN, NUM_DOF, typename FluidType::KernelWrapper, constitutive::ThermalCompressibleSinglePhaseFluid >;

        ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > dofNumberAccessor =
          elemManager.constructArrayViewAccessor< globalIndex, 1 >( dofKey );

        dofNumberAccessor.setName( solverName + "/accessors/" + dofKey );

        typename kernelType::SinglePhaseFlowAccessors singlePhaseFlowAccessors( elemManager, solverName );
        typename kernelType::ReactiveSinglePhaseFlowAccessors reactiveFlowAccessors( elemManager, solverName );
        typename kernelType::ThermalSinglePhaseFlowAccessors thermalSinglePhaseFlowAccessors( elemManager, solverName );
        typename kernelType::SinglePhaseFluidAccessors singlePhaseFluidAccessors( elemManager, solverName );
        typename kernelType::ReactiveSinglePhaseFluidAccessors reactiveFluidAccessors( elemManager, solverName );
        typename kernelType::ThermalReactiveSinglePhaseFluidAccessors thermalFluidAccessors( elemManager, solverName );
        typename kernelType::PermeabilityAccessors permeabilityAccessors( elemManager, solverName );
        typename kernelType::ThermalConductivityAccessors thermalConductivityAccessors( elemManager, solverName );

        kernelType kernel( rankOffset,
                           faceManager,
                           stencilWrapper,
                           fluidWrapper,
                           dofNumberAccessor,
                           singlePhaseFlowAccessors,
                           reactiveFlowAccessors,
                           thermalSinglePhaseFlowAccessors,
                           singlePhaseFluidAccessors,
                           reactiveFluidAccessors,
                           thermalFluidAccessors,
                           permeabilityAccessors,
                           thermalConductivityAccessors,
                           mobilePrimarySpeciesFlags,
                           dt,
                           localMatrix,
                           localRhs );

        kernelType::template launch< POLICY >( stencilWrapper.size(), kernel );
      } );
    } );
  }

};

} // namespace thermalSinglePhaseReactiveFVMKernels

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEREACTIVE_THERMALDIRICHLETFLUXCOMPUTEKERNEL_HPP
