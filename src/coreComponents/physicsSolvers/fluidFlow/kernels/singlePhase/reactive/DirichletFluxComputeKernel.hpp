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
 * @file DirichletFluxComputeKernel.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEREACTIVE_DIRICHLETFLUXCOMPUTEKERNEL_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEREACTIVE_DIRICHLETFLUXCOMPUTEKERNEL_HPP

#include "constitutive/fluid/reactivefluid/ReactiveSinglePhaseFluid.hpp"
#include "constitutive/fluid/reactivefluid/ReactiveFluidFields.hpp"
#include "constitutive/fluid/reactivefluid/ReactiveFluidSelector.hpp"
#include "physicsSolvers/fluidFlow/kernels/singlePhase/DirichletFluxComputeKernel.hpp"

namespace geos
{

namespace singlePhaseReactiveFVMKernels
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
class DirichletFluxComputeKernel : public singlePhaseFVMKernels::DirichletFluxComputeKernel< NUM_EQN, NUM_DOF,
                                                                                             FLUIDWRAPPER >
{
public:
  /// Compile time value for the number of primary species
  static constexpr integer numSpecies = NUM_SPECIES;

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
  using AbstractBase::m_dens;
  using AbstractBase::m_dDens;

  using Base = singlePhaseFVMKernels::DirichletFluxComputeKernel< NUM_EQN, NUM_DOF,
                                                                  FLUIDWRAPPER >;
  using Base::numDof;
  using Base::numEqn;
  using Base::m_stencilWrapper;
  using Base::m_seri;
  using Base::m_sesri;
  using Base::m_sei;

  using Base::m_facePres;
  using Base::m_fluidWrapper;

  using ReactiveSinglePhaseFlowAccessors =
    StencilAccessors< fields::flow::logPrimarySpeciesConcentration,
                      fields::flow::dMobility_dLogPrimaryConc >;

  using ReactiveSinglePhaseFluidAccessors =
    StencilMaterialAccessors< constitutive::reactivefluid::ReactiveSinglePhaseFluid< BASE_FLUID_TYPE >,
                              fields::reactivefluid::primarySpeciesMobileAggregateConcentration,
                              fields::reactivefluid::dPrimarySpeciesMobileAggregateConcentration_dLogPrimarySpeciesConcentrations >;

  /**
   * @brief Constructor for the kernel interface
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] faceManager the face manager
   * @param[in] stencilWrapper reference to the stencil wrapper
   * @param[in] fluidWrapper reference to the fluid wrapper
   * @param[in] dofNumberAccessor
   * @param[in] singlePhaseFlowAccessors
   * @param[in] singlePhaseFluidAccessors
   * @param[in] permeabilityAccessors
   * @param[in] dt time step size
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
                              SinglePhaseFluidAccessors const & singlePhaseFluidAccessors,
                              ReactiveSinglePhaseFluidAccessors const & reactiveSinglePhaseFluidAccessors,
                              PermeabilityAccessors const & permeabilityAccessors,
                              arrayView1d< integer const > const & mobilePrimarySpeciesFlags,
                              real64 const & solventDensity,
                              real64 const & dt,
                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                              arrayView1d< real64 > const & localRhs )
    : Base( rankOffset,
            faceManager,
            stencilWrapper,
            fluidWrapper,
            dofNumberAccessor,
            singlePhaseFlowAccessors,
            singlePhaseFluidAccessors,
            permeabilityAccessors,
            dt,
            localMatrix,
            localRhs ),
    m_dMob_dLogPrimaryConc( reactiveSinglePhaseFlowAccessors.get( fields::flow::dMobility_dLogPrimaryConc {} ) ),
    m_primarySpeciesMobileAggregateConc( reactiveSinglePhaseFluidAccessors.get( fields::reactivefluid::primarySpeciesMobileAggregateConcentration {} ) ),
    m_dPrimarySpeciesMobileAggregateConc_dLogPrimaryConc( reactiveSinglePhaseFluidAccessors.get(
                                                            fields::reactivefluid::dPrimarySpeciesMobileAggregateConcentration_dLogPrimarySpeciesConcentrations {} ) ),
    m_mobilePrimarySpeciesFlags( mobilePrimarySpeciesFlags ),
    m_solventDensity( solventDensity )
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

    using Base::StackVariables::transmissibility;
    using Base::StackVariables::dofColIndices;
    using Base::StackVariables::localFlux;
    using Base::StackVariables::localFluxJacobian;

  };

  /**
   * @brief Compute the local Dirichlet face flux contributions to the residual and Jacobian
   * @param[in] iconn the connection index
   * @param[inout] stack the stack variables
   * @param[in] compFluxKernelOp the function used to customize the computation of the component fluxes
   */
  template< typename FUNC = NoOpFunc >
  GEOS_HOST_DEVICE
  void computeFlux( localIndex const iconn,
                    StackVariables & stack,
                    FUNC && compFluxKernelOp = NoOpFunc{} ) const
  {
    Base::computeFlux( iconn, stack, [&] ( localIndex const seri,
                                           localIndex const sesri,
                                           localIndex const sei,
                                           localIndex const kf,
                                           real64 const & fluxVal,
                                           real64 const & dFlux_dP,
                                           real64 const & mobility_up,
                                           real64 const & dMobility_dP_up,
                                           real64 const & dens_up,
                                           real64 const & dDens_dP_up )
    {
      GEOS_UNUSED_VAR( kf );

      real64 speciesFlux[numSpecies]{};
      real64 dSpeciesFlux_dP[numSpecies]{};
      real64 dSpeciesFlux_dLogConc[numSpecies][numSpecies]{};

      for( integer is = 0; is < numSpecies; ++is )
      {
        real64 const aggregateConcMolarity_i = m_primarySpeciesMobileAggregateConc[seri][sesri][sei][0][is] * m_solventDensity;
        speciesFlux[is] = aggregateConcMolarity_i / dens_up * fluxVal * mobility_up;

        dSpeciesFlux_dP[is] = aggregateConcMolarity_i / dens_up * dFlux_dP * mobility_up
                              + aggregateConcMolarity_i / dens_up * fluxVal * dMobility_dP_up
                              - aggregateConcMolarity_i * fluxVal * mobility_up * dDens_dP_up / (dens_up * dens_up);

        for( integer js = 0; js < numSpecies; ++js )
        {
          real64 const dAggregateConcMolarity_i_dLogConc_j = m_dPrimarySpeciesMobileAggregateConc_dLogPrimaryConc[seri][sesri][sei][0][is][js] * m_solventDensity;
          dSpeciesFlux_dLogConc[is][js] += dAggregateConcMolarity_i_dLogConc_j / dens_up * fluxVal * mobility_up;
        }
      }

      /// populate local flux vector and derivatives
      for( integer is = 0; is < numSpecies; ++is )
      {
        stack.localFlux[numEqn - numSpecies + is] =  m_dt * speciesFlux[is] * m_mobilePrimarySpeciesFlags[is];
        stack.localFluxJacobian[numEqn - numSpecies + is][0] = m_dt * dSpeciesFlux_dP[is] * m_mobilePrimarySpeciesFlags[is];

        for( integer js = 0; js < numSpecies; ++js )
        {
          stack.localFluxJacobian[numEqn - numSpecies + is][numDof - numSpecies + js] += m_dt * dSpeciesFlux_dLogConc[is][js] * m_mobilePrimarySpeciesFlags[is];
        }
      }

      compFluxKernelOp( seri, sesri, sei, kf, fluxVal, dFlux_dP, mobility_up, dMobility_dP_up );
    } );
  }

  /**
   * @brief Performs the complete phase for the kernel.
   * @param[in] iconn the connection index
   * @param[inout] stack the stack variables
   */
  template< typename FUNC = NoOpFunc >
  GEOS_HOST_DEVICE
  void complete( localIndex const iconn,
                 StackVariables & stack,
                 FUNC && assemblyKernelOp = NoOpFunc{} ) const
  {
    Base::complete( iconn, stack, [&] ( localIndex const localRow )
    {
      for( integer is = 0; is < numSpecies; ++is )
      {
        RAJA::atomicAdd( parallelDeviceAtomic{}, &AbstractBase::m_localRhs[localRow + numEqn - numSpecies + is],
                         stack.localFlux[numEqn - numSpecies + is] );

        AbstractBase::m_localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >
          ( localRow + numEqn - numSpecies + is,
          stack.dofColIndices,
          stack.localFluxJacobian[numEqn - numSpecies + is],
          numDof );
      }

      assemblyKernelOp( localRow );
    } );
  }

protected:

  /// Views on derivatives of fluid mobilities
  ElementViewConst< arrayView2d< real64 const, compflow::USD_FLUID_DC > > const m_dMob_dLogPrimaryConc;

  /// Views on primary species aggregate concentration
  ElementViewConst< arrayView3d< real64 const, constitutive::reactivefluid::USD_SPECIES > > const m_primarySpeciesMobileAggregateConc;

  /// Views on the derivative of primary species mobile aggregate concentration wrt log of primary concentration
  ElementViewConst< arrayView4d< real64 const, constitutive::reactivefluid::USD_SPECIES_DC > > const m_dPrimarySpeciesMobileAggregateConc_dLogPrimaryConc;

  /// Array of flags to indicate mobile primary species
  arrayView1d< integer const > const m_mobilePrimarySpeciesFlags;

  /// Solvent density [kg/m³] used to convert molality [mol/kg] to molarity [mol/m³]
  real64 const m_solventDensity;

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
                   constitutive::reactivefluid::ReactiveCompressibleSinglePhaseFluid & reactiveFluid,
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
        integer constexpr NUM_DOF = 1+NS();
        integer constexpr NUM_EQN = 1+NS();

        using kernelType = DirichletFluxComputeKernel< NUM_SPECIES, NUM_EQN, NUM_DOF, typename FluidType::KernelWrapper,
                                                       constitutive::CompressibleSinglePhaseFluid >;

        ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > dofNumberAccessor =
          elemManager.constructArrayViewAccessor< globalIndex, 1 >( dofKey );

        dofNumberAccessor.setName( solverName + "/accessors/" + dofKey );

        typename kernelType::SinglePhaseFlowAccessors singlePhaseFlowAccessors( elemManager, solverName );
        typename kernelType::ReactiveSinglePhaseFlowAccessors reactiveFlowAccessors( elemManager, solverName );
        typename kernelType::SinglePhaseFluidAccessors singlePhaseFluidAccessors( elemManager, solverName );
        typename kernelType::ReactiveSinglePhaseFluidAccessors reactiveFluidAccessors( elemManager, solverName );
        typename kernelType::PermeabilityAccessors permeabilityAccessors( elemManager, solverName );

        kernelType kernel( rankOffset,
                           faceManager,
                           stencilWrapper,
                           fluidWrapper,
                           dofNumberAccessor,
                           singlePhaseFlowAccessors,
                           reactiveFlowAccessors,
                           singlePhaseFluidAccessors,
                           reactiveFluidAccessors,
                           permeabilityAccessors,
                           mobilePrimarySpeciesFlags,
                           reactiveFluid.solventDensity(),
                           dt,
                           localMatrix,
                           localRhs );

        kernelType::template launch< POLICY >( stencilWrapper.size(), kernel );
      } );
    } );
  }

};

} // namespace singlePhaseReactiveFVMKernels

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEREACTIVE_DIRICHLETFLUXCOMPUTEKERNEL_HPP
