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
 * @file ThermalSinglePhaseFVMKernels.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_THERMALSINGLEPHASEFVMKERNELS_HPP
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_THERMALSINGLEPHASEFVMKERNELS_HPP

#include "constitutive/thermalConductivity/SinglePhaseThermalConductivityBase.hpp"
#include "constitutive/thermalConductivity/ThermalConductivityExtrinsicData.hpp"
#include "constitutive/thermalConductivity/SinglePhaseThermalConductivityExtrinsicData.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseFVMKernels.hpp"

namespace geosx
{

namespace thermalSinglePhaseFVMKernels
{
using namespace constitutive;

/******************************** FaceBasedAssemblyKernel ********************************/

/**
 * @class FaceBasedAssemblyKernel
 * @tparam NUM_DOF number of degrees of freedom
 * @tparam STENCILWRAPPER the type of the stencil wrapper
 * @brief Define the interface for the assembly kernel in charge of flux terms
 */
template< integer NUM_DOF, typename STENCILWRAPPER >
class FaceBasedAssemblyKernel : public singlePhaseFVMKernels::FaceBasedAssemblyKernel< NUM_DOF, STENCILWRAPPER >
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

  using AbstractBase = singlePhaseFVMKernels::FaceBasedAssemblyKernelBase;
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

  using Base = singlePhaseFVMKernels::FaceBasedAssemblyKernel< NUM_DOF, STENCILWRAPPER >;
  using Base::numDof;
  using Base::numEqn;
  using Base::maxNumElems;
  using Base::maxNumConns;
  using Base::maxStencilSize;
  using Base::m_stencilWrapper;
  using Base::m_seri;
  using Base::m_sesri;
  using Base::m_sei;

  using ThermalSinglePhaseFlowAccessors =
    StencilAccessors< extrinsicMeshData::flow::temperature,
                      extrinsicMeshData::flow::dMobility_dTemperature >;

  using ThermalSinglePhaseFluidAccessors =
    StencilMaterialAccessors< SingleFluidBase,
                              extrinsicMeshData::singlefluid::dDensity_dTemperature,
                              extrinsicMeshData::singlefluid::enthalpy,
                              extrinsicMeshData::singlefluid::dEnthalpy_dPressure, 
                              extrinsicMeshData::singlefluid::dEnthalpy_dTemperature >;

  using ThermalConductivityAccessors =
    StencilMaterialAccessors< SinglePhaseThermalConductivityBase,
                              extrinsicMeshData::thermalconductivity::effectiveConductivity >;

  /**
   * @brief Constructor for the kernel interface
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] hasCapPressure flag specifying whether capillary pressure is used or not
   * @param[in] stencilWrapper reference to the stencil wrapper
   * @param[in] dofNumberAccessor accessor for the dofs numbers
   * @param[in] compFlowAccessor accessor for wrappers registered by the solver
   * @param[in] thermalCompFlowAccessors accessor for *thermal* wrappers registered by the solver
   * @param[in] multiFluidAccessor accessor for wrappers registered by the multifluid model
   * @param[in] thermalMultiFluidAccessors accessor for *thermal* wrappers registered by the multifluid model
   * @param[in] capPressureAccessors accessor for wrappers registered by the cap pressure model
   * @param[in] permeabilityAccessors accessor for wrappers registered by the permeability model
   * @param[in] thermalConductivityAccessors accessor for wrappers registered by the thermal conductivity model
   * @param[in] dt time step size
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  FaceBasedAssemblyKernel( globalIndex const rankOffset,
                           STENCILWRAPPER const & stencilWrapper,
                           DofNumberAccessor const & dofNumberAccessor,
                           SinglePhaseFlowAccessors const & singlePhaseFlowAccessors,
                           ThermalSinglePhaseFlowAccessors const & thermalSinglePhaseFlowAccessors,
                           SinglePhaseFluidAccessors const & singlePhaseFluidAccessors,
                           ThermalSinglePhaseFluidAccessors const & thermalSinglePhaseFluidAccessors,
                           PermeabilityAccessors const & permeabilityAccessors,
                           ThermalConductivityAccessors const & thermalConductivityAccessors,
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
    m_temp( thermalSinglePhaseFlowAccessors.get( extrinsicMeshData::flow::temperature {} ) ),
    m_dMob_dTemp( thermalSinglePhaseFlowAccessors.get( extrinsicMeshData::flow::dMobility_dTemperature {} ) ),
    m_dDens_dTemp( thermalSinglePhaseFluidAccessors.get( extrinsicMeshData::singlefluid::dDensity_dTemperature {} ) ),
    m_enthalpy( thermalSinglePhaseFluidAccessors.get( extrinsicMeshData::singlefluid::enthalpy {} ) ),
    m_dEnthalpy_dPres( thermalSinglePhaseFluidAccessors.get( extrinsicMeshData::singlefluid::dEnthalpy_dPressure {} ) ),
    m_dEnthalpy_dTemp( thermalSinglePhaseFluidAccessors.get( extrinsicMeshData::singlefluid::dEnthalpy_dTemperature {} ) ), 
    m_thermalConductivity( thermalConductivityAccessors.get( extrinsicMeshData::thermalconductivity::effectiveConductivity {} ) )
  {}

  struct StackVariables : public Base::StackVariables
  {
public:

    GEOSX_HOST_DEVICE
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

    // Thermal transmissibility (for now, no derivatives)

    real64 thermalTransmissibility[maxNumConns][2]{};

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
  GEOSX_HOST_DEVICE
  void computeFlux( localIndex const iconn,
                    StackVariables & stack ) const
  {
    // ***********************************************
    // First, we call the base computeFlux to compute:
    //  1) compFlux and its derivatives (including derivatives wrt temperature),
    //  2) enthalpy part of energyFlux  and its derivatives (including derivatives wrt temperature)
    //
    // Computing dCompFlux_dT and the enthalpy flux requires quantities already computed in the base computeFlux,
    // such as potGrad, phaseFlux, and the indices of the upwind cell
    // We use the lambda below (called **inside** the phase loop of the base computeFlux) to access these variables
    Base::computeFlux( iconn, stack, [&] ( localIndex const k_up,
                                           localIndex const er_up,
                                           localIndex const esr_up,
                                           localIndex const ei_up,
                                           real64 const & potGrad,
                                           real64 const & fluxVal,
                                           real64 const (&dFlux_dP)[maxStencilSize] )
    {
      // Step 1: compute the derivatives of the mean density at the interface wrt temperature

      stackArray1d< real64, maxNumElems > dDensMean_dT( stack.numFluxElems ); 

      for ( integer ke = 0; ke < stack.numFluxElems; ++ke )
      {
        localIndex const er  = m_seri( iconn, ke );
        localIndex const esr = m_sesri( iconn, ke );
        localIndex const ei  = m_sei( iconn, ke );

        real64 const dDens_dT = m_dDens_dTemp[er][esr][ei][0]; 
        dDensMean_dT[ke] = 0.5 * dDens_dT; 
      }

      // Step 2: compute the derivatives of the potential difference wrt temperature
      //***** calculation of flux *****

      stackArray1d< real64, maxStencilSize > dGravHead_dT( stack.numFluxElems );

      // compute potential difference 
      for( integer i = 0; i < stack.stencilSize; ++i )
      {
        localIndex const er  = m_seri( iconn, i );
        localIndex const esr = m_sesri( iconn, i );
        localIndex const ei  = m_sei( iconn, i );

        // compute derivative of gravity potential difference wrt temperature
        real64 const gravD = stack.transmissibility[0][i] * m_gravCoef[er][esr][ei];

        for( integer ke = 0; ke < stack.numFluxElems; ++ke )
        {
          dGravHead_dT[ke] += dDensMean_dT[ke] * gravD;
        }
      }

      // Step 3: compute the derivatives of the (upwinded) compFlux wrt temperature
      // *** upwinding ***

      real64 dFlux_dT[maxStencilSize]{};

      // Step 3.1: compute the derivative of flux wrt temperature
      for( integer ke = 0; ke < stack.numFluxElems; ++ke )
      {
        dFlux_dT[ke] -= dGravHead_dT[ke];
      }
      for( integer i = 0; i < stack.stencilSize; ++i )
      {
        dFlux_dT[i] *= m_mob[er_up][esr_up][ei_up];
      }

      dFlux_dT[k_up] += m_dMob_dTemp[er_up][esr_up][ei_up] * potGrad;

      // add dFlux_dTemp to localFluxJacobian
      for( integer i = 0; i < stack.stencilSize; ++i )
      {
        localIndex const localDofIndexTemp = i * numDof + numDof - 1; 
        stack.localFluxJacobian[0][localDofIndexTemp]      += m_dt * dFlux_dT[i]; 
        stack.localFluxJacobian[numEqn][localDofIndexTemp] -= m_dt * dFlux_dT[i]; 
      }

      // Step 4: compute the enthalpy flux

      real64 const enthalpy = m_enthalpy[er_up][esr_up][ei_up][0];
      stack.energyFlux += fluxVal * enthalpy;

      for( integer i = 0; i < stack.stencilSize; ++i )
      {
        stack.dEnergyFlux_dP[i] += dFlux_dP[i] * enthalpy;
        stack.dEnergyFlux_dT[i] += dFlux_dT[i] * enthalpy;
      }

      stack.dEnergyFlux_dP[k_up] += fluxVal * m_dEnthalpy_dPres[er_up][esr_up][ei_up][0];
      stack.dEnergyFlux_dT[k_up] += fluxVal * m_dEnthalpy_dTemp[er_up][esr_up][ei_up][0];
    }); 

    // *****************************************************
    // Computation of the conduction term in the energy flux
    // Note that the enthalpy term in the energy was computed above
    // Note that this term is computed using an explicit treatment of conductivity for now

    // Step 1: compute the thermal transmissibilities at this face
    // We follow how the thermal compositional multi-phase solver does to update the thermal transmissibility
    m_stencilWrapper.computeWeights( iconn,
                                     m_thermalConductivity,
                                     m_thermalConductivity, // we have to pass something here, so we just use thermal conductivity
                                     stack.thermalTransmissibility,
                                     stack.dTrans_dPres ); // again, we have to pass something here, but this is unused for now

    // Step 2: compute temperature difference at the interface
    for( integer i = 0; i < stack.stencilSize; ++i )
    {
      localIndex const er  = m_seri( iconn, i );
      localIndex const esr = m_sesri( iconn, i );
      localIndex const ei  = m_sei( iconn, i );

      stack.energyFlux += stack.thermalTransmissibility[0][i] * m_temp[er][esr][ei];
      stack.dEnergyFlux_dT[i] += stack.thermalTransmissibility[0][i];
    }

    // add energyFlux and its derivatives to localFlux and localFluxJacobian
    integer const localRowIndexEnergy = numEqn-1;
    stack.localFlux[localRowIndexEnergy]          += m_dt * stack.energyFlux;
    stack.localFlux[numEqn + localRowIndexEnergy] -= m_dt * stack.energyFlux;

    for( integer i = 0; i < stack.stencilSize; ++i )
    {
      integer const localDofIndexPres = i * numDof;
      stack.localFluxJacobian[localRowIndexEnergy][localDofIndexPres]          =  m_dt * stack.dEnergyFlux_dP[i];
      stack.localFluxJacobian[numEqn + localRowIndexEnergy][localDofIndexPres] = -m_dt * stack.dEnergyFlux_dP[i];
      integer const localDofIndexTemp = localDofIndexPres + numDof - 1;
      stack.localFluxJacobian[localRowIndexEnergy][localDofIndexTemp]          =  m_dt * stack.dEnergyFlux_dT[i];
      stack.localFluxJacobian[numEqn + localRowIndexEnergy][localDofIndexTemp] = -m_dt * stack.dEnergyFlux_dT[i];
    }
  }

  /**
   * @brief Performs the complete phase for the kernel.
   * @param[in] iconn the connection index
   * @param[inout] stack the stack variables
   */
  GEOSX_HOST_DEVICE
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
      RAJA::atomicAdd( parallelDeviceAtomic{}, &AbstractBase::m_localRhs[localRow + numEqn-1], stack.localFlux[i * numEqn + numEqn-1] ); 

      AbstractBase::m_localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( localRow + numEqn-1,
                                                                                        stack.dofColIndices.data(),
                                                                                        stack.localFluxJacobian[i * numEqn + numEqn-1].dataIfContiguous(),
                                                                                        stack.stencilSize * numDof );

    } );
  }

protected:

  /// Views on temperature
  ElementViewConst< arrayView1d< real64 const > > const m_temp;

  /// Views on derivatives of fluid mobilities
  ElementViewConst< arrayView1d< real64 const > > const m_dMob_dTemp;

  /// Views on derivatives of fluid densities
  ElementViewConst< arrayView2d< real64 const > > const m_dDens_dTemp;

  /// Views on phase enthalpies
  ElementViewConst< arrayView2d< real64 const > > const m_enthalpy;
  ElementViewConst< arrayView2d< real64 const > > const m_dEnthalpy_dPres;
  ElementViewConst< arrayView2d< real64 const > > const m_dEnthalpy_dTemp;

  /// View on thermal conductivity
  ElementViewConst< arrayView3d< real64 const > > m_thermalConductivity;

}; 

/**
 * @class FaceBasedAssemblyKernelFactory
 */
class FaceBasedAssemblyKernelFactory
{
public:

  /**
   * @brief Create a new kernel and launch
   * @tparam POLICY the policy used in the RAJA kernel
   * @tparam STENCILWRAPPER the type of the stencil wrapper
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
  createAndLaunch( globalIndex const rankOffset,
                   string const & dofKey,
                   string const & solverName,
                   ElementRegionManager const & elemManager,
                   STENCILWRAPPER const & stencilWrapper,
                   real64 const & dt,
                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                   arrayView1d< real64 > const & localRhs )
  {
    integer constexpr NUM_DOF = 2;

    ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > dofNumberAccessor =
      elemManager.constructArrayViewAccessor< globalIndex, 1 >( dofKey );
    dofNumberAccessor.setName( solverName + "/accessors/" + dofKey );

    using KernelType = FaceBasedAssemblyKernel< NUM_DOF, STENCILWRAPPER >;
    typename KernelType::SinglePhaseFlowAccessors flowAccessors( elemManager, solverName );
    typename KernelType::ThermalSinglePhaseFlowAccessors thermalFlowAccessors( elemManager, solverName );
    typename KernelType::SinglePhaseFluidAccessors fluidAccessors( elemManager, solverName );
    typename KernelType::ThermalSinglePhaseFluidAccessors thermalFluidAccessors( elemManager, solverName );
    typename KernelType::PermeabilityAccessors permAccessors( elemManager, solverName );
    typename KernelType::ThermalConductivityAccessors thermalConductivityAccessors( elemManager, solverName );

    KernelType kernel( rankOffset, stencilWrapper, dofNumberAccessor,
                       flowAccessors, thermalFlowAccessors, fluidAccessors, thermalFluidAccessors,
                       permAccessors, thermalConductivityAccessors,
                       dt, localMatrix, localRhs );
    KernelType::template launch< POLICY >( stencilWrapper.size(), kernel );
  }
};

} // namespace thermalSinglePhaseFVMKernels

} // namespace geosx

#endif //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_THERMALSINGLEPHASEFVMKERNELS_HPP
