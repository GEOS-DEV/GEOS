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
 * @file ThermalCompositionalMultiphaseFVMKernels.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_THERMALCOMPOSITIONALMULTIPHASEFVMKERNELS_HPP
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_THERMALCOMPOSITIONALMULTIPHASEFVMKERNELS_HPP

#include "constitutive/thermalConductivity/ThermalConductivityBase.hpp"
#include "constitutive/thermalConductivity/ThermalConductivityExtrinsicData.hpp"
#include "physicsSolvers/fluidFlow/IsothermalCompositionalMultiphaseFVMKernels.hpp"

namespace geosx
{

namespace thermalCompositionalMultiphaseFVMKernels
{

using namespace constitutive;

/******************************** PhaseMobilityKernel ********************************/

/**
 * @class PhaseMobilityKernel
 * @tparam NUM_COMP number of fluid components
 * @tparam NUM_PHASE number of fluid phases
 * @brief Define the interface for the property kernel in charge of computing the phase mobilities
 */
template< integer NUM_COMP, integer NUM_PHASE >
class PhaseMobilityKernel : public isothermalCompositionalMultiphaseFVMKernels::PhaseMobilityKernel< NUM_COMP, NUM_PHASE >
{
public:

  using Base = isothermalCompositionalMultiphaseFVMKernels::PhaseMobilityKernel< NUM_COMP, NUM_PHASE >;
  using Base::numPhase;
  using Base::m_phaseDens;
  using Base::m_dPhaseDens;
  using Base::m_phaseVisc;
  using Base::m_dPhaseVisc;
  using Base::m_dPhaseRelPerm_dPhaseVolFrac;

  /**
   * @brief Constructor
   * @param[in] subRegion the element subregion
   * @param[in] fluid the fluid model
   * @param[in] relperm the relperm model
   */
  PhaseMobilityKernel( ObjectManagerBase & subRegion,
                       MultiFluidBase const & fluid,
                       RelativePermeabilityBase const & relperm )
    : Base( subRegion, fluid, relperm ),
    m_dPhaseVolFrac_dTemp( subRegion.getExtrinsicData< extrinsicMeshData::flow::dPhaseVolumeFraction_dTemperature >() ),
    m_dPhaseMob_dTemp( subRegion.getExtrinsicData< extrinsicMeshData::flow::dPhaseMobility_dTemperature >() )
  {}

  /**
   * @brief Compute the phase mobilities in an element
   * @param[in] ei the element index
   */
  GEOSX_HOST_DEVICE
  void compute( localIndex const ei ) const
  {
    using Deriv = multifluid::DerivativeOffset;

    arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > const phaseDens = m_phaseDens[ei][0];
    arraySlice2d< real64 const, multifluid::USD_PHASE_DC - 2 > const dPhaseDens = m_dPhaseDens[ei][0];
    arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > const phaseVisc = m_phaseVisc[ei][0];
    arraySlice2d< real64 const, multifluid::USD_PHASE_DC - 2 > const dPhaseVisc = m_dPhaseVisc[ei][0];
    arraySlice2d< real64 const, relperm::USD_RELPERM_DS - 2 > const dPhaseRelPerm_dPhaseVolFrac = m_dPhaseRelPerm_dPhaseVolFrac[ei][0];
    arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const dPhaseVolFrac_dTemp = m_dPhaseVolFrac_dTemp[ei];

    arraySlice1d< real64, compflow::USD_PHASE - 1 > const dPhaseMob_dTemp = m_dPhaseMob_dTemp[ei];
    LvArray::forValuesInSlice( dPhaseMob_dTemp, []( real64 & val ){ val = 0.0; } );

    Base::compute( ei, [=] GEOSX_HOST_DEVICE ( localIndex const ip,
                                               real64 const & phaseMob,
                                               real64 const & GEOSX_UNUSED_PARAM( dPhaseMob_dPres ),
                                               arraySlice1d< real64 const, compflow::USD_PHASE_DC - 2 > const & GEOSX_UNUSED_PARAM( dPhaseMob_dComp ) )
    {
      // Step 1: compute the derivative of relPerm[ip] wrt temperature
      real64 dRelPerm_dT = 0.0;
      for( integer jp = 0; jp < numPhase; ++jp )
      {
        dRelPerm_dT += dPhaseRelPerm_dPhaseVolFrac[ip][jp] * dPhaseVolFrac_dTemp[jp];
      }

      // Step 2: compute the derivative of phaseMob[ip] wrt temperature
      dPhaseMob_dTemp[ip] = dRelPerm_dT * phaseDens[ip] / phaseVisc[ip]
                            + phaseMob * (dPhaseDens[ip][Deriv::dT] / phaseDens[ip] - dPhaseVisc[ip][Deriv::dT] / phaseVisc[ip] );
    } );
  }

protected:

  // inputs

  /// Views on thermal derivatives of phase volume fractions
  arrayView2d< real64 const, compflow::USD_PHASE > m_dPhaseVolFrac_dTemp;

  // outputs

  /// Views on thermal derivatives of phase mobilities
  arrayView2d< real64, compflow::USD_PHASE > m_dPhaseMob_dTemp;

};

/**
 * @class PhaseMobilityKernelFactory
 */
class PhaseMobilityKernelFactory
{
public:

  /**
   * @brief Create a new kernel and launch
   * @tparam POLICY the policy used in the RAJA kernel
   * @param[in] numComp the number of fluid components
   * @param[in] numPhase the number of fluid phases
   * @param[in] subRegion the element subregion
   * @param[in] fluid the fluid model
   * @param[in] relperm the relperm model
   */
  template< typename POLICY >
  static void
  createAndLaunch( integer const numComp,
                   integer const numPhase,
                   ObjectManagerBase & subRegion,
                   MultiFluidBase const & fluid,
                   RelativePermeabilityBase const & relperm )
  {
    if( numPhase == 2 )
    {
      isothermalCompositionalMultiphaseBaseKernels::
        internal::kernelLaunchSelectorCompSwitch( numComp, [&] ( auto NC )
      {
        integer constexpr NUM_COMP = NC();
        PhaseMobilityKernel< NUM_COMP, 2 > kernel( subRegion, fluid, relperm );
        PhaseMobilityKernel< NUM_COMP, 2 >::template launch< POLICY >( subRegion.size(), kernel );
      } );
    }
    else if( numPhase == 3 )
    {
      isothermalCompositionalMultiphaseBaseKernels::
        internal::kernelLaunchSelectorCompSwitch( numComp, [&] ( auto NC )
      {
        integer constexpr NUM_COMP = NC();
        PhaseMobilityKernel< NUM_COMP, 3 > kernel( subRegion, fluid, relperm );
        PhaseMobilityKernel< NUM_COMP, 3 >::template launch< POLICY >( subRegion.size(), kernel );
      } );
    }
  }
};


/******************************** FaceBasedAssemblyKernel ********************************/

/**
 * @class FaceBasedAssemblyKernel
 * @tparam NUM_COMP number of fluid components
 * @tparam NUM_DOF number of degrees of freedom
 * @tparam STENCILWRAPPER the type of the stencil wrapper
 * @brief Define the interface for the assembly kernel in charge of flux terms
 */
template< integer NUM_COMP, integer NUM_DOF, typename STENCILWRAPPER >
class FaceBasedAssemblyKernel : public isothermalCompositionalMultiphaseFVMKernels::FaceBasedAssemblyKernel< NUM_COMP, NUM_DOF, STENCILWRAPPER >
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

  using AbstractBase = isothermalCompositionalMultiphaseFVMKernels::FaceBasedAssemblyKernelBase;
  using DofNumberAccessor = AbstractBase::DofNumberAccessor;
  using CompFlowAccessors = AbstractBase::CompFlowAccessors;
  using MultiFluidAccessors = AbstractBase::MultiFluidAccessors;
  using CapPressureAccessors = AbstractBase::CapPressureAccessors;
  using PermeabilityAccessors = AbstractBase::PermeabilityAccessors;

  using AbstractBase::m_dt;
  using AbstractBase::m_numPhases;
  using AbstractBase::m_hasCapPressure;
  using AbstractBase::m_rankOffset;
  using AbstractBase::m_ghostRank;
  using AbstractBase::m_dofNumber;
  using AbstractBase::m_gravCoef;
  using AbstractBase::m_phaseMob;
  using AbstractBase::m_dPhaseMassDens;
  using AbstractBase::m_phaseCompFrac;
  using AbstractBase::m_dPhaseCompFrac;
  using AbstractBase::m_dCompFrac_dCompDens;
  using AbstractBase::m_dPhaseCapPressure_dPhaseVolFrac;

  using Base = isothermalCompositionalMultiphaseFVMKernels::FaceBasedAssemblyKernel< NUM_COMP, NUM_DOF, STENCILWRAPPER >;
  using Base::numComp;
  using Base::numDof;
  using Base::numEqn;
  using Base::maxNumElems;
  using Base::maxNumConns;
  using Base::maxStencilSize;
  using Base::m_stencilWrapper;
  using Base::m_seri;
  using Base::m_sesri;
  using Base::m_sei;

  using ThermalCompFlowAccessors =
    StencilAccessors< extrinsicMeshData::flow::temperature,
                      extrinsicMeshData::flow::deltaTemperature,
                      extrinsicMeshData::flow::dPhaseMobility_dTemperature,
                      extrinsicMeshData::flow::dPhaseVolumeFraction_dTemperature >;

  using ThermalMultiFluidAccessors =
    StencilMaterialAccessors< MultiFluidBase,
                              extrinsicMeshData::multifluid::phaseEnthalpy,
                              extrinsicMeshData::multifluid::dPhaseEnthalpy >;

  using ThermalConductivityAccessors =
    StencilMaterialAccessors< ThermalConductivityBase,
                              extrinsicMeshData::thermalconductivity::effectiveConductivity >;
  // for now, we treat thermal conductivity explicitly

  /**
   * @brief Constructor for the kernel interface
   * @param[in] numPhases the number of fluid phases
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] hasCapPressure flag specifying whether capillary pressure is used or not
   * @param[in] stencilWrapper reference to the stencil wrapper
   * @param[in] dofNumberAccessor
   * @param[in] compFlowAccessor
   * @param[in] thermalCompFlowAccessors
   * @param[in] multiFluidAccessor
   * @param[in] thermalMultiFluidAccessors
   * @param[in] capPressureAccessors
   * @param[in] permeabilityAccessors
   * @param[in] thermalConductivityAccessors
   * @param[in] dt time step size
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  FaceBasedAssemblyKernel( integer const numPhases,
                           globalIndex const rankOffset,
                           integer const hasCapPressure,
                           STENCILWRAPPER const & stencilWrapper,
                           DofNumberAccessor const & dofNumberAccessor,
                           CompFlowAccessors const & compFlowAccessors,
                           ThermalCompFlowAccessors const & thermalCompFlowAccessors,
                           MultiFluidAccessors const & multiFluidAccessors,
                           ThermalMultiFluidAccessors const & thermalMultiFluidAccessors,
                           CapPressureAccessors const & capPressureAccessors,
                           PermeabilityAccessors const & permeabilityAccessors,
                           ThermalConductivityAccessors const & thermalConductivityAccessors,
                           real64 const & dt,
                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                           arrayView1d< real64 > const & localRhs )
    : Base( numPhases,
            rankOffset,
            hasCapPressure,
            stencilWrapper,
            dofNumberAccessor,
            compFlowAccessors,
            multiFluidAccessors,
            capPressureAccessors,
            permeabilityAccessors,
            dt,
            localMatrix,
            localRhs ),
    m_temp( thermalCompFlowAccessors.get( extrinsicMeshData::flow::temperature {} ) ),
    m_dTemp( thermalCompFlowAccessors.get( extrinsicMeshData::flow::deltaTemperature {} ) ),
    m_dPhaseMob_dTemp( thermalCompFlowAccessors.get( extrinsicMeshData::flow::dPhaseMobility_dTemperature {} ) ),
    m_dPhaseVolFrac_dTemp( thermalCompFlowAccessors.get( extrinsicMeshData::flow::dPhaseVolumeFraction_dTemperature {} ) ),
    m_phaseEnthalpy( thermalMultiFluidAccessors.get( extrinsicMeshData::multifluid::phaseEnthalpy {} ) ),
    m_dPhaseEnthalpy( thermalMultiFluidAccessors.get( extrinsicMeshData::multifluid::dPhaseEnthalpy {} ) ),
    m_thermalConductivity( thermalConductivityAccessors.get( extrinsicMeshData::thermalconductivity::effectiveConductivity {} ) )
  {}

  struct StackVariables : public Base::StackVariables
  {
public:

    GEOSX_HOST_DEVICE
    StackVariables( localIndex const size, localIndex numElems )
      : Base::StackVariables( size, numElems ),
      dCompFlux_dT( size, numComp ),
      energyFlux( 0.0 ),
      dEnergyFlux_dP( size ),
      dEnergyFlux_dT( size ),
      dEnergyFlux_dC( size, numComp )
    {}

    using Base::StackVariables::stencilSize;
    using Base::StackVariables::numFluxElems;
    using Base::StackVariables::transmissibility;
    using Base::StackVariables::dTrans_dPres;
    using Base::StackVariables::dofColIndices;
    using Base::StackVariables::localFlux;
    using Base::StackVariables::localFluxJacobian;

    // Component fluxes and derivatives

    /// Derivatives of component fluxes wrt temperature
    stackArray2d< real64, maxStencilSize * numComp > dCompFlux_dT;

    // Thermal transmissibility (for now, no derivatives)

    real64 thermalTransmissibility[maxNumConns][2]{};

    // Energy fluxes and derivatives

    /// Energy fluxes
    real64 energyFlux;
    /// Derivatives of energy fluxes wrt pressure
    stackArray1d< real64, maxStencilSize > dEnergyFlux_dP;
    /// Derivatives of energy fluxes wrt temperature
    stackArray1d< real64, maxStencilSize > dEnergyFlux_dT;
    /// Derivatives of energy fluxes wrt component densities
    stackArray2d< real64, maxStencilSize * numComp > dEnergyFlux_dC;

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
    using Deriv = multifluid::DerivativeOffset;

    arraySlice2d< real64 > dCompFlux_dTSlice = stack.dCompFlux_dT.toSlice();

    // I did not find another way to pass the stack.energyFlux to the lambda...
    real64 * energyFluxPtr = &(stack.energyFlux);
    arraySlice1d< real64 > dEnergyFlux_dPSlice = stack.dEnergyFlux_dP.toSlice();
    arraySlice1d< real64 > dEnergyFlux_dTSlice = stack.dEnergyFlux_dT.toSlice();
    arraySlice2d< real64 > dEnergyFlux_dCSlice = stack.dEnergyFlux_dC.toSlice();

    // ***********************************************
    // First, we call the base computeFlux to compute:
    //  1) compFlux and its derivatives (including derivatives wrt temperature),
    //  2) enthalpyFlux and its derivatives (including derivatives wrt temperature
    //
    // Computing dCompFlux_dT and the enthalpyFlux requires quantities already computed in the base computeFlux,
    // such as potGrad, phaseFlux, and the indices of the upwind cell
    // We use the lambda below (called **inside** the phase loop of the base computeFlux) to access these variables
    Base::computeFlux( iconn, stack, [=] GEOSX_HOST_DEVICE ( integer const ip,
                                                             localIndex const k_up,
                                                             localIndex const er_up,
                                                             localIndex const esr_up,
                                                             localIndex const ei_up,
                                                             real64 const & potGrad,
                                                             real64 const & phaseFlux,
                                                             real64 const (&dPhaseFlux_dP)[maxStencilSize],
                                                             real64 const (&dPhaseFlux_dC)[maxStencilSize][numComp] )
    {
      GEOSX_UNUSED_VAR( dPhaseFlux_dP, dPhaseFlux_dC );

      // We are in the loop over phases, ip provides the current phase index.

      // Step 1: compute the derivatives of the mean density at the interface wrt temperature

      stackArray1d< real64, maxNumElems > dDensMean_dT( stack.numFluxElems );

      for( integer i = 0; i < stack.numFluxElems; ++i )
      {
        localIndex const er  = m_seri( iconn, i );
        localIndex const esr = m_sesri( iconn, i );
        localIndex const ei  = m_sei( iconn, i );

        real64 const dDens_dT = m_dPhaseMassDens[er][esr][ei][0][ip][Deriv::dT];
        dDensMean_dT[i] = 0.5 * dDens_dT;
      }

      // Step 2: compute the derivatives of the phase potential difference wrt temperature
      //***** calculation of flux *****

      stackArray1d< real64, maxStencilSize > dPresGrad_dT( stack.stencilSize );
      stackArray1d< real64, maxStencilSize > dGravHead_dT( stack.numFluxElems );

      // compute potential difference MPFA-style
      for( integer i = 0; i < stack.stencilSize; ++i )
      {
        localIndex const er  = m_seri( iconn, i );
        localIndex const esr = m_sesri( iconn, i );
        localIndex const ei  = m_sei( iconn, i );

        // Step 2.1: compute derivative of capillary pressure wrt temperature
        real64 dCapPressure_dT = 0.0;
        if( m_hasCapPressure )
        {
          for( integer jp = 0; jp < m_numPhases; ++jp )
          {
            real64 const dCapPressure_dS = m_dPhaseCapPressure_dPhaseVolFrac[er][esr][ei][0][ip][jp];
            dCapPressure_dT += dCapPressure_dS * m_dPhaseVolFrac_dTemp[er][esr][ei][jp];
          }
        }

        // Step 2.2: compute derivative of phase pressure difference wrt temperature
        dPresGrad_dT[i] += stack.transmissibility[0][i] * (1 - dCapPressure_dT);
        real64 const gravD = stack.transmissibility[0][i] * m_gravCoef[er][esr][ei];

        // Step 2.3: compute derivative of gravity potential difference wrt temperature
        for( integer j = 0; j < stack.numFluxElems; ++j )
        {
          dGravHead_dT[j] += dDensMean_dT[j] * gravD;
        }
      }

      // Step 3: compute the derivatives of the (upwinded) compFlux wrt temperature
      // *** upwinding ***

      // note: the upwinding is done in the base class, which is in charge of
      //       computing the following quantities: potGrad, phaseFlux, k_up, er_up, esr_up, ei_up

      real64 dPhaseFlux_dT[maxStencilSize]{};

      // Step 3.1: compute the derivative of phase flux wrt temperature
      for( integer ke = 0; ke < stack.stencilSize; ++ke )
      {
        dPhaseFlux_dT[ke] += dPresGrad_dT[ke];
      }
      for( integer ke = 0; ke < stack.numFluxElems; ++ke )
      {
        dPhaseFlux_dT[ke] -= dGravHead_dT[ke];
      }
      for( integer ke = 0; ke < stack.stencilSize; ++ke )
      {
        dPhaseFlux_dT[ke] *= m_phaseMob[er_up][esr_up][ei_up][ip];
      }
      dPhaseFlux_dT[k_up] += m_dPhaseMob_dTemp[er_up][esr_up][ei_up][ip] * potGrad;

      // Step 3.2: compute the derivative of component flux wrt temperature

      // slice some constitutive arrays to avoid too much indexing in component loop
      arraySlice1d< real64 const, multifluid::USD_PHASE_COMP - 3 > phaseCompFracSub =
        m_phaseCompFrac[er_up][esr_up][ei_up][0][ip];
      arraySlice2d< real64 const, multifluid::USD_PHASE_COMP_DC - 3 > dPhaseCompFracSub =
        m_dPhaseCompFrac[er_up][esr_up][ei_up][0][ip];

      for( integer ic = 0; ic < numComp; ++ic )
      {
        real64 const ycp = phaseCompFracSub[ic];
        for( integer ke = 0; ke < stack.stencilSize; ++ke )
        {
          dCompFlux_dTSlice[ke][ic] += dPhaseFlux_dT[ke] * ycp;
        }
        dCompFlux_dTSlice[k_up][ic] += phaseFlux * dPhaseCompFracSub[ic][Deriv::dT];
      }

      // Step 4: compute the enthalpy flux

      real64 const enthalpy = m_phaseEnthalpy[er_up][esr_up][ei_up][0][ip];
      *energyFluxPtr += phaseFlux * enthalpy;

      for( integer ke = 0; ke < stack.stencilSize; ++ke )
      {
        dEnergyFlux_dPSlice[ke] += dPhaseFlux_dP[ke] * enthalpy;
        dEnergyFlux_dTSlice[ke] += dPhaseFlux_dT[ke] * enthalpy;

        for( integer jc = 0; jc < numComp; ++jc )
        {
          dEnergyFlux_dCSlice[ke][jc] += dPhaseFlux_dC[ke][jc] * enthalpy;
        }
      }

      dEnergyFlux_dPSlice[k_up] += phaseFlux * m_dPhaseEnthalpy[er_up][esr_up][ei_up][0][ip][Deriv::dP];
      dEnergyFlux_dTSlice[k_up] += phaseFlux * m_dPhaseEnthalpy[er_up][esr_up][ei_up][0][ip][Deriv::dT];

      real64 dProp_dC[numComp]{};
      applyChainRule( numComp,
                      m_dCompFrac_dCompDens[er_up][esr_up][ei_up],
                      m_dPhaseEnthalpy[er_up][esr_up][ei_up][0][ip],
                      dProp_dC,
                      Deriv::dC );
      for( integer jc = 0; jc < numComp; ++jc )
      {
        dEnergyFlux_dCSlice[k_up][jc] += phaseFlux * dProp_dC[jc];
      }
    } );

    // *****************************************************
    // Computation of the conduction term in the energy flux
    // Note that this term is computed using an explicit treatment of conductivity for now
    // Step 1: compute the thermal transmissibilities at this face
    m_stencilWrapper.computeWeights( iconn,
                                     m_thermalConductivity,
                                     m_thermalConductivity, // unused for now
                                     stack.thermalTransmissibility,
                                     stack.dTrans_dPres ); // unused for now

    // Step 2: compute temperature difference at the interface
    for( integer i = 0; i < stack.stencilSize; ++i )
    {
      localIndex const er  = m_seri( iconn, i );
      localIndex const esr = m_sesri( iconn, i );
      localIndex const ei  = m_sei( iconn, i );

      stack.energyFlux += stack.thermalTransmissibility[0][i] * ( m_temp[er][esr][ei] + m_dTemp[er][esr][ei] );
      stack.dEnergyFlux_dT[i] += stack.thermalTransmissibility[0][i];
    }

    // **********************************************************************************
    // At this point, we have computed the energyFlux and the compFlux for all components
    // We have to do two things here:
    // 1) Add dCompFlux_dTemp to the localFluxJacobian of the component mass balance equations
    // 2) Add energyFlux and its derivatives to the localFlux(Jacobian) of the energy balance equation

    // Step 1: add dCompFlux_dTemp to localFluxJacobian
    for( integer ic = 0; ic < numComp; ++ic )
    {
      for( integer ke = 0; ke < stack.stencilSize; ++ke )
      {
        integer const localDofIndexTemp = ke * numDof + numDof - 1;
        stack.localFluxJacobian[ic][localDofIndexTemp]          =  m_dt * stack.dCompFlux_dT[ke][ic];
        stack.localFluxJacobian[numEqn + ic][localDofIndexTemp] = -m_dt * stack.dCompFlux_dT[ke][ic];
      }
    }

    // Step 2: add energyFlux and its derivatives to localFlux and localFluxJacobian
    integer const localRowIndexEnergy = numEqn-1;
    stack.localFlux[localRowIndexEnergy]          =  m_dt * stack.energyFlux;
    stack.localFlux[numEqn + localRowIndexEnergy] = -m_dt * stack.energyFlux;

    for( integer ke = 0; ke < stack.stencilSize; ++ke )
    {
      integer const localDofIndexPres = ke * numDof;
      stack.localFluxJacobian[localRowIndexEnergy][localDofIndexPres]          =  m_dt * stack.dEnergyFlux_dP[ke];
      stack.localFluxJacobian[numEqn + localRowIndexEnergy][localDofIndexPres] = -m_dt * stack.dEnergyFlux_dP[ke];
      integer const localDofIndexTemp = localDofIndexPres + numDof - 1;
      stack.localFluxJacobian[localRowIndexEnergy][localDofIndexTemp]          =  m_dt * stack.dEnergyFlux_dT[ke];
      stack.localFluxJacobian[numEqn + localRowIndexEnergy][localDofIndexTemp] = -m_dt * stack.dEnergyFlux_dT[ke];

      for( integer jc = 0; jc < numComp; ++jc )
      {
        integer const localDofIndexComp = localDofIndexPres + jc + 1;
        stack.localFluxJacobian[localRowIndexEnergy][localDofIndexComp]          =  m_dt * stack.dEnergyFlux_dC[ke][jc];
        stack.localFluxJacobian[numEqn + localRowIndexEnergy][localDofIndexComp] = -m_dt * stack.dEnergyFlux_dC[ke][jc];
      }
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
    // Step 1: assemble the component mass balance equations (i = 0 to i = numDof-2)
    Base::complete( iconn, stack );

    // Step 2: add contribution to residual and jacobian into the energy balance equation
    for( integer i = 0; i < stack.numFluxElems; ++i )
    {
      if( m_ghostRank[m_seri( iconn, i )][m_sesri( iconn, i )][m_sei( iconn, i )] < 0 )
      {
        globalIndex const globalRow = m_dofNumber[m_seri( iconn, i )][m_sesri( iconn, i )][m_sei( iconn, i )];
        localIndex const localRow = LvArray::integerConversion< localIndex >( globalRow - m_rankOffset );
        GEOSX_ASSERT_GE( localRow, 0 );
        GEOSX_ASSERT_GT( AbstractBase::m_localMatrix.numRows(), localRow + numEqn );

        // beware, there is  volume balance eqn in m_localRhs and m_localMatrix!
        RAJA::atomicAdd( parallelDeviceAtomic{}, &AbstractBase::m_localRhs[localRow + numEqn], stack.localFlux[i * numEqn + numEqn-1] );
        AbstractBase::m_localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >
          ( localRow + numEqn,
          stack.dofColIndices.data(),
          stack.localFluxJacobian[i * numEqn + numEqn-1].dataIfContiguous(),
          stack.stencilSize * numDof );
      }
    }
  }

protected:

  /// Views on temperature
  ElementViewConst< arrayView1d< real64 const > > const m_temp;
  ElementViewConst< arrayView1d< real64 const > > const m_dTemp;

  /// Views on derivatives of phase mobilities, volume fractions, mass densities and phase comp fractions
  ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const m_dPhaseMob_dTemp;
  ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const m_dPhaseVolFrac_dTemp;

  /// Views on phase enthalpies
  ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const m_phaseEnthalpy;
  ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const m_dPhaseEnthalpy;

  /// View on thermal conductivity
  ElementViewConst< arrayView3d< real64 const > > m_thermalConductivity;
  // for now, we treat thermal conductivity explicitly

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
   * @param[in] numComps the number of fluid components
   * @param[in] numPhases the number of fluid phases
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofKey string to get the element degrees of freedom numbers
   * @param[in] hasCapPressure flag specifying whether capillary pressure is used or not
   * @param[in] solverName name of the solver (to name accessors)
   * @param[in] elemManager reference to the element region manager
   * @param[in] stencilWrapper reference to the stencil wrapper
   * @param[in] dt time step size
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  template< typename POLICY, typename STENCILWRAPPER >
  static void
  createAndLaunch( integer const numComps,
                   integer const numPhases,
                   globalIndex const rankOffset,
                   string const & dofKey,
                   integer const hasCapPressure,
                   string const & solverName,
                   ElementRegionManager const & elemManager,
                   STENCILWRAPPER const & stencilWrapper,
                   real64 const & dt,
                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                   arrayView1d< real64 > const & localRhs )
  {
    isothermalCompositionalMultiphaseBaseKernels::
      internal::kernelLaunchSelectorCompSwitch( numComps, [&] ( auto NC )
    {
      integer constexpr NUM_COMP = NC();
      integer constexpr NUM_DOF = NC()+2;

      ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > dofNumberAccessor =
        elemManager.constructArrayViewAccessor< globalIndex, 1 >( dofKey );
      dofNumberAccessor.setName( solverName + "/accessors/" + dofKey );

      using KERNEL_TYPE = FaceBasedAssemblyKernel< NUM_COMP, NUM_DOF, STENCILWRAPPER >;
      typename KERNEL_TYPE::CompFlowAccessors compFlowAccessors( elemManager, solverName );
      typename KERNEL_TYPE::ThermalCompFlowAccessors thermalCompFlowAccessors( elemManager, solverName );
      typename KERNEL_TYPE::MultiFluidAccessors multiFluidAccessors( elemManager, solverName );
      typename KERNEL_TYPE::ThermalMultiFluidAccessors thermalMultiFluidAccessors( elemManager, solverName );
      typename KERNEL_TYPE::CapPressureAccessors capPressureAccessors( elemManager, solverName );
      typename KERNEL_TYPE::PermeabilityAccessors permeabilityAccessors( elemManager, solverName );
      typename KERNEL_TYPE::ThermalConductivityAccessors thermalConductivityAccessors( elemManager, solverName );

      KERNEL_TYPE kernel( numPhases, rankOffset, hasCapPressure, stencilWrapper, dofNumberAccessor,
                          compFlowAccessors, thermalCompFlowAccessors, multiFluidAccessors, thermalMultiFluidAccessors,
                          capPressureAccessors, permeabilityAccessors, thermalConductivityAccessors,
                          dt, localMatrix, localRhs );
      KERNEL_TYPE::template launch< POLICY >( stencilWrapper.size(), kernel );
    } );
  }
};


} // namespace thermalCompositionalMultiphaseFVMKernels

} // namespace geosx


#endif //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_THERMALCOMPOSITIONALMULTIPHASEFVMKERNELS_HPP
