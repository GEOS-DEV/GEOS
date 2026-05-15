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
 * @file ThermalDirichletFluxComputeKernel_impl.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_THERMALDIRICHLETFLUXCOMPUTEKERNEL_IMPL_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_THERMALDIRICHLETFLUXCOMPUTEKERNEL_IMPL_HPP

#include "ThermalDirichletFluxComputeKernel.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/DirichletFluxComputeKernel_impl.hpp"

namespace geos
{
namespace thermalCompositionalMultiphaseFVMKernels
{

template< typename FLUIDWRAPPER, integer NUM_COMP, integer NUM_DOF >
DirichletFluxComputeKernel< FLUIDWRAPPER, NUM_COMP, NUM_DOF >::
DirichletFluxComputeKernel( integer const numPhases,
                            globalIndex const rankOffset,
                            FaceManager const & faceManager,
                            BoundaryStencilWrapper const & stencilWrapper,
                            FLUIDWRAPPER const & fluidWrapper,
                            DofNumberAccessor const & dofNumberAccessor,
                            CompFlowAccessors const & compFlowAccessors,
                            ThermalCompFlowAccessors const & thermalCompFlowAccessors,
                            MultiFluidAccessors const & multiFluidAccessors,
                            ThermalMultiFluidAccessors const & thermalMultiFluidAccessors,
                            CapPressureAccessors const & capPressureAccessors,
                            PermeabilityAccessors const & permeabilityAccessors,
                            ThermalConductivityAccessors const & thermalConductivityAccessors,
                            real64 const dt,
                            CRSMatrixView< real64, globalIndex const > const & localMatrix,
                            arrayView1d< real64 > const & localRhs,
                            BitFlags< isothermalCompositionalMultiphaseFVMKernels::KernelFlags > kernelFlags )
  : Base( numPhases,
          rankOffset,
          faceManager,
          stencilWrapper,
          fluidWrapper,
          dofNumberAccessor,
          compFlowAccessors,
          multiFluidAccessors,
          capPressureAccessors,
          permeabilityAccessors,
          dt,
          localMatrix,
          localRhs,
          kernelFlags ),
  m_temp( thermalCompFlowAccessors.get( fields::flow::temperature {} ) ),
  m_phaseEnthalpy( thermalMultiFluidAccessors.get( fields::multifluid::phaseEnthalpy {} ) ),
  m_dPhaseEnthalpy( thermalMultiFluidAccessors.get( fields::multifluid::dPhaseEnthalpy {} ) ),
  m_thermalConductivity( thermalConductivityAccessors.get( fields::thermalconductivity::effectiveConductivity {} ) )
{}

template< typename FLUIDWRAPPER, integer NUM_COMP, integer NUM_DOF >
GEOS_HOST_DEVICE
void DirichletFluxComputeKernel< FLUIDWRAPPER, NUM_COMP, NUM_DOF >::
computeFlux( localIndex const iconn,
             StackVariables & stack ) const
{
  using Order = BoundaryStencil::Order;
  using Deriv = constitutive::multifluid::DerivativeOffset;

  // ***********************************************
  // First, we call the base computeFlux to compute:
  //  1) compFlux and its derivatives (including derivatives wrt temperature),
  //  2) enthalpy part of energyFlux  and its derivatives (including derivatives wrt temperature)
  //
  // Computing dCompFlux_dT and the enthalpy flux requires quantities already computed in the base computeFlux,
  // such as potGrad, phaseFlux, and the indices of the upwind cell
  // We use the lambda below (called **inside** the phase loop of the base computeFlux) to access these variables
  Base::computeFlux( iconn, stack, [&] ( integer const ip,
                                         localIndex const er,
                                         localIndex const esr,
                                         localIndex const ei,
                                         localIndex const kf,
                                         real64 const f,   // potGrad times trans
                                         real64 const facePhaseMob,
                                         arraySlice1d< const real64, constitutive::multifluid::USD_PHASE - 2 > const & facePhaseEnthalpy,
                                         arraySlice2d< const real64, constitutive::multifluid::USD_PHASE_COMP-2 > const & facePhaseCompFrac,
                                         real64 const phaseFlux,
                                         real64 const dPhaseFlux_dP,
                                         real64 const (&dPhaseFlux_dC)[numComp] )
  {
    // We are in the loop over phases, ip provides the current phase index.

    // Step 1: compute the derivatives of the mean density at the interface wrt temperature

    real64 const dDensMean_dT = 0.5 * m_dPhaseMassDens[er][esr][ei][0][ip][Deriv::dT];

    // Step 2: compute the derivatives of the phase potential difference wrt temperature
    //***** calculation of flux *****

    real64 const dF_dT = -stack.transmissibility * dDensMean_dT * ( m_gravCoef[er][esr][ei] - m_faceGravCoef[kf] );

    // Step 3: compute the derivatives of the (upwinded) compFlux wrt temperature
    // *** upwinding ***

    // note: the upwinding is done in the base class, which is in charge of
    //       computing the following quantities: potGrad, phaseFlux
    // It is easier to hard-code the if/else because it is difficult to address elem and face variables in a uniform way

    if( f >= 0 )   // the element is upstream
    {

      // Step 3.1.a: compute the derivative of phase flux wrt temperature
      real64 const dPhaseFlux_dT = m_phaseMob[er][esr][ei][ip] * dF_dT + m_dPhaseMob[er][esr][ei][ip][Deriv::dT] * f;

      // Step 3.2.a: compute the derivative of component flux wrt temperature

      // slice some constitutive arrays to avoid too much indexing in component loop
      arraySlice1d< real64 const, constitutive::multifluid::USD_PHASE_COMP - 3 > phaseCompFracSub =
        m_phaseCompFrac[er][esr][ei][0][ip];
      arraySlice2d< real64 const, constitutive::multifluid::USD_PHASE_COMP_DC - 3 > dPhaseCompFracSub =
        m_dPhaseCompFrac[er][esr][ei][0][ip];

      for( integer ic = 0; ic < numComp; ++ic )
      {
        real64 const ycp = phaseCompFracSub[ic];
        stack.dCompFlux_dT[ic] += dPhaseFlux_dT * ycp + phaseFlux * dPhaseCompFracSub[ic][Deriv::dT];
      }

      // Step 3.3.a: compute the enthalpy flux

      real64 const enthalpy = m_phaseEnthalpy[er][esr][ei][0][ip];
      stack.energyFlux += phaseFlux * enthalpy;
      stack.dEnergyFlux_dP += dPhaseFlux_dP * enthalpy + phaseFlux * m_dPhaseEnthalpy[er][esr][ei][0][ip][Deriv::dP];
      stack.dEnergyFlux_dT += dPhaseFlux_dT * enthalpy + phaseFlux * m_dPhaseEnthalpy[er][esr][ei][0][ip][Deriv::dT];

      real64 dProp_dC[numComp]{};
      applyChainRule( numComp,
                      m_dCompFrac_dCompDens[er][esr][ei],
                      m_dPhaseEnthalpy[er][esr][ei][0][ip],
                      dProp_dC,
                      Deriv::dC );
      for( integer jc = 0; jc < numComp; ++jc )
      {
        stack.dEnergyFlux_dC[jc] += dPhaseFlux_dC[jc] * enthalpy + phaseFlux * dProp_dC[jc];
      }

    }
    else   // the face is upstream
    {

      // Step 3.1.b: compute the derivative of phase flux wrt temperature
      real64 const dPhaseFlux_dT = facePhaseMob * dF_dT;

      // Step 3.2.b: compute the derivative of component flux wrt temperature

      for( integer ic = 0; ic < numComp; ++ic )
      {
        real64 const ycp = facePhaseCompFrac[ip][ic];
        stack.dCompFlux_dT[ic] += dPhaseFlux_dT * ycp;
      }

      // Step 3.3.b: compute the enthalpy flux

      real64 const enthalpy = facePhaseEnthalpy[ip];
      stack.energyFlux += phaseFlux * enthalpy;
      stack.dEnergyFlux_dP += dPhaseFlux_dP * enthalpy;
      stack.dEnergyFlux_dT += dPhaseFlux_dT * enthalpy;
      for( integer jc = 0; jc < numComp; ++jc )
      {
        stack.dEnergyFlux_dC[jc] += dPhaseFlux_dC[jc] * enthalpy;
      }

    }

  } );

  // *****************************************************
  // Computation of the conduction term in the energy flux
  // Note that the phase enthalpy term in the energy was computed above
  // Note that this term is computed using an explicit treatment of conductivity for now

  // Step 1: compute the thermal transmissibilities at this face
  // Below, the thermal conductivity used to compute (explicitly) the thermal conducivity
  // To avoid modifying the signature of the "computeWeights" function for now, we pass m_thermalConductivity twice
  // TODO: modify computeWeights to accomodate explicit coefficients
  real64 thermalTrans = 0.0;
  real64 dThermalTrans_dPerm[3]{};   // not used
  m_stencilWrapper.computeWeights( iconn,
                                   m_thermalConductivity,
                                   thermalTrans,
                                   dThermalTrans_dPerm );

  // Step 2: compute temperature difference at the interface
  stack.energyFlux += thermalTrans
                      * ( m_temp[m_seri( iconn, Order::ELEM )][m_sesri( iconn, Order::ELEM )][m_sei( iconn, Order::ELEM )] - m_faceTemp[m_sei( iconn, Order::FACE )] );
  stack.dEnergyFlux_dT += thermalTrans;


  // **********************************************************************************
  // At this point, we have computed the energyFlux and the compFlux for all components
  // We have to do two things here:
  // 1) Add dCompFlux_dTemp to the localFluxJacobian of the component mass balance equations
  // 2) Add energyFlux and its derivatives to the localFlux(Jacobian) of the energy balance equation

  // Step 1: add dCompFlux_dTemp to localFluxJacobian
  for( integer ic = 0; ic < numComp; ++ic )
  {
    stack.localFluxJacobian[ic][numDof-1] =  m_dt * stack.dCompFlux_dT[ic];
  }

  // Step 2: add energyFlux and its derivatives to localFlux and localFluxJacobian
  integer const localRowIndexEnergy = numEqn-1;
  stack.localFlux[localRowIndexEnergy] =  m_dt * stack.energyFlux;

  stack.localFluxJacobian[localRowIndexEnergy][0] =  m_dt * stack.dEnergyFlux_dP;
  stack.localFluxJacobian[localRowIndexEnergy][numDof-1] =  m_dt * stack.dEnergyFlux_dT;
  for( integer jc = 0; jc < numComp; ++jc )
  {
    stack.localFluxJacobian[localRowIndexEnergy][jc+1] =  m_dt * stack.dEnergyFlux_dC[jc];
  }
}

template< typename FLUIDWRAPPER, integer NUM_COMP, integer NUM_DOF >
GEOS_HOST_DEVICE
void DirichletFluxComputeKernel< FLUIDWRAPPER, NUM_COMP, NUM_DOF >::
complete( localIndex const iconn,
          StackVariables & stack ) const
{
  // Call Case::complete to assemble the component mass balance equations (i = 0 to i = numDof-2)
  // In the lambda, add contribution to residual and jacobian into the energy balance equation
  Base::complete( iconn, stack, [&] ( localIndex const localRow )
  {
    // beware, there is  volume balance eqn in m_localRhs and m_localMatrix!
    RAJA::atomicAdd( parallelDeviceAtomic{}, &AbstractBase::m_localRhs[localRow + numEqn], stack.localFlux[numEqn-1] );
    AbstractBase::m_localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >
      ( localRow + numEqn,
      stack.dofColIndices,
      stack.localFluxJacobian[numEqn-1],
      numDof );
  } );
}

} // namespace thermalCompositionalMultiphaseFVMKernels
} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_THERMALDIRICHLETFLUXCOMPUTEKERNEL_IMPL_HPP
