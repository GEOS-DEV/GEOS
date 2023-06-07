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
 * @file ElementBasedAssemblyKernel.hpp
 */

#ifndef GEOSX_ELEMENTBASEDASSEMBLYKERNEL_HPP
#define GEOSX_ELEMENTBASEDASSEMBLYKERNEL_HPP

#include "physicsSolvers/fluidFlow/kernels/compositionalMultiphase/ElementBasedAssemblyKernel.hpp"


namespace geos
{

namespace thermalCompositionalMultiphaseBaseKernels
{

using namespace constitutive;

/******************************** ElementBasedAssemblyKernel ********************************/

/**
 * @class ElementBasedAssemblyKernel
 * @tparam NUM_COMP number of fluid components
 * @tparam NUM_DOF number of degrees of freedom
 * @brief Define the interface for the assembly kernel in charge of thermal accumulation and volume balance
 */
template< localIndex NUM_COMP, localIndex NUM_DOF >
class ElementBasedAssemblyKernel : public isothermalCompositionalMultiphaseBaseKernels::ElementBasedAssemblyKernel< NUM_COMP, NUM_DOF >
{
public:

  using Base = isothermalCompositionalMultiphaseBaseKernels::ElementBasedAssemblyKernel< NUM_COMP, NUM_DOF >;
  using Base::numComp;
  using Base::numDof;
  using Base::numEqn;
  using Base::m_numPhases;
  using Base::m_rankOffset;
  using Base::m_dofNumber;
  using Base::m_elemGhostRank;
  using Base::m_volume;
  using Base::m_porosity_n;
  using Base::m_porosity;
  using Base::m_dPoro_dPres;
  using Base::m_dCompFrac_dCompDens;
  using Base::m_phaseVolFrac_n;
  using Base::m_phaseVolFrac;
  using Base::m_dPhaseVolFrac;
  using Base::m_phaseDens_n;
  using Base::m_phaseDens;
  using Base::m_dPhaseDens;
  using Base::m_phaseCompFrac_n;
  using Base::m_phaseCompFrac;
  using Base::m_dPhaseCompFrac;
  using Base::m_localMatrix;
  using Base::m_localRhs;

  /**
   * @brief Constructor
   * @param[in] numPhases the number of fluid phases
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofKey the string key to retrieve the degress of freedom numbers
   * @param[in] subRegion the element subregion
   * @param[in] fluid the fluid model
   * @param[in] solid the solid model
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  ElementBasedAssemblyKernel( geos::localIndex const numPhases,
                              geos::globalIndex const rankOffset,
                              geos::string const dofKey,
                              geos::ElementSubRegionBase const & subRegion,
                              geos::constitutive::MultiFluidBase const & fluid,
                              geos::constitutive::CoupledSolidBase const & solid,
                              geos::CRSMatrixView< geos::real64, geos::globalIndex const > const & localMatrix,
                              geos::arrayView1d< geos::real64 > const & localRhs )
    : Base( numPhases, rankOffset, dofKey, subRegion, fluid, solid, localMatrix, localRhs ),
    m_dPoro_dTemp( solid.getDporosity_dTemperature() ),
    m_phaseInternalEnergy_n( fluid.phaseInternalEnergy_n() ),
    m_phaseInternalEnergy( fluid.phaseInternalEnergy() ),
    m_dPhaseInternalEnergy( fluid.dPhaseInternalEnergy() ),
    m_rockInternalEnergy_n( solid.getInternalEnergy_n() ),
    m_rockInternalEnergy( solid.getInternalEnergy() ),
    m_dRockInternalEnergy_dTemp( solid.getDinternalEnergy_dTemperature() )
  { }

  struct StackVariables : public Base::StackVariables
  {
public:

    GEOS_HOST_DEVICE
    StackVariables()
      : Base::StackVariables()
    { }

    using Base::StackVariables::poreVolume;
    using Base::StackVariables::poreVolume_n;
    using Base::StackVariables::dPoreVolume_dPres;
    using Base::StackVariables::localRow;
    using Base::StackVariables::dofIndices;
    using Base::StackVariables::localResidual;
    using Base::StackVariables::localJacobian;

    // derivative of pore volume wrt temperature
    geos::real64 dPoreVolume_dTemp = 0.0;

    // Solid energy

    /// Solid energy at time n+1
    geos::real64 solidEnergy = 0.0;

    /// Solid energy at the previous converged time step
    geos::real64 solidEnergy_n = 0.0;

    /// Derivative of solid internal energy with respect to pressure
    geos::real64 dSolidEnergy_dPres = 0.0;

    /// Derivative of solid internal energy with respect to temperature
    geos::real64 dSolidEnergy_dTemp = 0.0;

  };

  /**
   * @brief Performs the setup phase for the kernel.
   * @param[in] ei the element index
   * @param[in] stack the stack variables
   */
  GEOS_HOST_DEVICE
  void setup( geos::localIndex const ei,
              StackVariables & stack ) const
  {
    ElementBasedAssemblyKernel< NUM_COMP, NUM_DOF >::setup( ei, stack );

    // derivative of pore volume wrt temperature
    stack.dPoreVolume_dTemp = m_volume[ei] * m_dPoro_dTemp[ei][0];

    // initialize the solid volume
    geos::real64 const solidVolume = m_volume[ei] * ( 1.0 - m_porosity[ei][0] );
    geos::real64 const solidVolume_n = m_volume[ei] * ( 1.0 - m_porosity_n[ei][0] );
    geos::real64 const dSolidVolume_dPres = -m_volume[ei] * m_dPoro_dPres[ei][0];
    geos::real64 const dSolidVolume_dTemp = -stack.dPoreVolume_dTemp;

    // initialize the solid internal energy
    stack.solidEnergy = solidVolume * m_rockInternalEnergy[ei][0];
    stack.solidEnergy_n = solidVolume_n * m_rockInternalEnergy_n[ei][0];
    stack.dSolidEnergy_dPres = dSolidVolume_dPres * m_rockInternalEnergy[ei][0];
    stack.dSolidEnergy_dTemp = solidVolume * m_dRockInternalEnergy_dTemp[ei][0]
                               + dSolidVolume_dTemp * m_rockInternalEnergy[ei][0];
  }

  /**
   * @brief Compute the local accumulation contributions to the residual and Jacobian
   * @tparam FUNC the type of the function that can be used to customize the kernel
   * @param[in] ei the element index
   * @param[inout] stack the stack variables
   */
  GEOS_HOST_DEVICE
  void computeAccumulation( geos::localIndex const ei,
                            StackVariables & stack ) const
  {
    using Deriv = geos::constitutive::multifluid::DerivativeOffset;

    Base::computeAccumulation( ei, stack, [&]( integer const ip,
                                               real64 const & phaseAmount,
                                               real64 const & phaseAmount_n,
                                               real64 const & dPhaseAmount_dP,
                                               real64 const (&dPhaseAmount_dC)[numComp] )
    {
      // We are in the loop over phases, ip provides the current phase index.
      // We have to do two things:
      //   1- Assemble the derivatives of the component mass balance equations with respect to temperature
      //   2- Assemble the phase-dependent part of the accumulation term of the energy equation

      real64 dPhaseInternalEnergy_dC[numComp]{};

      // construct the slices
      arraySlice2d< real64 const, compflow::USD_COMP_DC - 1 > dCompFrac_dCompDens = m_dCompFrac_dCompDens[ei];
      arraySlice1d< real64 const, compflow::USD_PHASE - 1 > phaseVolFrac = m_phaseVolFrac[ei];
      arraySlice2d< real64 const, compflow::USD_PHASE_DC - 1 > dPhaseVolFrac = m_dPhaseVolFrac[ei];
      arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > phaseDens = m_phaseDens[ei][0];
      arraySlice2d< real64 const, multifluid::USD_PHASE_DC - 2 > dPhaseDens = m_dPhaseDens[ei][0];
      arraySlice2d< real64 const, multifluid::USD_PHASE_COMP - 2 > phaseCompFrac = m_phaseCompFrac[ei][0];
      arraySlice3d< real64 const, multifluid::USD_PHASE_COMP_DC - 2 > dPhaseCompFrac = m_dPhaseCompFrac[ei][0];
      arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > phaseInternalEnergy_n = m_phaseInternalEnergy_n[ei][0];
      arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > phaseInternalEnergy = m_phaseInternalEnergy[ei][0];
      arraySlice2d< real64 const, multifluid::USD_PHASE_DC - 2 > dPhaseInternalEnergy = m_dPhaseInternalEnergy[ei][0];

      // Step 1: assemble the derivatives of the component mass balance equations with respect to temperature

      real64 const dPhaseAmount_dT = stack.dPoreVolume_dTemp * phaseVolFrac[ip] * phaseDens[ip]
                                     + stack.poreVolume * ( dPhaseVolFrac[ip][Deriv::dT] * phaseDens[ip] + phaseVolFrac[ip] * dPhaseDens[ip][Deriv::dT] );
      for( integer ic = 0; ic < numComp; ++ic )
      {
        stack.localJacobian[ic][numDof - 1] += dPhaseAmount_dT * phaseCompFrac[ip][ic]
                                               + phaseAmount * dPhaseCompFrac[ip][ic][Deriv::dT];
      }

      // Step 2: assemble the phase-dependent part of the accumulation term of the energy equation

      real64 const phaseEnergy = phaseAmount * phaseInternalEnergy[ip];
      real64 const phaseEnergy_n = phaseAmount_n * phaseInternalEnergy_n[ip];
      real64 const dPhaseEnergy_dP = dPhaseAmount_dP * phaseInternalEnergy[ip]
                                     + phaseAmount * dPhaseInternalEnergy[ip][Deriv::dP];
      real64 const dPhaseEnergy_dT = dPhaseAmount_dT * phaseInternalEnergy[ip]
                                     + phaseAmount * dPhaseInternalEnergy[ip][Deriv::dT];

      // local accumulation
      stack.localResidual[numEqn - 1] += phaseEnergy - phaseEnergy_n;

      // derivatives w.r.t. pressure and temperature
      stack.localJacobian[numEqn - 1][0] += dPhaseEnergy_dP;
      stack.localJacobian[numEqn - 1][numDof - 1] += dPhaseEnergy_dT;

      // derivatives w.r.t. component densities
      applyChainRule( numComp, dCompFrac_dCompDens, dPhaseInternalEnergy[ip], dPhaseInternalEnergy_dC, Deriv::dC );
      for( integer jc = 0; jc < numComp; ++jc )
      {
        stack.localJacobian[numEqn - 1][jc + 1] += phaseInternalEnergy[ip] * dPhaseAmount_dC[jc]
                                                   + dPhaseInternalEnergy_dC[jc] * phaseAmount;
      }
    } );

    // Step 3: assemble the solid part of the accumulation term

    // local accumulation and derivatives w.r.t. pressure and temperature
    stack.localResidual[numEqn - 1] += stack.solidEnergy - stack.solidEnergy_n;
    stack.localJacobian[numEqn - 1][0] += stack.dSolidEnergy_dPres;
    stack.localJacobian[numEqn - 1][numDof - 1] += stack.dSolidEnergy_dTemp;

  }

  /**
   * @brief Compute the local volume balance contributions to the residual and Jacobian
   * @tparam FUNC the type of the function that can be used to customize the kernel
   * @param[in] ei the element index
   * @param[inout] stack the stack variables
   */
  GEOS_HOST_DEVICE
  void computeVolumeBalance( localIndex const ei,
                             StackVariables & stack ) const
  {
    using Deriv = multifluid::DerivativeOffset;

    Base::computeVolumeBalance( ei, stack, [&]( real64 const & oneMinusPhaseVolFraction )
    {
      GEOS_UNUSED_VAR( oneMinusPhaseVolFraction );

      arraySlice2d< real64 const, compflow::USD_PHASE_DC - 1 > dPhaseVolFrac = m_dPhaseVolFrac[ei];

      for( integer ip = 0; ip < m_numPhases; ++ip )
      {
        stack.localJacobian[numEqn - 2][numDof - 1] -= dPhaseVolFrac[ip][Deriv::dT];
      }
    } );
  }

  GEOS_HOST_DEVICE
  void complete( localIndex const ei,
                 StackVariables & stack ) const
  {
    // Step 1: assemble the component mass balance equations and volume balance equations
    Base::complete( ei, stack );

    // Step 2: assemble the energy equation
    m_localRhs[stack.localRow + numEqn - 1] += stack.localResidual[numEqn - 1];
    m_localMatrix.template addToRow< serialAtomic >( stack.localRow + numEqn - 1,
                                                     stack.dofIndices,
                                                     stack.localJacobian[numEqn - 1],
                                                     numDof );
  }

protected:

  /// View on derivative of porosity w.r.t temperature
  arrayView2d< real64 const > const m_dPoro_dTemp;

  /// Views on phase internal energy
  arrayView3d< real64 const, multifluid::USD_PHASE > m_phaseInternalEnergy_n;
  arrayView3d< real64 const, multifluid::USD_PHASE > m_phaseInternalEnergy;
  arrayView4d< real64 const, multifluid::USD_PHASE_DC > m_dPhaseInternalEnergy;

  /// Views on rock internal energy
  arrayView2d< real64 const > m_rockInternalEnergy_n;
  arrayView2d< real64 const > m_rockInternalEnergy;
  arrayView2d< real64 const > m_dRockInternalEnergy_dTemp;

};

/**
 * @class ElementBasedAssemblyKernelFactory
 */
class ElementBasedAssemblyKernelFactory
{
public:

  /**
   * @brief Create a new kernel and launch
   * @tparam POLICY the policy used in the RAJA kernel
   * @param[in] numComps the number of fluid components
   * @param[in] numPhases the number of fluid phases
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofKey the string key to retrieve the degress of freedom numbers
   * @param[in] subRegion the element subregion
   * @param[in] fluid the fluid model
   * @param[in] solid the solid model
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  template< typename POLICY >
  static void
  createAndLaunch( localIndex const numComps,
                   localIndex const numPhases,
                   globalIndex const rankOffset,
                   string const dofKey,
                   ElementSubRegionBase const & subRegion,
                   MultiFluidBase const & fluid,
                   CoupledSolidBase const & solid,
                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                   arrayView1d< real64 > const & localRhs )
  {
    isothermalCompositionalMultiphaseBaseKernels::
      internal::kernelLaunchSelectorCompSwitch( numComps, [&]( auto NC )
    {
      localIndex constexpr NUM_COMP = NC();
      localIndex constexpr NUM_DOF = NC() + 2;
      ElementBasedAssemblyKernel< NUM_COMP, NUM_DOF >
      kernel( numPhases, rankOffset, dofKey, subRegion, fluid, solid, localMatrix, localRhs );
      ElementBasedAssemblyKernel< NUM_COMP, NUM_DOF >::template
      launch< POLICY, ElementBasedAssemblyKernel< NUM_COMP, NUM_DOF > >( subRegion.size(), kernel );
    } );
  }

};

}

}

#endif //GEOSX_ELEMENTBASEDASSEMBLYKERNEL_HPP
