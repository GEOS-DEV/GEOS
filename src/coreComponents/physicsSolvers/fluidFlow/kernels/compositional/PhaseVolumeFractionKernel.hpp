/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file PhaseVolumeFractionKernel.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_PHASEVOLUMEFRACTIONKERNEL_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_PHASEVOLUMEFRACTIONKERNEL_HPP

#include "physicsSolvers/fluidFlow/kernels/compositional/PropertyKernelBase.hpp"

namespace geos
{

namespace isothermalCompositionalMultiphaseBaseKernels
{

/******************************** PhaseVolumeFractionKernel ********************************/

/**
 * @class PhaseVolumeFractionKernel
 * @tparam NUM_COMP number of fluid components
 * @tparam NUM_PHASE number of fluid phases
 * @brief Define the interface for the property kernel in charge of computing the phase volume fractions
 */
template< integer NUM_COMP, integer NUM_PHASE >
class PhaseVolumeFractionKernel : public PropertyKernelBase< NUM_COMP >
{
public:

  using Base = PropertyKernelBase< NUM_COMP >;
  using Base::numComp;

  /// Compile time value for the number of phases
  static constexpr integer numPhase = NUM_PHASE;

  /**
   * @brief Constructor
   * @param[in] subRegion the element subregion
   * @param[in] fluid the fluid model
   */
  PhaseVolumeFractionKernel( ObjectManagerBase & subRegion,
                             constitutive::MultiFluidBase const & fluid )
    : Base(),
    m_phaseVolFrac( subRegion.getField< fields::flow::phaseVolumeFraction >() ),
    m_dPhaseVolFrac( subRegion.getField< fields::flow::dPhaseVolumeFraction >() ),
    m_compDens( subRegion.getField< fields::flow::globalCompDensity >() ),
    m_dCompFrac_dCompDens( subRegion.getField< fields::flow::dGlobalCompFraction_dGlobalCompDensity >() ),
    m_phaseFrac( fluid.phaseFraction() ),
    m_dPhaseFrac( fluid.dPhaseFraction() ),
    m_phaseDens( fluid.phaseDensity() ),
    m_dPhaseDens( fluid.dPhaseDensity() )
  {}

  /**
   * @brief Compute the phase volume fractions in an element
   * @tparam FUNC the type of the function that can be used to customize the kernel
   * @param[in] ei the element index
   * @param[in] phaseVolFractionKernelOp the function used to customize the kernel
   */
  template< typename FUNC = NoOpFunc >
  GEOS_HOST_DEVICE
  real64 compute( localIndex const ei,
                  FUNC && phaseVolFractionKernelOp = NoOpFunc{} ) const
  {
    using Deriv = constitutive::multifluid::DerivativeOffset;

    arraySlice1d< real64 const, compflow::USD_COMP - 1 > const compDens = m_compDens[ei];
    arraySlice2d< real64 const, compflow::USD_COMP_DC - 1 > const dCompFrac_dCompDens = m_dCompFrac_dCompDens[ei];
    arraySlice1d< real64 const, constitutive::multifluid::USD_PHASE - 2 > const phaseDens = m_phaseDens[ei][0];
    arraySlice2d< real64 const, constitutive::multifluid::USD_PHASE_DC - 2 > const dPhaseDens = m_dPhaseDens[ei][0];
    arraySlice1d< real64 const, constitutive::multifluid::USD_PHASE - 2 > const phaseFrac = m_phaseFrac[ei][0];
    arraySlice2d< real64 const, constitutive::multifluid::USD_PHASE_DC - 2 > const dPhaseFrac = m_dPhaseFrac[ei][0];
    arraySlice1d< real64, compflow::USD_PHASE - 1 > const phaseVolFrac = m_phaseVolFrac[ei];
    arraySlice2d< real64, compflow::USD_PHASE_DC - 1 > const dPhaseVolFrac = m_dPhaseVolFrac[ei];

    real64 work[numComp]{};

    // compute total density from component partial densities
    real64 totalDensity = 0.0;
    real64 const dTotalDens_dCompDens = 1.0;
    for( integer ic = 0; ic < numComp; ++ic )
    {
      totalDensity += compDens[ic];
    }

    real64 maxDeltaPhaseVolFrac = 0.0;

    for( integer ip = 0; ip < numPhase; ++ip )
    {

      // set the saturation to zero if the phase is absent
      bool const phaseExists = (phaseFrac[ip] > 0);
      if( !phaseExists )
      {
        phaseVolFrac[ip] = 0.;
        for( integer jc = 0; jc < numComp+2; ++jc )
        {
          dPhaseVolFrac[ip][jc] = 0.;
        }
        continue;
      }

      // Expression for volume fractions: S_p = (nu_p / rho_p) * rho_t
      real64 const phaseDensInv = 1.0 / phaseDens[ip];

      // store old saturation to compute change later
      real64 const satOld = phaseVolFrac[ip];

      // compute saturation and derivatives except multiplying by the total density
      phaseVolFrac[ip] = phaseFrac[ip] * phaseDensInv;

      dPhaseVolFrac[ip][Deriv::dP] =
        (dPhaseFrac[ip][Deriv::dP] - phaseVolFrac[ip] * dPhaseDens[ip][Deriv::dP]) * phaseDensInv;

      for( integer jc = 0; jc < numComp; ++jc )
      {
        dPhaseVolFrac[ip][Deriv::dC+jc] =
          (dPhaseFrac[ip][Deriv::dC+jc] - phaseVolFrac[ip] * dPhaseDens[ip][Deriv::dC+jc]) * phaseDensInv;
      }

      // apply chain rule to convert derivatives from global component fractions to densities
      applyChainRuleInPlace( numComp, dCompFrac_dCompDens, dPhaseVolFrac[ip], work, Deriv::dC );

      // call the lambda in the phase loop to allow the reuse of the phaseVolFrac and totalDensity
      // possible use: assemble the derivatives wrt temperature
      phaseVolFractionKernelOp( ip, phaseVolFrac[ip], phaseDensInv, totalDensity );

      // now finalize the computation by multiplying by total density
      for( integer jc = 0; jc < numComp; ++jc )
      {
        dPhaseVolFrac[ip][Deriv::dC+jc] *= totalDensity;
        dPhaseVolFrac[ip][Deriv::dC+jc] += phaseVolFrac[ip] * dTotalDens_dCompDens;
      }

      phaseVolFrac[ip] *= totalDensity;
      dPhaseVolFrac[ip][Deriv::dP] *= totalDensity;

      real64 const deltaPhaseVolFrac = LvArray::math::abs( phaseVolFrac[ip] - satOld );
      if( maxDeltaPhaseVolFrac < deltaPhaseVolFrac )
      {
        maxDeltaPhaseVolFrac = deltaPhaseVolFrac;
      }
    }
    return maxDeltaPhaseVolFrac;
  }

  /**
   * @brief Performs the kernel launch
   * @tparam POLICY the policy used in the RAJA kernels
   * @tparam KERNEL_TYPE the kernel type
   * @param[in] numElems the number of elements
   * @param[inout] kernelComponent the kernel component providing access to the compute function
   */
  template< typename POLICY, typename KERNEL_TYPE >
  static real64
  launch( localIndex const numElems,
          KERNEL_TYPE const & kernelComponent )
  {
    RAJA::ReduceMax< ReducePolicy< POLICY >, real64 > maxDeltaPhaseVolFrac( 0.0 );
    forAll< POLICY >( numElems, [=] GEOS_HOST_DEVICE ( localIndex const ei )
    {
      real64 const deltaPhaseVolFrac = kernelComponent.compute( ei );
      maxDeltaPhaseVolFrac.max( deltaPhaseVolFrac );
    } );
    return maxDeltaPhaseVolFrac.get();
  }

protected:

  // outputs

  /// Views on phase volume fractions
  arrayView2d< real64, compflow::USD_PHASE > m_phaseVolFrac;
  arrayView3d< real64, compflow::USD_PHASE_DC > m_dPhaseVolFrac;

  // inputs

  /// Views on component densities
  arrayView2d< real64 const, compflow::USD_COMP > m_compDens;
  arrayView3d< real64 const, compflow::USD_COMP_DC > m_dCompFrac_dCompDens;

  /// Views on phase fractions
  arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > m_phaseFrac;
  arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > m_dPhaseFrac;

  /// Views on phase densities
  arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > m_phaseDens;
  arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > m_dPhaseDens;

};

/**
 * @class PhaseVolumeFractionKernelFactory
 */
class PhaseVolumeFractionKernelFactory
{
public:

  /**
   * @brief Create a new kernel and launch
   * @tparam POLICY the policy used in the RAJA kernel
   * @param[in] numComp the number of fluid components
   * @param[in] numPhase the number of fluid phases
   * @param[in] subRegion the element subregion
   * @param[in] fluid the fluid model
   */
  template< typename POLICY >
  static real64
  createAndLaunch( integer const numComp,
                   integer const numPhase,
                   ObjectManagerBase & subRegion,
                   constitutive::MultiFluidBase const & fluid )
  {
    real64 maxDeltaPhaseVolFrac = 0.0;
    if( numPhase == 2 )
    {
      internal::kernelLaunchSelectorCompSwitch( numComp, [&] ( auto NC )
      {
        integer constexpr NUM_COMP = NC();
        PhaseVolumeFractionKernel< NUM_COMP, 2 > kernel( subRegion, fluid );
        maxDeltaPhaseVolFrac = PhaseVolumeFractionKernel< NUM_COMP, 2 >::template launch< POLICY >( subRegion.size(), kernel );
      } );
    }
    else if( numPhase == 3 )
    {
      internal::kernelLaunchSelectorCompSwitch( numComp, [&] ( auto NC )
      {
        integer constexpr NUM_COMP = NC();
        PhaseVolumeFractionKernel< NUM_COMP, 3 > kernel( subRegion, fluid );
        maxDeltaPhaseVolFrac = PhaseVolumeFractionKernel< NUM_COMP, 3 >::template launch< POLICY >( subRegion.size(), kernel );
      } );
    }
    return maxDeltaPhaseVolFrac;
  }
};

} // namespace isothermalCompositionalMultiphaseBaseKernels

} // namespace geos


#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_PHASEVOLUMEFRACTIONKERNEL_HPP
