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
 * @file CompositionalMultiphaseFluidUpdates.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_COMPOSITIONALMULTIPHASEFLUIDUPDATES_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_COMPOSITIONALMULTIPHASEFLUIDUPDATES_HPP_

#include "constitutive/fluid/multifluid/compositional/models/ComponentProperties.hpp"

#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
#include "constitutive/fluid/multifluid/MultiFluidUtils.hpp"

namespace geos
{
namespace constitutive
{

/**
 * @brief Kernel wrapper class for CompositionalMultiphaseFluid.
 * @tparam FLASH Class describing the phase equilibrium model
 * @tparam PHASE1 Class describing the phase property models for the first phase.
 * @tparam PHASE2 Class describing the phase property models for the second phase.
 * @tparam PHASE3 Class describing the phase property models for the possible third phase.
 */
template< typename FLASH, typename PHASE1, typename PHASE2, typename PHASE3 >
class CompositionalMultiphaseFluidUpdates final : public MultiFluidBase::KernelWrapper
{
public:
  CompositionalMultiphaseFluidUpdates( compositional::ComponentProperties const & componentProperties,
                                       FLASH const & flash,
                                       PHASE1 const & phase1,
                                       PHASE2 const & phase2,
                                       PHASE3 const & phase3,
                                       arrayView1d< real64 const > const & componentMolarWeight,
                                       bool const useMass,
                                       MultiFluidBase::PhaseProp::ViewType phaseFraction,
                                       MultiFluidBase::PhaseProp::ViewType phaseDensity,
                                       MultiFluidBase::PhaseProp::ViewType phaseMassDensity,
                                       MultiFluidBase::PhaseProp::ViewType phaseViscosity,
                                       MultiFluidBase::PhaseProp::ViewType phaseEnthalpy,
                                       MultiFluidBase::PhaseProp::ViewType phaseInternalEnergy,
                                       MultiFluidBase::PhaseComp::ViewType phaseCompFraction,
                                       MultiFluidBase::FluidProp::ViewType totalDensity );

  GEOS_HOST_DEVICE
  virtual void compute( real64 const pressure,
                        real64 const temperature,
                        arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition,
                        arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseFraction,
                        arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseDensity,
                        arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseMassDensity,
                        arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseEnthalpy,
                        arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseInternalEnergy,
                        arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseViscosity,
                        arraySlice2d< real64, multifluid::USD_PHASE_COMP-2 > const & phaseCompFraction,
                        real64 & totalDensity ) const override;

  GEOS_HOST_DEVICE
  virtual void compute( real64 const pressure,
                        real64 const temperature,
                        arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition,
                        MultiFluidBase::PhaseProp::SliceType const phaseFraction,
                        MultiFluidBase::PhaseProp::SliceType const phaseDensity,
                        MultiFluidBase::PhaseProp::SliceType const phaseMassDensity,
                        MultiFluidBase::PhaseProp::SliceType const phaseViscosity,
                        MultiFluidBase::PhaseProp::SliceType const phaseEnthalpy,
                        MultiFluidBase::PhaseProp::SliceType const phaseInternalEnergy,
                        MultiFluidBase::PhaseComp::SliceType const phaseCompFraction,
                        MultiFluidBase::FluidProp::SliceType const totalDensity ) const override;

  GEOS_HOST_DEVICE
  virtual void update( localIndex const k,
                       localIndex const q,
                       real64 const pressure,
                       real64 const temperature,
                       arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition ) const override;

private:
  // The component properties
  compositional::ComponentProperties::KernelWrapper m_componentProperties;

  // Flash kernel wrapper
  typename FLASH::KernelWrapper m_flash;

  // Phase model kernel wrappers
  typename PHASE1::KernelWrapper m_phase1;
  typename PHASE2::KernelWrapper m_phase2;
  typename PHASE3::KernelWrapper m_phase3;
};

template< typename FLASH, typename PHASE1, typename PHASE2, typename PHASE3 >
CompositionalMultiphaseFluidUpdates< FLASH, PHASE1, PHASE2, PHASE3 >::
CompositionalMultiphaseFluidUpdates( compositional::ComponentProperties const & componentProperties,
                                     FLASH const & flash,
                                     PHASE1 const & phase1,
                                     PHASE2 const & phase2,
                                     PHASE3 const & phase3,
                                     arrayView1d< real64 const > const & componentMolarWeight,
                                     bool const useMass,
                                     MultiFluidBase::PhaseProp::ViewType phaseFraction,
                                     MultiFluidBase::PhaseProp::ViewType phaseDensity,
                                     MultiFluidBase::PhaseProp::ViewType phaseMassDensity,
                                     MultiFluidBase::PhaseProp::ViewType phaseViscosity,
                                     MultiFluidBase::PhaseProp::ViewType phaseEnthalpy,
                                     MultiFluidBase::PhaseProp::ViewType phaseInternalEnergy,
                                     MultiFluidBase::PhaseComp::ViewType phaseCompFraction,
                                     MultiFluidBase::FluidProp::ViewType totalDensity ):
  MultiFluidBase::KernelWrapper( componentMolarWeight,
                                 useMass,
                                 std::move( phaseFraction ),
                                 std::move( phaseDensity ),
                                 std::move( phaseMassDensity ),
                                 std::move( phaseViscosity ),
                                 std::move( phaseEnthalpy ),
                                 std::move( phaseInternalEnergy ),
                                 std::move( phaseCompFraction ),
                                 std::move( totalDensity ) ),
  m_componentProperties( componentProperties.createKernelWrapper() ),
  m_flash( flash.createKernelWrapper() ),
  m_phase1( phase1.createKernelWrapper()),
  m_phase2( phase2.createKernelWrapper()),
  m_phase3( phase3.createKernelWrapper())
{}

template< typename FLASH, typename PHASE1, typename PHASE2, typename PHASE3 >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void
CompositionalMultiphaseFluidUpdates< FLASH, PHASE1, PHASE2, PHASE3 >::compute(
  real64 const pressure,
  real64 const temperature,
  arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition,
  arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseFrac,
  arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseDens,
  arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseMassDens,
  arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseVisc,
  arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseEnthalpy,
  arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseInternalEnergy,
  arraySlice2d< real64, multifluid::USD_PHASE_COMP - 2 > const & phaseCompFrac,
  real64 & totalDens ) const
{
  GEOS_UNUSED_VAR( phaseEnthalpy, phaseInternalEnergy );
  GEOS_UNUSED_VAR( pressure, temperature );

  integer constexpr maxNumComp = MultiFluidBase::MAX_NUM_COMPONENTS;
  integer constexpr maxNumPhase = MultiFluidBase::MAX_NUM_PHASES;
  integer const numComp = numComponents();
  integer const numPhase = numPhases();

  // 1. Convert input mass fractions to mole fractions and keep derivatives

  real64 compMoleFrac[maxNumComp]{};
  if( m_useMass )
  {
    convertToMoleFractions< maxNumComp >( composition, compMoleFrac );
  }
  else
  {
    for( integer ic = 0; ic < numComp; ++ic )
    {
      compMoleFrac[ic] = composition[ic];
    }
  }

  // 2. Calculate values

  for( integer ip = 0; ip < numPhase; ++ip )
  {
    phaseFrac[ip] = 1.0 / numPhase;     // TODO
    phaseDens[ip] = 40.0;               // TODO
    phaseMassDens[ip] = 1000.0;         // TODO
    phaseVisc[ip] = 0.001;              // TODO
    for( integer jc = 0; jc < numComp; ++jc )
    {
      phaseCompFrac[ip][jc] = compMoleFrac[jc];  // TODO
    }
  }

  // 4. if mass variables used instead of molar, perform the conversion

  if( m_useMass )
  {

    // unfortunately here, we have to copy the molecular weight coming from PVT package...
    real64 phaseMolecularWeight[maxNumPhase]{};
    for( integer ip = 0; ip < numPhase; ++ip )
    {
      phaseMolecularWeight[ip] = 1.0;
    }

    // convert mole fractions to mass fractions
    convertToMassFractions< maxNumComp >( phaseMolecularWeight,
                                          phaseFrac,
                                          phaseCompFrac );

  }

  // 5. Compute total fluid mass/molar density

  computeTotalDensity< maxNumComp, maxNumPhase >( phaseFrac,
                                                  phaseDens,
                                                  totalDens );
}

template< typename FLASH, typename PHASE1, typename PHASE2, typename PHASE3 >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void
CompositionalMultiphaseFluidUpdates< FLASH, PHASE1, PHASE2, PHASE3 >::compute(
  real64 const pressure,
  real64 const temperature,
  arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition,
  MultiFluidBase::PhaseProp::SliceType const phaseFraction,
  MultiFluidBase::PhaseProp::SliceType const phaseDensity,
  MultiFluidBase::PhaseProp::SliceType const phaseMassDensity,
  MultiFluidBase::PhaseProp::SliceType const phaseViscosity,
  MultiFluidBase::PhaseProp::SliceType const phaseEnthalpy,
  MultiFluidBase::PhaseProp::SliceType const phaseInternalEnergy,
  MultiFluidBase::PhaseComp::SliceType const phaseCompFraction,
  MultiFluidBase::FluidProp::SliceType const totalDensity ) const
{
  GEOS_UNUSED_VAR( phaseEnthalpy, phaseInternalEnergy );
  GEOS_UNUSED_VAR( pressure, temperature );

  using Deriv = multifluid::DerivativeOffset;

  integer constexpr maxNumComp = MultiFluidBase::MAX_NUM_COMPONENTS;
  integer constexpr maxNumPhase = MultiFluidBase::MAX_NUM_PHASES;
  integer const numComp = numComponents();
  integer const numPhase = numPhases();

  // 1. Convert input mass fractions to mole fractions and keep derivatives

  real64 compMoleFrac[maxNumComp]{};
  real64 dCompMoleFrac_dCompMassFrac[maxNumComp][maxNumComp]{};

  if( m_useMass )
  {
    // convert mass fractions to mole fractions
    convertToMoleFractions( composition,
                            compMoleFrac,
                            dCompMoleFrac_dCompMassFrac );
  }
  else
  {
    for( integer ic = 0; ic < numComp; ++ic )
    {
      compMoleFrac[ic] = composition[ic];
    }
  }

  // 2. Calculate values


  // 3. Extract phase split, phase properties and derivatives from result

  for( integer ip = 0; ip < numPhase; ++ip )
  {
    phaseFraction.value[ip] = 1.0 / numPhase;
    phaseFraction.derivs[ip][Deriv::dP] = 0.0;
    phaseFraction.derivs[ip][Deriv::dT] = 0.0;

    phaseDensity.value[ip] = 40.0;
    phaseDensity.derivs[ip][Deriv::dP] = 0.0;
    phaseDensity.derivs[ip][Deriv::dT] = 0.0;

    phaseMassDensity.value[ip] = 1000.0;
    phaseMassDensity.derivs[ip][Deriv::dP] = 0.0;
    phaseMassDensity.derivs[ip][Deriv::dT] = 0.0;

    phaseViscosity.value[ip] = 0.001;
    phaseViscosity.derivs[ip][Deriv::dP] = 0.0;
    phaseViscosity.derivs[ip][Deriv::dT] = 0.0;

    for( integer jc = 0; jc < numComp; ++jc )
    {
      phaseFraction.derivs[ip][Deriv::dC+jc] = 0.0;
      phaseDensity.derivs[ip][Deriv::dC+jc] = 0.0;
      phaseMassDensity.derivs[ip][Deriv::dC+jc] = 0.0;
      phaseViscosity.derivs[ip][Deriv::dC+jc] = 0.0;

      phaseCompFraction.value[ip][jc] = compMoleFrac[jc];
      phaseCompFraction.derivs[ip][jc][Deriv::dP] = 0.0;
      phaseCompFraction.derivs[ip][jc][Deriv::dT] = 0.0;

      for( integer ic = 0; ic < numComp; ++ic )
      {
        phaseCompFraction.derivs[ip][ic][Deriv::dC+jc] = 0.0;
      }
    }
  }

  // 4. if mass variables used instead of molar, perform the conversion
  if( m_useMass )
  {

    // unfortunately here, we have to copy the molecular weight coming from PVT package...
    real64 phaseMolecularWeight[maxNumPhase]{};
    real64 dPhaseMolecularWeight[maxNumPhase][maxNumComp+2]{};

    for( integer ip = 0; ip < numPhase; ++ip )
    {
      phaseMolecularWeight[ip] = 1.0;
      dPhaseMolecularWeight[ip][Deriv::dP] = 0.0;
      dPhaseMolecularWeight[ip][Deriv::dT] = 0.0;
      for( integer ic = 0; ic < numComp; ++ic )
      {
        dPhaseMolecularWeight[ip][Deriv::dC+ic] = 0.0;
      }
    }

    convertToMassFractions( dCompMoleFrac_dCompMassFrac,
                            phaseMolecularWeight,
                            dPhaseMolecularWeight,
                            phaseFraction,
                            phaseCompFraction,
                            phaseDensity.derivs,
                            phaseViscosity.derivs,
                            phaseEnthalpy.derivs,
                            phaseInternalEnergy.derivs );
  }

  // 5. Compute total fluid mass/molar density and derivatives

  computeTotalDensity( phaseFraction,
                       phaseDensity,
                       totalDensity );
}

template< typename FLASH, typename PHASE1, typename PHASE2, typename PHASE3 >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void
CompositionalMultiphaseFluidUpdates< FLASH, PHASE1, PHASE2, PHASE3 >::
update( localIndex const k,
        localIndex const q,
        real64 const pressure,
        real64 const temperature,
        arraySlice1d< geos::real64 const, compflow::USD_COMP - 1 > const & composition ) const
{
  compute( pressure,
           temperature,
           composition,
           m_phaseFraction( k, q ),
           m_phaseDensity( k, q ),
           m_phaseMassDensity( k, q ),
           m_phaseViscosity( k, q ),
           m_phaseEnthalpy( k, q ),
           m_phaseInternalEnergy( k, q ),
           m_phaseCompFraction( k, q ),
           m_totalDensity( k, q ) );
}

} /* namespace constitutive */

} /* namespace geos */

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_COMPOSITIONALMULTIPHASEFLUIDUPDATES_HPP_
