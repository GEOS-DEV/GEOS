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
                                       MultiFluidBase::PhaseProp::ViewType phaseFrac,
                                       MultiFluidBase::PhaseProp::ViewType phaseDens,
                                       MultiFluidBase::PhaseProp::ViewType phaseMassDensity,
                                       MultiFluidBase::PhaseProp::ViewType phaseVisc,
                                       MultiFluidBase::PhaseProp::ViewType phaseEnthalpy,
                                       MultiFluidBase::PhaseProp::ViewType phaseInternalEnergy,
                                       MultiFluidBase::PhaseComp::ViewType phaseCompFrac,
                                       MultiFluidBase::FluidProp::ViewType totalDensity );

  GEOS_HOST_DEVICE
  virtual void compute( real64 const pressure,
                        real64 const temperature,
                        arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition,
                        arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseFrac,
                        arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseDens,
                        arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseMassDensity,
                        arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseEnthalpy,
                        arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseInternalEnergy,
                        arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseVisc,
                        arraySlice2d< real64, multifluid::USD_PHASE_COMP-2 > const & phaseCompFrac,
                        real64 & totalDensity ) const override;

  GEOS_HOST_DEVICE
  virtual void compute( real64 const pressure,
                        real64 const temperature,
                        arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition,
                        MultiFluidBase::PhaseProp::SliceType const phaseFrac,
                        MultiFluidBase::PhaseProp::SliceType const phaseDens,
                        MultiFluidBase::PhaseProp::SliceType const phaseMassDensity,
                        MultiFluidBase::PhaseProp::SliceType const phaseVisc,
                        MultiFluidBase::PhaseProp::SliceType const phaseEnthalpy,
                        MultiFluidBase::PhaseProp::SliceType const phaseInternalEnergy,
                        MultiFluidBase::PhaseComp::SliceType const phaseCompFrac,
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
                                     MultiFluidBase::PhaseProp::ViewType phaseFrac,
                                     MultiFluidBase::PhaseProp::ViewType phaseDens,
                                     MultiFluidBase::PhaseProp::ViewType phaseMassDensity,
                                     MultiFluidBase::PhaseProp::ViewType phaseVisc,
                                     MultiFluidBase::PhaseProp::ViewType phaseEnthalpy,
                                     MultiFluidBase::PhaseProp::ViewType phaseInternalEnergy,
                                     MultiFluidBase::PhaseComp::ViewType phaseCompFrac,
                                     MultiFluidBase::FluidProp::ViewType totalDensity ):
  MultiFluidBase::KernelWrapper( componentMolarWeight,
                                 useMass,
                                 std::move( phaseFrac ),
                                 std::move( phaseDens ),
                                 std::move( phaseMassDensity ),
                                 std::move( phaseVisc ),
                                 std::move( phaseEnthalpy ),
                                 std::move( phaseInternalEnergy ),
                                 std::move( phaseCompFrac ),
                                 std::move( totalDensity ) ),
  m_componentProperties( componentProperties.createKernelWrapper() ),
  m_flash( flash.createKernelWrapper() ),
  m_phase1( phase1.createKernelWrapper() ),
  m_phase2( phase2.createKernelWrapper() )
{
  GEOS_UNUSED_VAR( phase3 );
}

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

  stackArray1d< real64, maxNumComp > compMoleFrac( numComp );
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

  // 2. Compute phase fractions and phase component fractions
  m_flash.compute( pressure,
                   temperature,
                   compMoleFrac.toSliceConst(),
                   phaseFrac,
                   phaseCompFrac );

  // 3. Calculate the phase densities
  m_phase1.density.compute( pressure,
                            temperature,
                            phaseCompFrac[0].toSliceConst(),
                            phaseDens[0],
                            phaseMassDens[0],
                            m_useMass );
  m_phase2.density.compute( pressure,
                            temperature,
                            phaseCompFrac[1].toSliceConst(),
                            phaseDens[1],
                            phaseMassDens[1],
                            m_useMass );

  // 4. Calculate the phase viscosities
  m_phase1.viscosity.compute( pressure,
                              temperature,
                              phaseCompFrac[0].toSliceConst(),
                              phaseMassDens[0],
                              phaseVisc[0],
                              m_useMass );
  m_phase2.viscosity.compute( pressure,
                              temperature,
                              phaseCompFrac[1].toSliceConst(),
                              phaseMassDens[1],
                              phaseVisc[1],
                              m_useMass );

  // 5. if mass variables used instead of molar, perform the conversion

  if( m_useMass )
  {

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

  // 6. Compute total fluid mass/molar density

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
  MultiFluidBase::PhaseProp::SliceType const phaseFrac,
  MultiFluidBase::PhaseProp::SliceType const phaseDens,
  MultiFluidBase::PhaseProp::SliceType const phaseMassDensity,
  MultiFluidBase::PhaseProp::SliceType const phaseVisc,
  MultiFluidBase::PhaseProp::SliceType const phaseEnthalpy,
  MultiFluidBase::PhaseProp::SliceType const phaseInternalEnergy,
  MultiFluidBase::PhaseComp::SliceType const phaseCompFrac,
  MultiFluidBase::FluidProp::SliceType const totalDensity ) const
{
  GEOS_UNUSED_VAR( phaseEnthalpy, phaseInternalEnergy );

  using Deriv = multifluid::DerivativeOffset;

  integer constexpr maxNumComp = MultiFluidBase::MAX_NUM_COMPONENTS;
  integer constexpr maxNumPhase = MultiFluidBase::MAX_NUM_PHASES;
  integer const numComp = numComponents();
  integer const numPhase = numPhases();

  // 1. Convert input mass fractions to mole fractions and keep derivatives

  stackArray1d< real64, maxNumComp > compMoleFrac( numComp );
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

  // 2. Compute phase fractions and phase component fractions
  m_flash.compute( pressure,
                   temperature,
                   compMoleFrac.toSliceConst(),
                   phaseFrac,
                   phaseCompFrac );

  // 3. Calculate the phase densities
  m_phase1.density.compute( pressure,
                            temperature,
                            phaseCompFrac.value[0].toSliceConst(),
                            phaseCompFrac.derivs[0].toSliceConst(),
                            phaseDens.value[0],
                            phaseDens.derivs[0],
                            phaseMassDensity.value[0],
                            phaseMassDensity.derivs[0],
                            m_useMass );
  m_phase2.density.compute( pressure,
                            temperature,
                            phaseCompFrac.value[1].toSliceConst(),
                            phaseCompFrac.derivs[1].toSliceConst(),
                            phaseDens.value[1],
                            phaseDens.derivs[1],
                            phaseMassDensity.value[1],
                            phaseMassDensity.derivs[1],
                            m_useMass );

  // 4. Calculate the phase viscosities
  m_phase1.viscosity.compute( pressure,
                              temperature,
                              phaseCompFrac.value[0].toSliceConst(),
                              phaseCompFrac.derivs[0].toSliceConst(),
                              phaseMassDensity.value[0],
                              phaseMassDensity.derivs[0].toSliceConst(),
                              phaseVisc.value[0],
                              phaseVisc.derivs[0],
                              m_useMass );
  m_phase2.viscosity.compute( pressure,
                              temperature,
                              phaseCompFrac.value[1].toSliceConst(),
                              phaseCompFrac.derivs[1].toSliceConst(),
                              phaseMassDensity.value[1],
                              phaseMassDensity.derivs[1].toSliceConst(),
                              phaseVisc.value[1],
                              phaseVisc.derivs[1],
                              m_useMass );

  // 5. if mass variables used instead of molar, perform the conversion
  if( m_useMass )
  {
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
                            phaseFrac,
                            phaseCompFrac,
                            phaseDens.derivs,
                            phaseVisc.derivs,
                            phaseEnthalpy.derivs,
                            phaseInternalEnergy.derivs );
  }

  // 6. Compute total fluid mass/molar density and derivatives

  computeTotalDensity( phaseFrac,
                       phaseDens,
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
