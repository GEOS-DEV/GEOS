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
 * @file CompositionalMultiphaseFluidUpdates.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_COMPOSITIONALMULTIPHASEFLUIDUPDATES_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_COMPOSITIONALMULTIPHASEFLUIDUPDATES_HPP_

#include "constitutive/fluid/multifluid/compositional/parameters/ComponentProperties.hpp"
#include "constitutive/fluid/multifluid/compositional/models/NullModel.hpp"

#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
#include "constitutive/fluid/multifluid/MultiFluidUtils.hpp"

namespace geos
{
namespace constitutive
{

/**
 * @brief Kernel wrapper class for CompositionalMultiphaseFluid.
 * @tparam FLASH Class describing the phase equilibrium model
 * @tparam PHASES Classes describing the phase property models for each phase.
 */
template< typename FLASH, typename ... PHASES >
class CompositionalMultiphaseFluidUpdates final : public MultiFluidBase::KernelWrapper
{
public:
  using FlashModel = FLASH;

  // Get the number of phases
  static constexpr integer NUM_PHASES = static_cast< integer >(sizeof...(PHASES));

  // Ensure that the number of phases matches the flash object
  static_assert( NUM_PHASES == FlashModel::KernelWrapper::getNumberOfPhases() );

public:
  CompositionalMultiphaseFluidUpdates( compositional::ComponentProperties const & componentProperties,
                                       FLASH const & flash,
                                       PHASES const & ... phases,
                                       arrayView1d< integer const > const & phaseOrder,
                                       arrayView1d< real64 const > const & componentMolarWeight,
                                       bool const useMass,
                                       MultiFluidBase::PhaseProp::ViewType phaseFraction,
                                       MultiFluidBase::PhaseProp::ViewType phaseDensity,
                                       MultiFluidBase::PhaseProp::ViewType phaseMassDensity,
                                       MultiFluidBase::PhaseProp::ViewType phaseViscosity,
                                       MultiFluidBase::PhaseProp::ViewType phaseEnthalpy,
                                       MultiFluidBase::PhaseProp::ViewType phaseInternalEnergy,
                                       MultiFluidBase::PhaseComp::ViewType phaseComponentFraction,
                                       MultiFluidBase::FluidProp::ViewType totalDensity,
                                       MultiFluidBase::PhaseComp::ViewValueType kValues );

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
                        MultiFluidBase::PhaseComp::SliceType const phaseComponentFraction,
                        MultiFluidBase::FluidProp::SliceType const totalDensity ) const override;

  GEOS_HOST_DEVICE
  virtual void update( localIndex const k,
                       localIndex const q,
                       real64 const pressure,
                       real64 const temperature,
                       arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition ) const override;

protected:
  GEOS_HOST_DEVICE
  void compute( real64 const pressure,
                real64 const temperature,
                arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition,
                MultiFluidBase::PhaseProp::SliceType const phaseFraction,
                MultiFluidBase::PhaseProp::SliceType const phaseDensity,
                MultiFluidBase::PhaseProp::SliceType const phaseMassDensity,
                MultiFluidBase::PhaseProp::SliceType const phaseViscosity,
                MultiFluidBase::PhaseProp::SliceType const phaseEnthalpy,
                MultiFluidBase::PhaseProp::SliceType const phaseInternalEnergy,
                MultiFluidBase::PhaseComp::SliceType const phaseComponentFraction,
                MultiFluidBase::FluidProp::SliceType const totalDensity,
                MultiFluidBase::PhaseComp::SliceType::ValueType const & kValues ) const;

  /**
   * @brief Helper to unpack phase property models and compute phase densities and viscosities.
   */
  template< std::size_t... Is >
  GEOS_HOST_DEVICE
  void computePhaseProperties( real64 const pressure,
                               real64 const temperature,
                               MultiFluidBase::PhaseComp::SliceType const phaseComponentFraction,
                               MultiFluidBase::PhaseProp::SliceType const phaseDensity,
                               MultiFluidBase::PhaseProp::SliceType const phaseMassDensity,
                               MultiFluidBase::PhaseProp::SliceType const phaseViscosity,
                               MultiFluidBase::PhaseProp::SliceType const phaseEnthalpy,
                               std::index_sequence< Is... > ) const;

  /**
   * @brief Convert derivatives from phase mole fraction to total mole fraction
   * @details Given property derivatives @c dProperty where composition derivatives are with
   *          respect to a phase compositions, this will transform that properties so that
   *          they the composition derivatives are with respect to total composition. The derivatives
   *          of the phase composition should be provided in @c dPhaseComposition.
   * @param[in] numComps The number of components
   * @param[in] dPhaseComposition Derivatives of the phase composition
   * @param[in,out] dProperty The derivatives of the property
   * @param[in] workSpace Temporary workspace
   */
  template< int USD1, int USD2 >
  GEOS_HOST_DEVICE
  static void convertDerivativesToTotalMoleFraction( integer const numComps,
                                                     arraySlice2d< real64 const, USD1 > const & dPhaseComposition,
                                                     arraySlice1d< real64, USD2 > const & dProperty,
                                                     arraySlice1d< real64 > const & workSpace );

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void setZero( real64 & val ){ val = 0.0; }

  GEOS_HOST_DEVICE
  static constexpr bool isThermalType()
  {
    return (!std::is_same_v< typename PHASES::Enthalpy, compositional::NullModel > && ...);
  }

private:
  // The component properties
  compositional::ComponentProperties::KernelWrapper m_componentProperties;

  // Flash kernel wrapper
  typename FLASH::KernelWrapper m_flash;

  // The ordering of phases
  arrayView1d< integer const > const m_phaseOrder;

  // Phase model kernel wrappers stored in a tuple
  camp::tuple< typename PHASES::KernelWrapper... > m_phases;

  // Backup variables
  MultiFluidBase::PhaseComp::ViewValueType m_kValues;
};

template< typename FLASH, typename ... PHASES >
CompositionalMultiphaseFluidUpdates< FLASH, PHASES... >::
CompositionalMultiphaseFluidUpdates( compositional::ComponentProperties const & componentProperties,
                                     FLASH const & flash,
                                     PHASES const &... phases,
                                     arrayView1d< integer const > const & phaseOrder,
                                     arrayView1d< real64 const > const & componentMolarWeight,
                                     bool const useMass,
                                     MultiFluidBase::PhaseProp::ViewType phaseFraction,
                                     MultiFluidBase::PhaseProp::ViewType phaseDensity,
                                     MultiFluidBase::PhaseProp::ViewType phaseMassDensity,
                                     MultiFluidBase::PhaseProp::ViewType phaseViscosity,
                                     MultiFluidBase::PhaseProp::ViewType phaseEnthalpy,
                                     MultiFluidBase::PhaseProp::ViewType phaseInternalEnergy,
                                     MultiFluidBase::PhaseComp::ViewType phaseComponentFraction,
                                     MultiFluidBase::FluidProp::ViewType totalDensity,
                                     MultiFluidBase::PhaseComp::ViewValueType kValues ):
  MultiFluidBase::KernelWrapper( componentMolarWeight,
                                 useMass,
                                 std::move( phaseFraction ),
                                 std::move( phaseDensity ),
                                 std::move( phaseMassDensity ),
                                 std::move( phaseViscosity ),
                                 std::move( phaseEnthalpy ),
                                 std::move( phaseInternalEnergy ),
                                 std::move( phaseComponentFraction ),
                                 std::move( totalDensity ) ),
  m_componentProperties( componentProperties.createKernelWrapper() ),
  m_flash( flash.createKernelWrapper() ),
  m_phaseOrder( phaseOrder ),
  m_phases( phases.createKernelWrapper()... ),
  m_kValues( kValues )
{}

template< typename FLASH, typename ... PHASES >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void
CompositionalMultiphaseFluidUpdates< FLASH, PHASES... >::compute(
  real64 const pressure,
  real64 const temperature,
  arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition,
  MultiFluidBase::PhaseProp::SliceType const phaseFraction,
  MultiFluidBase::PhaseProp::SliceType const phaseDensity,
  MultiFluidBase::PhaseProp::SliceType const phaseMassDensity,
  MultiFluidBase::PhaseProp::SliceType const phaseViscosity,
  MultiFluidBase::PhaseProp::SliceType const phaseEnthalpy,
  MultiFluidBase::PhaseProp::SliceType const phaseInternalEnergy,
  MultiFluidBase::PhaseComp::SliceType const phaseComponentFraction,
  MultiFluidBase::FluidProp::SliceType const totalDensity ) const
{
  integer constexpr maxNumComp = MultiFluidBase::MAX_NUM_COMPONENTS;
  integer constexpr maxNumPhase = MultiFluidBase::MAX_NUM_PHASES - 1;
  MultiFluidBase::PhaseComp::StackValueType< maxNumPhase *maxNumComp > kValues( 1, 1, numPhases() - 1, numComponents() );

  LvArray::forValuesInSlice( kValues[0][0], setZero );   // Force initialisation of k-Values

  compute( pressure,
           temperature,
           composition,
           phaseFraction,
           phaseDensity,
           phaseMassDensity,
           phaseViscosity,
           phaseEnthalpy,
           phaseInternalEnergy,
           phaseComponentFraction,
           totalDensity,
           kValues[0][0] );
}

template< typename FLASH, typename ... PHASES >
template< std::size_t... Is >
GEOS_HOST_DEVICE
void CompositionalMultiphaseFluidUpdates< FLASH, PHASES... >::computePhaseProperties(
  real64 const pressure,
  real64 const temperature,
  MultiFluidBase::PhaseComp::SliceType const phaseComponentFraction,
  MultiFluidBase::PhaseProp::SliceType const phaseDensity,
  MultiFluidBase::PhaseProp::SliceType const phaseMassDensity,
  MultiFluidBase::PhaseProp::SliceType const phaseViscosity,
  MultiFluidBase::PhaseProp::SliceType const phaseEnthalpy,
  std::index_sequence< Is... > ) const
{
  // Density computations
  ( camp::get< Is >( m_phases ).density.compute(
      m_componentProperties,
      pressure,
      temperature,
      phaseComponentFraction.value[m_phaseOrder[Is]].toSliceConst(),
      phaseDensity.value[m_phaseOrder[Is]],
      phaseDensity.derivs[m_phaseOrder[Is]],
      phaseMassDensity.value[m_phaseOrder[Is]],
      phaseMassDensity.derivs[m_phaseOrder[Is]],
      m_useMass ), ... );

  // Viscosity computations
  ( camp::get< Is >( m_phases ).viscosity.compute(
      m_componentProperties,
      pressure,
      temperature,
      phaseComponentFraction.value[m_phaseOrder[Is]].toSliceConst(),
      phaseMassDensity.value[m_phaseOrder[Is]],
      phaseMassDensity.derivs[m_phaseOrder[Is]].toSliceConst(),
      phaseViscosity.value[m_phaseOrder[Is]],
      phaseViscosity.derivs[m_phaseOrder[Is]],
      m_useMass ), ... );

  if constexpr (isThermalType())
  {
    // Enthalpy computations
    ( camp::get< Is >( m_phases ).enthalpy.compute(
        m_componentProperties,
        pressure,
        temperature,
        phaseComponentFraction.value[m_phaseOrder[Is]].toSliceConst(),
        phaseEnthalpy.value[m_phaseOrder[Is]],
        phaseEnthalpy.derivs[m_phaseOrder[Is]],
        m_useMass ), ... );
  }
}

template< typename FLASH, typename ... PHASES >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void
CompositionalMultiphaseFluidUpdates< FLASH, PHASES... >::compute(
  real64 const pressure,
  real64 const temperature,
  arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition,
  MultiFluidBase::PhaseProp::SliceType const phaseFraction,
  MultiFluidBase::PhaseProp::SliceType const phaseDensity,
  MultiFluidBase::PhaseProp::SliceType const phaseMassDensity,
  MultiFluidBase::PhaseProp::SliceType const phaseViscosity,
  MultiFluidBase::PhaseProp::SliceType const phaseEnthalpy,
  MultiFluidBase::PhaseProp::SliceType const phaseInternalEnergy,
  MultiFluidBase::PhaseComp::SliceType const phaseComponentFraction,
  MultiFluidBase::FluidProp::SliceType const totalDensity,
  MultiFluidBase::PhaseComp::SliceType::ValueType const & kValues ) const
{
  integer constexpr maxNumComp = MultiFluidBase::MAX_NUM_COMPONENTS;
  integer constexpr maxNumDof = MultiFluidBase::MAX_NUM_COMPONENTS + 2;
  integer const numComp = numComponents();
  integer const numDof = numComp + 2;

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
  m_flash.compute( m_componentProperties,
                   pressure,
                   temperature,
                   compMoleFrac.toSliceConst(),
                   kValues,
                   phaseFraction,
                   phaseComponentFraction );

  // 3. Calculate phase properties: density, viscosity and enthalpy
  computePhaseProperties( pressure,
                          temperature,
                          phaseComponentFraction,
                          phaseDensity,
                          phaseMassDensity,
                          phaseViscosity,
                          phaseEnthalpy,
                          std::make_index_sequence< NUM_PHASES >{} );

  // 5. Convert derivatives from phase composition to total composition
  stackArray1d< real64, maxNumDof > workSpace( numDof );
  for( integer ip = 0; ip < NUM_PHASES; ++ip )
  {
    convertDerivativesToTotalMoleFraction( numComp,
                                           phaseComponentFraction.derivs[ip].toSliceConst(),
                                           phaseDensity.derivs[ip],
                                           workSpace );
    convertDerivativesToTotalMoleFraction( numComp,
                                           phaseComponentFraction.derivs[ip].toSliceConst(),
                                           phaseMassDensity.derivs[ip],
                                           workSpace );
    convertDerivativesToTotalMoleFraction( numComp,
                                           phaseComponentFraction.derivs[ip].toSliceConst(),
                                           phaseViscosity.derivs[ip],
                                           workSpace );
    if constexpr (isThermalType())
    {
      convertDerivativesToTotalMoleFraction( numComp,
                                             phaseComponentFraction.derivs[ip].toSliceConst(),
                                             phaseEnthalpy.derivs[ip],
                                             workSpace );
    }
  }

  // 4. Calculate the internal energy
  if constexpr (isThermalType())
  {
    computeInternalEnergy( pressure,
                           phaseFraction,
                           phaseMassDensity,
                           phaseEnthalpy,
                           phaseInternalEnergy );
  }

  // 6. if mass variables used instead of molar, perform the conversion
  if( m_useMass )
  {
    real64 phaseMolecularWeight[NUM_PHASES]{};
    real64 dPhaseMolecularWeight[NUM_PHASES][maxNumDof]{};

    arrayView1d< real64 const > const & componentMolarWeight = m_componentProperties.m_componentMolarWeight;

    for( integer ip = 0; ip < NUM_PHASES; ++ip )
    {
      auto const & phaseComposition = phaseComponentFraction.value[ip].toSliceConst();
      auto const & dPhaseComposition = phaseComponentFraction.derivs[ip].toSliceConst();

      for( integer ic = 0; ic < numComp; ++ic )
      {
        phaseMolecularWeight[ip] += phaseComposition[ic] * componentMolarWeight[ic];
        for( integer kc = 0; kc < numDof; ++kc )
        {
          dPhaseMolecularWeight[ip][kc] += dPhaseComposition( ic, kc ) * componentMolarWeight[ic];
        }
      }
    }

    convertToMassFractions( dCompMoleFrac_dCompMassFrac,
                            phaseMolecularWeight,
                            dPhaseMolecularWeight,
                            phaseFraction,
                            phaseComponentFraction,
                            phaseMassDensity.derivs,
                            phaseViscosity.derivs,
                            phaseEnthalpy.derivs,
                            phaseInternalEnergy.derivs );

    // Molar density equals mass density
    for( integer ip = 0; ip < NUM_PHASES; ++ip )
    {
      phaseDensity.value[ip] = phaseMassDensity.value[ip];
      for( integer idof = 0; idof < numDof; ++idof )
      {
        phaseDensity.derivs( ip, idof ) = phaseMassDensity.derivs( ip, idof );
      }
    }
  }

  // 7. Compute total fluid mass/molar density and derivatives

  computeTotalDensity( phaseFraction,
                       phaseDensity,
                       totalDensity );
}

template< typename FLASH, typename ... PHASES >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void
CompositionalMultiphaseFluidUpdates< FLASH, PHASES... >::
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
           m_totalDensity( k, q ),
           m_kValues[k][q] );
}

template< typename FLASH, typename ... PHASES >
template< int USD1, int USD2 >
GEOS_HOST_DEVICE
void
CompositionalMultiphaseFluidUpdates< FLASH, PHASES... >::
convertDerivativesToTotalMoleFraction( integer const numComps,
                                       arraySlice2d< real64 const, USD1 > const & dPhaseComposition,
                                       arraySlice1d< real64, USD2 > const & dProperty,
                                       arraySlice1d< real64 > const & workSpace )
{
  using Deriv = constitutive::multifluid::DerivativeOffset;
  integer const numDofs = numComps + 2;
  for( integer kc = 0; kc < numDofs; ++kc )
  {
    workSpace[kc] = dProperty[kc];
  }
  for( integer ic = 0; ic < numComps; ++ic )
  {
    dProperty[Deriv::dC+ic] = 0.0;
  }
  for( integer kc = 0; kc < numDofs; ++kc )
  {
    for( integer ic = 0; ic < numComps; ++ic )
    {
      dProperty[kc] += (dPhaseComposition( ic, kc ) * workSpace[Deriv::dC+ic]);
    }
  }
}

} /* namespace constitutive */

} /* namespace geos */

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_COMPOSITIONALMULTIPHASEFLUIDUPDATES_HPP_
