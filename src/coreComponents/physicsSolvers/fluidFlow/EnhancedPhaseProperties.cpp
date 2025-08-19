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
 * @file EnhancedPhaseProperties.cpp
 * @brief Implementation of enhanced phase model interface with fractional flows and mobility functions
 */

#include "EnhancedPhaseProperties.hpp"
#include "common/GEOS_RAJA_Interface.hpp"

namespace geos
{

namespace phaseProperties
{

// ========================================
// EnhancedPhaseModel implementations
// ========================================

template< int NUM_PHASES >
GEOS_HOST_DEVICE
real64 EnhancedPhaseModel< NUM_PHASES >::computeFractionalFlow( integer const phaseIndex,
                                                                arraySlice1d< real64 const > const & mobilities )
{
  real64 totalMobility = computeTotalMobility( mobilities );
  
  // Avoid division by zero
  if( totalMobility < 1e-20 )
  {
    return 0.0;
  }
  
  return mobilities[phaseIndex] / totalMobility;
}

template< int NUM_PHASES >
GEOS_HOST_DEVICE
real64 EnhancedPhaseModel< NUM_PHASES >::computeTotalMobility( arraySlice1d< real64 const > const & mobilities )
{
  real64 totalMobility = 0.0;
  
  for( int phaseIndex = 0; phaseIndex < NUM_PHASES; ++phaseIndex )
  {
    totalMobility += mobilities[phaseIndex];
  }
  
  return totalMobility;
}

template< int NUM_PHASES >
GEOS_HOST_DEVICE
void EnhancedPhaseModel< NUM_PHASES >::computeFractionalFlowDerivatives( integer const phaseIndex,
                                                                         arraySlice1d< real64 const > const & mobilities,
                                                                         arraySlice2d< real64 const > const & dMobilities_dSaturations,
                                                                         arraySlice1d< real64 > const & dFractionalFlow_dSaturations )
{
  real64 const totalMobility = computeTotalMobility( mobilities );
  real64 const mobility_i = mobilities[phaseIndex];
  
  // Avoid division by zero
  if( totalMobility < 1e-20 )
  {
    for( int j = 0; j < NUM_PHASES; ++j )
    {
      dFractionalFlow_dSaturations[j] = 0.0;
    }
    return;
  }
  
  for( int j = 0; j < NUM_PHASES; ++j )
  {
    // Compute derivative of total mobility w.r.t. saturation j
    real64 dTotalMobility_dSj = 0.0;
    for( int k = 0; k < NUM_PHASES; ++k )
    {
      dTotalMobility_dSj += dMobilities_dSaturations[k][j];
    }
    
    // Apply quotient rule: d(f/g)/dx = (g*df/dx - f*dg/dx)/g^2
    real64 const numerator = totalMobility * dMobilities_dSaturations[phaseIndex][j] - 
                            mobility_i * dTotalMobility_dSj;
    dFractionalFlow_dSaturations[j] = numerator / (totalMobility * totalMobility);
  }
}

template< int NUM_PHASES >
GEOS_HOST_DEVICE
void EnhancedPhaseModel< NUM_PHASES >::computeTotalMobilityDerivatives( arraySlice2d< real64 const > const & dMobilities_dSaturations,
                                                                        arraySlice1d< real64 > const & dTotalMobility_dSaturations )
{
  for( int j = 0; j < NUM_PHASES; ++j )
  {
    dTotalMobility_dSaturations[j] = 0.0;
    for( int i = 0; i < NUM_PHASES; ++i )
    {
      dTotalMobility_dSaturations[j] += dMobilities_dSaturations[i][j];
    }
  }
}

// Specialization for two-phase systems
template<>
GEOS_HOST_DEVICE
real64 EnhancedPhaseModel< 2 >::computeFractionalFlow( integer const phaseIndex,
                                                       real64 const mobility1,
                                                       real64 const mobility2 )
{
  real64 const totalMobility = mobility1 + mobility2;
  
  if( totalMobility < 1e-20 )
  {
    return 0.0;
  }
  
  if( phaseIndex == 0 )
  {
    return mobility1 / totalMobility;
  }
  else
  {
    return mobility2 / totalMobility;
  }
}

template<>
GEOS_HOST_DEVICE
void EnhancedPhaseModel< 2 >::computeFractionalFlowDerivatives( integer const phaseIndex,
                                                                real64 const mobility1,
                                                                real64 const mobility2,
                                                                real64 const dMobility1_dS1,
                                                                real64 const dMobility2_dS1,
                                                                real64 & dFractionalFlow_dS1 )
{
  real64 const totalMobility = mobility1 + mobility2;
  
  if( totalMobility < 1e-20 )
  {
    dFractionalFlow_dS1 = 0.0;
    return;
  }
  
  real64 const dTotalMobility_dS1 = dMobility1_dS1 + dMobility2_dS1;
  
  if( phaseIndex == 0 )
  {
    // Phase 1 fractional flow derivative
    real64 const numerator = totalMobility * dMobility1_dS1 - mobility1 * dTotalMobility_dS1;
    dFractionalFlow_dS1 = numerator / (totalMobility * totalMobility);
  }
  else
  {
    // Phase 2 fractional flow derivative
    real64 const numerator = totalMobility * dMobility2_dS1 - mobility2 * dTotalMobility_dS1;
    dFractionalFlow_dS1 = numerator / (totalMobility * totalMobility);
  }
}

// ========================================
// MixtureProperties implementations
// ========================================

template< int NUM_PHASES, int NUM_COMPONENTS >
GEOS_HOST_DEVICE
real64 MixtureProperties< NUM_PHASES, NUM_COMPONENTS >::computeMixtureDensity( arraySlice1d< real64 const > const & phaseDensities,
                                                                               arraySlice1d< real64 const > const & saturations )
{
  real64 mixtureDensity = 0.0;
  
  for( int phaseIndex = 0; phaseIndex < NUM_PHASES; ++phaseIndex )
  {
    mixtureDensity += phaseDensities[phaseIndex] * saturations[phaseIndex];
  }
  
  return mixtureDensity;
}

template< int NUM_PHASES, int NUM_COMPONENTS >
GEOS_HOST_DEVICE
real64 MixtureProperties< NUM_PHASES, NUM_COMPONENTS >::computeMixtureEnthalpy( arraySlice1d< real64 const > const & phaseEnthalpies,
                                                                                arraySlice1d< real64 const > const & saturations,
                                                                                arraySlice1d< real64 const > const & phaseDensities )
{
  real64 mixtureEnthalpy = 0.0;
  real64 totalMass = 0.0;
  
  for( int phaseIndex = 0; phaseIndex < NUM_PHASES; ++phaseIndex )
  {
    real64 const phaseMass = phaseDensities[phaseIndex] * saturations[phaseIndex];
    mixtureEnthalpy += phaseEnthalpies[phaseIndex] * phaseMass;
    totalMass += phaseMass;
  }
  
  // Mass-weighted enthalpy
  return (totalMass > 1e-20) ? mixtureEnthalpy / totalMass : 0.0;
}

template< int NUM_PHASES, int NUM_COMPONENTS >
GEOS_HOST_DEVICE
real64 MixtureProperties< NUM_PHASES, NUM_COMPONENTS >::computeMixtureViscosity( arraySlice1d< real64 const > const & phaseViscosities,
                                                                                 arraySlice1d< real64 const > const & saturations )
{
  // Harmonic average for viscosity
  real64 inverseMixtureViscosity = 0.0;
  real64 totalSaturation = 0.0;
  
  for( int phaseIndex = 0; phaseIndex < NUM_PHASES; ++phaseIndex )
  {
    if( phaseViscosities[phaseIndex] > 1e-20 && saturations[phaseIndex] > 1e-20 )
    {
      inverseMixtureViscosity += saturations[phaseIndex] / phaseViscosities[phaseIndex];
      totalSaturation += saturations[phaseIndex];
    }
  }
  
  return (inverseMixtureViscosity > 1e-20) ? totalSaturation / inverseMixtureViscosity : 0.0;
}

template< int NUM_PHASES, int NUM_COMPONENTS >
GEOS_HOST_DEVICE
void MixtureProperties< NUM_PHASES, NUM_COMPONENTS >::computeMixtureDensityDerivatives( arraySlice1d< real64 const > const & phaseDensities,
                                                                                        arraySlice1d< real64 const > const & saturations,
                                                                                        arraySlice1d< real64 const > const & dDensities_dPressure,
                                                                                        arraySlice2d< real64 const > const & dDensities_dSaturations,
                                                                                        real64 & dMixtureDensity_dPressure,
                                                                                        arraySlice1d< real64 > const & dMixtureDensity_dSaturations )
{
  // Derivative w.r.t. pressure
  dMixtureDensity_dPressure = 0.0;
  for( int phaseIndex = 0; phaseIndex < NUM_PHASES; ++phaseIndex )
  {
    dMixtureDensity_dPressure += dDensities_dPressure[phaseIndex] * saturations[phaseIndex];
  }
  
  // Derivatives w.r.t. saturations
  for( int j = 0; j < NUM_PHASES; ++j )
  {
    dMixtureDensity_dSaturations[j] = phaseDensities[j]; // Direct contribution
    
    // Add contributions from density variations with saturation
    for( int i = 0; i < NUM_PHASES; ++i )
    {
      dMixtureDensity_dSaturations[j] += dDensities_dSaturations[i][j] * saturations[i];
    }
  }
}

// ========================================
// CapillaryPressureUtils implementations
// ========================================

template< int NUM_PHASES >
GEOS_HOST_DEVICE
real64 CapillaryPressureUtils< NUM_PHASES >::computeCapillaryPressure( integer const phase1Index,
                                                                       integer const phase2Index,
                                                                       arraySlice1d< real64 const > const & saturations )
{
  // This is a simplified implementation - actual capillary pressure models would be more complex
  // For example, Brooks-Corey or van Genuchten models
  
  GEOS_UNUSED_VAR( phase2Index );
  
  real64 const saturation = saturations[phase1Index];
  
  // Simple linear capillary pressure model for demonstration
  // Pc = Pc_entry * (1 - S_eff)
  real64 const Pc_entry = 1000.0; // Entry pressure (Pa)
  real64 const S_eff = LvArray::math::max( 0.0, LvArray::math::min( 1.0, saturation ) );
  
  return Pc_entry * (1.0 - S_eff);
}

template< int NUM_PHASES >
GEOS_HOST_DEVICE
void CapillaryPressureUtils< NUM_PHASES >::computeCapillaryPressureDerivatives( integer const phase1Index,
                                                                                integer const phase2Index,
                                                                                arraySlice1d< real64 const > const & saturations,
                                                                                arraySlice1d< real64 > const & dPc_dSaturations )
{
  GEOS_UNUSED_VAR( phase2Index );
  GEOS_UNUSED_VAR( saturations );
  
  // Initialize all derivatives to zero
  for( int i = 0; i < NUM_PHASES; ++i )
  {
    dPc_dSaturations[i] = 0.0;
  }
  
  // Simple linear model derivative
  real64 const Pc_entry = 1000.0;
  dPc_dSaturations[phase1Index] = -Pc_entry;
}

// Explicit template instantiations for common use cases
template struct EnhancedPhaseModel< 1 >;
template struct EnhancedPhaseModel< 2 >;
template struct EnhancedPhaseModel< 3 >;

template struct MixtureProperties< 2, 2 >;
template struct MixtureProperties< 2, 3 >;
template struct MixtureProperties< 3, 3 >;
template struct MixtureProperties< 3, 4 >;

template struct CapillaryPressureUtils< 2 >;
template struct CapillaryPressureUtils< 3 >;

} // namespace phaseProperties

} // namespace geos