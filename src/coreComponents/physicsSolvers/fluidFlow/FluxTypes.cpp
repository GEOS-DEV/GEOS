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
 * @file FluxTypes.cpp
 * @brief Implementation of flux computation interfaces
 */

#include "FluxTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"

namespace geos
{

namespace fluxTypes
{

// ========================================
// TotalMassFlux implementations
// ========================================

template< int NUM_PHASES >
GEOS_HOST_DEVICE
real64 TotalMassFlux< NUM_PHASES >::computeTotalMassFlux( real64 const pressure,
                                                          arraySlice1d< real64 const > const & densities,
                                                          arraySlice1d< real64 const > const & saturations,
                                                          arraySlice1d< real64 const > const & mobilities,
                                                          real64 const gravityTerm )
{
  real64 totalFlux = 0.0;
  
  for( int phaseIndex = 0; phaseIndex < NUM_PHASES; ++phaseIndex )
  {
    real64 const phaseFlux = densities[phaseIndex] * saturations[phaseIndex] * mobilities[phaseIndex] * 
                            (pressure + densities[phaseIndex] * gravityTerm);
    totalFlux += phaseFlux;
  }
  
  return totalFlux;
}

template< int NUM_PHASES >
GEOS_HOST_DEVICE
void TotalMassFlux< NUM_PHASES >::computeTotalMassFluxDerivatives( real64 const pressure,
                                                                   arraySlice1d< real64 const > const & densities,
                                                                   arraySlice1d< real64 const > const & saturations,
                                                                   arraySlice1d< real64 const > const & mobilities,
                                                                   real64 & dFlux_dPressure,
                                                                   arraySlice1d< real64 > const & dFlux_dSaturations )
{
  dFlux_dPressure = 0.0;
  
  for( int phaseIndex = 0; phaseIndex < NUM_PHASES; ++phaseIndex )
  {
    // Derivative w.r.t. pressure
    dFlux_dPressure += densities[phaseIndex] * saturations[phaseIndex] * mobilities[phaseIndex];
    
    // Derivatives w.r.t. saturations
    dFlux_dSaturations[phaseIndex] = densities[phaseIndex] * mobilities[phaseIndex] * 
                                    (pressure + densities[phaseIndex] * gravityTerm);
  }
}

// Specialization for single-phase flow
template<>
GEOS_HOST_DEVICE
real64 TotalMassFlux< 1 >::computeTotalMassFlux( real64 const pressure,
                                                 real64 const density,
                                                 real64 const mobility,
                                                 real64 const gravityTerm )
{
  return density * mobility * (pressure + density * gravityTerm);
}

template<>
GEOS_HOST_DEVICE
void TotalMassFlux< 1 >::computeTotalMassFluxDerivatives( real64 const pressure,
                                                          real64 const density,
                                                          real64 const mobility,
                                                          real64 & dFlux_dPressure )
{
  GEOS_UNUSED_VAR( pressure );
  dFlux_dPressure = density * mobility;
}

// Specialization for two-phase flow
template<>
GEOS_HOST_DEVICE
real64 TotalMassFlux< 2 >::computeTotalMassFlux( real64 const pressure,
                                                 real64 const density1,
                                                 real64 const density2,
                                                 real64 const saturation1,
                                                 real64 const saturation2,
                                                 real64 const mobility1,
                                                 real64 const mobility2,
                                                 real64 const gravityTerm )
{
  real64 const flux1 = density1 * saturation1 * mobility1 * (pressure + density1 * gravityTerm);
  real64 const flux2 = density2 * saturation2 * mobility2 * (pressure + density2 * gravityTerm);
  return flux1 + flux2;
}

// ========================================
// TotalCapillaryFlux implementations
// ========================================

template< int NUM_PHASES >
GEOS_HOST_DEVICE
real64 TotalCapillaryFlux< NUM_PHASES >::computeTotalCapillaryFlux( arraySlice1d< real64 const > const & capillaryPressures,
                                                                    arraySlice1d< real64 const > const & densities,
                                                                    arraySlice1d< real64 const > const & saturations,
                                                                    arraySlice1d< real64 const > const & capillaryMobilities )
{
  real64 totalCapillaryFlux = 0.0;
  
  for( int phaseIndex = 0; phaseIndex < NUM_PHASES - 1; ++phaseIndex )
  {
    real64 const capillaryFlux = densities[phaseIndex] * saturations[phaseIndex] * capillaryMobilities[phaseIndex] * 
                                capillaryPressures[phaseIndex];
    totalCapillaryFlux += capillaryFlux;
  }
  
  return totalCapillaryFlux;
}

// ========================================
// PhaseFlux implementations
// ========================================

template< int NUM_PHASES, int PHASE_INDEX >
GEOS_HOST_DEVICE
real64 PhaseFlux< NUM_PHASES, PHASE_INDEX >::computeMassicFlux( real64 const pressure,
                                                                real64 const density,
                                                                real64 const saturation,
                                                                real64 const mobility,
                                                                real64 const gravityTerm )
{
  return density * saturation * mobility * (pressure + density * gravityTerm);
}

template< int NUM_PHASES, int PHASE_INDEX >
GEOS_HOST_DEVICE
real64 PhaseFlux< NUM_PHASES, PHASE_INDEX >::computeCapillaryFlux( real64 const capillaryPressure,
                                                                   real64 const density,
                                                                   real64 const saturation,
                                                                   real64 const capillaryMobility )
{
  return density * saturation * capillaryMobility * capillaryPressure;
}

template< int NUM_PHASES, int PHASE_INDEX >
GEOS_HOST_DEVICE
real64 PhaseFlux< NUM_PHASES, PHASE_INDEX >::computeBuoyancyFlux( real64 const density,
                                                                  real64 const saturation,
                                                                  real64 const mobility,
                                                                  arraySlice1d< real64 const > const & gravityVector,
                                                                  arraySlice2d< real64 const > const & permeability )
{
  // Compute gravity flux contribution
  real64 gravityFlux = 0.0;
  
  // This is simplified - actual implementation would involve proper tensor operations
  for( int i = 0; i < 3; ++i )
  {
    for( int j = 0; j < 3; ++j )
    {
      gravityFlux += permeability[i][j] * gravityVector[j] * gravityVector[i];
    }
  }
  
  return density * density * saturation * mobility * gravityFlux;
}

// Explicit template instantiations for common use cases
template struct TotalMassFlux< 1 >;
template struct TotalMassFlux< 2 >;
template struct TotalMassFlux< 3 >;

template struct TotalCapillaryFlux< 2 >;
template struct TotalCapillaryFlux< 3 >;

template struct PhaseFlux< 2, 0 >;
template struct PhaseFlux< 2, 1 >;
template struct PhaseFlux< 3, 0 >;
template struct PhaseFlux< 3, 1 >;
template struct PhaseFlux< 3, 2 >;

} // namespace fluxTypes

} // namespace geos