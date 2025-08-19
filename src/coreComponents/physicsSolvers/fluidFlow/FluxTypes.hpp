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
 * @file FluxTypes.hpp
 * @brief Comprehensive interface for defining flux types in flow solvers
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_FLUXTYPES_HPP_
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_FLUXTYPES_HPP_

#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"

namespace geos
{

namespace fluxTypes
{

/**
 * @brief Enumeration for flux computation types
 */
enum class FluxComputationType
{
  HYPERBOLIC_SCALAR,           ///< No diffusion/conduction
  PARABOLIC_SCALAR,            ///< With diffusion/conduction
  HYPERBOLIC_PARABOLIC_SCALAR  ///< Mixed hyperbolic-parabolic
};

/**
 * @brief Base template for flux computations
 * @tparam NUM_PHASES Number of phases
 * @tparam FLUX_TYPE Type of flux computation
 */
template< int NUM_PHASES, FluxComputationType FLUX_TYPE >
struct FluxInterface
{
  static constexpr int numPhases() { return NUM_PHASES; }
  static constexpr FluxComputationType fluxType() { return FLUX_TYPE; }
};

/**
 * @brief Interface for total mass flux computations
 */
template< int NUM_PHASES >
struct TotalMassFlux : FluxInterface< NUM_PHASES, FluxComputationType::PARABOLIC_SCALAR >
{
  /**
   * @brief Compute total mass flux across all phases
   * @param[in] pressure Pressure field
   * @param[in] densities Phase densities
   * @param[in] saturations Phase saturations
   * @param[in] mobilities Phase mobilities
   * @param[in] gravityTerm Gravity contribution
   * @return Total mass flux
   */
  GEOS_HOST_DEVICE
  static real64 computeTotalMassFlux( real64 const pressure,
                                      arraySlice1d< real64 const > const & densities,
                                      arraySlice1d< real64 const > const & saturations,
                                      arraySlice1d< real64 const > const & mobilities,
                                      real64 const gravityTerm );

  /**
   * @brief Compute derivatives of total mass flux
   * @param[in] pressure Pressure field
   * @param[in] densities Phase densities
   * @param[in] saturations Phase saturations
   * @param[in] mobilities Phase mobilities
   * @param[out] dFlux_dPressure Derivative w.r.t. pressure
   * @param[out] dFlux_dSaturations Derivatives w.r.t. saturations
   */
  GEOS_HOST_DEVICE
  static void computeTotalMassFluxDerivatives( real64 const pressure,
                                               arraySlice1d< real64 const > const & densities,
                                               arraySlice1d< real64 const > const & saturations,
                                               arraySlice1d< real64 const > const & mobilities,
                                               real64 & dFlux_dPressure,
                                               arraySlice1d< real64 > const & dFlux_dSaturations );
};

/**
 * @brief Interface for total capillary flux computations
 */
template< int NUM_PHASES >
struct TotalCapillaryFlux : FluxInterface< NUM_PHASES, FluxComputationType::PARABOLIC_SCALAR >
{
  /**
   * @brief Compute total capillary flux
   * @param[in] capillaryPressures Capillary pressures between phases
   * @param[in] densities Phase densities
   * @param[in] saturations Phase saturations
   * @param[in] capillaryMobilities Capillary mobilities
   * @return Total capillary flux
   */
  GEOS_HOST_DEVICE
  static real64 computeTotalCapillaryFlux( arraySlice1d< real64 const > const & capillaryPressures,
                                           arraySlice1d< real64 const > const & densities,
                                           arraySlice1d< real64 const > const & saturations,
                                           arraySlice1d< real64 const > const & capillaryMobilities );
};

/**
 * @brief Interface for phase-specific flux computations
 */
template< int NUM_PHASES, int PHASE_INDEX >
struct PhaseFlux : FluxInterface< NUM_PHASES, FluxComputationType::PARABOLIC_SCALAR >
{
  static_assert( PHASE_INDEX < NUM_PHASES, "Phase index must be less than number of phases" );

  /**
   * @brief Compute massic flux for a specific phase
   * @param[in] pressure Pressure field
   * @param[in] density Phase density
   * @param[in] saturation Phase saturation
   * @param[in] mobility Phase mobility
   * @param[in] gravityTerm Gravity contribution
   * @return Phase massic flux
   */
  GEOS_HOST_DEVICE
  static real64 computeMassicFlux( real64 const pressure,
                                   real64 const density,
                                   real64 const saturation,
                                   real64 const mobility,
                                   real64 const gravityTerm );

  /**
   * @brief Compute capillary flux for a specific phase
   * @param[in] capillaryPressure Capillary pressure
   * @param[in] density Phase density
   * @param[in] saturation Phase saturation
   * @param[in] capillaryMobility Capillary mobility
   * @return Phase capillary flux
   */
  GEOS_HOST_DEVICE
  static real64 computeCapillaryFlux( real64 const capillaryPressure,
                                      real64 const density,
                                      real64 const saturation,
                                      real64 const capillaryMobility );

  /**
   * @brief Compute buoyancy flux for a specific phase
   * @param[in] density Phase density
   * @param[in] saturation Phase saturation
   * @param[in] mobility Phase mobility
   * @param[in] gravityVector Gravity vector
   * @param[in] permeability Permeability tensor
   * @return Phase buoyancy flux
   */
  GEOS_HOST_DEVICE
  static real64 computeBuoyancyFlux( real64 const density,
                                     real64 const saturation,
                                     real64 const mobility,
                                     arraySlice1d< real64 const > const & gravityVector,
                                     arraySlice2d< real64 const > const & permeability );
};

/**
 * @brief Specialization for single-phase flow (NUM_PHASES = 1)
 */
template<>
struct TotalMassFlux< 1 > : FluxInterface< 1, FluxComputationType::PARABOLIC_SCALAR >
{
  GEOS_HOST_DEVICE
  static real64 computeTotalMassFlux( real64 const pressure,
                                      real64 const density,
                                      real64 const mobility,
                                      real64 const gravityTerm );

  GEOS_HOST_DEVICE
  static void computeTotalMassFluxDerivatives( real64 const pressure,
                                               real64 const density,
                                               real64 const mobility,
                                               real64 & dFlux_dPressure );
};

/**
 * @brief Specialization for two-phase flow (NUM_PHASES = 2)
 */
template<>
struct TotalMassFlux< 2 > : FluxInterface< 2, FluxComputationType::PARABOLIC_SCALAR >
{
  GEOS_HOST_DEVICE
  static real64 computeTotalMassFlux( real64 const pressure,
                                      real64 const density1,
                                      real64 const density2,
                                      real64 const saturation1,
                                      real64 const saturation2,
                                      real64 const mobility1,
                                      real64 const mobility2,
                                      real64 const gravityTerm );
};

/**
 * @brief Template for hyperbolic flux computations (no capillarity/diffusion)
 */
template< int NUM_PHASES >
using HyperbolicScalar = FluxInterface< NUM_PHASES, FluxComputationType::HYPERBOLIC_SCALAR >;

/**
 * @brief Template for parabolic flux computations (with capillarity/diffusion)
 */
template< int NUM_PHASES >
using ParabolicScalar = FluxInterface< NUM_PHASES, FluxComputationType::PARABOLIC_SCALAR >;

/**
 * @brief Template for mixed hyperbolic-parabolic flux computations
 */
template< int NUM_PHASES >
using HyperbolicParabolicScalar = FluxInterface< NUM_PHASES, FluxComputationType::HYPERBOLIC_PARABOLIC_SCALAR >;

} // namespace fluxTypes

} // namespace geos

#endif /* GEOS_PHYSICSSOLVERS_FLUIDFLOW_FLUXTYPES_HPP_ */