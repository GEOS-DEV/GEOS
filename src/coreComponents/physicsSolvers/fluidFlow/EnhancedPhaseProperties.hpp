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
 * @file EnhancedPhaseProperties.hpp
 * @brief Enhanced phase model interface with fractional flows and mobility functions
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_ENHANCEDPHASEPROPERTIES_HPP_
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_ENHANCEDPHASEPROPERTIES_HPP_

#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"

namespace geos
{

namespace phaseProperties
{

/**
 * @brief Enhanced phase model interface for multiphase flow
 * @tparam NUM_PHASES Number of phases in the system
 */
template< int NUM_PHASES >
struct EnhancedPhaseModel
{
  static constexpr int numPhases() { return NUM_PHASES; }

  /**
   * @brief Compute fractional flow for a specific phase
   * @param[in] phaseIndex Index of the phase (0 to NUM_PHASES-1)
   * @param[in] mobilities Array of phase mobilities
   * @return Fractional flow for the specified phase
   */
  GEOS_HOST_DEVICE
  static real64 computeFractionalFlow( integer const phaseIndex,
                                       arraySlice1d< real64 const > const & mobilities );

  /**
   * @brief Compute total mobility across all phases
   * @param[in] mobilities Array of phase mobilities
   * @return Total mobility
   */
  GEOS_HOST_DEVICE
  static real64 computeTotalMobility( arraySlice1d< real64 const > const & mobilities );

  /**
   * @brief Compute derivatives of fractional flow with respect to saturations
   * @param[in] phaseIndex Index of the phase
   * @param[in] mobilities Array of phase mobilities
   * @param[in] dMobilities_dSaturations Derivatives of mobilities w.r.t. saturations
   * @param[out] dFractionalFlow_dSaturations Output derivatives
   */
  GEOS_HOST_DEVICE
  static void computeFractionalFlowDerivatives( integer const phaseIndex,
                                                arraySlice1d< real64 const > const & mobilities,
                                                arraySlice2d< real64 const > const & dMobilities_dSaturations,
                                                arraySlice1d< real64 > const & dFractionalFlow_dSaturations );

  /**
   * @brief Compute derivatives of total mobility with respect to saturations
   * @param[in] dMobilities_dSaturations Derivatives of mobilities w.r.t. saturations
   * @param[out] dTotalMobility_dSaturations Output derivatives
   */
  GEOS_HOST_DEVICE
  static void computeTotalMobilityDerivatives( arraySlice2d< real64 const > const & dMobilities_dSaturations,
                                               arraySlice1d< real64 > const & dTotalMobility_dSaturations );
};

/**
 * @brief Mixture properties interface for computing bulk mixture quantities
 * @tparam NUM_PHASES Number of phases
 * @tparam NUM_COMPONENTS Number of components
 */
template< int NUM_PHASES, int NUM_COMPONENTS >
struct MixtureProperties
{
  static constexpr int numPhases() { return NUM_PHASES; }
  static constexpr int numComponents() { return NUM_COMPONENTS; }

  /**
   * @brief Compute mixture density
   * @param[in] phaseDensities Array of phase densities
   * @param[in] saturations Array of phase saturations
   * @return Mixture density
   */
  GEOS_HOST_DEVICE
  static real64 computeMixtureDensity( arraySlice1d< real64 const > const & phaseDensities,
                                       arraySlice1d< real64 const > const & saturations );

  /**
   * @brief Compute mixture enthalpy
   * @param[in] phaseEnthalpies Array of phase enthalpies
   * @param[in] saturations Array of phase saturations
   * @param[in] phaseDensities Array of phase densities
   * @return Mixture enthalpy
   */
  GEOS_HOST_DEVICE
  static real64 computeMixtureEnthalpy( arraySlice1d< real64 const > const & phaseEnthalpies,
                                        arraySlice1d< real64 const > const & saturations,
                                        arraySlice1d< real64 const > const & phaseDensities );

  /**
   * @brief Compute mixture viscosity (harmonic average)
   * @param[in] phaseViscosities Array of phase viscosities
   * @param[in] saturations Array of phase saturations
   * @return Mixture viscosity
   */
  GEOS_HOST_DEVICE
  static real64 computeMixtureViscosity( arraySlice1d< real64 const > const & phaseViscosities,
                                         arraySlice1d< real64 const > const & saturations );

  /**
   * @brief Compute derivatives of mixture density
   * @param[in] phaseDensities Array of phase densities
   * @param[in] saturations Array of phase saturations
   * @param[in] dDensities_dPressure Derivatives of phase densities w.r.t. pressure
   * @param[in] dDensities_dSaturations Derivatives of phase densities w.r.t. saturations
   * @param[out] dMixtureDensity_dPressure Output derivative w.r.t. pressure
   * @param[out] dMixtureDensity_dSaturations Output derivatives w.r.t. saturations
   */
  GEOS_HOST_DEVICE
  static void computeMixtureDensityDerivatives( arraySlice1d< real64 const > const & phaseDensities,
                                                arraySlice1d< real64 const > const & saturations,
                                                arraySlice1d< real64 const > const & dDensities_dPressure,
                                                arraySlice2d< real64 const > const & dDensities_dSaturations,
                                                real64 & dMixtureDensity_dPressure,
                                                arraySlice1d< real64 > const & dMixtureDensity_dSaturations );

  /**
   * @brief Compute derivatives of mixture enthalpy
   * @param[in] phaseEnthalpies Array of phase enthalpies
   * @param[in] saturations Array of phase saturations
   * @param[in] phaseDensities Array of phase densities
   * @param[in] dEnthalpies_dPressure Derivatives of phase enthalpies w.r.t. pressure
   * @param[in] dEnthalpies_dTemperature Derivatives of phase enthalpies w.r.t. temperature
   * @param[out] dMixtureEnthalpy_dPressure Output derivative w.r.t. pressure
   * @param[out] dMixtureEnthalpy_dTemperature Output derivative w.r.t. temperature
   * @param[out] dMixtureEnthalpy_dSaturations Output derivatives w.r.t. saturations
   */
  GEOS_HOST_DEVICE
  static void computeMixtureEnthalpyDerivatives( arraySlice1d< real64 const > const & phaseEnthalpies,
                                                 arraySlice1d< real64 const > const & saturations,
                                                 arraySlice1d< real64 const > const & phaseDensities,
                                                 arraySlice1d< real64 const > const & dEnthalpies_dPressure,
                                                 arraySlice1d< real64 const > const & dEnthalpies_dTemperature,
                                                 real64 & dMixtureEnthalpy_dPressure,
                                                 real64 & dMixtureEnthalpy_dTemperature,
                                                 arraySlice1d< real64 > const & dMixtureEnthalpy_dSaturations );
};

/**
 * @brief Specialization for single-phase systems
 */
template<>
struct EnhancedPhaseModel< 1 >
{
  static constexpr int numPhases() { return 1; }

  GEOS_HOST_DEVICE
  static real64 computeFractionalFlow( integer const phaseIndex,
                                       real64 const mobility )
  {
    GEOS_UNUSED_VAR( phaseIndex );
    GEOS_UNUSED_VAR( mobility );
    return 1.0; // Single phase always has fractional flow = 1
  }

  GEOS_HOST_DEVICE
  static real64 computeTotalMobility( real64 const mobility )
  {
    return mobility;
  }
};

/**
 * @brief Specialization for two-phase systems
 */
template<>
struct EnhancedPhaseModel< 2 >
{
  static constexpr int numPhases() { return 2; }

  GEOS_HOST_DEVICE
  static real64 computeFractionalFlow( integer const phaseIndex,
                                       real64 const mobility1,
                                       real64 const mobility2 );

  GEOS_HOST_DEVICE
  static real64 computeTotalMobility( real64 const mobility1,
                                      real64 const mobility2 )
  {
    return mobility1 + mobility2;
  }

  GEOS_HOST_DEVICE
  static void computeFractionalFlowDerivatives( integer const phaseIndex,
                                                real64 const mobility1,
                                                real64 const mobility2,
                                                real64 const dMobility1_dS1,
                                                real64 const dMobility2_dS1,
                                                real64 & dFractionalFlow_dS1 );
};

/**
 * @brief Capillary pressure utilities
 * @tparam NUM_PHASES Number of phases
 */
template< int NUM_PHASES >
struct CapillaryPressureUtils
{
  /**
   * @brief Compute capillary pressure between phases
   * @param[in] phase1Index Index of first phase
   * @param[in] phase2Index Index of second phase
   * @param[in] saturations Array of phase saturations
   * @param[in] capillaryPressureFunction Function pointer for capillary pressure model
   * @return Capillary pressure between phases
   */
  GEOS_HOST_DEVICE
  static real64 computeCapillaryPressure( integer const phase1Index,
                                          integer const phase2Index,
                                          arraySlice1d< real64 const > const & saturations );

  /**
   * @brief Compute derivatives of capillary pressure
   * @param[in] phase1Index Index of first phase
   * @param[in] phase2Index Index of second phase
   * @param[in] saturations Array of phase saturations
   * @param[out] dPc_dSaturations Derivatives w.r.t. saturations
   */
  GEOS_HOST_DEVICE
  static void computeCapillaryPressureDerivatives( integer const phase1Index,
                                                   integer const phase2Index,
                                                   arraySlice1d< real64 const > const & saturations,
                                                   arraySlice1d< real64 > const & dPc_dSaturations );
};

} // namespace phaseProperties

} // namespace geos

#endif /* GEOS_PHYSICSSOLVERS_FLUIDFLOW_ENHANCEDPHASEPROPERTIES_HPP_ */