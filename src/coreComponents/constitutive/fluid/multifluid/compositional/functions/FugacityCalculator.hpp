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
 * @file FugacityCalculator.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_FUGACITYCALCULATOR_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_FUGACITYCALCULATOR_HPP_

#include "constitutive/fluid/multifluid/compositional/models/ComponentProperties.hpp"
#include "constitutive/fluid/multifluid/compositional/models/EquationOfState.hpp"
#include "constitutive/fluid/multifluid/compositional/functions/CubicEOSPhaseModel.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

struct FugacityCalculator
{
  /**
   * @brief Calculate the log fugacity for a phase
   * @param[in] numComps number of components
   * @param[in] pressure pressure
   * @param[in] temperature temperature
   * @param[in] composition composition of the phase
   * @param[in] componentProperties The compositional component properties
   * @param[in] equationOfState The equation of state parameters
   * @param[in] phaseIndex The index of the phase within the equation of state parameters
   * @param[out] logFugacity the calculated log fugacity
   */
  template< int USD >
  GEOS_HOST_DEVICE
  static void computeLogFugacityCoefficients( integer const numComps,
                                              real64 const pressure,
                                              real64 const temperature,
                                              arraySlice1d< real64 const, USD > const & composition,
                                              ComponentProperties::KernelWrapper const & componentProperties,
                                              EquationOfState::KernelWrapper const & equationOfState,
                                              integer const phaseIndex,
                                              arraySlice1d< real64 > const & logFugacity );

};

template< int USD >
GEOS_HOST_DEVICE
void FugacityCalculator::computeLogFugacityCoefficients( integer const numComps,
                                                         real64 const pressure,
                                                         real64 const temperature,
                                                         arraySlice1d< real64 const, USD > const & composition,
                                                         ComponentProperties::KernelWrapper const & componentProperties,
                                                         EquationOfState::KernelWrapper const & equationOfState,
                                                         integer const phaseIndex,
                                                         arraySlice1d< real64 > const & logFugacity )
{
  integer const eos_type = equationOfState.m_types[phaseIndex];
  if( EquationOfState::equals( eos_type, EquationOfStateType::PengRobinson ))
  {
    CubicEOSPhaseModel< PengRobinsonEOS >::
    computeLogFugacityCoefficients( numComps,
                                    pressure,
                                    temperature,
                                    composition,
                                    componentProperties,
                                    logFugacity );
  }
  else if( EquationOfState::equals( eos_type, EquationOfStateType::SoaveRedlichKwong ))
  {
    CubicEOSPhaseModel< SoaveRedlichKwongEOS >::
    computeLogFugacityCoefficients( numComps,
                                    pressure,
                                    temperature,
                                    composition,
                                    componentProperties,
                                    logFugacity );
  }
}

} // namespace compositional

} // namespace constitutive

} // namespace geos


#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_FUGACITYCALCULATOR_HPP_
