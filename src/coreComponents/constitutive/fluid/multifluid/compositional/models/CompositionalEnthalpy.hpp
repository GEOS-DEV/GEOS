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
 * @file CompositionalEnthalpy.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_COMPOSITIONALENTHALPY_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_COMPOSITIONALENTHALPY_HPP_

#include "FunctionBase.hpp"
#include "constitutive/fluid/multifluid/compositional/parameters/EquationOfState.hpp"
#include "constitutive/fluid/multifluid/compositional/parameters/PhaseType.hpp"

#include "constitutive/fluid/multifluid/MultiFluidUtils.hpp"
#include "constitutive/fluid/multifluid/MultiFluidConstants.hpp"
#include "constitutive/fluid/multifluid/compositional/functions/CompositionalProperties.hpp"
#include "constitutive/fluid/multifluid/compositional/functions/CubicEOSPhaseModel.hpp"
#include "constitutive/fluid/multifluid/compositional/functions/SoreideWhitsonEOSPhaseModel.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

class HeatCapacityCoefficients;

class CompositionalEnthalpyUpdate final : public FunctionBaseUpdate
{
  using Deriv = multifluid::DerivativeOffset;
public:
  CompositionalEnthalpyUpdate( PhaseType const phaseType,
                               EquationOfStateType const equationOfState,
                               arrayView1d< real64 const > const & referenceEnthalpy,
                               arrayView2d< real64 const > const & coefficients );

  template< integer USD1, integer USD2 >
  GEOS_HOST_DEVICE
  void compute( ComponentProperties::KernelWrapper const & componentProperties,
                real64 const & pressure,
                real64 const & temperature,
                arraySlice1d< real64 const, USD1 > const & phaseComposition,
                real64 & enthalpy,
                arraySlice1d< real64, USD2 > const & dEnthalpy,
                bool useMass ) const;

  /**
   * @brief Evaluates the Poling polynomial at a given temeperature given a list of coefficients
   * @param[in] T - the temperature
   * @param[in] a - the coefficients
   * @param[out] enthalpy - the enthalpy
   * @param[out] heatCapacity - the enthalpy derivative wrt temperature
   */
  GEOS_FORCE_INLINE
  GEOS_HOST_DEVICE
  static void evaluatePolynomial( real64 const & T,
                                  arraySlice1d< real64 const > const & a,
                                  real64 & enthalpy,
                                  real64 & heatCapacity )
  {
    real64 constexpr r1 = 1.0/2.0;
    real64 constexpr r2 = 1.0/3.0;
    real64 constexpr r3 = 1.0/4.0;
    real64 constexpr r4 = 1.0/5.0;
    enthalpy = ((((r4*a[4]*T + r3*a[3])*T + r2*a[2])*T + r1*a[1])*T + a[0])*T;
    heatCapacity = (((a[4]*T + a[3])*T + a[2])*T + a[1])*T + a[0];
  }

private:
  PhaseType const m_phaseType;
  EquationOfStateType const m_equationOfState;
  arrayView1d< real64 const > const m_referenceEnthalpy;
  arrayView2d< real64 const > const m_coefficients;
};

class CompositionalEnthalpy : public FunctionBase
{
public:
  CompositionalEnthalpy( string const & name,
                         ComponentProperties const & componentProperties,
                         integer const phaseIndex,
                         ModelParameters const & modelParameters );

  static string catalogName() { return "CompositionalEnthalpy"; }

  virtual FunctionType functionType() const override
  {
    return FunctionType::ENTHALPY;
  }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = CompositionalEnthalpyUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const;

  // Create parameters unique to this model
  static std::unique_ptr< ModelParameters > createParameters( std::unique_ptr< ModelParameters > parameters );

private:
  PhaseType m_phaseType;
  EquationOfStateType m_equationOfState;
  array1d< real64 > m_referenceEnthalpy;
  HeatCapacityCoefficients const * const m_heatCapacityCoefficients{nullptr};
};

template< integer USD1, integer USD2 >
GEOS_HOST_DEVICE
void CompositionalEnthalpyUpdate::compute(
  ComponentProperties::KernelWrapper const & componentProperties,
  real64 const & pressure,
  real64 const & temperature,
  arraySlice1d< real64 const, USD1 > const & phaseComposition,
  real64 & enthalpy,
  arraySlice1d< real64, USD2 > const & dEnthalpy,
  bool useMass ) const
{
  GEOS_UNUSED_VAR( useMass );
  GEOS_UNUSED_VAR( pressure );

  integer const numComps = componentProperties.m_componentMolarWeight.size();

  // 1. Calculate the ideal gas enthalpy
  real64 & hIdealGas = enthalpy;
  auto const & dhIdealGas = dEnthalpy;
  hIdealGas = 0.0;
  dhIdealGas[Deriv::dT] = 0.0;
  dhIdealGas[Deriv::dP] = 0.0;
  for( integer ic = 0; ic < numComps; ++ic )
  {
    real64 enthalpyI = 0.0;
    real64 heatCapacityI = 0.0;
    evaluatePolynomial( temperature, m_coefficients[ic], enthalpyI, heatCapacityI );
    enthalpyI += m_referenceEnthalpy[ic];
    hIdealGas += phaseComposition[ic] * enthalpyI;
    dhIdealGas[Deriv::dT] += phaseComposition[ic] * heatCapacityI;
    dhIdealGas[Deriv::dC+ic] = enthalpyI;
  }
}

} // end namespace compositional

} // end namespace constitutive

} // end namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_COMPOSITIONALENTHALPY_HPP_
