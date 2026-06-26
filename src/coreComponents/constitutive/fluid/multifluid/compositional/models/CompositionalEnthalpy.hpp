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
  CompositionalEnthalpyUpdate( integer const phaseIndex,
                               integer const vapourIndex,
                               PhaseType const phaseType,
                               EquationOfStateType const equationOfState,
                               real64 const refTemperature,
                               arrayView2d< real64 const > const & referenceEnthalpy,
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
   * @brief Evaluates the Poling polynomial at a given temperature given a list of coefficients
   * @param[in] dT - the temperature difference between the temperature and the reference
   * @param[in] a - the coefficients
   * @param[out] enthalpy - the enthalpy
   * @param[out] heatCapacity - the enthalpy derivative wrt temperature
   */
  GEOS_FORCE_INLINE
  GEOS_HOST_DEVICE
  static void evaluatePolynomial( real64 const & dT,
                                  arraySlice1d< real64 const > const & a,
                                  real64 & enthalpy,
                                  real64 & heatCapacity )
  {
    constexpr real64 r2 = 1.0 / 2.0;
    constexpr real64 r3 = 1.0 / 3.0;
    constexpr real64 r4 = 1.0 / 4.0;
    constexpr real64 r5 = 1.0 / 5.0;

    // c(T) = a[0] + a[1]*dT + a[2]*dT^2 + a[3]*dT^3 + a[4]*dT^4
    heatCapacity = a[0] + dT * (a[1] + dT * (a[2] + dT * (a[3] + dT * a[4])));

    // Evaluate enthalpy using the integral of c(s) ds from Tref to T
    // H(T) = a[0]*dT + (a[1]/2)*dT^2 + (a[2]/3)*dT^3 + (a[3]/4)*dT^4 + (a[4]/5)*dT^5
    enthalpy = dT * (a[0] + dT * (a[1] * r2 + dT * (a[2] * r3 + dT * (a[3] * r4 + dT * (a[4] * r5)))));
  }

  template< typename LAMBDA >
  GEOS_HOST_DEVICE
  static void selectEquationOfState( EquationOfStateType equationOfState, LAMBDA && lambda );

  template< typename EOS, integer USD1, integer USD2 >
  GEOS_HOST_DEVICE
  void calculateEquationOfStateEnthalpy( integer const numComps,
                                         ComponentProperties::KernelWrapper const & componentProperties,
                                         real64 const & pressure,
                                         real64 const & temperature,
                                         arraySlice1d< real64 const, USD1 > const & phaseComposition,
                                         real64 & enthalpy,
                                         arraySlice1d< real64, USD2 > const & dEnthalpy,
                                         PhaseType const & phaseType ) const;

private:
  integer const m_phaseIndex{-1};
  integer const m_vapourIndex{-1};
  PhaseType const m_phaseType;
  EquationOfStateType const m_equationOfState;
  real64 const m_refTemperature;
  arrayView2d< real64 const > const m_referenceEnthalpy;
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
  integer m_phaseIndex{-1};
  integer m_vapourIndex{-1};
  PhaseType m_phaseType;
  EquationOfStateType m_equationOfState;
  arrayView1d< real64 const > m_criticalTemperature;
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
  constexpr integer maxNumDof = MultiFluidConstants::MAX_NUM_COMPONENTS + 2;
  integer const numComps = componentProperties.m_componentMolarWeight.size();
  integer const numDofs = numComps + 2;

  auto const & componentMolarWeight = componentProperties.m_componentMolarWeight;
  auto const & criticalTemperature = componentProperties.m_componentCriticalTemperature;

  // 1. Calculate the polynomial (ideal) enthalpy
  real64 const dT = temperature - m_refTemperature;
  real64 & hIdealEnthalpy = enthalpy;
  auto const & dhIdealEnthalpy = dEnthalpy;
  hIdealEnthalpy = 0.0;
  dhIdealEnthalpy[Deriv::dT] = 0.0;
  dhIdealEnthalpy[Deriv::dP] = 0.0;
  for( integer ic = 0; ic < numComps; ++ic )
  {
    real64 enthalpyI = 0.0;
    real64 heatCapacityI = 0.0;
    evaluatePolynomial( dT, m_coefficients[ic], enthalpyI, heatCapacityI );
    // If temperature is greater that critical temperature, use the gas enthalpy if the gas phase exists
    if( criticalTemperature[ic] < temperature && 0 <= m_vapourIndex )
    {
      enthalpyI += m_referenceEnthalpy( m_vapourIndex, ic );
    }
    else
    {
      enthalpyI += m_referenceEnthalpy( m_phaseIndex, ic );
    }
    hIdealEnthalpy += phaseComposition[ic] * enthalpyI;
    dhIdealEnthalpy[Deriv::dT] += phaseComposition[ic] * heatCapacityI;
    dhIdealEnthalpy[Deriv::dC+ic] = enthalpyI;
  }

  // 2. Calculate the departure enthalpy from the EoS
  real64 hEosEnthalpy = 0.0;
  stackArray1d< real64, maxNumDof > dhEosEnthalpy( numDofs );
  LvArray::forValuesInSlice( dhEosEnthalpy.toSlice(), setZero );
  selectEquationOfState( m_equationOfState, [&]( auto eos )
  {
    using CubicModel = CubicEOSPhaseModel< decltype(eos) >;
    calculateEquationOfStateEnthalpy< CubicModel >( numComps,
                                                    componentProperties,
                                                    pressure,
                                                    temperature,
                                                    phaseComposition,
                                                    hEosEnthalpy,
                                                    dhEosEnthalpy.toSlice(),
                                                    m_phaseType );
  } );

  // 3. Combine the 2 values for the actual fluid enthalpy
  enthalpy += hEosEnthalpy;
  for( integer idof = 0; idof < numDofs; idof++ )
  {
    dEnthalpy[idof] += dhEosEnthalpy[idof];
  }

  // In the mass formulation, the enthalpy will be required to be (annoyingly) in J/kg. The
  // composition however is preconverted to molar composition before we get to this point.
  // The derivatives will also be output in molar composition and conversion is done at a
  // higher level within the fluid model itself. The only change here is to the "unit" of
  // enthalpy.
  // Convert from J/mol -> J/kg
  if( useMass )
  {
    real64 phaseMolarWeight = 0.0;
    for( integer ic = 0; ic < numComps; ++ic )
    {
      phaseMolarWeight += phaseComposition[ic] * componentMolarWeight[ic];
    }
    real64 const invPhaseMolarWeight = 1.0 / phaseMolarWeight;
    enthalpy *= invPhaseMolarWeight;
    dEnthalpy[Deriv::dP] *= invPhaseMolarWeight;
    dEnthalpy[Deriv::dT] *= invPhaseMolarWeight;
    for( integer ic = 0; ic < numComps; ++ic )
    {
      dEnthalpy[Deriv::dC+ic] = (dEnthalpy[Deriv::dC+ic] - componentMolarWeight[ic] * enthalpy)*invPhaseMolarWeight;
    }
  }
}

template< typename LAMBDA >
GEOS_HOST_DEVICE
void CompositionalEnthalpyUpdate::selectEquationOfState( EquationOfStateType equationOfState, LAMBDA && lambda )
{
  if( equationOfState == EquationOfStateType::SoaveRedlichKwong )
  {
    std::forward< LAMBDA >( lambda )( SoaveRedlichKwongEOS{} );
  }
  else if( equationOfState == EquationOfStateType::SoreideWhitson )
  {
    // Soreide-Whitson uses underying Peng-Robinson formulation
    std::forward< LAMBDA >( lambda )( PengRobinsonEOS{} );
  }
  else
  {
    std::forward< LAMBDA >( lambda )( PengRobinsonEOS{} );
  }
}

template< typename EOS, integer USD1, integer USD2 >
GEOS_HOST_DEVICE
void CompositionalEnthalpyUpdate::calculateEquationOfStateEnthalpy(
  integer const numComps,
  ComponentProperties::KernelWrapper const & componentProperties,
  real64 const & pressure,
  real64 const & temperature,
  arraySlice1d< real64 const, USD1 > const & phaseComposition,
  real64 & enthalpy,
  arraySlice1d< real64, USD2 > const & dEnthalpy,
  PhaseType const & phaseType ) const
{
  integer sizes[2] = {0, 0};
  arraySlice2d< real64 const > derivs( nullptr, sizes, sizes );
  typename EOS::template StackVariables< true > stack( numComps,
                                                       componentProperties.m_componentBinaryCoeff.toSliceConst(),
                                                       derivs.toSliceConst() );

  typename EOS::SelectedRoot const root = phaseType == PhaseType::LIQUID ?
                                          EOS::SelectedRoot::MINIMUM : EOS::SelectedRoot::MAXIMUM;

  EOS::initialiseStack( numComps, pressure, temperature, componentProperties, stack );
  EOS::computeMixtureCoefficients( numComps, phaseComposition, stack );
  EOS::template computeEnthalpy< USD1, true >( numComps,
                                               pressure,
                                               temperature,
                                               phaseComposition,
                                               componentProperties,
                                               stack,
                                               enthalpy,
                                               dEnthalpy,
                                               root );
}

} // end namespace compositional

} // end namespace constitutive

} // end namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_COMPOSITIONALENTHALPY_HPP_
