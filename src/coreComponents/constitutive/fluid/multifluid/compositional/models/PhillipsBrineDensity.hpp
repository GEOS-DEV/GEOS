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
 * @file PhillipsBrineDensity.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_PHILLIPSBRINEDENSITY_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_PHILLIPSBRINEDENSITY_HPP_

#include "FunctionBase.hpp"
#include "CompositionalDensity.hpp"

#include "constitutive/fluid/multifluid/Layouts.hpp"
#include "constitutive/fluid/multifluid/MultiFluidUtils.hpp"
#include "constitutive/fluid/multifluid/MultiFluidConstants.hpp"

#include "functions/TableFunction.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

class PhillipsBrineDensityUpdate final : public FunctionBaseUpdate
{
public:
  PhillipsBrineDensityUpdate( TableFunction const & brineVolumeShiftTable,
                              integer const waterIndex,
                              real64 const salinity,
                              real64 const brineMolarWeight,
                              EquationOfStateType const equationOfState );

  template< int USD1, int USD2 >
  GEOS_HOST_DEVICE
  void compute( ComponentProperties::KernelWrapper const & componentProperties,
                real64 const & pressure,
                real64 const & temperature,
                arraySlice1d< real64 const, USD1 > const & phaseComposition,
                real64 & molarDensity,
                arraySlice1d< real64, USD2 > const & dMolarDensity,
                real64 & massDensity,
                arraySlice1d< real64, USD2 > const & dMassDensity,
                bool useMass ) const;

  virtual void move( LvArray::MemorySpace const space, bool const touch ) override
  {
    FunctionBaseUpdate::move( space, touch );
    m_volumeShiftTable.move( space, touch );
  }

protected:
  /// Volume shift of the water component tabulated as a function (P,T)
  TableFunction::KernelWrapper m_volumeShiftTable;

  /// Index of the water component
  integer const m_waterIndex;

  /// The brine molecular weight
  real64 const m_brineMolarWeight;

  /// The salinity
  real64 const m_salinity{0.0};

  /// Equation of state for the density correction
  EquationOfStateType const m_equationOfState;
};

class PhillipsBrineDensity : public FunctionBase
{
public:
  PhillipsBrineDensity( string const & name,
                        ComponentProperties const & componentProperties,
                        integer const phaseIndex,
                        ModelParameters const & modelParameters );

  static string catalogName() { return "PhillipsBrineDensity"; }

  virtual FunctionType functionType() const override
  {
    return FunctionType::DENSITY;
  }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = PhillipsBrineDensityUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const;

  // Create parameters unique to this model
  static std::unique_ptr< ModelParameters > createParameters( std::unique_ptr< ModelParameters > parameters );

private:
  static real64 calculateBrineMolarWeight( real64 const & waterMolarWeight,
                                           real64 const & saltMolarWeight,
                                           real64 const & salinity );

  static std::pair< TableFunction const *, TableFunction const * > createSaturationTables();

  static void calculateBrineDensity( arraySlice1d< real64 const > const & pressureCoords,
                                     arraySlice1d< real64 const > const & temperatureCoords,
                                     real64 const & salinity,
                                     arraySlice1d< real64 > const & densities );

  static void calculatePureWaterDensity( arraySlice1d< real64 const > const & pressureCoords,
                                         arraySlice1d< real64 const > const & temperatureCoords,
                                         real64 const & compressibility,
                                         arraySlice1d< real64 > const & densities );

  static void calculateEosWaterMolarVolume( arraySlice1d< real64 const > const & pressureCoords,
                                            arraySlice1d< real64 const > const & temperatureCoords,
                                            ComponentProperties const & componentProperties,
                                            EquationOfStateType const equationOfState,
                                            real64 const salinity,
                                            integer const waterIndex,
                                            arraySlice1d< real64 > const & molarVolume );

  static TableFunction const * makeVolumeShiftTable( string const & name,
                                                     ComponentProperties const & componentProperties,
                                                     ModelParameters const & modelParameters,
                                                     EquationOfStateType const equationOfState,
                                                     real64 const brineMolarWeight,
                                                     integer const waterIndex );

private:
  /// Volume shift of the water component tabulated as a function (P,T)
  TableFunction const * m_volumeShiftTable;

  /// Index of the water phase
  integer m_waterIndex;

  /// Equation of state for the density correction
  EquationOfStateType m_equationOfState;

  /// The salinity
  real64 m_salinity{0.0};

  /// The brine molecular weight
  real64 m_brineMolarWeight;
};

template< int USD1, int USD2 >
GEOS_HOST_DEVICE
void PhillipsBrineDensityUpdate::compute(
  ComponentProperties::KernelWrapper const & componentProperties,
  real64 const & pressure,
  real64 const & temperature,
  arraySlice1d< real64 const, USD1 > const & phaseComposition,
  real64 & molarDensity,
  arraySlice1d< real64, USD2 > const & dMolarDensity,
  real64 & massDensity,
  arraySlice1d< real64, USD2 > const & dMassDensity,
  bool useMass ) const
{
  using Deriv = constitutive::multifluid::DerivativeOffset;
  GEOS_UNUSED_VAR( useMass );

  integer const numComps = componentProperties.m_componentMolarWeight.size();
  integer const numDofs = 2 + numComps;

  // Calculate the volume shift as a function of (P,T)
  real64 const input[2] = { pressure, temperature };
  real64 brineVolumeShiftDeriv[2]{};
  real64 const brineVolumeShift = m_volumeShiftTable.compute( input, brineVolumeShiftDeriv );

  // Calculate the compressibility factor of the mixture from the equation of state
  // Use molar density space for temporary derivatives
  real64 compressibilityFactor = 0.0;
  arraySlice1d< real64, USD2 > const & dCompressibilityFactor = dMolarDensity;
  CompositionalDensityUpdate::computeCompressibilityFactor( numComps,
                                                            pressure,
                                                            temperature,
                                                            phaseComposition,
                                                            componentProperties,
                                                            m_equationOfState,
                                                            m_salinity,
                                                            compressibilityFactor,
                                                            dCompressibilityFactor );

  // Convert to molar volume by scaling by (RT/P)
  // Scaling factor to convert compressibility factor (Z) to volume.
  real64 const idealGasVolume = constants::gasConstant * temperature  / pressure;
  real64 const dIdealGasVolume_dP = -constants::gasConstant * temperature  / (pressure * pressure);
  real64 const dIdealGasVolume_dT = constants::gasConstant / pressure;

  real64 molarVolume = idealGasVolume * compressibilityFactor;
  arraySlice1d< real64, USD2 > const & dMolarVolume = dMolarDensity;
  dMolarVolume[Deriv::dP] = idealGasVolume * dCompressibilityFactor[Deriv::dP] + dIdealGasVolume_dP * compressibilityFactor;
  dMolarVolume[Deriv::dT] = idealGasVolume * dCompressibilityFactor[Deriv::dT] + dIdealGasVolume_dT * compressibilityFactor;
  for( integer ic = 0; ic < numComps; ++ic )
  {
    dMolarVolume[Deriv::dC + ic] = idealGasVolume * dCompressibilityFactor[Deriv::dC + ic];
  }

  // Extract the water mole fraction
  real64 const x_h2o = phaseComposition[m_waterIndex];

  // Apply the volume shift
  molarVolume += x_h2o * brineVolumeShift;
  dMolarVolume[Deriv::dP] += x_h2o * brineVolumeShiftDeriv[0];
  dMolarVolume[Deriv::dT] += x_h2o * brineVolumeShiftDeriv[1];
  dMolarVolume[Deriv::dC + m_waterIndex] += brineVolumeShift;

  // Invert molar volume to get the molar density
  molarDensity = 1.0 / molarVolume;
  real64 const sqrMolarDen = molarDensity * molarDensity;
  for( integer idof = 0; idof < numDofs; ++idof )
  {
    dMolarDensity[idof] = -sqrMolarDen * dMolarVolume[idof];
  }

  // Calculate the mass density
  auto const & componentMolarWeight = componentProperties.m_componentMolarWeight;
  real64 phaseMolarWeight = 0.0;
  arraySlice1d< real64, USD2 > const & dPhaseMolarWeight = dMassDensity;
  phaseMolarWeight = m_brineMolarWeight * x_h2o;
  dPhaseMolarWeight[Deriv::dP] = 0.0;
  dPhaseMolarWeight[Deriv::dT] = 0.0;
  dPhaseMolarWeight[Deriv::dC + m_waterIndex] = m_brineMolarWeight;

  for( integer ic = 0; ic < m_waterIndex; ++ic )
  {
    phaseMolarWeight += componentMolarWeight[ic] * phaseComposition[ic];
    dPhaseMolarWeight[Deriv::dC + ic] = componentMolarWeight[ic];
  }
  for( integer ic = m_waterIndex+1; ic < numComps; ++ic )
  {
    phaseMolarWeight += componentMolarWeight[ic] * phaseComposition[ic];
    dPhaseMolarWeight[Deriv::dC + ic] = componentMolarWeight[ic];
  }

  // Multiply by the molar density
  massDensity = phaseMolarWeight * molarDensity;
  for( integer idof = 0; idof < numDofs; ++idof )
  {
    dMassDensity[idof] = phaseMolarWeight * dMolarDensity[idof] + dPhaseMolarWeight[idof] * molarDensity;
  }
}

} // end namespace compositional

} // end namespace constitutive

} // end namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_PHILLIPSBRINEDENSITY_HPP_
