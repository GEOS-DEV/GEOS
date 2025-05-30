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
 * @file PhillipsBrineDensity.cpp
 */

#include "PhillipsBrineDensity.hpp"

#include "constitutive/fluid/multifluid/compositional/parameters/EquationOfState.hpp"
#include "constitutive/fluid/multifluid/compositional/parameters/PressureTemperatureCoordinates.hpp"
#include "constitutive/fluid/multifluid/compositional/parameters/BrineSalinity.hpp"
#include "constitutive/fluid/multifluid/compositional/parameters/ImmiscibleWaterParameters.hpp"

#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
#include "constitutive/fluid/multifluid/CO2Brine/functions/PureWaterProperties.hpp"

#include "constitutive/ExponentialRelation.hpp"

#include "functions/FunctionManager.hpp"

#include "common/Units.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

PhillipsBrineDensityUpdate::PhillipsBrineDensityUpdate( TableFunction const & brineVolumeShiftTable,
                                                        integer const waterIndex,
                                                        real64 const salinity,
                                                        real64 const brineMolarWeight,
                                                        EquationOfStateType const equationOfState ):
  m_volumeShiftTable( brineVolumeShiftTable.createKernelWrapper() ),
  m_waterIndex( waterIndex ),
  m_brineMolarWeight( brineMolarWeight ),
  m_salinity( salinity ),
  m_equationOfState( equationOfState )
{}

std::unique_ptr< ModelParameters >
PhillipsBrineDensity::createParameters( std::unique_ptr< ModelParameters > parameters )
{
  std::unique_ptr< ModelParameters > params = std::move( parameters );
  params = EquationOfState::create( std::move( params ) );
  params = PressureTemperatureCoordinates::create( std::move( params ) );
  params = BrineSalinity::create( std::move( params ) );
  return params;
}

PhillipsBrineDensity::PhillipsBrineDensity( string const & name,
                                            ComponentProperties const & componentProperties,
                                            integer const phaseIndex,
                                            ModelParameters const & modelParameters ):
  FunctionBase( name, componentProperties )
{
  GEOS_UNUSED_VAR( phaseIndex );

  m_waterIndex = ImmiscibleWaterParameters::getWaterComponentIndex( componentProperties );
  GEOS_THROW_IF_LT_MSG( m_waterIndex, 0, "Water component not found", InputError );

  EquationOfState const * equationOfState = modelParameters.get< EquationOfState >();
  string const eosName = equationOfState->m_equationsOfStateNames[phaseIndex];
  m_equationOfState = EnumStrings< EquationOfStateType >::fromString( eosName );

  BrineSalinity const * brineSalinity = modelParameters.get< BrineSalinity >();
  m_salinity = brineSalinity->m_salinity;
  real64 const saltMolarWeight = brineSalinity->m_saltMolarWeight;
  real64 const waterMolarWeight = componentProperties.getComponentMolarWeight()[m_waterIndex];

  m_brineMolarWeight = calculateBrineMolarWeight( waterMolarWeight,
                                                  saltMolarWeight,
                                                  m_salinity );

  m_volumeShiftTable = makeVolumeShiftTable( name,
                                             componentProperties,
                                             modelParameters,
                                             m_equationOfState,
                                             m_brineMolarWeight,
                                             m_waterIndex );
}

PhillipsBrineDensity::KernelWrapper
PhillipsBrineDensity::createKernelWrapper() const
{
  return KernelWrapper( *m_volumeShiftTable,
                        m_waterIndex,
                        m_salinity,
                        m_brineMolarWeight,
                        m_equationOfState );
}

void PhillipsBrineDensity::calculateBrineDensity( arraySlice1d< real64 const > const & pressureCoords,
                                                  arraySlice1d< real64 const > const & temperatureCoords,
                                                  real64 const & salinity,
                                                  arraySlice1d< real64 > const & densities )
{
  // These coefficients come from Phillips et al. (1981), equations (4) and (5), pages 14 and 15
  constexpr real64 c1 = -9.9595;
  constexpr real64 c2 = 7.0845;
  constexpr real64 c3 = 3.9093;

  constexpr real64 a1 = -0.004539;
  constexpr real64 a2 = -0.0001638;
  constexpr real64 a3 = 0.00002551;

  constexpr real64 AA = -3.033405;
  constexpr real64 BB = 10.128163;
  constexpr real64 CC = -8.750567;
  constexpr real64 DD = 2.663107;

  localIndex const nPressures = pressureCoords.size();
  localIndex const nTemperatures = temperatureCoords.size();

  // Phillips correlation has pressure in bar, temperature in C and density in gcc
  constexpr real64 PA_2_BAR = 1.0e-5;
  constexpr real64 GCC_2_KGM3 = 1.0e3;

  real64 const x0 = c1 * exp( a1 * salinity );

  for( localIndex i = 0; i < nPressures; ++i )
  {
    real64 const pres_in_bar = pressureCoords[i] * PA_2_BAR;
    real64 const x1 = x0 + c3 * exp( a3 * pres_in_bar );

    for( localIndex j = 0; j < nTemperatures; ++j )
    {
      real64 const temperature = units::convertKToC( temperatureCoords[j] );

      real64 const x = x1 + c2 * exp( a2 * temperature );
      densities[j*nPressures+i] = (AA + x * (BB + x * (CC + x * DD))) * GCC_2_KGM3;
    }
  }
}

/**
   Here's what each part does:
   1. **`salinity * saltMolarWeight`**: This calculates the mass of salt in the brine.
   2. **`1.0 - salinity * saltMolarWeight`**: This calculates the mass of water in the brine
   by subtracting the mass of salt from the total mass of the brine (which is 1.0 kg).
   3. **`(1.0 - salinity * saltMolarWeight) / waterMolarWeight`**: This calculates the moles of
   water by dividing the mass of water by the molar weight of water.
   4. **`salinity + (1.0 - salinity * saltMolarWeight) / waterMolarWeight`**: This adds the
   moles of salt (which is `salinity`) to the moles of water.
   5. **`1.0 / (salinity + (1.0 - salinity * saltMolarWeight) / waterMolarWeight)`**: Finally,
   this calculates the molar weight of the brine by dividing the total mass of the brine
   (1.0 kg) by the total moles of salt and water.
 */
real64 PhillipsBrineDensity::calculateBrineMolarWeight( real64 const & waterMolarWeight,
                                                        real64 const & saltMolarWeight,
                                                        real64 const & salinity )
{
  // Salinity is in mol of salt per kg of brine
  return 1.0 / (salinity + (1.0 - salinity * saltMolarWeight) / waterMolarWeight);
}

std::pair< TableFunction const *, TableFunction const * >
PhillipsBrineDensity::createSaturationTables()
{
  string const functionName = "compositional_phillips_density";

  // Pure water saturation property tables are provided with temperature in C
  TableFunction const * waterSatDensityTable =
    PVTProps::PureWaterProperties::makeSaturationDensityTable( functionName, FunctionManager::getInstance() );
  TableFunction const * waterSatPressureTable =
    PVTProps::PureWaterProperties::makeSaturationPressureTable( functionName, FunctionManager::getInstance() );
  return {waterSatDensityTable, waterSatPressureTable};
}

void PhillipsBrineDensity::calculatePureWaterDensity( arraySlice1d< real64 const > const & pressureCoords,
                                                      arraySlice1d< real64 const > const & temperatureCoords,
                                                      real64 const & compressibility,
                                                      arraySlice1d< real64 > const & densities )
{
  // if no salinity, we fall back to the standard approach in three steps
  // 1- Get the saturation density as a function of temperature
  // 2- Get the saturation pressure as a function of temperature
  // 3- Get the pure water density

  using ExponentialCompute = detail::ExponentialCompute< real64, ExponentApproximationType::Full >;

  auto [waterSatDensityTable, waterSatPressureTable] = createSaturationTables();

  localIndex const nPressures = pressureCoords.size();
  localIndex const nTemperatures = temperatureCoords.size();

  for( localIndex j = 0; j < nTemperatures; ++j )
  {
    real64 const temperature = units::convertKToC( temperatureCoords[j] );

    // Step 1: get the saturation density
    real64 const waterSatDensity = waterSatDensityTable->evaluate( &temperature );
    // Step 2: get the saturation pressure
    real64 const waterSatPressure = waterSatPressureTable->evaluate( &temperature );

    for( localIndex i = 0; i < nPressures; ++i )
    {
      real64 const pressure = pressureCoords[i];
      // Step 3: get the pure water density
      real64 density = 0.0;
      ExponentialCompute::compute( waterSatPressure, waterSatDensity, compressibility, pressure, density );
      densities[j*nPressures+i] = density;
    }
  }
}

void PhillipsBrineDensity::calculateEosWaterMolarVolume( arraySlice1d< real64 const > const & pressureCoords,
                                                         arraySlice1d< real64 const > const & temperatureCoords,
                                                         ComponentProperties const & componentProperties,
                                                         EquationOfStateType const equationOfState,
                                                         real64 const salinity,
                                                         integer const waterIndex,
                                                         arraySlice1d< real64 > const & molarVolume )
{
  integer const numComps = componentProperties.getNumberOfComponents();
  integer const numDofs = 2 + numComps;

  // Create pure water composition
  stackArray1d< real64, MultiFluidConstants::MAX_NUM_COMPONENTS > waterComposition( numComps );
  stackArray1d< real64, 2 + MultiFluidConstants::MAX_NUM_COMPONENTS > tempDerivs( numDofs );
  for( integer ic = 0; ic < numComps; ++ic )
  {
    waterComposition[ic] = 0.0;
  }
  waterComposition[waterIndex] = 1.0;

  localIndex const nPressures = pressureCoords.size();
  localIndex const nTemperatures = temperatureCoords.size();

  auto const & componentPropertiesWrapper = componentProperties.createKernelWrapper();

  auto [waterSatDensityTable, waterSatPressureTable] = createSaturationTables();

  for( localIndex j = 0; j < nTemperatures; ++j )
  {
    real64 const temperature = temperatureCoords[j];

    // For each temperature, if the pressure is below the saturation pressure as calculated from the NIST
    // table, we will not call the EOS because it might give us the gas density. Instead we will replace
    // all values below the saturation pressure by the minimum value calculate at the specified temperature.
    real64 const temp_in_c = units::convertKToC( temperature );

    real64 const waterSatPressure = waterSatPressureTable->evaluate( &temp_in_c );

    real64 minMolarVolume = LvArray::NumericLimits< real64 >::max;
    for( localIndex i = 0; i < nPressures; ++i )
    {
      real64 const pressure = pressureCoords[i];
      if( pressure < waterSatPressure )
      {
        molarVolume[j*nPressures+i] = -1.0;
      }
      else
      {
        // Step 3: get the pure water density
        real64 compressibilityFactor = 0.0;
        CompositionalDensityUpdate::computeCompressibilityFactor( numComps,
                                                                  pressure,
                                                                  temperature,
                                                                  waterComposition.toSliceConst(),
                                                                  componentPropertiesWrapper,
                                                                  equationOfState,
                                                                  salinity,
                                                                  compressibilityFactor,
                                                                  tempDerivs.toSlice() );

        molarVolume[j*nPressures+i] = constants::gasConstant * temperature * compressibilityFactor / pressure;
        minMolarVolume = LvArray::math::min( minMolarVolume, molarVolume[j*nPressures+i] );
      }
    }
    GEOS_THROW_IF_GT_MSG( minMolarVolume, 1.0/MultiFluidConstants::minForSpeciesPresence,
                          GEOS_FMT( "Failed to calculate molar volume for undersaturated pressures below {} at {}.",
                                    units::formatValue( waterSatPressure, units::Unit::Pressure ),
                                    units::formatValue( temperature, units::Unit::Temperature )),
                          InputError );
    for( localIndex i = 0; i < nPressures; ++i )
    {
      if( molarVolume[j*nPressures+i] < 0.0 )
      {
        molarVolume[j*nPressures+i] = minMolarVolume;
      }
    }
  }
}

TableFunction const *
PhillipsBrineDensity::makeVolumeShiftTable( string const & name,
                                            ComponentProperties const & componentProperties,
                                            ModelParameters const & modelParameters,
                                            EquationOfStateType const equationOfState,
                                            real64 const brineMolarWeight,
                                            integer const waterIndex )
{
  FunctionManager & functionManager = FunctionManager::getInstance();

  string const tableName = name + "_brine_volume_shift_table";

  if( functionManager.hasGroup< TableFunction >( tableName ) )
  {
    TableFunction * const volumeShiftTable = functionManager.getGroupPointer< TableFunction >( tableName );
    volumeShiftTable->initializeFunction();
    volumeShiftTable->setDimUnits( { units::Pressure, units::Temperature } );
    volumeShiftTable->setValueUnits( units::MolarVolume );
    return volumeShiftTable;
  }

  PressureTemperatureCoordinates const * coordinates = modelParameters.get< PressureTemperatureCoordinates >();
  auto const pressureCoordinates = coordinates->m_pressureCoordinates.toSliceConst();
  localIndex const nPressures = pressureCoordinates.size();
  GEOS_THROW_IF_LE_MSG( nPressures, 0,
                        GEOS_FMT( "{}: Failed to determine pressure points for Phillips brine density interpolation. "
                                  "Provide values for {}.",
                                  name,
                                  PressureTemperatureCoordinates::viewKeyStruct::pressureCoordinatesString() ),
                        InputError );

  auto const temperatureCoordinates = coordinates->m_temperatureCoordinates.toSliceConst();
  localIndex const nTemperatures = temperatureCoordinates.size();
  GEOS_THROW_IF_LE_MSG( nTemperatures, 0,
                        GEOS_FMT( "{}: Failed to determine temperature points for Phillips brine density interpolation. "
                                  "Provide values for {}.",
                                  name,
                                  PressureTemperatureCoordinates::viewKeyStruct::temperatureCoordinatesString() ),
                        InputError );

  // Calculate the mass density of brine from the Phillips correlation
  array1d< real64 > brineDensity( nPressures*nTemperatures );
  BrineSalinity const * brineSalinity = modelParameters.get< BrineSalinity >();
  real64 const salinity = brineSalinity->m_salinity;
  calculateBrineDensity( pressureCoordinates, temperatureCoordinates, salinity, brineDensity );

  // The Phillips correlation is not accurate for zero or very low salinity
  // If salinity is less than 0.25 (magic value), we will linearly interpolate between
  // the pure water density and the Phillips density
  real64 constexpr salinityCuttoff = 0.25;
  if( salinity < salinityCuttoff )
  {
    array1d< real64 > pureWaterDensities( nPressures*nTemperatures );
    real64 const compressibility = brineSalinity->m_waterCompressibility;
    calculatePureWaterDensity( pressureCoordinates, temperatureCoordinates, compressibility, pureWaterDensities );

    real64 const weight = salinity / salinityCuttoff;
    for( localIndex k = 0; k < nPressures*nTemperatures; ++k )
    {
      brineDensity[k] = weight * brineDensity[k] + (1.0 - weight)*pureWaterDensities[k];
    }
  }

  // Calculate the molar volume of pure water from the equation of state
  array1d< real64 > eosMolarVolume( nPressures*nTemperatures );
  calculateEosWaterMolarVolume( pressureCoordinates,
                                temperatureCoordinates,
                                componentProperties,
                                equationOfState,
                                salinity,
                                waterIndex,
                                eosMolarVolume );

  // Calculate the volume shift by combining the brine density from the Phillips correlation and the molar
  // volume from the equation of state
  array1d< real64 > volumeShift( nPressures*nTemperatures );
  for( localIndex k = 0; k < nPressures*nTemperatures; ++k )
  {
    volumeShift[k] = brineMolarWeight / brineDensity[k] - eosMolarVolume[k];
  }

  // Construct the actual table object
  array1d< real64_array > tableCoords( 2 );
  tableCoords[0].resize( nPressures );
  tableCoords[1].resize( nTemperatures );
  for( localIndex k = 0; k < nPressures; ++k )
  {
    tableCoords[0][k] = pressureCoordinates[k];
  }
  for( localIndex k = 0; k < nTemperatures; ++k )
  {
    tableCoords[1][k] = temperatureCoordinates[k];
  }
  TableFunction * const volumeShiftTable = dynamicCast< TableFunction * >( functionManager.createChild( TableFunction::catalogName(), tableName ) );
  volumeShiftTable->setTableCoordinates( tableCoords, { units::Pressure, units::Temperature } );
  volumeShiftTable->setTableValues( volumeShift, units::MolarVolume );
  volumeShiftTable->setInterpolationMethod( TableFunction::InterpolationType::Linear );
  return volumeShiftTable;
}

} // namespace compositional

} // namespace constitutive

} // namespace geos
