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
#include "PressureTemperatureCoordinates.hpp"
#include "BrineSalinity.hpp"

#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
#include "constitutive/fluid/multifluid/CO2Brine/functions/PureWaterProperties.hpp"

#include "constitutive/ExponentialRelation.hpp"

#include "common/Units.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

PhillipsBrineDensityUpdate::PhillipsBrineDensityUpdate( TableFunction const & brineDensityTable,
                                                        integer const waterIndex ):
  m_brineDensityTable( brineDensityTable.createKernelWrapper() ),
  m_waterIndex( waterIndex )
{}

std::unique_ptr< ModelParameters >
PhillipsBrineDensity::createParameters( std::unique_ptr< ModelParameters > parameters )
{
  std::unique_ptr< ModelParameters > params = std::move( parameters );
  params = PressureTemperatureCoordinates::create( std::move( params ) );
  params = BrineSalinity::create( std::move( params ) );
  if( params && params->get< Parameters >() != nullptr )
  {
    return params;
  }
  params = std::make_unique< Parameters >( std::move( params ) );
  return params;
}

PhillipsBrineDensity::PhillipsBrineDensity( string const & name,
                                            ComponentProperties const & componentProperties,
                                            integer const phaseIndex,
                                            ModelParameters const & modelParameters ):
  FunctionBase( name, componentProperties )
{
  GEOS_UNUSED_VAR( phaseIndex );
  
  

}

PhillipsBrineDensity::KernelWrapper
PhillipsBrineDensity::createKernelWrapper() const
{
  return KernelWrapper( *m_brineDensityTable,
                        m_waterIndex );
}

PhillipsBrineDensity::Parameters::Parameters( std::unique_ptr< ModelParameters > parameters ):
  ModelParameters( std::move( parameters ) )
{}

void PhillipsBrineDensity::Parameters::registerParametersImpl( MultiFluidBase * fluid )
{
  fluid->registerWrapper( viewKeyStruct::waterCompressibilityString(), &m_waterCompressibility ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_waterCompressibility ).
    setDescription( "The water compressibility for the Phillips correlation" );
}

void PhillipsBrineDensity::Parameters::postInputInitializationImpl( MultiFluidBase const * fluid,
                                                                    ComponentProperties const & componentProperties )
{
  GEOS_UNUSED_VAR( fluid );
  GEOS_UNUSED_VAR( componentProperties );

  real64 constexpr epsilon = MultiFluidConstants::epsilon;

  GEOS_THROW_IF_LT_MSG( m_waterCompressibility, epsilon,
                        GEOS_FMT( "{}: invalid salinity value provided in '{}'. Compressibility should be positive.",
                                  fluid->getFullName(),
                                  viewKeyStruct::waterCompressibilityString() ),
                        InputError );
}

void PhillipsBrineDensity::calculateBrineDensity( arraySlice1d< real64 const > const & pressureCoords,
                                                  arraySlice1d< real64 const > const & temperatureCoords,
                                                  real64 const & salinity,
                                                  arraySlice1d< real64 > const & densities )
{
  // these coefficients come from Phillips et al. (1981), equations (4) and (5), pages 14 and 15
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

  for( localIndex i = 0; i < nPressures; ++i )
  {
    real64 const pres_in_bar = pressureCoords[i] * PA_2_BAR;

    for( localIndex j = 0; j < nTemperatures; ++j )
    {
      real64 const temperature = units::convertKToC( temperatureCoords[j] );

      real64 const x = c1 * exp( a1 * salinity ) + c2 * exp( a2 * temperature ) + c3 * exp( a3 * pres_in_bar );
      densities[j*nPressures+i] = (AA + x * (BB + x * (CC + x * DD))) * GCC_2_KGM3;
    }
  }
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

  string const functionName = "compositional_phillips_density";

  // Pure water saturation property tables are provided with temperature in C
  TableFunction const * waterSatDensityTable =
    PVTProps::PureWaterProperties::makeSaturationDensityTable( functionName, FunctionManager::getInstance() );
  TableFunction const * waterSatPressureTable =
    PVTProps::PureWaterProperties::makeSaturationPressureTable( functionName, FunctionManager::getInstance() );

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

} // namespace compositional

} // namespace constitutive

} // namespace geos
