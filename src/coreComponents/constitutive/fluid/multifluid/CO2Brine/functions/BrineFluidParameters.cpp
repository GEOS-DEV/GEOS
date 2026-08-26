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
 * @file BrineFluidParameters.cpp
 */

#include "BrineFluidParameters.hpp"
#include "constitutive/fluid/multifluid/CO2Brine/functions/PVTFunctionHelpers.hpp"
#include "constitutive/fluid/multifluid/MultiFluidConstants.hpp"
#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"

#include "dataRepository/InputFlags.hpp"
#include "dataRepository/RestartFlags.hpp"
#include "common/logger/Logger.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{
namespace PVTProps
{

template< bool FLASH, bool EZROKHI_DENSITY, bool EZROKHI_VISCOSITY >
void BrineFluidParameters::registerOnFluid( MultiFluidBase * fluid )
{
  if constexpr (FLASH)
  {
    fluid->registerWrapper( viewKeyStruct::solubilityModelString(), &m_solubilityModel ).
      setInputFlag( InputFlags::OPTIONAL ).
      setRestartFlags( RestartFlags::NO_WRITE ).
      setDescription( "The solubility model for calculating phase equilibrium. "
                      "Available options are: " + EnumStrings< SolubilityModel >::concat( "\n* " ) ).
      setDefaultValue( m_solubilityModel );

    fluid->registerWrapper( viewKeyStruct::solubilityTablesString(), &m_solubilityTables ).
      setInputFlag( InputFlags::OPTIONAL ).
      setRestartFlags( RestartFlags::NO_WRITE ).
      setDescription( "Names of solubility tables for each phase" );
  }

  fluid->registerWrapper( viewKeyStruct::pressureCoordinatesString(), &m_pressureCoordinates ).
    setInputFlag( InputFlags::REQUIRED ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "List of pressure values for interpolation of function values." );

  fluid->registerWrapper( viewKeyStruct::pressureIntervalString(), &m_pressureInterval ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( GEOS_FMT( "Step size used to evenly discretize the range pressure range given by "
                              "{0} into a set of evenly spaced values. If a positive is givem, it defines the "
                              "spacing between generated points across the full range from the first to the "
                              "last input value. If zero or a negative number is given then the points in {0} "
                              " are left as is.", viewKeyStruct::pressureIntervalString() )).
    setDefaultValue( m_pressureInterval );

  fluid->registerWrapper( viewKeyStruct::temperatureCoordinatesString(), &m_temperatureCoordinates ).
    setInputFlag( InputFlags::REQUIRED ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "List of temperature values for interpolation of function values." );

  fluid->registerWrapper( viewKeyStruct::temperatureIntervalString(), &m_temperatureInterval ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( GEOS_FMT( "Step size used to evenly discretize the range temperature range given by "
                              "{0} into a set of evenly spaced values. If a positive is givem, it defines the "
                              "spacing between generated points across the full range from the first to the "
                              "last input value. If zero or a negative number is given then the points in {0} "
                              " are left as is.", viewKeyStruct::temperatureIntervalString() )).
    setDefaultValue( m_temperatureInterval );

  fluid->registerWrapper( viewKeyStruct::waterCompressibilityString(), &m_waterCompressibility ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "The compressibility of pure water" ).
    setDefaultValue( m_waterCompressibility );

  fluid->registerWrapper( viewKeyStruct::salinityString(), &m_salinity ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "The salinity of brine" ).
    setDefaultValue( m_salinity );

  fluid->registerWrapper( viewKeyStruct::toleranceString(), &m_tolerance ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Tolerance to used for fluid property calculation" ).
    setDefaultValue( m_tolerance );

  // Check if we have the Ezrokhi model
  if constexpr ( EZROKHI_DENSITY )
  {
    fluid->registerWrapper( viewKeyStruct::ezrokhiDensityCoefficientsString(), &m_ezrokhiDensityCoefficients ).
      setInputFlag( InputFlags::OPTIONAL ).
      setRestartFlags( RestartFlags::NO_WRITE ).
      setDescription( "Ezrokhi correlation coefficients for brine density (see :ref:`CO2-EOS` for details)." );
  }

  if constexpr ( EZROKHI_VISCOSITY )
  {
    fluid->registerWrapper( viewKeyStruct::ezrokhiViscosityCoefficientsString(), &m_ezrokhiViscosityCoefficients ).
      setInputFlag( InputFlags::OPTIONAL ).
      setRestartFlags( RestartFlags::NO_WRITE ).
      setDescription( "Ezrokhi correlation coefficients for brine viscosity (see :ref:`CO2-EOS` for details)." );
  }
}

template< bool FLASH, bool EZROKHI_DENSITY, bool EZROKHI_VISCOSITY >
void BrineFluidParameters::postInputInitialization( MultiFluidBase * fluid )
{
  string const fullName = fluid->getFullName();

  if constexpr (FLASH)
  {
    // Check if we need solubility tables
    bool const hasNoTables = m_solubilityTables.empty();
    bool const needsTables = m_solubilityModel == SolubilityModel::Tables;
    GEOS_THROW_IF( needsTables && hasNoTables,
                   GEOS_FMT( "{}: If the solubility model is set to {} then the {} field should be provided.",
                             fullName,
                             m_solubilityModel,
                             viewKeyStruct::solubilityTablesString() ),
                   InputError );

    if( needsTables )
    {
      // The user must provide 1 or 2 tables.
      GEOS_THROW_IF( m_solubilityTables.size() != 1 && m_solubilityTables.size() != 2,
                     GEOS_FMT( "{}: The number of table names in {} must be 1 or 2", fullName, viewKeyStruct::solubilityTablesString() ),
                     InputError );
    }
  }

  // Water compressibility must be positive
  GEOS_THROW_IF_LT_MSG( m_waterCompressibility, MultiFluidConstants::epsilon,
                        GEOS_FMT( "{}: invalid water compressibility {}. "
                                  "Value must be positive", fullName, viewKeyStruct::waterCompressibilityString() ),
                        InputError );

  // Salinity must not be negative
  GEOS_THROW_IF_LT_MSG( m_salinity, 0.0,
                        GEOS_FMT( "{}: invalid salinity {}. "
                                  "Value must not be negative", fullName, viewKeyStruct::salinityString() ),
                        InputError );

  // Flash tolerance must be positive
  GEOS_THROW_IF_LT_MSG( m_tolerance, MultiFluidConstants::epsilon,
                        GEOS_FMT( "{}: invalid flash tolerance {}. "
                                  "Value must be positive", fullName, viewKeyStruct::toleranceString() ),
                        InputError );

  // Coordinates values must be at least 2
  GEOS_THROW_IF_LT_MSG( m_pressureCoordinates.size(), 2,
                        GEOS_FMT( "{}: invalid number of pressure coordinates provided in {}. "
                                  "At least 2 values must be provided", fullName, viewKeyStruct::pressureCoordinatesString() ),
                        InputError );

  // Values must be strictly increasing
  GEOS_THROW_IF( !isStrictlyIncreasing( m_pressureCoordinates.toSliceConst()),
                 GEOS_FMT( "{}: invalid values of pressure coordinates provided in {}. "
                           "Values must be strictly increasing.", fullName, viewKeyStruct::pressureCoordinatesString() ),
                 InputError );

  // Coordinates values must be at least 2
  GEOS_THROW_IF_LT_MSG( m_temperatureCoordinates.size(), 2,
                        GEOS_FMT( "{}: invalid number of temperature coordinates provided in {}. "
                                  "At least 2 values must be provided", fullName, viewKeyStruct::temperatureCoordinatesString() ),
                        InputError );

  // Values must be strictly increasing
  GEOS_THROW_IF( !isStrictlyIncreasing( m_temperatureCoordinates.toSliceConst()),
                 GEOS_FMT( "{}: invalid values of temperature coordinates provided in {}. "
                           "Values must be strictly increasing.", fullName, viewKeyStruct::temperatureCoordinatesString() ),
                 InputError );

  real64 const minTemp = m_temperatureCoordinates[0];
  real64 const maxTemp = m_temperatureCoordinates[m_temperatureCoordinates.size()-1];
  real64 const minTempInK = units::convertCToK( minimumTemperature );
  real64 const maxTempInK = units::convertCToK( maximumTemperature );
  GEOS_THROW_IF_LT_MSG( minTemp, minTempInK,
                        GEOS_FMT( "{}: Minimum temperature must be at least {}K ({} in C). "
                                  "The lowest value provided in {} is {}K", fullName,
                                  minTempInK, minimumTemperature,
                                  viewKeyStruct::temperatureCoordinatesString(), minTemp ),
                        InputError );
  GEOS_THROW_IF_GT_MSG( maxTemp, maxTempInK,
                        GEOS_FMT( "{}: Maximum temperature must be at most {}K ({} in C). "
                                  "The highest value provided in {} is {}K", fullName,
                                  maxTempInK, minimumTemperature,
                                  viewKeyStruct::temperatureCoordinatesString(), maxTemp ),
                        InputError );

  if constexpr ( EZROKHI_DENSITY )
  {
    if( m_ezrokhiDensityCoefficients.empty())
    {
      // Default values for Ezrokhi density coefficients
      m_ezrokhiDensityCoefficients.emplace_back( 0.1003 );
      m_ezrokhiDensityCoefficients.emplace_back( -2.2991e-5 );
      m_ezrokhiDensityCoefficients.emplace_back( -2.3658e-6 );
    }
    GEOS_THROW_IF_NE_MSG( m_ezrokhiDensityCoefficients.size(), 3,
                          GEOS_FMT( "{}: invalid number of Ezrokhi density coefficients provided in {}. "
                                    "Exactly 3 values must be provided.", fullName, viewKeyStruct::ezrokhiDensityCoefficientsString() ),
                          InputError );
  }

  if constexpr ( EZROKHI_VISCOSITY )
  {
    if( m_ezrokhiViscosityCoefficients.empty())
    {
      // Default to zero Ezrokhi viscosity coefficients
      m_ezrokhiViscosityCoefficients.emplace_back( 0.0 );
      m_ezrokhiViscosityCoefficients.emplace_back( 0.0 );
      m_ezrokhiViscosityCoefficients.emplace_back( 0.0 );
    }
    GEOS_THROW_IF_NE_MSG( m_ezrokhiViscosityCoefficients.size(), 3,
                          GEOS_FMT( "{}: invalid number of Ezrokhi viscosity coefficients provided in {}. "
                                    "Exactly 3 values must be provided.", fullName, viewKeyStruct::ezrokhiViscosityCoefficientsString() ),
                          InputError );
  }
}

void BrineFluidParameters::initializePropertyTable( BrineFluidParameters const & fluidParameters,
                                                    PTTableCoordinates & tableCoords )
{
  real64 const dP = fluidParameters.m_pressureInterval;
  if( MultiFluidConstants::epsilon < dP )
  {
    integer const n = fluidParameters.m_pressureCoordinates.size();
    real64 const startPressure = fluidParameters.m_pressureCoordinates[0];
    real64 const endPressure = fluidParameters.m_pressureCoordinates[n-1];
    for( real64 pressure = startPressure; pressure <= endPressure; pressure += dP )
    {
      tableCoords.appendPressure( pressure );
    }
  }
  else
  {
    for( real64 const pressure : fluidParameters.m_pressureCoordinates )
    {
      tableCoords.appendPressure( pressure );
    }
  }

  real64 const dT = fluidParameters.m_temperatureInterval;
  if( MultiFluidConstants::epsilon < dT )
  {
    integer const n = fluidParameters.m_temperatureCoordinates.size();
    real64 const startTemperature = fluidParameters.m_temperatureCoordinates[0];
    real64 const endTemperature = fluidParameters.m_temperatureCoordinates[n-1];
    for( real64 temperature = startTemperature; temperature <= endTemperature; temperature += dT )
    {
      tableCoords.appendTemperature( units::convertKToC( temperature ) );
    }
  }
  else
  {
    for( real64 const temperature : fluidParameters.m_temperatureCoordinates )
    {
      tableCoords.appendTemperature( units::convertKToC( temperature ) );
    }
  }
}

bool BrineFluidParameters::isStrictlyIncreasing( arraySlice1d< real64 const > const & array )
{
  localIndex const size = array.size();
  GEOS_ASSERT( 1 < size );
  real64 constexpr epsilon = MultiFluidConstants::epsilon;
  for( localIndex i = 1; i < size; ++i )
  {
    if( array[i] - array[i-1] < epsilon )
    {
      return false;
    }
  }
  return true;
}

// Explicitly instantiate
template void BrineFluidParameters::registerOnFluid< false, false, false >( MultiFluidBase * );
template void BrineFluidParameters::registerOnFluid< true, false, false >( MultiFluidBase * );
template void BrineFluidParameters::registerOnFluid< true, true, true >( MultiFluidBase * );
template void BrineFluidParameters::postInputInitialization< false, false, false >( MultiFluidBase * );
template void BrineFluidParameters::postInputInitialization< true, false, false >( MultiFluidBase * );
template void BrineFluidParameters::postInputInitialization< true, true, true >( MultiFluidBase * );

} // end namespace PVTProps

} // end namespace constitutive

} // end namespace geos
