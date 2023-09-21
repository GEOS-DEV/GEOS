/*
 * -------------------------------------------------------------------------------------
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
 * @file Units.hpp
 * @brief Enumerates the Units that are in use in GEOS and regroups useful conversion and formatting functions.
 */
#ifndef GEOS_MATH_UNITS_HPP_
#define GEOS_MATH_UNITS_HPP_

#include "common/DataTypes.hpp"
#include "common/PhysicsConstants.hpp"

namespace geos
{

namespace units
{


/**
 * @return the input Kelvin degrees converted in Celsius
 * @param kelvin degrees input
 */
inline constexpr double convertKToC( double kelvin )
{ return kelvin - constants::zeroDegreesCelsiusInKelvin; }
/**
 * @return the input Celsius degrees converted in Kelvin
 * @param celsius degrees input
 */
inline constexpr double convertCToK( double celsius )
{ return celsius + constants::zeroDegreesCelsiusInKelvin; }


/**
 * @brief Enumerator of available unit types. Units are in SI by default.
 */
enum Unit : integer
{
  /// Default label when a unit is not known for a value
  Unknown,

  /// Label to use when a value has not physical dimension (ratio values, propotions...)
  Dimensionless,

  /// Pressure in Pascal
  Pressure,

  /// Temperature in Kelvin
  Temperature,

  /// Temperature in Celcius
  TemperatureInC,

  /// Distance in meter
  Distance,

  /// Time in seconds
  Time,

  /// Viscosity in Pa*s
  Viscosity,

  /// Density in kg/m³
  Density,

  /// Enthalpy in J/kg
  Enthalpy,

  /// Solubility in g/L
  Solubility,
};


/**
 * @param unit The unit we want the information.
 * @return The name of the specified unit followed by its symbol in brackets.
 */
constexpr inline std::string_view getDescription( Unit unit )
{
  switch( unit )
  {
    default:              return "unknown [?]";
    case Dimensionless:   return "dimensionless [1]";
    case Pressure:        return "pressure [Pa]";
    case Temperature:     return "temperature [K]";
    case TemperatureInC:  return "temperature [C]";
    case Distance:        return "distance [m]";
    case Time:            return "time [s]";
  }
}

/**
 * @param unit The unit we want the information.
 * @return The symbol of the specified unit.
 */
constexpr inline std::string_view getSymbol( Unit unit )
{
  switch( unit )
  {
    default:              return "?";
    case Dimensionless:   return "1";
    case Pressure:        return "Pa";
    case Temperature:     return "K";
    case TemperatureInC:  return "C";
    case Distance:        return "m";
    case Time:            return "s";
  }
}

/**
 * @brief Format the specified value coherently with the specified unit.
 * @param value The value to format.
 * @param unit The unit of the specified value.
 * @return A string that can be easily be integrated in a sentence and takes the form 
 * "nature" of "value" ["unit symbol"]
 */
inline string formatValue( real64 value, Unit unit )
{
  switch( unit )
  {
    default:              return GEOS_FMT( "value of {} [?]", value );
    case Dimensionless:   return GEOS_FMT( "value of {} [1]", value );
    case Pressure:        return GEOS_FMT( "pressure of {} [Pa]", value );
    case Temperature:     return GEOS_FMT( "temperature of {} [K]", value );
    case TemperatureInC:  return GEOS_FMT( "temperature of {} [K]", convertCToK( value ) );
    case Distance:        return GEOS_FMT( "distance of {} [s]", value );
    case Time:            return GEOS_FMT( "time of {} [s]", value );
  }
}


} // end namespace units

} // end namespace geos

#endif //GEOS_MATH_PHYSICSCONSTANTS_HPP_
