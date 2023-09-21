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
 * @brief Regroups useful constants and functions relative to Units Of Measures
 * that are globally used for math and physics computations.
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


/// Enumerator of available unit types. Units are in SI by default.
enum Unit : integer
{
  Unknown,
  Dimensionless,
  Pressure,
  Temperature,
  TemperatureInC,
  Distance,
  Time,
};

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

inline string formatValue( real64 value, Unit dim )
{
  switch( dim )
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
