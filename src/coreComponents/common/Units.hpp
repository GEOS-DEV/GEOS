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
 * "characteristic" of "value" ["unit symbol"]
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


using SystemClock = std::chrono::system_clock;

// One year = 365.2425 days (= 146097 / 400).
using YearDaysRatio = std::ratio< 146097, 400 >;
// C++20 equivalent types
using Days = std::chrono::duration< int64_t, std::ratio_multiply< std::ratio< 24 >, std::chrono::hours::period > >;
using Years = std::chrono::duration< int64_t, std::ratio_multiply< YearDaysRatio, Days::period > >;

static constexpr double YearDays = ( double )YearDaysRatio::num / YearDaysRatio::den;
static constexpr double MinuteSeconds = 60.0;
static constexpr double HourSeconds = 60.0 * MinuteSeconds;
static constexpr double DaySeconds = 24.0 * HourSeconds;
static constexpr double YearSeconds = YearDays * DaySeconds;



struct TimeFormatInfo
{
  double m_totalSeconds = 0.0;
  int m_years = 0;
  int m_days = 0;
  int m_hours = 0;
  int m_minutes = 0;
  int m_seconds = 0;


  TimeFormatInfo( double totalSeconds, int years, int days, int hours, int minutes, int seconds );
  static TimeFormatInfo fromSeconds( double const seconds );
  template< typename Duration > static TimeFormatInfo fromDuration( Duration duration );

  friend std::ostream & operator<<( std::ostream & os, TimeFormatInfo const & ctx );

  string toString() const;
};


} // end namespace units

} // end namespace geos


/**
 * @brief Formatter to be able to directly use a DurationInfo as a GEOS_FMT() argument.
 */
template<>
struct GEOS_FMT_NS::formatter< geos::units::TimeFormatInfo > : GEOS_FMT_NS::formatter< std::string >
{
  /**
   * @brief Format the specified TimeFormatInfo to a string.
   * @param durationData the TimeFormatInfo object to format
   * @param ctx formatting state consisting of the formatting arguments and the output iterator
   * @return iterator to the output buffer
   */
  auto format( geos::units::TimeFormatInfo const & durationData, format_context & ctx )
  {
    return GEOS_FMT_NS::formatter< std::string >::format( durationData.toString(), ctx );
  }
};

#endif //GEOS_MATH_PHYSICSCONSTANTS_HPP_
