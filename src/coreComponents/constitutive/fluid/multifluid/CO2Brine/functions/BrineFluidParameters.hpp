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
 * @file BrineFluidParameters.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_CO2BRINE_FUNCTIONS_BRINEFLUIDPARAMETERS_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_CO2BRINE_FUNCTIONS_BRINEFLUIDPARAMETERS_HPP_

#include "common/DataTypes.hpp"
#include "common/format/EnumStrings.hpp"

namespace geos
{
namespace constitutive
{
class MultiFluidBase;

namespace PVTProps
{

class PTTableCoordinates;

/// A structure to contain the properties required to create a brine fluid model
struct BrineFluidParameters
{
  // Temperature limits for the correlation used (in C)
  static real64 constexpr minimumTemperature = 10.0;
  static real64 constexpr maximumTemperature = 200.0;

  // Solubility models
  enum class SolubilityModel : integer
  {
    DuanSun,
    SpycherPruess,
    Tables
  };

  /// Solubility model
  SolubilityModel m_solubilityModel{SolubilityModel::DuanSun};

  /// Pressure discretisation points
  array1d< real64 > m_pressureCoordinates;

  /// Interval for pressure discretisation
  real64 m_pressureInterval{0.0};

  /// Temperature discretisation points
  array1d< real64 > m_temperatureCoordinates;

  /// Interval for temperature discretisation
  real64 m_temperatureInterval{0.0};

  /// Names of solubility tables for each phase
  string_array m_solubilityTables;

  /// Water compressibility
  real64 m_waterCompressibility{4.5e-10};

  /// Flash tolerance used to solve any EOS states
  real64 m_tolerance{1.0e-9};

  /// Salinity
  real64 m_salinity{0.0};

  /// The Ezrokhi density Coefficients if required
  array1d< real64 > m_ezrokhiDensityCoefficients;

  /// The Ezrokhi viscosity Coefficients if required
  array1d< real64 > m_ezrokhiViscosityCoefficients;

  struct viewKeyStruct
  {
    static constexpr char const * solubilityModelString() { return "solubilityModel"; }
    static constexpr char const * pressureCoordinatesString() { return "pressureCoordinates"; }
    static constexpr char const * pressureIntervalString() { return "pressureInterval"; }
    static constexpr char const * temperatureCoordinatesString() { return "temperatureCoordinates"; }
    static constexpr char const * temperatureIntervalString() { return "temperatureInterval"; }
    static constexpr char const * salinityString() { return "salinity"; }
    static constexpr char const * waterCompressibilityString() { return "waterCompressibility"; }
    static constexpr char const * solubilityTablesString() { return "solubilityTableNames"; }
    static constexpr char const * toleranceString() { return "tolerance"; }
    static constexpr char const * ezrokhiDensityCoefficientsString() { return "ezrokhiDensityCoefficients"; }
    static constexpr char const * ezrokhiViscosityCoefficientsString() { return "ezrokhiViscosityCoefficients"; }
  };

  /**
   * @brief Attaches the fluid properties to a fluid Group
   * @param[in/out] fluid The fluid object
   */
  template< bool FLASH, bool EZROKHI_DENSITY, bool EZROKHI_VISCOSITY >
  void registerOnFluid( MultiFluidBase * fluid );

  /**
   * @brief Validates the user input for the fluid properties
   * @param[in] fluid The fluid object
   */
  template< bool FLASH, bool EZROKHI_DENSITY, bool EZROKHI_VISCOSITY >
  void postInputInitialization( MultiFluidBase * fluid );

  /**
   * @brief Populate the coordinate table with pressure and temperature
   * @param[in] fluidProperties the user provided properties
   * @param[out] tableCoords the (p,T) coordinates of the table
   * @note This will output temperatures in C
   */
  static void initializePropertyTable( BrineFluidParameters const & fluidParameters,
                                       PTTableCoordinates & tableCoords );

  /**
   * @brief Check if an array has strictly increasing values
   * @param[in] array The array to be checked
   * @return @c true if the array has values that are strictly increasing
   */
  static bool isStrictlyIncreasing( arraySlice1d< real64 const > const & array );
};

ENUM_STRINGS( BrineFluidParameters::SolubilityModel,
              "DuanSun",
              "SpycherPruess",
              "Tables" );

} // end namespace PVTProps

} // end namespace constitutive

} // end namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_CO2BRINE_FUNCTIONS_BRINEFLUIDPARAMETERS_HPP_
