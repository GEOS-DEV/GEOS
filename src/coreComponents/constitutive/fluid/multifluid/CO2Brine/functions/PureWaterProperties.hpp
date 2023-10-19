/*
 * ------------------------------------------------------------------------------------------------------------
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
 * @file PureWaterProperties.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_CO2BRINE_FUNCTIONS_PUREWATERPROPERTIES_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_CO2BRINE_FUNCTIONS_PUREWATERPROPERTIES_HPP_

#include "PVTFunctionBase.hpp"

#include "functions/FunctionManager.hpp"
#include "functions/TableFunction.hpp"

namespace geos
{

namespace constitutive
{

namespace PVTProps
{

struct PureWaterProperties
{

  /**
   * @brief Creates a table of pure water viscosity [Pa.s] as a function of temperature [degC]
   * @param[in] functionName the name of the PVT function calling this function
   * @param[in] functionManager a reference to the FunctionManager
   * @return a pointer to the newly created TableFunction
   */
  static
  TableFunction const * makeSaturationViscosityTable( string const & functionName,
                                                      FunctionManager & functionManager );

  /**
   * @brief Creates a table of pure water saturation density [kg/m^3] as a function of temperature [degC]
   * @param[in] functionName the name of the PVT function calling this function
   * @param[in] functionManager a reference to the FunctionManager
   * @return a pointer to the newly created TableFunction
   */
  static
  TableFunction const * makeSaturationDensityTable( string const & functionName,
                                                    FunctionManager & functionManager );

  /**
   * @brief Creates a table of pure water saturation pressure [Pa] as a function of temperature [degC]
   * @param[in] functionName the name of the PVT function calling this function
   * @param[in] functionManager a reference to the FunctionManager
   * @return a pointer to the newly created TableFunction
   */
  static
  TableFunction const * makeSaturationPressureTable( string const & functionName,
                                                     FunctionManager & functionManager );

  /// Water molecular weight in kg/mol
  static constexpr real64 MOLECULAR_WEIGHT = 18e-3;
};

} // end namespace PVTProps

} // end namespace constitutive

} // end namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_PVTFUNCTIONS_PUREWATERPROPERTIES_HPP_
