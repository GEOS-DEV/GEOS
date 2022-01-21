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

#ifndef GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_PUREWATERPROPERTIES_HPP_
#define GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_PUREWATERPROPERTIES_HPP_

#include "PVTFunctionBase.hpp"

#include "functions/FunctionManager.hpp"
#include "functions/TableFunction.hpp"

namespace geosx
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
  TableFunction const * makeViscosityTable( string const & functionName,
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

};

} // end namespace PVTProps

} // end namespace constitutive

} // end namespace geosx

#endif //GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_PUREWATERPROPERTIES_HPP_
