/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file PressureTemperatureCoordinates.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_PARAMETERS_PRESSURETEMPERATURECOORDINATES_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_PARAMETERS_PRESSURETEMPERATURECOORDINATES_HPP_

#include "ModelParameters.hpp"
#include "common/DataTypes.hpp"

namespace geos
{

class FunctionBase;

namespace constitutive
{

namespace compositional
{

/**
   Parameter container for pressure and temperature interpolation points
   These are used by functions that are tabulated functions of pressure and temperature. These
   user inputs allow the user to explicitly select the points at which the functions should be
   tabulated.
 */
class PressureTemperatureCoordinates : public ModelParameters
{
public:
  PressureTemperatureCoordinates( std::unique_ptr< ModelParameters > parameters );
  ~PressureTemperatureCoordinates() override = default;

  static std::unique_ptr< ModelParameters > create( std::unique_ptr< ModelParameters > parameters );

  array1d< real64 > m_pressureCoordinates;
  array1d< real64 > m_temperatureCoordinates;

  struct viewKeyStruct
  {
    static constexpr char const * pressureCoordinatesString() { return "pressureCoordinates"; }
    static constexpr char const * temperatureCoordinatesString() { return "temperatureCoordinates"; }
  };

protected:
  void registerParametersImpl( MultiFluidBase * fluid ) override;
  void postInputInitializationImpl( MultiFluidBase const * fluid, ComponentProperties const & componentProperties ) override;

protected:
  /**
   * @brief Check if an array has strictly increasing values
   * @param[in] array The array to be checked
   * @return @c true if the array has values that are strictly increasing
   */
  static bool isStrictlyIncreasing( arraySlice1d< real64 const > const & array );

  /**
   * @brief Get the ordering of input variables on a user supplied function
   * @details Checks a table (or function) of 2 variables pressure and temperature and determines if
              pressure is the first variable. This will return a pair of integers {pIndex, tIndex} which
              will be {0,1} if pressure is the first input variable on the function and {1,0} if
              temperature is the first input variable.
   * @param[in] function The function to be checked
   * @return @c {1,0} if temperature is the first input variable otherwise {0,1}
   */
  static std::pair< integer, integer > getVariableIndices( FunctionBase const * function );
};

} // end namespace compositional

} // end namespace constitutive

} // end namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_PARAMETERS_PRESSURETEMPERATURECOORDINATES_HPP_
