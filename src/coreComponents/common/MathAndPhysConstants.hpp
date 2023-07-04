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
 * @file MathAndPhysConstants.hpp
 * @brief Regroups constant values, as well as conversion and calculation functions that are
 * globally used for math and physics computations.
 */
#ifndef GEOS_COMMON_MATH_AND_PHYS_CONSTANTS_HPP_
#define GEOS_COMMON_MATH_AND_PHYS_CONSTANTS_HPP_

namespace geos
{

namespace constants
{

/// @brief Zero degree Celsius in Kelvin
constexpr double zeroCInK = 273.15;

/**
 * @return the input Kelvin degrees converted in Celsius
 * @param kelvin degrees input
 */
inline constexpr double convertKToC( double kelvin )
{ return degrees - zeroKInC; }
/**
 * @return the input Celsius degrees converted in Kelvin
 * @param kelvin degrees input
 */
inline constexpr double convertCToK( double celsius )
{ return celsius + zeroKInC; }

}

}

#endif //GEOS_COMMON_MATH_AND_PHYS_CONSTANTS_HPP_
