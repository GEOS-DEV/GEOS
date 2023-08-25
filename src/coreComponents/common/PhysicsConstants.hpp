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
 * @file PhysicsConstants.hpp
 * @brief Regroups useful constants and functions relative to Units Of Measures
 * that are globally used for math and physics computations.
 */
#ifndef GEOS_MATH_PHYSICSCONSTANTS_HPP_
#define GEOS_MATH_PHYSICSCONSTANTS_HPP_

namespace geos
{

namespace constants
{

/// @brief Zero degree Celsius in Kelvin
constexpr double zeroDegreesCelsiusInKelvin = 273.15;

/**
 * @return the input Kelvin degrees converted in Celsius
 * @param kelvin degrees input
 */
inline constexpr double convertKToC( double kelvin )
{ return kelvin - zeroDegreesCelsiusInKelvin; }
/**
 * @return the input Celsius degrees converted in Kelvin
 * @param celsius degrees input
 */
inline constexpr double convertCToK( double celsius )
{ return celsius + zeroDegreesCelsiusInKelvin; }

}

}

#endif //GEOS_MATH_PHYSICSCONSTANTS_HPP_
