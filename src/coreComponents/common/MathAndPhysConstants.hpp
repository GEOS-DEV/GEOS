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

constexpr double zeroKInC = 273.15;

inline constexpr double convertCToK( double degrees )
{ return degrees - zeroKInC; }
inline constexpr double convertKToC( double degrees )
{ return degrees + zeroKInC; }

}

}

#endif //GEOS_COMMON_MATH_AND_PHYS_CONSTANTS_HPP_
