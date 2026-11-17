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
 * @file FlowSolverUtilities.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_FLOWSOLVERUTILITIES_HPP_
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_FLOWSOLVERUTILITIES_HPP_

#include "common/DataTypes.hpp"
#include "common/format/EnumStrings.hpp"

namespace geos
{

namespace flowSolverUtilities
{

/**
 * @brief Solution scaling type
 */
enum class ScalingType : integer
{
  Global,         ///< Scale the Newton update with a unique scaling factor
  Local            ///< Scale the Newton update locally (modifies the Newton direction)
};

ENUM_STRINGS( ScalingType,
              "Global",
              "Local" );

} // namespace flowSolverUtilities

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_FLOWSOLVERUTILITIES_HPP_
