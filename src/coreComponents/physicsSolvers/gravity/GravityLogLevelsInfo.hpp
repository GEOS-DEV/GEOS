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
 * @file GravityLogLevels.hpp
 */



#ifndef GEOS_SOLVERS_GRAVITY_GRAVITYLOGLEVELSINFO_HPP
#define GEOS_SOLVERS_GRAVITY_GRAVITYLOGLEVELSINFO_HPP

#include "common/DataTypes.hpp"

namespace geos
{
namespace gravity
{
/// Log basic solver status and configuration.
struct GravitySolverStatus
{
  static constexpr integer getMinLogLevel() { return 1; }
  static constexpr std::string_view getDescription()
  {
    return "Log basic Gravity solver status and configuration on initialization.";
  }
};

/// Log the computed vertical gravity component (Gz).
struct GravityComponentDebug
{
  static constexpr integer getMinLogLevel() { return 2; }
  static constexpr std::string_view getDescription()
  {
    return "Log the calculated vertical gravity component (Gz) for debugging.";
  }
};

/// Log adjoint-related computations.
struct GravityAdjointDebug
{
  static constexpr integer getMinLogLevel() { return 3; }
  static constexpr std::string_view getDescription()
  {
    return "Log the adjoint calculations for debugging.";
  }
};

/// Log detailed, per-element properties like density.
struct GravityPropertiesDebug
{
  static constexpr integer getMinLogLevel() { return 4; }
  static constexpr std::string_view getDescription()
  {
    return "Log detailed per-element properties (e.g., density).";
  }
};

} // namespace gravity
} // namespace geos

#endif // GEOS_SOLVERS_GRAVITY_GRAVITYLOGLEVELSINFO_HPP
