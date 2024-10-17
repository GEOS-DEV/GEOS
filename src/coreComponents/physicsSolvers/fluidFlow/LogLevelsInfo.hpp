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
 * @file LogLevelsInfo.hpp
 * This file contains log level informations for flow solvers
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_LOGLEVELSINFO_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_LOGLEVELSINFO_HPP

#include "common/DataTypes.hpp"

namespace geos
{

namespace logInfo
{

/**
 * @name Wells LogLevels info structures. They must comply with the `is_log_level_info` trait.
 */
///@{

/// @cond DO_NOT_DOCUMENT

struct CFL
{
  static constexpr int getMinLogLevel() { return 1; }
  static constexpr std::string_view getDescription() { return "CFL information"; }
};

struct Statistics
{
  static constexpr int getMinLogLevel() { return 1; }
  static constexpr std::string_view getDescription() { return "Print statistics"; }
};

/// @endcond
///@}

}

}

#endif // GEOS_PHYSICSSOLVERS_FLUIDFLOW_LOGLEVELSINFO_HPP