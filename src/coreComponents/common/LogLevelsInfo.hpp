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
 * @file LogLevels.hpp
 */
#ifndef GEOS_COMMON_LOGLEVELSINFO_HPP
#define GEOS_COMMON_LOGLEVELSINFO_HPP

#include "DataTypes.hpp"

namespace geos
{

namespace logInfo
{

/**
 * @brief Trait used to check whether a LOG_LEVEL_INFO structure is valid.
 * @tparam LOG_LEVEL_INFO The log level structure to check.
 */
template< typename LOG_LEVEL_INFO >
static constexpr bool is_log_level_info =
  std::is_same_v< int, decltype(LOG_LEVEL_INFO::getMinLogLevel()) > &&
  std::is_same_v< std::string_view, decltype(LOG_LEVEL_INFO::getDescription()) >;

/** This 3 method would replace the ones in Logger.hpp  */
/**
 * @brief Output messages based on current Group's log level.
 * @param[in] logInfo Strut containing log level desscription
 * @param[in] msg a message to log (any expression that can be stream inserted)
 */
#define GEOS_LOG_INFO_LEVEL( logInfo, msg ) GEOS_INFO_IF( isLogLevelActive< logInfo >( this->getLogLevel() ), msg );

/**
 * @brief Output messages (only on rank 0) based on current Group's log level.
 * @param[in] logInfo Strut containing log level desscription
 * @param[in] msg a message to log (any expression that can be stream inserted)
 */
#define GEOS_LOG_LEVEL_INFO_RANK_0( logInfo, msg ) GEOS_INFO_IF( isLogLevelActive< logInfo >( this->getLogLevel() ), msg );

/**
 * @brief Output messages (with one line per rank) based on current Group's log level.
 * @param[in] logInfo Strut containing log level desscription
 * @param[in] msg a message to log (any expression that can be stream inserted)
 */
#define GEOS_LOG_LEVEL_INFO_BY_RANK( logInfo, msg ) GEOS_INFO_IF( isLogLevelActive< logInfo >( this->getLogLevel() ), msg );


/// @cond DO_NOT_DOCUMENT
struct LineSearch
{
  static constexpr int getMinLogLevel() { return 1; }
  static constexpr std::string_view getDescription() { return "Information on line search"; }
};

/// @cond DO_NOT_DOCUMENT
struct LineSearchFailed
{
  static constexpr int getMinLogLevel() { return 1; }
  static constexpr std::string_view getDescription() { return "On Incorrect solution, Information on failed line search"; }
};

/// @cond DO_NOT_DOCUMENT
struct ScalingFactor
{
  static constexpr int getMinLogLevel() { return 1; }
  static constexpr std::string_view getDescription() { return "Information on global solution scaling factor"; }
};

/// @cond DO_NOT_DOCUMENT
struct TimeStep
{
  static constexpr int getMinLogLevel() { return 1; }
  static constexpr std::string_view getDescription() { return "Information on the timestep"; }
};

/// @cond DO_NOT_DOCUMENT
struct SolverTimers
{
  static constexpr int getMinLogLevel() { return 1; }
  static constexpr std::string_view getDescription() { return "Information on solver timers"; }
};

/// @cond DO_NOT_DOCUMENT
struct ScreenLinearSystem
{
  static constexpr int getMinLogLevel() { return 2; }
  static constexpr std::string_view getDescription() { return "Output to screen the assembled linear system and solutions (matrices and vectors)"; }

};

/// @cond DO_NOT_DOCUMENT
struct FileLinearSystem
{
  static constexpr int getMinLogLevel() { return 3; }
  static constexpr std::string_view getDescription() { return "Output to file the assembled linear system and solutions (matrices and vectors)"; }
};

/// @cond DO_NOT_DOCUMENT
struct SolverConfig
{
  static constexpr int getMinLogLevel() { return 1; }
  static constexpr std::string_view getDescription() { return "On non convergance, Information about testing new configuration and print the time step"; }
};

/// @cond DO_NOT_DOCUMENT
struct ResidualNorm
{
  static constexpr int getMinLogLevel() { return 1; }
  static constexpr std::string_view getDescription() { return "Print residual norm"; }
};

/// @cond DO_NOT_DOCUMENT
struct ResidualValues
{
  static constexpr int getMinLogLevel() { return 1; }
  static constexpr std::string_view getDescription() { return "Print the residual values"; }
};

/// @cond DO_NOT_DOCUMENT
struct LinearSystem
{
  static constexpr int getMinLogLevel() { return 1; }
  static constexpr std::string_view getDescription() { return "Information oon linear system"; }
};

/// @cond DO_NOT_DOCUMENT
struct CrossflowWarning
{
  static constexpr int getMinLogLevel() { return 1; }
  static constexpr std::string_view getDescription() { return "If the well is injector and crossflow enabled, display informations about crossflow for injectors"; }
};

/// @cond DO_NOT_DOCUMENT
struct HydraulicAperture
{
  static constexpr int getMinLogLevel() { return 1; }
  static constexpr std::string_view getDescription() { return "Information on aperture and hydraulic aperture"; }
};

/// @cond DO_NOT_DOCUMENT
struct SolverTimeStep
{
  static constexpr int getMinLogLevel() { return 3; }
  static constexpr std::string_view getDescription() { return "Informations on solver time step"; }
};

/// @cond DO_NOT_DOCUMENT
struct Dof
{
  static constexpr int getMinLogLevel() { return 3; }
  static constexpr std::string_view getDescription() { return "The summary of declared fields and coupling"; }
};

/// @cond DO_NOT_DOCUMENT
struct PoromechanicsPhaseFraction
{
  static constexpr int getMinLogLevel() { return 1; }
  static constexpr std::string_view getDescription() { return "Print phase volume fraction"; }
};

}

}

#endif
