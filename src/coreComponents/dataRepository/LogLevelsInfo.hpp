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
 * This file contains all log level informations and the mecanism to ensure LOG_LEVEL_INFO structure is valid
 */
#ifndef GEOS_COMMON_LOGLEVELSINFO_HPP
#define GEOS_COMMON_LOGLEVELSINFO_HPP

#include "common/DataTypes.hpp"
#include "common/format/Format.hpp"

namespace geos
{

namespace logInfo
{

/**
 * @name Common LogLevels info structures. They must comply with the `is_log_level_info` trait.
 */
///@{

/// @cond DO_NOT_DOCUMENT
struct LineSearch
{
  static constexpr int getMinLogLevel() { return 1; }
  static constexpr std::string_view getDescription() { return "Information on line search"; }
};

struct LineSearchFailed
{
  static constexpr int getMinLogLevel() { return 1; }
  static constexpr std::string_view getDescription() { return "On Incorrect solution, Information on failed line search"; }
};

struct ScalingFactor
{
  static constexpr int getMinLogLevel() { return 1; }
  static constexpr std::string_view getDescription() { return "Information on global solution scaling factor"; }
};

struct TimeStep
{
  static constexpr int getMinLogLevel() { return 1; }
  static constexpr std::string_view getDescription() { return "Information on the timestep"; }
};

struct SolverTimers
{
  static constexpr int getMinLogLevel() { return 1; }
  static constexpr std::string_view getDescription() { return "Information on solver timers"; }
};

struct ScreenLinearSystem
{
  static constexpr int getMinLogLevel() { return 1; }
  static constexpr std::string_view getDescription() { return "Output to screen the assembled linear system and solutions (matrices and vectors)"; }

};

struct FileLinearSystem
{
  static constexpr int getMinLogLevel() { return 2; }
  static constexpr std::string_view getDescription() { return "Output to file the assembled linear system and solutions (matrices and vectors)"; }
};

struct SolverConfig
{
  static constexpr int getMinLogLevel() { return 1; }
  static constexpr std::string_view getDescription() { return "On non convergance, Information about testing new configuration and print the time step"; }
};

struct ResidualNorm
{
  static constexpr int getMinLogLevel() { return 1; }
  static constexpr std::string_view getDescription() { return "Print residual norm"; }
};

struct ResidualValues
{
  static constexpr int getMinLogLevel() { return 1; }
  static constexpr std::string_view getDescription() { return "Print the residual values"; }
};

struct LinearSystem
{
  static constexpr int getMinLogLevel() { return 1; }
  static constexpr std::string_view getDescription() { return "Information on linear system"; }
};

struct CrossflowWarning
{
  static constexpr int getMinLogLevel() { return 1; }
  static constexpr std::string_view getDescription() { return "If the well is injector and crossflow enabled, display informations about crossflow for injectors"; }
};

struct HydraulicAperture
{
  static constexpr int getMinLogLevel() { return 1; }
  static constexpr std::string_view getDescription() { return "Information on aperture and hydraulic aperture"; }
};

struct SolverTimeStep
{
  static constexpr int getMinLogLevel() { return 1; }
  static constexpr std::string_view getDescription() { return "Informations on solver time step"; }
};

struct SolverNextDt
{
  static constexpr int getMinLogLevel() { return 3; }
  static constexpr std::string_view getDescription() { return "Informations on changing DT"; }
};

struct Dof
{
  static constexpr int getMinLogLevel() { return 2; }
  static constexpr std::string_view getDescription() { return "The summary of declared fields and coupling"; }
};

struct PoromechanicsPhaseFraction
{
  static constexpr int getMinLogLevel() { return 1; }
  static constexpr std::string_view getDescription() { return "Print poromechanics phase volume fraction"; }
};
/// @endcond
///@}

}


/**
 * @brief Trait used to check whether a LOG_LEVEL_INFO structure is valid.
 * @tparam LOG_LEVEL_INFO The log level structure to check.
 */
template< typename LOG_LEVEL_INFO >
static constexpr bool is_log_level_info =
  std::is_same_v< int, decltype(LOG_LEVEL_INFO::getMinLogLevel()) > &&
  std::is_same_v< std::string_view, decltype(LOG_LEVEL_INFO::getDescription()) >;

/**
 * @brief Verify if a log level is active
 * @tparam LOG_LEVEL_INFO The structure containing log level information.
 * @param level Log level to be checked.
 * @return `true` if the log level is active, `false` otherwise.
 * @pre `LOG_LEVEL_INFO` must satisfy `logInfo::is_log_level_info`.
 *
 */
template< typename LOG_LEVEL_INFO >
std::enable_if_t< geos::is_log_level_info< LOG_LEVEL_INFO >, bool >
isLogLevelActive( int level )
{
  return level >= LOG_LEVEL_INFO::getMinLogLevel();
}

/** ThOSE 3 macros would replace the ones in Logger.hpp  */
/**
 * @brief Output messages based on current Group's log level.
 * @param[in] logInfoStruct Strut containing log level desscription
 * @param[in] msg a message to log (any expression that can be stream inserted)
 */
#define GEOS_LOG_LEVEL_INFO( logInfoStruct, msg ) GEOS_LOG_IF( isLogLevelActive< logInfoStruct >( this->getLogLevel() ), msg );

/**
 * @brief Output messages (only on rank 0) based on current Group's log level.
 * @param[in] logInfoStruct Strut containing log level desscription
 * @param[in] msg a message to log (any expression that can be stream inserted)
 */
#define GEOS_LOG_LEVEL_INFO_RANK_0( logInfoStruct, msg ) GEOS_LOG_RANK_0_IF( isLogLevelActive< logInfoStruct >( this->getLogLevel() ), msg );

/**
 * @brief Output messages (with one line per rank) based on current Group's log level.
 * @param[in] logInfoStruct Strut containing log level desscription
 * @param[in] msg a message to log (any expression that can be stream inserted)
 */
#define GEOS_LOG_LEVEL_INFO_BY_RANK( logInfoStruct, msg ) GEOS_LOG_RANK_IF( isLogLevelActive< logInfoStruct >( this->getLogLevel() ), msg );

}

#endif
