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

/** ThOSE 3 macros would replace the ones in Logger.hpp  */
/**
 * @brief Output messages based on current Group's log level.
 * @param[in] logInfo Strut containing log level desscription
 * @param[in] msg a message to log (any expression that can be stream inserted)
 */
#define GEOS_LOG_LEVEL_INFO( logInfo, msg ) GEOS_INFO_IF( isLogLevelActive< logInfo >( this->getLogLevel() ), msg );

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
  static constexpr std::string_view getDescription() { return "Information oon linear system"; }
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
  static constexpr std::string_view getDescription() { return "Print phase volume fraction"; }
};
/// @endcond

}

}

#endif
