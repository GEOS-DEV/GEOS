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

template< typename LOG_LEVEL_INFO >
static constexpr bool is_log_level_info =
  std::is_same_v< int, decltype(LOG_LEVEL_INFO::getMinLogLevel()) > &&
  std::is_same_v< string_view, decltype(LOG_LEVEL_INFO::getDescription()) >;

struct LineSearchLogLevel
{
  static constexpr int getMinLogLevel() { return 1; }
  static constexpr string_view getDescription() { return "Information on line search"; }
};

struct LineSearchFailedLogLevel
{
  static constexpr int getMinLogLevel() { return 1; }
  static constexpr string_view getDescription() { return "On Incorrect solution, Information on failed line search"; }
};

struct ScalingFactorLogLevel
{
  static constexpr int getMinLogLevel() { return 1; }
  static constexpr string_view getDescription() { return "Information on global solution scaling factor"; }
};

struct TimeStepLogLevel
{
  static constexpr int getMinLogLevel() { return 1; }
  static constexpr string_view getDescription() { return "Information on the timestep"; }
};

struct SolverTimersLogLevel
{
  static constexpr int getMinLogLevel() { return 1; }
  static constexpr string_view getDescription() { return "Information on solver timers"; }
};

struct ScreenLinearSystemLogLevel
{
  static constexpr int getMinLogLevel() { return 2; }
  static constexpr string_view getDescription() { return "Output to screen the assembled linear system and solutions (matrices and vectors)"; }

};

struct FileLinearSystemLogLevel
{
  static constexpr int getMinLogLevel() { return 3; }
  static constexpr string_view getDescription() { return "Output to file the assembled linear system and solutions (matrices and vectors)"; }
};

struct SolverConfigLogLevel
{
  static constexpr int getMinLogLevel() { return 1; }
  static constexpr string_view getDescription() { return "On non convergance, Information about testing new configuration and print the time step"; }
};

struct ResidualNormLogLevel
{
  static constexpr int getMinLogLevel() { return 1; }
  static constexpr string_view getDescription() { return "Print residual norm"; }
};

struct LinearSystemLogLevel
{
  static constexpr int getMinLogLevel() { return 1; }
  static constexpr string_view getDescription() { return "Information oon linear system"; }
};






