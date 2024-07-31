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
 * @file LogMessage.hpp
 *
 * Collection of utilities to facilitate I/O of enumeration types.
 * Provides a macro definition that allows associating string names
 * with enumeration constants and a set of functions that make use
 * of these strings, like stream insertion/extraction operators.
 */

#ifndef GEOS_COMMON_LOGMESSAGE_HPP
#define GEOS_COMMON_LOGMESSAGE_HPP

#include "common/DataTypes.hpp"
#include "codingUtilities/EnumBimap.hpp"
#include "codingUtilities/EnumStrings.hpp"
#include "codingUtilities/SourceCodeLocation.hpp"


namespace geos
{ // TODO : document


/// @brief Enumerate the logging levels from the most importants messages to the one that contain the more details.
enum class LogLevel
{
  Silent = -1,    // Useful as a globalLogLevel for silencing the logger. shouldn't be used as a message level.
  Important = 0,  // Application level messages (help, almost blocking warnings, application informations & phases)
  Progress = 1,   // Time step attempts, newton loop progress of Solver objects
  Detailed = 2,   // More detailed info and stats that affect progress, e.g. residual norms.
  Trace = 3,      // This is intended as a user-level debugging tool. detailed info trace what each component is doing. e.g. a line for
                  // every part of the assembly process, every boundary condition that's being applied by a physics solver, etc.
  Debug = 4,      // Information that are only relevant for a developper in a debugging context.
  DebugTrace = 5, // This level is useful to have deeper debugging information (which can be potencially heavy to log).
};
ENUM_BIMAP( LogLevel,
            { LogLevel::Silent, "Silent" },
            { LogLevel::Important, "Important" },
            { LogLevel::Progress, "Progress" },
            { LogLevel::Detailed, "Detailed" },
            { LogLevel::Trace, "Trace" },
            { LogLevel::Debug, "Debug" },
            { LogLevel::DebugTrace, "DebugTrace" } );

struct LogMsgType
{
  Info,
  Warning,
  Error,
};

struct LogMsgCallParams
{
  SourceCodeLocation m_location;
  LogMsgType m_type;
  LogLevel m_logLevel;
};

#define GEOS_LOG_MSG_CTX( type, logLevel ) LogMsgParams( GEOS_SRCLOC(), type, logLevel )

#define GEOS_INFO_IMPORTANT       GEOS_LOG_MSG_PARAMS( LogMsgType::Info, LogLevel::Important )
#define GEOS_INFO_PROGRESS        GEOS_LOG_MSG_PARAMS( LogMsgType::Info, LogLevel::Progress )
#define GEOS_INFO_DETAILED        GEOS_LOG_MSG_PARAMS( LogMsgType::Info, LogLevel::Detailed )
#define GEOS_INFO_TRACE           GEOS_LOG_MSG_PARAMS( LogMsgType::Info, LogLevel::Trace )
#define GEOS_INFO_DEBUG           GEOS_LOG_MSG_PARAMS( LogMsgType::Info, LogLevel::Debug )
#define GEOS_INFO_DEBUGTRACE      GEOS_LOG_MSG_PARAMS( LogMsgType::Info, LogLevel::Debugtrace )

#define GEOS_WARNING_IMPORTANT    GEOS_LOG_MSG_PARAMS( LogMsgType::Warning, LogLevel::Important )
#define GEOS_WARNING_PROGRESS     GEOS_LOG_MSG_PARAMS( LogMsgType::Warning, LogLevel::Progress )
#define GEOS_WARNING_DETAILED     GEOS_LOG_MSG_PARAMS( LogMsgType::Warning, LogLevel::Detailed )
#define GEOS_WARNING_TRACE        GEOS_LOG_MSG_PARAMS( LogMsgType::Warning, LogLevel::Trace )
#define GEOS_WARNING_DEBUG        GEOS_LOG_MSG_PARAMS( LogMsgType::Warning, LogLevel::Debug )
#define GEOS_WARNING_DEBUGTRACE   GEOS_LOG_MSG_PARAMS( LogMsgType::Warning, LogLevel::Debugtrace )

#define GEOS_ERROR_IMPORTANT      GEOS_LOG_MSG_PARAMS( LogMsgType::Error, LogLevel::Important )
#define GEOS_ERROR_PROGRESS       GEOS_LOG_MSG_PARAMS( LogMsgType::Error, LogLevel::Progress )
#define GEOS_ERROR_DETAILED       GEOS_LOG_MSG_PARAMS( LogMsgType::Error, LogLevel::Detailed )
#define GEOS_ERROR_TRACE          GEOS_LOG_MSG_PARAMS( LogMsgType::Error, LogLevel::Trace )
#define GEOS_ERROR_DEBUG          GEOS_LOG_MSG_PARAMS( LogMsgType::Error, LogLevel::Debug )
#define GEOS_ERROR_DEBUGTRACE     GEOS_LOG_MSG_PARAMS( LogMsgType::Error, LogLevel::Debugtrace )



struct LogMsgGeneralContext {
  integer m_rank;
  SystemClock m_rankTimeStamp;
  string m_logSectionTitle;
  real64 m_timeStepStart;
};

struct LogMsgTargetContext {
  string m_targetName;
  string m_targetDataContext;
};


struct LogMsgContext
{
  LogMsgParams m_params;
  LogMsgGeneralContext m_generalContext;
  LogMsgTargetContext m_targetContext;
};


struct LogMsg
{
  string m_text;
  LogMsgContext m_context;
};


}

#endif /* GEOS_COMMON_LOGMESSAGE_HPP */
