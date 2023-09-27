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
 * @file Logger.hpp
 */

#ifndef GEOS_COMMON_LOGGER_HPP
#define GEOS_COMMON_LOGGER_HPP

// Source incldes
#include "common/GeosxConfig.hpp"
#include "common/GeosxMacros.hpp"
#include "common/Format.hpp"
#include "codingUtilities/EnumBimap.hpp"
#include "LvArray/src/Macros.hpp"

// System includes
#include <stdexcept>
#include <memory>

#if defined(GEOSX_USE_MPI)
  #include <mpi.h>
#endif

/**
 * @brief Log a message on screen.
 * @details The expression to log must evaluate something that can be stream inserted.
 */
//deprecated// #define GEOS_LOG( ... ) LVARRAY_LOG( __VA_ARGS__ )

/**
 * @brief Log an expression and its value on screen.
 * @details The expression to log must evaluate something that can be stream inserted.
 */
#define GEOS_LOG_VAR( ... ) logger.stdLog( STRINGIZE( __VA_ARGS__ ) " = ", __VA_ARGS__ )

/**
 * @brief Conditionally log a message on screen on rank 0.
 * @param EXP an expression that will be evaluated as a predicate
 * @param msg a message to log (any expression that can be stream inserted)
 */
//deprecated// #define GEOS_LOG_RANK_0_IF( EXP, msg ) \/
// do { \/
//   if( ::geos::logger.rank == 0 && EXP ) \/
//   { \/
//     std::ostringstream oss; \/
//     oss << msg; \/
//     std::cout << oss.str() << std::endl; \/
//   } \/
// } while( false )

/**
 * @brief Log a message on screen on rank 0.
 * @param msg a message to log (any expression that can be stream inserted)
 */
//deprecated// #define GEOS_LOG_RANK_0( msg ) GEOS_LOG_RANK_0_IF( true, msg )

/**
 * @brief Conditionally log a message to the rank output stream.
 * @param EXP an expression that will be evaluated as a predicate
 * @param msg a message to log (any expression that can be stream inserted)
 */
// #if defined(GEOS_DEVICE_COMPILE)
// //deprecated// #define GEOS_LOG_RANK_IF( EXP, msg )
// #else
// //deprecated// #define GEOS_LOG_RANK_IF( EXP, msg ) \/
// do { \/
//   if( EXP ) \/
//   { \/
//     std::ostringstream oss; \/
//     oss << ::geos::logger.rankMsgPrefix << msg; \/
//     *logger.outStream << oss.str() << std::endl; \/
//   } \/
// } while( false )
// #endif

/**
 * @brief Log a message to the rank output stream.
 * @param msg a message to log (any expression that can be stream inserted)
 */
//deprecated// #define GEOS_LOG_RANK( msg ) GEOS_LOG_RANK_IF( true, msg )

/**
 * @brief Log a variable/expression name and value on screen to the rank output stream.
 * @param var a variable or expression accessible from current scope that can be stream inserted
 */
#define GEOS_LOG_RANK_VAR( ... ) logger.rankLog( STRINGIZE( __VA_ARGS__ ) " = ", __VA_ARGS__ )//GEOS_LOG_RANK( #var " = " << var )

/**
 * @brief Conditionally raise a hard error and terminate the program.
 * @param EXP an expression that will be evaluated as a predicate
 * @param msg a message to log (any expression that can be stream inserted)
 */
#if defined(GEOS_DEVICE_COMPILE)
#define GEOS_ERROR_IF( EXP, msg ) LVARRAY_ERROR_IF( EXP, msg )
#else
#define GEOS_ERROR_IF( EXP, msg ) LVARRAY_ERROR_IF( EXP, "***** " << ::geos::logger.rankMsgPrefix << msg )
#endif

/**
 * @brief Conditionally throw an exception.
 * @param EXP an expression that will be evaluated as a predicate
 * @param msg a message to log (any expression that can be stream inserted)
 * @param TYPE the type of exception to throw
 */
#define GEOS_THROW_IF( EXP, msg, TYPE ) LVARRAY_THROW_IF( EXP, "***** " << ::geos::logger.rankMsgPrefix << msg, TYPE )

/**
 * @brief Raise a hard error and terminate the program.
 * @param msg a message to log (any expression that can be stream inserted)
 */
#define GEOS_ERROR( msg ) GEOS_ERROR_IF( true, msg )

/**
 * @brief Throw an exception.
 * @param msg a message to log (any expression that can be stream inserted)
 * @param TYPE the type of exception to throw
 */
#define GEOS_THROW( msg, TYPE ) GEOS_THROW_IF( true, msg, TYPE )

/**
 * @brief Assert a condition in debug builds.
 * @param EXP an expression that will be evaluated as a predicate
 * @param msg a message to log (any expression that can be stream inserted)
 */
#define GEOS_ASSERT_MSG( EXP, msg ) LVARRAY_ASSERT_MSG( EXP, "***** " << ::geos::logger.rankMsgPrefix << msg )

/**
 * @brief Assert a condition in debug builds.
 * @param EXP an expression that will be evaluated as a predicate
 */
#define GEOS_ASSERT( EXP ) GEOS_ASSERT_MSG( EXP, "" )

/**
 * @brief Conditionally report a warning.
 * @param EXP an expression that will be evaluated as a predicate
 * @param msg a message to log (any expression that can be stream inserted)
 */
#define GEOS_WARNING_IF( EXP, msg ) LVARRAY_WARNING_IF( EXP, msg )

/**
 * @brief Report a warning.
 * @param msg a message to log (any expression that can be stream inserted)
 */
#define GEOS_WARNING( msg ) LVARRAY_WARNING( msg )

/**
 * @brief Conditionally log an info message.
 * @param EXP an expression that will be evaluated as a predicate
 * @param msg a message to log (any expression that can be stream inserted)
 */
//deprecated// #define GEOS_INFO_IF( EXP, msg ) LVARRAY_INFO_IF( EXP, msg )

/**
 * @brief Log an info message.
 * @param msg a message to log (any expression that can be stream inserted)
 */
//deprecated// #define GEOS_INFO( msg ) LVARRAY_INFO( msg )

/**
 * @brief Raise a hard error if two values are equal.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param msg a message to log (any expression that can be stream inserted)
 */
#define GEOS_ERROR_IF_EQ_MSG( lhs, rhs, msg ) LVARRAY_ERROR_IF_EQ_MSG( lhs, rhs, "***** " << ::geos::logger.rankMsgPrefix << msg )

/**
 * @brief Raise a hard error if two values are equal.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param msg a message to log (any expression that can be stream inserted)
 * @param TYPE the type of exception to throw
 */
#define GEOS_THROW_IF_EQ_MSG( lhs, rhs, msg, TYPE ) LVARRAY_THROW_IF_EQ_MSG( lhs, rhs, "***** " << ::geos::logger.rankMsgPrefix << msg, TYPE )

/**
 * @brief Raise a hard error if two values are equal.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 */
#define GEOS_ERROR_IF_EQ( lhs, rhs ) GEOS_ERROR_IF_EQ_MSG( lhs, rhs, "" )

/**
 * @brief Raise a hard error if two values are equal.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param TYPE the type of exception to throw
 */
#define GEOS_THROW_IF_EQ( lhs, rhs, TYPE ) GEOS_THROW_IF_EQ_MSG( lhs, rhs, "", TYPE )

/**
 * @brief Raise a hard error if two values are not equal.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param msg a message to log (any expression that can be stream inserted)
 */
#define GEOS_ERROR_IF_NE_MSG( lhs, rhs, msg ) LVARRAY_ERROR_IF_NE_MSG( lhs, rhs, "***** " << ::geos::logger.rankMsgPrefix << msg )

/**
 * @brief Throw an exception if two values are not equal.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param msg a message to log (any expression that can be stream inserted)
 * @param TYPE the type of exception to throw
 */
#define GEOS_THROW_IF_NE_MSG( lhs, rhs, msg, TYPE ) LVARRAY_THROW_IF_NE_MSG( lhs, rhs, "***** " << ::geos::logger.rankMsgPrefix << msg, TYPE )

/**
 * @brief Raise a hard error if two values are not equal.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 */
#define GEOS_ERROR_IF_NE( lhs, rhs ) GEOS_ERROR_IF_NE_MSG( lhs, rhs, "" )

/**
 * @brief Throw an exception if two values are not equal.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param TYPE the type of exception to throw
 */
#define GEOS_THROW_IF_NE( lhs, rhs, TYPE ) GEOS_THROW_IF_NE_MSG( lhs, rhs, "", TYPE )

/**
 * @brief Raise a hard error if one value compares greater than the other.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param msg a message to log (any expression that can be stream inserted)
 */
#define GEOS_ERROR_IF_GT_MSG( lhs, rhs, msg ) LVARRAY_ERROR_IF_GT_MSG( lhs, rhs, "***** " << ::geos::logger.rankMsgPrefix << msg )

/**
 * @brief Throw an exception if one value compares greater than the other.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param msg a message to log (any expression that can be stream inserted)
 * @param TYPE the type of exception to throw
 */
#define GEOS_THROW_IF_GT_MSG( lhs, rhs, msg, TYPE ) LVARRAY_THROW_IF_GT_MSG( lhs, rhs, "***** " << ::geos::logger.rankMsgPrefix << msg, TYPE )

/**
 * @brief Raise a hard error if one value compares greater than the other.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 */
#define GEOS_ERROR_IF_GT( lhs, rhs ) GEOS_ERROR_IF_GT_MSG( lhs, rhs, "" )

/**
 * @brief Throw an exception if one value compares greater than the other.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param TYPE the type of exception to throw
 */
#define GEOS_THROW_IF_GT( lhs, rhs, TYPE ) GEOS_ERROR_IF_GT_MSG( lhs, rhs, "", TYPE )

/**
 * @brief Raise a hard error if one value compares greater than or equal to the other.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param msg a message to log (any expression that can be stream inserted)
 */
#define GEOS_ERROR_IF_GE_MSG( lhs, rhs, msg ) LVARRAY_ERROR_IF_GE_MSG( lhs, rhs, "***** " << ::geos::logger.rankMsgPrefix << msg )

/**
 * @brief Throw an exception if one value compares greater than or equal to the other.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param msg a message to log (any expression that can be stream inserted)
 * @param TYPE the type of exception to throw
 */
#define GEOS_THROW_IF_GE_MSG( lhs, rhs, msg, TYPE ) LVARRAY_THROW_IF_GE_MSG( lhs, rhs, "***** " << ::geos::logger.rankMsgPrefix << msg, TYPE )

/**
 * @brief Raise a hard error if one value compares greater than or equal to the other.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 */
#define GEOS_ERROR_IF_GE( lhs, rhs ) GEOS_ERROR_IF_GE_MSG( lhs, rhs, "" )

/**
 * @brief Throw an exception if one value compares greater than or equal to the other.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param TYPE the type of exception to throw
 */
#define GEOS_THROW_IF_GE( lhs, rhs, TYPE ) GEOS_ERROR_IF_GE_MSG( lhs, rhs, "", TYPE )

/**
 * @brief Raise a hard error if one value compares less than the other.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param msg a message to log (any expression that can be stream inserted)
 */
#define GEOS_ERROR_IF_LT_MSG( lhs, rhs, msg ) LVARRAY_ERROR_IF_LT_MSG( lhs, rhs, "***** " << ::geos::logger.rankMsgPrefix << msg )

/**
 * @brief Throw an exception if one value compares less than the other.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param msg a message to log (any expression that can be stream inserted)
 * @param TYPE the type of exception to throw
 */
#define GEOS_THROW_IF_LT_MSG( lhs, rhs, msg, TYPE ) LVARRAY_THROW_IF_LT_MSG( lhs, rhs, "***** " << ::geos::logger.rankMsgPrefix << msg, TYPE )

/**
 * @brief Raise a hard error if one value compares less than the other.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 */
#define GEOS_ERROR_IF_LT( lhs, rhs ) GEOS_ERROR_IF_LT_MSG( lhs, rhs, "" )

/**
 * @brief Throw an exception if one value compares less than the other.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param TYPE the type of exception to throw
 */
#define GEOS_THROW_IF_LT( lhs, rhs, TYPE ) GEOS_ERROR_IF_LT_MSG( lhs, rhs, "", TYPE )

/**
 * @brief Raise a hard error if one value compares less than or equal to the other.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param msg a message to log (any expression that can be stream inserted)
 */
#define GEOS_ERROR_IF_LE_MSG( lhs, rhs, msg ) LVARRAY_ERROR_IF_LE_MSG( lhs, rhs, "***** " << ::geos::logger.rankMsgPrefix << msg )

/**
 * @brief Throw an exception if one value compares less than or equal to the other.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param msg a message to log (any expression that can be stream inserted)
 * @param TYPE the type of exception to throw
 */
#define GEOS_THROW_IF_LE_MSG( lhs, rhs, msg, TYPE ) LVARRAY_THROW_IF_LE_MSG( lhs, rhs, "***** " << ::geos::logger.rankMsgPrefix << msg, TYPE )

/**
 * @brief Raise a hard error if one value compares less than or equal to the other.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 */
#define GEOS_ERROR_IF_LE( lhs, rhs ) GEOS_ERROR_IF_LE_MSG( lhs, rhs, "" )

/**
 * @brief Throw an exception if one value compares less than or equal to the other.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param TYPE the type of exception to throw
 */
#define GEOS_THROW_IF_LE( lhs, rhs, TYPE ) GEOS_ERROR_IF_LE_MSG( lhs, rhs, "", TYPE )

/**
 * @brief Assert that two values compare equal in debug builds.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param msg a message to log (any expression that can be stream inserted)
 */
#define GEOS_ASSERT_EQ_MSG( lhs, rhs, msg ) LVARRAY_ASSERT_EQ_MSG( lhs, rhs, "***** " << ::geos::logger.rankMsgPrefix << msg )

/**
 * @brief Assert that two values compare equal in debug builds.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 */
#define GEOS_ASSERT_EQ( lhs, rhs ) GEOS_ASSERT_EQ_MSG( lhs, rhs, "" )

/**
 * @brief Assert that two values compare not equal in debug builds.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param msg a message to log (any expression that can be stream inserted)
 */
#define GEOS_ASSERT_NE_MSG( lhs, rhs, msg ) LVARRAY_ASSERT_NE_MSG( lhs, rhs, msg )

/**
 * @brief Assert that two values compare not equal in debug builds.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 */
#define GEOS_ASSERT_NE( lhs, rhs ) LVARRAY_ASSERT_NE( lhs, rhs )

/**
 * @brief Assert that one value compares greater than the other in debug builds.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param msg a message to log (any expression that can be stream inserted)
 */
#define GEOS_ASSERT_GT_MSG( lhs, rhs, msg ) LVARRAY_ASSERT_GT_MSG( lhs, rhs, "***** " << ::geos::logger.rankMsgPrefix << msg )

/**
 * @brief Assert that one value compares greater than the other in debug builds.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 */
#define GEOS_ASSERT_GT( lhs, rhs ) GEOS_ASSERT_GT_MSG( lhs, rhs, "" )

/**
 * @brief Assert that one value compares greater than or equal to the other in debug builds.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param msg a message to log (any expression that can be stream inserted)
 */
#define GEOS_ASSERT_GE_MSG( lhs, rhs, msg ) LVARRAY_ASSERT_GE_MSG( lhs, rhs, "***** " << ::geos::logger.rankMsgPrefix << msg )

/**
 * @brief Assert that one value compares greater than or equal to the other in debug builds.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 */
#define GEOS_ASSERT_GE( lhs, rhs ) GEOS_ASSERT_GE_MSG( lhs, rhs, "" )

/**
 * @brief Macro used to turn on/off a function based on the log level.
 * @param[in] minLevel Minimum log level
 * @param[in] fn Function to filter
 */
// //deprecated// #define GEOS_LOG_LEVEL_FN( minLevel, fn )                                      \/
//   do {                                                                         \/
//     if( this->getLogLevel() >= minLevel )                                      \/
//     {                                                                          \/
//       fn;                                                                      \/
//     }                                                                          \/
//   } while( false )

/**
 * @brief Output messages based on current Group's log level.
 * @param[in] minLevel minimum log level
 * @param[in] msg a message to log (any expression that can be stream inserted)
 */
//deprecated// #define GEOS_LOG_LEVEL( minLevel, ... ) logger.stdLogIf( this->getLogLevel() >= minLevel, __VA_ARGS__ )//GEOS_INFO_IF(
// this->getLogLevel()
// >= minLevel, msg );

/**
 * @brief Output messages (only on rank 0) based on current Group's log level.
 * @param[in] minLevel minimum log level
 * @param[in] msg a message to log (any expression that can be stream inserted)
 */
//deprecated// #define GEOS_LOG_LEVEL_RANK_0( minLevel, ... ) logger.rank0LogIf( this->getLogLevel() >= minLevel, __VA_ARGS__
// )//GEOS_LOG_RANK_0_IF(
// this->getLogLevel() >=
// minLevel, msg )

/**
 * @brief Output messages (with one line per rank) based on current Group's log level.
 * @param[in] minLevel minimum log level
 * @param[in] msg a message to log (any expression that can be stream inserted)
 */
//deprecated// #define GEOS_LOG_LEVEL_BY_RANK( minLevel, ... ) logger.rankLogIf( this->getLogLevel() >= minLevel, __VA_ARGS__
// )//GEOS_LOG_RANK_IF(
// this->getLogLevel() >=
// minLevel, msg )


namespace geos
{


/// @brief Enumerate the logging levels from the most importants messages to the one that contain the more details.
enum class LogLevel : int16_t
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


/// @brief Class used to log messages in the standard output and rank log files.
/// A "logger" extern global is available for general use, but other instances can be created where
/// a specialized logger is needed.
/// Doesn't have any GPU logging capacities for now.
class Logger
{
public:

  /// @brief The default message log level when logging directly with the Logger
  static constexpr LogLevel defaultMsgLogLevel = LogLevel::Debug;
  // TODO :
  static constexpr LogLevel defaultProblemLogLevel = LogLevel::Detailed;// TODO Logger: mettre ça à  Progress? Detailed? faire en f° du tri des messages
  // TODO :
  static constexpr LogLevel defaultGroupLogLevel = LogLevel::Detailed;// TODO Logger: mettre ça à  Progress? Detailed? faire en f° du tri des messages

  /**
   * @brief Construct the logger. Default output is the standard one.
   */
  Logger();
  /**
   * @brief Re-initialize the logger configuration. Close the potential file streams.
   */
  void reset();
  /**
   * @brief Initialize the logger in a parallel build. Must be called at least once, after the
   * MPI_Init() call. Must not be called in a serial build.
   * @param mpiComm global MPI communicator
   */
  void initMpi( MPI_Comm mpiComm );

  /**
   * @brief Change the logger output to the specified stream.
   * Also close the potential previously opened file output stream.
   * @param stream the stream to use when calling rankLog()
   */
  void setRankOutputToStream( std::ostream & stream );
  /**
   * @brief Change the logger output to a automatically created file in the specified folder.
   * The file name is based on the rank number. The file extension is a ".out".
   * @param output_dir optional output directory for rank log files to stream to when calling rankLog()
   */
  void setRankOutputToFile( const std::string & output_dir = "" );

  /**
   * @param param set the log level for global messages (ie. applied to simulation progress messages).
   */
  void setProblemLogLevel( LogLevel param );
  /**
   * @param param value to override all groups LogLevel and the problem LogLevel.
   */
  void enableLogLevelOverride( LogLevel param );
  /**
   * @param param set the maximal log level allowed for all messages. With a value of 2, all
   * messages from level 3 and more will not be output. With a value of -1, the logger become silent.
   * This value is prioritized over setMinLogLevel.
   */
  void disableLogLevelOverride();


  /**
   * @brief log one or more inputs, only from the rank 0, to the standard output (and to
   * the rank 0 file stream if used).
   * If ProblemLogLevel is at least set to Detailed, a prefix est added to know which rank is streaming
   * the message. As an exemple: "Rank 0: Hello World!"
   * @param MSG_LEVEL the level of the message to log: The message will be ignored if the MSG_LEVEL
   * is strictly higher than the current globalLogLevel value.
   * @param inputs the inputs to log.
   * @tparam INPUTS types of the inputs.
   */
  template< LogLevel MSG_LEVEL = defaultMsgLogLevel, typename ... INPUTS >
  void logRank0( INPUTS ... inputs );

  template< LogLevel MSG_LEVEL = defaultMsgLogLevel, typename ... INPUTS >
  void logRank0If( bool cond, INPUTS ... inputs );


  /**
   * @brief log one or more inputs to the rank file stream if used, or to the standard
   * output (the message is streamed repeatedly if it comes from multiple ranks).
   * If ProblemLogLevel is at least set to Detailed, a prefix est added to know which rank is streaming
   * the message. As an exemple: "Rank 54: Hello World!"
   * @param MSG_LEVEL the level of the message to log: The message will be ignored if the MSG_LEVEL
   * is strictly higher than the current globalLogLevel value.
   * @param input the inputs to log.
   * @tparam INPUTS types of the inputs.
   */
  template< LogLevel MSG_LEVEL = defaultMsgLogLevel, typename ... INPUTS >
  void log( INPUTS ... inputs );

  template< LogLevel MSG_LEVEL = defaultMsgLogLevel, typename ... INPUTS >
  void logIf( bool cond, INPUTS ... inputs );


  // /**
  //  * @brief Log a message to the rank file stream if used, or log it to the std::cout only once per
  //  * form (and therefore does a MPI barrier).
  //  * If ProblemLogLevel is at least set to Detailed, a prefix est added to know which rank is streaming
  //  * the message. As an exemple, for a message that can be "Hello World!" or "Hello Folks!":
  //  * Rank 0->54, 56, 60, 67->127: Hello World!
  //  * Rank 55, 57->59, 61->66: Hello Folks!
  //  * @param MSG_LEVEL the level of the message to log: The message will be ignored if the MSG_LEVEL
  //  * is strictly higher than the current globalLogLevel value.
  //  * @param input the inputs to log.
  //  * @tparam INPUTS types of the inputs.
  //  */
  // template< LogLevel MSG_LEVEL = defaultMsgLogLevel, typename ... INPUTS >
  // void mergedLog( INPUTS ... inputs );
  //
  // template< LogLevel MSG_LEVEL = defaultMsgLogLevel, typename ... INPUTS >
  // void mergedLogIf( bool cond, INPUTS ... inputs );


  //TODO logDevice()


public://todo: private

  /// @see setProblemLogLevel( LogLevel param )
  LogLevel m_problemLogLevel;
  /// @see enableLogLevelOverride( LogLevel param ) and disableLogLevelOverride()
  std::optional< LogLevel > m_overrideLogLevel;

  /// @brief MPI communicator
  MPI_Comm m_comm;
  /// @brief MPI rank id
  int m_rank;
  /// @brief MPI total ranks count
  int m_ranksCount;
  /// @brief prefix to add before a message that is specific to a rank. Empty in a serial build.
  std::string m_rankMsgPrefix;

  /// @brief Pointer to the rank output stream
  std::ostream * m_outStream;
  /// @brief Smart pointer to the active file output stream. Equals to nullptr if not outputing to a file.
  std::unique_ptr< std::ostream > m_fileOutStream;

  // todo comment
  template< typename INPUT >
  static void streamLog( std::ostream & stream, INPUT input );
  // todo comment
  template< typename INPUT, typename ... MORE_INPUTS >
  static void streamLog( std::ostream & stream, INPUT input, MORE_INPUTS ... moreInputs );

};
extern Logger logger;


struct LogSource
{

  LogLevel m_logLevel;


  LogSource():
    m_logLevel( Logger::defaultGroupLogLevel )
  {}
  LogSource( LogLevel logLevel ):
    m_logLevel( logLevel )
  {}


  //TODO implementations
  template< LogLevel MSG_LEVEL = defaultMsgLogLevel, typename ... INPUTS >
  void log( INPUTS ... inputs );

  template< LogLevel MSG_LEVEL = defaultMsgLogLevel, typename ... INPUTS >
  void logIf( bool cond, INPUTS ... inputs );

  template< LogLevel MSG_LEVEL = defaultMsgLogLevel, typename ... INPUTS >
  void mergedLog( INPUTS ... inputs );

  template< LogLevel MSG_LEVEL = defaultMsgLogLevel, typename ... INPUTS >
  void mergedLogIf( bool cond, INPUTS ... inputs );

  template< LogLevel MSG_LEVEL = defaultMsgLogLevel, typename ... INPUTS >
  void rank0Log( INPUTS ... inputs );

  template< LogLevel MSG_LEVEL = defaultMsgLogLevel, typename ... INPUTS >
  void rank0LogIf( bool cond, INPUTS ... inputs );

};


/**
 * @brief Exception class used to report errors in user input.
 */
struct InputError : public std::runtime_error
{
  /**
   * @brief Constructor
   * @param what the error message
   */
  InputError( std::string const & what ):
    std::runtime_error( what )
  {}

  /**
   * @brief Constructor
   * @param what the error message
   */
  InputError( char const * const what ):
    std::runtime_error( what )
  {}

  /**
   * @brief Construct an InputError from an underlying exception.
   * @param subException An exception to base this new one on.
   * @param msgToInsert The error message. It will be inserted into the one inside of subException.
   */
  InputError( std::exception const & subException, std::string const & msgToInsert );
};

/**
 * @brief Exception class used for special control flow.
 */
class NotAnError : public std::exception
{};


template< typename INPUT >
void Logger::streamLog( std::ostream & stream, INPUT input )
{
  stream << input << std::endl;
}

template< typename INPUT, typename ... MORE_INPUTS >
void Logger::streamLog( std::ostream & stream, INPUT input, MORE_INPUTS ... moreInputs )
{
  stream << input;
  streamLog( stream, moreInputs ... );
}


//WIP section...
template< LogLevel MSG_LEVEL, typename ... INPUTS >
void Logger::rankLog( INPUTS ... inputs )
{
  if( MSG_LEVEL <= m_problemLogLevel )
  {
    streamLog( *m_outStream, m_rankMsgPrefix, inputs ... );
  }
}

template< LogLevel MSG_LEVEL, typename ... INPUTS >
void Logger::log( INPUTS ... inputs )
{
  if( MSG_LEVEL <= m_problemLogLevel )
  {
    if( m_outStream != &std::cout )
    {
      streamLog( *m_outStream, m_rankMsgPrefix, inputs ... );
    }
    else
    {
      std::ostringstream oss;
      streamLog( *m_outStream, m_rankMsgPrefix, inputs ... );

    }
  }
}

template< LogLevel MSG_LEVEL, typename ... INPUTS >
void Logger::rank0Log( INPUTS ... inputs )
{
  if( m_rank == 0 && MSG_LEVEL <= m_problemLogLevel )
  {
    streamLog( *m_outStream, m_rankMsgPrefix, inputs ... );
  }
}


template< LogLevel MSG_LEVEL, typename ... INPUTS >
void Logger::stdLogIf( bool condition, INPUTS ... inputs )
{
  if( condition )
  {
    stdLog< MSG_LEVEL >( inputs ... );
  }
}

template< LogLevel MSG_LEVEL, typename ... INPUTS >
void Logger::mergedLogIf( bool condition, INPUTS ... inputs )
{
  if( condition )
  {
    rankLog< MSG_LEVEL >( inputs ... );
  }
}

template< LogLevel MSG_LEVEL, typename ... INPUTS >
void Logger::rank0LogIf( bool condition, INPUTS ... inputs )
{
  if( condition )
  {
    rank0Log< MSG_LEVEL >( inputs ... );
  }
}

// TODO Logger: decide to use this version or not (streams on std output + rank 0 file stream)
// template< typename ... INPUTS >
// void Logger::rank0log( INPUTS ... inputs )
// {
//   if( rank == 0 )
//   {
//     if( fileOutStream!=nullptr )
//     {
//       streamLog( &std::cout, rankMsgPrefix, inputs ... );
//     }
//     else
//     {
//       std::ostringstream oss;
//       streamLog( oss, rankMsgPrefix, inputs ... );
//       string const str( oss.str() );
//       std::cout << str;
//       fileOutStream << str;
//     }
//   }
// }

} // namespace geos

#endif /* GEOS_COMMON_LOGGER_HPP */
