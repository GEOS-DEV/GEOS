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
#define GEOS_LOG( ... ) LVARRAY_LOG( __VA_ARGS__ )

/**
 * @brief Log an expression and its value on screen.
 * @details The expression to log must evaluate something that can be stream inserted.
 */
#define GEOS_LOG_VAR( ... ) LVARRAY_LOG_VAR( __VA_ARGS__ )

/**
 * @brief Conditionally log a message on screen on rank 0.
 * @param EXP an expression that will be evaluated as a predicate
 * @param msg a message to log (any expression that can be stream inserted)
 */
#define GEOS_LOG_RANK_0_IF( EXP, msg ) \
  do { \
    if( ::geos::logger.rank == 0 && EXP ) \
    { \
      std::ostringstream oss; \
      oss << msg; \
      std::cout << oss.str() << std::endl; \
    } \
  } while( false )

/**
 * @brief Log a message on screen on rank 0.
 * @param msg a message to log (any expression that can be stream inserted)
 */
#define GEOS_LOG_RANK_0( msg ) GEOS_LOG_RANK_0_IF( true, msg )

/**
 * @brief Conditionally log a message to the rank output stream.
 * @param EXP an expression that will be evaluated as a predicate
 * @param msg a message to log (any expression that can be stream inserted)
 */
#if defined(GEOS_DEVICE_COMPILE)
#define GEOS_LOG_RANK_IF( EXP, msg )
#else
#define GEOS_LOG_RANK_IF( EXP, msg ) \
  do { \
    if( EXP ) \
    { \
      std::ostringstream oss; \
      oss << ::geos::logger.rankMsgPrefix << msg; \
      *logger.outStream << oss.str() << std::endl; \
    } \
  } while( false )
#endif

/**
 * @brief Log a message to the rank output stream.
 * @param msg a message to log (any expression that can be stream inserted)
 */
#define GEOS_LOG_RANK( msg ) GEOS_LOG_RANK_IF( true, msg )

/**
 * @brief Log a variable/expression name and value on screen to the rank output stream.
 * @param var a variable or expression accessible from current scope that can be stream inserted
 */
#define GEOS_LOG_RANK_VAR( var ) GEOS_LOG_RANK( #var " = " << var )

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
#define GEOS_INFO_IF( EXP, msg ) LVARRAY_INFO_IF( EXP, msg )

/**
 * @brief Log an info message.
 * @param msg a message to log (any expression that can be stream inserted)
 */
#define GEOS_INFO( msg ) LVARRAY_INFO( msg )

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
#define GEOS_LOG_LEVEL_FN( minLevel, fn )                                      \
  do {                                                                         \
    if( this->getLogLevel() >= minLevel )                                      \
    {                                                                          \
      fn;                                                                      \
    }                                                                          \
  } while( false )

/**
 * @brief Output messages based on current Group's log level.
 * @param[in] minLevel minimum log level
 * @param[in] msg a message to log (any expression that can be stream inserted)
 */
#define GEOS_LOG_LEVEL( minLevel, msg ) GEOS_INFO_IF( this->getLogLevel() >= minLevel, msg );

/**
 * @brief Output messages (only on rank 0) based on current Group's log level.
 * @param[in] minLevel minimum log level
 * @param[in] msg a message to log (any expression that can be stream inserted)
 */
#define GEOS_LOG_LEVEL_RANK_0( minLevel, msg ) GEOS_LOG_RANK_0_IF( this->getLogLevel() >= minLevel, msg )

/**
 * @brief Output messages (with one line per rank) based on current Group's log level.
 * @param[in] minLevel minimum log level
 * @param[in] msg a message to log (any expression that can be stream inserted)
 */
#define GEOS_LOG_LEVEL_BY_RANK( minLevel, msg ) GEOS_LOG_RANK_IF( this->getLogLevel() >= minLevel, msg )


namespace geos
{


class Logger
{
public:

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
   * @brief Change the logger output to a automatically created file in the specified folder.
   * The file name is based on the rank number. The file extension is a ".out".
   * @param output_dir optional output directory for rank log files
   */
  void setOutputToFile( const std::string & output_dir = "" );
  /**
   * @brief Change the logger output to the standard output (std::cout).
   * Also close the potential file output stream.
   */
  void setOutputToStd();

  /**
   * @param value set the log level for global messages (ie. applied to simulation progress messages).
   */
  void setGlobalLogLevel( int value );
  /**
   * @param value set the minimal log level requiered for all messages. With a value of 2, all
   * messages from level 0, 1 and 2 will be output.
   */
  void setMinLogLevel( int value );
  /**
   * @param value set the maximal log level allowed for all messages. With a value of 2, all
   * messages from level 3 and more will not be output.
   */
  void setMaxLogLevel( int value );

public://todo: private
  /// @see setGlobalLogLevel( int value )
  int globalLogLevel;
  /// @see setMinLogLevel( int value )
  int minLogLevel;
  /// @see setMaxLogLevel( int value )
  int maxLogLevel;

  /// @brief MPI communicator
  MPI_Comm comm;
  /// @brief MPI rank id
  int rank;
  /// @brief MPI total ranks count
  int ranksCount;
  /// @brief prefix to add before a message that is specific to a rank. Empty in a serial build.
  std::string rankMsgPrefix;

  /// @brief Pointer to the selected output stream
  std::ostream * outStream;
  /// @brief Smart pointer to the active file output stream. Equals to nullptr if not outputing to a file.
  std::unique_ptr< std::ostream > fileOutStream;
};
extern Logger logger;


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


} // namespace geos

#endif /* GEOS_COMMON_LOGGER_HPP */
