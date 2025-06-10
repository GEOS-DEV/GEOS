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
 * @file Logger.hpp
 */

#ifndef GEOS_COMMON_LOGGER_HPP
#define GEOS_COMMON_LOGGER_HPP

// Source incldes
#include "common/GeosxConfig.hpp"
#include "common/GeosxMacros.hpp"
#include "common/format/Format.hpp"
#include "LvArray/src/Macros.hpp"
#include "common/logger/ErrorHandling.hpp"

// System includes
#include <stdexcept>

#if defined(GEOS_USE_MPI)
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
 * @brief Conditionally log a message.
 * @param EXP an expression that will be evaluated as a predicate
 * @param msg a message to log (any expression that can be stream inserted)
 */
#if defined(GEOS_DEVICE_COMPILE)
#define GEOS_LOG_IF( EXP, msg )
#else
#define GEOS_LOG_IF( EXP, msg ) \
  do { \
    if( EXP ) \
    { \
      std::cout<< msg << std::endl; \
    } \
  } while( false )
#endif


/**
 * @brief Conditionally log a message on screen on rank 0.
 * @param EXP an expression that will be evaluated as a predicate
 * @param msg a message to log (any expression that can be stream inserted)
 */
#define GEOS_LOG_RANK_0_IF( EXP, msg ) \
  do { \
    if( ::geos::logger::internal::rank == 0 && EXP ) \
    { \
      std::ostringstream oss; \
      oss << msg; \
      std::cout << oss.str() << std::endl; \
    } \
  } while( false )

/**
 * @brief Conditionally log a message on screen on rank 0 without line breaking.
 * @param EXP an expression that will be evaluated as a predicate
 * @param msg a message to log (any expression that can be stream inserted)
 */
#define GEOS_LOG_RANK_0_IF_NLR( EXP, msg ) \
  do { \
    if( ::geos::logger::internal::rank == 0 && EXP ) \
    { \
      std::ostringstream oss; \
      oss << msg; \
      std::cout << oss.str(); \
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
      oss << "Rank " << ::geos::logger::internal::rankString << ": " << msg; \
      *logger::internal::rankStream << oss.str() << std::endl; \
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
 * @param MSG a message to log (any expression that can be stream inserted)
 */
#define GEOS_ERROR_IF_IMPL( EXP, MSG ) \
  do \
  { \
    if( EXP ) \
    { \
      std::ostringstream __oss; \
      __oss << "***** ERROR\n"; \
      __oss << "***** LOCATION: " LOCATION "\n"; \
      __oss << "***** Controlling expression (should be false): " STRINGIZE( EXP ) "\n"; \
      __oss << MSG << "\n"; \
      std::string stackHistory = LvArray::system::stackTrace( true ); \
      __oss << stackHistory; \
      std::cout << __oss.str() << std::endl; \
      std::ostringstream __msgoss; \
      __msgoss << MSG; \
      if( errorLogger.writeFile() ) \
      { \
        ErrorLogger::ErrorMsg msgStruct( ErrorLogger::MsgType::Error, \
                                         __msgoss.str(), \
                                         __FILE__, \
                                         __LINE__ ); \
        msgStruct.setRank( ::geos::logger::internal::rank ); \
        msgStruct.addCallStackInfo( stackHistory ); \
        errorLogger.write( msgStruct ); \
      } \
      LvArray::system::callErrorHandler(); \
    } \
  } while( false )

/**
 * @brief Conditionally raise a hard error and terminate the program.
 * @param EXP an expression that will be evaluated as a predicate
 * @param MSG a message to log (any expression that can be stream inserted)
 * @param ... One or more DataContext (current error context information) 
 */
#define GEOS_ERROR_CTX_IF( EXP, MSG, ... ) \
  do \
  { \
    if( EXP ) \
    { \
      std::ostringstream __oss; \
      __oss << "***** ERROR\n"; \
      __oss << "***** LOCATION: " LOCATION "\n"; \
      __oss << "***** Controlling expression (should be false): " STRINGIZE( EXP ) "\n"; \
      __oss << MSG << "\n"; \
      std::string stackHistory = LvArray::system::stackTrace( true ); \
      __oss << stackHistory; \
      std::cout << __oss.str() << std::endl; \
      std::ostringstream __msgoss; \
      __msgoss << MSG; \
      if( errorLogger.writeFile() ) \
      { \
        ErrorLogger::ErrorMsg msgStruct( ErrorLogger::MsgType::Error, \
                                         __msgoss.str(), \
                                         __FILE__, \
                                         __LINE__ ); \
        msgStruct.setRank( ::geos::logger::internal::rank ); \
        msgStruct.addContextInfo( __VA_ARGS__ ); \
        msgStruct.addCallStackInfo( stackHistory ); \
        errorLogger.write( msgStruct ); \
      } \
      LvArray::system::callErrorHandler(); \
    } \
  } while( false )

/**
 * @brief Conditionally raise a hard error and terminate the program.
 * @param EXP an expression that will be evaluated as a predicate
 * @param msg a message to log (any expression that can be stream inserted)
 */
#if defined(GEOS_DEVICE_COMPILE)
#define GEOS_ERROR_IF( EXP, msg ) GEOS_ERROR_IF_IMPL( EXP, msg )
#else
#define GEOS_ERROR_IF( EXP, msg ) GEOS_ERROR_IF_IMPL( EXP, "***** Rank " << ::geos::logger::internal::rankString << ": " << msg )
#endif

/**
 * @brief Raise a hard error and terminate the program.
 * @param msg a message to log (any expression that can be stream inserted)
 */
#define GEOS_ERROR( msg ) GEOS_ERROR_IF( true, msg )

/**
 * @brief Conditionally throw an exception.
 * @param EXP an expression that will be evaluated as a predicate
 * @param MSG a message to log (any expression that can be stream inserted)
 * @param TYPE the type of exception to throw
 */
#define GEOS_THROW_IF_IMPL( EXP, MSG, EXCEPTIONTYPE ) \
  do \
  { \
    if( EXP ) \
    { \
      std::ostringstream __oss; \
      __oss << "\n"; \
      __oss << "***** LOCATION: " LOCATION "\n"; \
      __oss << "***** Controlling expression (should be false): " STRINGIZE( EXP ) "\n"; \
      __oss << MSG << "\n"; \
      std::string stackHistory = LvArray::system::stackTrace( true ); \
      __oss << stackHistory; \
      std::cout << __oss.str() << std::endl; \
      std::ostringstream __msgoss; \
      __msgoss << MSG; \
      if( errorLogger.writeFile() ) \
      { \
        errorLogger.currentErrorMsg() \
          .setType( ErrorLogger::MsgType::Exception ) \
          .setCodeLocation( __FILE__, __LINE__ ) \
          .addToMsg( __msgoss.str() ) \
          .setRank( ::geos::logger::internal::rank ) \
          .addCallStackInfo( stackHistory ); \
      } \
      throw EXCEPTIONTYPE( __oss.str() ); \
    } \
  } while( false )

/**
 * @brief Conditionally throw an exception.
 * @param EXP an expression that will be evaluated as a predicate
 * @param MSG a message to log (any expression that can be stream inserted)
 * @param TYPE the type of exception to throw
 * @param ... One or more DataContext (current error context information) 
 */
#define GEOS_THROW_CTX_IF( EXP, MSG, EXCEPTIONTYPE, ... ) \
  do \
  { \
    if( EXP ) \
    { \
      std::ostringstream __oss; \
      __oss << "\n"; \
      __oss << "***** LOCATION: " LOCATION "\n"; \
      __oss << "***** Controlling expression (should be false): " STRINGIZE( EXP ) "\n"; \
      __oss << MSG << "\n"; \
      std::string stackHistory = LvArray::system::stackTrace( true ); \
      __oss << stackHistory; \
      std::cout << __oss.str() << std::endl; \
      std::ostringstream __msgoss; \
      __msgoss << MSG; \
      if( errorLogger.writeFile() ) \
      { \
        errorLogger.currentErrorMsg() \
          .setType( ErrorLogger::MsgType::Exception ) \
          .setCodeLocation( __FILE__, __LINE__ ) \
          .addToMsg( __msgoss.str() ) \
          .setRank( ::geos::logger::internal::rank ) \
          .addCallStackInfo( stackHistory ) \
          .addContextInfo( __VA_ARGS__ ); \
      } \
      throw EXCEPTIONTYPE( __oss.str() ); \
    } \
  } while( false )

/**
 * @brief Conditionally throw an exception.
 * @param EXP an expression that will be evaluated as a predicate
 * @param msg a message to log (any expression that can be stream inserted)
 * @param TYPE the type of exception to throw
 */
#define GEOS_THROW_IF( EXP, msg, TYPE ) GEOS_THROW_IF_IMPL( EXP, "***** Rank " << ::geos::logger::internal::rankString << ": " << msg, TYPE )

/**
 * @brief Throw an exception.
 * @param msg a message to log (any expression that can be stream inserted)
 * @param TYPE the type of exception to throw
 */
#define GEOS_THROW( msg, TYPE ) GEOS_THROW_IF( true, msg, TYPE )

/**
 * @brief Conditionally report a warning
 * @param EXP an expression that will be evaluated as a predicate
 * @param MSG a message to log (any expression that can be stream inserted)
 */
#define GEOS_WARNING_IF_IMPL( EXP, MSG ) \
  do \
  { \
    if( EXP ) \
    { \
      std::ostringstream __oss; \
      __oss << "***** WARNING\n"; \
      __oss << "***** LOCATION: " LOCATION "\n"; \
      __oss << "***** Controlling expression (should be false): " STRINGIZE( EXP ) "\n"; \
      __oss << MSG; \
      std::cout << __oss.str() << std::endl; \
      std::ostringstream __msgoss; \
      __msgoss << MSG; \
      if( errorLogger.writeFile() ) \
      { \
        ErrorLogger::ErrorMsg msgStruct( ErrorLogger::MsgType::Warning, \
                                         __msgoss.str(), \
                                         __FILE__, \
                                         __LINE__ ); \
        msgStruct.setRank( ::geos::logger::internal::rank ); \
        msgStruct.addCallStackInfo( LvArray::system::stackTrace( true ) ); \
        errorLogger.write( msgStruct ); \
      } \
    } \
  } while( false )

/**
 * @brief Conditionally report a warning
 * @param EXP an expression that will be evaluated as a predicate
 * @param MSG a message to log (any expression that can be stream inserted)
 * @param ... One or more DataContext (current error context information) 
 */
#define GEOS_WARNING_CTX_IF( EXP, MSG, ... ) \
  do \
  { \
    if( EXP ) \
    { \
      std::ostringstream __oss; \
      __oss << "***** WARNING\n"; \
      __oss << "***** LOCATION: " LOCATION "\n"; \
      __oss << "***** Controlling expression (should be false): " STRINGIZE( EXP ) "\n"; \
      __oss << MSG; \
      std::cout << __oss.str() << std::endl; \
      std::ostringstream __msgoss; \
      __msgoss << MSG; \
      if( errorLogger.writeFile() ) \
      { \
        ErrorLogger::ErrorMsg msgStruct( ErrorLogger::MsgType::Warning, \
                                         __msgoss.str(), \
                                         __FILE__, \
                                         __LINE__ ); \
        msgStruct.setRank( ::geos::logger::internal::rank ); \
        msgStruct.addContextInfo( __VA_ARGS__ ); \
        msgStruct.addCallStackInfo( LvArray::system::stackTrace( true ) ); \
        errorLogger.write( msgStruct ); \
      } \
    } \
  } while( false )

/**
 * @brief Conditionally report a warning.
 * @param EXP an expression that will be evaluated as a predicate
 * @param msg a message to log (any expression that can be stream inserted)
 */
#define GEOS_WARNING_IF( EXP, msg ) GEOS_WARNING_IF_IMPL( EXP, msg )

/**
 * @brief Report a warning.
 * @param msg a message to log (any expression that can be stream inserted)
 */
#define GEOS_WARNING( msg ) GEOS_WARNING_IF( true, msg )

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
 * @brief Abort execution if @p lhs @p OP @p rhs.
 * @param lhs The left side of the operation.
 * @param OP The operation to apply.
 * @param NOP The opposite of @p OP, used in the message.
 * @param rhs The right side of the operation.
 * @param msg The message to diplay.
 */
#define GEOS_ERROR_IF_OP_MSG( lhs, OP, NOP, rhs, msg ) \
  GEOS_ERROR_IF_IMPL( lhs OP rhs, \
                    msg << "\n" << \
                    "Expected " << #lhs << " " << #NOP << " " << #rhs << "\n" << \
                    "  " << #lhs << " = " << lhs << "\n" << \
                    "  " << #rhs << " = " << rhs << "\n" )

/**
 * @brief Raise a hard error if two values are equal.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param msg a message to log (any expression that can be stream inserted)
 */
#define LOG_GEOS_ERROR_IF_EQ_MSG( lhs, rhs, msg ) GEOS_ERROR_IF_OP_MSG( lhs, ==, !=, rhs, msg )

/**
 * @brief Raise a hard error if two values are equal.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param msg a message to log (any expression that can be stream inserted)
 */
#define GEOS_ERROR_IF_EQ_MSG( lhs, rhs, msg ) LOG_GEOS_ERROR_IF_EQ_MSG( lhs, rhs, "***** Rank " << ::geos::logger::internal::rankString << ": " << msg )

/**
 * @brief Throw an exception if @p lhs @p OP @p rhs.
 * @param lhs The left side of the operation.
 * @param OP The operation to apply.
 * @param NOP The opposite of @p OP, used in the message.
 * @param rhs The right side of the operation.
 * @param msg The message to diplay.
 * @param TYPE the type of exception to throw.
 */
#define GEOS_THROW_IF_OP_MSG( lhs, OP, NOP, rhs, msg, TYPE ) \
  GEOS_THROW_IF_IMPL( lhs OP rhs, \
                    msg << "\n" << \
                    "Expected " << #lhs << " " << #NOP << " " << #rhs << "\n" << \
                    "  " << #lhs << " = " << lhs << "\n" << \
                    "  " << #rhs << " = " << rhs << "\n", TYPE )

/**
 * @brief Throw an exception if two values are equal.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param msg a message to log (any expression that can be stream inserted)
 * @param TYPE the type of exception to throw
 */
#define LOG_GEOS_THROW_IF_EQ_MSG( lhs, rhs, msg, TYPE ) GEOS_THROW_IF_OP_MSG( lhs, ==, !=, rhs, msg, TYPE )

/**
 * @brief Raise a hard error if two values are equal.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param msg a message to log (any expression that can be stream inserted)
 * @param TYPE the type of exception to throw
 */
#define GEOS_THROW_IF_EQ_MSG( lhs, rhs, msg, TYPE ) LOG_GEOS_THROW_IF_EQ_MSG( lhs, rhs, "***** Rank " << ::geos::logger::internal::rankString << ": " << msg, TYPE )

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
#define LOG_GEOS_ERROR_IF_NE_MSG( lhs, rhs, msg ) GEOS_ERROR_IF_OP_MSG( lhs, !=, ==, rhs, msg )

/**
 * @brief Raise a hard error if two values are not equal.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param msg a message to log (any expression that can be stream inserted)
 */
#define GEOS_ERROR_IF_NE_MSG( lhs, rhs, msg ) LOG_GEOS_ERROR_IF_NE_MSG( lhs, rhs, "***** Rank " << ::geos::logger::internal::rankString << ": " << msg )

/**
 * @brief Throw an exception if two values are not equal.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param msg a message to log (any expression that can be stream inserted)
 * @param TYPE the type of exception to throw
 */
#define LOG_GEOS_THROW_IF_NE_MSG( lhs, rhs, msg, TYPE ) GEOS_THROW_IF_OP_MSG( lhs, !=, ==, rhs, msg, TYPE )

/**
 * @brief Throw an exception if two values are not equal.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param msg a message to log (any expression that can be stream inserted)
 * @param TYPE the type of exception to throw
 */
#define GEOS_THROW_IF_NE_MSG( lhs, rhs, msg, TYPE ) LOG_GEOS_THROW_IF_NE_MSG( lhs, rhs, "***** Rank " << ::geos::logger::internal::rankString << ": " << msg, TYPE )

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
#define LOG_GEOS_ERROR_IF_GT_MSG( lhs, rhs, msg ) GEOS_ERROR_IF_OP_MSG( lhs, >, <=, rhs, msg )

/**
 * @brief Raise a hard error if one value compares greater than the other.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param msg a message to log (any expression that can be stream inserted)
 */
#define GEOS_ERROR_IF_GT_MSG( lhs, rhs, msg ) LOG_GEOS_ERROR_IF_GT_MSG( lhs, rhs, "***** Rank " << ::geos::logger::internal::rankString << ": " << msg )

/**
 * @brief Throw an exception if one value compares greater than the other.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param msg a message to log (any expression that can be stream inserted)
 * @param TYPE the type of exception to throw
 */
#define LOG_GEOS_THROW_IF_GT_MSG( lhs, rhs, msg, TYPE ) GEOS_THROW_IF_OP_MSG( lhs, >, <=, rhs, msg, TYPE )

/**
 * @brief Throw an exception if one value compares greater than the other.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param msg a message to log (any expression that can be stream inserted)
 * @param TYPE the type of exception to throw
 */
#define GEOS_THROW_IF_GT_MSG( lhs, rhs, msg, TYPE ) LOG_GEOS_THROW_IF_GT_MSG( lhs, rhs, "***** Rank " << ::geos::logger::internal::rankString << ": " << msg, TYPE )

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
#define LOG_GEOS_ERROR_IF_GE_MSG( lhs, rhs, msg ) GEOS_ERROR_IF_OP_MSG( lhs, >=, <, rhs, msg )


/**
 * @brief Raise a hard error if one value compares greater than or equal to the other.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param msg a message to log (any expression that can be stream inserted)
 */
#define GEOS_ERROR_IF_GE_MSG( lhs, rhs, msg ) LOG_GEOS_ERROR_IF_GE_MSG( lhs, rhs, "***** Rank " << ::geos::logger::internal::rankString << ": " << msg )

/**
 * @brief Throw an exception if one value compares greater than or equal to the other.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param msg a message to log (any expression that can be stream inserted)
 * @param TYPE the type of exception to throw
 */
#define LOG_GEOS_THROW_IF_GE_MSG( lhs, rhs, msg, TYPE ) GEOS_THROW_IF_OP_MSG( lhs, >=, <, rhs, msg, TYPE )

/**
 * @brief Throw an exception if one value compares greater than or equal to the other.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param msg a message to log (any expression that can be stream inserted)
 * @param TYPE the type of exception to throw
 */
#define GEOS_THROW_IF_GE_MSG( lhs, rhs, msg, TYPE ) LOG_GEOS_THROW_IF_GE_MSG( lhs, rhs, "***** Rank " << ::geos::logger::internal::rankString << ": " << msg, TYPE )

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
#define LOG_GEOS_ERROR_IF_LT_MSG( lhs, rhs, msg ) GEOS_ERROR_IF_OP_MSG( lhs, <, >=, rhs, msg )

/**
 * @brief Raise a hard error if one value compares less than the other.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param msg a message to log (any expression that can be stream inserted)
 */
#define GEOS_ERROR_IF_LT_MSG( lhs, rhs, msg ) LOG_GEOS_ERROR_IF_LT_MSG( lhs, rhs, "***** Rank " << ::geos::logger::internal::rankString << ": " << msg )

/**
 * @brief Throw an exception if one value compares less than the other.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param msg a message to log (any expression that can be stream inserted)
 * @param TYPE the type of exception to throw
 */
#define LOG_GEOS_THROW_IF_LT_MSG( lhs, rhs, msg, TYPE ) GEOS_THROW_IF_OP_MSG( lhs, <, >=, rhs, msg, TYPE )

/**
 * @brief Throw an exception if one value compares less than the other.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param msg a message to log (any expression that can be stream inserted)
 * @param TYPE the type of exception to throw
 */
#define GEOS_THROW_IF_LT_MSG( lhs, rhs, msg, TYPE ) LOG_GEOS_THROW_IF_LT_MSG( lhs, rhs, "***** Rank " << ::geos::logger::internal::rankString << ": " << msg, TYPE )

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
#define LOG_GEOS_ERROR_IF_LE_MSG( lhs, rhs, msg ) GEOS_ERROR_IF_OP_MSG( lhs, <=, >, rhs, msg )

/**
 * @brief Raise a hard error if one value compares less than or equal to the other.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param msg a message to log (any expression that can be stream inserted)
 */
#define GEOS_ERROR_IF_LE_MSG( lhs, rhs, msg ) LOG_GEOS_ERROR_IF_LE_MSG( lhs, rhs, "***** Rank " << ::geos::logger::internal::rankString << ": " << msg )

/**
 * @brief Throw an exception if one value compares less than or equal to the other.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param msg a message to log (any expression that can be stream inserted)
 * @param TYPE the type of exception to throw
 */
#define LOG_GEOS_THROW_IF_LE_MSG( lhs, rhs, msg, TYPE ) GEOS_THROW_IF_OP_MSG( lhs, <=, >, rhs, msg, TYPE )

/**
 * @brief Throw an exception if one value compares less than or equal to the other.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param msg a message to log (any expression that can be stream inserted)
 * @param TYPE the type of exception to throw
 */
#define GEOS_THROW_IF_LE_MSG( lhs, rhs, msg, TYPE ) LOG_GEOS_THROW_IF_LE_MSG( lhs, rhs, "***** Rank " << ::geos::logger::internal::rankString << ": " << msg, TYPE )

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
 * @brief Abort execution if @p EXP is false but only when
 *        NDEBUG is not defined..
 * @param EXP The expression to check.
 * @param MSG The message to associate with the error, can be anything streamable to a std::ostream.
 * @note This macro can be used in both host and device code.
 * @note Tries to provide as much information about the location of the error
 *       as possible. On host this should result in the file and line of the error
 *       and a stack trace along with the provided message. On device none of this is
 *       guaranteed. In fact it is only guaranteed to abort the current kernel.
 */
#if !defined(NDEBUG)
#define GEOS_ASSERT_MSG_IF( EXP, MSG ) GEOS_ERROR_IF_IMPL( !(EXP), MSG )
#else
#define GEOS_ASSERT_MSG_IF( EXP, MSG ) ((void) 0)
#endif

/**
 * @brief Assert a condition in debug builds.
 * @param EXP an expression that will be evaluated as a predicate
 * @param msg a message to log (any expression that can be stream inserted)
 */
#define GEOS_ASSERT_MSG( EXP, msg ) GEOS_ASSERT_MSG_IF( EXP, "***** Rank " << ::geos::logger::internal::rankString << ": " << msg )

/**
 * @brief Assert a condition in debug builds.
 * @param EXP an expression that will be evaluated as a predicate
 */
#define GEOS_ASSERT( EXP ) GEOS_ASSERT_MSG( EXP, "" )

/**
 * @brief Abort execution if @p lhs @p OP @p rhs is false.
 * @param lhs The left side of the operation.
 * @param OP The operation to apply.
 * @param rhs The right side of the operation.
 * @param msg The message to diplay.
 */
#define GEOS_ASSERT_OP_MSG( lhs, OP, rhs, msg ) \
  GEOS_ASSERT_MSG_IF( lhs OP rhs, \
                      msg << "\n" << \
                      "  " << #lhs << " = " << lhs << "\n" << \
                      "  " << #rhs << " = " << rhs << "\n" )

/**
 * @brief Assert that two values compare equal in debug builds.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param msg a message to log (any expression that can be stream inserted)
 */
#define GEOS_ASSERT_EQ_MSG_IF( lhs, rhs, msg ) GEOS_ASSERT_OP_MSG( lhs, ==, rhs, msg )

/**
 * @brief Assert that two values compare equal in debug builds.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param msg a message to log (any expression that can be stream inserted)
 */
#define GEOS_ASSERT_EQ_MSG( lhs, rhs, msg ) GEOS_ASSERT_EQ_MSG_IF( lhs, rhs, "***** Rank " << ::geos::logger::internal::rankString << ": " << msg )

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
#define GEOS_ASSERT_NE_MSG_IF( lhs, rhs, msg ) GEOS_ASSERT_OP_MSG( lhs, !=, rhs, msg )

/**
 * @brief Assert that two values compare not equal in debug builds.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param msg a message to log (any expression that can be stream inserted)
 */
#define GEOS_ASSERT_NE_MSG( lhs, rhs, msg ) GEOS_ASSERT_NE_MSG_IF( lhs, rhs, msg )

/**
 * @brief Assert that two values compare not equal in debug builds.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 */
#define GEOS_ASSERT_NE_IF( lhs, rhs ) GEOS_ASSERT_NE_MSG( lhs, rhs, "" )

/**
 * @brief Assert that two values compare not equal in debug builds.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 */
#define GEOS_ASSERT_NE( lhs, rhs ) GEOS_ASSERT_NE_IF( lhs, rhs )

/**
 * @brief Assert that one value compares greater than the other in debug builds.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param msg a message to log (any expression that can be stream inserted)
 */
#define GEOS_ASSERT_GT_MSG_IF( lhs, rhs, msg ) GEOS_ASSERT_OP_MSG( lhs, >, rhs, msg )

/**
 * @brief Assert that one value compares greater than the other in debug builds.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param msg a message to log (any expression that can be stream inserted)
 */
#define GEOS_ASSERT_GT_MSG( lhs, rhs, msg ) GEOS_ASSERT_GT_MSG_IF( lhs, rhs, "***** Rank " << ::geos::logger::internal::rankString << ": " << msg )

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
#define GEOS_ASSERT_GE_MSG_IF( lhs, rhs, msg ) GEOS_ASSERT_OP_MSG( lhs, >=, rhs, msg )

/**
 * @brief Assert that one value compares greater than or equal to the other in debug builds.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param msg a message to log (any expression that can be stream inserted)
 */
#define GEOS_ASSERT_GE_MSG( lhs, rhs, msg ) GEOS_ASSERT_GE_MSG_IF( lhs, rhs, "***** Rank " << ::geos::logger::internal::rankString << ": " << msg )

/**
 * @brief Assert that one value compares greater than or equal to the other in debug builds.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 */
#define GEOS_ASSERT_GE( lhs, rhs ) GEOS_ASSERT_GE_MSG( lhs, rhs, "" )

namespace geos
{

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
   * @brief Constructs an InputError from an underlying exception.
   * @param subException The exception on which the created one is based.
   * @param msgToInsert The error message that will be inserted in the subException error message.
   */
  InputError( std::exception const & subException, std::string const & msgToInsert );
};

/**
 * @brief Exception class used to report errors in user input.
 */
struct SimulationError : public std::runtime_error
{
  /**
   * @brief Constructor
   * @param what the error message
   */
  SimulationError( std::string const & what ):
    std::runtime_error( what )
  {}

  /**
   * @brief Constructor
   * @param what the error message
   */
  SimulationError( char const * const what ):
    std::runtime_error( what )
  {}

  /**
   * @brief Construct a SimulationError from an underlying exception.
   * @param subException An exception to base this new one on.
   * @param msgToInsert The error message.
   * It will be inserted before the error message inside of subException.
   */
  SimulationError( std::exception const & subException, std::string const & msgToInsert );
};

/**
 * @brief Exception class used to report errors from type conversion
 * @todo (ErrorManager EPIC #2940) Consider adding a way to precise custom exception parameters, to add
 * expected & encountered typeid for this one (in order to manage the exception output more precisely).
 * We could also manage this by having:    BadTypeErrorABC <|--- BadTypeError< T >      /!\ compilation time
 */
struct BadTypeError : public std::runtime_error
{
  /**
   * @brief Constructor
   * @param what the error message
   */
  BadTypeError( std::string const & what ):
    std::runtime_error( what )
  {}
};

/**
 * @brief Exception class used for special control flow.
 */
class NotAnError : public std::exception
{};

namespace logger
{

namespace internal
{

extern int rank;

extern std::string rankString;

extern int n_ranks;

extern std::ostream * rankStream;

#if defined(GEOS_USE_MPI)
extern MPI_Comm comm;
#endif
} // namespace internal

#if defined(GEOS_USE_MPI)
/**
 * @brief Initialize the logger in a parallel build.
 * @param comm global MPI communicator
 * @param rank_output_dir output directory for rank log files
 */
void InitializeLogger( MPI_Comm comm, const std::string & rank_output_dir="" );
#endif

/**
 * @brief Initialize the logger in a serial build.
 * @param rank_output_dir output directory for rank log files
 */
void InitializeLogger( const std::string & rank_output_dir="" );

/**
 * @brief Finalize the logger and close the rank streams.
 */
void FinalizeLogger();

} // namespace logger

} // namespace geos

#endif /* GEOS_COMMON_LOGGER_HPP */
