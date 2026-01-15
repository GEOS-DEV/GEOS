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

// Source includes
#include "LvArray/src/Macros.hpp"
#include "common/logger/GeosExceptions.hpp"

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
    if( ::geos::logger::internal::g_rank == 0 && EXP ) \
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
    if( ::geos::logger::internal::g_rank == 0 && EXP ) \
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
      oss << "Rank " << ::geos::logger::internal::g_rankString << ": " << msg; \
      *logger::internal::g_rankStream << oss.str() << std::endl; \
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
 * @brief Error logger instance to use in GEOS_ERROR*, GEOS_ASSERT*, GEOS_THROW*, GEOS_WARNING* macros.
 * @note - Currently not available on GPU.
 *       - Possible to pre-define it in any source file (e.g. for unit tests)
 */
#if !defined(GEOS_DEVICE_COMPILE) && !defined(GEOS_GLOBAL_LOGGER)
#define GEOS_GLOBAL_LOGGER ErrorLogger::global()
#endif

/**
 * @brief Conditionally raise a hard error and terminate the program.
 *        Implementation of GEOS_ERROR_* and GEOS_ASSERT_* macros.
 * @param COND A condition that causes the error if true.
 * @param CAUSE_MESSAGE The condition that caused the error, in a readable text format for the user.
 * @param ... Variable arguments with the following structure:
 *            - Mandatory first parameter, the message to log (must be streamable)
 *            - Optional following parameters, context information on the current error (DataContext)
 */
#if !defined(GEOS_DEVICE_COMPILE)
#define GEOS_ERROR_IF_CAUSE( COND, CAUSE_MESSAGE, ... ) \
  do \
  { \
    if( COND ) \
    { \
      std::ostringstream __msgoss; \
      __msgoss << GEOS_DETAIL_FIRST_ARG( __VA_ARGS__ ); \
      std::ostringstream __causemsgsoss; \
      __causemsgsoss << CAUSE_MESSAGE; \
      GEOS_GLOBAL_LOGGER.initCurrentExceptionMessage( MsgType::Error, __msgoss.str(), \
                                                      ::geos::logger::internal::g_rank ) \
        .setCodeLocation( __FILE__, __LINE__ ) \
        .setCause( __causemsgsoss.str() ) \
        .addCallStackInfo( LvArray::system::stackTrace( true ) ) \
        .addContextInfo( GEOS_DETAIL_REST_ARGS( __VA_ARGS__ )) \
      GEOS_GLOBAL_LOGGER.flushCurrentExceptionMessage(); \
      LvArray::system::callErrorHandler(); \
    } \
  }while( false )
  #elif __CUDA_ARCH__
#define GEOS_ERROR_IF_CAUSE( COND, CAUSE_MESSAGE, ... ) \
  do \
  { \
    if( COND ) \
    { \
      constexpr char const * formatString = "***** ERROR\n" \
                                            "***** LOCATION" LOCATION "\n" \
                                                                      "***** BLOCK:  [%u, %u, %u]\n" \
                                                                      "***** THREAD: [%u, %u, %u]\n" \
                                                                      "***** " STRINGIZE( CAUSE_MESSAGE ) "\n" \
                                                                                                          "***** " STRINGIZE( GEOS_DETAIL_FIRST_ARG( __VA_ARGS__ ) ) "\n\n"; \
      printf( formatString, blockIdx.x, blockIdx.y, blockIdx.z, threadIdx.x, threadIdx.y, threadIdx.z ); \
      asm ( "trap;" ); \
    } \
  } while( false )
#endif

/**
 * @brief Conditionally raise a hard error and terminate the program.
 * @param COND A condition that causes the error if true.
 * @param ... Variable arguments with the following structure:
 *            - Mandatory first parameter, the message to log (must be streamable)
 *            - Optional following parameters, context information on the current error (DataContext)
 */
#define GEOS_ERROR_IF( COND, ... ) \
  GEOS_ERROR_IF_CAUSE( COND, "Error cause: " STRINGIZE( COND ), __VA_ARGS__ )

/**
 * @brief Raise a hard error and terminate the program.
 * @param ... Variable arguments with the following structure:
 *            - Mandatory first parameter, the message to log (must be streamable)
 *            - Optional following parameters, context information on the current error (DataContext)
 */
#define GEOS_ERROR( ... ) GEOS_ERROR_IF_CAUSE( true, "", __VA_ARGS__ )

/**
 * @brief Conditionally throw an exception.
 * @param COND an expression that will be evaluated as a predicate
 * @param CAUSE_MESSAGE The condition that caused the error, in a readable text format for the user.
 * @param MSG a message to log (any expression that can be stream inserted)
 * @param ... Variable arguments with the following structure:
 *            - Mandatory first parameter, the type of the exception to throw
 *            - Optional following parameters, context information on the current error (DataContext)
 */
 #if !defined(GEOS_DEVICE_COMPILE)
#define GEOS_THROW_IF_CAUSE( COND, CAUSE_MESSAGE, MSG, ... ) \
  do \
  { \
    if( COND ) \
    { \
      std::ostringstream __msgoss; \
      __msgoss << MSG; \
      std::ostringstream __causemsgsoss; \
      __causemsgsoss << CAUSE_MESSAGE; \
      DiagnosticMsg exceptionMsg =  GEOS_GLOBAL_LOGGER.initCurrentExceptionMessage( MsgType::Exception, __msgoss.str(), \
                                                                                    ::geos::logger::internal::g_rank ) \
                                     .setCodeLocation( __FILE__, __LINE__ ) \
                                     .setCause( __causemsgsoss.str() ) \
                                     .addCallStackInfo( LvArray::system::stackTrace( true ) ) \
                                     .addContextInfo( GEOS_DETAIL_REST_ARGS( __VA_ARGS__ )) \
                                     .getDiagnosticMsg(); \
      auto ex = GEOS_DETAIL_FIRST_ARG( __VA_ARGS__ )(); \
      ex.prepareWhat( exceptionMsg ); \
      throw ex; \
    } \
  }while( false )
  #elif __CUDA_ARCH__
#define GEOS_THROW_IF_CAUSE( COND, CAUSE_MESSAGE, MSG, ... ) \
  do \
  { \
    if( COND ) \
    { \
      static char const formatString[] = "***** ERROR\n" \
                                         "***** LOCATION" LOCATION "\n" \
                                                                   "***** BLOCK:  [%u, %u, %u]\n" \
                                                                   "***** THREAD: [%u, %u, %u]\n" \
                                                                   "***** " STRINGIZE( CAUSE_MESSAGE ) "\n" \
                                                                                                       "***** " STRINGIZE( GEOS_DETAIL_FIRST_ARG( __VA_ARGS__ ) ) "\n\n"; \
      printf( formatString, blockIdx.x, blockIdx.y, blockIdx.z, threadIdx.x, threadIdx.y, threadIdx.z ); \
      asm ( "trap;" ); \
    } \
  } while( false )
#endif

/**
 * @brief Conditionally raise a hard error and terminate the program.
 * @param COND A condition that causes the error if true.
 * @param MSG a message to log (any expression that can be stream inserted)
 * @param ... Variable arguments with the following structure:
 *            - Mandatory first parameter, the type of the exception to throw
 *            - Optional following parameters, context information on the current error (DataContext)
 */
#define GEOS_THROW_IF( COND, MSG, ... ) \
  GEOS_THROW_IF_CAUSE( COND, "Error cause: " STRINGIZE( COND ), MSG, __VA_ARGS__ )

/**
 * @brief Conditionally raise a hard error and terminate the program.
 * @param MSG a message to log (any expression that can be stream inserted)
 * @param ... Variable arguments with the following structure:
 *            - Mandatory first parameter, the type of the exception to throw
 *            - Optional following parameters, context information on the current error (DataContext)
 */
#define GEOS_THROW( MSG, ... ) GEOS_THROW_IF_CAUSE( true, "", MSG, __VA_ARGS__ )

/**
 * @brief Conditionally report a warning
 * @param COND A condition that causes the error if true.
 * @param CAUSE_MESSAGE The condition that caused the error, in a readable text format for the user.
 * @param ... Variable arguments with the following structure:
 *            - Mandatory first parameter, the message to log (must be streamable)
 *            - Optional following parameters, context information on the current error (DataContext)
 */
#if !defined(GEOS_DEVICE_COMPILE)
#define GEOS_WARNING_IF_CAUSE( COND, CAUSE_MESSAGE, ... ) \
  do \
  { \
    if( COND ) \
    { \
      std::ostringstream __msgoss; \
      __msgoss << GEOS_DETAIL_FIRST_ARG( __VA_ARGS__ ); \
      std::ostringstream __causemsgsoss; \
      __causemsgsoss << CAUSE_MESSAGE; \
      DiagnosticMsg __warningMsg; \
      GEOS_GLOBAL_LOGGER.flushErrorMsg( DiagnosticMsgBuilder::init( __warningMsg, \
                                                                    MsgType::Warning, __msgoss.str(), \
                                                                    ::geos::logger::internal::g_rank ) \
                                          .setCodeLocation( __FILE__, __LINE__ ) \
                                          .setCause( __causemsgsoss.str() ) \
                                          .addCallStackInfo( LvArray::system::stackTrace( true ) ) \
                                          .addContextInfo( GEOS_DETAIL_REST_ARGS( __VA_ARGS__ )) \
                                          .getDiagnosticMsg() ); \
    } \
  }while( false )
#elif __CUDA_ARCH__
#define GEOS_WARNING_IF_CAUSE( COND, CAUSE_MESSAGE, ... ) \
  do \
  { \
    if( COND ) \
    { \
      static char const formatString[] = "***** WARNING\n" \
                                         "***** LOCATION" LOCATION "\n" \
                                                                   "***** BLOCK:  [%u, %u, %u]\n" \
                                                                   "***** THREAD: [%u, %u, %u]\n" \
                                                                   "***** " STRINGIZE( CAUSE_MESSAGE ) "\n" \
                                                                                                       "***** " STRINGIZE( GEOS_DETAIL_FIRST_ARG( __VA_ARGS__ ) ) "\n\n"; \
      printf( formatString, blockIdx.x, blockIdx.y, blockIdx.z, threadIdx.x, threadIdx.y, threadIdx.z ); \
      asm ( "trap;" ); \
    } \
  } while( false )
#endif

/**
 * @brief Conditionally report a warning.
 * @param COND an expression that will be evaluated as a predicate
 * @param ... Variable arguments with the following structure:
 *            - Mandatory first parameter, the message to log (must be streamable)
 *            - Optional following parameters, context information on the current error (DataContext)
 */
#define GEOS_WARNING_IF( COND, ... ) \
  GEOS_WARNING_IF_CAUSE( COND, "Warning cause: " STRINGIZE( COND ), __VA_ARGS__ )

/**
 * @brief Report a warning.
 * @param ... Variable arguments with the following structure:
 *            - Mandatory first parameter, the message to log (must be streamable)
 *            - Optional following parameters, context information on the current error (DataContext)
 */
#define GEOS_WARNING( ... ) GEOS_WARNING_IF_CAUSE( true, "", __VA_ARGS__ )

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
 * @brief Declares variables for "assertion" evaluation only on CPU; no-op on GPU to avoid device compilation errors.
 * @param lhs The left side of the operation.
 * @param rhs The right side of the operation.
 */
#define GEOS_ERROR_LHS_RHS_DECLS( lhs, rhs ) \
  GEOS_MAYBE_UNUSED auto const lhsResult = (lhs); \
  GEOS_MAYBE_UNUSED auto const rhsResult = (rhs)

/**
 * @brief Abort execution if @p lhs @p OP @p rhs.
 * @param lhs The left side of the operation.
 * @param OP The operation to apply.
 * @param NOP The operation that caused the error, used in the message (typically opposite of @p OP).
 * @param rhs The right side of the operation.
 * @param ... Variable arguments with the following structure:
 *            - Mandatory first parameter, the message to log (must be streamable)
 *            - Optional following parameters, context information on the current error (DataContext)
 */
#define GEOS_ERROR_IF_OP_MSG( lhs, OP, NOP, rhs, ... ) \
  do { \
    GEOS_ERROR_LHS_RHS_DECLS( lhs, rhs ); \
    GEOS_ERROR_IF_CAUSE( lhsResult OP rhsResult, \
                         "Expected: " #lhs " " #NOP " " #rhs "\n* " #lhs " = " << lhsResult << "\n* " #rhs " = " << rhsResult << "\n", \
                         __VA_ARGS__ ); \
  } while(false)


/**
 * @brief Raise a hard error if two values are equal.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param ... Variable arguments with the following structure:
 *            - Mandatory first parameter, the message to log (must be streamable)
 *            - Optional following parameters, context information on the current error (DataContext)
 */
#define GEOS_ERROR_IF_EQ_MSG( lhs, rhs, ... ) GEOS_ERROR_IF_OP_MSG( lhs, ==, !=, rhs, __VA_ARGS__ )

/**
 * @brief Raise a hard error if two values are equal.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 */
#define GEOS_ERROR_IF_EQ( lhs, rhs ) GEOS_ERROR_IF_EQ_MSG( lhs, rhs, "" )

/**
 * @brief Raise a hard error if two values are not equal.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param ... Variable arguments with the following structure:
 *            - Mandatory first parameter, the message to log (must be streamable)
 *            - Optional following parameters, context information on the current error (DataContext)
 */
#define GEOS_ERROR_IF_NE_MSG( lhs, rhs, ... ) GEOS_ERROR_IF_OP_MSG( lhs, !=, ==, rhs, __VA_ARGS__ )

/**
 * @brief Raise a hard error if two values are not equal.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 */
#define GEOS_ERROR_IF_NE( lhs, rhs ) GEOS_ERROR_IF_NE_MSG( lhs, rhs, "" )

/**
 * @brief Raise a hard error if one value compares greater than the other.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param ... Variable arguments with the following structure:
 *            - Mandatory first parameter, the message to log (must be streamable)
 *            - Optional following parameters, context information on the current error (DataContext)
 */
#define GEOS_ERROR_IF_GT_MSG( lhs, rhs, ... ) GEOS_ERROR_IF_OP_MSG( lhs, >, <=, rhs, __VA_ARGS__ )

/**
 * @brief Raise a hard error if one value compares greater than the other.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 */
#define GEOS_ERROR_IF_GT( lhs, rhs ) GEOS_ERROR_IF_GT_MSG( lhs, rhs, "" )

/**
 * @brief Raise a hard error if one value compares greater than or equal to the other.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param ... Variable arguments with the following structure:
 *            - Mandatory first parameter, the message to log (must be streamable)
 *            - Optional following parameters, context information on the current error (DataContext)
 */
#define GEOS_ERROR_IF_GE_MSG( lhs, rhs, ... ) GEOS_ERROR_IF_OP_MSG( lhs, >=, <, rhs, __VA_ARGS__ )

/**
 * @brief Raise a hard error if one value compares greater than or equal to the other.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 */
#define GEOS_ERROR_IF_GE( lhs, rhs ) GEOS_ERROR_IF_GE_MSG( lhs, rhs, "" )

/**
 * @brief Raise a hard error if one value compares less than the other.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param ... Variable arguments with the following structure:
 *            - Mandatory first parameter, the message to log (must be streamable)
 *            - Optional following parameters, context information on the current error (DataContext)
 */
#define GEOS_ERROR_IF_LT_MSG( lhs, rhs, ... ) GEOS_ERROR_IF_OP_MSG( lhs, <, >=, rhs, __VA_ARGS__ )

/**
 * @brief Raise a hard error if one value compares less than the other.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 */
#define GEOS_ERROR_IF_LT( lhs, rhs ) GEOS_ERROR_IF_LT_MSG( lhs, rhs, "" )

/**
 * @brief Raise a hard error if one value compares less than or equal to the other.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param ... Variable arguments with the following structure:
 *            - Mandatory first parameter, the message to log (must be streamable)
 *            - Optional following parameters, context information on the current error (DataContext)
 */
#define GEOS_ERROR_IF_LE_MSG( lhs, rhs, ... ) GEOS_ERROR_IF_OP_MSG( lhs, <=, >, rhs, __VA_ARGS__ )


/**
 * @brief Raise a hard error if one value compares less than or equal to the other.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 */
#define GEOS_ERROR_IF_LE( lhs, rhs ) GEOS_ERROR_IF_LE_MSG( lhs, rhs, "" )

/**
 * @brief Log a warning if @p lhs @p OP @p rhs.
 * @param lhs The left side of the operation.
 * @param OP The operation to apply.
 * @param NOP The operation that caused the error, used in the message (typically opposite of @p OP).
 * @param rhs The right side of the operation.
 * @param ... Variable arguments with the following structure:
 *            - Mandatory first parameter, the message to log (must be streamable)
 *            - Optional following parameters, context information on the current error (DataContext)
 */
#define GEOS_WARNING_IF_OP_MSG( lhs, OP, NOP, rhs, ... ) \
  do { \
    GEOS_ERROR_LHS_RHS_DECLS( lhs, rhs ); \
    GEOS_WARNING_IF_CAUSE( lhsResult OP rhsResult, \
                           "Expected: " #lhs " " #NOP " " #rhs "\n* " #lhs " = " << lhsResult << "\n* " #rhs " = " << rhsResult << "\n", \
                           __VA_ARGS__ ); \
  } while(false)


/**
 * @brief Log a warning if two values are equal.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param ... Variable arguments with the following structure:
 *            - Mandatory first parameter, the message to log (must be streamable)
 *            - Optional following parameters, context information on the current error (DataContext)
 */
#define GEOS_WARNING_IF_EQ_MSG( lhs, rhs, ... ) GEOS_WARNING_IF_OP_MSG( lhs, ==, !=, rhs, __VA_ARGS__ )

/**
 * @brief Log a warning if two values are equal.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 */
#define GEOS_WARNING_IF_EQ( lhs, rhs ) GEOS_WARNING_IF_EQ_MSG( lhs, rhs, "" )

/**
 * @brief Log a warning if two values are not equal.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param ... Variable arguments with the following structure:
 *            - Mandatory first parameter, the message to log (must be streamable)
 *            - Optional following parameters, context information on the current error (DataContext)
 */
#define GEOS_WARNING_IF_NE_MSG( lhs, rhs, ... ) GEOS_WARNING_IF_OP_MSG( lhs, !=, ==, rhs, __VA_ARGS__ )

/**
 * @brief Log a warning if two values are not equal.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 */
#define GEOS_WARNING_IF_NE( lhs, rhs ) GEOS_WARNING_IF_NE_MSG( lhs, rhs, "" )

/**
 * @brief Log a warning if one value compares greater than the other.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param ... Variable arguments with the following structure:
 *            - Mandatory first parameter, the message to log (must be streamable)
 *            - Optional following parameters, context information on the current error (DataContext)
 */
#define GEOS_WARNING_IF_GT_MSG( lhs, rhs, ... ) GEOS_WARNING_IF_OP_MSG( lhs, >, <=, rhs, __VA_ARGS__ )

/**
 * @brief Log a warning if one value compares greater than the other.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 */
#define GEOS_WARNING_IF_GT( lhs, rhs ) GEOS_WARNING_IF_GT_MSG( lhs, rhs, "" )

/**
 * @brief Log a warning if one value compares greater than or equal to the other.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param ... Variable arguments with the following structure:
 *            - Mandatory first parameter, the message to log (must be streamable)
 *            - Optional following parameters, context information on the current error (DataContext)
 */
#define GEOS_WARNING_IF_GE_MSG( lhs, rhs, ... ) GEOS_WARNING_IF_OP_MSG( lhs, >=, <, rhs, __VA_ARGS__ )

/**
 * @brief Log a warning if one value compares greater than or equal to the other.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 */
#define GEOS_WARNING_IF_GE( lhs, rhs ) GEOS_WARNING_IF_GE_MSG( lhs, rhs, "" )

/**
 * @brief Log a warning if one value compares less than the other.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param ... Variable arguments with the following structure:
 *            - Mandatory first parameter, the message to log (must be streamable)
 *            - Optional following parameters, context information on the current error (DataContext)
 */
#define GEOS_WARNING_IF_LT_MSG( lhs, rhs, ... ) GEOS_WARNING_IF_OP_MSG( lhs, <, >=, rhs, __VA_ARGS__ )

/**
 * @brief Log a warning if one value compares less than the other.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 */
#define GEOS_WARNING_IF_LT( lhs, rhs ) GEOS_WARNING_IF_LT_MSG( lhs, rhs, "" )

/**
 * @brief Log a warning if one value compares less than or equal to the other.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param ... Variable arguments with the following structure:
 *            - Mandatory first parameter, the message to log (must be streamable)
 *            - Optional following parameters, context information on the current error (DataContext)
 */
#define GEOS_WARNING_IF_LE_MSG( lhs, rhs, ... ) GEOS_WARNING_IF_OP_MSG( lhs, <=, >, rhs, __VA_ARGS__ )

/**
 * @brief Log a warning if one value compares less than or equal to the other.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 */
#define GEOS_WARNING_IF_LE( lhs, rhs ) GEOS_WARNING_IF_LE_MSG( lhs, rhs, "" )

/**
 * @brief Throw an exception if @p lhs @p OP @p rhs.
 * @param lhs The left side of the operation.
 * @param OP The operation to apply.
 * @param NOP The operation that caused the error, used in the message (typically opposite of @p OP).
 * @param rhs The right side of the operation.
 * @param MSG a message to log (any expression that can be stream inserted)
 * @param ... Variable arguments with the following structure:
 *            - Mandatory first parameter, the type of the exception to throw
 *            - Optional following parameters, context information on the current error (DataContext)
 */
#define GEOS_THROW_IF_OP_MSG( lhs, OP, NOP, rhs, MSG, ... ) \
  do { \
    GEOS_ERROR_LHS_RHS_DECLS( lhs, rhs ); \
    GEOS_THROW_IF_CAUSE( lhsResult OP rhsResult, \
                         "Expected: " #lhs " " #NOP " " #rhs "\n* " #lhs " = " << lhsResult << "\n* " #rhs " = " << rhsResult << "\n", \
                         MSG, __VA_ARGS__ ); \
  } while(false)



/**
 * @brief Raise a hard error if two values are equal.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param MSG a message to log (any expression that can be stream inserted)
 * @param ... Variable arguments with the following structure:
 *            - Mandatory first parameter, the type of the exception to throw
 *            - Optional following parameters, context information on the current error (DataContext)
 */
#define GEOS_THROW_IF_EQ_MSG( lhs, rhs, MSG, ... ) GEOS_THROW_IF_OP_MSG( lhs, ==, !=, rhs, MSG, __VA_ARGS__ )

/**
 * @brief Raise a hard error if two values are equal.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param ... Variable arguments with the following structure:
 *            - Mandatory first parameter, the type of the exception to throw
 *            - Optional following parameters, context information on the current error (DataContext)
 */
#define GEOS_THROW_IF_EQ( lhs, rhs, ... ) GEOS_THROW_IF_EQ_MSG( lhs, rhs, "", __VA_ARGS__ )

/**
 * @brief Throw an exception if two values are not equal.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param MSG a message to log (any expression that can be stream inserted)
 * @param ... Variable arguments with the following structure:
 *            - Mandatory first parameter, the type of the exception to throw
 *            - Optional following parameters, context information on the current error (DataContext)
 */
#define GEOS_THROW_IF_NE_MSG( lhs, rhs, MSG, ... ) GEOS_THROW_IF_OP_MSG( lhs, !=, ==, rhs, MSG, __VA_ARGS__ )

/**
 * @brief Throw an exception if two values are not equal.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param ... Variable arguments with the following structure:
 *            - Mandatory first parameter, the type of the exception to throw
 *            - Optional following parameters, context information on the current error (DataContext)
 */
#define GEOS_THROW_IF_NE( lhs, rhs, ... ) GEOS_THROW_IF_NE_MSG( lhs, rhs, "", __VA_ARGS__ )

/**
 * @brief Throw an exception if one value compares greater than the other.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param MSG a message to log (any expression that can be stream inserted)
 * @param ... Variable arguments with the following structure:
 *            - Mandatory first parameter, the type of the exception to throw
 *            - Optional following parameters, context information on the current error (DataContext)
 */
#define GEOS_THROW_IF_GT_MSG( lhs, rhs, MSG, ... ) GEOS_THROW_IF_OP_MSG( lhs, >, <=, rhs, MSG, __VA_ARGS__ )

/**
 * @brief Throw an exception if one value compares greater than the other.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param ... Variable arguments with the following structure:
 *            - Mandatory first parameter, the type of the exception to throw
 *            - Optional following parameters, context information on the current error (DataContext)
 */
#define GEOS_THROW_IF_GT( lhs, rhs, ... ) GEOS_THROW_IF_GT_MSG( lhs, rhs, "", __VA_ARGS__ )

/**
 * @brief Throw an exception if one value compares greater than or equal to the other.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param MSG a message to log (any expression that can be stream inserted)
 * @param ... Variable arguments with the following structure:
 *            - Mandatory first parameter, the type of the exception to throw
 *            - Optional following parameters, context information on the current error (DataContext)
 */
#define GEOS_THROW_IF_GE_MSG( lhs, rhs, MSG, ... ) GEOS_THROW_IF_OP_MSG( lhs, >=, <, rhs, MSG, __VA_ARGS__ )

/**
 * @brief Throw an exception if one value compares greater than or equal to the other.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param ... Variable arguments with the following structure:
 *            - Mandatory first parameter, the type of the exception to throw
 *            - Optional following parameters, context information on the current error (DataContext)
 */
#define GEOS_THROW_IF_GE( lhs, rhs, ... ) GEOS_THROW_IF_GE_MSG( lhs, rhs, "", __VA_ARGS__ )

/**
 * @brief Throw an exception if one value compares less than the other.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param MSG a message to log (any expression that can be stream inserted)
 * @param ... Variable arguments with the following structure:
 *            - Mandatory first parameter, the type of the exception to throw
 *            - Optional following parameters, context information on the current error (DataContext)
 */
#define GEOS_THROW_IF_LT_MSG( lhs, rhs, MSG, ... ) GEOS_THROW_IF_OP_MSG( lhs, <, >=, rhs, MSG, __VA_ARGS__ )

/**
 * @brief Throw an exception if one value compares less than the other.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param ... Variable arguments with the following structure:
 *            - Mandatory first parameter, the type of the exception to throw
 *            - Optional following parameters, context information on the current error (DataContext)
 */
#define GEOS_THROW_IF_LT( lhs, rhs, ... ) GEOS_THROW_IF_LT_MSG( lhs, rhs, "", __VA_ARGS__ )

/**
 * @brief Throw an exception if one value compares less than or equal to the other.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param MSG a message to log (any expression that can be stream inserted)
 * @param ... Variable arguments with the following structure:
 *            - Mandatory first parameter, the type of the exception to throw
 *            - Optional following parameters, context information on the current error (DataContext)
 */
#define GEOS_THROW_IF_LE_MSG( lhs, rhs, MSG, ... ) GEOS_THROW_IF_OP_MSG( lhs, <=, >, rhs, MSG, __VA_ARGS__ )

/**
 * @brief Throw an exception if one value compares less than or equal to the other.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param ... Variable arguments with the following structure:
 *            - Mandatory first parameter, the type of the exception to throw
 *            - Optional following parameters, context information on the current error (DataContext)
 */
#define GEOS_THROW_IF_LE( lhs, rhs, ... ) GEOS_THROW_IF_LE_MSG( lhs, rhs, "", __VA_ARGS__ )

#if !defined(NDEBUG) || defined(GEOS_ASSERT_ENABLED)

/**
 * @brief Enables assertion macros (GEOS_ASSERT*) when NDEBUG is not defined or previously explicitly enabled.
 */
#define GEOS_ASSERT_ENABLED

/**
 * @brief Abort execution if @p COND is false but only when NDEBUG is not defined..
 * @param COND The condition to check, causes an error if false.
 * @param ... Variable arguments with the following structure:
 *            - Mandatory first parameter, the message to log (must be streamable)
 *            - Optional following parameters, context information on the current error (DataContext)
 * @note This macro can be used in both host and device code.
 * @note Tries to provide as much information about the location of the error
 *       as possible. On host this should result in the file and line of the error
 *       and a stack trace along with the provided message. On device none of this is
 *       guaranteed. In fact it is only guaranteed to abort the current kernel.
 */
#define  GEOS_ASSERT_MSG( COND, ... ) \
  GEOS_ERROR_IF_CAUSE( !( COND ), "Expected: " STRINGIZE( COND ), __VA_ARGS__ )

/**
 * @brief Abort execution if @p lhs @p OP @p rhs is false.
 * @param lhs The left side of the operation.
 * @param OP The operation to apply.
 * @param rhs The right side of the operation.
 * @param ... Variable arguments with the following structure:
 *            - Mandatory first parameter, the message to log (must be streamable)
 *            - Optional following parameters, context information on the current error (DataContext)
 */
#define GEOS_ASSERT_OP_MSG( lhs, OP, rhs, ... ) \
  { \
    GEOS_ERROR_LHS_RHS_DECLS( lhs, rhs ); \
    GEOS_ERROR_IF_CAUSE( !( lhsResult OP rhsResult ), \
                         "Expected: " #lhs " " #OP " " #rhs "\n* " #lhs " = " << lhsResult << "\n* " #rhs " = " << rhsResult << "\n", \
                         __VA_ARGS__ ); \
  }

#else

#define GEOS_ASSERT_MSG( ... ) ((void) 0)
#define GEOS_ASSERT_OP_MSG( ... ) ((void) 0)

#endif

/**
 * @brief Assert a condition in debug builds.
 * @param COND The condition to check, causes an error if false.
 */
#define GEOS_ASSERT( COND ) GEOS_ASSERT_MSG( COND, "" )

/**
 * @brief Assert that two values compare equal in debug builds.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param ... Variable arguments with the following structure:
 *            - Mandatory first parameter, the message to log (must be streamable)
 *            - Optional following parameters, context information on the current error (DataContext)
 */
#define GEOS_ASSERT_EQ_MSG( lhs, rhs, ... ) GEOS_ASSERT_OP_MSG( lhs, ==, rhs, __VA_ARGS__ )

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
 * @param ... Variable arguments with the following structure:
 *            - Mandatory first parameter, the message to log (must be streamable)
 *            - Optional following parameters, context information on the current error (DataContext)
 */
#define GEOS_ASSERT_NE_MSG( lhs, rhs, ... ) GEOS_ASSERT_OP_MSG( lhs, !=, rhs, __VA_ARGS__ )

/**
 * @brief Assert that two values compare not equal in debug builds.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 */
#define GEOS_ASSERT_NE( lhs, rhs ) GEOS_ASSERT_NE_MSG( lhs, rhs, "" )

/**
 * @brief Assert that one value compares greater than the other in debug builds.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param ... Variable arguments with the following structure:
 *            - Mandatory first parameter, the message to log (must be streamable)
 *            - Optional following parameters, context information on the current error (DataContext)
 */
#define GEOS_ASSERT_GT_MSG( lhs, rhs, ... ) GEOS_ASSERT_OP_MSG( lhs, >, rhs, __VA_ARGS__ )

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
 * @param ... Variable arguments with the following structure:
 *            - Mandatory first parameter, the message to log (must be streamable)
 *            - Optional following parameters, context information on the current error (DataContext)
 */
#define GEOS_ASSERT_GE_MSG( lhs, rhs, ... ) GEOS_ASSERT_OP_MSG( lhs, >=, rhs, __VA_ARGS__ )

/**
 * @brief Assert that one value compares greater than or equal to the other in debug builds.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 */
#define GEOS_ASSERT_GE( lhs, rhs ) GEOS_ASSERT_GE_MSG( lhs, rhs, "" )

/**
 * @brief Assert that one value compares greater than the other in debug builds.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param ... Variable arguments with the following structure:
 *            - Mandatory first parameter, the message to log (must be streamable)
 *            - Optional following parameters, context information on the current error (DataContext)
 */
#define GEOS_ASSERT_LT_MSG( lhs, rhs, ... ) GEOS_ASSERT_OP_MSG( lhs, <, rhs, __VA_ARGS__ )

/**
 * @brief Assert that one value compares greater than the other in debug builds.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 */
#define GEOS_ASSERT_LT( lhs, rhs ) GEOS_ASSERT_LT_MSG( lhs, rhs, "" )

/**
 * @brief Assert that one value compares greater than or equal to the other in debug builds.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param ... Variable arguments with the following structure:
 *            - Mandatory first parameter, the message to log (must be streamable)
 *            - Optional following parameters, context information on the current error (DataContext)
 */
#define GEOS_ASSERT_LE_MSG( lhs, rhs, ... ) GEOS_ASSERT_OP_MSG( lhs, <=, rhs, __VA_ARGS__ )

/**
 * @brief Assert that one value compares greater than or equal to the other in debug builds.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 */
#define GEOS_ASSERT_LE( lhs, rhs ) GEOS_ASSERT_LE_MSG( lhs, rhs, "" )

namespace geos
{

namespace logger
{

namespace internal
{

extern int g_rank;

extern std::string g_rankString;

extern std::ostream * g_rankStream;

}       // namespace internal

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

}     // namespace logger

}   // namespace geos

#endif /* GEOS_COMMON_LOGGER_HPP */
