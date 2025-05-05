/*
 * Copyright (c) 2021, Lawrence Livermore National Security, LLC and LvArray contributors.
 * All rights reserved.
 * See the LICENSE file for details.
 * SPDX-License-Identifier: (BSD-3-Clause)
 */

/**
 * @file Macros.hpp
 * @brief Contains a bunch of macro definitions.
 */

#pragma once

// Source includes
#include "LvArrayConfig.hpp"
#include "system.hpp"

// System includes
#include <fstream>
#include <sstream>
#include <iostream>
#include <type_traits>


#if defined(LVARRAY_USE_CUDA) || defined(LVARRAY_USE_HIP)
/// Macro defined when using a device.
#define LVARRAY_USE_DEVICE
#endif

#if defined(LVARRAY_USE_CUDA)
#define LVARRAY_DEFAULT_DEVICE_SPACE MemorySpace::cuda
#elif defined(LVARRAY_USE_HIP)
#define LVARRAY_DEFAULT_DEVICE_SPACE MemorySpace::hip
#endif

#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
/// Macro defined when currently compiling on device (only defined in the device context).
#define LVARRAY_DEVICE_COMPILE
/// Marks a function/lambda for inlining
#define LVARRAY_FORCE_INLINE __forceinline__
#else
/// Marks a function/lambda for inlining
#define LVARRAY_FORCE_INLINE inline
#endif

#if defined(__CUDACC__) || defined(__HIPCC__)
// Denotes whether to define decorator macros later in this file.
#define LVARRAY_DECORATE
#endif


//#if !defined(NDEBUG) && defined(LVARRAY_DEVICE_COMPILE)
  #include <cassert>
//#endif

/**
 * @brief Convert @p A into a string.
 * @param A the token to convert to a string.
 */
#define STRINGIZE_NX( A ) #A

/**
 * @brief Convert the macro expansion of @p A into a string.
 * @param A the token to convert to a string.
 */
#define STRINGIZE( A ) STRINGIZE_NX( A )

/**
 * @brief Mark @p X as an unused argument, used to silence compiler warnings.
 * @param X the unused argument.
 */
#define LVARRAY_UNUSED_ARG( X )

/**
 * @brief Mark @p X as an unused variable, used to silence compiler warnings.
 * @param X the unused variable.
 */
#define LVARRAY_UNUSED_VARIABLE( X ) ( ( void ) X )

/**
 * @brief Mark @p X as an debug variable, used to silence compiler warnings.
 * @param X the debug variable.
 */
#define LVARRAY_DEBUG_VAR( X ) LVARRAY_UNUSED_VARIABLE( X )

/// Expands to a string representing the current file and line.
#define LOCATION __FILE__ ":" STRINGIZE( __LINE__ )

/**
 * @brief Given an expression @p X that evaluates to a pointer, expands to the type pointed to.
 * @param X The expression to evaluate.
 */
#define TYPEOFPTR( X ) std::remove_pointer_t< decltype( X ) >

/**
 * @brief Given an expression @p X that evaluates to a reference, expands to the type referred to.
 * @param X The expression to evaluate.
 */
#define TYPEOFREF( X ) std::remove_reference_t< decltype( X ) >

/**
 * @brief Print the expression.
 */
#define LVARRAY_LOG( ... ) std::cout << __VA_ARGS__ << std::endl

/**
 * @brief Print the expression string along with its value.
 */
#define LVARRAY_LOG_VAR( ... ) LVARRAY_LOG( STRINGIZE( __VA_ARGS__ ) << " = " << __VA_ARGS__ )

/**
 * @brief Abort execution if @p EXP is true.
 * @param EXP The expression to check.
 * @param MSG The message to associate with the error, can be anything streamable to a std::ostream.
 * @note This macro can be used in both host and device code.
 * @note Tries to provide as much information about the location of the error
 *       as possible. On host this should result in the file and line of the error
 *       and a stack trace along with the provided message. On device none of this is
 *       guaranteed. In fact it is only guaranteed to abort the current kernel.
 */

#if defined(LVARRAY_DEVICE_COMPILE)
//   #if defined(__HIP_DEVICE_COMPILE__)
// // empty impl to avoid the possibility of printfs in device code
// //   on AMD, which can cause performance degradation just by being present
// #define LVARRAY_ERROR_IF( EXP, MSG )
  #if (!defined(NDEBUG)) || defined(__HIP_DEVICE_COMPILE__)
#define LVARRAY_ERROR_IF( EXP, MSG ) \
  do \
  { \
    if( EXP ) \
    { \
      assert( false && "EXP = " STRINGIZE( EXP ) "MSG = " STRINGIZE( MSG ) ); \
    } \
  } while( false )
  #else
#define LVARRAY_ERROR_IF( EXP, MSG ) \
  do \
  { \
    if( EXP ) \
    { \
      constexpr char const * formatString = "***** ERROR\n" \
                                            "***** LOCATION: " LOCATION "\n" \
                                                                        "***** Block: [%u, %u, %u]\n" \
                                                                        "***** Thread: [%u, %u, %u]\n" \
                                                                        "***** Controlling expression (should be false): " STRINGIZE( EXP ) "\n" \
                                                                                                                                            "***** MSG: " STRINGIZE( MSG ) "\n\n"; \
      printf( formatString, blockIdx.x, blockIdx.y, blockIdx.z, threadIdx.x, threadIdx.y, threadIdx.z ); \
      asm ( "trap;" ); \
    } \
  } while( false )
  #endif
#else
#define LVARRAY_ERROR_IF( EXP, MSG ) \
  do \
  { \
    if( EXP ) \
    { \
      std::ostringstream __oss; \
      __oss << "***** ERROR\n"; \
      __oss << "***** LOCATION: " LOCATION "\n"; \
      __oss << "***** Controlling expression (should be false): " STRINGIZE( EXP ) "\n"; \
      __oss << MSG << "\n"; \
      __oss << LvArray::system::stackTrace( true ); \
      std::cout << __oss.str() << std::endl; \
      LvArray::system::callErrorHandler(); \
    } \
  } while( false )
#endif

/**
 * @brief Abort execution.
 * @param MSG The message to associate with the error, can be anything streamable to a std::ostream.
 */
#define LVARRAY_ERROR( MSG ) LVARRAY_ERROR_IF( true, MSG )

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
#define LVARRAY_ASSERT_MSG( EXP, MSG ) LVARRAY_ERROR_IF( !(EXP), MSG )
#else
#define LVARRAY_ASSERT_MSG( EXP, MSG ) ((void) 0)
#endif

/**
 * @brief Conditionally throw an exception.
 * @param EXP an expression that will be evaluated as a predicate
 * @param MSG a message to log (any expression that can be stream inserted)
 * @param TYPE the type of exception to throw
 */
#define LVARRAY_THROW_IF( EXP, MSG, TYPE ) \
  do \
  { \
    if( EXP ) \
    { \
      std::ostringstream __oss; \
      __oss << "\n"; \
      __oss << "***** LOCATION: " LOCATION "\n"; \
      __oss << "***** Controlling expression (should be false): " STRINGIZE( EXP ) "\n"; \
      __oss << MSG << "\n"; \
      __oss << LvArray::system::stackTrace( true ); \
      throw TYPE( __oss.str() ); \
    } \
  } while( false )

/**
 * @brief Throw an exception.
 * @param MSG The message to associate with the error, can be anything streamable to a std::ostream.
 */
#define LVARRAY_THROW( MSG, TYPE ) LVARRAY_THROW_IF( true, MSG, TYPE )

/// Assert @p EXP is true with no message.
#define LVARRAY_ASSERT( EXP ) LVARRAY_ASSERT_MSG( EXP, "" )

/**
 * @brief Print a warning if @p EXP is true.
 * @param EXP The expression to check.
 * @param MSG The message to associate with the warning, can be anything streamable to a std::ostream.
 */
#define LVARRAY_WARNING_IF( EXP, MSG ) \
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
    } \
  } while( false )

/**
 * @brief Print a warning with a message.
 * @param MSG The message to print.
 */
#define LVARRAY_WARNING( MSG ) LVARRAY_WARNING_IF( true, MSG )

/**
 * @brief Print @p msg along with the location if @p EXP is true.
 * @param EXP The expression to test.
 * @param MSG The message to print.
 */
#define LVARRAY_INFO_IF( EXP, MSG ) \
  do \
  { \
    if( EXP ) \
    { \
      std::ostringstream __oss; \
      __oss << "***** INFO\n"; \
      __oss << "***** LOCATION: " LOCATION "\n"; \
      __oss << "***** Controlling expression: " STRINGIZE( EXP ) "\n"; \
      __oss << MSG; \
      std::cout << __oss.str() << std::endl; \
    } \
  } while( false )

/**
 * @brief Print @p msg along with the location.
 * @param msg The message to print.
 */
#define LVARRAY_INFO( msg ) LVARRAY_INFO_IF( true, msg )

/**
 * @brief Abort execution if @p lhs @p OP @p rhs.
 * @param lhs The left side of the operation.
 * @param OP The operation to apply.
 * @param NOP The opposite of @p OP, used in the message.
 * @param rhs The right side of the operation.
 * @param msg The message to diplay.
 */
#define LVARRAY_ERROR_IF_OP_MSG( lhs, OP, NOP, rhs, msg ) \
  LVARRAY_ERROR_IF( lhs OP rhs, \
                    msg << "\n" << \
                    "Expected " << #lhs << " " << #NOP << " " << #rhs << "\n" << \
                    "  " << #lhs << " = " << lhs << "\n" << \
                    "  " << #rhs << " = " << rhs << "\n" )

/**
 * @brief Throw an exception if @p lhs @p OP @p rhs.
 * @param lhs The left side of the operation.
 * @param OP The operation to apply.
 * @param NOP The opposite of @p OP, used in the message.
 * @param rhs The right side of the operation.
 * @param msg The message to diplay.
 * @param TYPE the type of exception to throw.
 */
#define LVARRAY_THROW_IF_OP_MSG( lhs, OP, NOP, rhs, msg, TYPE ) \
  LVARRAY_THROW_IF( lhs OP rhs, \
                    msg << "\n" << \
                    "Expected " << #lhs << " " << #NOP << " " << #rhs << "\n" << \
                    "  " << #lhs << " = " << lhs << "\n" << \
                    "  " << #rhs << " = " << rhs << "\n", TYPE )

/**
 * @brief Raise a hard error if two values are equal.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param msg a message to log (any expression that can be stream inserted)
 */
#define LVARRAY_ERROR_IF_EQ_MSG( lhs, rhs, msg ) LVARRAY_ERROR_IF_OP_MSG( lhs, ==, !=, rhs, msg )

/**
 * @brief Throw an exception if two values are equal.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param msg a message to log (any expression that can be stream inserted)
 * @param TYPE the type of exception to throw
 */
#define LVARRAY_THROW_IF_EQ_MSG( lhs, rhs, msg, TYPE ) LVARRAY_THROW_IF_OP_MSG( lhs, ==, !=, rhs, msg, TYPE )

/**
 * @brief Raise a hard error if two values are equal.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 */
#define LVARRAY_ERROR_IF_EQ( lhs, rhs ) LVARRAY_ERROR_IF_EQ_MSG( lhs, rhs, "" )

/**
 * @brief Throw an exception if two values are equal.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param TYPE the type of exception to throw
 */
#define LVARRAY_THROW_IF_EQ( lhs, rhs, TYPE ) LVARRAY_THROW_IF_EQ_MSG( lhs, rhs, "", TYPE )

/**
 * @brief Raise a hard error if two values are not equal.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param msg a message to log (any expression that can be stream inserted)
 */
#define LVARRAY_ERROR_IF_NE_MSG( lhs, rhs, msg ) LVARRAY_ERROR_IF_OP_MSG( lhs, !=, ==, rhs, msg )

/**
 * @brief Throw an exception if two values are not equal.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param msg a message to log (any expression that can be stream inserted)
 * @param TYPE the type of exception to throw
 */
#define LVARRAY_THROW_IF_NE_MSG( lhs, rhs, msg, TYPE ) LVARRAY_THROW_IF_OP_MSG( lhs, !=, ==, rhs, msg, TYPE )

/**
 * @brief Raise a hard error if two values are not equal.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 */
#define LVARRAY_ERROR_IF_NE( lhs, rhs ) LVARRAY_ERROR_IF_NE_MSG( lhs, rhs, "" )

/**
 * @brief Throw an exception if two values are not equal.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param TYPE the type of exception to throw
 */
#define LVARRAY_THROW_IF_NE( lhs, rhs, TYPE ) LVARRAY_THROW_IF_NE_MSG( lhs, rhs, "", TYPE )

/**
 * @brief Raise a hard error if one value compares greater than the other.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param msg a message to log (any expression that can be stream inserted)
 */
#define LVARRAY_ERROR_IF_GT_MSG( lhs, rhs, msg ) LVARRAY_ERROR_IF_OP_MSG( lhs, >, <=, rhs, msg )

/**
 * @brief Throw an exception if one value compares greater than the other.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param msg a message to log (any expression that can be stream inserted)
 * @param TYPE the type of exception to throw
 */
#define LVARRAY_THROW_IF_GT_MSG( lhs, rhs, msg, TYPE ) LVARRAY_THROW_IF_OP_MSG( lhs, >, <=, rhs, msg, TYPE )

/**
 * @brief Raise a hard error if one value compares greater than the other.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 */
#define LVARRAY_ERROR_IF_GT( lhs, rhs ) LVARRAY_ERROR_IF_GT_MSG( lhs, rhs, "" )

/**
 * @brief Throw an exception if one value compares greater than the other.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param TYPE the type of exception to throw
 */
#define LVARRAY_THROW_IF_GT( lhs, rhs, TYPE ) LVARRAY_THROW_IF_GT_MSG( lhs, rhs, "", TYPE )

/**
 * @brief Raise a hard error if one value compares greater than or equal to the other.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param msg a message to log (any expression that can be stream inserted)
 */
#define LVARRAY_ERROR_IF_GE_MSG( lhs, rhs, msg ) LVARRAY_ERROR_IF_OP_MSG( lhs, >=, <, rhs, msg )

/**
 * @brief Throw an exception if one value compares greater than or equal to the other.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param msg a message to log (any expression that can be stream inserted)
 * @param TYPE the type of exception to throw
 */
#define LVARRAY_THROW_IF_GE_MSG( lhs, rhs, msg, TYPE ) LVARRAY_THROW_IF_OP_MSG( lhs, >=, <, rhs, msg, TYPE )

/**
 * @brief Raise a hard error if one value compares greater than or equal to the other.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 */
#define LVARRAY_ERROR_IF_GE( lhs, rhs ) LVARRAY_ERROR_IF_GE_MSG( lhs, rhs, "" )

/**
 * @brief Throw an exception if one value compares greater than or equal to the other.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param TYPE the type of exception to throw
 */
#define LVARRAY_THROW_IF_GE( lhs, rhs, TYPE ) LVARRAY_THROW_IF_GE_MSG( lhs, rhs, "", TYPE )

/**
 * @brief Raise a hard error if one value compares less than the other.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param msg a message to log (any expression that can be stream inserted)
 */
#define LVARRAY_ERROR_IF_LT_MSG( lhs, rhs, msg ) LVARRAY_ERROR_IF_OP_MSG( lhs, <, >=, rhs, msg )

/**
 * @brief Throw an exception if one value compares less than the other.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param msg a message to log (any expression that can be stream inserted)
 * @param TYPE the type of exception to throw
 */
#define LVARRAY_THROW_IF_LT_MSG( lhs, rhs, msg, TYPE ) LVARRAY_THROW_IF_OP_MSG( lhs, <, >=, rhs, msg, TYPE )

/**
 * @brief Raise a hard error if one value compares less than the other.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 */
#define LVARRAY_ERROR_IF_LT( lhs, rhs ) LVARRAY_ERROR_IF_LT_MSG( lhs, rhs, "" )

/**
 * @brief Throw an exception if one value compares less than the other.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param TYPE the type of exception to throw
 */
#define LVARRAY_THROW_IF_LT( lhs, rhs, TYPE ) LVARRAY_THROW_IF_LT_MSG( lhs, rhs, "", TYPE )

/**
 * @brief Raise a hard error if one value compares less than or equal to the other.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param msg a message to log (any expression that can be stream inserted)
 */
#define LVARRAY_ERROR_IF_LE_MSG( lhs, rhs, msg ) LVARRAY_ERROR_IF_OP_MSG( lhs, <=, >, rhs, msg )

/**
 * @brief Throw an exception if one value compares less than or equal to the other.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param msg a message to log (any expression that can be stream inserted)
 * @param TYPE the type of exception to throw
 */
#define LVARRAY_THROW_IF_LE_MSG( lhs, rhs, msg, TYPE ) LVARRAY_THROW_IF_OP_MSG( lhs, <=, >, rhs, msg, TYPE )

/**
 * @brief Raise a hard error if one value compares less than or equal to the other.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 */
#define LVARRAY_ERROR_IF_LE( lhs, rhs ) LVARRAY_ERROR_IF_GE_MSG( lhs, rhs, "" )

/**
 * @brief Throw an exception if one value compares less than or equal to the other.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param TYPE the type of exception to throw
 */
#define LVARRAY_THROW_IF_LE( lhs, rhs, TYPE ) LVARRAY_THROW_IF_GE_MSG( lhs, rhs, "", TYPE )

/**
 * @brief Abort execution if @p lhs @p OP @p rhs is false.
 * @param lhs The left side of the operation.
 * @param OP The operation to apply.
 * @param rhs The right side of the operation.
 * @param msg The message to diplay.
 */
#define LVARRAY_ASSERT_OP_MSG( lhs, OP, rhs, msg ) \
  LVARRAY_ASSERT_MSG( lhs OP rhs, \
                      msg << "\n" << \
                      "  " << #lhs << " = " << lhs << "\n" << \
                      "  " << #rhs << " = " << rhs << "\n" )

/**
 * @brief Assert that two values compare equal in debug builds.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param msg a message to log (any expression that can be stream inserted)
 */
#define LVARRAY_ASSERT_EQ_MSG( lhs, rhs, msg ) LVARRAY_ASSERT_OP_MSG( lhs, ==, rhs, msg )

/**
 * @brief Assert that two values compare equal in debug builds.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 */
#define LVARRAY_ASSERT_EQ( lhs, rhs ) LVARRAY_ASSERT_EQ_MSG( lhs, rhs, "" )

/**
 * @brief Assert that two values compare not equal in debug builds.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param msg a message to log (any expression that can be stream inserted)
 */
#define LVARRAY_ASSERT_NE_MSG( lhs, rhs, msg ) LVARRAY_ASSERT_OP_MSG( lhs, !=, rhs, msg )

/**
 * @brief Assert that two values compare not equal in debug builds.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 */
#define LVARRAY_ASSERT_NE( lhs, rhs ) LVARRAY_ASSERT_NE_MSG( lhs, rhs, "" )

/**
 * @brief Assert that one value compares greater than the other in debug builds.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param msg a message to log (any expression that can be stream inserted)
 */
#define LVARRAY_ASSERT_GT_MSG( lhs, rhs, msg ) LVARRAY_ASSERT_OP_MSG( lhs, >, rhs, msg )

/**
 * @brief Assert that one value compares greater than the other in debug builds.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 */
#define LVARRAY_ASSERT_GT( lhs, rhs ) LVARRAY_ASSERT_GT_MSG( lhs, rhs, "" )

/**
 * @brief Assert that one value compares greater than or equal to the other in debug builds.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param msg a message to log (any expression that can be stream inserted)
 */
#define LVARRAY_ASSERT_GE_MSG( lhs, rhs, msg ) LVARRAY_ASSERT_OP_MSG( lhs, >=, rhs, msg )

/**
 * @brief Assert that one value compares greater than or equal to the other in debug builds.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 */
#define LVARRAY_ASSERT_GE( lhs, rhs ) LVARRAY_ASSERT_GE_MSG( lhs, rhs, "" )

#if defined(LVARRAY_DECORATE)
/// Mark a function for both host and device usage.
#define LVARRAY_HOST_DEVICE __host__ __device__

#if defined( LVARRAY_USE_HIP )
/// Mark a function for both host and device usage when using HIP only.
#define LVARRAY_HOST_DEVICE_HIP __host__ __device__
#else
/// Mark a function for both host and device usage when using HIP only.
#define LVARRAY_HOST_DEVICE_HIP
#endif

/// Mark a function for only device usage.
#define LVARRAY_DEVICE __device__

/**
 * @brief Disable host device warnings.
 * @details This pragma disables nvcc warnings about calling a host function from a host-device
 *   function. This is used on templated host-device functions where some template instantiations
 *   call host only code. This is safe as long as the host only instantiations are only called on
 *   the host. To use place directly above a the template.
 */
#if defined(LVARRAY_USE_CUDA)
#define DISABLE_HD_WARNING _Pragma("hd_warning_disable")
#else
#define DISABLE_HD_WARNING
#endif
#else
/// Mark a function for both host and device usage.
#define LVARRAY_HOST_DEVICE
/// Mark a function for both host and device usage when using HIP only.
#define LVARRAY_HOST_DEVICE_HIP

/// Mark a function for only device usage.
#define LVARRAY_DEVICE

/**
 * @brief Disable host device warnings.
 * @details This pragma disables nvcc warnings about calling a host function from a host-device
 *   function. This is used on templated host-device functions where some template instantiations
 *   call host only code. This is safe as long as the host only instantiations are only called on
 *   the host. To use place directly above a the template.
 */
#define DISABLE_HD_WARNING
#endif


#if defined(__clang__)
#define LVARRAY_RESTRICT __restrict__
#define LVARRAY_RESTRICT_REF __restrict__
#define LVARRAY_INTEL_CONSTEXPR constexpr
#elif defined(__GNUC__)
  #if defined(__INTEL_COMPILER)
#define LVARRAY_RESTRICT __restrict__
#define LVARRAY_RESTRICT_REF __restrict__
#define LVARRAY_INTEL_CONSTEXPR
  #else
#define LVARRAY_RESTRICT __restrict__
#define LVARRAY_RESTRICT_REF __restrict__
#define LVARRAY_INTEL_CONSTEXPR constexpr
  #endif
#endif

#if !defined(LVARRAY_BOUNDS_CHECK)
/**
 * @brief Expands to constexpr when array bound checking is disabled.
 */
#define CONSTEXPR_WITHOUT_BOUNDS_CHECK constexpr
#else
/**
 * @brief Expands to constexpr when array bound checking is disabled.
 */
#define CONSTEXPR_WITHOUT_BOUNDS_CHECK
#endif

#if defined(NDEBUG)
/**
 * @brief Expands to constexpr in release builds (when NDEBUG is defined).
 */
#define CONSTEXPR_WITH_NDEBUG constexpr
#else
/**
 * @brief Expands to constexpr in release builds (when NDEBUG is defined).
 */
#define CONSTEXPR_WITH_NDEBUG
#endif

#if !defined(LVARRAY_BOUNDS_CHECK)
/**
 * @brief Expands to constexpr when array bound checking is disabled.
 */
#define CONSTEXPR_WITHOUT_BOUNDS_CHECK constexpr
#else
/**
 * @brief Expands to constexpr when array bound checking is disabled.
 */
#define CONSTEXPR_WITHOUT_BOUNDS_CHECK
#endif

#if defined(NDEBUG)
/**
 * @brief Expands to constexpr in release builds (when NDEBUG is defined).
 */
#define CONSTEXPR_WITH_NDEBUG constexpr
#else
/**
 * @brief Expands to constexpr in release builds (when NDEBUG is defined).
 */
#define CONSTEXPR_WITH_NDEBUG
#endif

// TPL includes
#include <RAJA/RAJA.hpp>

template< typename >
struct RAJAHelper
{};

using serialPolicy = RAJA::seq_exec;

template<>
struct RAJAHelper< serialPolicy >
{
  using ReducePolicy = RAJA::seq_reduce;
  using AtomicPolicy = RAJA::seq_atomic;
};

#if defined(RAJA_ENABLE_OPENMP)

using parallelHostPolicy = RAJA::omp_parallel_for_exec;

template<>
struct RAJAHelper< parallelHostPolicy >
{
  using ReducePolicy = RAJA::omp_reduce;
  using AtomicPolicy = RAJA::omp_atomic;
};

#endif

#if defined(LVARRAY_USE_CUDA)

template< unsigned long THREADS_PER_BLOCK >
using parallelDevicePolicy = RAJA::cuda_exec< THREADS_PER_BLOCK >;

template< unsigned long N >
struct RAJAHelper< RAJA::cuda_exec< N > >
{
  using ReducePolicy = RAJA::cuda_reduce;
  using AtomicPolicy = RAJA::cuda_atomic;
};

#elif defined(LVARRAY_USE_HIP)

template< unsigned long THREADS_PER_BLOCK >
using parallelDevicePolicy = RAJA::hip_exec< THREADS_PER_BLOCK >;

template< unsigned long N >
struct RAJAHelper< RAJA::hip_exec< N > >
{
  using ReducePolicy = RAJA::hip_reduce;
  using AtomicPolicy = RAJA::hip_atomic;
};

#endif
