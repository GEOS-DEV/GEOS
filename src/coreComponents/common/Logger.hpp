/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file Logger.hpp
 */

#ifndef GEOSX_COMMON_LOGGER_HPP
#define GEOSX_COMMON_LOGGER_HPP

// Source incldes
#include "common/GeosxConfig.hpp"
#include "cxx-utilities/src/Macros.hpp"

// System includes
#if defined(GEOSX_USE_MPI)
  #include <mpi.h>
#endif

#define GEOSX_LOG( msg ) LVARRAY_LOG( msg )

#define GEOSX_LOG_VAR( var ) GEOSX_LOG( #var << " = " << var )

#define GEOSX_LOG_RANK_0_IF( EXP, msg ) \
  do { \
    if( logger::internal::rank == 0 && EXP ) \
    { \
      std::ostringstream oss; \
      oss << msg; \
      std::cout << oss.str() << std::endl; \
    } \
  } while( false )

#define GEOSX_LOG_RANK_0( msg ) GEOSX_LOG_RANK_0_IF( true, msg )

#define GEOSX_LOG_RANK_IF( EXP, msg ) \
  do { \
    if( EXP ) \
    { \
      std::ostringstream oss; \
      oss << "Rank " << logger::internal::rank << ": " << msg; \
      *logger::internal::rankStream << oss.str() << std::endl; \
    } \
  } while( false )

#define GEOSX_LOG_RANK( msg ) GEOSX_LOG_RANK_IF( true, msg )

#define GEOSX_LOG_RANK_VAR( var ) GEOSX_LOG_RANK( #var " = " << var )

#define GEOSX_ERROR_IF( EXP, msg ) LVARRAY_ERROR_IF( EXP, msg )
#define GEOSX_ERROR( msg ) LVARRAY_ERROR( msg )

#define GEOSX_ASSERT_MSG( EXP, msg ) LVARRAY_ASSERT_MSG( EXP, msg )
#define GEOSX_ASSERT( EXP ) LVARRAY_ASSERT( EXP )

#define GEOSX_WARNING_IF( EXP, msg ) LVARRAY_WARNING_IF( EXP, msg )
#define GEOSX_WARNING( msg ) LVARRAY_WARNING( msg )

#define GEOSX_INFO_IF( EXP, msg ) LVARRAY_INFO_IF( EXP, msg )
#define GEOSX_INFO( msg ) LVARRAY_INFO( msg )

#define GEOSX_CHECK( EXP, msg ) LVARRAY_CHECK( EXP, msg )

#define GEOSX_ERROR_IF_EQ_MSG( lhs, rhs, msg ) LVARRAY_ERROR_IF_EQ_MSG( lhs, rhs, msg )
#define GEOSX_ERROR_IF_EQ( lhs, rhs ) GEOSX_ERROR_IF_EQ( lhs, rhs )

#define GEOSX_ERROR_IF_NE_MSG( lhs, rhs, msg ) LVARRAY_ERROR_IF_NE_MSG( lhs, rhs, msg )
#define GEOSX_ERROR_IF_NE( lhs, rhs ) LVARRAY_ERROR_IF_NE( lhs, rhs )

#define GEOSX_ERROR_IF_GT_MSG( lhs, rhs, msg ) LVARRAY_ERROR_IF_GT_MSG( lhs, rhs, msg )
#define GEOSX_ERROR_IF_GT( lhs, rhs ) LVARRAY_ERROR_IF_GT( lhs, rhs )

#define GEOSX_ERROR_IF_GE_MSG( lhs, rhs, msg ) LVARRAY_ERROR_IF_GE_MSG( lhs, rhs, msg )
#define GEOSX_ERROR_IF_GE( lhs, rhs ) LVARRAY_ERROR_IF_GE( lhs, rhs )

#define GEOSX_ERROR_IF_LT_MSG( lhs, rhs, msg ) LVARRAY_ERROR_IF_LT_MSG( lhs, rhs, msg )
#define GEOSX_ERROR_IF_LT( lhs, rhs ) LVARRAY_ERROR_IF_LT( lhs, rhs )

#define GEOSX_ERROR_IF_LE_MSG( lhs, rhs, msg ) LVARRAY_ERROR_IF_LE_MSG( lhs, rhs, msg )
#define GEOSX_ERROR_IF_LE( lhs, rhs ) LVARRAY_ERROR_IF_LE( lhs, rhs )

#define GEOSX_ASSERT_EQ_MSG( lhs, rhs, msg ) LVARRAY_ASSERT_EQ_MSG( lhs, rhs, msg )
#define GEOSX_ASSERT_EQ( lhs, rhs ) LVARRAY_ASSERT_EQ( lhs, rhs )

#define GEOSX_ASSERT_GT_MSG( lhs, rhs, msg ) LVARRAY_ASSERT_GT_MSG( lhs, rhs, msg )
#define GEOSX_ASSERT_GT( lhs, rhs ) LVARRAY_ASSERT_GT( lhs, rhs )

#define GEOSX_ASSERT_GE_MSG( lhs, rhs, msg ) LVARRAY_ASSERT_GE_MSG( lhs, rhs, msg )
#define GEOSX_ASSERT_GE( lhs, rhs ) LVARRAY_ASSERT_GE( lhs, rhs )


/**
 * @brief Macro used to turn on/off a function based on the log level.
 * @param[in] minLevel Minimum log level
 * @param[in] fn Function to filter
 */
#define GEOSX_LOG_LEVEL_FN( minLevel, fn )                                      \
  do {                                                                         \
    if( this->getLogLevel() >= minLevel )                                      \
    {                                                                          \
      fn;                                                                      \
    }                                                                          \
  } while( false )

/**
 * @brief Macro used to output messages based on the log level.
 * @param[in] minLevel Minimum log level
 * @param[in] msg Log message
 */
#define GEOSX_LOG_LEVEL( minLevel, msg ) GEOSX_INFO_IF( this->getLogLevel() >= minLevel, msg );

/**
 * @brief Macro used to output messages (only on rank 0) based on the log level.
 * @param[in] minLevel Minimum log level
 * @param[in] msg Log message
 */
#define GEOSX_LOG_LEVEL_RANK_0( minLevel, msg ) GEOSX_LOG_RANK_0_IF( this->getLogLevel() >= minLevel, msg )

/**
 * @brief Macro used to output messages (with one line per rank) based on the log level.
 * @param[in] minLevel Minimum log level
 * @param[in] msg Log message
 */
#define GEOSX_LOG_LEVEL_BY_RANK( minLevel, msg ) GEOSX_LOG_RANK_IF( this->getLogLevel() >= minLevel, msg )


namespace geosx
{

namespace logger
{

namespace internal
{

extern int rank;

extern int n_ranks;

extern std::ostream * rankStream;

#if defined(GEOSX_USE_MPI)
extern MPI_Comm comm;
#endif
} // namespace internal

#if defined(GEOSX_USE_MPI)
void InitializeLogger( MPI_Comm comm, const std::string & rank_output_dir="" );
#endif

void InitializeLogger( const std::string & rank_output_dir="" );

void FinalizeLogger();

} // namespace logger

} // namespace geosx

#endif /* GEOSX_COMMON_LOGGER_HPP */
