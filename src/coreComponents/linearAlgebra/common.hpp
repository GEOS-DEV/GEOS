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
 * @file common.hpp
 */
#ifndef GEOSX_LINEARALGEBRA_INTERFACES_COMMON_HPP_
#define GEOSX_LINEARALGEBRA_INTERFACES_COMMON_HPP_

#include "common/DataTypes.hpp"

/// Whether to check preconditions at runtime in LAI functions
#define GEOSX_LAI_RUNTIME_ASSERT 1

#if GEOSX_LAI_RUNTIME_ASSERT
#define GEOSX_LAI_ASSERT( expr ) GEOSX_ERROR_IF( !(expr), "" )
#define GEOSX_LAI_ASSERT_MSG( expr, msg ) GEOSX_ERROR_IF( !(expr), msg )
#define GEOSX_LAI_ASSERT_EQ( lhs, rhs ) GEOSX_ERROR_IF_NE( lhs, rhs )
#define GEOSX_LAI_ASSERT_NE( lhs, rhs ) GEOSX_ERROR_IF_EQ( lhs, rhs )
#define GEOSX_LAI_ASSERT_GT( lhs, rhs ) GEOSX_ERROR_IF_GE( rhs, lhs )
#define GEOSX_LAI_ASSERT_GE( lhs, rhs ) GEOSX_ERROR_IF_GT( rhs, lhs )
#else
#define GEOSX_LAI_ASSERT( expr ) GEOSX_ASSERT( expr )
#define GEOSX_LAI_ASSERT_MSG( expr, msg ) GEOSX_ASSERT_MSG( expr, msg )
#define GEOSX_LAI_ASSERT_EQ( lhs, rhs ) GEOSX_ASSERT_EQ( lhs, rhs )
#define GEOSX_LAI_ASSERT_NE( lhs, rhs ) GEOSX_ASSERT_NE( lhs, rhs )
#define GEOSX_LAI_ASSERT_GT( lhs, rhs ) GEOSX_ASSERT_GT( lhs, rhs )
#define GEOSX_LAI_ASSERT_GE( lhs, rhs ) GEOSX_ASSERT_GE( lhs, rhs )
#endif

/// Macro for checking and reporting error codes from TPL packages
#define GEOSX_LAI_CHECK_ERROR( call ) \
  do { \
    auto const ierr = call; \
    GEOSX_ERROR_IF_NE_MSG( ierr, 0, "Error in call to " << #call ); \
  } while( false )

#define GEOSX_LAI_CHECK_ERROR_NNEG( call ) \
  do { \
    auto const ierr = call; \
    GEOSX_ERROR_IF_GT_MSG( 0, ierr, "Error in call to " << #call ); \
  } while( false )

namespace geosx
{

enum class LAIOutputFormat
{
  NATIVE_ASCII,
  NATIVE_BINARY,
  MATLAB_ASCII,
  MATLAB_BINARY,
  MATRIX_MARKET
};

struct MatrixLayout
{
  using ROW_MAJOR_PERM = RAJA::PERM_IJ;
  using COL_MAJOR_PERM = RAJA::PERM_JI;

  constexpr static int const ROW_MAJOR = LvArray::getStrideOneDimension( ROW_MAJOR_PERM {} );
  constexpr static int const COL_MAJOR = LvArray::getStrideOneDimension( COL_MAJOR_PERM {} );
};

}

#endif //GEOSX_LINEARALGEBRA_INTERFACES_COMMON_HPP_
