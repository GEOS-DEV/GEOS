/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * @file testDenseLAOperations.cpp
 */

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wglobal-constructors"
#pragma clang diagnostic ignored "-Wexit-time-destructors"
#if __clang_major__ >= 5 && !defined(__APPLE__)
#pragma clang diagnostic ignored "-Wzero-as-null-pointer-constant"
#endif
#endif

#include <gtest/gtest.h>

#include "common/DataTypes.hpp"

#include "LapackInterface.hpp"


/**
 * \file testDenseLAOperations.cpp
 * \brief This test file is part of the ctest suite and tests Lapack- and
 * MAGMA-based dense linear algebra operations.
 */

using namespace geosx;

/*! @name Ctest tests.
 * @brief Runs similar testing functions using different Linear Algebra Interfaces (LAIs).
 */

/*! @name Utility functions.
 * @brief Functions used to construct useful matrices in the test files.
 */
//@{

/**
 * @brief Compute an identity matrix
 *
 * \param comm MPI communicator.
 * \param N global size of the square identity matrix.
 */

// BEGIN_RST_NARRATIVE testDenseLAOperations.rst

// ==============================
// Test Linear Algebra Operations
// ==============================
// In these 3 functions we test the linear algebra operations, the native solvers from the
// libraries as well as the re-implemented GEOSX solvers for CG and BiCGSTAB. We run these
// on both monolithic and block matrices.

// -------------------------------------
// Test libraries operations and solvers
// -------------------------------------
// We start by testing the linear algebra operations. We fill two matrices (one will be a
// preconditioner) and make sure the sparse storage is behaving properly. We then test the
// iterative and direct solvers available.

//------------------------------
// Test matrix constructor
//------------------------------
template< typename LAI >
void testMatrixConstructors()
{
  using DenseMatrix = typename LAI::DenseMatrix;

  // Empty constructor

  DenseMatrix A(2,2);
  DenseMatrix B(2,2);

//  SerialDenseMatrix A5(3,3);
//  SerialDenseMatrix B(3,2);
//  SerialDenseMatrix C(3,2);
//  SerialDenseMatrix D(2,2);
//  SerialDenseMatrix E(2,3);
//  SerialDenseMatrix F(1,1);

  EXPECT_EQ( A.get_nRows(), 2);
  EXPECT_EQ( B.get_nRows(), 3);
}


// END_RST_NARRATIVE

//@}


#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wglobal-constructors"
#pragma clang diagnostic ignored "-Wexit-time-destructors"
#if __clang_major__ >= 5 && !defined(__APPLE__)
#pragma clang diagnostic ignored "-Wzero-as-null-pointer-constant"
#endif
#endif

/*! @name Ctest tests.
 * @brief Runs similar testing functions using different Linear Algebra Interfaces (LAIs).
 */
//@{



/*! @function testLapackDenseLAOperations.
 * @brief Runs all tests using the Lapack interface.
 */
TEST(testDenseLAOperations,testLapackDenseLAOperations)
{
  testMatrixConstructors<LapackInterface>();
}

//@}

#ifdef __clang__
#pragma clang diagnostic pop
#endif
