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

#include <BlasLapackLA.hpp>

#include <random>

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
template<typename LAI>
void testArray1dLA()
{

  BlasLapackLA denseLA;

  real64 machinePrecision = 10. * std::numeric_limits<real64>::epsilon();

  // Repeat the following step for vectors of increasing size:
  //
  // a. Resize vectors v1, v2, and v3 to be of equal size
  //    (not initialized to 0)
  // b. Fill v1 with random coefficient
  // c. Initialize v2 to have zero coefficients
  // d. Compute v2 = v2 + alfa*v1
  // e. Compute v1 = alfa*v1
  // f. Copy v1 in v3
  // g. Compute v1 = v1 - v2
  // h. Check that (norm(v1, type) / norm(v2, type)) < epsilon
  //    with type = \{1-norm, 2-norm, Infinity-norm)
  // i. check that ( (norm(v2, norm-2) - sqrt(v2_dot_v2)) / norm(v2, norm-2)) < epsilon

  array1d<real64> v1, v2, v3;
  std::default_random_engine generator;
  std::uniform_real_distribution< real64 > distribution( 0.0, 1.0 );
  real64 alfa = distribution( generator );

  for( int iSize = 1 ; iSize <= 10 ; ++iSize )
  {
    // a.
    v1.resize( iSize );
    v2.resize( iSize );
    v3.resize( iSize );

    // b.
    for (int i = 0; i < iSize; ++i)
    {
      v1[i] = distribution(generator);
    }

    // c.
    v2 = 0.0;

    // d.
    denseLA.vectorVectorAdd(v1, v2, alfa );

    // e.
    denseLA.vectorScale(v1, alfa );

    // f.
    denseLA.vectorCopy(v1, v3);

    // g.
    denseLA.vectorVectorAdd(v2, v1, -1.0 );

    // h.
    EXPECT_NEAR( denseLA.vectorNorm1(v1),
                 denseLA.vectorNorm1(v2),
                 machinePrecision*denseLA.vectorNorm1(v2));
    EXPECT_NEAR( denseLA.vectorNorm2(v1),
                 denseLA.vectorNorm2(v2),
                 machinePrecision*denseLA.vectorNorm2(v2) );
    EXPECT_NEAR( denseLA.vectorNormInf(v1),
                 denseLA.vectorNormInf(v2),
                 machinePrecision*denseLA.vectorNormInf(v2) );

    // i.
    real64 beta = denseLA.vectorDot(v2, v2);
    EXPECT_NEAR( std::sqrt( denseLA.vectorDot(v2, v2)),
                 denseLA.vectorNorm2(v2),
                 machinePrecision*denseLA.vectorNorm2(v2) );

  }

}

template<typename LAI>
void testArray2dLA()
{
  BlasLapackLA denseLA;

  real64 machinePrecision = 10. * std::numeric_limits<real64>::epsilon();

  array2d<real64> A, B, C, D, LHS, RHS;
  std::default_random_engine generator;
  std::uniform_real_distribution<double> distribution( 0.0, 1.0 );
  real64 alfa = distribution( generator );
  real64 beta = distribution( generator );

  localIndex MA, NA, MB, NB, ND, K;

  // Test 2: repeat the following step for vectors of increasing size:
  //
  // a. compute LHS = ( alfa*A*B + beta*C ) * D
  // b. compute RHS = alfa*A*(B*D) + beta*C*D
  // c. check that (norm(LHS - RHS, type) / norm(LHS, type)) < epsilon
  //    with type = \{Infinity-norm, 1-norm, Frobenius-norm, )

  MA = 10;
  NB = 10;
  ND = 10;
  K = 20;

  for( localIndex mA = 1 ; mA <= MA ; ++mA )
    for( localIndex nB = 1 ; nB <= NB ; ++nB )
      for( localIndex nD = 1 ; nD <= ND ; ++nD )
        for( localIndex k = 1 ; k <= K ; ++k )
        {
          // Resize matrix operators
          A.resize( mA, k );
          B.resize( k, nB );
          C.resize( A.size(0), B.size(1) );
          D.resize( B.size(1), nD );
          LHS.resize( A.size(0), D.size(1) );
          RHS.resize( A.size(0), D.size(1) );

          // Populate A, B, C, and D with uniformly distributed random
          // coefficients
//          A.rand();
//          B.rand();
//          C.rand();
//          D.rand();
//
//          // Compute tmp1 = ( alfa*A*B + beta*C )
//          SerialMatrix tmp1( C );
//          A.matrixMultiply( B, tmp1, alfa, beta );
//
//          // Compute LHS = tmp*D
//          tmp1.matrixMultiply( D, LHS );
//
//          // Compute tmp2 = B*D
//          SerialMatrix tmp2( B.getNumRows(), D.getNumCols() );
//          B.matrixMultiply( D, tmp2 );
//
//          // Compute RHS = alfa*A*tmp2
//          A.matrixMultiply( tmp2, RHS, alfa );
//
//          // Compute RHS = RHS + beta*C*D
//          C.matrixMultiply( D, RHS, beta, 1. );
//
//          // Check norms
//          RHS.matrixAdd( LHS, -1. );
//          EXPECT_LT( RHS.normInf() / LHS.normInf(), machinePrecision );
//          EXPECT_LT( RHS.norm1() / LHS.norm1(), machinePrecision );
//          EXPECT_LT( RHS.normFrobenius() / LHS.normFrobenius(), machinePrecision );
//
        }
//
//  // Test 3: repeat the following step for vectors of increasing size:
//  //
//  // a. compute LHS = ( alfa*A^T*B + beta*C ) * D
//  // b. compute RHS = alfa*A^T*(B*D) + beta*C*D
//  // c. check that (norm(LHS - RHS, type) / norm(LHS, type)) < epsilon
//  //    with type = \{Infinity-norm, 1-norm, Frobenius-norm, )
//
//  NA = 10;
//  NB = 10;
//  ND = 10;
//  K = 20;
//
//  for( localIndex nA = 1 ; nA <= NA ; ++nA )
//    for( localIndex nB = 1 ; nB <= NB ; ++nB )
//      for( localIndex nD = 1 ; nD <= ND ; ++nD )
//        for( localIndex k = 1 ; k <= K ; ++k )
//        {
//
//          // Resize matrix operators
//          A.resize( k, nA );
//          B.resize( k, nB );
//          C.resize( A.getNumCols(), B.getNumCols() );
//          D.resize( B.getNumCols(), nD );
//          LHS.resize( A.getNumCols(), D.getNumCols() );
//          RHS.resize( A.getNumCols(), D.getNumCols() );
//
//          // Populate A, B, C, and D with uniformly distributed random
//          // coefficients
//          A.rand();
//          B.rand();
//          C.rand();
//          D.rand();
//
//          // Compute tmp1 = ( alfa*A^T*B + beta*C )
//          SerialMatrix tmp1( C );
//          A.TmatrixMultiply( B, tmp1, alfa, beta );
//
//          // Compute LHS = tmp*D
//          tmp1.matrixMultiply( D, LHS );
//
//          // Compute tmp2 = B*D
//          SerialMatrix tmp2( B.getNumRows(), D.getNumCols() );
//          B.matrixMultiply( D, tmp2 );
//
//          // Compute RHS = alfa*A^T*tmp2
//          A.TmatrixMultiply( tmp2, RHS, alfa );
//
//          // Compute RHS = RHS + beta*C*D
//          C.matrixMultiply( D, RHS, beta, 1. );
//
//          // Check norms
//
//          RHS.matrixAdd( LHS, -1. );
//          EXPECT_LT( RHS.normInf() / LHS.normInf(), machinePrecision );
//          EXPECT_LT( RHS.norm1() / LHS.norm1(), machinePrecision );
//          EXPECT_LT( RHS.normFrobenius() / LHS.normFrobenius(), machinePrecision );
//
//        }
//
//  // Test 3: repeat the following step for vectors of increasing size:
//  //
//  // a. compute LHS = ( alfa*A*B^T + beta*C ) * D
//  // b. compute RHS = alfa*A*(B^T*D) + beta*C*D
//  // c. check that (norm(LHS - RHS, type) / norm(LHS, type)) < epsilon
//  //    with type = \{Infinity-norm, 1-norm, Frobenius-norm, )
//
//  MA = 10;
//  MB = 10;
//  ND = 10;
//  K = 20;
//
//  for( localIndex mA = 1 ; mA <= MA ; ++mA )
//    for( localIndex mB = 1 ; mB <= MB ; ++mB )
//      for( localIndex nD = 1 ; nD <= ND ; ++nD )
//        for( localIndex k = 1 ; k <= K ; ++k )
//        {
//
//          // Resize matrix operators
//          A.resize( mA, k );
//          B.resize( mB, k );
//          C.resize( A.getNumRows(), B.getNumRows() );
//          D.resize( B.getNumRows(), nD );
//          LHS.resize( A.getNumRows(), D.getNumCols() );
//          RHS.resize( A.getNumRows(), D.getNumCols() );
//
//          // Populate A, B, C, and D with uniformly distributed random
//          // coefficients
//          A.rand();
//          B.rand();
//          C.rand();
//          D.rand();
//
//          // Compute tmp1 = ( alfa*A*B^T + beta*C )
//          SerialMatrix tmp1( C );
//          A.matrixTMultiply( B, tmp1, alfa, beta );
//
//          // Compute LHS = tmp*D
//          tmp1.matrixMultiply( D, LHS );
//
//          // Compute tmp2 = B^T*D
//          SerialMatrix tmp2( B.getNumCols(), D.getNumCols() );
//          B.TmatrixMultiply( D, tmp2 );
//
//          // Compute RHS = alfa*A*tmp2
//          A.matrixMultiply( tmp2, RHS, alfa );
//
//          // Compute RHS = RHS + beta*C*D
//          C.matrixMultiply( D, RHS, beta, 1. );
//
//          // Check norms
//          RHS.matrixAdd( LHS, -1. );
//          EXPECT_LT( RHS.normInf() / LHS.normInf(), machinePrecision );
//          EXPECT_LT( RHS.norm1() / LHS.norm1(), machinePrecision );
//          EXPECT_LT( RHS.normFrobenius() / LHS.normFrobenius(), machinePrecision );
//
//        }
//
//  // Test 4: repeat the following step for vectors of increasing size:
//  //
//  // a. compute LHS = ( alfa*A^T*B^T + beta*C ) * D
//  // b. compute RHS = alfa*A^T*(B^T*D) + beta*C*D
//  // c. check that (norm(LHS - RHS, type) / norm(LHS, type)) < epsilon
//  //    with type = \{Infinity-norm, 1-norm, Frobenius-norm, )
//
//  NA = 10;
//  MB = 10;
//  ND = 10;
//  K = 20;
//
//  for( localIndex nA = 1 ; nA <= NA ; ++nA )
//    for( localIndex mB = 1 ; mB <= MB ; ++mB )
//      for( localIndex nD = 1 ; nD <= ND ; ++nD )
//        for( localIndex k = 1 ; k <= K ; ++k )
//        {
//
//          // Resize matrix operators
//          A.resize( k, nA );
//          B.resize( mB, k );
//          C.resize( A.getNumCols(), B.getNumRows() );
//          D.resize( B.getNumRows(), nD );
//          LHS.resize( A.getNumCols(), D.getNumCols() );
//          RHS.resize( A.getNumCols(), D.getNumCols() );
//
//          // Populate A, B, C, and D with uniformly distributed random
//          // coefficients
//          A.rand();
//          B.rand();
//          C.rand();
//          D.rand();
//
//          // Compute tmp1 = ( alfa*A^T*B^T + beta*C )
//          SerialMatrix tmp1( C );
//          A.TmatrixTMultiply( B, tmp1, alfa, beta );
//
//          // Compute LHS = tmp*D
//          tmp1.matrixMultiply( D, LHS );
//
//          // Compute tmp2 = B^T*D
//          SerialMatrix tmp2( B.getNumCols(), D.getNumCols() );
//          B.TmatrixMultiply( D, tmp2 );
//
//          // Compute RHS = alfa*A^T*tmp2
//          A.TmatrixMultiply( tmp2, RHS, alfa );
//
//          // Compute RHS = RHS + beta*C*D
//          C.matrixMultiply( D, RHS, beta, 1. );
//
//          // Check norms
//          RHS.matrixAdd( LHS, -1. );
//          EXPECT_LT( RHS.normInf() / LHS.normInf(), machinePrecision );
//          EXPECT_LT( RHS.norm1() / LHS.norm1(), machinePrecision );
//          EXPECT_LT( RHS.normFrobenius() / LHS.normFrobenius(), machinePrecision );
//
//        }
//
}

//template<typename LAI>
//void testSerialMatrixVector()
//{
//
//  using SerialMatrix = typename LAI::SerialMatrix;
//  using SerialVector = typename LAI::SerialVector;
//
//  real64 machinePrecision = 10. * std::numeric_limits<real64>::epsilon();
//  SerialMatrix A, yT, tmp;
//  SerialVector x, y, yTT;
//
//  // Test 5: repeat the following step for vectors of increasing size:
//  //
//  // a. compute y = A*x
//  // b. compute compute yT = x^T * A^T
//  // c. check that (norm(y, 2-norm) / norm(yT, Frobenius-norm)) < epsilon
//  // d. create vector yTT = yT^T
//  // e. compute beta = sqrt(yTT dot y);
//  // f. check that ( (norm(y, norm-2) - sqrt(yTT_dot_y)) / norm(y, norm-2)) < epsilon
//
//  localIndex MA = 10;
//  localIndex NA = 10;
//  real64 alfa, beta;
//
//  for( localIndex mA = 1 ; mA <= MA ; ++mA )
//    for( localIndex nA = 1 ; nA <= NA ; ++nA )
//    {
//      // Resize matrices and vectors
//      A.resize( mA, nA );
//      x.resize( nA);
//      y.resize( mA);
//
//      // Populate matrix A and vector x with uniformly distributed random
//      // coefficients
//      A.rand();
//      x.rand();
//
//      // a.
//      A.vectorMultiply(x, y);
//
//      // b.
//      // --- construct tmp = x^T
//      tmp.resize(1, x.getSize());
//      for (localIndex i = 0; i < x.getSize(); ++i)
//        tmp(0,i) = x(i);
//
//      // -- compute yT = tmp * A^T
//      yT.resize(1, A.getNumRows());
//      tmp.matrixTMultiply(A, yT);
//
//      // c.
//      alfa = y.norm2();
//      beta = std::fabs(alfa - yT.normFrobenius());
//      EXPECT_LT( beta/alfa, machinePrecision );
//
//      // d.
//      yTT.resize(yT.getNumCols());
//      for (localIndex i = 0; i < yT.getNumCols(); ++i)
//        yTT(i) = yT(0,i);
//
//      // e.
//      beta = std::sqrt( y.dot(yTT) );
//
//      // f. check that ( (norm(y, norm-2) - sqrt(yTT_dot_y)) / norm(y, norm-2)) < epsilon
//      beta = alfa - beta;
//      EXPECT_LT( beta/alfa, machinePrecision );
//
//    }
//}

//template<typename LAI>
//void testSerialMatrixInverse()
//{
//
//  using SerialMatrix = typename LAI::SerialMatrix;
//
//  real64 machinePrecision = 10. * std::numeric_limits<real64>::epsilon();
//
//  SerialMatrix E;
//  SerialMatrix Einv;
//  SerialMatrix EinvXE;
//
//  // Test 5: repeat the following step for matrices of increasing size:
//  // a. Construct matrix E (1d discrete Laplacian)
//  // b. Compute Einv = E^-1
//  // c. Compute EinvXE = Einv*E
//  // d. Check that det(EinvXE) = 1.
//
//  real64 det;
//  localIndex max_dim = 10;
//
//  for( localIndex order = 1 ; order <= max_dim ; ++order )
//  {
//    // a.
//    E.resize( order );
//    for( localIndex i = 0 ; i < E.getNumCols() ; ++i )
//      for( localIndex j = 0 ; j < E.getNumRows() ; ++j )
//      {
//        if( i == j )
//          E( i, i ) = 2;
//        else if( abs( i - j ) == 1 )
//          E( i, j ) = -1;
//      }
//
//    // b.
//    E.computeInverse( Einv, det );
//
//    // c.
//    EinvXE.resize( E.getNumRows() );
//    Einv.matrixMultiply( E, EinvXE );
//
//    // d.
//    EXPECT_LT( std::fabs( EinvXE.determinant() - 1. ), machinePrecision );
//  }
//}

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

  testArray1dLA<BlasLapackLA>();
  testArray2dLA<BlasLapackLA>();
//  testSerialMatrixVector<LapackSuiteInterface>();
//  testSerialMatrixInverse<LapackSuiteInterface>();

}

//@}

#ifdef __clang__
#pragma clang diagnostic pop
#endif
