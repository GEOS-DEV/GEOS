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

#include <LapackSuiteInterface.hpp>

//#include <random>

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
void testSerialVector()
{

  using SerialVector = typename LAI::SerialVector;

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
  // i. Random forward permutations of v2 and v3
  // j. Check that abs(norm(v2, type) - norm(v3, type)) / norm(v2, type)) < epsilon
  //    with type = \{1-norm, 2-norm, Infinity-norm)
  // k. Apply inverse permutation to v2 and v3
  // l. Copy v3 in v1
  // m. Compute v1 = v1 - v2
  // n. Check that (norm(v1, type) / norm(v2, type)) < epsilon
  //    with type = \{1-norm, 2-norm, Infinity-norm)

  SerialVector v1, v2, v3;
  std::default_random_engine generator;
  std::uniform_real_distribution<double> distribution( 0.0, 1.0 );
  real64 alfa = distribution( generator );

  for( int iSize = 1 ; iSize <= 10 ; ++iSize )
  {
    // a.
    v1.resize( iSize );
    v2.resize( iSize );
    v3.resize( iSize );

    // b.
    v1.rand();

    // c.
    v2.zero();

    // d.
    v2.vectorAdd( v1, alfa );

    // e.
    v1.scale( alfa );

    // f.
    v3.copy( v1 );

    // g.
    v1.vectorAdd( v2, -1. );

    // h.
    EXPECT_LT( v1.norm1() / v2.norm1(), machinePrecision );
    EXPECT_LT( v1.norm2() / v2.norm2(), machinePrecision );
    EXPECT_LT( v1.normInf() / v2.normInf(), machinePrecision );

    // i.
    array1d<int> v2Perm( iSize ); //start counting from 1
    for( int i = 0 ; i < iSize ; ++i )
      v2Perm( i ) = i + 1;
    array1d<int> v3Perm( v2Perm );
    std::random_shuffle( v2Perm.begin(), v2Perm.end() );
    std::random_shuffle( v3Perm.begin(), v3Perm.end() );
    v2.permute( v2Perm );
    v3.permute( v3Perm );

    // j.
    EXPECT_LT( std::fabs( v2.norm1() - v3.norm1() ) / v2.norm1(), machinePrecision );
    EXPECT_LT( std::fabs( v2.norm2() - v3.norm2() ) / v2.norm2(), machinePrecision );
    EXPECT_LT( std::fabs( v2.normInf() - v3.normInf() ) / v2.normInf(), machinePrecision );

    // k.
    v2.permute( v2Perm, false );
    v3.permute( v3Perm, false );

    // l.
    v1.copy( v3 );

    // m .
    v1.vectorAdd( v2, -1. );

    // n.
    EXPECT_LT( v1.norm1() / v2.norm1(), machinePrecision );
    EXPECT_LT( v1.norm2() / v2.norm2(), machinePrecision );
    EXPECT_LT( v1.normInf() / v2.normInf(), machinePrecision );

  }

}

template<typename LAI>
void testSerialMatrix()
{

  using SerialMatrix = typename LAI::SerialMatrix;
//  using SerialVector = typename LAI::SerialVector;

  real64 machinePrecision = 10. * std::numeric_limits<real64>::epsilon();
  SerialMatrix A, B, C, D, LHS, RHS;
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
          C.resize( A.getNumRows(), B.getNumCols() );
          D.resize( B.getNumCols(), nD );
          LHS.resize( A.getNumRows(), D.getNumCols() );
          RHS.resize( A.getNumRows(), D.getNumCols() );

          // Populate A, B, C, and D with uniformly distributed random
          // coefficients
          A.rand();
          B.rand();
          C.rand();
          D.rand();

          // Compute tmp1 = ( alfa*A*B + beta*C )
          SerialMatrix tmp1( C );
          A.matrixMultiply( B, tmp1, alfa, beta );

          // Compute LHS = tmp*D
          tmp1.matrixMultiply( D, LHS );

          // Compute tmp2 = B*D
          SerialMatrix tmp2( B.getNumRows(), D.getNumCols() );
          B.matrixMultiply( D, tmp2 );

          // Compute RHS = alfa*A*tmp2
          A.matrixMultiply( tmp2, RHS, alfa );

          // Compute RHS = RHS + beta*C*D
          C.matrixMultiply( D, RHS, beta, 1. );

          // Check norms
          RHS.matrixAdd( LHS, -1. );
          EXPECT_LT( RHS.normInf() / LHS.normInf(), machinePrecision );
          EXPECT_LT( RHS.norm1() / LHS.norm1(), machinePrecision );
          EXPECT_LT( RHS.normFrobenius() / LHS.normFrobenius(), machinePrecision );

        }

  // Test 3: repeat the following step for vectors of increasing size:
  //
  // a. compute LHS = ( alfa*A^T*B + beta*C ) * D
  // b. compute RHS = alfa*A^T*(B*D) + beta*C*D
  // c. check that (norm(LHS - RHS, type) / norm(LHS, type)) < epsilon
  //    with type = \{Infinity-norm, 1-norm, Frobenius-norm, )

  NA = 10;
  NB = 10;
  ND = 10;
  K = 20;

  for( localIndex nA = 1 ; nA <= NA ; ++nA )
    for( localIndex nB = 1 ; nB <= NB ; ++nB )
      for( localIndex nD = 1 ; nD <= ND ; ++nD )
        for( localIndex k = 1 ; k <= K ; ++k )
        {

          // Resize matrix operators
          A.resize( k, nA );
          B.resize( k, nB );
          C.resize( A.getNumCols(), B.getNumCols() );
          D.resize( B.getNumCols(), nD );
          LHS.resize( A.getNumCols(), D.getNumCols() );
          RHS.resize( A.getNumCols(), D.getNumCols() );

          // Populate A, B, C, and D with uniformly distributed random
          // coefficients
          A.rand();
          B.rand();
          C.rand();
          D.rand();

          // Compute tmp1 = ( alfa*A^T*B + beta*C )
          SerialMatrix tmp1( C );
          A.TmatrixMultiply( B, tmp1, alfa, beta );

          // Compute LHS = tmp*D
          tmp1.matrixMultiply( D, LHS );

          // Compute tmp2 = B*D
          SerialMatrix tmp2( B.getNumRows(), D.getNumCols() );
          B.matrixMultiply( D, tmp2 );

          // Compute RHS = alfa*A^T*tmp2
          A.TmatrixMultiply( tmp2, RHS, alfa );

          // Compute RHS = RHS + beta*C*D
          C.matrixMultiply( D, RHS, beta, 1. );

          // Check norms

          RHS.matrixAdd( LHS, -1. );
          EXPECT_LT( RHS.normInf() / LHS.normInf(), machinePrecision );
          EXPECT_LT( RHS.norm1() / LHS.norm1(), machinePrecision );
          EXPECT_LT( RHS.normFrobenius() / LHS.normFrobenius(), machinePrecision );

        }

  // Test 3: repeat the following step for vectors of increasing size:
  //
  // a. compute LHS = ( alfa*A*B^T + beta*C ) * D
  // b. compute RHS = alfa*A*(B^T*D) + beta*C*D
  // c. check that (norm(LHS - RHS, type) / norm(LHS, type)) < epsilon
  //    with type = \{Infinity-norm, 1-norm, Frobenius-norm, )

  MA = 10;
  MB = 10;
  ND = 10;
  K = 20;

  for( localIndex mA = 1 ; mA <= MA ; ++mA )
    for( localIndex mB = 1 ; mB <= MB ; ++mB )
      for( localIndex nD = 1 ; nD <= ND ; ++nD )
        for( localIndex k = 1 ; k <= K ; ++k )
        {

          // Resize matrix operators
          A.resize( mA, k );
          B.resize( mB, k );
          C.resize( A.getNumRows(), B.getNumRows() );
          D.resize( B.getNumRows(), nD );
          LHS.resize( A.getNumRows(), D.getNumCols() );
          RHS.resize( A.getNumRows(), D.getNumCols() );

          // Populate A, B, C, and D with uniformly distributed random
          // coefficients
          A.rand();
          B.rand();
          C.rand();
          D.rand();

          // Compute tmp1 = ( alfa*A*B^T + beta*C )
          SerialMatrix tmp1( C );
          A.matrixTMultiply( B, tmp1, alfa, beta );

          // Compute LHS = tmp*D
          tmp1.matrixMultiply( D, LHS );

          // Compute tmp2 = B^T*D
          SerialMatrix tmp2( B.getNumCols(), D.getNumCols() );
          B.TmatrixMultiply( D, tmp2 );

          // Compute RHS = alfa*A*tmp2
          A.matrixMultiply( tmp2, RHS, alfa );

          // Compute RHS = RHS + beta*C*D
          C.matrixMultiply( D, RHS, beta, 1. );

          // Check norms
          RHS.matrixAdd( LHS, -1. );
          EXPECT_LT( RHS.normInf() / LHS.normInf(), machinePrecision );
          EXPECT_LT( RHS.norm1() / LHS.norm1(), machinePrecision );
          EXPECT_LT( RHS.normFrobenius() / LHS.normFrobenius(), machinePrecision );

        }

  // Test 4: repeat the following step for vectors of increasing size:
  //
  // a. compute LHS = ( alfa*A^T*B^T + beta*C ) * D
  // b. compute RHS = alfa*A^T*(B^T*D) + beta*C*D
  // c. check that (norm(LHS - RHS, type) / norm(LHS, type)) < epsilon
  //    with type = \{Infinity-norm, 1-norm, Frobenius-norm, )

  NA = 10;
  MB = 10;
  ND = 10;
  K = 20;

  for( localIndex nA = 1 ; nA <= NA ; ++nA )
    for( localIndex mB = 1 ; mB <= MB ; ++mB )
      for( localIndex nD = 1 ; nD <= ND ; ++nD )
        for( localIndex k = 1 ; k <= K ; ++k )
        {

          // Resize matrix operators
          A.resize( k, nA );
          B.resize( mB, k );
          C.resize( A.getNumCols(), B.getNumRows() );
          D.resize( B.getNumRows(), nD );
          LHS.resize( A.getNumCols(), D.getNumCols() );
          RHS.resize( A.getNumCols(), D.getNumCols() );

          // Populate A, B, C, and D with uniformly distributed random
          // coefficients
          A.rand();
          B.rand();
          C.rand();
          D.rand();

          // Compute tmp1 = ( alfa*A^T*B^T + beta*C )
          SerialMatrix tmp1( C );
          A.TmatrixTMultiply( B, tmp1, alfa, beta );

          // Compute LHS = tmp*D
          tmp1.matrixMultiply( D, LHS );

          // Compute tmp2 = B^T*D
          SerialMatrix tmp2( B.getNumCols(), D.getNumCols() );
          B.TmatrixMultiply( D, tmp2 );

          // Compute RHS = alfa*A^T*tmp2
          A.TmatrixMultiply( tmp2, RHS, alfa );

          // Compute RHS = RHS + beta*C*D
          C.matrixMultiply( D, RHS, beta, 1. );

          // Check norms
          RHS.matrixAdd( LHS, -1. );
          EXPECT_LT( RHS.normInf() / LHS.normInf(), machinePrecision );
          EXPECT_LT( RHS.norm1() / LHS.norm1(), machinePrecision );
          EXPECT_LT( RHS.normFrobenius() / LHS.normFrobenius(), machinePrecision );

        }

//  // Empty constructor
//
//  SerialMatrix A(2,2), B(1,1), C(4,4);
//  SerialMatrix Cinv;
//  array1d<localIndex> v_tmp;
//
//  std::cout << "default length: " << v_tmp.size() << std::endl;
//
//  // Assign values to A
//  A(0,0) = 3;
//  A(1,0) = 1;
//  A(0,1) = 1;
//  A(1,1) = 3;
//
//  C(0,0) = 3;
//  C(1,1) = 1;
//  C(2,2) = 0.5;
//  C(3,3) = -1;
//
//  std::cout << "C : \n";
//  C.print();
//  std::cout << "det(C) : " << C.determinant() << "\n";
//
//
//
//  C.computeInverse(Cinv);
//
//
//  SerialMatrix mat1(2,3);
//  SerialMatrix mat2(4,2);
//  mat1(0,0) = 1;
//  mat1(0,1) = 2;
//  mat1(0,2) = 3;
//  mat1(1,0) = 4;
//  mat1(1,1) = 5;
//  mat1(1,2) = 6;
//
//  mat2(0,0) = 1;
//  mat2(0,1) = 2;
//  mat2(1,0) = 3;
//  mat2(1,1) = 4;
//  mat2(2,0) = 5;
//  mat2(2,1) = 6;
//  mat2(3,0) = 7;
//  mat2(3,1) = 8;
//
//  std::cout << "normInf(mat2): " << mat2.normInf() << "\n";
//  std::cout << "  norm1(mat2): " << mat2.norm1() << "\n";
//  std::cout << "  normF(mat2): " << mat2.normFrobenius() << "\n";
//
//  SerialMatrix mat3(mat2.getNumRows(), mat1.getNumCols());
//
//  std::cout << "mat1 : \n";
//  mat1.print();
//  std::cout << "mat2 : \n";
//  mat2.print();
//
//  mat2.matrixMultiply(mat1, mat3);
//  std::cout << "mat2^T * mat1^T: " << std::endl;
//  mat3.print();
//
//  mat3.resize(mat1.getNumCols(),mat2.getNumRows());
//  mat1.TmatrixTMultiply(mat2, mat3);
//  std::cout << "mat2^T * mat1^T: " << std::endl;
//  mat3.print();
//
//  SerialMatrix mat4(4,3);
//  mat4(0,0) = 1;
//  mat4(0,1) = 2;
//  mat4(0,2) = 3;
//  mat4(1,0) = 4;
//  mat4(1,1) = 5;
//  mat4(1,2) = 6;
//  mat4(2,0) = 7;
//  mat4(2,1) = 8;
//  mat4(2,2) = 9;
//  mat4(3,0) = 10;
//  mat4(3,1) = 11;
//  mat4(3,2) = 12;
//
//  mat3.resize(mat1.getNumRows(),mat4.getNumRows());
//  mat1.matrixTMultiply(mat4, mat3);
//  std::cout << "mat2^T * mat1^T: " << std::endl;
//  mat3.print();
//
//
//  mat3.resize(mat4.getNumRows(),mat1.getNumRows());
//  mat4.matrixTMultiply(mat1, mat3);
//  std::cout << "mat2^T * mat1^T: " << std::endl;
//  mat3.print();
//
////  A.print();
////  std::cout << "Determinant of A: " << A.determinant() << std::endl;
////  std::cout << "Determinant of B: " << B.determinant() << std::endl;
////  std::cout << "Matrix C" << std::endl;
////  C.print();
////  std::cout << "Matrix C^-1" << std::endl;
////  Cinv.print();
////  BlasMatrix CinvxC(4,4);
////  CinvxC.GEMM(Cinv,C);
////  std::cout << "Matrix C^-1*C" << std::endl;
//////  CinvxC.print();
////
////  CinvxC.GEMM(C,Cinv);
////  std::cout << "Matrix C*C^-1" << std::endl;
////  CinvxC.print();
//
//  SerialMatrix D(C);
//  D.matrixAdd(Cinv);
//  std::cout << "Matrix D" << std::endl;
//  D.print();
//
//
//  SerialMatrix E;
//  SerialMatrix Einv;
//  SerialMatrix EinvXE;
//  real64 det;
//  localIndex max_dim = 8;
//  for (localIndex order = 1; order <= max_dim; ++order )
//  {
//    E.resize(order);
//    std::cout << "Matrix E" << std::endl;
//    E.print();
//    for (localIndex i = 0; i < E.getNumCols(); ++i)
//      for (localIndex j = 0; j < E.getNumRows(); ++j)
//      {
//        if ( i == j)
//        {
//          E(i,i) = 2;
//        }
//        else if (abs(i - j) == 1)
//        {
//          std::cout << "i: " << i << "j: " << j << "; abs(i - j): " << abs(i - j) << std::endl;
//          E(i,j) = -1;
//        }
//
//      }
//
//    std::cout << "Matrix E" << std::endl;
//    E.print();
//
//    E.computeInverse(Einv, det);
//    std::cout << "Matrix Einv" << std::endl;
//    Einv.print();
//    std::cout << "\n\ndet(E): " << det << "; det(Einv): " << Einv.determinant() << "; product: " << det*Einv.determinant() << "\n\n";
//
//    EinvXE.resize(E.getNumRows());
//    E.matrixMultiply(Einv, EinvXE);
//    std::cout << "Matrix Einv*E" << std::endl;
//    EinvXE.print();
//    std::cout << "Determinant(Einv*E): " << EinvXE.determinant() << std::endl;
//    std::cout << "****************************\n\n";
//  }
//
//  SerialMatrix matrixPermute(5,7);
//
//  for (localIndex j = 0; j < matrixPermute.getNumCols(); ++j)
//    for (localIndex i = 0; i < matrixPermute.getNumRows(); ++i)
//      matrixPermute(i,j) = static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
//
//  SerialMatrix matrixOrig(matrixPermute);
//
//  array1d<int> rowsPermutation,
//               colsPermutation; //start counting from 1
//  for(int i = 1; i <= matrixPermute.getNumRows(); ++i)
//  {
//    rowsPermutation.push_back(i);
//    std::cout << rowsPermutation[i-1] << " ";
//  }
//  std::cout << std::endl;
//  std::random_shuffle ( rowsPermutation.begin(), rowsPermutation.end() );
//  for(int i = 0; i < matrixPermute.getNumRows(); ++i)
//    std::cout << rowsPermutation[i] << " ";
//  std::cout << std::endl;
//
//  for(int i = 1; i <= matrixPermute.getNumCols(); ++i)
//  {
//    colsPermutation.push_back(i);
//    std::cout << colsPermutation[i-1] << " ";
//  }
//  std::cout << std::endl;
//  std::random_shuffle ( colsPermutation.begin(), colsPermutation.end() );
//  for(int i = 0; i < matrixPermute.getNumCols(); ++i)
//    std::cout << colsPermutation[i] << " ";
//  std::cout << std::endl;
//
//  std::cout << "\n\n\n\n\nMatrixPermute" << std::endl;
//  matrixPermute.print();
//  std::cout << "Row permuted matrixPermute" << std::endl;
//  matrixPermute.permuteRows(rowsPermutation);
//  matrixPermute.print();
//  std::cout << "Col permuted matrixPermute" << std::endl;
//  matrixPermute.permuteCols(colsPermutation);
//  matrixPermute.print();
//  std::cout << "Col permuted matrixPermute" << std::endl;
//  matrixPermute.permuteCols(colsPermutation, false);
//  matrixPermute.print();
//  std::cout << "Row permuted matrixPermute" << std::endl;
//  matrixPermute.permuteRows(rowsPermutation, false);
//  matrixPermute.print();
//
//  std::cout << "Row orig" << std::endl;
//  matrixOrig.print();
//  matrixOrig.matrixAdd(matrixPermute, -1.0);
//  std::cout << "Row orig" << std::endl;
//  matrixOrig.print();
//
//  SerialVector src_vec(matrixPermute.getNumCols()), dst_vec(matrixPermute.getNumRows());
//  for (localIndex i = 0; i < matrixPermute.getNumCols(); ++i)
//    src_vec(i) = static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
//
//
//
//  matrixPermute.vectorMultiply(src_vec, dst_vec);
//  std::cout << "dst_vec" << std::endl;
//  dst_vec.print();
//
//  matrixPermute.vectorMultiply(src_vec, dst_vec,1, -1);
//  std::cout << "dst_vec" << std::endl;
//  dst_vec.print();
//
//
////  SerialDenseMatrix A5(3,3);
////  SerialDenseMatrix B(3,2);
////  SerialDenseMatrix C(3,2);
////  SerialDenseMatrix D(2,2);
////  SerialDenseMatrix E(2,3);
////  SerialDenseMatrix F(1,1);S
//
//  EXPECT_EQ( A.getNumRows(), 2);
//  EXPECT_EQ( B.getNumRows(), 3);
}

template<typename LAI>
void testSerialMatrixVector()
{

  using SerialMatrix = typename LAI::SerialMatrix;
  using SerialVector = typename LAI::SerialVector;

  real64 machinePrecision = 10. * std::numeric_limits<real64>::epsilon();
  SerialMatrix A, yT, tmp;
  SerialVector x, y, yTT;

  // Test 5: repeat the following step for vectors of increasing size:
  //
  // a. compute y = A*x
  // b. compute compute yT = x^T * A^T
  // c. check that (norm(y, 2-norm) / norm(yT, Frobenius-norm)) < epsilon
  // d. create vector yTT = yT^T
  // e. compute beta = sqrt(yTT dot y);
  // f. check that ( (norm(y, norm-2) - sqrt(yTT_dot_y)) / norm(y, norm-2)) < epsilon

  localIndex MA = 10;
  localIndex NA = 10;
  real64 alfa, beta;

  for( localIndex mA = 1 ; mA <= MA ; ++mA )
    for( localIndex nA = 1 ; nA <= NA ; ++nA )
    {
      // Resize matrices and vectors
      A.resize( mA, nA );
      x.resize( nA);
      y.resize( mA);

      // Populate matrix A and vector x with uniformly distributed random
      // coefficients
      A.rand();
      x.rand();

      // a.
      A.vectorMultiply(x, y);

      // b.
      // --- construct tmp = x^T
      tmp.resize(1, x.getSize());
      for (localIndex i = 0; i < x.getSize(); ++i)
        tmp(0,i) = x(i);

      // -- compute yT = tmp * A^T
      yT.resize(1, A.getNumRows());
      tmp.matrixTMultiply(A, yT);

      // c.
      alfa = y.norm2();
      beta = std::fabs(alfa - yT.normFrobenius());
      EXPECT_LT( beta/alfa, machinePrecision );

      // d.
      yTT.resize(yT.getNumCols());
      for (localIndex i = 0; i < yT.getNumCols(); ++i)
        yTT(i) = yT(0,i);

      // e.
      beta = std::sqrt( y.dot(yTT) );

      // f. check that ( (norm(y, norm-2) - sqrt(yTT_dot_y)) / norm(y, norm-2)) < epsilon
      beta = alfa - beta;
      EXPECT_LT( beta/alfa, machinePrecision );

    }
}

template<typename LAI>
void testSerialMatrixInverse()
{

  using SerialMatrix = typename LAI::SerialMatrix;

  real64 machinePrecision = 10. * std::numeric_limits<real64>::epsilon();

  SerialMatrix E;
  SerialMatrix Einv;
  SerialMatrix EinvXE;

  // Test 5: repeat the following step for matrices of increasing size:
  // a. Construct matrix E (1d discrete Laplacian)
  // b. Compute Einv = E^-1
  // c. Compute EinvXE = Einv*E
  // d. Check that det(EinvXE) = 1.

  real64 det;
  localIndex max_dim = 10;

  for( localIndex order = 1 ; order <= max_dim ; ++order )
  {
    // a.
    E.resize( order );
    for( localIndex i = 0 ; i < E.getNumCols() ; ++i )
      for( localIndex j = 0 ; j < E.getNumRows() ; ++j )
      {
        if( i == j )
          E( i, i ) = 2;
        else if( abs( i - j ) == 1 )
          E( i, j ) = -1;
      }

    // b.
    E.computeInverse( Einv, det );

    // c.
    EinvXE.resize( E.getNumRows() );
    Einv.matrixMultiply( E, EinvXE );

    // d.
    EXPECT_LT( std::fabs( EinvXE.determinant() - 1. ), machinePrecision );
  }
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

  testSerialVector<LapackSuiteInterface>();
  testSerialMatrix<LapackSuiteInterface>();
  testSerialMatrixVector<LapackSuiteInterface>();
  testSerialMatrixInverse<LapackSuiteInterface>();

}

//@}

#ifdef __clang__
#pragma clang diagnostic pop
#endif
