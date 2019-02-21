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
  using SerialMatrix = typename LAI::SerialMatrix;
  using SerialVector = typename LAI::SerialVector;

  // Empty constructor

  SerialMatrix A(2,2), B(1,1), C(4,4);
  SerialMatrix Cinv;
  array1d<localIndex> v_tmp;

  std::cout << "default length: " << v_tmp.size() << std::endl;

  // Assign values to A
  A(0,0) = 3;
  A(1,0) = 1;
  A(0,1) = 1;
  A(1,1) = 3;

  C(0,0) = 3;
  C(1,1) = 1;
  C(2,2) = 0.5;
  C(3,3) = -1;

  std::cout << "C : \n";
  C.print();
  std::cout << "det(C) : " << C.determinant() << "\n";



  C.computeInverse(Cinv);


  SerialMatrix mat1(2,3);
  SerialMatrix mat2(4,2);
  mat1(0,0) = 1;
  mat1(0,1) = 2;
  mat1(0,2) = 3;
  mat1(1,0) = 4;
  mat1(1,1) = 5;
  mat1(1,2) = 6;

  mat2(0,0) = 1;
  mat2(0,1) = 2;
  mat2(1,0) = 3;
  mat2(1,1) = 4;
  mat2(2,0) = 5;
  mat2(2,1) = 6;
  mat2(3,0) = 7;
  mat2(3,1) = 8;

  std::cout << "normInf(mat2): " << mat2.normInf() << "\n";
  std::cout << "  norm1(mat2): " << mat2.norm1() << "\n";
  std::cout << "  normF(mat2): " << mat2.normFrobenius() << "\n";

  SerialMatrix mat3(mat2.getNumRows(), mat1.getNumCols());

  std::cout << "mat1 : \n";
  mat1.print();
  std::cout << "mat2 : \n";
  mat2.print();

  mat2.matrixMultiply(mat1, mat3);
  std::cout << "mat2^T * mat1^T: " << std::endl;
  mat3.print();

  mat3.resize(mat1.getNumCols(),mat2.getNumRows());
  mat1.TmatrixTMultiply(mat2, mat3);
  std::cout << "mat2^T * mat1^T: " << std::endl;
  mat3.print();

  SerialMatrix mat4(4,3);
  mat4(0,0) = 1;
  mat4(0,1) = 2;
  mat4(0,2) = 3;
  mat4(1,0) = 4;
  mat4(1,1) = 5;
  mat4(1,2) = 6;
  mat4(2,0) = 7;
  mat4(2,1) = 8;
  mat4(2,2) = 9;
  mat4(3,0) = 10;
  mat4(3,1) = 11;
  mat4(3,2) = 12;

  mat3.resize(mat1.getNumRows(),mat4.getNumRows());
  mat1.matrixTMultiply(mat4, mat3);
  std::cout << "mat2^T * mat1^T: " << std::endl;
  mat3.print();


  mat3.resize(mat4.getNumRows(),mat1.getNumRows());
  mat4.matrixTMultiply(mat1, mat3);
  std::cout << "mat2^T * mat1^T: " << std::endl;
  mat3.print();

//  A.print();
//  std::cout << "Determinant of A: " << A.determinant() << std::endl;
//  std::cout << "Determinant of B: " << B.determinant() << std::endl;
//  std::cout << "Matrix C" << std::endl;
//  C.print();
//  std::cout << "Matrix C^-1" << std::endl;
//  Cinv.print();
//  BlasMatrix CinvxC(4,4);
//  CinvxC.GEMM(Cinv,C);
//  std::cout << "Matrix C^-1*C" << std::endl;
////  CinvxC.print();
//
//  CinvxC.GEMM(C,Cinv);
//  std::cout << "Matrix C*C^-1" << std::endl;
//  CinvxC.print();

  SerialMatrix D(C);
  D.matrixAdd(Cinv);
  std::cout << "Matrix D" << std::endl;
  D.print();


  SerialMatrix E;
  SerialMatrix Einv;
  SerialMatrix EinvXE;
  real64 det;
  localIndex max_dim = 8;
  for (localIndex order = 1; order <= max_dim; ++order )
  {
    E.resize(order);
    std::cout << "Matrix E" << std::endl;
    E.print();
    for (localIndex i = 0; i < E.getNumCols(); ++i)
      for (localIndex j = 0; j < E.getNumRows(); ++j)
      {
        if ( i == j)
        {
          E(i,i) = 2;
        }
        else if (abs(i - j) == 1)
        {
          std::cout << "i: " << i << "j: " << j << "; abs(i - j): " << abs(i - j) << std::endl;
          E(i,j) = -1;
        }

      }

    std::cout << "Matrix E" << std::endl;
    E.print();

    E.computeInverse(Einv, det);
    std::cout << "Matrix Einv" << std::endl;
    Einv.print();
    std::cout << "\n\ndet(E): " << det << "; det(Einv): " << Einv.determinant() << "; product: " << det*Einv.determinant() << "\n\n";

    EinvXE.resize(E.getNumRows());
    E.matrixMultiply(Einv, EinvXE);
    std::cout << "Matrix Einv*E" << std::endl;
    EinvXE.print();
    std::cout << "Determinant(Einv*E): " << EinvXE.determinant() << std::endl;
    std::cout << "****************************\n\n";
  }

  SerialMatrix matrixPermute(5,7);

  for (localIndex j = 0; j < matrixPermute.getNumCols(); ++j)
    for (localIndex i = 0; i < matrixPermute.getNumRows(); ++i)
      matrixPermute(i,j) = static_cast<double>(rand()) / static_cast<double>(RAND_MAX);

  SerialMatrix matrixOrig(matrixPermute);

  array1d<int> rowsPermutation,
               colsPermutation; //start counting from 1
  for(int i = 1; i <= matrixPermute.getNumRows(); ++i)
  {
    rowsPermutation.push_back(i);
    std::cout << rowsPermutation[i-1] << " ";
  }
  std::cout << std::endl;
  std::random_shuffle ( rowsPermutation.begin(), rowsPermutation.end() );
  for(int i = 0; i < matrixPermute.getNumRows(); ++i)
    std::cout << rowsPermutation[i] << " ";
  std::cout << std::endl;

  for(int i = 1; i <= matrixPermute.getNumCols(); ++i)
  {
    colsPermutation.push_back(i);
    std::cout << colsPermutation[i-1] << " ";
  }
  std::cout << std::endl;
  std::random_shuffle ( colsPermutation.begin(), colsPermutation.end() );
  for(int i = 0; i < matrixPermute.getNumCols(); ++i)
    std::cout << colsPermutation[i] << " ";
  std::cout << std::endl;

  std::cout << "\n\n\n\n\nMatrixPermute" << std::endl;
  matrixPermute.print();
  std::cout << "Row permuted matrixPermute" << std::endl;
  matrixPermute.permuteRows(rowsPermutation);
  matrixPermute.print();
  std::cout << "Col permuted matrixPermute" << std::endl;
  matrixPermute.permuteCols(colsPermutation);
  matrixPermute.print();
  std::cout << "Col permuted matrixPermute" << std::endl;
  matrixPermute.permuteCols(colsPermutation, false);
  matrixPermute.print();
  std::cout << "Row permuted matrixPermute" << std::endl;
  matrixPermute.permuteRows(rowsPermutation, false);
  matrixPermute.print();

  std::cout << "Row orig" << std::endl;
  matrixOrig.print();
  matrixOrig.matrixAdd(matrixPermute, -1.0);
  std::cout << "Row orig" << std::endl;
  matrixOrig.print();

  Vector src_vec(matrixPermute.getNumCols()), dst_vec(matrixPermute.getNumRows());
  for (localIndex i = 0; i < matrixPermute.getNumCols(); ++i)
    src_vec(i) = static_cast<double>(rand()) / static_cast<double>(RAND_MAX);



  matrixPermute.vectorMultiply(src_vec, dst_vec);
  std::cout << "dst_vec" << std::endl;
  dst_vec.print();

  matrixPermute.vectorMultiply(src_vec, dst_vec,1, -1);
  std::cout << "dst_vec" << std::endl;
  dst_vec.print();


//  SerialDenseMatrix A5(3,3);
//  SerialDenseMatrix B(3,2);
//  SerialDenseMatrix C(3,2);
//  SerialDenseMatrix D(2,2);
//  SerialDenseMatrix E(2,3);
//  SerialDenseMatrix F(1,1);S

  EXPECT_EQ( A.getNumRows(), 2);
  EXPECT_EQ( B.getNumRows(), 3);
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
  testMatrixConstructors<LapackSuiteInterface>();
}

//@}

#ifdef __clang__
#pragma clang diagnostic pop
#endif
