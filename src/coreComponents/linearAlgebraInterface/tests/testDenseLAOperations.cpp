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
  using DenseMatrix = typename LAI::DenseMatrix;

  // Empty constructor

  DenseMatrix A(2,2), B(1,1), C(4,4);
  DenseMatrix Cinv;
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

  C.computeInverse(Cinv);

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

  BlasMatrix D(C);
  D.MatAdd(Cinv);
  std::cout << "Matrix D" << std::endl;
  D.print();


  BlasMatrix E;
  localIndex max_dim = 6;
  for (localIndex order = 1; order <= max_dim; ++order )
  {
    E.resize(order);
    for (localIndex i = 0; i < E.getNumCols(); ++i)
      for (localIndex j = 0; j < E.getNumRows(); ++j)
      {
        if ( i == j)
        {
          E(i,i) = 2;
        }
        else if (abs(i - j) == 1)
        {
          E(i,j) = -1;
        }

      }

    std::cout << "Matrix E" << std::endl;
    E.print();

    BlasMatrix Einv;
    E.computeInverse(Einv);
    std::cout << "Matrix Einv" << std::endl;
    Einv.print();

    BlasMatrix EinvXE(E.getNumRows());
    EinvXE.GEMM(Einv, E);
    std::cout << "Matrix Einv*E" << std::endl;
    EinvXE.print();
  }
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
