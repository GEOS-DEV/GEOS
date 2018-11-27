/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
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
 * @file testLAOperations.cpp
 *
 *  Created on: Sep 19, 2018
 *      Author: Matthias Cremon
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

#include <iostream>
#include <vector>
#include <mpi.h>

#include "common/DataTypes.hpp"

#include "TrilinosInterface.hpp"
//#include "HypreInterface.hpp"
#include "NativeCG.hpp"
#include "NativeBICGSTAB.hpp"
#include "BlockMatrixView.hpp"


/**
 * \file testLAOperations.cpp
 * \brief This test file is part of the ctest suite and tests the Trilinos based solvers
 * as well as block GEOSX solvers. It mainly uses dummy 2D Laplace operator matrices
 * to test the solvers along with a simple (block) identity preconditioner.
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

// BEGIN_RST_NARRATIVE testLAOperations.rst
// ==============================
// Compute Identity
// ==============================
// This function computes the identity matrix. It can be used to generate a dummy
// preconditioner.

template< typename LAI >
typename LAI::ParallelMatrix computeIdentity( MPI_Comm comm,
                                              typename LAI::gid N )
{
  // Declare matrix
  typename LAI::ParallelMatrix I;

  // Create a matrix of size N with 1 non-zero per row
  I.createWithGlobalSize(N,1,comm);

  // Loop over rows to fill the matrix
  for ( typename LAI::gid i = I.ilower(); i < I.iupper(); i++ )
  {
    // Set the value for element (i,i) to 1
    I.insert( i, i, 1.0);
  }

  // Close the matrix (make data contiguous in memory)
  I.close();

  // Return the matrix.
  return I;
}


/**
 * @brief Compute the 2D Laplace operator
 *
 * \param comm MPI communicator.
 * \param global size N of the square 2D Laplace operator matrix.
 */

// ==============================
// Compute 2D Laplace Operator
// ==============================

// This function computes the matrix corresponding to a 2D Laplace operator. These
// matrices arise from a classical finite volume formulation on a cartesian mesh
// (5-point stencil).

template< typename LAI >
typename LAI::ParallelMatrix compute2DLaplaceOperator( MPI_Comm comm,
                                                       typename LAI::gid N )
{
  // Declare matrix
  typename LAI::ParallelMatrix laplace2D;

  // Create a matrix of global size N with 5 non-zeros per row
  laplace2D.createWithGlobalSize(N,5,comm);

  // Get the size of the dummy mesh back to be able to put values in the correct
  // diagonals.
  typename LAI::gid n = std::sqrt( N );

  // Allocate arrays to fill the matrix (values and columns)
  real64 values[5];
  typename LAI::gid cols[5];

  // Loop over rows to fill the matrix
  for ( typename LAI::gid i = laplace2D.ilower(); i < laplace2D.iupper(); i++ )
  {
    // Re-set the number of non-zeros for row i to 0.
    typename LAI::lid nnz = 0;

    // The left -n: position i-n
    if ( i-n >= 0 )
    {
      cols[nnz] = i-n;
      values[nnz] = -1.0;
      nnz++;
    }

    // The left -1: position i-1
    if ( i-1 >= 0 )
    {
      cols[nnz] = i-1;
      values[nnz] = -1.0;
      nnz++;
    }

    // Set the diagonal: position i
    cols[nnz] = i;
    values[nnz] = 4.0;
    nnz++;

    // The right +1: position i+1
    if ( i+1 < N )
    {
      cols[nnz] = i+1;
      values[nnz] = -1.0;
      nnz++;
    }

    // The right +n: position i+n
    if ( i+n < N )
    {
      cols[nnz] = i+n;
      values[nnz] = -1.0;
      nnz++;
    }

    // Set the values for row i
    laplace2D.insert( i, cols, values, nnz );
  }

  // Close the matrix (make data contiguous in memory)
  laplace2D.close();

  // Return the matrix
  return laplace2D;
}

//@}


/*! @name Test functions.
 * @brief Templated functions to test the linear solvers.
 */
//@{
/**
 * @function testInterfaceSolvers
 *
 * @brief Test the packaged solvers from the LAI as well as basic linear algebra operations,
 * such as matrix-vector products, dot products, norms and residuals.
 */

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

template< typename LAI >
void testInterfaceSolvers()
{
  // Define aliases templated on the Linear Algebra Interface (LAI).
  using ParallelMatrix = typename LAI::ParallelMatrix;
  using ParallelVector = typename LAI::ParallelVector;
  using LinearSolver = typename LAI::LinearSolver;
  using LAIgid = typename LAI::gid;

  // Get the MPI rank
  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  // Use an nxn cartesian mesh to generate the Laplace 2D operator.
  LAIgid n = 100; 
  LAIgid N = n*n;

  // Compute a 2D Laplace operator
  ParallelMatrix matrix = compute2DLaplaceOperator<LAI>( MPI_COMM_WORLD, N );

  // Define some vectors
  ParallelVector x_true, 
                 x_comp, 
                 b;

  x_true.createWithGlobalSize( N, MPI_COMM_WORLD );
  x_comp.createWithGlobalSize( N, MPI_COMM_WORLD );
  b.createWithGlobalSize( N, MPI_COMM_WORLD );

  // We have some simple initialization options for vectors:
  x_true.rand(); // random
  x_comp.zero(); // zero 
  b.set(1.0);    // ones

  // Also define a residual vector, this time using the copy constructor
  ParallelVector r( b );

  // Test dot product: r.b = b.b = N
  real64 dotTest = r.dot(b);
  EXPECT_DOUBLE_EQ( dotTest, N );

  // Test various norms
  real64 norm1 = b.norm1();
  real64 norm2 = b.norm2();
  real64 normInf = b.normInf();

  EXPECT_DOUBLE_EQ( norm1, N );
  EXPECT_DOUBLE_EQ( norm2, n );
  EXPECT_DOUBLE_EQ( normInf, 1. );

  // Compute the matrix/vector multiplication. We compute b as Ax and will aim to get x
  // back from the solvers.
  matrix.multiply( x_true, b );

  // Test the residual function by computing r = b - Ax = 0
  matrix.residual( x_true, b, r );
  real64 normRes = r.normInf();
  EXPECT_DOUBLE_EQ( normRes, 0. );

  // Now declare a solver object with desired options
  //LinearSolverParameters parameters;
  //                       parameters.krylov_tolerance = 1e-8;
  //                       parameters.max_krylov_iterations = 1000;
  LinearSolver solver = LinearSolver();
  //solver.setParameters(parameters);

  // Solve using the iterative solver and compare norms with true solution
  solver.solve( matrix, x_comp, b, 1000, 1e-8);
  real64 norm_comp = x_comp.norm2();
  real64 norm_true = x_true.norm2();
  EXPECT_LT( std::fabs( norm_comp/norm_true - 1. ), 1e-6 );

  // We also check a couple random elements to see if they are equal.
  //JAW  EXPECT_LT( ( std::fabs( solIterative.getElement( n ) - x.getElement( n ))/std::fabs( x.getElement( n ) ) ), 1e-4 );
  //JAW  EXPECT_LT( ( std::fabs( solIterative.getElement( 3*n ) - x.getElement( 3*n ))/std::fabs( x.getElement( 3*n ) ) ), 1e-4 );

  // We now do the same using a direct solver.
  // Again the norm should be the norm of x. We use a tougher tolerance on the test
  // compared to the iterative solution. This should be accurate to machine precision
  // and some round off error. We (arbitrarily) chose 1e-12 as a good guess.
  //parameters.type = "direct";
  //solver.setParameters(parameters);
  x_comp.zero();
  solver.dsolve( matrix, x_comp, b );
  norm_comp = x_comp.norm2();
  EXPECT_LT( std::fabs( norm_comp/norm_true - 1. ), 1e-12 );

  // Again we check a couple elements
  //JAW  EXPECT_LT( ( std::fabs( solDirect.getElement( n ) - x.getElement( n ))/std::fabs( x.getElement( n ) ) ), 1e-12 );
  //JAW  EXPECT_LT( ( std::fabs( solDirect.getElement( 3*n ) - x.getElement( 3*n ))/std::fabs( x.getElement( 3*n ) ) ), 1e-12 );

#if 0
  // Declare integers to run tests
  integer numValRow0,numValRow1,numValRown;

  // Declare vectors to run tests
  std::vector<real64> vecValuesRow0( 5 ),vecValuesRow1( 5 ),vecValuesRown( 5 );
  std::vector<LAI::lid> vecIndicesRow0( 5 ),vecIndicesRow1( 5 ),vecIndicesRown( 5 );

  // Get values and columns in specific rows
  matrix.getLocalRow( 0, numValRow0, vecValuesRow0, vecIndicesRow0 );
  matrix.getLocalRow( 1, numValRow1, vecValuesRow1, vecIndicesRow1 );
  matrix.getLocalRow( static_cast<LAI::lid>( n ), numValRown, vecValuesRown, vecIndicesRown );

  // Run checks on rank 0 to see if the matrix was correctly constructed
  if (rank == 0)
  {
    // Check number of values per row
    EXPECT_DOUBLE_EQ( numValRow0, 3 );
    EXPECT_DOUBLE_EQ( numValRow1, 4 );
    EXPECT_DOUBLE_EQ( numValRown, 5 );

    // Check column indices for row 0
    EXPECT_DOUBLE_EQ( vecIndicesRow0[0], 0 );
    EXPECT_DOUBLE_EQ( vecIndicesRow0[1], 1 );
    EXPECT_DOUBLE_EQ( vecIndicesRow0[2], n );

    // Check column indices for row 1
    EXPECT_DOUBLE_EQ( vecIndicesRow1[0], 0 );
    EXPECT_DOUBLE_EQ( vecIndicesRow1[1], 1 );
    EXPECT_DOUBLE_EQ( vecIndicesRow1[2], 2 );
    EXPECT_DOUBLE_EQ( vecIndicesRow1[3], n+1 );

    // Check column indices for row n
    EXPECT_DOUBLE_EQ( vecIndicesRown[0], 0 );
    EXPECT_DOUBLE_EQ( vecIndicesRown[1], n-1 );
    EXPECT_DOUBLE_EQ( vecIndicesRown[2], n );
    EXPECT_DOUBLE_EQ( vecIndicesRown[3], n+1 );
    EXPECT_DOUBLE_EQ( vecIndicesRown[4], 2*n );

    // Check values for row 0
    EXPECT_DOUBLE_EQ( vecValuesRow0[0], 4. );
    EXPECT_DOUBLE_EQ( vecValuesRow0[1], -1. );
    EXPECT_DOUBLE_EQ( vecValuesRow0[2], -1. );

    // Check values for row 1
    EXPECT_DOUBLE_EQ( vecValuesRow1[0], -1. );
    EXPECT_DOUBLE_EQ( vecValuesRow1[1], 4. );
    EXPECT_DOUBLE_EQ( vecValuesRow1[2], -1. );
    EXPECT_DOUBLE_EQ( vecValuesRow1[3], -1. );

    // Check values for row n
    EXPECT_DOUBLE_EQ( vecValuesRown[0], -1. );
    EXPECT_DOUBLE_EQ( vecValuesRown[1], -1. );
    EXPECT_DOUBLE_EQ( vecValuesRown[2], 4. );
    EXPECT_DOUBLE_EQ( vecValuesRown[3], -1. );
    EXPECT_DOUBLE_EQ( vecValuesRown[4], -1. );
  }

  // We now test the linear algebra operations, such as matrix-vector product and
  // compute residuals.

  // For that we first fill some standard vectors
  std::vector<real64> ones, zer, random, random2;
  for (integer j = 0; j < N; j++)
  {
    // Vector of zeros
    zer.push_back( 0 );
    // Vector of ones
    ones.push_back( 1 );
    // Vector of random integers
    random.push_back( rand() % 20 - 10 );
    // Vector of random integers (different range to make sure the two random vectors
    // are not the same, which would defeat the purpose).
    random2.push_back( rand() % 30 - 15 );
  }

  // Test the clearRow function (this has to be done after the solves so that we can
  // use the same matrix when we are done with it).
  // We clear the row and multiply the diagonal value by 2.
  matrix.clearRow( 2*N/4+n, 2.0 );
  matrix.getLocalRow(static_cast<LAI::lid>( n ),numValRown,vecValuesRown,vecIndicesRown);
  if ( rank == 2 )
  {
    EXPECT_DOUBLE_EQ( vecValuesRown[2], 8.0 );
  }
#endif

}

/**
 * @function testGEOSXSolvers
 *
 * @brief Test the GEOSX solvers for monolithic matrices by solving a system with
 * a Laplace operator and the identity matrix as a (dummy) preconditioner.
 */

// -----------------------------------------
// Test GEOSX solvers on monolithic matrices
// -----------------------------------------
// We now test the GEOSX implementation of the Conjugate Gradient (CG) and BiCGSTAB algorithms
// on monolithic matrices.

template< typename LAI >
void testGEOSXSolvers()
{
  // Define aliases templated on the Linear Algebra Interface (LAI).
  using ParallelMatrix = typename LAI::ParallelMatrix;
  using ParallelVector = typename LAI::ParallelVector;
  using LAIgid = typename LAI::gid;

  // Get the MPI rank
  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  // Use nxn cartesian mesh to generate the Laplace 2D operator.
  LAIgid n = 100;
  LAIgid N = n*n;

  // Compute a 2D Laplace operator and identity matrix
  ParallelMatrix matrix = compute2DLaplaceOperator<LAI>( MPI_COMM_WORLD, N );
  ParallelMatrix identity = computeIdentity<LAI>( MPI_COMM_WORLD, N );

  // Define vectors 
  ParallelVector x_true,
                 x_comp,
                 b;

  x_true.createWithGlobalSize( N, MPI_COMM_WORLD);
  x_comp.createWithGlobalSize( N, MPI_COMM_WORLD);
  b.createWithGlobalSize( N, MPI_COMM_WORLD);

  x_true.rand();
  x_comp.zero();

  matrix.multiply( x_true, b );

  // Test with CG solver
  CGsolver<LAI> testCG;
                testCG.solve( matrix, x_comp, b, identity );

  real64 norm_true = x_true.norm2();
  real64 norm_comp = x_comp.norm2();
  EXPECT_LT( std::fabs( norm_comp/norm_true - 1. ), 5e-6 );

  // Test with BiCGSTAB solver 
  x_comp.zero();

  BiCGSTABsolver<LAI> testBiCGSTAB;
  testBiCGSTAB.solve( matrix, x_comp, b, identity );

  norm_true = x_true.norm2();
  norm_comp = x_comp.norm2();
  EXPECT_LT( std::fabs( norm_comp/norm_true - 1. ), 5e-6 );
}


/**
 * @function testGEOSXBlockSolvers
 *
 * @brief Test the GEOSX block solvers by solving a system with a block matrix
 * made of tiled Laplace operators and using the identity as a preconditioner.
 */

// -----------------------------------------
// Test GEOSX solvers on block matrices
// -----------------------------------------
// We finish by testing the GEOSX implementation of the Conjugate Gradient (CG) and BiCGSTAB algorithms
// on block matrices.

template< typename LAI >
void testGEOSXBlockSolvers()
{
  // The usual typenames
  using ParallelMatrix = typename LAI::ParallelMatrix;
  using ParallelVector = typename LAI::ParallelVector;
  using LAIgid = typename LAI::gid;

  // Get the MPI rank
  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  // We are going to assembly the following dummy system
  // [L L] [x_true] = [b_0]
  // [L L] [x_true] = [b_1]

  LAIgid n = 100;
  LAIgid N = n*n;

  ParallelMatrix matrix   = compute2DLaplaceOperator<LAI>( MPI_COMM_WORLD, N );
  ParallelMatrix identity = computeIdentity<LAI>( MPI_COMM_WORLD, N );

  ParallelVector x_true,
                 x_comp_0,
                 x_comp_1,
                 b_0,
                 b_1;

  x_true.createWithGlobalSize( N, MPI_COMM_WORLD);
  x_comp_0.createWithGlobalSize( N, MPI_COMM_WORLD);
  x_comp_1.createWithGlobalSize( N, MPI_COMM_WORLD);
  b_0.createWithGlobalSize( N, MPI_COMM_WORLD);
  b_1.createWithGlobalSize( N, MPI_COMM_WORLD);

  x_true.rand();
  x_comp_0.zero();
  x_comp_1.zero();

  // Size of the block system
  integer nRows = 2;
  integer nCols = 2;

  // Declare and allocate block matrices/vectors
  // System matrix

  BlockMatrixView<LAI> block_matrix( nRows, nCols );
  BlockMatrixView<LAI> block_precon( nRows, nCols );
  BlockVectorView<LAI> block_x_true( nCols );
  BlockVectorView<LAI> block_x_comp( nCols );
  BlockVectorView<LAI> block_rhs   ( nRows );

  // In this test we simply tile the laplace operator, so we assign a duplicate
  // of the monolithic matrix to every block of the matrix.

  block_matrix.set( 0, 0, matrix );
  block_matrix.set( 0, 1, matrix );
  block_matrix.set( 1, 0, matrix );
  block_matrix.set( 1, 1, matrix );

  // We do the same for the preconditioner to get the block identity matrix.
  // We ignore the off-diagonal blocks and leave them as null-pointers.

  block_precon.set( 0, 0, identity );
  block_precon.set( 1, 1, identity );

  // true solution
  block_x_true.set( 0, x_true );
  block_x_true.set( 1, x_true );

  // Set initial guess blocks (here we need multiple objects since we cannot
  // have them point to the same memory location). These objects are initial
  // guesses as input but solution vectors as output.

  block_x_comp.set( 0, x_comp_0 );
  block_x_comp.set( 1, x_comp_1 );

  // Set right hand side blocks.
  block_rhs.set( 0, b_0 );
  block_rhs.set( 1, b_1 );
  block_matrix.multiply( block_x_true, block_rhs);

//TODO: Need to refactor Native block solvers.  Disable this for now.
#if 0
  // Create block CG solver object and solve
  CGsolver<LAI> testCG;
  testCG.solve( block_matrix, block_x_comp, block_rhs, block_precon );

  // The true solution is the vector x, so we check if the norm are equal.
  // Note: the tolerance is higher that in the previous cases because the matrix is
  // twice as big and the condition number is higher. See details on the error
  // bounds of Krylov methods wrt the condition number.
  real64 norm_x_true = block_x_true.norm2();
  real64 norm_x_comp = block_x_comp.norm2();
  EXPECT_LT( std::fabs( norm_x_comp/norm_x_true - 1. ), 5e-6 );

  // now try out the BiCGstab solver

  x_comp_0.zero(); // TODO: fix block zero()
  x_comp_1.zero();

  BiCGSTABsolver<LAI> testBiCGSTAB;
  testBiCGSTAB.solve( block_matrix, block_x_comp, block_rhs, block_precon );

  norm_x_true = block_x_true.norm2();
  norm_x_comp = block_x_comp.norm2();
  EXPECT_LT( std::fabs( norm_x_comp/norm_x_true - 1. ), 5e-6 );
#endif
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



/*! @function testEpetraLAOperations.
 * @brief Runs all tests using the Trilinos interface.
 */

TEST(testLAOperations,testEpetraLAOperations)
{
  MPI_Init( nullptr, nullptr );
  testInterfaceSolvers<TrilinosInterface>();
  testGEOSXSolvers<TrilinosInterface>();
  testGEOSXBlockSolvers<TrilinosInterface>();
  MPI_Finalize();
}

/*! @function testHypreLAOperations.
 * @brief Runs all tests using the Hypre interface.
 */
TEST(testLAOperations,testHypreLAOperations)
{
  //MPI_Init( nullptr, nullptr );
  //testInterfaceSolvers<HypreInterface>();
  //testGEOSXSolvers<HypreInterface>();
  //testGEOSXBlockSolvers<HypreInterface>();
  //MPI_Finalize();
}

/*! @function testPETScLAOperations.
 * @brief Runs all tests using the PETSc interface.
 */
TEST(testLAOperations,testPETScLAOperations)
{

  //MPI_Init( nullptr, nullptr );
  //testInterfaceSolvers<PETScInterface>();
  //testGEOSXSolvers<PETScInterface>();
  //testGEOSXBlockSolvers<PETScInterface>();
  //MPI_Finalize();
}

//@}

#ifdef __clang__
#pragma clang diagnostic pop
#endif
