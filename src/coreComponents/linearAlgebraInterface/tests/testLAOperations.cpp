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

#include "TrilinosInterface.hpp"
//#include "HypreInterface.hpp"
#include "BlockMatrixView.hpp"
#include "../src/NativeCG.hpp"
#include "../src/NativeBICGSTAB.hpp"

#include "common/DataTypes.hpp"

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
  integer n = integer( std::sqrt( N ) );

  // Allocate arrays to fill the matrix (values and columns)

  real64 values[5];
  typename LAI::gid cols[5];

  // Loop over rows to fill the matrix
  for ( typename LAI::gid i = laplace2D.ilower(); i < laplace2D.iupper(); i++ )
  {
    // Re-set the number of non-zeros for row i to 0.
    integer nnz = 0;

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
 * @function testNativeSolvers
 *
 * @brief Test the native solvers from the LAI as well as basic linear algebra operations,
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
void testNativeSolvers()
{
  // Define aliases templated on the Linear Algebra Interface (LAI).
  // These objects can use all of the available libraries (trilinos and
  // soon HYPRE and PETSc).
  using ParallelMatrix = typename LAI::ParallelMatrix;
  using ParallelVector = typename LAI::ParallelVector;
  using LinearSolver = typename LAI::LinearSolver;
  using LAIgid = typename LAI::gid;
 
  // Initialize MPI
  MPI_Init( nullptr, nullptr );

  // Get the MPI rank
  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  // Set the MPI communicator
  // MPI_Comm comm = MPI_COMM_WORLD;

  MPI_Comm test_comm;
  MPI_Comm_dup( MPI_COMM_WORLD, &test_comm );
  MPI_Comm comm = test_comm;

  // Size of the dummy cartesian mesh we use to generate the Laplace 2D operator.
  LAIgid n = 100;
  // Number of degrees of freedom
  LAIgid N = n*n;

  // Compute a 2D Laplace operator
  ParallelMatrix testMatrix = compute2DLaplaceOperator<LAI>( comm, N );

  // Compute a dummy preconditioner (identity matrix)
  ParallelMatrix preconditioner = computeIdentity<LAI>( comm, N );


#if 0
  // Declare integers to run tests
  integer numValRow0,numValRow1,numValRown;

  // Declare vectors to run tests
  std::vector<real64> vecValuesRow0( 5 ),vecValuesRow1( 5 ),vecValuesRown( 5 );
  std::vector<LAI::lid> vecIndicesRow0( 5 ),vecIndicesRow1( 5 ),vecIndicesRown( 5 );

  // Get values and columns in specific rows
  testMatrix.getLocalRow( 0, numValRow0, vecValuesRow0, vecIndicesRow0 );
  testMatrix.getLocalRow( 1, numValRow1, vecValuesRow1, vecIndicesRow1 );
  testMatrix.getLocalRow( static_cast<LAI::lid>( n ), numValRown, vecValuesRown, vecIndicesRown );

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
#endif

  // Define x, b and x0 guess vectors
  ParallelVector x, b, x0;

  // Place holder for b (will be recomputed a few lines down). We do use the vector to test
  // the linear algebra operations
  b.createWithGlobalSize( N, comm);
  b.set(1.0);

  // Random vector for x (the true solution of Ax = b).
  x.createWithGlobalSize( N, comm );
  x.rand();

  // Get the inf-norm of the true solution (x).
  real64 norm2x = x.norm2();

  // Initial guess x0 = 0
  x0.createWithGlobalSize( N, comm );
  x0.zero();

  // Fill initial guess for the direct solver
  ParallelVector solDirect( x0 );

  // Fill initial guess for the iterative solver
  ParallelVector solIterative( x0 );

  // Residual vector (vector of zeros)
  ParallelVector r( b );

  // Test dot product
  real64 dotTest = b.dot(b);
  EXPECT_DOUBLE_EQ( dotTest, N );

  // Test 1-norm
  real64 norm1 = b.norm1();
  // The 1-norm should be equal to N
  EXPECT_DOUBLE_EQ( norm1, N );

  // Test 2-norm
  real64 norm2 = b.norm2();
  // The 2-norm should be equal to n
  EXPECT_DOUBLE_EQ( norm2, n );

  // Test inf-norm
  real64 norminf = b.normInf();
  // The inf-norm should be equal to 1
  EXPECT_DOUBLE_EQ( norminf, 1. );

  // Compute the matrix/vector multiplication. We compute b as Ax and will aim to get x
  // back from the solvers.
  testMatrix.multiply( x, b );

  // First we test the residual function by computing r = b - Ax.
  testMatrix.residual( x, b, r );

  // Compute the norm of the residual
  real64 normRes = r.normInf();
  // The inf-norm should be equal to 0 (all elements are 0 in exact algebra).
  EXPECT_DOUBLE_EQ( normRes, 0. );

  // We now test the linear solvers from the libraries.

  // Declare solver object.
  LinearSolver solver = LinearSolver();

  // Run the iterative solver (TODO remove hard coded options, currently
  // using ILUT preconditioned GMRES).
  if (rank == 0)
    std::cout << std::endl << "AztecOO iterative solver:";
  solver.solve( testMatrix, solIterative, b, 1000, 1e-8);

  // Get the inf-norm of the iterative solution.
  real64 normIterativeSol = solIterative.norm2();
  // The true solution is the vector x, so we check if the norm are equal.
  EXPECT_LT( std::fabs( normIterativeSol/norm2x - 1. ), 1e-6 );

  // We also check a couple random elements to see if they are equal.
//JAW  EXPECT_LT( ( std::fabs( solIterative.getElement( n ) - x.getElement( n ))/std::fabs( x.getElement( n ) ) ), 1e-4 );
//JAW  EXPECT_LT( ( std::fabs( solIterative.getElement( 3*n ) - x.getElement( 3*n ))/std::fabs( x.getElement( 3*n ) ) ), 1e-4 );

  // We now do the same using a direct solver from Amesos (Klu)
  if (rank == 0)
    std::cout << std::endl << "Amesos direct solver:" << std::endl << std::endl;
  solver.dsolve( testMatrix, solDirect, b );

  // Again the norm should be the norm of x. We use a tougher tolerance on the test
  // compared to the iterative solution. This should be accurate to machine precision
  // and some round off error. We (arbitrarily) chose 1e-12 as a good guess.
  real64 normDirectSol = solDirect.norm2();
  EXPECT_LT( std::fabs( normDirectSol/norm2x - 1 ), 1e-12 );
  // Again we check a couple elements
//JAW  EXPECT_LT( ( std::fabs( solDirect.getElement( n ) - x.getElement( n ))/std::fabs( x.getElement( n ) ) ), 1e-12 );
//JAW  EXPECT_LT( ( std::fabs( solDirect.getElement( 3*n ) - x.getElement( 3*n ))/std::fabs( x.getElement( 3*n ) ) ), 1e-12 );

  // Test the clearRow function (this has to be done after the solves so that we can
  // use the same matrix when we are done with it).

  // We clear the row and multiply the diagonal value by 2.
#if 0
  testMatrix.clearRow( 2*N/4+n, 2.0 );
  testMatrix.getLocalRow(static_cast<LAI::lid>( n ),numValRown,vecValuesRown,vecIndicesRown);
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
  // These objects can use all of the available libraries (trilinos and
  // soon HYPRE and PETSc).
  using ParallelMatrix = typename LAI::ParallelMatrix;
  using ParallelVector = typename LAI::ParallelVector;
  using LAIgid = typename LAI::gid;

  // Get the MPI rank
  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  // Set the MPI communicator
  // MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm test_comm;
  MPI_Comm_dup( MPI_COMM_WORLD, &test_comm );
  MPI_Comm comm = test_comm;

  // Size of the dummy cartesian mesh we use to generate the Laplace 2D operator.
  LAIgid n = 100;
  // Number of degrees of freedom
  LAIgid N = n*n;

  // Compute a 2D Laplace operator
  ParallelMatrix testMatrix = compute2DLaplaceOperator<LAI>( comm, N );

  // Compute a dummy preconditioner (identity matrix)
  ParallelMatrix preconditioner = computeIdentity<LAI>( comm, N );

  /*
  // We first fill some standard vectors
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
  */

  // Define vectors for x, b and x0
  ParallelVector x, b, x0;

  // Place holder for b (will be recomputed a few lines down). We do use the vector to test
  // the linear algebra operations
  b.createWithGlobalSize( N, comm);
  //b.ones();

  // Random vector for x (the true solution of Ax = b).
  x.createWithGlobalSize( N, comm );
  x.rand();

  // Initial guess x0 = 0
  x0.createWithGlobalSize( N, comm );
  x0.zero();

  // Get the norm of the true solution
  real64 norm2x = x.norm2();

  // Perform a matrix/vector multiplication to compute b
  testMatrix.multiply( x, b );

  // Initial guess for CG
  ParallelVector solCG( x0 );

  // Initial guess for BiCGSTAB
  ParallelVector solBiCGSTAB( x0 );

  // Create a CG solver object
  CGsolver<LAI> testCG;

  // Solve the left-preconditioned system with CG
  testCG.solve( testMatrix, solCG, b, preconditioner );

  // The true solution is the vector x, so we check if the norm are equal.
  real64 normCG = solCG.norm2();
  EXPECT_LT( std::fabs( normCG/norm2x - 1. ), 5e-6 );

  // Create a BiCGSTAB solver object
  BiCGSTABsolver<LAI> testBiCGSTAB;

  // Solve the left-preconditioned system with BiCGSTAB.
  testBiCGSTAB.solve( testMatrix, solBiCGSTAB, b, preconditioner );

  // The true solution is the vector x, so we check if the norm are equal.
  real64 normBiCGSTAB = solBiCGSTAB.norm2();
  EXPECT_LT( std::fabs( normBiCGSTAB/norm2x - 1. ), 5e-6 );

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

  // Define aliases templated on the Linear Algebra Interface (LAI).
  // These objects can use all of the available libraries (trilinos and
  // soon HYPRE and PETSc).
  using ParallelMatrix = typename LAI::ParallelMatrix;
  using ParallelVector = typename LAI::ParallelVector;

  // Get the MPI rank
  typename LAI::LAI::lid rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  // Set the MPI communicator
  // MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm test_comm;
  MPI_Comm_dup( MPI_COMM_WORLD, &test_comm );
  MPI_Comm comm = test_comm;

  // Size of the dummy cartesian mesh we use to generate the Laplace 2D operator.
  LAI::gid n = 100;
  // Number of degrees of freedom
  LAI::gid N = n*n;

  // Compute a 2D Laplace operator
  ParallelMatrix matrix00 = compute2DLaplaceOperator<LAI>( comm, N );

  // Compute a dummy preconditioner (identity matrix)
  ParallelMatrix preconditioner00 = computeIdentity<LAI>( comm, N );

#if 0

  // We first fill some standard vectors
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

  // Define vectors for x, b and x0
  ParallelVector x, b, x0;

  // Right hand side for multiplication (b)
  b.create( zer );

  // Vector of ones for multiplication (x)
  x.create( random );

  // Random vector for initial guess (x0)
  x0.create( random2 );

//  // Matrix-vector product
//  matrix00.multiply( x, b );

  // Create vector for the true solution
  ParallelVector trueSolution( x );

  // Set the initial guesses for CG and BiCGSTAB
  ParallelVector solutionCG0( x0 );
  ParallelVector solutionBiCGSTAB0( x0 );
  ParallelVector solutionCG1( x0 );
  ParallelVector solutionBiCGSTAB1( x0 );

  // Set the right hand side vectors guesses for CG and BiCGSTAB
  ParallelVector rhs0( b );
  ParallelVector rhs1( b );

  // Size of the block system
  integer nRows = 2;
  integer nCols = 2;

  // Declare and allocate block matrices/vectors
  // System matrix
  BlockMatrixView<LAI> blockMatrix( nRows, nCols );

  // Preconditioner
  BlockMatrixView<LAI> blockPreconditioner( nRows, nCols );

  // True solution
  BlockVectorView<LAI> blockTrueSolution( nCols );

  // Solution from CG
  BlockVectorView<LAI> blockSolutionCG( nCols );

  // Solution from BiCGSTAB
  BlockVectorView<LAI> blockSolutionBiCGSTAB( nCols );

  // Right hand side (we only need one, it is a const input)
  BlockVectorView<LAI> blockRhs( nRows );

  // In this test we simply tile the laplace operator, so we assign a duplicate
  // of the monolithic matrix to every block of the matrix.
  // Block (0,0)
  blockMatrix.setBlock( 0, 0, matrix00 );

  // Block (0,1)
  blockMatrix.setBlock( 0, 1, matrix00 );

  // Block (1,0)
  blockMatrix.setBlock( 1, 0, matrix00 );

  // Block (1,1)
  blockMatrix.setBlock( 1, 1, matrix00 );

  // We do the same for the preconditioner to get the block identity matrix.
  // We ignore the off-diagonal blocks and leave them as null-pointers.
  // Block (0,0)
  blockPreconditioner.setBlock( 0, 0, preconditioner00 );
  // Block (1,1)
  blockPreconditioner.setBlock( 1, 1, preconditioner00 );

  // Construct the true solution
  blockTrueSolution.setBlock( 0, trueSolution );
  blockTrueSolution.setBlock( 1, trueSolution );

  // Set initial guess blocks (here we need multiple objects since we cannot
  // have them point to the same memory location). These objects are initial
  // guesses as input but solution vectors as output.
  // CG vector
  blockSolutionCG.setBlock( 0, solutionCG0 );
  blockSolutionCG.setBlock( 1, solutionCG1 );

  // BiCGSTAB vector
  blockSolutionBiCGSTAB.setBlock( 0, solutionBiCGSTAB0 );
  blockSolutionBiCGSTAB.setBlock( 1, solutionBiCGSTAB1 );

  // Set right hand side blocks (this can be one object, it is passed as const).
  // We keep 2 of them for potential tests with more flexibility (rhs0 != rhs1).
  blockRhs.setBlock( 0, rhs0 );
  blockRhs.setBlock( 1, rhs1 );

  blockMatrix.multiply( blockTrueSolution, blockRhs);

  // Get the norm of the true solution
  real64 norm2trueSolution;
  blockTrueSolution.norm2( norm2trueSolution );

  // Create block CG solver object
  CGsolver<LAI> testCG;
  // Solve the left-preconditioned block system with CG
  testCG.solve( blockMatrix, blockSolutionCG, blockRhs, blockPreconditioner );

  // The true solution is the vector x, so we check if the norm are equal.
  // Note: the tolerance is higher that in the previous cases because the matrix is
  // twice as big and the condition number is higher. See details on the error
  // bounds of Krylov methods wrt the condition number.
  real64 normCG;
  blockSolutionCG.norm2( normCG );
  EXPECT_LT( std::fabs( normCG/norm2trueSolution - 1. ), 5e-6 );

  // Create block BiCGSTAB solver object
  BiCGSTABsolver<LAI> testBiCGSTAB;
  // Solve the left-preconditioned block system with BiCGSTAB
  testBiCGSTAB.solve( blockMatrix, blockSolutionBiCGSTAB, blockRhs, blockPreconditioner );

  // The true solution is the vector x, so we check if the norm are equal.
  real64 normBiCGSTAB;
  blockSolutionBiCGSTAB.norm2( normBiCGSTAB );
  EXPECT_LT( std::fabs( normBiCGSTAB/norm2trueSolution - 1. ), 5e-6 );

#endif
  // Finalize MPI
  MPI_Finalize();
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

  testNativeSolvers<TrilinosInterface>();
  testGEOSXSolvers<TrilinosInterface>();
  testGEOSXBlockSolvers<TrilinosInterface>();

}

/*! @function testHypreLAOperations.
 * @brief Runs all tests using the Hypre interface.
 */
TEST(testLAOperations,testHypreLAOperations)
{
  //testLaplaceOperator<HypreInterface>();
  //testGEOSXSolvers<HypreInterface>();
  //testGEOSXBlockSolvers<HypreInterface>();
}

/*! @function testPETScLAOperations.
 * @brief Runs all tests using the PETSc interface.
 */
TEST(testLAOperations,testPETScLAOperations)
{

  //testLaplaceOperator<PETScInterface>();
  //testGEOSXSolvers<PETScInterface>();
  //testGEOSXBlockSolvers<PETScInterface>();
}

//@}

#ifdef __clang__
#pragma clang diagnostic pop
#endif
