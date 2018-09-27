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
 * @file testLAOperations.hpp
 *
 *  Created on: Sep 19, 2018
 *      Author: Matthias Cremon
 */

#include "gtest/gtest.h"

#include <iostream>
#include <vector>
#include <mpi.h>

#include "TrilinosInterface.hpp"
//#include "HypreInterface.hpp"
#include "BlockMatrixView.hpp"
#include "CGsolver.hpp"
#include "BiCGSTABsolver.hpp"

#include "common/DataTypes.hpp"

/**
 * \file testLAOperations.cpp
 * \brief This test file is part of the ctest suite and tests the Trilinos based solvers
 * as well as block GEOSX solvers. It mainly uses dummy 2D Laplace operator matrices
 * to test the solvers along with a simple (block) identity precondioner.
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
 * @brief Compute the identity matrix of size <tt>N</tt>.
 *
 * \param comm MPI communicator.
 * \param N size of the squared identity matrix.
 */
template< typename LAI >
typename LAI::ParallelMatrix computeIdentity( MPI_Comm comm,
                                              typename LAI::laiGID N )
{
  // Declare matrix
  typename LAI::ParallelMatrix I;
  // Create a matrix of size N with 1 non-zero per row
  I.create( comm, N, 1 );

  // Loop over rows to fill the matrix
  for ( typename LAI::laiGID i = I.ilower(); i < I.iupper(); i++ )
  {
    // Placeholder for the 1 value (we need to take the memory address)
    real64 temp = 1;
    // Set the value for element (i,i) to 1
    I.insert( i, 1, &temp, &i );
  }

  // Close the matrix (make data contiguous in memory)
  I.close();

  return I;
}

/**
 * @brief Compute the 2D Laplace operator of size <tt>N</tt>.
 *
 * \param comm MPI communicator.
 * \param N size of the squared 2D Laplace operator matrix.
 */
template< typename LAI >
typename LAI::ParallelMatrix compute2DLaplaceOperator( MPI_Comm comm,
                                                       typename LAI::laiGID N )
{
  // Declare matrix
  typename LAI::ParallelMatrix laplace2D;
  // Create a matrix of size N with 5 non-zeros per row
  laplace2D.create( comm, N, 5 );

  // Get the size of the dummy mesh back to be able to put values in the correct
  // diagonals.
  integer n = integer( std::sqrt( N ) );

  // Allocate arrays to fill the matrix (values and columns)
  real64 values[5];
  typename LAI::laiGID cols[5];

  // Loop over rows to fill the matrix
  for ( typename LAI::laiGID i = laplace2D.ilower(); i < laplace2D.iupper(); i++ )
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

    // The left -n: position i-n
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

    // The right +n: position i-n
    if ( i+n < N )
    {
      cols[nnz] = i+n;
      values[nnz] = -1.0;
      nnz++;
    }

    // Set the values for row i
    laplace2D.insert( i, nnz, values, cols );

  }

  // Close the matrix
  laplace2D.close();

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
template< typename LAI >
void testNativeSolvers()
{
  // Define aliases templated on the Linear Algebra Interface (LAI).
  // These objects can use all of the available libraries (trilinos and
  // soon HYPRE and PETSc).
  using ParallelMatrix = typename LAI::ParallelMatrix;
  using ParallelVector = typename LAI::ParallelVector;
  using LinearSolver = typename LAI::LinearSolver;
  using laiLID = typename LAI::laiLID;
  using laiGID = typename LAI::laiGID;

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
  laiGID n = 100;
  // Number of degrees of freedom
  laiGID N = n*n;

  // Compute a 2D Laplace operator
  ParallelMatrix testMatrix = compute2DLaplaceOperator<LAI>( comm, N );

  // Compute a dummy preconditioner (identity matrix)
  ParallelMatrix preconditioner = computeIdentity<LAI>( comm, N );

  // Declare integers to run tests
  integer numValRow0,numValRow1,numValRown;

  // Declare vectors to run tests
  std::vector<real64> vecValuesRow0( 5 ),vecValuesRow1( 5 ),vecValuesRown( 5 );
  std::vector<laiLID> vecIndicesRow0( 5 ),vecIndicesRow1( 5 ),vecIndicesRown( 5 );

  // Get values and columns in specific rows
  testMatrix.getLocalRow( 0, numValRow0, vecValuesRow0, vecIndicesRow0 );
  testMatrix.getLocalRow( 1, numValRow1, vecValuesRow1, vecIndicesRow1 );
  testMatrix.getLocalRow( static_cast<laiLID>( n ), numValRown, vecValuesRown, vecIndicesRown );

  // Run checks on rank 0 to see if the matrix was correctly constructed
  if (rank == 0)
  {
    // Check number of values per row
    EXPECT_TRUE( numValRow0 == 3 );
    EXPECT_TRUE( numValRow1 == 4 );
    EXPECT_TRUE( numValRown == 5 );

    // Check column indices for row 0
    EXPECT_TRUE( vecIndicesRow0[0] == 0 );
    EXPECT_TRUE( vecIndicesRow0[1] == 1 );
    EXPECT_TRUE( vecIndicesRow0[2] == n );

    // Check column indices for row 1
    EXPECT_TRUE( vecIndicesRow1[0] == 0 );
    EXPECT_TRUE( vecIndicesRow1[1] == 1 );
    EXPECT_TRUE( vecIndicesRow1[2] == 2 );
    EXPECT_TRUE( vecIndicesRow1[3] == n+1 );

    // Check column indices for row n
    EXPECT_TRUE( vecIndicesRown[0] == 0 );
    EXPECT_TRUE( vecIndicesRown[1] == n-1 );
    EXPECT_TRUE( vecIndicesRown[2] == n );
    EXPECT_TRUE( vecIndicesRown[3] == n+1 );
    EXPECT_TRUE( vecIndicesRown[4] == 2*n );
  }

  // We now test the linear algebra operations, such as matrix-vector product and
  // compute residuals.

  // For that we first fill some standard vectors
  std::vector<real64> ones, zer, random;
  for (integer j = 0; j < N; j++)
  {
    // Vector of zeros
    zer.push_back( 0 );
    // Vector of ones
    ones.push_back( 1 );
    // Vector of random integers
    random.push_back( rand() % 20 - 10 );
  }

  // Define x, b and x0 guess vectors
  ParallelVector x, b, x0;

  // Right hand side for multiplication (b)
  b.create( zer );

  // Vector of ones for multiplication (x)
  x.create( ones );

  // Random vector for the initial guess.
  x0.create( random );

  // Fill initial guess for the direct solver
  ParallelVector solDirect( x0 );

  // Fill initial guess for the iterative solver
  ParallelVector solIterative( x0 );

  // Residual vector (vector of zeros)
  ParallelVector r( b );

  // Compute the matrix/vector multiplication.
  testMatrix.multiply( x, b );

  // Test dot product
  real64 dotTest;
  x.dot( x, dotTest );

  // Test 1-norm
  real64 norm1;
  x.norm1( norm1 );
  // The 1-norm should be equal to N
  EXPECT_TRUE( std::fabs( norm1 - N ) <= 1e-6 );

  // Test 2-norm
  real64 norm2;
  x.norm2( norm2 );
  // The 2-norm should be equal to n
  EXPECT_TRUE( std::fabs( norm2 - n ) <= 1e-6 );

  // Test inf-norm
  real64 norminf;
  x.normInf( norminf );
  // The inf-norm should be equal to 1
  EXPECT_TRUE( std::fabs( norminf - 1 ) <= 1e-6 );

  // Test residual function by computing r = b - Ax.
  testMatrix.residual( x, b, r );
  real64 normRes;
  r.normInf( normRes );
  // The inf-norm should be equal to 0 (all elements are 0).
  EXPECT_TRUE( std::fabs( normRes ) <= 1e-6 );

  // We now test the linear solvers from the libraries.

  // Declare solver object.
  LinearSolver solver = LinearSolver();

  // Run the iterative solver (TODO remove hard coded options, currently
  // using ILUT preconditioned GMRES).
  if (rank == 0)
    std::cout << std::endl << "AztecOO iterative solver:";
  solver.solve( testMatrix, solIterative, b, 1000, 1e-8);

  // Get the inf-norm of the solution.
  real64 normIterativeSol;
  solIterative.normInf( normIterativeSol );

  // The true solution is a vector of ones, so the inf-norm should be 1.
  EXPECT_TRUE( std::fabs( normIterativeSol - 1 ) <= 5e-5 );

  // We now do the same using a direct solver from Amesos (Klu)
  if (rank == 0)
    std::cout << std::endl << "Amesos direct solver:" << std::endl << std::endl;
  solver.dsolve( testMatrix, solDirect, b );

  // Again the inf-norm should be 1
  real64 normDirectSol;
  solDirect.normInf( normDirectSol );
  EXPECT_TRUE( std::fabs( normDirectSol - 1 ) <= 1e-8 );

  // Test the clearRow function (this has to be after the solves so that we do not
  // use a different matrix.

  // We clear the row and multiply the diagonal value by 2.
  testMatrix.clearRow( 2*N/4+n, 2.0 );
  testMatrix.getLocalRow(static_cast<laiLID>( n ),numValRown,vecValuesRown,vecIndicesRown);
  if ( rank == 2 )
  {
    EXPECT_TRUE( std::fabs( vecValuesRown[2] - 8.0 ) <= 1e-6 );
  }

  //MPI_Finalize();

}

/**
 * @function testGEOSXSolvers
 *
 * @brief Test the GEOSX solvers for monolithic matrices by solving a system with
 * a Laplace operator and the identity matrix as a preconditioner.
 */
template< typename LAI >
void testGEOSXSolvers()
{
  // Define aliases templated on the Linear Algebra Interface (LAI).
  // These objects can use all of the available libraries (trilinos and
  // soon HYPRE and PETSc).
  using ParallelMatrix = typename LAI::ParallelMatrix;
  using ParallelVector = typename LAI::ParallelVector;
  using laiGID = typename LAI::laiGID;

  // Get the MPI rank
  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  // Set the MPI communicator
  // MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm test_comm;
  MPI_Comm_dup( MPI_COMM_WORLD, &test_comm );
  MPI_Comm comm = test_comm;

  // Size of the dummy cartesian mesh we use to generate the Laplace 2D operator.
  laiGID n = 100;
  // Number of degrees of freedom
  laiGID N = n*n;

  // Compute a 2D Laplace operator
  ParallelMatrix testMatrix = compute2DLaplaceOperator<LAI>( comm, N );

  // Compute a dummy preconditioner (identity matrix)
  ParallelMatrix preconditioner = computeIdentity<LAI>( comm, N );

  // We first fill some standard vectors
  std::vector<real64> ones, zer, random;
  for (integer j = 0; j < N; j++)
  {
    // Vector of zeros
    zer.push_back( 0 );
    // Vector of ones
    ones.push_back( 1 );
    // Vector of random integers
    random.push_back( rand() % 20 - 10 );
  }

  // Define vectors for x, b and x0
  ParallelVector x, b, x0;

  // Right hand side for multiplication (b)
  b.create( zer );

  // Vector of ones for multiplication (x)
  x.create( ones );

  // Vector of random values for initial guesses
  x0.create( random );

  // Residual vector for CG
  ParallelVector rCG( b );

  // Residual vector for BiCGSTAB
  ParallelVector rBiCGSTAB( b );

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

  // Check if the inf norm is 1 (expected solution is a vector of 1).
  real64 normCG;
  solCG.normInf( normCG );
  EXPECT_TRUE( std::fabs( normCG - 1 ) <= 1e-6 );

  // Create a BiCGSTAB solver object
  BiCGSTABsolver<LAI> testBiCGSTAB;

  // Solve the left-preconditioned system with BiCGSTAB.
  testBiCGSTAB.solve( testMatrix, solBiCGSTAB, b, preconditioner );

  // Check if the inf norm is 1 (expected solution is a vector of 1).
  real64 normBiCGSTAB;
  solBiCGSTAB.normInf( normBiCGSTAB );
  EXPECT_TRUE( std::fabs( normBiCGSTAB - 1 ) <= 1e-6 );

  //MPI_Finalize();

}

/**
 * @function testGEOSXBlockSolvers
 *
 * @brief Test the GEOSX block solvers by solving a system with a block matrix
 * made of tiled Laplace operators and using the identity as a preconditioner.
 */
template< typename LAI >
void testGEOSXBlockSolvers()
{
  // Define aliases templated on the Linear Algebra Interface (LAI).
  // These objects can use all of the available libraries (trilinos and
  // soon HYPRE and PETSc).
  using ParallelMatrix = typename LAI::ParallelMatrix;
  using ParallelVector = typename LAI::ParallelVector;
  using laiGID = typename LAI::laiGID;

  // MPI_Init(nullptr,nullptr);

  // Get the MPI rank
  typename LAI::laiLID rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  // Set the MPI communicator
  // MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm test_comm;
  MPI_Comm_dup( MPI_COMM_WORLD, &test_comm );
  MPI_Comm comm = test_comm;

  // Size of the dummy cartesian mesh we use to generate the Laplace 2D operator.
  laiGID n = 100;
  // Number of degrees of freedom
  laiGID N = n*n;

  // Compute a 2D Laplace operator
  ParallelMatrix matrix00 = compute2DLaplaceOperator<LAI>( comm, N );

  // Compute a dummy preconditioner (identity matrix)
  ParallelMatrix preconditioner00 = computeIdentity<LAI>( comm, N );

  // We first fill some standard vectors
  std::vector<real64> ones, zer, random;
  for (integer j = 0; j < N; j++)
  {
    // Vector of zeros
    zer.push_back( 0 );
    // Vector of ones
    ones.push_back( 1 );
    // Vector of random integers
    random.push_back( rand() % 20 - 10 );
  }

  // Define vectors for x, b and x0
  ParallelVector x, b, x0;

  // Right hand side for multiplication (b)
  b.create( zer );

  // Vector of ones for multiplication (x)
  x.create( ones );

  // Random vector for initial guess (x0)
  x0.create( random );

  // Matrix-vector product
  matrix00.multiply( x, b );

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
  blockPreconditioner.setBlock( 0, 0, preconditioner00 );
  blockPreconditioner.setBlock( 1, 1, preconditioner00 );

  // Set initial guess blocks (here we need multiple objects since we cannot
  // have them point to the same memory location. These objects are initial
  // guesses as input but solution vectors as output.
  // CG vector
  blockSolutionCG.setBlock( 0, solutionCG0 );
  blockSolutionCG.setBlock( 1, solutionCG1 );

  // BiCGSTAB vector
  blockSolutionBiCGSTAB.setBlock( 0, solutionBiCGSTAB0 );
  blockSolutionBiCGSTAB.setBlock( 1, solutionBiCGSTAB1 );

  // Set right hand side blocks (this can be one object, it is passed as const.
  // We keep 2 of them for potential tests with more flexibility (rhs0 != rhs1).
  blockRhs.setBlock( 0, rhs0 );
  blockRhs.setBlock( 1, rhs1 );

  // Create block CG solver object
  CGsolver<LAI> testCG;
  // Solve the left-preconditioned block system with CG
  testCG.solve( blockMatrix, blockSolutionCG, blockRhs, blockPreconditioner );

  // Check if the inf norm is 0.5 (the expected solution is a vector of 0.5).
  real64 normCG;
  blockSolutionCG.normInf( normCG );
  EXPECT_TRUE( std::fabs( normCG - 0.5 ) <= 1e-6 );

  // Create block BiCGSTAB solver object
  BiCGSTABsolver<LAI> testBiCGSTAB;
  // Solve the left-preconditioned block system with BiCGSTAB
  testBiCGSTAB.solve( blockMatrix, blockSolutionBiCGSTAB, blockRhs, blockPreconditioner );

  // Check if the inf norm is 0.5 (the expected solution is a vector of 0.5).
  real64 normBiCGSTAB;
  blockSolutionBiCGSTAB.normInf( normBiCGSTAB );
  EXPECT_TRUE( std::fabs( normBiCGSTAB - 0.5 ) <= 1e-6 );

  MPI_Finalize();

}
//@}

/*! @name Ctest tests.
 * @brief Runs similar testing functions using different Linear Algebra Interfaces (LAIs).
 */
//@{

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wglobal-constructors"
#pragma clang diagnostic ignored "-Wexit-time-destructors"
#pragma clang diagnostic ignored "-Wzero-as-null-pointer-constant"
#endif

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

}

//@}

#ifdef __clang__
#pragma clang diagnostic pop
#endif
