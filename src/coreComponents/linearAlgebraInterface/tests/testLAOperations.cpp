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

/*
 * testLAOperations.hpp
 *
 *  Created on: Sep 19, 2018
 *      Author: Matthias
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
                                              typename LAI::laiGID N)
{
  // Declare matrix
  typename LAI::ParallelMatrix I;
  I.create( comm, N, 1 );

  // Construct a dummy Laplace matrix (5 points stencil)
  for (typename LAI::laiGID i = I.ilower(); i < I.iupper(); i++)
  {
    // Set the diagonal value to 1
    real64 temp = 1;
    // Set the values for row i
    I.insert( i, 1, &temp, &i );
  }

  // Close the matrix
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
                                                       typename LAI::laiGID N)
{
  typename LAI::ParallelMatrix laplace2D;
  laplace2D.create( comm, N, 5 );

  integer n = integer(std::sqrt(N));

  // Allocate arrays to fill dummy 2D Laplace (cartesian) matrix
  real64 values[5];
  typename LAI::laiGID cols[5];

  // Construct a dummy Laplace matrix (5 points stencil)
  for (typename LAI::laiGID i = laplace2D.ilower(); i < laplace2D.iupper(); i++)
  {
    integer nnz = 0;
    /* The left -n: position i-n */
    if (i-n >= 0)
    {
      cols[nnz] = i-n;
      values[nnz] = -1.0;
      nnz++;
    }
    /* The left -n: position i-n */
    if (i-1 >= 0)
    {
      cols[nnz] = i-1;
      values[nnz] = -1.0;
      nnz++;
    }
    /* Set the diagonal: position i */
    cols[nnz] = i;
    values[nnz] = 4.0;
    nnz++;
    /* The right +1: position i+1 */
    if (i+1 < N)
    {
      cols[nnz] = i+1;
      values[nnz] = -1.0;
      nnz++;
    }
    /* The right +n: position i-n */
    if (i+n < N)
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
 * @brief Test the native solvers from the LAI.
 */
template< typename LAI >
void testNativeSolvers()
{

  using ParallelMatrix = typename LAI::ParallelMatrix;
  using ParallelVector = typename LAI::ParallelVector;
  using LinearSolver = typename LAI::LinearSolver;
  using laiLID = typename LAI::laiLID;
  using laiGID = typename LAI::laiGID;

  MPI_Init( nullptr, nullptr );

  // Get the MPI rank
  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  // Set the MPI communicator
  // MPI_Comm comm = MPI_COMM_WORLD;

  MPI_Comm test_comm;
  MPI_Comm_dup( MPI_COMM_WORLD, &test_comm );
  MPI_Comm comm = test_comm;

  // Create Dummy Laplace matrix (5 points stencil)
  // Size of the mesh
  laiGID n = 100;
  // Size of the problem
  laiGID N = n*n;

  // Compute a 2D Laplace operator (symmetric) for testing purposes
  ParallelMatrix testMatrix = compute2DLaplaceOperator<TrilinosInterface>( comm, N );

  // Compute the identity to test the preconditioners
  ParallelMatrix preconditioner = computeIdentity<TrilinosInterface>( comm, N );

  integer numValRow0,numValRow1,numValRown;
  std::vector<real64> vecValuesRow0( 5 ),vecValuesRow1( 5 ),vecValuesRown( 5 );
  std::vector<laiLID> vecIndicesRow0( 5 ),vecIndicesRow1( 5 ),vecIndicesRown( 5 );
  testMatrix.getLocalRow( 0, numValRow0, vecValuesRow0, vecIndicesRow0 );
  testMatrix.getLocalRow( 1, numValRow1, vecValuesRow1, vecIndicesRow1 );
  testMatrix.getLocalRow( static_cast<laiLID>( n ), numValRown, vecValuesRown, vecIndicesRown );

  if (rank == 0)
  {
    // Check number of values per row
    EXPECT_TRUE( numValRow0 == 3 );
    EXPECT_TRUE( numValRow1 == 4 );
    EXPECT_TRUE( numValRown == 5 );

    // Check actual values TODO modify with better check for doubles
    //  EXPECT_TRUE( vecValuesRow0[0] == 4 );
    //  EXPECT_TRUE( vecValuesRow0[1] == -1 );
    //  EXPECT_TRUE( vecValuesRow0[2] == -1 );
    //
    //  EXPECT_TRUE( vecValuesRow1[0] == -1 );
    //  EXPECT_TRUE( vecValuesRow1[1] == 4 );
    //  EXPECT_TRUE( vecValuesRow1[2] == -1 );
    //  EXPECT_TRUE( vecValuesRow1[3] == -1 );
    //
    //  EXPECT_TRUE( vecValuesRown[0] == -1 );
    //  EXPECT_TRUE( vecValuesRown[1] == -1 );
    //  EXPECT_TRUE( vecValuesRown[2] == 4 );
    //  EXPECT_TRUE( vecValuesRown[3] == -1 );
    //  EXPECT_TRUE( vecValuesRown[4] == -1 );

    // Check column indices
    EXPECT_TRUE( vecIndicesRow0[0] == 0 );
    EXPECT_TRUE( vecIndicesRow0[1] == 1 );
    EXPECT_TRUE( vecIndicesRow0[2] == n );

    EXPECT_TRUE( vecIndicesRow1[0] == 0 );
    EXPECT_TRUE( vecIndicesRow1[1] == 1 );
    EXPECT_TRUE( vecIndicesRow1[2] == 2 );
    EXPECT_TRUE( vecIndicesRow1[3] == n+1 );

    EXPECT_TRUE( vecIndicesRown[0] == 0 );
    EXPECT_TRUE( vecIndicesRown[1] == n-1 );
    EXPECT_TRUE( vecIndicesRown[2] == n );
    EXPECT_TRUE( vecIndicesRown[3] == n+1 );
    EXPECT_TRUE( vecIndicesRown[4] == 2*n );
  }
  // Fill standard vectors
  std::vector<real64> ones, zer, random;
  for (integer j = 0; j < N; j++)
  {
    zer.push_back( 0 );
    ones.push_back( 1 );
    random.push_back( rand() % 20 - 10 );
  }

  // Define vectors
  ParallelVector x, b, init;
  // Right hand side for multiplication (b)
  b.create( zer );
  // Vector of ones for multiplication (x)
  x.create( ones );
  // Random vector for an initial guess.
  init.create( random );

  // Fill initial guesses
  ParallelVector solDirect( init );
  ParallelVector solIterative( init );
  ParallelVector solIterativeML( init );

  // Residual vector (zeros)
  ParallelVector r0( b );
  ParallelVector r1( b );

  // Matrix/vector multiplication to compute a known right hand side.
  testMatrix.multiply( x, b );

  // Test dot product
  real64 dotTest;
  x.dot( x, dotTest );

  // Test norms
  real64 norm1;
  x.norm1( norm1 );
  EXPECT_TRUE( std::fabs( norm1 - N ) <= 1e-6 );
  real64 norm2;
  x.norm2( norm2 );
  EXPECT_TRUE( std::fabs( norm2 - n ) <= 1e-6 );
  real64 norminf;
  x.normInf( norminf );
  EXPECT_TRUE( std::fabs( norminf - 1 ) <= 1e-6 );

  // Test residual
  testMatrix.residual( x, b, r0 );
  real64 normRes;
  r0.normInf( normRes );
  EXPECT_TRUE( std::fabs( normRes ) <= 1e-6 );

  // Test solvers
  LinearSolver solver = LinearSolver();

  // AztecOO iterative solver (TODO remove hard coded options)
  solver.solve( testMatrix, solIterative, b, 1000, 1e-8);
  real64 normIterativeSol;
  solIterative.normInf( normIterativeSol );
  EXPECT_TRUE( std::fabs( normIterativeSol - 1 ) <= 1e-5 );

  real64 norm2ItSol, norm2InitGuess = 0;
  solIterative.norm2( norm2ItSol );

  // AztecOO iterative solver (TODO remove hard coded options)
  // with ML_Epetra preconditioner (TODO remove hard coded options)
//  solver.ml_solve( testMatrix, solIterativeML, b, 500, 1e-8 );
//
//  real64 normIterativeSolML;
//  solIterativeML.normInf( normIterativeSolML );
//  EXPECT_TRUE( std::fabs( normIterativeSolML - 1 ) <= 1e-5 );

  // Amesos (Klu) direct solver
  solver.dsolve( testMatrix, solDirect, b );

  real64 normDirectSol;
  solDirect.normInf( normDirectSol );
  EXPECT_TRUE( std::fabs( normDirectSol - 1 ) <= 1e-8 );

  // Test clearRow
  testMatrix.clearRow( 2*N/4+n, 2.0 );
  testMatrix.getLocalRow(static_cast<laiLID>( n ),numValRown,vecValuesRown,vecIndicesRown);
  if ( rank == 2 )
  {
    EXPECT_TRUE( std::fabs( vecValuesRown[2] - 8.0 ) <= 1e-8 );
  }

  //MPI_Finalize();

}

/**
 * @function testGEOSXSolvers
 *
 * @brief Test the GEOSX solvers for monolithic matrices (mostly for debugging purposes).
 */
template< typename LAI >
void testGEOSXSolvers()
{

  using ParallelMatrix = typename LAI::ParallelMatrix;
  using ParallelVector = typename LAI::ParallelVector;
  using laiGID = typename LAI::laiGID;

//  MPI_Init( nullptr, nullptr );

  // Get the MPI rank
  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  // Set the MPI communicator
  // MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm test_comm;
  MPI_Comm_dup( MPI_COMM_WORLD, &test_comm );
  MPI_Comm comm = test_comm;

  // Create Dummy Laplace matrix (5 points stencil)
  // Size of the mesh
  laiGID n = 200;
  // Size of the problem
  laiGID N = n*n;

  // Compute a 2D Laplace operator (symmetric) for testing purposes
  ParallelMatrix testMatrix = compute2DLaplaceOperator<TrilinosInterface>( comm, N );

  // Compute the identity to test the preconditioners
  ParallelMatrix preconditioner = computeIdentity<TrilinosInterface>( comm, N );

  // Fill standard vectors
  std::vector<real64> ones, zer, random;
  for (integer j = 0; j < N; j++)
  {
    zer.push_back( 0 );
    ones.push_back( 1 );
    random.push_back( rand() % 20 - 10 );
  }

  // Define vectors
  ParallelVector x, b, init;
  // Right hand side for multiplication (b)
  b.create( zer );
  // Vector of ones for multiplication (x)
  x.create( ones );
  // Vector of random values for initial guess
  init.create( random );

  // Residual vectors
  ParallelVector rCG( b );
  ParallelVector rBiCGSTAB( b );

  // Matrix/vector multiplication
  testMatrix.multiply( x, b );

  // Initial guesses for CG and BiCGSTAB
  ParallelVector solCG( init );
  ParallelVector solBiCGSTAB( init );

  // Create CG solver object
  CGsolver<TrilinosInterface> testCG;
  // Solve the left-preconditioned system with CG
  testCG.solve( testMatrix, solCG, b, preconditioner );

  // Create BiCGSTAB solver object
  BiCGSTABsolver<TrilinosInterface> testBiCGSTAB;
  // Solve the left-preconditioned system with BiCGSTAB
  testBiCGSTAB.solve( testMatrix, solBiCGSTAB, b, preconditioner );

  //MPI_Finalize();

}

/**
 * @function testGEOSXBlockSolvers
 *
 * @brief Test the GEOSX block solvers.
 */
template< typename LAI >
void testGEOSXBlockSolvers()
{

  using ParallelMatrix = typename LAI::ParallelMatrix;
  using ParallelVector = typename LAI::ParallelVector;

  // MPI_Init(nullptr,nullptr);

  // Get the MPI rank
  typename LAI::laiLID rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  // Set the MPI communicator
  // MPI_Comm comm = MPI_COMM_WORLD;

  MPI_Comm test_comm;
  MPI_Comm_dup( MPI_COMM_WORLD, &test_comm );
  MPI_Comm comm = test_comm;

  // Create Dummy Laplace matrix (5 points stencil)
  typename LAI::laiGID n = 200;
  typename LAI::laiGID N = n*n;

  // Compute a 2D Laplace operator (symmetric) for testing purposes
  ParallelMatrix matrix00 = compute2DLaplaceOperator<TrilinosInterface>( comm, N );

  // Compute the identity to test the preconditioners
  ParallelMatrix preconditioner00 = computeIdentity<TrilinosInterface>( comm, N );

  // Fill standard vectors
  std::vector<real64> ones, zer, random;
  for (integer j = 0; j < N; j++)
  {
    zer.push_back( 0 );
    ones.push_back( 1 );
    random.push_back( rand() % 20 - 10 );
  }

  // Define vectors
  ParallelVector x, b, init;
  // Right hand side for multiplication (b)
  b.create( zer );
  // Vector of ones for multiplication (x)
  x.create( ones );
  init.create( random );

  // Matrix-vector product
  matrix00.multiply( x, b );

  // Set the initial guesses for CG and BiCGSTAB
  ParallelVector solutionCG0( init );
  ParallelVector solutionBiCGSTAB0( init );
  ParallelVector solutionCG1( init );
  ParallelVector solutionBiCGSTAB1( init );

  // Set the right hand side vectors guesses for CG and BiCGSTAB
  ParallelVector rhs0( b );
  ParallelVector rhs1( b );

  // Size of the block system
  integer nRows = 2;
  integer nCols = 2;

  // Declare and allocate block matrices/vectors
  BlockMatrixView<TrilinosInterface> blockMatrix( nRows, nCols );
  BlockMatrixView<TrilinosInterface> blockPreconditioner( nRows, nCols );
  BlockVectorView<TrilinosInterface> blockSolutionCG( nCols );
  BlockVectorView<TrilinosInterface> blockSolutionBiCGSTAB( nCols );
  BlockVectorView<TrilinosInterface> blockRhs( nRows );

  // Duplicate matrices if needed
  ParallelMatrix matrix01( matrix00 );
  ParallelMatrix matrix10( matrix00 );
  ParallelMatrix matrix11( matrix00 );
  ParallelMatrix preconditioner11( preconditioner00 );

  // Set matrix blocks (tiled Laplace operators)
  blockMatrix.setBlock( 0, 0, matrix00 );
  blockMatrix.setBlock( 0, 1, matrix01 );
  blockMatrix.setBlock( 1, 0, matrix10 );
  blockMatrix.setBlock( 1, 1, matrix11 );

  // Set preconditioner blocks (diagonal identity)
  blockPreconditioner.setBlock( 0, 0, preconditioner00 );
  blockPreconditioner.setBlock( 1, 1, preconditioner11 );

  // Set initial guess blocks
  blockSolutionCG.setBlock( 0, solutionCG0 );
  blockSolutionCG.setBlock( 1, solutionCG1 );
  blockSolutionBiCGSTAB.setBlock( 0, solutionBiCGSTAB0 );
  blockSolutionBiCGSTAB.setBlock( 1, solutionBiCGSTAB1 );

  // Set right hand side blocks
  blockRhs.setBlock( 0, rhs0 );
  blockRhs.setBlock( 1, rhs1 );

  // Create block CG solver object
  CGsolver<TrilinosInterface> testCG;
  // Solve the left-preconditioned block system with CG
  testCG.solve( blockMatrix, blockSolutionCG, blockRhs, blockPreconditioner );

  // Create block BiCGSTAB solver object
  BiCGSTABsolver<TrilinosInterface> testBiCGSTAB;
  // Solve the left-preconditioned block system with BiCGSTAB
  testBiCGSTAB.solve( blockMatrix, blockSolutionBiCGSTAB, blockRhs, blockPreconditioner );

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
