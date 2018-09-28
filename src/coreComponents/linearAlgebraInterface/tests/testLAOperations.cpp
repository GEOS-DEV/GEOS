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

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wglobal-constructors"
#pragma clang diagnostic ignored "-Wexit-time-destructors"
#if __clang_major__ >= 5
#pragma clang diagnostic ignored "-Wzero-as-null-pointer-constant"
#endif
#endif

#include <gtest/gtest.h>


#include <iostream>
#include <vector>
#include <mpi.h>

#include "TrilinosInterface.hpp"
//#include "HypreInterface.hpp"

#include "common/DataTypes.hpp"

using namespace geosx;

template< typename LAI >
void testLaplaceOperator()
{

  using ParallelMatrix = typename LAI::ParallelMatrix;
  using ParallelVector = typename LAI::ParallelVector;
  using LinearSolver = typename LAI::LinearSolver;
  using laiLID = typename LAI::laiLID;
  using laiGID = typename LAI::laiGID;

  MPI_Init(nullptr,nullptr);

  // Get the MPI rank
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // Set the MPI communicator
  // MPI_Comm comm = MPI_COMM_WORLD;

  MPI_Comm test_comm;
  MPI_Comm_dup(MPI_COMM_WORLD,&test_comm);
  MPI_Comm comm = test_comm;

  // Create Dummy Laplace matrix (5 points stencil)
  laiGID n = 100;
  laiGID N = n*n;

  ParallelMatrix testMatrix;
  testMatrix.create(comm,N,5);

  ParallelMatrix preconditioner;
  preconditioner.create(comm,N,1);

  // Allocate arrays to fill dummy 2D Laplace (cartesian) matrix
  real64 values[5];
  laiGID cols[5];

  // Construct a dummy Laplace matrix (5 points stencil)
  for (laiGID i = testMatrix.ilower(); i < testMatrix.iupper(); i++)
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

    real64 temp = 0.25;

    // Set the values for row i
    testMatrix.insert(i,nnz,values,cols);
    preconditioner.insert(i,1,&temp,&i);

  }

  testMatrix.close();
  preconditioner.close();

  integer numValRow0,numValRow1,numValRown;
  std::vector<real64> vecValuesRow0(5),vecValuesRow1(5),vecValuesRown(5);
  std::vector<laiLID> vecIndicesRow0(5),vecIndicesRow1(5),vecIndicesRown(5);
  testMatrix.getLocalRow(0,numValRow0,vecValuesRow0,vecIndicesRow0);
  testMatrix.getLocalRow(1,numValRow1,vecValuesRow1,vecIndicesRow1);
  testMatrix.getLocalRow(static_cast<laiLID>(n),numValRown,vecValuesRown,vecIndicesRown);

  if (rank == 0)
  {
    // Check number of values per row
    EXPECT_TRUE( numValRow0 == 3 );
    EXPECT_TRUE( numValRow1 == 4 );
    EXPECT_TRUE( numValRown == 5 );

    // Check actual values TODO modify with better double check
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
  std::vector<real64> ones, zer;
  for (integer j = 0; j < N; j++)
  {
    zer.push_back(0);
    ones.push_back(1);
  }

  // Define vectors
  ParallelVector x, b;
  // Right hand side for multiplication (b)
  b.create(zer);
  // Vector of ones for multiplication (x)
  x.create(ones);
  // Vector of zeros for iterative and direct solutions
  ParallelVector solIterativeML(b);
  ParallelVector solDirect(b);
  // Residual vector
  ParallelVector r(b);

  // Matrix/vector multiplication
  testMatrix.multiply(x, b);

  ParallelVector solIterative(b);

  // Test dot product
  real64 dotTest;
  x.dot(x,&dotTest);

  // test norms
  real64 norm1;
  x.norm1(norm1);
  EXPECT_TRUE( std::fabs(norm1 - N) <= 1e-6 );

  real64 norm2;
  x.norm2(norm2);
  EXPECT_TRUE( std::fabs(norm2 - n) <= 1e-6 );

  real64 norminf;
  x.normInf(norminf);
  EXPECT_TRUE( std::fabs(norminf - 1) <= 1e-6 );

  testMatrix.residual(x,b,r);
  real64 normRes;
  r.normInf(normRes);
  EXPECT_TRUE( std::fabs(normRes) <= 1e-6 );

  // Test solvers
  LinearSolver solver = LinearSolver();

  // AztecOO iterative solver
  solver.solve(testMatrix,solIterative,b,1000,1e-8);
  real64 normIterativeSol;
  solIterative.normInf(normIterativeSol);
  EXPECT_TRUE( std::fabs(normIterativeSol - 1) <= 1e-5 );

  real64 norm2ItSol, norm2InitGuess = 0;
  solIterative.norm2(norm2ItSol);

  solver.ml_solve(testMatrix,solIterativeML,b,500,1e-8);

  real64 normIterativeSolML;
  solIterativeML.normInf(normIterativeSolML);
  EXPECT_TRUE( std::fabs(normIterativeSolML - 1) <= 1e-5 );

  // Amesos (Klu) direct solver
  solver.dsolve(testMatrix,solDirect,b);

  real64 normDirectSol;
  solDirect.normInf(normDirectSol);
  EXPECT_TRUE( std::fabs(normDirectSol - 1) <= 1e-8 );

  // Test clearRow
  testMatrix.clearRow(2*N/4+n,2.0);
  testMatrix.getLocalRow(static_cast<laiLID>(n),numValRown,vecValuesRown,vecIndicesRown);
  if (rank == 2)
  {
    EXPECT_TRUE( std::fabs(vecValuesRown[2] - 8.0) <= 1e-8 );
  }

  MPI_Finalize();

}

TEST(testLAOperations,testEpetraLAOperations)
{

  testLaplaceOperator<TrilinosInterface>();

}

TEST(testLAOperations,testHypreLAOperations)
{

  //testLaplaceOperator<HypreInterface>();

}

#ifdef __clang__
#pragma clang diagnostic pop
#endif
