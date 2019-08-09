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
 * @file testLAOperations.cpp
 */

#include <gtest/gtest.h>

#include <iostream>
#include <mpi.h>

#include "common/DataTypes.hpp"
#include "managers/initialization.hpp"

#include "TrilinosInterface.hpp"
#include "PetscInterface.hpp"
//#include "HypreInterface.hpp"
#include "LinearSolverParameters.hpp"

#include "NativeCG.hpp"
#include "NativeBICGSTAB.hpp"
#include "BlockMatrixView.hpp"
#include "BlockVectorView.hpp"

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
template<typename LAI>
typename LAI::ParallelMatrix computeIdentity( MPI_Comm comm,
                                              globalIndex N )
{
  // Declare matrix
  typename LAI::ParallelMatrix I;

  // Create a matrix of size N with 1 non-zero per row
  I.createWithGlobalSize( N, 1, comm );

  // Loop over rows to fill the matrix
  for( globalIndex i = I.ilower() ; i < I.iupper() ; i++ )
  {
    // Set the value for element (i,i) to 1
    I.insert( i, i, 1.0 );
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
 * \param n size of the nxn mesh for the square 2D Laplace operator matrix. Matrix size will be N=n^2.
 */

// ==============================
// Compute 2D Laplace Operator
// ==============================
// This function computes the matrix corresponding to a 2D Laplace operator. These
// matrices arise from a classical finite volume formulation on a cartesian mesh
// (5-point stencil).  Input is the mesh size, n, from which the total dofs is N = n^2;
template<typename LAI>
typename LAI::ParallelMatrix compute2DLaplaceOperator( MPI_Comm comm,
                                                       globalIndex n )
{
  // total dofs = n^2
  globalIndex N = n * n;

  // Declare matrix
  typename LAI::ParallelMatrix laplace2D;

  // Create a matrix of global size N with 5 non-zeros per row
  laplace2D.createWithGlobalSize( N, 5, comm );

  // Allocate arrays to fill the matrix (values and columns)
  real64 values[5];
  globalIndex cols[5];

  // Loop over rows to fill the matrix
  for( globalIndex i = laplace2D.ilower() ; i < laplace2D.iupper() ; i++ )
  {
    // Re-set the number of non-zeros for row i to 0.
    localIndex nnz = 0;

    // The left -n: position i-n
    if( i - n >= 0 )
    {
      cols[nnz] = i - n;
      values[nnz] = -1.0;
      nnz++;
    }

    // The left -1: position i-1
    if( i - 1 >= 0 )
    {
      cols[nnz] = i - 1;
      values[nnz] = -1.0;
      nnz++;
    }

    // Set the diagonal: position i
    cols[nnz] = i;
    values[nnz] = 4.0;
    nnz++;

    // The right +1: position i+1
    if( i + 1 < N )
    {
      cols[nnz] = i + 1;
      values[nnz] = -1.0;
      nnz++;
    }

    // The right +n: position i+n
    if( i + n < N )
    {
      cols[nnz] = i + n;
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
 * @function testVectorFunction
 *
 * @brief Test vector functions including create, add/set, accessors,
 * and linear algebra operations.
 */
// -----------------------------------------
// Test vector functions
// -----------------------------------------
template<typename LAI>
void testVectorFunctions()
{
  // Define alias
  using ParallelVector = typename LAI::ParallelVector;

  // Get the MPI rank
  int rank, numRanks;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  MPI_Comm_size( MPI_COMM_WORLD, &numRanks );

  ParallelVector x;
  localIndex const localSize = 3;
  globalIndex const globalSize = localSize * numRanks;
  globalIndex const offset = rank * localSize;

  // Testing createWithLocalSize
  x.createWithLocalSize( localSize, MPI_COMM_WORLD );
  EXPECT_EQ( x.localSize(), localSize );
  EXPECT_EQ( x.globalSize(), globalSize );

  // Testing iupper/ilower
  EXPECT_EQ( x.ilower(), offset );
  EXPECT_EQ( x.iupper(), offset + localSize );

  // Testing setting/getting values locally
  for( globalIndex i = x.ilower(); i < x.iupper(); ++i )
  {
    x.set( i, 2 * i );
  }
  x.close();
  for( globalIndex i = x.ilower(); i < x.iupper(); ++i )
  {
    EXPECT_DOUBLE_EQ( x.get( i ), 2 * i );
  }

  // Testing createWithGlobalSize
  x.createWithGlobalSize( globalSize, MPI_COMM_WORLD );
  EXPECT_EQ( x.localSize(), localSize );
  EXPECT_EQ( x.globalSize(), globalSize );

  // Testing setting global values on rank 0 and getting locally
  if( rank == 0 )
  {
    for( globalIndex i = 0; i < x.globalSize(); ++i )
    {
      x.set( i, 2 * i );
    }
  }
  x.close();
  for( globalIndex i = x.ilower(); i < x.iupper(); ++i )
  {
    EXPECT_EQ( x.get( i ), 2 * i );
  }

  // Testing getLocalRowID
  for( globalIndex i = x.ilower(); i < x.iupper(); ++i )
  {
    EXPECT_EQ( x.getLocalRowID( i ), i % localSize );
  }

  // Testing create with array1d
  array1d< real64 > localVals( localSize );
  for( localIndex i = 0; i < localSize; ++i )
  {
    localVals[i] = real64( i + rank * localSize );
  }

  ParallelVector v;
  v.create( localVals, MPI_COMM_WORLD );
  for( globalIndex i = v.ilower(); i < v.iupper(); ++i )
  {
    EXPECT_EQ( v.get( i ), localVals[v.getLocalRowID( i )] );
  }

  // Testing copy constructor, create with ParallelVector, 
  // get element
  ParallelVector y( x );
  ParallelVector z;
  z.create( x );
  for( globalIndex i = x.ilower(); i < x.iupper(); ++i )
  {
    EXPECT_EQ( x.get( i ), y.get( i ) );
    EXPECT_EQ( x.get( i ), z.get( i ) );
  }

  // Testing zero
  z.zero();
  for( globalIndex i = y.ilower(); i < y.iupper(); ++i )
  {
    EXPECT_EQ( z.get( i ), 0 );
  }

  // Testing copy
  z.copy( x );
  for( globalIndex i = y.ilower(); i < y.iupper(); ++i )
  {
    EXPECT_EQ( x.get( i ), z.get( i ) );
  }

  // Testing scale, z = x
  z.scale( 4.0 );
  for( globalIndex i = y.ilower(); i < y.iupper(); ++i )
  {
    EXPECT_EQ( 4.0 * x.get( i ), z.get( i ) );
  }

  // Testing add/set single element
  x.set( offset, -1 );
  x.close(); // set/add can't be interchanged
  x.add( offset + 1, 10 );
  x.close();
  EXPECT_DOUBLE_EQ( x.get( offset ), -1 );
  EXPECT_DOUBLE_EQ( x.get( offset + 1 ), offset * 2 + 12 );

  // Testing add/set c-style
  {
    globalIndex const inds[3] = { offset, offset + 1, offset + 2 };
    real64 const vals[3] = { -5.0, -6.0, 0.0 };
    y.zero();
    y.set( inds, vals, 2 );
    y.close();
    z.set( 1.0 );
    z.add( inds, vals, 2 );
    z.close();
    for( localIndex i = 0; i < 3; ++i )
    {
      EXPECT_DOUBLE_EQ( y.get( inds[i] ), vals[i] );
      EXPECT_DOUBLE_EQ( z.get( inds[i] ), vals[i] + 1.0 );
    }
  }

  // Testing add/set array1d-style
  {
    array1d< globalIndex > inds( 3 );
    inds[0] = offset; inds[1] = offset + 1; inds[2] = offset + 2;
    array1d< real64 > vals( 3 );
    vals[0] = -5.0; vals[1] = -6.0; vals[2] = 0.0;
    y.zero();
    y.set( inds, vals );
    y.close();
    z.set( 1.0 );
    z.add( inds, vals );
    z.close();
    for( localIndex i = 0; i < 3; ++i )
    {
      EXPECT_DOUBLE_EQ( y.get( inds[i] ), vals[i] );
      EXPECT_DOUBLE_EQ( z.get( inds[i] ), vals[i] + 1.0 );
    }
  }

  // Testing dot, axpy, axpby
  x.set( 1.0 );
  y.set( 2.0 );
  z.set( 3.0 );

  real64 const dotprod = x.dot( y );
  EXPECT_EQ( dotprod, 2 * y.globalSize() ); // sum_size 2

  y.axpy( 2.0, x );
  for( globalIndex i = y.ilower(); i < y.iupper(); ++i )
  {
    EXPECT_EQ( y.get( i ), 4.0 ); // 2*1 + 2
  }

  z.axpby( 2.0, x, 3.0 );
  for( globalIndex i = z.ilower(); i < z.iupper(); ++i )
  {
    EXPECT_EQ( z.get( i ), 11.0 ); // 2*1 + 3*3
  }

  // Testing norms
  x.zero();
  if( rank == 0 )
  {
    globalIndex const inds2[2] = { 0, 1 };
    real64 const vals2[2] = { 3.0, -4.0 };
    x.set( inds2, vals2, 2 ); // 3, -4, 0
  }
  x.close();
  EXPECT_EQ( x.norm1(), 7.0 );
  EXPECT_EQ( x.norm2(), 5.0 );
  EXPECT_EQ( x.normInf(), 4.0 );

  // Testing extractLocalVector
  real64 * localVec;
  x.extractLocalVector( &localVec );
  for( globalIndex i = x.ilower(); i < x.iupper(); ++i )
  {
    EXPECT_EQ( localVec[x.getLocalRowID( i )], x.get( i ) );
  }
}
/**
 * @function testMatrixFunctions
 *
 * @brief Test matrix functions including create, add/set, accessors,
 * and linear algebra operations.
 */
// -----------------------------------------
// Test matrix functions
// -----------------------------------------
template<typename LAI>
void testMatrixFunctions()
{
  // Define aliases
  using ParallelVector = typename LAI::ParallelVector;
  using ParallelMatrix = typename LAI::ParallelMatrix;

  // Get the MPI rank
  int rank;
  int numranks;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  MPI_Comm_size( MPI_COMM_WORLD, &numranks );

  // Define some vectors, matrices
  ParallelVector vec1, vec2, vec3;
  ParallelMatrix mat1, mat2, mat3, mat4;

  mat1.createWithLocalSize( 2, 2, MPI_COMM_WORLD ); // 2*numranks x 2*numranks 
  mat2.createWithGlobalSize( 2, 2, MPI_COMM_WORLD ); // 2x2
  mat3.createWithLocalSize( 2, 3, 3, MPI_COMM_WORLD ); // 2*numranks x 3*numranks 
  mat4.createWithGlobalSize( 3, 4, 3, MPI_COMM_WORLD ); // 3x4 

  // Testing create, globalRows, globalCols
  localIndex rows1 = mat1.globalRows();
  localIndex cols1 = mat1.globalCols();
  localIndex rows2 = mat2.globalRows();
  localIndex cols2 = mat2.globalCols();
  localIndex rows3 = mat3.globalRows();
  localIndex cols3 = mat3.globalCols();
  localIndex rows4 = mat4.globalRows();
  localIndex cols4 = mat4.globalCols();
  EXPECT_EQ( rows1, 2*numranks );
  EXPECT_EQ( cols1, 2*numranks );
  EXPECT_EQ( rows2, 2 );
  EXPECT_EQ( cols2, 2 );
  EXPECT_EQ( rows3, 2*numranks );
  EXPECT_EQ( cols3, 3*numranks );
  EXPECT_EQ( rows4, 3 );
  EXPECT_EQ( cols4, 4 );

  // Testing add/set/insert element
  mat1.insert( 1, 0, .5 );
  mat1.close();
  mat1.set( 1, 0, 5 );
  mat1.close();
  mat1.add( 1, 0, 1 );
  mat1.add( 1, 0, 2 );
  mat1.close();

  // Testing add/set/insert c-style, getRowCopy
  globalIndex inds1[2] = {0, 2};
  globalIndex inds2[1] = {0};
  globalIndex inds3[3] = {0, 1, 2};
  real64 vals1[2] = {5, 10};
  real64 vals2[1] = {1};
  real64 vals3[3] = {.5, 1, 2};
  mat4.insert( 1, inds3, vals3, 3 );
  mat4.close();
  mat4.set( 1, inds1, vals1, 2 );
  mat4.close();
  mat4.add( 1, inds2, vals2, 1 );
  mat4.close();

  array1d<real64> colvals;
  array1d<globalIndex> colinds;
  mat4.getRowCopy( 1, colinds, colvals );
  EXPECT_EQ( colinds.size(), 3 );
  EXPECT_DOUBLE_EQ( colvals[0], 6 );
  EXPECT_DOUBLE_EQ( colvals[1], 1 );
  EXPECT_DOUBLE_EQ( colvals[2], 10 );

  // Testing add/set/insert array1d
  ParallelMatrix mat6;
  mat6.createWithGlobalSize( 4, 4, MPI_COMM_WORLD );
  array1d<real64> vals6( 3 );
  array1d<real64> vals7( 3 );
  array1d<globalIndex> inds6( 3 );
  vals6[0] = 1;
  vals6[1] = .5;
  vals6[2] = -3;
  vals7[0] = 1;
  vals7[1] = 1;
  vals7[2] = 1;
  inds6[0] = 0;
  inds6[1] = 1;
  inds6[2] = 3; 
  mat6.insert(0, inds6, vals6);
  mat6.close();
  mat6.set(0, inds6, vals7);
  mat6.close();
  mat6.add(0, inds6, vals6);
  mat6.close();

  // Testing add/set/insert array2d
  ParallelMatrix mat7;
  mat7.createWithGlobalSize( 4, 4, MPI_COMM_WORLD );
  array1d<globalIndex> rows( 2 );
  array1d<globalIndex> cols( 2 );
  array2d<real64> vals8( 2, 2 );
  rows[0] = 0;
  rows[1] = 2;
  cols[0] = 1;
  cols[1] = 3;
  vals8[0][0] = 1;
  vals8[0][1] = 2;
  vals8[1][0] = 3;
  vals8[1][1] = 4;
  mat7.insert( rows, cols, vals8 );
  mat7.close();
  mat7.add( rows, cols, vals8 );
  mat7.close();

  // Testing set and zero
  mat7.set( 2 );
  mat7.close();
  mat1.zero();
  mat1.close();

  // Testing vector multiply, matrix multiply, MatrixMatrixMultiply
  vec1.createWithGlobalSize( 2, MPI_COMM_WORLD );
  vec2.createWithGlobalSize( 2, MPI_COMM_WORLD );
  vec1.set( 1 );
  vec1.close();
  globalIndex inds4[2] = {0, 1};
  real64 vals4[2] = {1, 3};
  real64 vals5[2] = {2, 1};
  mat2.insert( 0, inds4, vals4, 2 );
  mat2.insert( 1, inds4, vals5, 2 );
  mat2.close();
  mat2.multiply(vec1, vec2);
  EXPECT_DOUBLE_EQ( vec2.get(0), 4 );
  EXPECT_DOUBLE_EQ( vec2.get(1), 3 );
  mat2.multiply(mat2, mat1);
  array1d<real64> colvals2;
  array1d<globalIndex> colinds2;
  mat1.getRowCopy( 0, colinds2, colvals2 );
  EXPECT_DOUBLE_EQ( colvals2[0], 7 );
  EXPECT_DOUBLE_EQ( colvals2[1], 6 );
  ParallelMatrix mat8;
  mat8.createWithGlobalSize( 2, 2, MPI_COMM_WORLD );
  mat1.MatrixMatrixMultiply( false, mat2, false, mat8, false );
  mat1.MatrixMatrixMultiply( true, mat2, false, mat8, false );
  mat1.MatrixMatrixMultiply( true, mat2, true, mat8, false );

  // Testing residual, gemv
  vec3.createWithGlobalSize( 2, MPI_COMM_WORLD );
  mat2.residual(vec1, vec2, vec3); 
  EXPECT_DOUBLE_EQ( vec3.get(0), 0 ); // mat2*vec1 = vec2
  EXPECT_DOUBLE_EQ( vec3.get(1), 0 );
  vec3.copy(vec2);
  mat2.gemv( 2, vec1, .5, vec2 );
  EXPECT_DOUBLE_EQ( vec2.get(0), 10 );
  EXPECT_DOUBLE_EQ( vec2.get(1), 7.5 );
  mat2.gemv( 2, vec1, .5, vec3, true );
  EXPECT_DOUBLE_EQ( vec3.get(0), 8 );
  EXPECT_DOUBLE_EQ( vec3.get(1), 9.5 );

  // Testing scale, left/right scaling
  mat2.scale( 2 );
  array1d<real64> colvals4;
  array1d<globalIndex> colinds4;
  mat2.getRowCopy( 0, colinds4, colvals4 );
  EXPECT_DOUBLE_EQ( colvals4[0], 2 );
  EXPECT_DOUBLE_EQ( colvals4[1], 6 );
  vec1.set( 0, 2 );
  vec1.set( 1, 0.5 );
  vec1.close();
  vec2.set( 0, .5 );
  vec2.set( 1, 3 );
  vec2.close();
  mat2.leftScale( vec1 );
  mat2.rightScale( vec2 );
  mat2.leftRightScale( vec2, vec1 );
  array1d<real64> colvals5;
  array1d<globalIndex> colinds5;
  mat2.getRowCopy( 1, colinds5, colvals5 );
  EXPECT_DOUBLE_EQ( colvals5[0], 6 );
  EXPECT_DOUBLE_EQ( colvals5[1], 4.5 );

  // Testing clearRow
  ParallelMatrix mat5 = compute2DLaplaceOperator<LAI>( MPI_COMM_WORLD, 2 );
  mat5.clearRow(0, 5);
  array1d<real64> colvals3;
  array1d<globalIndex> colinds3;
  mat5.getRowCopy( 0, colinds3, colvals3 );
  EXPECT_DOUBLE_EQ( colvals3[0], 5 );

  // Testing print, write, printParallelMatrix
  // mat5.print(); 
  // mat5.write("matout.mtx", true);
  // mat5.write("matout.m", false);
  // mat5.printParallelMatrix("matout2.mtx");

  // Testing ilower, iupper, getLocalRowID, getGlobalRowID, numMyCols
  if (numranks == 1) {
    EXPECT_EQ( mat5.ilower(), 0 );
    EXPECT_EQ( mat5.iupper(), 4 ); 
  }
  EXPECT_EQ( mat5.globalRows(), 4 );
  EXPECT_EQ( mat5.globalCols(), 4 ); 
  EXPECT_EQ( mat5.getLocalRowID( 0 ), 0 );
  EXPECT_EQ( mat5.getLocalRowID( 2 ), 2 );
  EXPECT_EQ( mat5.getGlobalRowID( 1 ), 1 );
  EXPECT_EQ( mat5.getGlobalRowID( 3 ), 3 );  
  // EXPECT_EQ( mat5.numMyCols(), 4 ); 

}

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
template<typename LAI>
void testInterfaceSolvers()
{
  // Define aliases templated on the Linear Algebra Interface (LAI).
  using ParallelMatrix = typename LAI::ParallelMatrix;
  using ParallelVector = typename LAI::ParallelVector;
  using LinearSolver = typename LAI::LinearSolver;

  // Get the MPI rank
  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  // Use an nxn cartesian mesh to generate the Laplace 2D operator.
  globalIndex n = 100;
  globalIndex N = n * n;

  // Compute a 2D Laplace operator
  ParallelMatrix matrix = compute2DLaplaceOperator<LAI>( MPI_COMM_WORLD, n );

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
  b.set( 1.0 ); // ones

  // Also define a residual vector, this time using the copy constructor
  ParallelVector r( b );

  // Test dot product: r.b = b.b = N
  real64 dotTest = r.dot( b );
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

  // Now create a solver parameter list and solver
  LinearSolverParameters parameters;
  LinearSolver solver( parameters );

  // Set basic options
  parameters.verbosity = 0;
  parameters.solverType = "cg";

  parameters.krylov.tolerance = 1e-8;
  parameters.krylov.maxIterations = 250;
  parameters.preconditionerType = "amg";
  parameters.amg.smootherType = "gaussSeidel";
  parameters.amg.coarseType = "direct";

  // Solve using the iterative solver and compare norms with true solution
  solver.solve( matrix, x_comp, b );
  real64 norm_comp = x_comp.norm2();
  real64 norm_true = x_true.norm2();
  EXPECT_LT( std::fabs( norm_comp / norm_true - 1. ), 1e-6 );

  // We now do the same using a direct solver.
  // Again the norm should be the norm of x. We use a tougher tolerance on the test
  // compared to the iterative solution. This should be accurate to machine precision
  // and some round off error. We (arbitrarily) chose 1e-12 as a good guess.
  x_comp.zero();
  parameters.solverType = "direct";
  solver.solve( matrix, x_comp, b );
  norm_comp = x_comp.norm2();
  EXPECT_LT( std::fabs( norm_comp / norm_true - 1. ), 1e-12 );

  // Option to write files (for direct comparison)
  // matrix.write("matrix.dat");
  // x_true.write("x_true.dat");
  // x_comp.write("x_comp.dat");

  // Try getting access to matrix entries

  array1d<real64> col_values;
  array1d<globalIndex> col_indices;

  if( rank == 0 )
  {
    matrix.getRowCopy( 0, col_indices, col_values );
    EXPECT_EQ( col_indices.size(), 3 );
    matrix.getRowCopy( 1, col_indices, col_values );
    EXPECT_EQ( col_indices.size(), 4 );
    matrix.getRowCopy( n + 1, col_indices, col_values );
    EXPECT_EQ( col_indices.size(), 5 );
  }
  // Try clearing rows and setting diagonal value

  double diagValue = 100.0;
  globalIndex firstRow = matrix.ilower();

  matrix.clearRow( firstRow, diagValue );
  matrix.close();

  matrix.getRowCopy( firstRow, col_indices, col_values );
  for( localIndex i = 0 ; i < col_indices.size() ; ++i )
  {
    if( firstRow == col_indices[i] )
      EXPECT_DOUBLE_EQ( col_values[i], diagValue );
    else
      EXPECT_DOUBLE_EQ( col_values[i], 0.0 );
  }
  EXPECT_DOUBLE_EQ( matrix.getDiagValue( firstRow ), diagValue );
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
template<typename LAI>
void testGEOSXSolvers()
{
  // Define aliases templated on the Linear Algebra Interface (LAI).
  using ParallelMatrix = typename LAI::ParallelMatrix;
  using ParallelVector = typename LAI::ParallelVector;

  // Get the MPI rank
  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  // Use nxn cartesian mesh to generate the Laplace 2D operator.
  globalIndex n = 100;
  globalIndex N = n * n;

  // Compute a 2D Laplace operator and identity matrix
  ParallelMatrix matrix = compute2DLaplaceOperator<LAI>( MPI_COMM_WORLD, n );
  ParallelMatrix identity = computeIdentity<LAI>( MPI_COMM_WORLD, N );

  // Define vectors 
  ParallelVector x_true,
      x_comp,
      b;

  x_true.createWithGlobalSize( N, MPI_COMM_WORLD );
  x_comp.createWithGlobalSize( N, MPI_COMM_WORLD );
  b.createWithGlobalSize( N, MPI_COMM_WORLD );

  x_true.rand();
  x_comp.zero();

  matrix.multiply( x_true, b );

  // Test with CG solver
  CGsolver<LAI> testCG;
  testCG.solve( matrix, x_comp, b, identity );

  real64 norm_true = x_true.norm2();
  real64 norm_comp = x_comp.norm2();
  EXPECT_LT( std::fabs( norm_comp / norm_true - 1. ), 5e-6 );

  // Test with BiCGSTAB solver 
  x_comp.zero();

  BiCGSTABsolver<LAI> testBiCGSTAB;
  testBiCGSTAB.solve( matrix, x_comp, b, identity );

  norm_true = x_true.norm2();
  norm_comp = x_comp.norm2();
  EXPECT_LT( std::fabs( norm_comp / norm_true - 1. ), 5e-6 );

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
template<typename LAI>
void testGEOSXBlockSolvers()
{
  // The usual typenames
  using ParallelMatrix = typename LAI::ParallelMatrix;
  using ParallelVector = typename LAI::ParallelVector;

  // Get the MPI rank
  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  // We are going to assembly the following dummy system
  // [L L] [x_true] = [b_0]
  // [L L] [x_true] = [b_1]

  globalIndex n = 100;
  globalIndex N = n * n;

  ParallelMatrix matrix = compute2DLaplaceOperator<LAI>( MPI_COMM_WORLD, n );
  ParallelMatrix identity = computeIdentity<LAI>( MPI_COMM_WORLD, N );

  ParallelVector x_true,
      x_comp_0,
      x_comp_1,
      b_0,
      b_1;

  x_true.createWithGlobalSize( N, MPI_COMM_WORLD );
  x_comp_0.createWithGlobalSize( N, MPI_COMM_WORLD );
  x_comp_1.createWithGlobalSize( N, MPI_COMM_WORLD );
  b_0.createWithGlobalSize( N, MPI_COMM_WORLD );
  b_1.createWithGlobalSize( N, MPI_COMM_WORLD );

  x_true.rand();
  x_comp_0.zero();
  x_comp_1.zero();

  // Size of the block system
  localIndex nRows = 2;
  localIndex nCols = 2;

  // Declare and allocate block matrices/vectors
  // System matrix

  

  BlockMatrixView<LAI> block_matrix( nRows, nCols );
  BlockMatrixView<LAI> block_precon( nRows, nCols );
  BlockVectorView<LAI> block_x_true( nCols );
  BlockVectorView<LAI> block_x_comp( nCols );
  BlockVectorView<LAI> block_rhs( nRows );

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
  block_matrix.multiply( block_x_true, block_rhs );

//TODO: Need to refactor Native block solvers.  Disable this testing section for now.
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

  x_comp_0.zero();// TODO: fix block zero()
  x_comp_1.zero();

  BiCGSTABsolver<LAI> testBiCGSTAB;
  testBiCGSTAB.solve( block_matrix, block_x_comp, block_rhs, block_precon );

  norm_x_true = block_x_true.norm2();
  norm_x_comp = block_x_comp.norm2();
  EXPECT_LT( std::fabs( norm_x_comp/norm_x_true - 1. ), 5e-6 );
#endif
}

//------------------------------
// Test matrix-matrix operations
//------------------------------
// Currently just test matrix-matrix multiply, but eventually
// should include add and other level-III operations
template<typename LAI>
void testMatrixMatrixOperations()
{
  using ParallelMatrix = typename LAI::ParallelMatrix;

  globalIndex n = 100;
  globalIndex N = n * n;

  ParallelMatrix A = compute2DLaplaceOperator<LAI>( MPI_COMM_WORLD, n );

  ParallelMatrix A_squared;
  A_squared.createWithGlobalSize( N, 1, MPI_COMM_WORLD );

  A.multiply( A, A_squared );

  real64 a = A.normInf();
  real64 b = A_squared.normInf();

  EXPECT_DOUBLE_EQ( a * a, b );
}

//-----------------------------------
// Test rectangular matrix operations
//-----------------------------------
// Create a rectangular matrix and check its sizes and norms
template<typename LAI>
void testRectangularMatrixOperations()
{
  using ParallelMatrix = typename LAI::ParallelMatrix;

  int mpiSize;
  MPI_Comm_size( MPI_COMM_WORLD, &mpiSize );

  // Set a size that allows to run with arbitrary number of processes
  globalIndex const nRows = std::max( 100, mpiSize );
  globalIndex const nCols = 2 * nRows;

  ParallelMatrix A;
  A.createWithGlobalSize( nRows, nCols, 2, MPI_COMM_WORLD );

  for( globalIndex i = A.ilower() ; i < A.iupper() ; ++i )
  {
    real64 const entry = static_cast<real64>( i + 1 );
    A.insert( i, 2 * i, entry );
    A.insert( i, 2 * i + 1, -entry );
  }
  A.close();

  // Check on sizes
  EXPECT_EQ( A.globalRows(), nRows );
  EXPECT_EQ( A.globalCols(), nCols );

  // Check on norms
  real64 const a = A.norm1();
  real64 const b = A.normInf();
  real64 const c = A.normFrobenius();

  EXPECT_DOUBLE_EQ( a, static_cast<real64>( nRows ) );
  EXPECT_DOUBLE_EQ( b, static_cast<real64>( nCols ) );
  EXPECT_DOUBLE_EQ( c, sqrt( static_cast<real64>( nRows * ( nRows + 1 ) * ( 2 * nRows + 1 ) ) / 3 ) );
}

// END_RST_NARRATIVE

//@}

/*! @name Ctest tests.
 * @brief Runs similar testing functions using different Linear Algebra Interfaces (LAIs).
 */
//@{

/*! @function testEpetraLAOperations.
 * @brief Runs all tests using the Trilinos interface.
 */
TEST( testLAOperations, testEpetraLAOperations )
{
  testVectorFunctions< TrilinosInterface >();
  // testMatrixFunctions<TrilinosInterface>();
  testInterfaceSolvers< TrilinosInterface >();
  testGEOSXSolvers< TrilinosInterface >();
  testGEOSXBlockSolvers< TrilinosInterface >();
  testMatrixMatrixOperations< TrilinosInterface >();
  testRectangularMatrixOperations< TrilinosInterface >();
}

/*! @function testHypreLAOperations.
 * @brief Runs all tests using the Hypre interface.
 */
//TEST(testLAOperations,testHypreLAOperations)
//{
  //testVectorFunctions<HypreInterface>();
  //testMatrixFunctions<HypreInterface>();
  //testInterfaceSolvers<HypreInterface>();
  //testGEOSXSolvers<HypreInterface>();
  //testGEOSXBlockSolvers<HypreInterface>();
//}

/*! @function testPETScLAOperations.
 * @brief Runs all tests using the PETSc interface.
 */
TEST( testLAOperations, testPETScLAOperations )
{
  PetscInterface();
  PetscInitializeNoArguments(); // set up PETSc directly in tests
  testVectorFunctions< PetscInterface >();
  // testMatrixFunctions<PetscInterface>();
  testInterfaceSolvers< PetscInterface >();
  testGEOSXSolvers< PetscInterface >();
  // testGEOSXBlockSolvers<PetscInterface>(); // doesn't build
  testMatrixMatrixOperations< PetscInterface >();
  testRectangularMatrixOperations< PetscInterface >();
  PetscFinalize();
}
//@}

int main( int argc, char ** argv )
 {
   ::testing::InitGoogleTest( &argc, argv );
   setupMPI( argc, argv );
   setupOpenMP();
   setupMKL();
   logger::InitializeLogger( MPI_COMM_GEOSX );

   int const result = RUN_ALL_TESTS();

   geosx::basicCleanup();
   return result;
 }
