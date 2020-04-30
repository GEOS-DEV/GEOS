/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file testLAOperations.cpp
 */

#include <gtest/gtest.h>
#include "codingUtilities/UnitTestUtilities.hpp"
#include "testLinearAlgebraUtils.hpp"
#include "common/DataTypes.hpp"
#include "managers/initialization.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "linearAlgebra/utilities/LinearSolverParameters.hpp"

/**
 * \file testLAOperations.cpp
 * \brief This test file is part of the ctest suite and tests linear algebra interfaces.
 * It mainly uses dummy 2D Laplace operator matrices
 * to test the solvers along with a simple (block) identity preconditioner.
 */

using namespace geosx;

static real64 const machinePrecision = 20.0 * std::numeric_limits< real64 >::epsilon();

// BEGIN_RST_NARRATIVE testLAOperations.rst

/**
 * @brief Test case class for LA operations tests
 * @tparam LAI Linear Algebra Interface
 */
template< typename LAI >
class LAOperationsTest : public ::testing::Test
{};

TYPED_TEST_SUITE_P( LAOperationsTest );

/*! @name Test functions.
 * @brief Templated functions to test the linear solvers.
 */
//@{

// ==============================
// Test Linear Algebra Operations
// ==============================
// In these functions we test the linear algebra (vector and matrix) operations

/**
 * @function testVectorFunction
 *
 * @brief Test vector functions including create, add/set, accessors,
 * and linear algebra operations.
 */
// -----------------------------------------
// Test vector functions
// -----------------------------------------
TYPED_TEST_P( LAOperationsTest, VectorFunctions )
{
  // Define alias
  using Vector = typename TypeParam::ParallelVector;

  // Get the MPI rank
  int const rank = MpiWrapper::Comm_rank( MPI_COMM_GEOSX );
  int const numRanks = MpiWrapper::Comm_size( MPI_COMM_GEOSX );

  Vector x;
  localIndex const localSize = 3;
  globalIndex const globalSize = localSize * numRanks;
  globalIndex const offset = rank * localSize;

  // Testing createWithLocalSize
  x.createWithLocalSize( localSize, MPI_COMM_GEOSX );
  EXPECT_EQ( x.localSize(), localSize );
  EXPECT_EQ( x.globalSize(), globalSize );

  // Testing iupper/ilower
  EXPECT_EQ( x.ilower(), offset );
  EXPECT_EQ( x.iupper(), offset + localSize );

  // Testing setting/getting values locally
  x.open();
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
  x.createWithGlobalSize( globalSize, MPI_COMM_GEOSX );
  EXPECT_EQ( x.localSize(), localSize );
  EXPECT_EQ( x.globalSize(), globalSize );

  // Testing adding global values on rank 0 and getting locally
  x.open();
  if( rank == 0 )
  {
    for( globalIndex i = 0; i < x.globalSize(); ++i )
    {
      x.add( i, 2 * i );
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

  Vector v;
  v.create( localVals, MPI_COMM_GEOSX );
  for( globalIndex i = v.ilower(); i < v.iupper(); ++i )
  {
    EXPECT_EQ( v.get( i ), localVals[v.getLocalRowID( i )] );
  }

  // Testing copy constructor, assignment, get element
  Vector y( x );
  Vector z;
  z = x;
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
  x.open();
  x.set( offset, -1 );
  x.close(); // set/add can't be interchanged
  x.open();
  x.add( offset + 1, 10 );
  x.close();
  EXPECT_DOUBLE_EQ( x.get( offset ), -1 );
  EXPECT_DOUBLE_EQ( x.get( offset + 1 ), offset * 2 + 12 );

  // Testing add/set c-style
  {
    globalIndex const inds[3] = { offset, offset + 1, offset + 2 };
    real64 const vals[3] = { -5.0, -6.0, 0.0 };
    y.zero();
    y.open();
    y.set( inds, vals, 2 );
    y.close();
    z.set( 1.0 );
    z.open();
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
    inds[0] = offset;
    inds[1] = offset + 1;
    inds[2] = offset + 2;
    array1d< real64 > vals( 3 );
    vals[0] = -5.0;
    vals[1] = -6.0;
    vals[2] = 0.0;
    y.zero();
    y.open();
    y.set( inds, vals );
    y.close();
    z.set( 1.0 );
    z.open();
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
  x.open();
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
  real64 const * localVec = x.extractLocalVector();
  for( globalIndex i = x.ilower(); i < x.iupper(); ++i )
  {
    EXPECT_EQ( localVec[x.getLocalRowID( i )], x.get( i ) );
  }
}

#if 0
/**
 * @function testMatrixFunctions
 *
 * @brief Test matrix functions including create, add/set, accessors,
 * and linear algebra operations.
 */
// -----------------------------------------
// Test matrix functions
// -----------------------------------------
TYPED_TEST_P( LAOperationsTest, MatrixFunctions )
{
  // Define aliases
  using Vector = typename TypeParam::ParallelVector;
  using Matrix = typename TypeParam::ParallelMatrix;

  // Get the MPI rank
  int numranks = MpiWrapper::Comm_size( MPI_COMM_GEOSX );
  int rank = MpiWrapper::Comm_rank( MPI_COMM_GEOSX );

  std::cout << "*** Rank: " << rank << std::endl;

  // Dummy vector and Matrix
  Matrix C;
  Matrix D;
  {
    // Test matrix-matrix product: C = A*B
    Matrix A;
    compute2DLaplaceOperator( MPI_COMM_GEOSX, 2 * numranks, A );
    Matrix B( A );

    A.multiply( B, C );
    A.leftMultiplyTranspose( A, D );
  }

  // Define some vectors, matrices
  Vector vec1, vec2, vec3;
  Matrix mat1, mat2, mat3, mat4;
  mat1.createWithLocalSize( 2, 2, MPI_COMM_GEOSX ); // 2*numranks x 2*numranks
  mat2.createWithGlobalSize( 2, 2, MPI_COMM_GEOSX ); // 2x2
  mat3.createWithLocalSize( 2, 3, 3, MPI_COMM_GEOSX ); // 2*numranks x 3*numranks
  mat4.createWithGlobalSize( 3, 4, 3, MPI_COMM_GEOSX ); // 3x4

  // Testing create, globalRows, globalCols
  EXPECT_EQ( mat1.numGlobalRows(), 2 * numranks );
  EXPECT_EQ( mat1.numGlobalCols(), 2 * numranks );
  EXPECT_EQ( mat2.numGlobalRows(), 2 );
  EXPECT_EQ( mat2.numGlobalCols(), 2 );
  EXPECT_EQ( mat3.numGlobalRows(), 2 * numranks );
  EXPECT_EQ( mat3.numGlobalCols(), 3 * numranks );
  EXPECT_EQ( mat4.numGlobalRows(), 3 );
  EXPECT_EQ( mat4.numGlobalCols(), 4 );

  // Testing add/set/insert element
  //  mat1.insert( 1, 0, .5 );
  //  mat1.close();
  //  mat1.set( 1, 0, 5 );
  //  mat1.close();
  mat1.open();
  mat1.add( 1, 0, 1 );
  mat1.add( 1, 0, 2 );
  mat1.close();

  // Testing add/set/insert c-style, getRowCopy
  globalIndex inds1[2] = { 0, 2 };
  globalIndex inds2[1] = { 0 };
  globalIndex inds3[3] = { 0, 1, 2 };
  real64 vals1[2] = { 5, 10 };
  real64 vals2[1] = { 1 };
  real64 vals3[3] = { .5, 1, 2 };

  globalIndex iRow = 1;

  if( ( mat4.ilower() <= iRow ) && ( iRow < mat4.iupper() ) )
  {
    mat4.insert( iRow, inds3, vals3, 3 );
    //    mat4.close();
    //    mat4.open();
    mat4.set( iRow, inds1, vals1, 2 );
    //    mat4.close();
    //    mat4.open();
    mat4.add( iRow, inds2, vals2, 1 );
    //    mat4.close();
  }
  mat4.close();

  array1d< real64 > colvals_CHECK( 3 );
  colvals_CHECK( 0 ) = 6;
  colvals_CHECK( 1 ) = 1;
  colvals_CHECK( 2 ) = 10;

  if( ( mat4.ilower() <= iRow ) && ( iRow < mat4.iupper() ) )
  {
    localIndex const rowLength = mat.getGlobalRowLength( iRow );
    EXPECT_EQ( rowLength, colvals_CHECK.size() );
    array1d< real64 > colvals( rowLength );
    array1d< globalIndex > colinds( rowLength );
    mat4.getRowCopy( iRow, colinds, colvals );
    for( int i = 0; i < 3; ++i )
    {
      EXPECT_DOUBLE_EQ( colvals( colinds[i] ), colvals_CHECK( i ) ); //HYPRE does not return sorted cols!
    }
  }
  // Testing add/set/insert array1d
  Matrix mat6;
  mat6.createWithGlobalSize( 4, 4, MPI_COMM_GEOSX );
  array1d< real64 > vals6( 3 );
  array1d< real64 > vals7( 3 );
  array1d< globalIndex > inds6( 3 );
  vals6[0] = 1;
  vals6[1] = .5;
  vals6[2] = -3;
  vals7[0] = 1;
  vals7[1] = 1;
  vals7[2] = 1;
  inds6[0] = 0;
  inds6[1] = 1;
  inds6[2] = 3;

  iRow = 0;
  if( ( mat6.ilower() <= iRow ) && ( iRow < mat6.iupper() ) )
  {
    mat6.insert( iRow, inds6, vals6 );
    //	  mat6.close();
    mat6.set( iRow, inds6, vals7 );
    //	  mat6.close();
    mat6.add( iRow, inds6, vals6 );
  }
  mat6.close();

  // Testing add/set/insert array2d
  Matrix mat7;
  mat7.createWithGlobalSize( 4, 4, MPI_COMM_GEOSX );
  array1d< globalIndex > rows( 2 );
  array1d< globalIndex > cols( 2 );
  array2d< real64 > vals8( 2, 2 );
  rows[0] = 0;
  rows[1] = 2;
  cols[0] = 1;
  cols[1] = 3;
  vals8[0][0] = 1;
  vals8[0][1] = 2;
  vals8[1][0] = 3;
  vals8[1][1] = 4;
  if( ( mat7.ilower() <= *std::min_element( rows.data(), rows.data() + rows.size() ) ) &&
      ( *std::max_element( rows.data(), rows.data() + rows.size() ) < mat7.iupper() ) )
  {
    mat7.insert( rows, cols, vals8 );
    //    mat7.close();
    mat7.add( rows, cols, vals8 );
    //    mat7.close();
  }
  mat7.close();

  // Testing set and zero
  mat7.open();
  mat7.set( 2 );
  mat7.close();

  mat7.open();
  mat1.zero();
  mat1.close();

  // Testing vector multiply, matrix multiply, MatrixMatrixMultiply
  vec1.createWithGlobalSize( 2, MPI_COMM_GEOSX );
  vec2.createWithGlobalSize( 2, MPI_COMM_GEOSX );
  vec1.set( 1 );
  vec1.close();
  globalIndex inds4[2] = { 0, 1 };
  real64 vals4[2] = { 1, 3 };
  real64 vals5[2] = { 2, 1 };

  if( ( mat2.ilower() <= 0 ) && ( 0 < mat2.iupper() ) )
  {
    mat2.insert( 0, inds4, vals4, 2 );
  }
  if( ( mat2.ilower() <= 1 ) && ( 1 < mat2.iupper() ) )
  {
    mat2.insert( 1, inds4, vals5, 2 );
  }
  mat2.close();

  mat2.multiply( vec1, vec2 );

  if( ( vec2.ilower() <= 0 ) && ( 0 < vec2.iupper() ) )
  {
    EXPECT_DOUBLE_EQ( vec2.get( 0 ), 4 );
  }
  if( ( vec2.ilower() <= 1 ) && ( 1 < vec2.iupper() ) )
  {
    EXPECT_DOUBLE_EQ( vec2.get( 1 ), 3 );
  }

  // Matrix-Matrix multiply
  Matrix mat2mat2;

  {
    Matrix mat22( mat2 );

    mat22.multiply( mat22, mat2mat2 );
  }
}
#endif

// -------------------------------------
// Test libraries operations and solvers
// -------------------------------------
// We start by testing the linear algebra operations. We fill two matrices (one will be a
// preconditioner) and make sure the sparse storage is behaving properly. We then test the
// iterative and direct solvers available.

/**
 * @function testInterfaceSolvers
 *
 * @brief Test the packaged solvers from the LAI as well as basic linear algebra operations,
 * such as matrix-vector products, dot products, norms and residuals.
 */
TYPED_TEST_P( LAOperationsTest, InterfaceSolvers )
{
#ifdef GEOSX_USE_PETSC
  SKIP_TEST_IF( (std::is_same< TypeParam, PetscInterface >::value), "https://github.com/GEOSX/GEOSX/issues/790" );
#endif

  // Define aliases templated on the Linear Algebra Interface (LAI).
  using Matrix = typename TypeParam::ParallelMatrix;
  using Vector = typename TypeParam::ParallelVector;
  using Solver = typename TypeParam::LinearSolver;

  // Use an nxn cartesian mesh to generate the Laplace 2D operator.
  globalIndex n = 100;
  globalIndex N = n * n;

  // Compute a 2D Laplace operator
  Matrix matrix;
  compute2DLaplaceOperator( MPI_COMM_GEOSX, n, matrix );

  // Define some vectors
  Vector x_true;
  Vector x_comp;
  Vector b;

  x_true.createWithGlobalSize( N, MPI_COMM_GEOSX );
  x_comp.createWithGlobalSize( N, MPI_COMM_GEOSX );
  b.createWithGlobalSize( N, MPI_COMM_GEOSX );

  // We have some simple initialization options for vectors:
  x_true.rand(); // random
  x_comp.zero(); // zero
  b.set( 1.0 ); // ones

  // Also define a residual vector, this time using the copy constructor
  Vector r( b );

  // Test dot product: r.b = b.b = N
  real64 dotTest = r.dot( b );
  EXPECT_DOUBLE_EQ( dotTest, N );

  // Test various norms
  real64 norm1 = b.norm1();
  real64 norm2 = b.norm2();
  real64 normInf = b.normInf();

  EXPECT_NEAR( norm1, N, N * machinePrecision );
  EXPECT_NEAR( norm2, n, n * machinePrecision );
  EXPECT_NEAR( normInf, 1., machinePrecision );

  // Compute the matrix/vector multiplication. We compute b as Ax and will aim to get x
  // back from the solvers.
  matrix.apply( x_true, b );

  // Test the residual function by computing r = b - Ax = 0
  matrix.residual( x_true, b, r );
  real64 normRes = r.normInf();
  EXPECT_NEAR( normRes, 0., machinePrecision );

  // Now create a solver parameter list and solver
  LinearSolverParameters parameters;
  Solver solver( parameters );

  // Set basic options
  parameters.logLevel = 0;
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
}


//------------------------------
// Test matrix-matrix operations
//------------------------------
// Currently just test matrix-matrix multiply, but eventually
// should include add and other level-III operations
TYPED_TEST_P( LAOperationsTest, MatrixMatrixOperations )
{
  using Matrix = typename TypeParam::ParallelMatrix;

  globalIndex const n = 100;
  Matrix A;
  compute2DLaplaceOperator( MPI_COMM_GEOSX, n, A );

  Matrix A_squared;
  A.multiply( A, A_squared );

  real64 const a = A.normInf();
  real64 const b = A_squared.normInf();

  EXPECT_DOUBLE_EQ( a, 8.0 );
  EXPECT_DOUBLE_EQ( b, 64.0 );
}

//-----------------------------------
// Test rectangular matrix operations
//-----------------------------------
// Create a rectangular matrix and check its sizes and norms
TYPED_TEST_P( LAOperationsTest, RectangularMatrixOperations )
{
  using Matrix = typename TypeParam::ParallelMatrix;

  int mpiSize = MpiWrapper::Comm_size( MPI_COMM_GEOSX );

  // Set a size that allows to run with arbitrary number of processes
  globalIndex const nRows = std::max( 100, mpiSize );
  globalIndex const nCols = 2 * nRows;

  Matrix A;
  A.createWithGlobalSize( nRows, nCols, 2, MPI_COMM_GEOSX );

  A.open();
  for( globalIndex i = A.ilower(); i < A.iupper(); ++i )
  {
    real64 const entry = static_cast< real64 >( i + 1 );
    A.insert( i, 2 * i, entry );
    A.insert( i, 2 * i + 1, -entry );
  }
  A.close();

  // Check on sizes
  EXPECT_EQ( A.numGlobalRows(), nRows );
  EXPECT_EQ( A.numGlobalCols(), nCols );

  // Check on norms
  real64 const a = A.norm1();
  real64 const b = A.normInf();
  real64 const c = A.normFrobenius();

  EXPECT_DOUBLE_EQ( a, static_cast< real64 >( nRows ) );
  EXPECT_DOUBLE_EQ( b, static_cast< real64 >( nCols ) );
  EXPECT_DOUBLE_EQ( c, std::sqrt( static_cast< real64 >( nRows * ( nRows + 1 ) * ( 2 * nRows + 1 ) ) / 3.0 ) );
}

// END_RST_NARRATIVE

//@}

REGISTER_TYPED_TEST_SUITE_P( LAOperationsTest,
                             VectorFunctions,
                             //MatrixFunctions,
                             InterfaceSolvers,
                             MatrixMatrixOperations,
                             RectangularMatrixOperations );

#ifdef GEOSX_USE_TRILINOS
INSTANTIATE_TYPED_TEST_SUITE_P( Trilinos, LAOperationsTest, TrilinosInterface, );
#endif

#ifdef GEOSX_USE_HYPRE
INSTANTIATE_TYPED_TEST_SUITE_P( Hypre, LAOperationsTest, HypreInterface, );
#endif

#ifdef GEOSX_USE_PETSC
INSTANTIATE_TYPED_TEST_SUITE_P( Petsc, LAOperationsTest, PetscInterface, );
#endif

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );

  // Avoid setting up signal handlers, due to mysterious ML FPE crashes on Mac
  //geosx::basicSetup( argc, argv );
  setupMPI( argc, argv );
  setupLogger();
  setupOpenMP();
  setupMKL();
  setupLogger();
  setupLAI( argc, argv );

  int const result = RUN_ALL_TESTS();
  geosx::basicCleanup();
  return result;
}
