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
 * @file testKrylovSolvers.cpp
 */

#include <gtest/gtest.h>

#include "testLinearAlgebraUtils.hpp"

#include "common/DataTypes.hpp"
#include "managers/initialization.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "linearAlgebra/utilities/BlockOperatorWrapper.hpp"
#include "linearAlgebra/utilities/BlockOperator.hpp"
#include "linearAlgebra/utilities/BlockVectorWrapper.hpp"
#include "linearAlgebra/solvers/PreconditionerIdentity.hpp"
#include "linearAlgebra/solvers/CGsolver.hpp"
#include "linearAlgebra/solvers/BiCGSTABsolver.hpp"

/**
 * \file testKrylovSolvers.cpp
 * \brief This test file is part of the ctest suite and tests native GEOSX  Krylov solvers.
 * It mainly uses dummy 2D Laplace operator matrices to test the solvers along with a
 * simple (block) identity preconditioner.
 */

using namespace geosx;

// BEGIN_RST_NARRATIVE testKrylovSolvers.rst

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
template< template< typename T > class SOLVER, typename LAI >
void testGEOSXSolvers()
{
  // Define aliases templated on the Linear Algebra Interface (LAI).
  using Matrix = typename LAI::ParallelMatrix;
  using Vector = typename LAI::ParallelVector;

  // Use nxn cartesian mesh to generate the Laplace 2D operator.
  globalIndex constexpr n = 100;
  globalIndex constexpr N = n * n;

  // Compute a 2D Laplace operator and identity matrix
  Matrix matrix;
  Matrix identity;
  compute2DLaplaceOperator( MPI_COMM_GEOSX, n, matrix );
  computeIdentity( MPI_COMM_GEOSX, N, identity );

  // Define vectors
  Vector x_true, x_comp, b;
  x_true.createWithGlobalSize( N, MPI_COMM_GEOSX );
  x_comp.createWithGlobalSize( N, MPI_COMM_GEOSX );
  b.createWithGlobalSize( N, MPI_COMM_GEOSX );

  x_true.rand();
  x_comp.zero();

  matrix.apply( x_true, b );

  // Solve
  SOLVER< Vector > solver( matrix, identity, 1e-8, 300 );
  solver.solve( b, x_comp );
  /////////////////////////////////////////
//  GEOSX_LOG_RANK_0("Iterations: " + std::to_string(solver.numIterations()));
//  if ( MpiWrapper::Comm_rank( MPI_COMM_GEOSX ) == 0)
//  {
//    GEOSX_LOG_RANK_VAR( solver.relativeResidual() );
//  }
  /////////////////////////////////////////

  real64 const norm_true = x_true.norm2();
  real64 const norm_comp = x_comp.norm2();
  EXPECT_LT( std::fabs( norm_comp / norm_true - 1. ), 5e-6 );

  PreconditionerIdentity< LAI > preconIdentity;
  preconIdentity.compute( matrix );
  SOLVER< Vector > solver2( matrix, preconIdentity, 1e-8, 300 );
  x_comp.zero();

  solver2.solve( b, x_comp );
  /////////////////////////////////////////
//  GEOSX_LOG_RANK_0("Iterations: " + std::to_string(solver.numIterations()));
//  if ( MpiWrapper::Comm_rank( MPI_COMM_GEOSX ) == 0)
//  {
//    GEOSX_LOG_RANK_VAR( solver.residualNormVector() );
//  }
  /////////////////////////////////////////


  real64 const norm_true2 = x_true.norm2();
  real64 const norm_comp2 = x_comp.norm2();
  EXPECT_LT( std::fabs( norm_comp2 / norm_true2 - 1. ), 5e-6 );
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
template< template< typename T > class SOLVER, typename LAI >
void testGEOSXBlockSolvers()
{
  // The usual typenames
  using Matrix = typename LAI::ParallelMatrix;
  using Vector = typename LAI::ParallelVector;

  // We are going to assembly the following dummy system
  // [L L] [x_true] = [b_0]
  // [0 L] [x_true] = [b_1]

  globalIndex constexpr n = 100;
  globalIndex constexpr N = n * n;

  // Declare and allocate block matrices/vectors

  // Set up a block operator consisting of two laplacian blocks
  // Here we demonstrate creating a thin wrapper pointing to the same matrix twice
  Matrix matrix;
  compute2DLaplaceOperator( MPI_COMM_GEOSX, n, matrix );
  BlockOperatorWrapper< Vector, Matrix > block_matrix( 2, 2 );
  block_matrix.set( 0, 0, matrix );
  block_matrix.set( 1, 1, matrix );

  // Set up block identity preconditioner
  // Here we instantiate a concrete block matrix with two independent blocks
  BlockOperator< Vector, Matrix > block_precon( 2, 2 );
  computeIdentity( MPI_COMM_GEOSX, N, block_precon.block( 0, 0 ) );
  computeIdentity( MPI_COMM_GEOSX, N, block_precon.block( 1, 1 ) );
  computeZero( MPI_COMM_GEOSX, N, block_precon.block( 0, 1 ) );
  computeZero( MPI_COMM_GEOSX, N, block_precon.block( 1, 0 ) );

  // Create rhs and solution vectors
  BlockVector< Vector > x_true( 2 );
  BlockVector< Vector > x_comp( 2 );
  BlockVector< Vector > b( 2 );

  x_true.block( 0 ).createWithGlobalSize( N, MPI_COMM_GEOSX );
  x_true.block( 1 ).createWithGlobalSize( N, MPI_COMM_GEOSX );
  x_comp.block( 0 ).createWithGlobalSize( N, MPI_COMM_GEOSX );
  x_comp.block( 1 ).createWithGlobalSize( N, MPI_COMM_GEOSX );
  b.block( 0 ).createWithGlobalSize( N, MPI_COMM_GEOSX );
  b.block( 1 ).createWithGlobalSize( N, MPI_COMM_GEOSX );

  x_true.rand();
  x_comp.zero();

  // Set right hand side.
  block_matrix.apply( x_true, b );

  // Create block CG solver object and solve
  SOLVER< BlockVectorView< Vector > > solver( block_matrix, block_precon, 1e-8, 300 );
  solver.solve( b, x_comp );

  // The true solution is the vector x, so we check if the norm are equal.
  // Note: the tolerance is higher that in the previous cases because the matrix is
  // twice as big and the condition number is higher. See details on the error
  // bounds of Krylov methods wrt the condition number.
  real64 const norm_x_true = x_true.norm2();
  real64 const norm_x = x_comp.norm2();
  EXPECT_LT( std::fabs( norm_x / norm_x_true - 1. ), 5e-6 );
}

// END_RST_NARRATIVE

template< typename LAI >
class KrylovSolverTest : public ::testing::Test
{};

TYPED_TEST_SUITE_P( KrylovSolverTest );

TYPED_TEST_P( KrylovSolverTest, CG )
{
  testGEOSXSolvers< CGsolver, TypeParam >();
}

TYPED_TEST_P( KrylovSolverTest, BiCGSTAB )
{
  testGEOSXSolvers< BiCGSTABsolver, TypeParam >();
}

TYPED_TEST_P( KrylovSolverTest, CG_block )
{
  testGEOSXBlockSolvers< CGsolver, TypeParam >();
}

TYPED_TEST_P( KrylovSolverTest, BiCGSTAB_block )
{
  testGEOSXBlockSolvers< BiCGSTABsolver, TypeParam >();
}

REGISTER_TYPED_TEST_SUITE_P( KrylovSolverTest,
                             CG,
                             BiCGSTAB,
                             CG_block,
                             BiCGSTAB_block );

#ifdef GEOSX_USE_TRILINOS
INSTANTIATE_TYPED_TEST_SUITE_P( Trilinos, KrylovSolverTest, TrilinosInterface, );
#endif

#ifdef GEOSX_USE_HYPRE
INSTANTIATE_TYPED_TEST_SUITE_P( Hypre, KrylovSolverTest, HypreInterface, );
#endif

#ifdef GEOSX_USE_PETSC
INSTANTIATE_TYPED_TEST_SUITE_P( Petsc, KrylovSolverTest, PetscInterface, );
#endif


int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  geosx::basicSetup( argc, argv );
  int const result = RUN_ALL_TESTS();
  geosx::basicCleanup();
  return result;
}
