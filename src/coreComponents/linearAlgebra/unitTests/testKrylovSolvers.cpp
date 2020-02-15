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
#include "linearAlgebra/utilities/BlockOperatorView.hpp"
#include "linearAlgebra/utilities/BlockVectorView.hpp"
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
template<template<typename T> class SOLVER, typename LAI>
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
  compute2DLaplaceOperator<LAI>( MPI_COMM_WORLD, n, matrix );
  computeIdentity<LAI>( MPI_COMM_WORLD, n, identity );

  // Define vectors
  Vector x_true, x_comp, b;
  x_true.createWithGlobalSize( N, MPI_COMM_WORLD );
  x_comp.createWithGlobalSize( N, MPI_COMM_WORLD );
  b.createWithGlobalSize( N, MPI_COMM_WORLD );

  x_true.rand();
  x_comp.zero();

  matrix.multiply( x_true, b );

  // Solve
  SOLVER<Vector> solver( matrix, identity, 1e-8, -1 );
  solver.solve( b, x_comp );

  real64 const norm_true = x_true.norm2();
  real64 const norm_comp = x_comp.norm2();
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
template<template<typename T> class SOLVER, typename LAI>
void testGEOSXBlockSolvers()
{
  // The usual typenames
  using Matrix = typename LAI::ParallelMatrix;
  using Vector = typename LAI::ParallelVector;

  // We are going to assembly the following dummy system
  // [L L] [x_true] = [b_0]
  // [L L] [x_true] = [b_1]

  globalIndex constexpr n = 100;
  globalIndex constexpr N = n * n;

  Matrix matrix;
  Matrix identity;
  compute2DLaplaceOperator<LAI>( MPI_COMM_WORLD, n, matrix );
  computeIdentity<LAI>( MPI_COMM_WORLD, n, identity );

  BlockVector<Vector> x_true(2);
  BlockVector<Vector> x_comp( 2);
  BlockVector<Vector> b(2);

  x_true.block( 0 ).createWithGlobalSize( N, MPI_COMM_WORLD );
  x_true.block( 1 ).createWithGlobalSize( N, MPI_COMM_WORLD );
  x_comp.block( 0 ).createWithGlobalSize( N, MPI_COMM_WORLD );
  x_comp.block( 1 ).createWithGlobalSize( N, MPI_COMM_WORLD );
  b.block( 0 ).createWithGlobalSize( N, MPI_COMM_WORLD );
  b.block( 1 ).createWithGlobalSize( N, MPI_COMM_WORLD );

  x_true.rand();
  x_comp.zero();

  // Declare and allocate block matrices/vectors
  // System matrix
  BlockOperatorView<Vector> block_matrix( 2, 2 );
  BlockOperatorView<Vector> block_precon( 2, 2 );

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

  // Set right hand side.
  block_matrix.multiply( x_true, b );

  // Create block CG solver object and solve
  SOLVER< BlockVectorView< Vector > > solver( block_matrix, block_precon, 1e-8, -1 );
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

template<typename LAI>
class KrylovSolverTest : public ::testing::Test
{

};

TYPED_TEST_CASE_P( KrylovSolverTest );

TYPED_TEST_P( KrylovSolverTest, CG )
{
  testGEOSXSolvers<CGsolver, TypeParam>();
}

TYPED_TEST_P( KrylovSolverTest, BiCGSTAB )
{
  testGEOSXSolvers<BiCGSTABsolver, TypeParam>();
}

TYPED_TEST_P( KrylovSolverTest, CG_block )
{
  testGEOSXBlockSolvers<CGsolver, TypeParam>();
}

TYPED_TEST_P( KrylovSolverTest, BiCGSTAB_block )
{
  testGEOSXBlockSolvers<BiCGSTABsolver, TypeParam>();
}

REGISTER_TYPED_TEST_CASE_P( KrylovSolverTest,
                            CG,
                            BiCGSTAB,
                            CG_block,
                            BiCGSTAB_block );

#ifdef GEOSX_USE_TRILINOS
INSTANTIATE_TYPED_TEST_CASE_P( Trilinos, KrylovSolverTest, TrilinosInterface );
#endif

#ifdef GEOSX_USE_HYPRE
//INSTANTIATE_TYPED_TEST_CASE_P( Hypre, KrylovSolverTest, HypreInterface );
#endif

#ifdef GEOSX_USE_PETSC
INSTANTIATE_TYPED_TEST_CASE_P( Petsc, KrylovSolverTest, PetscInterface );
#endif


int main( int argc, char ** argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  geosx::basicSetup( argc, argv );
  int const result = RUN_ALL_TESTS();
  geosx::basicCleanup();
  return result;
}
