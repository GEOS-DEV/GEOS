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
#include "linearAlgebra/utilities/BlockMatrixView.hpp"
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
template<typename LAI>
void testGEOSXSolvers()
{
  // Define aliases templated on the Linear Algebra Interface (LAI).
  using Matrix = typename LAI::ParallelMatrix;
  using Vector = typename LAI::ParallelVector;


  // Use nxn cartesian mesh to generate the Laplace 2D operator.
  globalIndex n = 100;
  globalIndex N = n * n;

  // Compute a 2D Laplace operator and identity matrix
  Matrix matrix;
  Matrix identity;
  compute2DLaplaceOperator<LAI>( MPI_COMM_WORLD,
                                 n,
                                 matrix );
  computeIdentity<LAI>( MPI_COMM_WORLD,
                        n,
                        identity );

  // Define vectors
  Vector x_true,
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
  using Matrix = typename LAI::ParallelMatrix;
  using Vector = typename LAI::ParallelVector;


  // We are going to assembly the following dummy system
  // [L L] [x_true] = [b_0]
  // [L L] [x_true] = [b_1]

  globalIndex n = 100;
  globalIndex N = n * n;

  Matrix matrix;
  Matrix identity;
  compute2DLaplaceOperator<LAI>( MPI_COMM_WORLD,
                                 n,
                                 matrix );
  computeIdentity<LAI>( MPI_COMM_WORLD,
                        n,
                        identity );

  Vector x_true,
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

// END_RST_NARRATIVE

#ifdef GEOSX_USE_TRILINOS
TEST( testKrylovSolvers, withTrilinos )
{
  testGEOSXSolvers<TrilinosInterface>();
  testGEOSXBlockSolvers<TrilinosInterface>();
}
#endif

#ifdef GEOSX_USE_HYPRE
TEST( testKrylovSolvers, withHypre )
{
//  testGEOSXSolvers<HypreInterface>();
//  testGEOSXBlockSolvers<HypreInterface>();
}
#endif

#ifdef GEOSX_USE_PETSC
TEST( testKrylovSolvers, withPetsc )
{
  testGEOSXSolvers<PetscInterface>();
  testGEOSXBlockSolvers<PetscInterface>();
}
#endif


int main( int argc, char ** argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  geosx::basicSetup( argc, argv );

  int const result = RUN_ALL_TESTS();

  geosx::basicCleanup();

  return result;
}
