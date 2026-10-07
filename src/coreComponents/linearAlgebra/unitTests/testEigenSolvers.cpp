/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file testEigenSolvers.cpp
 */

#include "linearAlgebra/unitTests/testLinearAlgebraUtils.hpp"
#include "linearAlgebra/solvers/EigenSolverBase.hpp"
#include "linearAlgebra/utilities/DiagonalOperator.hpp"
#include "common/GEOS_RAJA_Interface.hpp"

#include <gtest/gtest.h>

#include <vector>

using namespace geos;

namespace
{

/// Number of unknowns, divisible by the numbers of ranks of the tests
integer constexpr numUnknowns = 60;

/// Number of repeated zero eigenvalues of the test pencil
integer constexpr numZeros = 6;

/**
 * @brief Fill the diagonals of the pencil K = diag( 0 (numZeros times), 1, 2, ..., numUnknowns - numZeros ), M = I,
 *        or K = diag( 1 (multiplicity times), 2 (multiplicity times), ... ), and of (K - sigma M)^{-1}.
 * @tparam VEC type of the vectors
 * @param[in] sigma the shift
 * @param[in] multiplicity if positive, every eigenvalue 1, 2, 3, ... is repeated this number of times, instead of
 *            the six zeros and 1, 2, ...
 * @param[in,out] k the diagonal of K
 * @param[in,out] m the diagonal of M
 * @param[in,out] inverse the diagonal of (K - sigma M)^{-1}
 *
 * This is a free function because nvcc does not accept extended lambdas in constructors.
 */
template< typename VEC >
void fillDiagonals( real64 const sigma, integer const multiplicity, VEC & k, VEC & m, VEC & inverse )
{
  int const commSize = MpiWrapper::commSize( MPI_COMM_GEOS );
  GEOS_ERROR_IF( numUnknowns % commSize != 0, "The number of unknowns must be divisible by the number of ranks" );
  localIndex const localSize = numUnknowns / commSize;
  k.create( localSize, MPI_COMM_GEOS );
  m.create( localSize, MPI_COMM_GEOS );
  inverse.create( localSize, MPI_COMM_GEOS );
  globalIndex const offset = k.ilower();

  arrayView1d< real64 > const kView = k.open();
  arrayView1d< real64 > const mView = m.open();
  arrayView1d< real64 > const iView = inverse.open();
  forAll< geos::parallelDevicePolicy<> >( localSize, [=] GEOS_HOST_DEVICE ( localIndex const i )
  {
    globalIndex const row = offset + i;
    // Six zeros, then 1, 2, 3, ...; or the eigenvalues 1, 2, 3, ... with the given multiplicity
    real64 const value = multiplicity > 0 ? static_cast< real64 >( row / multiplicity + 1 )
                                          : ( row < numZeros ? 0.0 : static_cast< real64 >( row - numZeros + 1 ) );
    kView[i] = value;
    mView[i] = 1.0;
    iView[i] = 1.0 / ( value - sigma );
  } );
  k.close();
  m.close();
  inverse.close();
}

/// Vectors of the diagonal pencil of fillDiagonals
template< typename VEC >
struct DiagonalPencil
{
  /// Create the diagonals, for a shift sigma
  explicit DiagonalPencil( real64 const sigma, integer const multiplicity = 0 )
  {
    fillDiagonals( sigma, multiplicity, k, m, inverse );
  }

  /// Diagonal of K
  VEC k;
  /// Diagonal of M
  VEC m;
  /// Diagonal of (K - sigma M)^{-1}
  VEC inverse;
};

} // namespace

template< typename LAI >
class EigenSolversTest : public ::testing::Test
{
public:
  using Vector = typename LAI::ParallelVector;

  /// Solve a diagonal problem with Arnoldi and return the result
  static EigenSolverResult solveArnoldi( real64 const sigma, integer const multiplicity, integer const numModes, integer const blockSize,
                                         integer const subspaceSize = 0 )
  {
    DiagonalPencil< Vector > const pencil( sigma, multiplicity );
    DiagonalOperator< Vector > const stiffness( pencil.k );
    DiagonalOperator< Vector > const mass( pencil.m );
    DiagonalOperator< Vector > const inverse( pencil.inverse );

    EigenSolverParameters params;
    params.solverType = EigenSolverParameters::SolverType::arnoldi;
    params.numEigenvalues = numModes;
    params.shift = sigma;
    params.blockSize = blockSize;
    params.subspaceSize = subspaceSize;
    params.tolerance = 1.0e-10;

    typename GeneralizedEigenSolver< Vector >::Problem problem{ stiffness, mass, &inverse };
    std::vector< Vector > modes;
    return GeneralizedEigenSolver< Vector >::create( params )->solve( problem, pencil.k, modes );
  }

  /// Solve the pencil with the six zeros and check the returned eigenvalues against the exact ones
  static void checkArnoldi( integer const numModes, integer const blockSize )
  {
    EigenSolverResult const result = solveArnoldi( -1.0, 0, numModes, blockSize );
    EXPECT_TRUE( result.converged );
    ASSERT_EQ( result.eigenvalues.size(), numModes );
    for( integer i = 0; i < numModes; ++i )
    {
      real64 const expected = i < numZeros ? 0.0 : static_cast< real64 >( i - numZeros + 1 );
      EXPECT_NEAR( result.eigenvalues[i], expected, 1.0e-8 ) << "eigenvalue " << i << ", block size " << blockSize;
    }
  }

  /// Solve the pencil with repeated eigenvalues 1, 2, 3, ... and check the returned eigenvalues against the exact ones
  static void checkArnoldiRepeated( real64 const sigma, integer const multiplicity, integer const numModes, integer const blockSize,
                                    integer const subspaceSize = 0 )
  {
    EigenSolverResult const result = solveArnoldi( sigma, multiplicity, numModes, blockSize, subspaceSize );
    EXPECT_TRUE( result.converged );
    ASSERT_EQ( result.eigenvalues.size(), numModes );
    for( integer i = 0; i < numModes; ++i )
    {
      EXPECT_NEAR( result.eigenvalues[i], static_cast< real64 >( i / multiplicity + 1 ), 1.0e-7 )
        << "eigenvalue " << i << ", multiplicity " << multiplicity << ", shift " << sigma << ", block size " << blockSize
        << ", subspace " << subspaceSize;
    }
  }
};

TYPED_TEST_SUITE_P( EigenSolversTest );

// With a single vector, the completeness check finds the copies of the zero eigenvalue a few at a time
TYPED_TEST_P( EigenSolversTest, arnoldiRepeatedEigenvaluesSingleVector )
{
  TestFixture::checkArnoldi( 10, 1 );
}

TYPED_TEST_P( EigenSolversTest, arnoldiRepeatedEigenvaluesBlockOfMultiplicity )
{
  TestFixture::checkArnoldi( 10, 6 );
}

TYPED_TEST_P( EigenSolversTest, arnoldiRepeatedEigenvaluesSmallBlock )
{
  TestFixture::checkArnoldi( 10, 2 );
}

TYPED_TEST_P( EigenSolversTest, arnoldiOnlyRepeatedEigenvalues )
{
  TestFixture::checkArnoldi( 6, 1 );
}

// The eigenvalues are close to each other compared with the distance to the shift, so that the Ritz values of
// the check cycle are clustered. A check that looks only at the best new Ritz value can miss the second copy
// of a double eigenvalue.
TYPED_TEST_P( EigenSolversTest, arnoldiClusteredDoubleEigenvalues )
{
  for( integer const numModes : { 2, 3, 4, 5, 7 } )
  {
    SCOPED_TRACE( numModes );
    TestFixture::checkArnoldiRepeated( -1.0e3, 2, numModes, 1 );
  }
  TestFixture::checkArnoldiRepeated( -1.0e5, 2, 5, 1 );
  TestFixture::checkArnoldiRepeated( -1.0e3, 2, 5, 1, 30 );
}

// Multiplicities of three and four, found one copy at a time with a single vector
TYPED_TEST_P( EigenSolversTest, arnoldiClusteredHigherMultiplicity )
{
  for( integer const multiplicity : { 3, 4 } )
  {
    for( real64 const sigma : { -1.0e2, -1.0e4 } )
    {
      for( integer const numModes : { 2, 4, 5, 8, 10 } )
      {
        SCOPED_TRACE( ::testing::Message() << "multiplicity " << multiplicity << ", shift " << sigma << ", modes " << numModes );
        TestFixture::checkArnoldiRepeated( sigma, multiplicity, numModes, 1 );
      }
    }
  }
}

// Few modes with the default basis size must converge
TYPED_TEST_P( EigenSolversTest, arnoldiFewModesDefaultBasis )
{
  for( integer const numModes : { 1, 2, 3 } )
  {
    SCOPED_TRACE( numModes );
    TestFixture::checkArnoldiRepeated( -1.0e3, 2, numModes, 1 );
  }
}

REGISTER_TYPED_TEST_SUITE_P( EigenSolversTest,
                             arnoldiClusteredDoubleEigenvalues,
                             arnoldiClusteredHigherMultiplicity,
                             arnoldiFewModesDefaultBasis,
                             arnoldiRepeatedEigenvaluesSingleVector,
                             arnoldiRepeatedEigenvaluesBlockOfMultiplicity,
                             arnoldiRepeatedEigenvaluesSmallBlock,
                             arnoldiOnlyRepeatedEigenvalues );

#ifdef GEOS_USE_TRILINOS
INSTANTIATE_TYPED_TEST_SUITE_P( Trilinos, EigenSolversTest, TrilinosInterface, );
#endif

#ifdef GEOS_USE_HYPRE
INSTANTIATE_TYPED_TEST_SUITE_P( Hypre, EigenSolversTest, HypreInterface, );
#endif

#ifdef GEOS_USE_PETSC
INSTANTIATE_TYPED_TEST_SUITE_P( Petsc, EigenSolversTest, PetscInterface, );
#endif

int main( int argc, char * * argv )
{
  geos::testing::LinearAlgebraTestScope scope( argc, argv );
  return RUN_ALL_TESTS();
}
