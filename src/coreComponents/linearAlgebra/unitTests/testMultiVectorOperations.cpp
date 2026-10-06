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
 * @file testMultiVectorOperations.cpp
 */

#include "linearAlgebra/unitTests/testLinearAlgebraUtils.hpp"
#include "linearAlgebra/utilities/MultiVectorOperations.hpp"
#include "common/GEOS_RAJA_Interface.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <vector>

using namespace geos;

/** ---------------------- Helpers ---------------------- **/

/// Fill a vector with values that differ for every vector index and every global row.
/// The local size is different on every rank, to test a non-uniform distribution.
template< typename VEC >
void createAndFill( localIndex const startSize, integer const index, VEC & x )
{
  int const rank = MpiWrapper::commRank( MPI_COMM_GEOS );
  localIndex const localSize = rank + startSize;
  globalIndex const rankOffset = rank * ( rank - 1 ) / 2 + rank * startSize;

  x.create( localSize, MPI_COMM_GEOS );
  arrayView1d< real64 > const values = x.open();
  forAll< geos::parallelDevicePolicy<> >( localSize, [=] GEOS_HOST_DEVICE ( localIndex const i )
  {
    real64 const row = static_cast< real64 >( rankOffset + i );
    values[i] = std::sin( 0.37 * ( index + 1 ) + 0.113 * row ) + 0.01 * ( index + 1 );
  } );
  x.close();
}

template< typename VEC >
std::vector< VEC > createBlock( localIndex const startSize, integer const count, integer const first )
{
  std::vector< VEC > block( count );
  for( integer k = 0; k < count; ++k )
  {
    createAndFill( startSize, first + k, block[k] );
  }
  return block;
}

template< typename VEC >
std::vector< VEC const * > pointers( std::vector< VEC > const & block )
{
  std::vector< VEC const * > result;
  for( VEC const & v : block )
  {
    result.push_back( &v );
  }
  return result;
}

/** ---------------------- Tests ---------------------- **/

template< typename LAI >
class MultiVectorOperationsTest : public ::testing::Test
{
public:
  using Vector = typename LAI::ParallelVector;
};

TYPED_TEST_SUITE_P( MultiVectorOperationsTest );

// More vectors than the tile size, so that several tiles are needed in both directions
TYPED_TEST_P( MultiVectorOperationsTest, dots )
{
  using Vector = typename TestFixture::Vector;
  localIndex constexpr startSize = 3000;
  std::vector< Vector > const X = createBlock< Vector >( startSize, 19, 0 );
  std::vector< Vector > const Y = createBlock< Vector >( startSize, 17, 40 );

  array2d< real64 > products;
  multiVectorOperations::dots( pointers( X ), pointers( Y ), products );

  ASSERT_EQ( products.size( 0 ), 19 );
  ASSERT_EQ( products.size( 1 ), 17 );
  for( integer i = 0; i < 19; ++i )
  {
    for( integer j = 0; j < 17; ++j )
    {
      real64 const expected = X[i].dot( Y[j] );
      EXPECT_NEAR( products( i, j ), expected, 1.0e-11 * ( 1.0 + std::fabs( expected ) ) ) << "pair " << i << ", " << j;
    }
  }
}

TYPED_TEST_P( MultiVectorOperationsTest, dotsOfOneVector )
{
  using Vector = typename TestFixture::Vector;
  std::vector< Vector > const X = createBlock< Vector >( 100, 5, 0 );

  array2d< real64 > products;
  multiVectorOperations::dots( pointers( X ), pointers( X ), products );
  for( integer i = 0; i < 5; ++i )
  {
    EXPECT_NEAR( products( i, i ), X[i].dot( X[i] ), 1.0e-11 * X[i].dot( X[i] ) );
    for( integer j = 0; j < 5; ++j )
    {
      EXPECT_NEAR( products( i, j ), products( j, i ), 1.0e-11 * std::fabs( products( i, j ) ) );
    }
  }
}

TYPED_TEST_P( MultiVectorOperationsTest, combine )
{
  using Vector = typename TestFixture::Vector;
  localIndex constexpr startSize = 2500;
  integer constexpr count = 21;
  std::vector< Vector > const V = createBlock< Vector >( startSize, count, 0 );

  std::vector< real64 > coefficients( count );
  for( integer k = 0; k < count; ++k )
  {
    coefficients[k] = std::cos( 0.9 * k ) - 0.2;
  }

  // Reference: a sequence of axpy
  Vector expected;
  createAndFill( startSize, 100, expected );
  Vector reference( expected );
  reference.zero();
  for( integer k = 0; k < count; ++k )
  {
    reference.axpy( coefficients[k], V[k] );
  }

  // Overwrite
  Vector result;
  createAndFill( startSize, 100, result );
  multiVectorOperations::combine( pointers( V ), coefficients, result, false );
  Vector difference( result );
  difference.axpy( -1.0, reference );
  EXPECT_LT( difference.normInf(), 1.0e-12 * ( 1.0 + reference.normInf() ) );

  // Accumulate: the previous content is kept
  Vector accumulated;
  createAndFill( startSize, 100, accumulated );
  multiVectorOperations::combine( pointers( V ), coefficients, accumulated, true );
  Vector accumulatedReference( reference );
  accumulatedReference.axpy( 1.0, expected );
  accumulated.axpy( -1.0, accumulatedReference );
  EXPECT_LT( accumulated.normInf(), 1.0e-12 * ( 1.0 + accumulatedReference.normInf() ) );
}

TYPED_TEST_P( MultiVectorOperationsTest, combineNothing )
{
  using Vector = typename TestFixture::Vector;
  Vector x;
  createAndFill( 50, 3, x );
  Vector kept( x );
  multiVectorOperations::combine( std::vector< Vector const * >{}, std::vector< real64 >{}, x, true );
  kept.axpy( -1.0, x );
  EXPECT_DOUBLE_EQ( kept.normInf(), 0.0 );
  multiVectorOperations::combine( std::vector< Vector const * >{}, std::vector< real64 >{}, x, false );
  EXPECT_DOUBLE_EQ( x.normInf(), 0.0 );
}

REGISTER_TYPED_TEST_SUITE_P( MultiVectorOperationsTest,
                             dots,
                             dotsOfOneVector,
                             combine,
                             combineNothing );

#ifdef GEOS_USE_TRILINOS
INSTANTIATE_TYPED_TEST_SUITE_P( Trilinos, MultiVectorOperationsTest, TrilinosInterface, );
#endif

#ifdef GEOS_USE_HYPRE
INSTANTIATE_TYPED_TEST_SUITE_P( Hypre, MultiVectorOperationsTest, HypreInterface, );
#endif

#ifdef GEOS_USE_PETSC
INSTANTIATE_TYPED_TEST_SUITE_P( Petsc, MultiVectorOperationsTest, PetscInterface, );
#endif

int main( int argc, char * * argv )
{
  geos::testing::LinearAlgebraTestScope scope( argc, argv );
  return RUN_ALL_TESTS();
}
