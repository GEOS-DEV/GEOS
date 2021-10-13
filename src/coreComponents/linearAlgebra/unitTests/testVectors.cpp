/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file testVectors.cpp
 */

#include "linearAlgebra/unitTests/testLinearAlgebraUtils.hpp"

#include <gtest/gtest.h>

using namespace geosx;

/** ---------------------- Helpers ---------------------- **/

array1d< real64 > makeLocalValues( localIndex const localSize, globalIndex const rankOffset )
{
  array1d< real64 > localValues( localSize );
  for( localIndex i = 0; i < localSize; ++i )
  {
    globalIndex const row = rankOffset + i;
    localValues[i] = std::pow( -1.0, row ) * ( 1.0 + row );
  }
  return localValues;
}

array1d< real64 > makeLocalValuesUniform( localIndex const localSize )
{
  int const rank = MpiWrapper::commRank( MPI_COMM_GEOSX );
  localIndex const rankOffset = localSize * rank;

  return makeLocalValues( localSize, rankOffset );
}

array1d< real64 > makeLocalValuesNonUniform( localIndex const startSize )
{
  int const rank = MpiWrapper::commRank( MPI_COMM_GEOSX );
  localIndex const localSize = rank + startSize;
  globalIndex const rankOffset = ( rank + startSize ) * ( rank + startSize - 1 ) / 2;

  return makeLocalValues( localSize, rankOffset );
}

namespace ops
{

auto const identity = []( real64 const x ){ return x; };

struct multiply
{
  real64 factor;
  real64 operator()( real64 const x ) const { return factor * x; }
};

auto const reciprocal = []( real64 const x ){ return 1.0 / x; };

}

template< typename OP_X = decltype( ops::identity ), typename OP_Y = decltype( ops::identity ) >
void compareValues( arrayView1d< real64 > const & x,
                    arrayView1d< real64 > const & y,
                    bool exact = true,
                    OP_X const op_x = ops::identity,
                    OP_Y const op_y = ops::identity )
{
  EXPECT_EQ( x.size(), y.size() );
  x.move( LvArray::MemorySpace::host, false );
  y.move( LvArray::MemorySpace::host, false );
  for( localIndex i = 0; i < x.size(); ++i )
  {
    if( exact )
    {
      EXPECT_EQ( op_x( x[i] ), op_y( y[i] ) );
    }
    else
    {
      EXPECT_DOUBLE_EQ( op_x( x[i] ), op_y( y[i] ) );
    }
  }
}

template< typename VEC, typename OP_X = decltype( ops::identity ), typename OP_Y = decltype( ops::identity ) >
void compareVectors( VEC const & x,
                     VEC const & y,
                     bool exact = true,
                     OP_X const op_x = ops::identity,
                     OP_Y const op_y = ops::identity )
{
  EXPECT_TRUE( MpiWrapper::commCompare( y.getComm(), x.getComm() ) );
  EXPECT_EQ( y.localSize(), x.localSize() );
  EXPECT_EQ( y.globalSize(), x.globalSize() );

  array1d< real64 > xval( x.localSize() );
  x.extract( xval );

  array1d< real64 > yval( y.localSize() );
  y.extract( yval );

  compareValues( xval, yval, exact, op_x, op_y );
}

/** ---------------------- Tests ---------------------- **/

template< typename LAI >
class VectorTest : public ::testing::Test
{};

TYPED_TEST_SUITE_P( VectorTest );

TYPED_TEST_P( VectorTest, createWithLocalSize )
{
  using Vector = typename TypeParam::ParallelVector;

  MPI_Comm const comm = MPI_COMM_GEOSX;
  int const rank = MpiWrapper::commRank( comm );
  int const nproc = MpiWrapper::commSize( comm );

  localIndex const localSize = rank + 1;
  globalIndex const globalSize = ( nproc + 1 ) * nproc / 2;
  globalIndex const rankOffset = ( rank + 1 ) * rank / 2;

  Vector x;
  x.createWithLocalSize( localSize, comm );

  EXPECT_TRUE( x.ready() );
  EXPECT_TRUE( MpiWrapper::commCompare( x.getComm(), comm ) );
  EXPECT_EQ( x.localSize(), localSize );
  EXPECT_EQ( x.globalSize(), globalSize );
  EXPECT_EQ( x.ilower(), rankOffset );
  EXPECT_EQ( x.iupper(), rankOffset + localSize );
}

TYPED_TEST_P( VectorTest, createWithGlobalSize )
{
  using Vector = typename TypeParam::ParallelVector;

  MPI_Comm const comm = MPI_COMM_GEOSX;
  int const rank = MpiWrapper::commRank( comm );
  int const nproc = MpiWrapper::commSize( comm );

  localIndex const localSize = 3;
  globalIndex const globalSize = localSize * nproc;
  globalIndex const rankOffset = localSize * rank;

  Vector x;
  x.createWithGlobalSize( globalSize, comm );

  EXPECT_TRUE( x.ready() );
  EXPECT_TRUE( MpiWrapper::commCompare( x.getComm(), comm ) );
  EXPECT_EQ( x.localSize(), localSize );
  EXPECT_EQ( x.globalSize(), globalSize );
  EXPECT_EQ( x.ilower(), rankOffset );
  EXPECT_EQ( x.iupper(), rankOffset + localSize );
}

TYPED_TEST_P( VectorTest, create )
{
  using Vector = typename TypeParam::ParallelVector;

  array1d< real64 > const valuesInitial = makeLocalValuesNonUniform( 1 );
  localIndex const localSize = valuesInitial.size();
  globalIndex const globalSize = MpiWrapper::sum( localSize, MPI_COMM_GEOSX );

  Vector x;
  x.create( valuesInitial, MPI_COMM_GEOSX );

  EXPECT_TRUE( x.ready() );
  EXPECT_TRUE( MpiWrapper::commCompare( x.getComm(), MPI_COMM_GEOSX ) );
  EXPECT_EQ( x.localSize(), localSize );
  EXPECT_EQ( x.globalSize(), globalSize );

  array1d< real64 > const valuesExtracted( localSize );
  x.extract( valuesExtracted );

  compareValues( valuesExtracted, valuesInitial );
}

TYPED_TEST_P( VectorTest, copyConstruction )
{
  using Vector = typename TypeParam::ParallelVector;

  array1d< real64 > const localValues = makeLocalValuesUniform( 3 );

  Vector x;
  x.create( localValues, MPI_COMM_GEOSX );

  Vector y( x );

  // Test that values are equal after copy contruction
  EXPECT_TRUE( y.ready() );
  compareVectors( x, y );

  // Test that vectors don't share data after copy construction
  real64 const factor = 0.5;
  x.scale( factor );

  compareVectors( x, y, true, ops::identity, ops::multiply{factor} );
}

TYPED_TEST_P( VectorTest, moveConstruction )
{
  using Vector = typename TypeParam::ParallelVector;

  array1d< real64 > const valuesInitial = makeLocalValuesUniform( 3 );
  localIndex const localSize = valuesInitial.size();
  globalIndex const globalSize = MpiWrapper::sum( localSize );

  Vector x;
  x.create( valuesInitial, MPI_COMM_GEOSX );

  Vector y( std::move( x ) );

  EXPECT_TRUE( y.ready() );
  EXPECT_TRUE( MpiWrapper::commCompare( y.getComm(), MPI_COMM_GEOSX ) );
  EXPECT_EQ( y.localSize(), localSize );
  EXPECT_EQ( y.globalSize(), globalSize );

  array1d< real64 > const valuesExtracted( localSize );
  y.extract( valuesExtracted );

  compareValues( valuesExtracted, valuesInitial );
}

TYPED_TEST_P( VectorTest, copy )
{
  using Vector = typename TypeParam::ParallelVector;

  localIndex const localSize = 3;

  Vector x;
  x.createWithLocalSize( localSize, MPI_COMM_GEOSX );
  x.rand();

  Vector y;
  y.createWithLocalSize( localSize, MPI_COMM_GEOSX );
  y.copy( x );

  compareVectors( x, y );
}

TYPED_TEST_P( VectorTest, extract )
{
  using Vector = typename TypeParam::ParallelVector;

  array1d< real64 > const valuesInitial = makeLocalValuesUniform( 3 );

  Vector x;
  x.create( valuesInitial, MPI_COMM_GEOSX );

  array1d< real64 > const valuesExtracted( valuesInitial.size() );
  x.extract( valuesExtracted );

  compareValues( valuesExtracted, valuesInitial );
}

TYPED_TEST_P( VectorTest, localGlobalRowID )
{
  using Vector = typename TypeParam::ParallelVector;

  int const rank = MpiWrapper::commRank( MPI_COMM_GEOSX );

  localIndex const localSize = rank + 1;
  globalIndex const offset = ( rank + 1 ) * rank / 2;

  Vector x;
  x.createWithLocalSize( localSize, MPI_COMM_GEOSX );

  for( localIndex i = 0; i < localSize; ++i )
  {
    EXPECT_EQ( x.getGlobalRowID( i ), offset + i );
    EXPECT_EQ( x.getLocalRowID( offset + i ), i );
  }
}

TYPED_TEST_P( VectorTest, getSingleValue )
{
  using Vector = typename TypeParam::ParallelVector;

  localIndex const localSize = 3;
  array1d< real64 > const localValues = makeLocalValuesUniform( 3 );
  array1d< real64 > const localValuesCopy = localValues;

  Vector x;
  x.create( localValues, MPI_COMM_GEOSX );

  for( localIndex i = 0; i < localSize; ++i )
  {
    EXPECT_EQ( x.get( x.getGlobalRowID( i ) ), localValuesCopy[i] );
  }
}

TYPED_TEST_P( VectorTest, getMultipleValues )
{
  using Vector = typename TypeParam::ParallelVector;

  localIndex const localSize = 3;
  array1d< real64 > const valuesInitial = makeLocalValuesUniform( 3 );

  Vector x;
  x.create( valuesInitial, MPI_COMM_GEOSX );

  array1d< globalIndex > const rowIndices( localSize );
  array1d< real64 > const valuesExtracted( localSize );

  for( localIndex i = 0; i < localSize; ++i )
  {
    rowIndices[i] = x.getGlobalRowID( i );
  }

  x.get( rowIndices, valuesExtracted );
  compareValues( valuesExtracted, valuesInitial );
}

TYPED_TEST_P( VectorTest, setSingleValue )
{
  using Vector = typename TypeParam::ParallelVector;

  localIndex const localSize = 3;

  Vector x;
  x.createWithLocalSize( localSize, MPI_COMM_GEOSX );

  x.open();
  for( localIndex i = 0; i < localSize; ++i )
  {
    x.set( x.getGlobalRowID( i ), x.ilower() + i );
  }
  x.close();

  for( localIndex i = 0; i < localSize; ++i )
  {
    EXPECT_EQ( x.get( x.getGlobalRowID( i ) ), x.ilower() + i );
  }
}

TYPED_TEST_P( VectorTest, addSingleValue )
{
  using Vector = typename TypeParam::ParallelVector;

  localIndex const localSize = 3;

  Vector x;
  x.createWithLocalSize( localSize, MPI_COMM_GEOSX );
  x.set( 1.0 );

  x.open();
  for( localIndex i = 0; i < localSize; ++i )
  {
    x.add( x.getGlobalRowID( i ), x.ilower() + i );
  }
  x.close();

  for( localIndex i = 0; i < localSize; ++i )
  {
    EXPECT_EQ( x.get( x.getGlobalRowID( i ) ), 1.0 + x.ilower() + i );
  }
}

TYPED_TEST_P( VectorTest, setAllValues )
{
  using Vector = typename TypeParam::ParallelVector;

  localIndex const localSize = 3;
  real64 const value = 1.23;

  Vector x;
  x.createWithLocalSize( localSize, MPI_COMM_GEOSX );
  x.set( value );

  array1d< real64 > const valuesExtracted( localSize );
  x.extract( valuesExtracted );
  valuesExtracted.move( LvArray::MemorySpace::host, false );

  for( localIndex i = 0; i < localSize; ++i )
  {
    EXPECT_EQ( valuesExtracted[i], value );
  }
}

TYPED_TEST_P( VectorTest, zeroAllValues )
{
  using Vector = typename TypeParam::ParallelVector;

  localIndex const localSize = 3;

  Vector x;
  x.createWithLocalSize( localSize, MPI_COMM_GEOSX );
  x.zero();

  array1d< real64 > const valuesExtracted( localSize );
  x.extract( valuesExtracted );
  valuesExtracted.move( LvArray::MemorySpace::host, false );

  for( localIndex i = 0; i < localSize; ++i )
  {
    EXPECT_EQ( valuesExtracted[i], 0.0 );
  }
}

TYPED_TEST_P( VectorTest, scaleValues )
{
  using Vector = typename TypeParam::ParallelVector;

  localIndex const localSize = 3;
  array1d< real64 > const valuesInitial = makeLocalValuesUniform( 3 );
  array1d< real64 > const valuesCopy = valuesInitial;

  Vector x;
  x.create( valuesInitial, MPI_COMM_GEOSX );

  real64 const factor = 0.5;
  x.scale( factor );

  array1d< real64 > const valuesExtracted( localSize );
  x.extract( valuesExtracted );

  compareValues( valuesExtracted, valuesCopy, true, ops::identity, ops::multiply{factor} );
}

TYPED_TEST_P( VectorTest, reciprocal )
{
  using Vector = typename TypeParam::ParallelVector;

  localIndex const localSize = 3;
  array1d< real64 > const valuesInitial = makeLocalValuesUniform( 3 );
  array1d< real64 > const valuesCopy = valuesInitial;

  Vector x;
  x.create( valuesInitial, MPI_COMM_GEOSX );

  x.reciprocal();

  array1d< real64 > const valuesExtracted( localSize );
  x.extract( valuesExtracted );

  compareValues( valuesExtracted, valuesCopy, true, ops::reciprocal );
}

TYPED_TEST_P( VectorTest, dotProduct )
{
  using Vector = typename TypeParam::ParallelVector;

  array1d< real64 > const xValues = makeLocalValuesUniform( 3 );
  array1d< real64 > const yValues = xValues;

  Vector x;
  x.create( xValues, MPI_COMM_GEOSX );

  Vector y;
  y.create( yValues, MPI_COMM_GEOSX );
  y.reciprocal();

  real64 const dp = x.dot( y );
  EXPECT_DOUBLE_EQ( dp, x.globalSize() );
}

TYPED_TEST_P( VectorTest, axpy )
{
  using Vector = typename TypeParam::ParallelVector;

  real64 const alpha = 2.0;

  array1d< real64 > const xValues = makeLocalValuesUniform( 3 );
  array1d< real64 > const yValues = xValues;

  Vector x;
  x.create( xValues, MPI_COMM_GEOSX );

  Vector y;
  y.create( yValues, MPI_COMM_GEOSX );

  x.axpy( alpha, y );

  compareVectors( y, x, false, ops::multiply{ 1.0 + alpha } );
}

TYPED_TEST_P( VectorTest, axpby )
{
  using Vector = typename TypeParam::ParallelVector;

  real64 const alpha = 2.0;
  real64 const beta = 3.0;

  array1d< real64 > const xValues = makeLocalValuesUniform( 3 );
  array1d< real64 > const yValues = xValues;

  Vector x;
  x.create( xValues, MPI_COMM_GEOSX );

  Vector y;
  y.create( yValues, MPI_COMM_GEOSX );

  x.axpby( alpha, y, beta );

  compareVectors( y, x, false, ops::multiply{ alpha + beta } );
}

TYPED_TEST_P( VectorTest, norm1 )
{
  using Vector = typename TypeParam::ParallelVector;

  array1d< real64 > const xValues = makeLocalValuesUniform( 3 );

  Vector x;
  x.create( xValues, MPI_COMM_GEOSX );
  x.scale( -1.0 );

  real64 const normTrue = ( x.globalSize() + 1 ) * x.globalSize() / 2;

  EXPECT_DOUBLE_EQ( x.norm1(), normTrue );
}

TYPED_TEST_P( VectorTest, norm2 )
{
  using Vector = typename TypeParam::ParallelVector;

  array1d< real64 > const xValues = makeLocalValuesUniform( 3 );

  Vector x;
  x.create( xValues, MPI_COMM_GEOSX );
  x.scale( -1.0 );

  real64 const normTrue = std::sqrt( ( 2 * x.globalSize() + 1 ) * ( x.globalSize() + 1 ) * x.globalSize() / 6 );

  EXPECT_DOUBLE_EQ( x.norm2(), normTrue );
}

TYPED_TEST_P( VectorTest, normInf )
{
  using Vector = typename TypeParam::ParallelVector;

  array1d< real64 > const xValues = makeLocalValuesUniform( 3 );

  Vector x;
  x.create( xValues, MPI_COMM_GEOSX );
  x.scale( -1.0 );

  real64 const normTrue = x.globalSize();

  EXPECT_DOUBLE_EQ( x.normInf(), normTrue );
}

REGISTER_TYPED_TEST_SUITE_P( VectorTest,
                             createWithLocalSize,
                             createWithGlobalSize,
                             create,
                             copyConstruction,
                             moveConstruction,
                             copy,
                             extract,
                             localGlobalRowID,
                             getSingleValue,
                             getMultipleValues,
                             setSingleValue,
                             addSingleValue,
                             setAllValues,
                             zeroAllValues,
                             scaleValues,
                             reciprocal,
                             dotProduct,
                             axpy,
                             axpby,
                             norm1,
                             norm2,
                             normInf );

#ifdef GEOSX_USE_TRILINOS
INSTANTIATE_TYPED_TEST_SUITE_P( Trilinos, VectorTest, TrilinosInterface, );
#endif

#ifdef GEOSX_USE_HYPRE
INSTANTIATE_TYPED_TEST_SUITE_P( Hypre, VectorTest, HypreInterface, );
#endif

#ifdef GEOSX_USE_PETSC
INSTANTIATE_TYPED_TEST_SUITE_P( Petsc, VectorTest, PetscInterface, );
#endif

int main( int argc, char * * argv )
{
  geosx::testing::LinearAlgebraTestScope scope( argc, argv );
  return RUN_ALL_TESTS();
}
