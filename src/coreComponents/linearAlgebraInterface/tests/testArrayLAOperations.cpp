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

//#ifdef __clang__
//#pragma clang diagnostic push
//#pragma clang diagnostic ignored "-Wglobal-constructors"
//#pragma clang diagnostic ignored "-Wexit-time-destructors"
//#if __clang_major__ >= 5 && !defined(__APPLE__)
//#pragma clang diagnostic ignored "-Wzero-as-null-pointer-constant"
//#endif
//#endif
//
//#include <gtest/gtest.h>

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wglobal-constructors"
#pragma clang diagnostic ignored "-Wexit-time-destructors"
#pragma clang diagnostic ignored "-Wused-but-marked-unused"
#endif

#include "gtest/gtest.h"

#ifdef __clang__
#pragma clang diagnostic push
#define __null nullptr
#endif

#include "common/DataTypes.hpp"
#include <random>
#include <BlasLapackLA.hpp>

using namespace geosx;

extern "C"
{
  #include "BlasLapackFunctions.hpp"
}

using INDEX_TYPE = std::ptrdiff_t;

// Generate a random GEOSX array1d
void vectorRand( array1d<real64> & x,
                 std::uniform_real_distribution<real64> & distribution,
                 std::default_random_engine & generator );

// Generate a random GEOSX array2d
void matrixRand( array2d<real64> & A,
                 std::uniform_real_distribution<real64> & distribution,
                 std::default_random_engine & generator );

static real64 machinePrecision = 10. * std::numeric_limits<real64>::epsilon();

//////////////////////////////////////////////////////////

template<typename LAI>
void VectorNorm1_Test()
{
  array1d<real64> v;
  INDEX_TYPE N = 24;
  v.resize( N );

  // Populate vector with random coefficients
  int IDIST = 2;
  LAI::vectorRand(v, IDIST);

  // Compute norm1
  real64 norm1 = 0.0;
  for( INDEX_TYPE i = 0 ; i < N ; ++i )
  {
    norm1 += std::abs(v(i));
  }

  // Check
  EXPECT_NEAR( norm1,
               LAI::vectorNorm1( v ),
               norm1 * machinePrecision );
}

template<typename LAI>
void VectorNorm2_Test()
{
  array1d<real64> v;
  INDEX_TYPE N = 24;
  v.resize( N );

  // Populate vector with random coefficients
  int IDIST = 2;
  LAI::vectorRand(v, IDIST);

  // Compute norm2
  real64 norm2 = 0.0;
  for( INDEX_TYPE i = 0 ; i < N ; ++i )
  {
    norm2 += v(i)*v(i);
  }
  norm2 = std::sqrt(norm2);

  // Check
  EXPECT_NEAR( norm2,
               LAI::vectorNorm2( v ),
               norm2 * machinePrecision );
}

template<typename LAI>
void VectorNormInf_Test()
{
  array1d<real64> v;
  INDEX_TYPE N = 24;
  v.resize( N );

  // Populate vector with random coefficients
  int IDIST = 2;
  LAI::vectorRand(v, IDIST);

  // Compute normInf
  real64 normInf = std::abs(v(0));

  if (N > 1)
  {
    for( INDEX_TYPE i = 1 ; i < N ; ++i )
    {
      normInf = std::max( normInf, std::abs(v(i)) );
    }
  }

  // Check
  EXPECT_NEAR( normInf,
               LAI::vectorNormInf( v ),
               normInf * machinePrecision );
}

template<typename LAI>
void Determinant_Test()
{
  array2d<real64> Laplacian1d;
  real64 determinant;

  for (INDEX_TYPE N = 1; N <= 6; ++ N)
  {
    // Construct 1d Laplacian
    Laplacian1d.resize( N, N );
    Laplacian1d = 0;
    for( INDEX_TYPE i = 0 ; i < N ; ++i )
    {
      for( INDEX_TYPE j = 0 ; j < N ; ++j )
      {
        if( i == j )
        {
          Laplacian1d( i, i ) = 2;
        }
        else if( std::abs( i - j ) == 1 )
        {
          Laplacian1d( i, j ) = -1;
        }
      }
    }

    // Check
    determinant = static_cast<real64>(N+1);
    EXPECT_NEAR( determinant,
                 LAI::determinant( Laplacian1d ),
                 determinant * machinePrecision );
  }
}

template<typename LAI>
void MatrixNormInf_Test()
{
  array2d<real64> A;
  INDEX_TYPE M = 6;
  INDEX_TYPE N = 24;
  A.resize( M, N );

  // Populate matrix with random coefficients
  int IDIST = 2;
  LAI::matrixRand(A, IDIST);

  // Compute normInf
  real64 normInf = 0.0;
  for( INDEX_TYPE i = 0 ; i < M ; ++i )
  {
    real64 tmp = 0.0;
    for( INDEX_TYPE j = 0 ; j < N ; ++j )
    {
      tmp += std::abs( A(i,j) );
    }
    normInf = std::max( normInf, tmp);
  }

  // Check
  EXPECT_NEAR( normInf,
               LAI::matrixNormInf( A ),
               normInf * machinePrecision );
}

template<typename LAI>
void MatrixNorm1_Test()
{
  array2d<real64> A;
  INDEX_TYPE M = 6;
  INDEX_TYPE N = 24;
  A.resize( M, N );

  // Populate matrix with random coefficients
  int IDIST = 2;
  LAI::matrixRand(A, IDIST);

  // Compute norm1
  array1d<real64> tmp(N);
  tmp = 0;
  for( INDEX_TYPE i = 0 ; i < M ; ++i )
  {
    for( INDEX_TYPE j = 0 ; j < N ; ++j )
    {
      tmp(j) += std::abs( A(i,j) );
    }
  }
  real64 norm1 = tmp(0);
  for( INDEX_TYPE j = 0 ; j < N ; ++j )
  {
    norm1 = std::max( norm1, tmp(j)) ;
  }


  // Check
  EXPECT_NEAR( norm1,
               LAI::matrixNorm1( A ),
               norm1 * machinePrecision );
}

template<typename LAI>
void MatrixNormFrobenius_Test()
{
  array2d<real64> A;
  INDEX_TYPE M = 6;
  INDEX_TYPE N = 24;
  A.resize( M, N );

  // Populate matrix with random coefficients
  int IDIST = 2;
  LAI::matrixRand(A, IDIST);

  // Compute normFrobeniusm
  real64 normFrobenius = 0.0;
  for( INDEX_TYPE i = 0 ; i < M ; ++i )
  {
    for( INDEX_TYPE j = 0 ; j < N ; ++j )
    {
      normFrobenius += A(i,j)*A(i,j);
    }
  }
  normFrobenius = std::sqrt(normFrobenius);

  // Check
  EXPECT_NEAR( normFrobenius,
               LAI::matrixNormFrobenius( A ),
               normFrobenius * machinePrecision );
}

template<typename LAI>
void VectorVectorAdd_Test()
{
  array1d<real64> v1;
  array1d<real64> v2;
  INDEX_TYPE N = 24;
  v.resize( N );

  // Populate vector with random coefficients
  int IDIST = 2;
  LAI::vectorRand(v, IDIST);

  // Compute normInf
  real64 normInf = std::abs(v(0));

  if (N > 1)
  {
    for( INDEX_TYPE i = 1 ; i < N ; ++i )
    {
      normInf = std::max( normInf, std::abs(v(i)) );
    }
  }

  // Check
  EXPECT_NEAR( normInf,
               LAI::vectorNormInf( v ),
               normInf * machinePrecision );
}

TEST( Array1D, VectorNorm1)
{
  VectorNorm1_Test<BlasLapackLA>();
}

TEST( Array1D, VectorNorm2)
{
  VectorNorm2_Test<BlasLapackLA>();
}

TEST( Array1D, VectorNormInf)
{
  VectorNormInf_Test<BlasLapackLA>();
}

TEST( Array2D, MatrixNormInf)
{
  MatrixNormInf_Test<BlasLapackLA>();
}

TEST( Array2D, determinant)
{
  Determinant_Test<BlasLapackLA>();
}

TEST( Array2D, MatrixNorm1)
{
  MatrixNorm1_Test<BlasLapackLA>();
}

TEST( Array2D, MatrixNormFrobenius)
{
  MatrixNormFrobenius_Test<BlasLapackLA>();
}

TEST( Array1D, VectorVectorAdd)
{
  VectorVectorAdd_Test<BlasLapackLA>();
}

/////////////////////////////////////////////////////////



template<typename LAI>
void testArray1dLA()
{

  // Repeat the following step for vectors of increasing size:
  //
  // a. Resize vectors v1, v2, and v3 to be of equal size
  //    (not initialized to 0)
  // b. Fill v1 with random coefficient
  // c. Initialize v2 to have zero coefficients
  // d. Compute v2 = v2 + alfa*v1
  // e. Compute v1 = alfa*v1
  // f. Copy v1 in v3
  // g. Check that (norm(v2, type) == norm(v3, type)), type = \{1-norm, 2-norm, Infinity-norm)
  // h. check that ( (norm(v2, norm-2) - sqrt(v2_dot_v2)) / norm(v2, norm-2)) < epsilon

  array1d<real64> v1, v2, v3;
  std::default_random_engine generator;
  std::uniform_real_distribution<real64> distribution( 0.0, 1.0 );
  real64 alfa = distribution( generator );

  for( int iSize = 1 ; iSize <= 10 ; ++iSize )
  {
    // a.
    v1.resize( iSize );
    v2.resize( iSize );
    v3.resize( iSize );

    // b.
    for( int i = 0 ; i < iSize ; ++i )
    {
      v1[i] = distribution( generator );
    }

    // c.
    v2 = 0.0;

    // d.
    LAI::vectorVectorAdd( v1, v2, alfa );

    // e.
    LAI::vectorScale( alfa, v1 );

    // f.
    LAI::vectorCopy( v1, v3 );

    // g.
    EXPECT_NEAR( LAI::vectorNorm1( v2 ),
                 LAI::vectorNorm1( v3 ),
                 machinePrecision * LAI::vectorNorm1( v2 ) );
    EXPECT_NEAR( LAI::vectorNorm2( v2 ),
                 LAI::vectorNorm2( v3 ),
                 machinePrecision * LAI::vectorNorm2( v2 ) );
    EXPECT_NEAR( LAI::vectorNormInf( v2 ),
                 LAI::vectorNormInf( v3 ),
                 machinePrecision * LAI::vectorNormInf( v2 ) );

    // i.
    real64 beta = LAI::vectorDot( v2, v2 );
    EXPECT_NEAR( std::sqrt( LAI::vectorDot( v2, v2 ) ),
                 LAI::vectorNorm2( v2 ),
                 machinePrecision * LAI::vectorNorm2( v2 ) );

  }

}

template<typename LAI>
void testArray2dLA()
{

  array2d<real64> A, B, C, D, LHS, RHS;
  std::default_random_engine generator;
  std::uniform_real_distribution<real64> distribution( 0.0, 1.0 );
  real64 alfa = distribution( generator );
  real64 beta = distribution( generator );

  localIndex MA, NA, MB, NB, ND, K;

  // Repeat the following step for vectors of increasing size:
  //
  // a. compute LHS = ( alfa*A*B + beta*C ) * D
  // b. compute RHS = alfa*A*(B*D) + beta*C*D
  // c. check that (norm(LHS) == norm(RHS, type)), with type = \{Infinity-norm, 1-norm, Frobenius-norm, )

  MA = 10;
  NB = 10;
  ND = 10;
  K = 20;

  for( localIndex mA = 1 ; mA <= MA ; ++mA )
  {
    for( localIndex nB = 1 ; nB <= NB ; ++nB )
    {
      for( localIndex nD = 1 ; nD <= ND ; ++nD )
      {
        for( localIndex k = 1 ; k <= K ; ++k )
        {
          // Resize matrix operators
          A.resize( mA, k );
          B.resize( k, nB );
          C.resize( A.size( 0 ), B.size( 1 ) );
          D.resize( B.size( 1 ), nD );
          LHS.resize( A.size( 0 ), D.size( 1 ) );
          RHS.resize( A.size( 0 ), D.size( 1 ) );

          // Populate A, B, C, and D with uniformly distributed random
          // coefficients
          matrixRand( A,
                      distribution,
                      generator );
          matrixRand( B,
                      distribution,
                      generator );
          matrixRand( C,
                      distribution,
                      generator );
          matrixRand( D,
                      distribution,
                      generator );

          // Compute tmp1 = ( alfa*A*B + beta*C )
          array2d<real64> tmp1( C );
          LAI::matrixMatrixMultiply( A, B, tmp1, alfa, beta );

          // Compute LHS = tmp1*D
          LAI::matrixMatrixMultiply( tmp1, D, LHS );

          // Compute tmp2 = B*D
          array2d<real64> tmp2( B.size( 0 ), D.size( 1 ) );
          LAI::matrixMatrixMultiply( B, D, tmp2 );

          // Compute RHS = alfa*A*tmp2
          LAI::matrixMatrixMultiply( A, tmp2, RHS, alfa );

          // Compute RHS = RHS + beta*C*D
          LAI::matrixMatrixMultiply( C, D, RHS, beta, 1. );

          // Check norms
          EXPECT_NEAR( LAI::matrixNormInf( LHS ),
                       LAI::matrixNormInf( RHS ),
                       LAI::matrixNormInf( LHS ) * machinePrecision );
          EXPECT_NEAR( LAI::matrixNorm1( LHS ),
                       LAI::matrixNorm1( RHS ),
                       LAI::matrixNorm1( LHS ) * machinePrecision );
          EXPECT_NEAR( LAI::matrixNormFrobenius( LHS ),
                       LAI::matrixNormFrobenius( RHS ),
                       LAI::matrixNormFrobenius( LHS ) * machinePrecision );
        }
      }
    }
  }
  // Repeat the following step for vectors of increasing size:
  //
  // a. compute LHS = ( alfa*A^T*B + beta*C ) * D
  // b. compute RHS = alfa*A^T*(B*D) + beta*C*D
  // c. check that (norm(LHS - RHS, type) / norm(LHS, type)) < epsilon
  //    with type = \{Infinity-norm, 1-norm, Frobenius-norm, )

  NA = 10;
  NB = 10;
  ND = 10;
  K = 20;

  for( localIndex nA = 1 ; nA <= NA ; ++nA )
  {
    for( localIndex nB = 1 ; nB <= NB ; ++nB )
    {
      for( localIndex nD = 1 ; nD <= ND ; ++nD )
      {
        for( localIndex k = 1 ; k <= K ; ++k )
        {

          // Resize matrix operators
          A.resize( k, nA );
          B.resize( k, nB );
          C.resize( A.size( 1 ), B.size( 1 ) );
          D.resize( B.size( 1 ), nD );
          LHS.resize( A.size( 1 ), D.size( 1 ) );
          RHS.resize( A.size( 1 ), D.size( 1 ) );

          // Populate A, B, C, and D with uniformly distributed random
          // coefficients
          // Populate A, B, C, and D with uniformly distributed random
          // coefficients
          matrixRand( A,
                      distribution,
                      generator );
          matrixRand( B,
                      distribution,
                      generator );
          matrixRand( C,
                      distribution,
                      generator );
          matrixRand( D,
                      distribution,
                      generator );

          // Compute tmp1 = ( alfa*A^T*B + beta*C )
          array2d<real64> tmp1( C );
          LAI::matrixTMatrixMultiply( A, B, tmp1, alfa, beta );

          // Compute LHS = tmp*D
          LAI::matrixMatrixMultiply( tmp1, D, LHS );

          // Compute tmp2 = B*D
          array2d<real64> tmp2( B.size( 0 ), D.size( 1 ) );
          LAI::matrixMatrixMultiply( B, D, tmp2 );

          // Compute RHS = alfa*A^T*tmp2
          LAI::matrixTMatrixMultiply( A, tmp2, RHS, alfa );

          // Compute RHS = RHS + beta*C*D
          LAI::matrixMatrixMultiply( C, D, RHS, beta, 1. );

          // Check norms
          EXPECT_NEAR( LAI::matrixNormInf( LHS ),
                       LAI::matrixNormInf( RHS ),
                       LAI::matrixNormInf( LHS ) * machinePrecision );
          EXPECT_NEAR( LAI::matrixNorm1( LHS ),
                       LAI::matrixNorm1( RHS ),
                       LAI::matrixNorm1( LHS ) * machinePrecision );
          EXPECT_NEAR( LAI::matrixNormFrobenius( LHS ),
                       LAI::matrixNormFrobenius( RHS ),
                       LAI::matrixNormFrobenius( LHS ) * machinePrecision );

        }
      }
    }
  }

  // Repeat the following step for vectors of increasing size:
  //
  // a. compute LHS = ( alfa*A*B^T + beta*C ) * D
  // b. compute RHS = alfa*A*(B^T*D) + beta*C*D
  // c. check that (norm(LHS - RHS, type) / norm(LHS, type)) < epsilon
  //    with type = \{Infinity-norm, 1-norm, Frobenius-norm, )

  MA = 10;
  MB = 10;
  ND = 10;
  K = 20;

  for( localIndex mA = 1 ; mA <= MA ; ++mA )
  {
    for( localIndex mB = 1 ; mB <= MB ; ++mB )
    {
      for( localIndex nD = 1 ; nD <= ND ; ++nD )
      {
        for( localIndex k = 1 ; k <= K ; ++k )
        {

          // Resize matrix operators
          A.resize( mA, k );
          B.resize( mB, k );
          C.resize( A.size( 0 ), B.size( 0 ) );
          D.resize( B.size( 0 ), nD );
          LHS.resize( A.size( 0 ), D.size( 1 ) );
          RHS.resize( A.size( 0 ), D.size( 1 ) );

          // Populate A, B, C, and D with uniformly distributed random
          // coefficients
          matrixRand( A,
                      distribution,
                      generator );
          matrixRand( B,
                      distribution,
                      generator );
          matrixRand( C,
                      distribution,
                      generator );
          matrixRand( D,
                      distribution,
                      generator );

          // Compute tmp1 = ( alfa*A*B^T + beta*C )
          array2d<real64> tmp1( C );
          LAI::matrixMatrixTMultiply( A, B, tmp1, alfa, beta );

          // Compute LHS = tmp*D
          LAI::matrixMatrixMultiply( tmp1, D, LHS );

          // Compute tmp2 = B^T*D
          array2d<real64> tmp2( B.size( 1 ), D.size( 1 ) );
          LAI::matrixTMatrixMultiply( B, D, tmp2 );

          // Compute RHS = alfa*A*tmp2
          LAI::matrixMatrixMultiply( A, tmp2, RHS, alfa );

          // Compute RHS = RHS + beta*C*D
          LAI::matrixMatrixMultiply( C, D, RHS, beta, 1. );

          // Check norms
          EXPECT_NEAR( LAI::matrixNormInf( LHS ),
                       LAI::matrixNormInf( RHS ),
                       LAI::matrixNormInf( LHS ) * machinePrecision );
          EXPECT_NEAR( LAI::matrixNorm1( LHS ),
                       LAI::matrixNorm1( RHS ),
                       LAI::matrixNorm1( LHS ) * machinePrecision );
          EXPECT_NEAR( LAI::matrixNormFrobenius( LHS ),
                       LAI::matrixNormFrobenius( RHS ),
                       LAI::matrixNormFrobenius( LHS ) * machinePrecision );

        }
      }
    }
  }

  // Repeat the following step for vectors of increasing size:
  //
  // a. compute LHS = ( alfa*A^T*B^T + beta*C ) * D
  // b. compute RHS = alfa*A^T*(B^T*D) + beta*C*D
  // c. check that (norm(LHS - RHS, type) / norm(LHS, type)) < epsilon
  //    with type = \{Infinity-norm, 1-norm, Frobenius-norm, )

  NA = 10;
  MB = 10;
  ND = 10;
  K = 20;

  for( localIndex nA = 1 ; nA <= NA ; ++nA )
  {
    for( localIndex mB = 1 ; mB <= MB ; ++mB )
    {
      for( localIndex nD = 1 ; nD <= ND ; ++nD )
      {
        for( localIndex k = 1 ; k <= K ; ++k )
        {

          // Resize matrix operators
          A.resize( k, nA );
          B.resize( mB, k );
          C.resize( A.size( 1 ), B.size( 0 ) );
          D.resize( B.size( 0 ), nD );
          LHS.resize( A.size( 1 ), D.size( 1 ) );
          RHS.resize( A.size( 1 ), D.size( 1 ) );

          // Populate A, B, C, and D with uniformly distributed random
          // coefficients
          matrixRand( A,
                      distribution,
                      generator );
          matrixRand( B,
                      distribution,
                      generator );
          matrixRand( C,
                      distribution,
                      generator );
          matrixRand( D,
                      distribution,
                      generator );

          // Compute tmp1 = ( alfa*A^T*B^T + beta*C )
          array2d<real64> tmp1( C );
          LAI::matrixTMatrixTMultiply( A, B, tmp1, alfa, beta );

          // Compute LHS = tmp*D
          LAI::matrixMatrixMultiply( tmp1, D, LHS );

          // Compute tmp2 = B^T*D
          array2d<real64> tmp2( B.size( 1 ), D.size( 1 ) );
          LAI::matrixTMatrixMultiply( B, D, tmp2 );

          // Compute RHS = alfa*A^T*tmp2
          LAI::matrixTMatrixMultiply( A, tmp2, RHS, alfa );

          // Compute RHS = RHS + beta*C*D
          LAI::matrixMatrixMultiply( C, D, RHS, beta, 1. );

          // Check norms
          EXPECT_NEAR( LAI::matrixNormInf( LHS ),
                       LAI::matrixNormInf( RHS ),
                       LAI::matrixNormInf( LHS ) * machinePrecision );
          EXPECT_NEAR( LAI::matrixNorm1( LHS ),
                       LAI::matrixNorm1( RHS ),
                       LAI::matrixNorm1( LHS ) * machinePrecision );
          EXPECT_NEAR( LAI::matrixNormFrobenius( LHS ),
                       LAI::matrixNormFrobenius( RHS ),
                       LAI::matrixNormFrobenius( LHS ) * machinePrecision );

        }
      }
    }
  }

}

template<typename LAI>
void testArray2dArray1dLA()
{

  std::default_random_engine generator;
  std::uniform_real_distribution<real64> distribution( 0.0, 1.0 );

  array2d<real64> A, yT, tmp;
  array1d<real64> x, y, yTT;

  // Repeat the following step for vectors of increasing size:
  //
  // a. compute y = A*x
  // b. compute compute yT = x^T * A^T
  // c. check that (norm(y, 2-norm) / norm(yT, Frobenius-norm)) < epsilon
  // d. create vector yTT = yT^T
  // e. compute beta = sqrt(yTT dot y);
  // f. check that ( (norm(y, norm-2) - sqrt(yTT_dot_y)) / norm(y, norm-2)) < epsilon

  localIndex MA = 10;
  localIndex NA = 10;
  real64 alfa, beta;

  for( localIndex mA = 1 ; mA <= MA ; ++mA )
  {
    for( localIndex nA = 1 ; nA <= NA ; ++nA )
    {
      // Resize matrices and vectors
      A.resize( mA, nA );
      x.resize( nA );
      y.resize( mA );

      // Populate matrix A and vector x with uniformly distributed random
      // coefficients
      matrixRand( A,
                  distribution,
                  generator );
      vectorRand( x,
                  distribution,
                  generator );

      // a.
      LAI::matrixVectorMultiply( A, x, y );

      // b.
      // --- construct tmp = x^T
      tmp.resize( 1, x.size() );
      for( int i = 0 ; i < x.size() ; ++i )
      {
        tmp( 0, i ) = x( i );
      }

      // -- compute yT = tmp * A^T
      yT.resize( 1, A.size( 0 ) );
      LAI::matrixMatrixTMultiply( tmp, A, yT );

      // c.
      alfa = LAI::vectorNorm2( y );
      beta = LAI::matrixNormFrobenius( yT );
      EXPECT_NEAR( alfa,
                   beta,
                   alfa * machinePrecision );

      // d.
      yTT.resize( yT.size( 1 ) );
      for( int i = 0 ; i < yT.size( 1 ) ; ++i )
      {
        yTT( i ) = yT( 0, i );
      }

      // e.
      beta = std::sqrt( LAI::vectorDot( y, yTT ) );

      // f. check that ( (norm(y, norm-2) == sqrt(yTT_dot_y))
      EXPECT_NEAR( alfa,
                   beta,
                   alfa * machinePrecision );
    }
  }
}

template<typename LAI>
void testArray2dInverseLA()
{

  array2d<real64> E;
  array2d<real64> Einv;
  array2d<real64> EinvXE;

  // Repeat the following step for matrices of increasing size:
  // a. Construct matrix E (1d discrete Laplacian)
  // b. Compute Einv = E^-1
  // c. Compute EinvXE = Einv*E
  // d. Check that det(EinvXE) = 1.

  real64 det;
  localIndex max_dim = 10;

  for( int order = 1 ; order <= max_dim ; ++order )
  {
    // a.
    E.resize( order, order );
    E = 0;
    for( int i = 0 ; i < E.size(0) ; ++i )
    {
      for( int j = 0 ; j < E.size(1) ; ++j )
      {
        if( i == j )
        {
          E( i, i ) = 2;
        }
        else if( std::abs( i - j ) == 1 )
        {
          E( i, j ) = -1;
        }
      }
    }

    // b.
    Einv.resize( E.size(0), E.size(1) );
    LAI::matrixInverse(E, Einv, det );

    // c.
    EinvXE.resize( E.size(0), E.size(1) );
    LAI::matrixMatrixMultiply( Einv, E, EinvXE );

    // d.
    EXPECT_NEAR( LAI::determinant(EinvXE),
                 1.,
                 machinePrecision );
  }

}

// -----------------------------------------------------------------------------


// --- Aggregate
TEST(BlasLapackLA,testArray1dLA)
{
  testArray1dLA<BlasLapackLA>();
}

TEST(BlasLapackLA,testArray2dLA)
{
  testArray2dLA<BlasLapackLA>();
}

TEST(BlasLapackLA,testArray2dArray1dLA)
{
  testArray2dArray1dLA<BlasLapackLA>();
}

TEST(BlasLapackLA,testArray2dInverseLA)
{
  testArray2dInverseLA<BlasLapackLA>();
}

#ifdef __clang__
#pragma clang diagnostic pop
#endif

// Generate a random GEOSX array1d
void vectorRand( array1d<real64> & x,
                 std::uniform_real_distribution<real64> & distribution,
                 std::default_random_engine & generator )
{
  for( int i = 0 ; i < x.size() ; ++i )
  {
    x( i ) = distribution( generator );
  }
  return;
}

// Generate a random GEOSX array2d
void matrixRand( array2d<real64> & A,
                 std::uniform_real_distribution<real64> & distribution,
                 std::default_random_engine & generator )
{
  for( int i = 0 ; i < A.size( 0 ) ; ++i )
  {
    for( int j = 0 ; j < A.size( 1 ) ; ++j )
    {
      A( i, j ) = distribution( generator );
    }
  }
  return;
}

