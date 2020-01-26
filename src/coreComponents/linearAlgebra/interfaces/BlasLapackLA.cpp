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
 * @file BlasLapack.cpp
 */

// Include the corresponding header file.
#include "BlasLapackLA.hpp"

#include <random>

// BLAS and LAPACK function declaration
#include "BlasLapackFunctions.h"

// Put everything under the geosx namespace.
namespace geosx
{

// Random device and random number generator seed integer array used
// to populate a vector/matrix with random coefficients
static std::random_device rd;
static std::mt19937 gen(rd());
static std::uniform_int_distribution<int> dis(0, 4095);
static std::uniform_int_distribution<int> disOdd(0, 2047);
static int ISEED[] = {dis(gen), dis(gen), dis(gen), disOdd(gen)*2 + 1};
  
real64 BlasLapackLA::vectorNorm1( array1d<real64> const & X )
{
  int const INCX = 1;
  int const N = integer_conversion<int>( X.size() );
  return GEOSX_dasum( &N, X.data(), &INCX);
}

real64 BlasLapackLA::vectorNorm2( array1d<real64> const & X )
{
  int const INCX = 1;
  int const N = integer_conversion<int>( X.size() );
  return GEOSX_dnrm2( &N, X.data(), &INCX);
}

real64 BlasLapackLA::vectorNormInf( array1d<real64> const & X )
                                    {
  int const INCX = 1;
  int const N = integer_conversion<int>( X.size() );
  int ind = GEOSX_idamax( &N, X.data(), &INCX);
  ind -= 1; // Fortran convention, subtract 1
  return std::abs( X(ind) );
}

real64 BlasLapackLA::determinant( array2d<real64> const & A )
                                  {
  // --- check that matrix is square
  GEOSX_ASSERT_MSG( A.size( 0 ) == A.size( 1 ) &&
                   A.size( 0 ) > 0,
                   "Matrix must be square with order greater than zero" );

  real64 det;
  switch( A.size( 0 ) )
  {
    case 1:
    {
      det =  A( 0, 0 );
      break;
    }
    case 2:
    {
      det = A( 0, 0 ) * A( 1, 1 ) - A( 0, 1 ) * A( 1, 0 );
      break;
    }
    case 3:
    {
      det =
      A( 0, 0 ) * ( A( 1, 1 ) * A( 2, 2 ) - A( 1, 2 ) * A( 2, 1 ) ) +
      A( 0, 1 ) * ( A( 1, 2 ) * A( 2, 0 ) - A( 1, 0 ) * A( 2, 2 ) ) +
      A( 0, 2 ) * ( A( 1, 0 ) * A( 2, 1 ) - A( 1, 1 ) * A( 2, 0 ) );
      break;
    }
    case 4:
    {
      det =
      A( 0, 0 ) * ( A( 1, 1 ) * ( A( 2, 2 ) * A( 3, 3 ) - A( 3, 2 ) * A( 2, 3 ) ) -
                    A( 1, 2 ) * ( A( 2, 1 ) * A( 3, 3 ) - A( 3, 1 ) * A( 2, 3 ) ) +
                    A( 1, 3 ) * ( A( 2, 1 ) * A( 3, 2 ) - A( 3, 1 ) * A( 2, 2 ) )
                  ) -
      A( 0, 1 ) * ( A( 1, 0 ) * ( A( 2, 2 ) * A( 3, 3 ) - A( 3, 2 ) * A( 2, 3 ) ) -
                    A( 1, 2 ) * ( A( 2, 0 ) * A( 3, 3 ) - A( 3, 0 ) * A( 2, 3 ) ) +
                    A( 1, 3 ) * ( A( 2, 0 ) * A( 3, 2 ) - A( 3, 0 ) * A( 2, 2 ) )
                  ) +
      A( 0, 2 ) * ( A( 1, 0 ) * ( A( 2, 1 ) * A( 3, 3 ) - A( 3, 1 ) * A( 2, 3 ) ) -
                    A( 1, 1 ) * ( A( 2, 0 ) * A( 3, 3 ) - A( 3, 0 ) * A( 2, 3 ) ) +
                    A( 1, 3 ) * ( A( 2, 0 ) * A( 3, 1 ) - A( 3, 0 ) * A( 2, 1 ) )
                  ) -
      A( 0, 3 ) * ( A( 1, 0 ) * ( A( 2, 1 ) * A( 3, 2 ) - A( 3, 1 ) * A( 2, 2 ) ) -
                    A( 1, 1 ) * ( A( 2, 0 ) * A( 3, 2 ) - A( 3, 0 ) * A( 2, 2 ) ) +
                    A( 1, 2 ) * ( A( 2, 0 ) * A( 3, 1 ) - A( 3, 0 ) * A( 2, 1 ) )
                  );
      break;
    }
    default:
    {
      // Compute the determinant via LU factorization
      int const NN = integer_conversion<int>( A.size(0) );
      int INFO;
      array1d<int> IPIV( NN );
      array2d<double> LUFactor( A );

      // We compute the LU factors for the transpose matrix, i.e. choosing the
      // LAPACK_COL_MAJOR ordering, to avoid transposition/copy requires for
      // LAPACK_ROW_MAJOR ordering.
      GEOSX_dgetrf( &NN, &NN,  LUFactor.data(), &NN, IPIV.data(), &INFO);

      GEOSX_ASSERT_MSG( INFO == 0, "LAPACK dgetrf error code: " << INFO );

      det = 1.0;
      for( int i = 0 ; i < NN ; ++i )
      {
        if( IPIV[i] != i + 1 ) //IPIV is based on Fortran convention (counting from 1)
        {
          det *= -LUFactor( i, i);
        }
        else
        {
          det *= LUFactor( i, i);
        }
      }

      break;
    }
  }

  return det;

}

real64 BlasLapackLA::matrixNormInf( array2d<real64> const & A )
{
  // Computed as one-norm of the transpose matrix
  char const NORM = '1';
  int const  M = integer_conversion<int>( A.size(0) );
  int const  N = integer_conversion<int>( A.size(1) );
  return GEOSX_dlange(&NORM,  &N, &M, A.data(), &N, nullptr);

}

real64 BlasLapackLA::matrixNorm1( array2d<real64> const & A )
{
  // Computed as infinity-norm of the transpose matrix
  char const NORM = 'I';
  int const M = integer_conversion<int>( A.size(0) );
  int const N = integer_conversion<int>( A.size(1) );
  array1d<double> WORK(N);
  return  GEOSX_dlange( &NORM,
                   &N,
                   &M,
                   A.data(),
                   &N,
                   WORK.data());
}

real64 BlasLapackLA::matrixNormFrobenius( array2d<real64> const & A )
{
  // Computed using the transpose matrix
  char const NORM = 'F';
  int const M = integer_conversion<int>( A.size(0) );
  int const N = integer_conversion<int>( A.size(1) );
  return GEOSX_dlange(&NORM, &N, &M, A.data(), &N, nullptr);

}

void BlasLapackLA::vectorVectorAdd( array1d<real64> const & X,
                                    array1d<real64> & Y,
                                    real64 const alpha )
{

  GEOSX_ASSERT_MSG( X.size() == Y.size(),
                   "Vector dimensions not compatible for sum" );

  int const INCX = 1;
  int const INCY = 1;
  int const N = static_cast<int>( X.size() );
  GEOSX_daxpy( &N, &alpha, X.data(), &INCX, Y.data(), &INCY );

  return;
}

void BlasLapackLA::matrixMatrixAdd( array2d<real64> const & A,
                                    array2d<real64> & B,
                                    real64 const alpha )
{

  GEOSX_ASSERT_MSG( A.size( 0 ) == B.size( 0 ) &&
                   A.size( 1 ) == B.size( 1 ),
                   "Matrix dimensions not compatible for sum" );

  int const INCX = 1;
  int const INCY = 1;
  int const N = static_cast<int>( A.size() );
  GEOSX_daxpy( &N, &alpha, A.data(), &INCX, B.data(), &INCY );

  return;
}

void BlasLapackLA::vectorScale( real64 alpha,
                                array1d<real64> & X )
{

  int const INCX = 1;
  int const N = static_cast<int>( X.size() );
  GEOSX_dscal( &N, &alpha, X.data(), &INCX );

  return;
}

void BlasLapackLA::matrixScale( real64 alpha,
                                array2d<real64> & A )
{

  int const INCX = 1;
  int const N = static_cast<int>( A.size() );
  GEOSX_dscal( &N, &alpha, A.data(), &INCX );

  return;
}

real64 BlasLapackLA::vectorDot( array1d<real64> const & X,
                                array1d<real64> const & Y )
{
  GEOSX_ASSERT_MSG( X.size() == Y.size(),
                   "Vector dimensions not compatible for dot product" );

  int const INCX = 1;
  int const INCY = 1;
  int const N = static_cast<int>( X.size() );
  return GEOSX_ddot( &N, X.data(), &INCX, Y.data(), &INCY );

}

void BlasLapackLA::matrixVectorMultiply( array2d<real64> const & A,
                                         array1d<real64> const & X,
                                         array1d<real64> & Y,
                                         real64 const alpha,
                                         real64 const beta )
{
  GEOSX_ASSERT_MSG( A.size( 1 ) == X.size() &&
                   A.size( 0 ) == Y.size(),
                   "Matrix, source vector and destination vector not compatible" );

  int const M = integer_conversion<int>( A.size( 0 ) );
  int const N = 1;
  int const K = integer_conversion<int>( A.size( 1 ) );

  // A*X = Y is computed as X^T * A^T = Y^T, i.e. accessing the transpose
  // matrix using a column-major layout
  char const TRANS1 = 'N';
  char const TRANS2 = 'N';

  GEOSX_dgemm(&TRANS1, &TRANS2, &N, &M, &K, &alpha, X.data(),&N, A.data(), &K, &beta, Y.data(), &N);

  return;
}

void BlasLapackLA::matrixTVectorMultiply( array2d<real64> const & A,
                                          array1d<real64> const & X,
                                          array1d<real64> & Y,
                                          real64 const alpha,
                                          real64 const beta )
{
  GEOSX_ASSERT_MSG( A.size( 0 ) == X.size() &&
                       A.size( 1 ) == Y.size(),
                   "Matrix, source vector and destination vector not compatible" );

  int const M = integer_conversion<int>( A.size( 1 ) );
  int const N = 1;
  int const K = integer_conversion<int>( A.size( 0 ) );

  // A^T*X = Y is computed as X^T * A = Y^T, i.e. accessing the transpose
  // matrix using a column-major layout
  char const TRANS1 = 'N';
  char const TRANS2 = 'T';

  GEOSX_dgemm(&TRANS1, &TRANS2, &N, &M, &K, &alpha, X.data(), &N, A.data(), &M, &beta, Y.data(), &N);

  return;
}

void BlasLapackLA::matrixMatrixMultiply( array2d<real64> const & A,
                                         array2d<real64> const & B,
                                         array2d<real64> & C,
                                         real64 const alpha,
                                         real64 const beta )
{

  GEOSX_ASSERT_MSG( C.size( 0 ) == A.size( 0 ) &&
                   C.size( 1 ) == B.size( 1 ) &&
                   A.size( 1 ) == B.size( 0 ),
                   "Matrix dimensions not compatible for product" );

  int const M = integer_conversion<int>( A.size( 0 ) );
  int const N = integer_conversion<int>( B.size( 1 ) );
  int const K = integer_conversion<int>( A.size( 1 ) );

  // A*B = C is computed as B^T * A^T = C^T, i.e. accessing the transpose
  // matrices using a column-major layout
  char const TRANS1 = 'N';
  char const TRANS2 = 'N';

  GEOSX_dgemm(&TRANS1, &TRANS2, &N, &M, &K, &alpha, B.data(), &N, A.data(), &K, &beta, C.data(), &N);

  return;
}

void BlasLapackLA::matrixTMatrixMultiply( array2d<real64> const & A,
                                          array2d<real64> const & B,
                                          array2d<real64> & C,
                                          real64 const alpha,
                                          real64 const beta )
{

  GEOSX_ASSERT_MSG( C.size( 0 ) == A.size( 1 ) &&
                   C.size( 1 ) == B.size( 1 ) &&
                   A.size( 0 ) == B.size( 0 ),
                   "Matrix dimensions not compatible for product" );

  int const M = integer_conversion<int>( A.size( 1 ) );
  int const N = integer_conversion<int>( B.size( 1 ) );
  int const K = integer_conversion<int>( A.size( 0 ) );

  // A^T*B = C is computed as B^T * A = C^T, i.e. accessing the transpose
  // matrices using a column-major layout

  char const TRANS1 = 'N';
  char const TRANS2 = 'T';

  GEOSX_dgemm(&TRANS1, &TRANS2, &N, &M, &K, &alpha, B.data(), &N, A.data(), &M, &beta, C.data(), &N);

  return;
}

void BlasLapackLA::matrixMatrixTMultiply( array2d<real64> const & A,
                                          array2d<real64> const & B,
                                          array2d<real64> & C,
                                          real64 const alpha,
                                          real64 const beta )
{

  GEOSX_ASSERT_MSG( C.size( 0 ) == A.size( 0 ) &&
                   C.size( 1 ) == B.size( 0 ) &&
                   A.size( 1 ) == B.size( 1 ),
                   "Matrix dimensions not compatible for product" );

  int const M = integer_conversion<int>( A.size( 0 ) );
  int const N = integer_conversion<int>( B.size( 0 ) );
  int const K = integer_conversion<int>( A.size( 1 ) );

  // A*B^T = C is computed as B * A^T = C^T, i.e. accessing the transpose
  // matrices using a column-major layout

  char const TRANS1 = 'T';
  char const TRANS2 = 'N';

  GEOSX_dgemm(&TRANS1, &TRANS2, &N, &M, &K, &alpha, B.data(), &K, A.data(), &K, &beta, C.data(), &N);

  return;
}

void BlasLapackLA::matrixTMatrixTMultiply( array2d<real64> const & A,
                                           array2d<real64> const & B,
                                           array2d<real64> & C,
                                           real64 const alpha,
                                           real64 const beta )
{

  GEOSX_ASSERT_MSG( C.size( 0 ) == A.size( 1 ) &&
                       C.size( 1 ) == B.size( 0 ) &&
                       A.size( 0 ) == B.size( 1 ),
                   "Matrix dimensions not compatible for product" );

  int const M = integer_conversion<int>( A.size( 1 ) );
  int const N = integer_conversion<int>( B.size( 0 ) );
  int const K = integer_conversion<int>( A.size( 0 ) );

  // A^T*B^T = C is computed as B * A = C^T, i.e. accessing the transpose
  // matrices using a column-major layout

  char const TRANS1 = 'T';
  char const TRANS2 = 'T';

  GEOSX_dgemm(&TRANS1, &TRANS2, &N, &M, &K, &alpha, B.data(), &K, A.data(), &M, &beta, C.data(), &N);

  return;
}

void BlasLapackLA::matrixInverse( array2d<real64> const & A,
                                  array2d<real64> & Ainv )
{
  real64 detA;
  matrixInverse( A, Ainv, detA );
}

void BlasLapackLA::matrixInverse( array2d<real64> const & A,
                                  array2d<real64> & Ainv,
                                  real64 & detA )
{
  // --- Check that source matrix is square
  int const NN = integer_conversion<int>(A.size( 0 ));
  GEOSX_ASSERT_MSG( NN > 0 &&
                   NN == A.size( 1 ),
                   "Matrix must be square" );

  // --- Check that inverse matrix has appropriate dimension
  GEOSX_ASSERT_MSG( Ainv.size( 0 ) == NN &&
                   Ainv.size( 1 ) == NN,
                   "Inverse matrix has wrong dimensions" );

  // --- Check if matrix is singular by computing the determinant
  //     note: if order greater than 3 we compute the determinant by
  //           first constructing the LU factors, later reused for calculating
  //           the inverse.
  array1d<int> IPIV;
  array1d<double> INV_WORK;

  if (NN <= 3)
  {
    detA = determinant(A);
  }
  else
  {
    // Copy A in Ainv
    matrixCopy(A, Ainv);

    // Declare workspace for permutations and scratch array
    IPIV.resize(NN);
    INV_WORK.resize(NN);

    // Compute determinant (not done calling directly the function determinant
    // (avoid computing twice LUFactors, currently stored in Ainv, needed for
    // computing the inverse). Again we compute the LU factors for the
    // transpose matrix, i.e. choosing the LAPACK_COL_MAJOR ordering, to
    // avoid transposition/copy requires for LAPACK_ROW_MAJOR ordering.
    int INFO;
    GEOSX_dgetrf( &NN, &NN, Ainv.data(), &NN, IPIV.data(), &INFO);

    GEOSX_ASSERT_MSG( INFO == 0, "LAPACK dgetrf error code: " << INFO );

    detA = 1.0;
    for( int i = 0 ; i < NN ; ++i )
    {
      if( IPIV[i] != i + 1 ) //IPIV is based on Fortran convention (counting from 1)
      {
        detA *= -Ainv(i,i);
      }
      else
      {
        detA *= Ainv(i,i);
      }
    }
  }

  // Check if matrix is singular
  GEOSX_ASSERT_MSG( std::abs(detA) >
                   std::numeric_limits<real64>::epsilon() *
                   matrixNormFrobenius(A),
                   "Matrix is singular" );
  real64 oneOverDetA = 1. / detA;

  // --- Compute inverse
  switch( NN )
    {
      case 1:
      {
        Ainv( 0, 0 ) = oneOverDetA;
        break;
      }

      case 2:
      {
        Ainv( 0, 0 ) =  A( 1, 1 ) * oneOverDetA;
        Ainv( 0, 1 ) = -A( 0, 1 ) * oneOverDetA;
        Ainv( 1, 0 ) = -A( 1, 0 ) * oneOverDetA;
        Ainv( 1, 1 ) =  A( 0, 0 ) * oneOverDetA;
        break;
      }

      case 3:
      {
        Ainv( 0, 0 ) = ( A( 1, 1 ) * A( 2, 2 ) -
                         A( 1, 2 ) * A( 2, 1 ) ) * oneOverDetA;
        Ainv( 0, 1 ) = ( A( 0, 2 ) * A( 2, 1 ) -
                         A( 0, 1 ) * A( 2, 2 ) ) * oneOverDetA;
        Ainv( 0, 2 ) = ( A( 0, 1 ) * A( 1, 2 ) -
                         A( 0, 2 ) * A( 1, 1 ) ) * oneOverDetA;
        Ainv( 1, 0 ) = ( A( 1, 2 ) * A( 2, 0 ) -
                         A( 1, 0 ) * A( 2, 2 ) ) * oneOverDetA;
        Ainv( 1, 1 ) = ( A( 0, 0 ) * A( 2, 2 ) -
                         A( 0, 2 ) * A( 2, 0 ) ) * oneOverDetA;
        Ainv( 1, 2 ) = ( A( 0, 2 ) * A( 1, 0 ) -
                         A( 0, 0 ) * A( 1, 2 ) ) * oneOverDetA;
        Ainv( 2, 0 ) = ( A( 1, 0 ) * A( 2, 1 ) -
                         A( 1, 1 ) * A( 2, 0 ) ) * oneOverDetA;
        Ainv( 2, 1 ) = ( A( 0, 1 ) * A( 2, 0 ) -
                         A( 0, 0 ) * A( 2, 1 ) ) * oneOverDetA;
        Ainv( 2, 2 ) = ( A( 0, 0 ) * A( 1, 1 ) -
                         A( 0, 1 ) * A( 1, 0 ) ) * oneOverDetA;
        break;
      }
      default:
    {
    // Invert (LAPACK function DGETRI). The LU factors computed for the
    // transpose matrix stored in Ainv are used.
    int INFO;
    GEOSX_dgetri( &NN, Ainv.data(), &NN, IPIV.data(), INV_WORK.data(), &NN, &INFO);

    GEOSX_ASSERT_MSG( INFO == 0, "LAPACK dgetri error code: " << INFO );

    break;
    }
  }
  return;
}

void BlasLapackLA::vectorCopy( array1d<real64> const & X,
                               array1d<real64> & Y )
{
  GEOSX_ASSERT_MSG( X.size() == Y.size(),
                   "Vector dimensions not compatible for copying" );

  int const INCX = 1;
  int const INCY = 1;
  int const N = static_cast<int>( X.size() );
  GEOSX_dcopy( &N, X.data(), &INCX, Y.data(), &INCY );

  return;
}

void BlasLapackLA::matrixCopy( array2d<real64> const & A,
                               array2d<real64> & B )
{
  GEOSX_ASSERT_MSG( A.size(0) == B.size(0) &&
                   A.size(1) == B.size(1),
                   "Matrix dimensions not compatible for copying" );

  int const INCX = 1;
  int const INCY = 1;
  int const N = static_cast<int>( A.size() );
  GEOSX_dcopy( &N, A.data(), &INCX, B.data(), &INCY );

  return;
}

void BlasLapackLA::setRandomNumberGeneratorSeed( array1d<int> const & seed)
{
  // Error checking
  GEOSX_ASSERT_MSG( seed.size() >= 4,
                   "Seed array must have size at least four");

  GEOSX_ASSERT_MSG( 0 <= seed(0) && seed(0) <= 4095 &&
                   0 <= seed(1) && seed(1) <= 4095 &&
                   0 <= seed(2) && seed(2) <= 4095 &&
                   0 <= seed(3) && seed(3) <= 4095,
                  "Seed array integer entries must be in interval [0,4095]");

  GEOSX_ASSERT_MSG( seed(3) % 2 > 0,
                   "Seed array 4th element must be odd");

  for (int i = 0; i < 4; ++i)
  {
    ISEED[i] = seed[i];
  }
}

void BlasLapackLA::getRandomNumberGeneratorSeed( array1d<int> & seed)
{
  // Error checking
  GEOSX_ASSERT_MSG( seed.size() >= 4,
                   "Seed array must have size at least four");

  for (int i = 0; i < 4; ++i)
  {
    seed[i] = ISEED[i];
  }
}

void BlasLapackLA::vectorRand( array1d<real64> & X,
                               RandomNumberDistribution const & idist)
{

  int IDIST = static_cast<int>(idist);
  int const N = static_cast<int>( X.size() );

  GEOSX_ASSERT_MSG( N > 0,
                   "The vector cannot be empty");

  GEOSX_dlarnv( &IDIST, ISEED, &N, X.data());

  return;
}

void BlasLapackLA::matrixRand( array2d<real64> & A,
                               RandomNumberDistribution const & idist)
{

  int IDIST = static_cast<int>(idist);
  int const NN = static_cast<int>( A.size() );

  GEOSX_ASSERT_MSG( NN > 0,
                   "The matrix cannot be empty");

  GEOSX_dlarnv( &IDIST, ISEED, &NN, A.data());

  return;
}

void BlasLapackLA::matrixSVD( array2d<real64> const & A,
		              array2d<real64> & U,
			      array1d<real64> & S,
			      array2d<real64> & VT )
{
  int const minDim = (A.size(0) < A.size(1))
                   ? static_cast<int>( A.size(0) )
                   : static_cast<int>( A.size(1) );

  GEOSX_ASSERT_MSG( A.size(0) == U.size(0) && minDim == U.size(1),
                   "The matrices A and U have an incompatible size");

  GEOSX_ASSERT_MSG( minDim == VT.size(0) && A.size(1) == VT.size(1),
                   "The matrices A and V have an incompatible size");

  GEOSX_ASSERT_MSG( S.size() == minDim,
                   "The matrix A and vector S have an incompatible size");

  
  // note: this is *not* an optimized implementation
  // this is just a temporary solution to have something that works

  array2d<real64> AT( A.size(1), A.size(0) );
  array2d<real64> UT( U.size(1), U.size(0) );
  array2d<real64> V( VT.size(1), VT.size(0) ); 

  // fill the transpose of A
  // this is done because Lapack expects column-major
  for (int i = 0; i < A.size(0); ++i)
  {
    for (int j = 0; j < A.size(1); ++j)
    {
      AT(j,i) = A(i,j); 
    }
  }

  // define the arguments of dgesvd
  int M     = static_cast<int>( A.size(0) );
  int N     = static_cast<int>( A.size(1) ); 
  int LDA   = M;
  int LDU   = M;
  int LDVT  = minDim;
  int LWORK = 0;
  int INFO  = 0;
  double WKOPT = 0.0;

  // 1) query and allocate the optimal workspace
  LWORK = -1;
  GEOSX_dgesvd( "S", "S",
		&M, &N, AT.data(), &LDA,
		S.data(), UT.data(), &LDU, V.data(), &LDVT,
		&WKOPT, &LWORK, &INFO );

  LWORK = static_cast<int>( WKOPT );
  array1d<real64> WORK( LWORK );
  
  // 2) compute svd 
  GEOSX_dgesvd( "S", "S",
		&M, &N, AT.data(), &LDA,
		S.data(), UT.data(), &LDU, V.data(), &LDVT,
		WORK.data(), &LWORK, &INFO );
  
  GEOSX_ASSERT_MSG( INFO == 0, 
                    "The algorithm computing SVD failed to converge." );
  
  // fill the transpose of U and VT
  // this is done because Lapack uses column-major
  for (int i = 0; i < UT.size(0); ++i)
  {
    for (int j = 0; j < UT.size(1); ++j)
    {
      U(j,i) = UT(i,j); 
    }
  }
  for (int i = 0; i < V.size(0); ++i)
  {
    for (int j = 0; j < V.size(1); ++j)
    {
      VT(j,i) = V(i,j); 
    }
  }
}
  
} // end geosx namespace
