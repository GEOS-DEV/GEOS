/*
 * LapackLA.cpp
 *
 *  Created on: Apr 8, 2019
 *      Author: castelletto1
 */

// Include the corresponding header file.
#include "BlasLapackLA.hpp"

// Put everything under the geosx namespace.
namespace geosx
{

//-------------------------------------------------------Mathematical methods---

real64 BlasLapackLA::vectorNorm1( array1d<real64> const & X ) const
                                  {
  return cblas_dasum( integer_conversion<int>( X.size() ),
                      X.data(),
                      1 );
}

real64 BlasLapackLA::vectorNorm2( array1d<real64> const & X ) const
                                  {
  return cblas_dnrm2( integer_conversion<int>( X.size() ),
                      X.data(),
                      1 );
}

real64 BlasLapackLA::vectorNormInf( array1d<real64> const & X ) const
                                    {
  int ind = cblas_idamax( integer_conversion<int>( X.size() ),
                          X.data(),
                          1 );
  return std::abs( X[ind] );
}

real64 BlasLapackLA::determinant( array2d<real64> const & A ) const
                                  {
  // --- check that matrix is square
  GEOS_ASSERT_MSG( A.size( 0 ) == A.size( 1 ) &&
                   A.size( 0 ) > 0,
                   "Matrix must be square with order greater than zero" );

  switch( A.size( 0 ) )
  {
    case 1:
      {
      return A( 0, 0 );
    }
    case 2:
      {
      return A( 0, 0 ) * A( 1, 1 ) - A( 0, 1 ) * A( 1, 0 );
    }
    case 3:
      {
      return
      A( 0, 0 ) * ( A( 1, 1 ) * A( 2, 2 ) - A( 1, 2 ) * A( 2, 1 ) ) +
      A( 0, 1 ) * ( A( 1, 2 ) * A( 2, 0 ) - A( 1, 0 ) * A( 2, 2 ) ) +
      A( 0, 2 ) * ( A( 1, 0 ) * A( 2, 1 ) - A( 1, 1 ) * A( 2, 0 ) );
    }
    case 4:
      {
      return
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
    }
    default:

      // Compute the determinant via LU factorization
      lapack_int NN = integer_conversion<lapack_int>( A.size(0) );
      array1d<lapack_int> IPIV( NN );
      array2d<double> LUfactor( A );

      // We compute the LU factors for the transpose matrix, i.e. choosing the
      // LAPACK_COL_MAJOR ordering, to avoid transposition/copy requires for
      // LAPACK_ROW_MAJOR ordering.
      lapack_int INFO =  LAPACKE_dgetrf( LAPACK_COL_MAJOR,
                                         NN,
                                         NN,
                                         LUfactor.data(),
                                         NN,
                                         IPIV.data() );

      GEOS_ASSERT_MSG( INFO == 0, "LAPACKE_dgetrf error code: " << INFO );

      real64 det = 1.0;
      for( int i = 0 ; i < NN ; ++i )
      {
        if( IPIV[i] != i + 1 ) //IPIV is based on Fortran convention (counting from 1)
        {
          det *= -LUfactor( i, i);
        }
        else
        {
          det *= LUfactor( i, i);
        }
      }

      return det;

  }

}

real64 BlasLapackLA::matrixNormInf( array2d<real64> const & A ) const
{
  // Computed as one-norm of the transpose matrix, i.e. assuming
  // column major ordering
  return LAPACKE_dlange( LAPACK_COL_MAJOR,
                         '1',
                         integer_conversion<lapack_int>( A.size( 1 ) ),
                         integer_conversion<lapack_int>( A.size( 0 ) ),
                         A.data(),
                         integer_conversion<lapack_int>( A.size( 1 ) ) );

}

real64 BlasLapackLA::matrixNorm1( array2d<real64> const & A ) const
{
  // Computed as infinity-norm of the transpose matrix, i.e. assuming
  // column major ordering
  return LAPACKE_dlange( LAPACK_COL_MAJOR,
                         'I',
                         integer_conversion<lapack_int>( A.size( 1 ) ),
                         integer_conversion<lapack_int>( A.size( 0 ) ),
                         A.data(),
                         integer_conversion<lapack_int>( A.size( 1 ) ) );

}

real64 BlasLapackLA::matrixNormFrobenius( array2d<real64> const & A ) const
{
  // Computed using the transpose matrix, i.e. assuming column major ordering
  return LAPACKE_dlange( LAPACK_COL_MAJOR,
                         'F',
                         integer_conversion<lapack_int>( A.size( 1 ) ),
                         integer_conversion<lapack_int>( A.size( 0 ) ),
                         A.data(),
                         integer_conversion<lapack_int>( A.size( 1 ) ) );

}

void BlasLapackLA::vectorVectorAdd( array1d<real64> const & X,
                                    array1d<real64> & Y,
                                    real64 const alpha )
{

  GEOS_ASSERT_MSG( X.size() == Y.size(),
                   "Vector dimensions not compatible for sum" );

  cblas_daxpy( integer_conversion<int>( X.size() ),
               alpha,
               X.data(),
               1,
               Y.data(),
               1 );
  return;
}

void BlasLapackLA::matrixMatrixAdd( array2d<real64> const & A,
                                    array2d<real64> & B,
                                    real64 const alpha )
{

  GEOS_ASSERT_MSG( A.size( 0 ) == B.size( 0 ) &&
                       A.size( 1 ) == B.size( 1 ),
                   "Matrix dimensions not compatible for sum" );

  cblas_daxpy( static_cast<int>( A.size( 0 ) * A.size( 1 ) ),
               alpha,
               A.data(),
               1,
               B.data(),
               1 );
  return;
}

void BlasLapackLA::vectorScale( array1d<real64> & X,
                                real64 alpha )
{

  cblas_dscal( integer_conversion<int>( X.size() ),
               alpha,
               X.data(),
               1 );
  return;
}

void BlasLapackLA::matrixScale( array2d<real64> & A,
                                real64 alpha )
{
  cblas_dscal( integer_conversion<int>( A.size( 0 ) * A.size( 1 ) ),
               alpha,
               A.data(),
               1 );
  return;
}

real64 BlasLapackLA::vectorDot( array1d<real64> const & X,
                                array1d<real64> const & Y )
{
  GEOS_ASSERT_MSG( X.size() == Y.size(),
                   "Vector dimensions not compatible for dot product" );

  return cblas_ddot( integer_conversion<int>( X.size() ),
                     X.data(),
                     1,
                     Y.data(),
                     1 );
}

void BlasLapackLA::matrixVectorMultiply( array2d<real64> const & A,
                                         array1d<real64> const & X,
                                         array1d<real64> & Y,
                                         real64 const alpha,
                                         real64 const beta )
{
  GEOS_ASSERT_MSG( A.size( 1 ) == X.size() &&
                   A.size( 0 ) == Y.size(),
                   "Matrix, source vector and destination vector not compatible" );

//  int M = integer_conversion<int>( A.size( 0 ) );
//  int N = integer_conversion<int>( A.size( 1 ) );
//
//  cblas_dgemv( CblasRowMajor,
//               CblasNoTrans,
//               M,
//               N,
//               alpha,
//               A.data(),
//               N,
//               X.data(),
//               1,
//               beta,
//               Y.data(),
//               1 );

  int M = integer_conversion<int>( A.size( 0 ) );
  int N = 1;
  int K = integer_conversion<int>( A.size( 1 ) );

  // CblasColMajor used for better performance. A*X = Y is computed
  // as X^T * A^T = Y^T
  cblas_dgemm( CblasColMajor,
               CblasNoTrans,
               CblasNoTrans,
               N,
               M,
               K,
               alpha,
               X.data(),
               N,
               A.data(),
               K,
               beta,
               Y.data(),
               N );

  return;
}

void BlasLapackLA::matrixTVectorMultiply( array2d<real64> const & A,
                                          array1d<real64> const & X,
                                          array1d<real64> & Y,
                                          real64 const alpha,
                                          real64 const beta )
{
  GEOS_ASSERT_MSG( A.size( 0 ) == X.size() &&
                       A.size( 1 ) == Y.size(),
                   "Matrix, source vector and destination vector not compatible" );

//  int M = integer_conversion<int>( A.size( 0 ) );
//  int N = integer_conversion<int>( A.size( 1 ) );
//
//  cblas_dgemv( CblasRowMajor,
//               CblasTrans,
//               M,
//               N,
//               alpha,
//               A.data(),
//               N,
//               X.data(),
//               1,
//               beta,
//               Y.data(),
//               1 );

  int M = integer_conversion<int>( A.size( 1 ) );
  int N = 1;
  int K = integer_conversion<int>( A.size( 0 ) );

  // CblasColMajor used for better performance. A^T*X = Y is computed
  // as X^T * A = Y^T
  cblas_dgemm( CblasColMajor,
               CblasNoTrans,
               CblasTrans,
               N,
               M,
               K,
               alpha,
               X.data(),
               N,
               A.data(),
               M,
               beta,
               Y.data(),
               N );

  return;
}

void BlasLapackLA::matrixMatrixMultiply( array2d<real64> const & A,
                                         array2d<real64> const & B,
                                         array2d<real64> & C,
                                         real64 const alpha,
                                         real64 const beta )
{

  GEOS_ASSERT_MSG( C.size( 0 ) == A.size( 0 ) &&
                       C.size( 1 ) == B.size( 1 ) &&
                       A.size( 1 ) == B.size( 0 ),
                   "Matrix dimensions not compatible for product" );

  int M = integer_conversion<int>( A.size( 0 ) );
  int N = integer_conversion<int>( B.size( 1 ) );
  int K = integer_conversion<int>( A.size( 1 ) );

//  cblas_dgemm( CblasRowMajor,
//               CblasNoTrans,
//               CblasNoTrans,
//               M,
//               N,
//               K,
//               alpha,
//               A.data(),
//               K,
//               B.data(),
//               N,
//               beta,
//               C.data(),
//               N );
  // CblasColMajor used for better performance. A*B = C is computed
  // as B^T * A^T = C^T
  cblas_dgemm( CblasColMajor,
               CblasNoTrans,
               CblasNoTrans,
               N,
               M,
               K,
               alpha,
               B.data(),
               N,
               A.data(),
               K,
               beta,
               C.data(),
               N );
  return;
}

void BlasLapackLA::matrixTMatrixMultiply( array2d<real64> const & A,
                                          array2d<real64> const & B,
                                          array2d<real64> & C,
                                          real64 const alpha,
                                          real64 const beta )
{

  GEOS_ASSERT_MSG( C.size( 0 ) == A.size( 1 ) &&
                       C.size( 1 ) == B.size( 1 ) &&
                       A.size( 0 ) == B.size( 0 ),
                   "Matrix dimensions not compatible for product" );

  int M = integer_conversion<int>( A.size( 1 ) );
  int N = integer_conversion<int>( B.size( 1 ) );
  int K = integer_conversion<int>( A.size( 0 ) );

//  cblas_dgemm( CblasRowMajor,
//               CblasTrans,
//               CblasNoTrans,
//               M,
//               N,
//               K,
//               alpha,
//               A.data(),
//               M,
//               B.data(),
//               N,
//               beta,
//               C.data(),
//               N );
  // CblasColMajor used for better performance. A^T*B = C is computed
  // as B^T * A = C^T
  cblas_dgemm( CblasColMajor,
               CblasNoTrans,
               CblasTrans,
               N,
               M,
               K,
               alpha,
               B.data(),
               N,
               A.data(),
               M,
               beta,
               C.data(),
               N );
  return;
}

void BlasLapackLA::matrixMatrixTMultiply( array2d<real64> const & A,
                                          array2d<real64> const & B,
                                          array2d<real64> & C,
                                          real64 const alpha,
                                          real64 const beta )
{

  GEOS_ASSERT_MSG( C.size( 0 ) == A.size( 0 ) &&
                       C.size( 1 ) == B.size( 0 ) &&
                       A.size( 1 ) == B.size( 1 ),
                   "Matrix dimensions not compatible for product" );

  int M = integer_conversion<int>( A.size( 0 ) );
  int N = integer_conversion<int>( B.size( 0 ) );
  int K = integer_conversion<int>( A.size( 1 ) );

//  cblas_dgemm( CblasRowMajor,
//               CblasNoTrans,
//               CblasTrans,
//               M,
//               N,
//               K,
//               alpha,
//               A.data(),
//               K,
//               B.data(),
//               K,
//               beta,
//               C.data(),
//               N );
  // CblasColMajor used for better performance. A*B^T = C is computed
  // as B * A^T = C^T
  cblas_dgemm( CblasColMajor,
               CblasTrans,
               CblasNoTrans,
               N,
               M,
               K,
               alpha,
               B.data(),
               K,
               A.data(),
               K,
               beta,
               C.data(),
               N );
  return;
}

void BlasLapackLA::matrixTMatrixTMultiply( array2d<real64> const & A,
                                           array2d<real64> const & B,
                                           array2d<real64> & C,
                                           real64 const alpha,
                                           real64 const beta )
{

  GEOS_ASSERT_MSG( C.size( 0 ) == A.size( 1 ) &&
                       C.size( 1 ) == B.size( 0 ) &&
                       A.size( 0 ) == B.size( 1 ),
                   "Matrix dimensions not compatible for product" );

  int M = integer_conversion<int>( A.size( 1 ) );
  int N = integer_conversion<int>( B.size( 0 ) );
  int K = integer_conversion<int>( A.size( 0 ) );

//  cblas_dgemm( CblasRowMajor,
//               CblasTrans,
//               CblasTrans,
//               M,
//               N,
//               K,
//               alpha,
//               A.data(),
//               M,
//               B.data(),
//               K,
//               beta,
//               C.data(),
//               N );
  // CblasColMajor used for better performance. A^T*B^T = C is computed
  // as B * A = C^T
  cblas_dgemm( CblasColMajor,
               CblasTrans,
               CblasTrans,
               N,
               M,
               K,
               alpha,
               B.data(),
               K,
               A.data(),
               M,
               beta,
               C.data(),
               N );
  return;
}

void BlasLapackLA::matrixInverse( array2d<real64> const & A,
                                  array2d<real64> & Ainv )
{
  real64 detA;
  matrixInverse( A,
                 Ainv,
                 detA );
}

void BlasLapackLA::matrixInverse( array2d<real64> const & A,
                                  array2d<real64> & Ainv,
                                  real64 & detA )
{
  // --- Check that source matrix is square
  int order = integer_conversion<lapack_int>(A.size( 0 ));
  GEOS_ASSERT_MSG( order > 0 &&
                   order == A.size( 1 ),
                   "Matrix must be square" );

  // --- Check that inverse matrix has appropriate dimension
  GEOS_ASSERT_MSG( Ainv.size( 0 ) == order &&
                   Ainv.size( 1 ) == order,
                   "Inverse matrix has wrong dimensions" );

  // --- Check if matrix is singular by computing the determinant
  //     note: if order greater than 3 we compute the determinant by
  //           first constructing the LU factors, later reused for calculating
  //           the inverse.
  lapack_int NN = order;
  array1d<lapack_int> IPIV;
  array1d<double> INV_WORK;

  if (order <= 3)
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
    // (avoid computing twice LUfactors, currently stored in Ainv, needed for
    // computing the inverse). Again we compute the LU factors for the
    // transpose matrix, i.e. choosing the LAPACK_COL_MAJOR ordering, to
    // avoid transposition/copy requires for LAPACK_ROW_MAJOR ordering.
    lapack_int INFO = LAPACKE_dgetrf( LAPACK_COL_MAJOR,
                                      NN,
                                      NN,
                                      Ainv.data(),
                                      NN,
                                      IPIV.data() );
    GEOS_ASSERT_MSG( INFO == 0, "LAPACKE_dgetrf error code: " << INFO );

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
  GEOS_ASSERT_MSG( std::fabs(detA) >
                   std::numeric_limits<real64>::epsilon() *
                   matrixNormFrobenius(A),
                   "Matrix is singular" );
  real64 oneOverDetA = 1. / detA;

  // --- Compute inverse
  switch( order )
    {
      case 1:
        Ainv( 0, 0 ) = oneOverDetA;
        return;

        // Case 2 to 4 copied from deal.ii full_matrix.templates.h (Maple generated)
      case 2:
        {
        Ainv( 0, 0 ) =  A( 1, 1 ) * oneOverDetA;
        Ainv( 0, 1 ) = -A( 0, 1 ) * oneOverDetA;
        Ainv( 1, 0 ) = -A( 1, 0 ) * oneOverDetA;
        Ainv( 1, 1 ) =  A( 0, 0 ) * oneOverDetA;
        return;
      }
        ;

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
        return;
      }
      default:
    {
    // Invert (LAPACK function DGETRI). The LU factors computed for the
    // transpose matrix stored in Ainv are used.
    lapack_int INFO = LAPACKE_dgetri( LAPACK_COL_MAJOR,
                           NN,
                           Ainv.data(),
                           NN,
                           IPIV.data() );
    GEOS_ASSERT_MSG( INFO == 0, "LAPACKE_dgetri error code: " << INFO );

    return;
    }
  }
}

void BlasLapackLA::vectorCopy( array1d<real64> const & X,
                               array1d<real64> & Y )
{
  GEOS_ASSERT_MSG( X.size() == Y.size(),
                   "Vector dimensions not compatible for copying" );

  // Call to BLAS using CBLAS interface
  cblas_dcopy( integer_conversion<int>( X.size() ),
               X.data(),
               1,
               Y.data(),
               1 );
  return;
}

void BlasLapackLA::matrixCopy( array2d<real64> const & A,
                               array2d<real64> & B )
{
  GEOS_ASSERT_MSG( A.size(0) == B.size(0) &&
                   A.size(1) == B.size(1),
                   "Matrix dimensions not compatible for copying" );

  // Call to BLAS using CBLAS interface
  cblas_dcopy( integer_conversion<int>( A.size(0)*A.size(1) ),
               A.data(),
               1,
               B.data(),
               1 );
  return;
}

//----------------------------------------------------------------I/O methods---
// vector nice output
void BlasLapackLA::printVector( array1d<real64> const & X )
{
  for( int i = 0 ; i < X.size() ; ++i )
  {
    printf( "%10.2e ", X[i] );
    printf( "\n" );
  }
}

// vector nice output
void BlasLapackLA::printMatrix( array2d<real64> const & X )
{
  for( int i = 0 ; i < X.size( 0 ) ; ++i )
  {
    for( int j = 0 ; j < X.size( 1 ) ; ++j )
    {
      printf( "%10.2e ", X( i, j ) );
    }
    printf( "\n" );
  }
}

} // end geosx namespace
