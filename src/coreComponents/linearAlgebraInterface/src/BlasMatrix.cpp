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

/**
 * @file LapackMatrix.cpp
 */

// Include the corresponding header file.
#include <BlasMatrix.hpp>

// Put everything under the geosx namespace.
namespace geosx
{

//----------------------------------------------Constructor/destructor methods---
BlasMatrix::BlasMatrix()
{

}

BlasMatrix::BlasMatrix( localIndex nRows, localIndex nCols )
{
  this->resize( nRows, nCols );
}

BlasMatrix::BlasMatrix( localIndex order )
{
  this->resize( order, order );
}

BlasMatrix::BlasMatrix( BlasMatrix const & src )
:
    m_numRows( src.m_numRows ),
    m_numCols( src.m_numCols ),
    m_values( src.m_values )
{

}

BlasMatrix::~BlasMatrix()
{

}

//-----------------------------------------------------Shaping/sizing methods---
void BlasMatrix::resize( localIndex nRows,
                         localIndex nCols )
{
  GEOS_ASSERT_MSG( nRows > 0, "Matrix number rows must be > 0" );
  GEOS_ASSERT_MSG( nCols > 0, "Matrix number of columns must be > 0" );
  m_numRows = nRows;
  m_numCols = nCols;
  m_values.resizeDefault( m_numRows * m_numCols );
  this->zero();
}

void BlasMatrix::resize( localIndex order )
{
  this->resize( order, order );
}

void BlasMatrix::zero()
{
  m_values = 0;
}

BlasMatrix &BlasMatrix::operator=(double value)
{
   m_values = value;
   return *this;
}

void BlasMatrix::permuteRows( array1d<int> permVector,
                              const bool forwardPermutation)
{
  // --- check that permutation vector has correct size
  GEOS_ASSERT_MSG( permVector.size() == m_numRows,
                   "Row permutation vector not consistent with matrix size" );
  LAPACKE_dlapmr( LAPACK_COL_MAJOR,
                  forwardPermutation,
                  integer_conversion<lapack_int>( m_numRows ),
                  integer_conversion<lapack_int >( m_numCols ),
                  m_values.data(),
                  integer_conversion<lapack_int>( m_numRows ),
                  permVector.data());

  return;
}

void BlasMatrix::permuteCols( array1d<int> permVector,
                              const bool forwardPermutation)
{
  // --- check that permutation vector has correct size
  GEOS_ASSERT_MSG( permVector.size() == m_numCols,
                   "Cols permutation vector not consistent with matrix size" );
  LAPACKE_dlapmt( LAPACK_COL_MAJOR,
                  forwardPermutation,
                  integer_conversion<lapack_int>( m_numRows ),
                  integer_conversion<lapack_int >( m_numCols ),
                  m_values.data(),
                  integer_conversion<lapack_int>( m_numRows ),
                  permVector.data());

  return;
}

//-------------------------------------------------------Mathematical methods---
// determinant calculation
real64 BlasMatrix::determinant() const
{
  // --- check that matrix is square
  GEOS_ASSERT_MSG( m_numRows == m_numCols && m_numRows > 0,
                   "Matrix must be square with order greater than zero" );

  switch( m_numRows )
  {
    case 1:
      return m_values[0];
    case 2:
      return m_values[0] * m_values[3] - m_values[1] * m_values[2];
    case 3:
      {
      real64 const * const v = m_values.data();
      return
          v[0] * ( v[4] * v[8] - v[5] * v[7] ) +
          v[3] * ( v[2] * v[7] - v[1] * v[8] ) +
          v[6] * ( v[1] * v[5] - v[2] * v[4] );
    }
    case 4:
      {
      real64 const * const v = m_values.data();
      return
      v[0] * ( v[5] * ( v[10] * v[15] - v[11] * v[14] ) -
               v[9] * ( v[6] * v[15] - v[7] * v[14] ) +
               v[13] * ( v[6] * v[11] - v[7] * v[10] )
             ) -
      v[4] * ( v[1] * ( v[10] * v[15] - v[11] * v[14] ) -
               v[9] * ( v[2] * v[15] - v[3] * v[14] ) +
               v[13] * ( v[2] * v[11] - v[3] * v[10] )
              ) +
      v[8] * ( v[1] * ( v[6] * v[15] - v[7] * v[14] ) -
               v[5] * ( v[2] * v[15] - v[3] * v[14] ) +
               v[13] * ( v[2] * v[7] - v[3] * v[6] )
          ) -
      v[12] * ( v[1] * ( v[6] * v[11] - v[7] * v[10] ) -
                v[5] * ( v[2] * v[11] - v[3] * v[10] ) +
                v[9] * ( v[2] * v[7] - v[3] * v[6] )
          );
    }
    default:

      // Compute the determinant via LU factorization
      lapack_int INFO;
      lapack_int NN = integer_conversion<lapack_int>( this->getNumRows() );
      array1d<lapack_int> IPIV( NN );
      array1d<double> LUfactor( m_values );

      INFO = LAPACKE_dgetrf( LAPACK_COL_MAJOR,
                             NN,
                             NN,
                             LUfactor.data(),
                             NN,
                             IPIV.data() );
      GEOS_ASSERT_MSG( INFO == 0,
                       "LAPACKE_dgetrf error code: " + std::to_string( INFO ) );

      real64 det = 1.0;
      for( int i = 0 ; i < NN ; ++i )
      {
        if( IPIV[i] != i + 1 ) //IPIV is based on Fortran convention (counting from 1)
        {
          det *= -LUfactor[NN * i + i];
        }
        else
        {
          det *= LUfactor[NN * i + i];
        }
      }
      return det;
  }
}

// computes inverse matrix
void BlasMatrix::computeInverse( BlasMatrix & dst )
{
  real64 det;
  this->computeInverse( dst, det );
}

// computes inverse matrix
void BlasMatrix::computeInverse( BlasMatrix & dst, real64& det )
{
  // --- Check that source matrix is square
  localIndex order = this->getNumRows();
  GEOS_ASSERT_MSG( order > 0 && order == this->getNumCols(),
                   "Matrix must be square" );

  // --- Initialize the inverse matrix to the appropriate dimension
  dst.resize( order,
              order );

  // --- Check if matrix is singular by computing the determinant
  //     note: if order greater than 3 we compute the determinant by
  //           first constructing the LU factors, later reused for calculating
  //           the inverse.
  lapack_int NN;
  array1d<lapack_int> IPIV;
  lapack_int INFO;
  array1d<double> INV_WORK;

  if (order <= 3)
  {
    det = this->determinant();
    real64 oneOverDet = 1. / this->determinant();
  }
  else
  {
    // Copy this in dst
    dst.m_values = this->m_values;

    // Declare workspace for permutations and scratch array
    NN = integer_conversion<lapack_int>( this->getNumCols() );
    IPIV.resize(NN);
    INV_WORK.resize(NN);

    // Call to LAPACK using LAPACKE
    // --- Compute LU factorization (LAPACK function DGETRF)
    INFO = LAPACKE_dgetrf( LAPACK_COL_MAJOR,
                           NN,
                           NN,
                           dst.m_values.data(),
                           NN,
                           IPIV.data() );
    GEOS_ASSERT_MSG( INFO == 0,
                     "LAPACKE_dgetrf error code: " + std::to_string( INFO ) );

    // --- Compute determinant (not done calling directly the function determinant
    det = 1.0;
    for( int i = 0 ; i < NN ; ++i )
    {
      if( IPIV[i] != i + 1 ) //IPIV is based on Fortran convention (counting from 1)
      {
        det *= -dst.m_values[NN * i + i];
      }
      else
      {
        det *= dst.m_values[NN * i + i];
      }
    }
  }

  real64 oneOverDet = 1. / this->determinant();
  GEOS_ASSERT_MSG( !( std::isinf( oneOverDet ) ), "Matrix is singular" );

  // --- Compute inverse
  switch( order )
    {
      case 1:
        dst( 0, 0 ) = oneOverDet;
        return;

        // Case 2 to 4 copied from deal.ii full_matrix.templates.h (Maple generated)
      case 2:
        {
        dst( 0, 0 ) = ( *this )( 1, 1 ) * oneOverDet;
        dst( 0, 1 ) = -( *this )( 0, 1 ) * oneOverDet;
        dst( 1, 0 ) = -( *this )( 1, 0 ) * oneOverDet;
        dst( 1, 1 ) = ( *this )( 0, 0 ) * oneOverDet;
        return;
      }
        ;

      case 3:
        {
        dst( 0, 0 ) = ( ( *this )( 1, 1 ) * ( *this )( 2, 2 ) -
                        ( *this )( 1, 2 ) * ( *this )( 2, 1 ) ) * oneOverDet;
        dst( 0, 1 ) = ( ( *this )( 0, 2 ) * ( *this )( 2, 1 ) -
                        ( *this )( 0, 1 ) * ( *this )( 2, 2 ) ) * oneOverDet;
        dst( 0, 2 ) = ( ( *this )( 0, 1 ) * ( *this )( 1, 2 ) -
                        ( *this )( 0, 2 ) * ( *this )( 1, 1 ) ) * oneOverDet;
        dst( 1, 0 ) = ( ( *this )( 1, 2 ) * ( *this )( 2, 0 ) -
                        ( *this )( 1, 0 ) * ( *this )( 2, 2 ) ) * oneOverDet;
        dst( 1, 1 ) = ( ( *this )( 0, 0 ) * ( *this )( 2, 2 ) -
                        ( *this )( 0, 2 ) * ( *this )( 2, 0 ) ) * oneOverDet;
        dst( 1, 2 ) = ( ( *this )( 0, 2 ) * ( *this )( 1, 0 ) -
                        ( *this )( 0, 0 ) * ( *this )( 1, 2 ) ) * oneOverDet;
        dst( 2, 0 ) = ( ( *this )( 1, 0 ) * ( *this )( 2, 1 ) -
                        ( *this )( 1, 1 ) * ( *this )( 2, 0 ) ) * oneOverDet;
        dst( 2, 1 ) = ( ( *this )( 0, 1 ) * ( *this )( 2, 0 ) -
                        ( *this )( 0, 0 ) * ( *this )( 2, 1 ) ) * oneOverDet;
        dst( 2, 2 ) = ( ( *this )( 0, 0 ) * ( *this )( 1, 1 ) -
                        ( *this )( 0, 1 ) * ( *this )( 1, 0 ) ) * oneOverDet;
        return;
      }
      default:
    {
    // --- Invert (LAPACK function DGETRI)
    INFO = LAPACKE_dgetri( LAPACK_COL_MAJOR,
                           NN,
                           dst.m_values.data(),
                           NN,
                           IPIV.data() );
    GEOS_ASSERT_MSG( INFO == 0,
                     "LAPACKE_dgetri error code: " + std::to_string( INFO ) );

    return;
    }
  }
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Inf-norm.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Returns the infinity norm of the matrix.
real64 BlasMatrix::normInf() const
{
  return LAPACKE_dlange( LAPACK_COL_MAJOR,
                         'I',
                         integer_conversion<lapack_int>(m_numRows),
                         integer_conversion<lapack_int>(m_numCols),
                         m_values.data(),
                         integer_conversion<lapack_int>(m_numRows));

}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// 1-norm.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Returns the one norm of the matrix.
real64 BlasMatrix::norm1() const
{
  return LAPACKE_dlange( LAPACK_COL_MAJOR,
                         '1',
                         integer_conversion<lapack_int>(m_numRows),
                         integer_conversion<lapack_int>(m_numCols),
                         m_values.data(),
                         integer_conversion<lapack_int>(m_numRows));

}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Frobenius-norm.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Returns the Frobenius norm of the matrix.
real64 BlasMatrix::normFrobenius() const
{
  return LAPACKE_dlange( LAPACK_COL_MAJOR,
                         'F',
                         integer_conversion<lapack_int>(m_numRows),
                         integer_conversion<lapack_int>(m_numCols),
                         m_values.data(),
                         integer_conversion<lapack_int>(m_numRows));

}

// matrix-matrix sum (optional scaling)
void BlasMatrix::matrixAdd( BlasMatrix const & A,
                            real64 const scalarA )
{

  GEOS_ASSERT_MSG( A.getNumRows() == m_numRows && A.getNumCols() == m_numCols,
                   "Matrix dimensions not compatible for sum" );

  cblas_daxpy( static_cast<int>( m_values.size() ),
               scalarA,
               A.m_values.data(),
               1,
               m_values.data(),
               1 );
  return;
}

// in-place scalar-matrix product
void BlasMatrix::scale(real64 scalarThis)
{
  // Call to BLAS using CBLAS interface
  cblas_dscal(integer_conversion<int>(m_values.size()),
              scalarThis,
              m_values.data(),
              1);
  return;
}

// matrix-matrix multiplication (optional scaling/accumulation)
void BlasMatrix::matrixMultiply( BlasMatrix const &src,
                                 BlasMatrix &dst,
                                 real64 const scalarThisSrc,
                                 real64 const scalarDst )
{

  GEOS_ASSERT_MSG( dst.getNumRows() == m_numRows &&
                       dst.getNumCols() == src.getNumCols() &&
                       m_numCols == src.getNumRows(),
                   "Matrix dimensions not compatible for product" );

  int M = integer_conversion<int>( m_numRows );
  int N = integer_conversion<int>( src.getNumCols() );
  int K = integer_conversion<int>( m_numCols );

  cblas_dgemm( CblasColMajor,
               CblasNoTrans,
               CblasNoTrans,
               M,
               N,
               K,
               scalarThisSrc,
               m_values.data(),
               M,
               src.m_values.data(),
               K,
               scalarDst,
               dst.m_values.data(),
               M );
}

// transpose(matrix)-matrix multiplication (optional scaling/accumulation)
void BlasMatrix::TmatrixMultiply( BlasMatrix const &src,
                                  BlasMatrix &dst,
                                  real64 const scalarThisSrc,
                                  real64 const scalarDst )
{

  GEOS_ASSERT_MSG( dst.getNumRows() == m_numCols &&
                       dst.getNumCols() == src.getNumCols() &&
                       m_numRows == src.getNumRows(),
                   "Matrix dimensions not compatible for product" );

  int M = integer_conversion<int>( m_numCols );
  int N = integer_conversion<int>( src.getNumCols() );
  int K = integer_conversion<int>( m_numRows );

  cblas_dgemm( CblasColMajor,
               CblasTrans,
               CblasNoTrans,
               M,
               N,
               K,
               scalarThisSrc,
               m_values.data(),
               K,
               src.m_values.data(),
               K,
               scalarDst,
               dst.m_values.data(),
               M );
}

// transpose(matrix)-matrix multiplication (optional scaling/accumulation)
void BlasMatrix::matrixTMultiply( BlasMatrix const &src,
                                  BlasMatrix &dst,
                                  real64 const scalarThisSrc,
                                  real64 const scalarDst )
{

  GEOS_ASSERT_MSG( dst.getNumRows() == m_numRows &&
                       dst.getNumCols() == src.getNumRows() &&
                       m_numCols == src.getNumCols(),
                   "Matrix dimensions not compatible for product" );

  int M = integer_conversion<int>( m_numRows );
  int N = integer_conversion<int>( src.getNumRows() );
  int K = integer_conversion<int>( m_numCols );

  cblas_dgemm( CblasColMajor,
               CblasNoTrans,
               CblasTrans,
               M,
               N,
               K,
               scalarThisSrc,
               m_values.data(),
               M,
               src.m_values.data(),
               N,
               scalarDst,
               dst.m_values.data(),
               M );
}

// transpose(matrix)-transpose(matrix) multiplication
// (optional scaling/accumulation)
void BlasMatrix::TmatrixTMultiply( BlasMatrix const &src,
                                   BlasMatrix &dst,
                                   real64 const scalarThisSrc,
                                   real64 const scalarDst )
{

  GEOS_ASSERT_MSG( dst.getNumRows() == m_numCols &&
                       dst.getNumCols() == src.getNumRows() &&
                       m_numRows == src.getNumCols(),
                   "Matrix dimensions not compatible for product" );

  int M = integer_conversion<int>( m_numCols );
  int N = integer_conversion<int>( src.getNumRows() );
  int K = integer_conversion<int>( m_numRows );

  cblas_dgemm( CblasColMajor,
               CblasTrans,
               CblasTrans,
               M,
               N,
               K,
               scalarThisSrc,
               m_values.data(),
               K,
               src.m_values.data(),
               N,
               scalarDst,
               dst.m_values.data(),
               M );
}

// matrix-vector multiplication
void BlasMatrix::vectorMultiply(BlasMatrix const &src,
                                BlasMatrix &dst,
                                real64 const scalarThisSrc,
                                real64 const scalarDst)
{
  GEOS_ASSERT_MSG(m_numCols == src.getNumRows() && m_numRows == dst.getNumRows(),
                  "Matrix, source vector and destination vector not compatible");

  this->matrixMultiply(src,
                       dst,
                       scalarThisSrc,
                       scalarDst);

  return;
}

// transpose(matrix)-vector multiplication
void BlasMatrix::TvectorMultiply(BlasMatrix const &src,
                                 BlasMatrix &dst,
                                 real64 const scalarThisSrc,
                                 real64 const scalarDst)
{
  GEOS_ASSERT_MSG(m_numRows == src.getNumRows() && m_numCols == dst.getNumRows(),
                  "Matrix, source vector and destination vector not compatible");

  this->TmatrixMultiply(src,
                        dst,
                        scalarThisSrc,
                        scalarDst);

  return;
}


//
//------------------------------------------------------Data Accessor methods---

// Returns number of matrix rows.
localIndex BlasMatrix::getNumRows() const
{
  return m_numRows;
}

// Returns number of matrix columns.
localIndex BlasMatrix::getNumCols() const
{
  return m_numCols;
}

//----------------------------------------------------------------I/O methods---
// matrix nice output
void BlasMatrix::print()
{
  for( localIndex i = 0 ; i < m_numRows ; ++i )
  {
    for( localIndex j = 0 ; j < m_numCols ; ++j )
      printf( "%10.2e ", m_values[j * m_numRows + i] );
    printf( "\n" );
  }
}

} // end geosx namespace
