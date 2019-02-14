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
  m_values.resizeDefault( m_numRows * m_numCols);
  m_values = 0.0;
}

void BlasMatrix::resize( localIndex order )
{
  this->resize( order, order );
}

//-------------------------------------------------------Mathematical methods---
// determinant calculation
real64 BlasMatrix::determinant()
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
      v[ 0] * ( v[ 5] * ( v[10] * v[15] - v[11] * v[14] ) -
                v[ 9] * ( v[ 6] * v[15] - v[ 7] * v[14] ) +
                v[13] * ( v[ 6] * v[11] - v[ 7] * v[10] )
              ) -
      v[ 4] * ( v[ 1] * ( v[10] * v[15] - v[11] * v[14] ) -
                v[ 9] * ( v[ 2] * v[15] - v[ 3] * v[14] ) +
                v[13] * ( v[ 2] * v[11] - v[ 3] * v[10] )
              ) +
      v[ 8] * ( v[ 1] * ( v[ 6] * v[15] - v [7] * v[14] ) -
                v[ 5] * ( v[ 2] * v[15] - v[ 3] * v[14] ) +
                v[13] * ( v[ 2] * v[ 7] - v[ 3] * v[ 6] )
              ) -
      v[12] * ( v[ 1] * ( v[ 6] * v[11] - v[ 7] * v[10] ) -
                v[ 5] * ( v[ 2] * v[11] - v[ 3] * v[10] ) +
                v[ 9] * ( v[ 2] * v[ 7] - v[ 3] * v[ 6] )
              );
     }
    default:

      // Compute the determinant via LU factorization
      lapack_int INFO;
      lapack_int NN = integer_conversion<lapack_int>( this->getNumRows() );
      array1d<lapack_int> IPIV( NN );
      array1d<double> LUfactor(m_values);

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
        if( IPIV[i] != i+1 ) //IPIV is based on Fortran convention (counting from 1)
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
void BlasMatrix::computeInverse( BlasMatrix & dst)
{
  real64 det;
  this->computeInverse( dst, det);
}

// computes inverse matrix
void BlasMatrix::computeInverse( BlasMatrix & dst, real64& det )
{
  // --- Check that source matrix is square
  localIndex order = this->getNumRows();
  GEOS_ASSERT_MSG( order == this->getNumCols(), "Matrix must be square" );

  // --- Initialize the inverse matrix to the appropriate dimension
  dst.resize( order,
              order );

  det = this->determinant();
  real64 oneOverDet = 1. / this->determinant();
  GEOS_ASSERT_MSG( !(std::isinf( oneOverDet )) , "Matrix is singular" );

  switch( m_numCols )
  {
    case 1:
      dst( 0, 0 ) = oneOverDet;
      return;

      // Case 2 to 4 copied from deal.ii full_matrix.templates.h (Maple generated)
    case 2:
      {
      const double t4 = 1. / ( ( *this )( 0, 0 ) * ( *this )( 1, 1 ) - ( *this )( 0, 1 ) * ( *this )( 1, 0 ) );
      dst( 0, 0 ) = ( *this )( 1, 1 ) * t4;
      dst( 0, 1 ) = -( *this )( 0, 1 ) * t4;
      dst( 1, 0 ) = -( *this )( 1, 0 ) * t4;
      dst( 1, 1 ) = ( *this )( 0, 0 ) * t4;
      return;
    }
      ;

    case 3:
      {
      const double t4 = ( *this )( 0, 0 ) * ( *this )( 1, 1 ),
          t6 = ( *this )( 0, 0 ) * ( *this )( 1, 2 ),
          t8 = ( *this )( 0, 1 ) * ( *this )( 1, 0 ),
          t00 = ( *this )( 0, 2 ) * ( *this )( 1, 0 ),
          t01 = ( *this )( 0, 1 ) * ( *this )( 2, 0 ),
          t04 = ( *this )( 0, 2 ) * ( *this )( 2, 0 ),
          t07 = 1. / ( t4 * ( *this )( 2, 2 ) - t6 * ( *this )( 2, 1 ) - t8 * ( *this )( 2, 2 ) +
              t00 * ( *this )( 2, 1 ) + t01 * ( *this )( 1, 2 ) - t04 * ( *this )( 1, 1 ) );
      dst( 0, 0 ) = ( ( *this )( 1, 1 ) * ( *this )( 2, 2 ) - ( *this )( 1, 2 ) * ( *this )( 2, 1 ) ) * t07;
      dst( 0, 1 ) = -( ( *this )( 0, 1 ) * ( *this )( 2, 2 ) - ( *this )( 0, 2 ) * ( *this )( 2, 1 ) ) * t07;
      dst( 0, 2 ) = -( -( *this )( 0, 1 ) * ( *this )( 1, 2 ) + ( *this )( 0, 2 ) * ( *this )( 1, 1 ) ) * t07;
      dst( 1, 0 ) = -( ( *this )( 1, 0 ) * ( *this )( 2, 2 ) - ( *this )( 1, 2 ) * ( *this )( 2, 0 ) ) * t07;
      dst( 1, 1 ) = ( ( *this )( 0, 0 ) * ( *this )( 2, 2 ) - t04 ) * t07;
      dst( 1, 2 ) = -( t6 - t00 ) * t07;
      dst( 2, 0 ) = -( -( *this )( 1, 0 ) * ( *this )( 2, 1 ) + ( *this )( 1, 1 ) * ( *this )( 2, 0 ) ) * t07;
      dst( 2, 1 ) = -( ( *this )( 0, 0 ) * ( *this )( 2, 1 ) - t01 ) * t07;
      dst( 2, 2 ) = ( t4 - t8 ) * t07;
      return;
    }
      ;

    case 4:
      {
      const double t14 = ( *this )( 0, 0 ) * ( *this )( 1, 1 );
      const double t15 = ( *this )( 2, 2 ) * ( *this )( 3, 3 );
      const double t17 = ( *this )( 2, 3 ) * ( *this )( 3, 2 );
      const double t19 = ( *this )( 0, 0 ) * ( *this )( 2, 1 );
      const double t20 = ( *this )( 1, 2 ) * ( *this )( 3, 3 );
      const double t22 = ( *this )( 1, 3 ) * ( *this )( 3, 2 );
      const double t24 = ( *this )( 0, 0 ) * ( *this )( 3, 1 );
      const double t25 = ( *this )( 1, 2 ) * ( *this )( 2, 3 );
      const double t27 = ( *this )( 1, 3 ) * ( *this )( 2, 2 );
      const double t29 = ( *this )( 1, 0 ) * ( *this )( 0, 1 );
      const double t32 = ( *this )( 1, 0 ) * ( *this )( 2, 1 );
      const double t33 = ( *this )( 0, 2 ) * ( *this )( 3, 3 );
      const double t35 = ( *this )( 0, 3 ) * ( *this )( 3, 2 );
      const double t37 = ( *this )( 1, 0 ) * ( *this )( 3, 1 );
      const double t38 = ( *this )( 0, 2 ) * ( *this )( 2, 3 );
      const double t40 = ( *this )( 0, 3 ) * ( *this )( 2, 2 );
      const double t42 = t14 * t15 - t14 * t17 - t19 * t20 + t19 * t22 +
          t24 * t25 - t24 * t27 - t29 * t15 + t29 * t17 +
          t32 * t33 - t32 * t35 - t37 * t38 + t37 * t40;
      const double t43 = ( *this )( 2, 0 ) * ( *this )( 0, 1 );
      const double t46 = ( *this )( 2, 0 ) * ( *this )( 1, 1 );
      const double t49 = ( *this )( 2, 0 ) * ( *this )( 3, 1 );
      const double t50 = ( *this )( 0, 2 ) * ( *this )( 1, 3 );
      const double t52 = ( *this )( 0, 3 ) * ( *this )( 1, 2 );
      const double t54 = ( *this )( 3, 0 ) * ( *this )( 0, 1 );
      const double t57 = ( *this )( 3, 0 ) * ( *this )( 1, 1 );
      const double t60 = ( *this )( 3, 0 ) * ( *this )( 2, 1 );
      const double t63 = t43 * t20 - t43 * t22 - t46 * t33 + t46 * t35 +
          t49 * t50 - t49 * t52 - t54 * t25 + t54 * t27 +
          t57 * t38 - t57 * t40 - t60 * t50 + t60 * t52;
      const double t65 = 1. / ( t42 + t63 );
      const double t71 = ( *this )( 0, 2 ) * ( *this )( 2, 1 );
      const double t73 = ( *this )( 0, 3 ) * ( *this )( 2, 1 );
      const double t75 = ( *this )( 0, 2 ) * ( *this )( 3, 1 );
      const double t77 = ( *this )( 0, 3 ) * ( *this )( 3, 1 );
      const double t81 = ( *this )( 0, 1 ) * ( *this )( 1, 2 );
      const double t83 = ( *this )( 0, 1 ) * ( *this )( 1, 3 );
      const double t85 = ( *this )( 0, 2 ) * ( *this )( 1, 1 );
      const double t87 = ( *this )( 0, 3 ) * ( *this )( 1, 1 );
      const double t101 = ( *this )( 1, 0 ) * ( *this )( 2, 2 );
      const double t103 = ( *this )( 1, 0 ) * ( *this )( 2, 3 );
      const double t105 = ( *this )( 2, 0 ) * ( *this )( 1, 2 );
      const double t107 = ( *this )( 2, 0 ) * ( *this )( 1, 3 );
      const double t109 = ( *this )( 3, 0 ) * ( *this )( 1, 2 );
      const double t111 = ( *this )( 3, 0 ) * ( *this )( 1, 3 );
      const double t115 = ( *this )( 0, 0 ) * ( *this )( 2, 2 );
      const double t117 = ( *this )( 0, 0 ) * ( *this )( 2, 3 );
      const double t119 = ( *this )( 2, 0 ) * ( *this )( 0, 2 );
      const double t121 = ( *this )( 2, 0 ) * ( *this )( 0, 3 );
      const double t123 = ( *this )( 3, 0 ) * ( *this )( 0, 2 );
      const double t125 = ( *this )( 3, 0 ) * ( *this )( 0, 3 );
      const double t129 = ( *this )( 0, 0 ) * ( *this )( 1, 2 );
      const double t131 = ( *this )( 0, 0 ) * ( *this )( 1, 3 );
      const double t133 = ( *this )( 1, 0 ) * ( *this )( 0, 2 );
      const double t135 = ( *this )( 1, 0 ) * ( *this )( 0, 3 );
      dst( 0, 0 ) = ( ( *this )( 1, 1 ) * ( *this )( 2, 2 ) * ( *this )( 3, 3 ) - ( *this )( 1, 1 ) * ( *this )( 2, 3 ) * ( *this )(
          3, 2 ) -
          ( *this )( 2, 1 ) * ( *this )( 1, 2 ) * ( *this )( 3, 3 ) + ( *this )( 2, 1 ) * ( *this )( 1, 3 ) * ( *this )(
          3, 2 ) +
          ( *this )( 3, 1 ) * ( *this )( 1, 2 ) * ( *this )( 2, 3 ) - ( *this )( 3, 1 ) * ( *this )( 1, 3 ) * ( *this )(
          2, 2 ) ) * t65;
      dst( 0, 1 ) = -( ( *this )( 0, 1 ) * ( *this )( 2, 2 ) * ( *this )( 3, 3 ) - ( *this )( 0, 1 ) * ( *this )( 2, 3 ) * ( *this )(
          3, 2 ) -
          t71 * ( *this )( 3, 3 ) + t73 * ( *this )( 3, 2 ) + t75 * ( *this )( 2, 3 ) - t77 * ( *this )( 2, 2 ) ) * t65;
      dst( 0, 2 ) = ( t81 * ( *this )( 3, 3 ) - t83 * ( *this )( 3, 2 ) - t85 * ( *this )( 3, 3 ) + t87 * ( *this )( 3,
                                                                                                                     2 ) +
          t75 * ( *this )( 1, 3 ) - t77 * ( *this )( 1, 2 ) ) * t65;
      dst( 0, 3 ) = -( t81 * ( *this )( 2, 3 ) - t83 * ( *this )( 2, 2 ) - t85 * ( *this )( 2, 3 ) + t87 * ( *this )(
          2, 2 ) +
          t71 * ( *this )( 1, 3 ) - t73 * ( *this )( 1, 2 ) ) * t65;
      dst( 1, 0 ) = -( t101 * ( *this )( 3, 3 ) - t103 * ( *this )( 3, 2 ) - t105 * ( *this )( 3, 3 ) + t107 * ( *this )(
          3, 2 ) +
          t109 * ( *this )( 2, 3 ) - t111 * ( *this )( 2, 2 ) ) * t65;
      dst( 1, 1 ) = ( t115 * ( *this )( 3, 3 ) - t117 * ( *this )( 3, 2 ) - t119 * ( *this )( 3, 3 ) + t121 * ( *this )(
          3, 2 ) +
          t123 * ( *this )( 2, 3 ) - t125 * ( *this )( 2, 2 ) ) * t65;
      dst( 1, 2 ) = -( t129 * ( *this )( 3, 3 ) - t131 * ( *this )( 3, 2 ) - t133 * ( *this )( 3, 3 ) + t135 * ( *this )(
          3, 2 ) +
          t123 * ( *this )( 1, 3 ) - t125 * ( *this )( 1, 2 ) ) * t65;
      dst( 1, 3 ) = ( t129 * ( *this )( 2, 3 ) - t131 * ( *this )( 2, 2 ) - t133 * ( *this )( 2, 3 ) + t135 * ( *this )(
          2, 2 ) +
          t119 * ( *this )( 1, 3 ) - t121 * ( *this )( 1, 2 ) ) * t65;
      dst( 2, 0 ) = ( t32 * ( *this )( 3, 3 ) - t103 * ( *this )( 3, 1 ) - t46 * ( *this )( 3, 3 ) + t107 * ( *this )(
          3, 1 ) +
          t57 * ( *this )( 2, 3 ) - t111 * ( *this )( 2, 1 ) ) * t65;
      dst( 2, 1 ) = -( t19 * ( *this )( 3, 3 ) - t117 * ( *this )( 3, 1 ) - t43 * ( *this )( 3, 3 ) + t121 * ( *this )(
          3, 1 ) +
          t54 * ( *this )( 2, 3 ) - t125 * ( *this )( 2, 1 ) ) * t65;
      dst( 2, 2 ) = ( t14 * ( *this )( 3, 3 ) - t131 * ( *this )( 3, 1 ) - t29 * ( *this )( 3, 3 ) + t135 * ( *this )(
          3, 1 ) +
          t54 * ( *this )( 1, 3 ) - t125 * ( *this )( 1, 1 ) ) * t65;
      dst( 2, 3 ) = -( t14 * ( *this )( 2, 3 ) - t131 * ( *this )( 2, 1 ) - t29 * ( *this )( 2, 3 ) + t135 * ( *this )(
          2, 1 ) +
          t43 * ( *this )( 1, 3 ) - t121 * ( *this )( 1, 1 ) ) * t65;
      dst( 3, 0 ) = -( t32 * ( *this )( 3, 2 ) - t101 * ( *this )( 3, 1 ) - t46 * ( *this )( 3, 2 ) + t105 * ( *this )(
          3, 1 ) +
          t57 * ( *this )( 2, 2 ) - t109 * ( *this )( 2, 1 ) ) * t65;
      dst( 3, 1 ) = ( t19 * ( *this )( 3, 2 ) - t115 * ( *this )( 3, 1 ) - t43 * ( *this )( 3, 2 ) + t119 * ( *this )(
          3, 1 ) +
          t54 * ( *this )( 2, 2 ) - t123 * ( *this )( 2, 1 ) ) * t65;
      dst( 3, 2 ) = -( t14 * ( *this )( 3, 2 ) - t129 * ( *this )( 3, 1 ) - t29 * ( *this )( 3, 2 ) + t133 * ( *this )(
          3, 1 ) +
          t54 * ( *this )( 1, 2 ) - t123 * ( *this )( 1, 1 ) ) * t65;
      dst( 3, 3 ) = ( t14 * ( *this )( 2, 2 ) - t129 * ( *this )( 2, 1 ) - t29 * ( *this )( 2, 2 ) + t133 * ( *this )(
          2, 1 ) +
          t43 * ( *this )( 1, 2 ) - t119 * ( *this )( 1, 1 ) ) * t65;

      return;
    }

    default:
      {
      // Copy M in this
      dst.m_values = this->m_values;

      // Declare workspace for permutations and scratch array
      lapack_int NN = integer_conversion<lapack_int>( this->getNumCols() );
      array1d<lapack_int> IPIV( NN );
      lapack_int INFO;
      array1d<double> INV_WORK( NN );

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

// matrix-matrix sum (optional scaling)
void BlasMatrix::MatAdd( BlasMatrix const & A,
                         real64 const scalarA )
{
  unsigned int nRowA = A.getNumRows();
  unsigned int nColA = A.getNumCols();

  GEOS_ASSERT_MSG( A.getNumRows() == m_numRows,
                   "Matrix dimensions not compatible for sum" );
  GEOS_ASSERT_MSG( A.getNumCols() == m_numCols,
                   "Matrix dimensions not compatible for sum" );

  cblas_daxpy( static_cast<int>( m_values.size() ),
               scalarA,
               A.m_values.data(),
               1,
               m_values.data(),
               1 );
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

//// matrix-vector multiplication
//void SerialDenseMatrix::MatVecMult(SerialDenseMatrix& x,
//                                   SerialDenseMatrix& y)
//{
//  ASSERT(this->m_nRows == y.get_nRows(), "Matrix and destination vector "
//                                         "dimensions not compatible for "
//                                         "product");
//  ASSERT(this->m_nCols == x.get_nRows(), "Matrix and destination vector "
//                                         "dimensions not compatible for "
//                                         "product");
//
//  y.MatMatMult(*this, x);
//
//  return;
//}
//
//// matrix-vector multiplication
//void SerialDenseMatrix::MatTVecMult(SerialDenseMatrix& x,
//                                    SerialDenseMatrix& y)
//{
//  ASSERT(this->m_nCols == y.get_nRows(), "Matrix and destination vector "
//                                         "dimensions not compatible for "
//                                         "product");
//  ASSERT(this->m_nRows == x.get_nRows(), "Matrix and destination vector "
//                                         "dimensions not compatible for "
//                                         "product");
//
//  y.MatTMatMult(*this, x);
//
//  return;
//}
//
//// in-place scalar-matrix product
//void SerialDenseMatrix::Scale(double scalarThis)
//{
//
//#ifdef WITH_TRILINOS
//  // Call to BLAS using Epetra BLAS Wrapper Class (Epetra_BLAS)
//  Epetra_BLAS Blas;
//  Blas.SCAL(static_cast<int>(m_data.size()), scalarThis, &(*this)(0,0), 1);
//  return;
//#else
//  // Data accessing by column
//  for (unsigned int j = 0; j < m_nCols; ++j)
//    for (unsigned int i = 0; i  <m_nRows; ++i)
//       (*this)(i,j) = scalarThis*(*this)(i,j);
//#endif
//
//}
//
//------------------------------------------------------Data Accessor methods---

//// matrix nice output
//void SerialDenseMatrix::assign_value(double value)
//{
//  std::fill(m_data.begin(), m_data.end(), value);
//  return;
//};

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
