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
    m_height( src.m_height ),
    m_width( src.m_width ),
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
  m_height = nRows;
  m_width = nCols;
  m_values.resizeDefault( m_height * m_width, 0.0 );
}

void BlasMatrix::resize( localIndex order )
{
  this->resize( order, order );
}

//-------------------------------------------------------Mathematical methods---
// determinant calculation
double BlasMatrix::determinant()
{
  // --- check that matrix is square
  GEOS_ASSERT_MSG( m_height == m_width, "Matrix must be square" );

  switch( m_width )
  {
    case 1:
      return m_values[0];
    case 2:
      return m_values[0] * m_values[3] - m_values[2] * m_values[1];
    case 3:
      return m_values[0] *
          ( m_values[4] * m_values[8] - m_values[7] * m_values[5] )
          - m_values[1] *
              ( m_values[3] * m_values[8] - m_values[6] * m_values[5] )
          + m_values[2] *
              ( m_values[3] * m_values[7] - m_values[6] * m_values[4] );
    default:
      GEOS_ERROR_IF( true, "Determinant computation implemented up to matrix 3x3" );
      return 1;
  }
}

// computes inverse matrix
void BlasMatrix::invert( BlasMatrix & src )
{
  // --- Check that source matrix is square
  localIndex order = src.get_nRows();
  GEOS_ASSERT_MSG( order == src.get_nCols(), "Source matrix must be square" );

  // --- Initialize the inverse matrix to the appropriate dimension
  this->resize( order,
                order );

  switch( m_width )
  {
    case 1:
      ( *this )( 0, 0 ) = 1. / src( 0, 0 );
      return;

      // Case 2 to 4 copied from deal.ii full_matrix.templates.h (Maple generated)
    case 2:
      {
      const double t4 = 1. / ( src( 0, 0 ) * src( 1, 1 ) - src( 0, 1 ) * src( 1, 0 ) );
      ( *this )( 0, 0 ) = src( 1, 1 ) * t4;
      ( *this )( 0, 1 ) = -src( 0, 1 ) * t4;
      ( *this )( 1, 0 ) = -src( 1, 0 ) * t4;
      ( *this )( 1, 1 ) = src( 0, 0 ) * t4;
      return;
    }
      ;

    case 3:
      {
      const double t4 = src( 0, 0 ) * src( 1, 1 ),
          t6 = src( 0, 0 ) * src( 1, 2 ),
          t8 = src( 0, 1 ) * src( 1, 0 ),
          t00 = src( 0, 2 ) * src( 1, 0 ),
          t01 = src( 0, 1 ) * src( 2, 0 ),
          t04 = src( 0, 2 ) * src( 2, 0 ),
          t07 = 1. / ( t4 * src( 2, 2 ) - t6 * src( 2, 1 ) - t8 * src( 2, 2 ) +
              t00 * src( 2, 1 ) + t01 * src( 1, 2 ) - t04 * src( 1, 1 ) );
      ( *this )( 0, 0 ) = ( src( 1, 1 ) * src( 2, 2 ) - src( 1, 2 ) * src( 2, 1 ) ) * t07;
      ( *this )( 0, 1 ) = -( src( 0, 1 ) * src( 2, 2 ) - src( 0, 2 ) * src( 2, 1 ) ) * t07;
      ( *this )( 0, 2 ) = -( -src( 0, 1 ) * src( 1, 2 ) + src( 0, 2 ) * src( 1, 1 ) ) * t07;
      ( *this )( 1, 0 ) = -( src( 1, 0 ) * src( 2, 2 ) - src( 1, 2 ) * src( 2, 0 ) ) * t07;
      ( *this )( 1, 1 ) = ( src( 0, 0 ) * src( 2, 2 ) - t04 ) * t07;
      ( *this )( 1, 2 ) = -( t6 - t00 ) * t07;
      ( *this )( 2, 0 ) = -( -src( 1, 0 ) * src( 2, 1 ) + src( 1, 1 ) * src( 2, 0 ) ) * t07;
      ( *this )( 2, 1 ) = -( src( 0, 0 ) * src( 2, 1 ) - t01 ) * t07;
      ( *this )( 2, 2 ) = ( t4 - t8 ) * t07;
      return;
    }
      ;

    case 4:
      {
      const double t14 = src( 0, 0 ) * src( 1, 1 );
      const double t15 = src( 2, 2 ) * src( 3, 3 );
      const double t17 = src( 2, 3 ) * src( 3, 2 );
      const double t19 = src( 0, 0 ) * src( 2, 1 );
      const double t20 = src( 1, 2 ) * src( 3, 3 );
      const double t22 = src( 1, 3 ) * src( 3, 2 );
      const double t24 = src( 0, 0 ) * src( 3, 1 );
      const double t25 = src( 1, 2 ) * src( 2, 3 );
      const double t27 = src( 1, 3 ) * src( 2, 2 );
      const double t29 = src( 1, 0 ) * src( 0, 1 );
      const double t32 = src( 1, 0 ) * src( 2, 1 );
      const double t33 = src( 0, 2 ) * src( 3, 3 );
      const double t35 = src( 0, 3 ) * src( 3, 2 );
      const double t37 = src( 1, 0 ) * src( 3, 1 );
      const double t38 = src( 0, 2 ) * src( 2, 3 );
      const double t40 = src( 0, 3 ) * src( 2, 2 );
      const double t42 = t14 * t15 - t14 * t17 - t19 * t20 + t19 * t22 +
          t24 * t25 - t24 * t27 - t29 * t15 + t29 * t17 +
          t32 * t33 - t32 * t35 - t37 * t38 + t37 * t40;
      const double t43 = src( 2, 0 ) * src( 0, 1 );
      const double t46 = src( 2, 0 ) * src( 1, 1 );
      const double t49 = src( 2, 0 ) * src( 3, 1 );
      const double t50 = src( 0, 2 ) * src( 1, 3 );
      const double t52 = src( 0, 3 ) * src( 1, 2 );
      const double t54 = src( 3, 0 ) * src( 0, 1 );
      const double t57 = src( 3, 0 ) * src( 1, 1 );
      const double t60 = src( 3, 0 ) * src( 2, 1 );
      const double t63 = t43 * t20 - t43 * t22 - t46 * t33 + t46 * t35 +
          t49 * t50 - t49 * t52 - t54 * t25 + t54 * t27 +
          t57 * t38 - t57 * t40 - t60 * t50 + t60 * t52;
      const double t65 = 1. / ( t42 + t63 );
      const double t71 = src( 0, 2 ) * src( 2, 1 );
      const double t73 = src( 0, 3 ) * src( 2, 1 );
      const double t75 = src( 0, 2 ) * src( 3, 1 );
      const double t77 = src( 0, 3 ) * src( 3, 1 );
      const double t81 = src( 0, 1 ) * src( 1, 2 );
      const double t83 = src( 0, 1 ) * src( 1, 3 );
      const double t85 = src( 0, 2 ) * src( 1, 1 );
      const double t87 = src( 0, 3 ) * src( 1, 1 );
      const double t101 = src( 1, 0 ) * src( 2, 2 );
      const double t103 = src( 1, 0 ) * src( 2, 3 );
      const double t105 = src( 2, 0 ) * src( 1, 2 );
      const double t107 = src( 2, 0 ) * src( 1, 3 );
      const double t109 = src( 3, 0 ) * src( 1, 2 );
      const double t111 = src( 3, 0 ) * src( 1, 3 );
      const double t115 = src( 0, 0 ) * src( 2, 2 );
      const double t117 = src( 0, 0 ) * src( 2, 3 );
      const double t119 = src( 2, 0 ) * src( 0, 2 );
      const double t121 = src( 2, 0 ) * src( 0, 3 );
      const double t123 = src( 3, 0 ) * src( 0, 2 );
      const double t125 = src( 3, 0 ) * src( 0, 3 );
      const double t129 = src( 0, 0 ) * src( 1, 2 );
      const double t131 = src( 0, 0 ) * src( 1, 3 );
      const double t133 = src( 1, 0 ) * src( 0, 2 );
      const double t135 = src( 1, 0 ) * src( 0, 3 );
      ( *this )( 0, 0 ) = ( src( 1, 1 ) * src( 2, 2 ) * src( 3, 3 ) - src( 1, 1 ) * src( 2, 3 ) * src( 3, 2 ) -
          src( 2, 1 ) * src( 1, 2 ) * src( 3, 3 ) + src( 2, 1 ) * src( 1, 3 ) * src( 3, 2 ) +
          src( 3, 1 ) * src( 1, 2 ) * src( 2, 3 ) - src( 3, 1 ) * src( 1, 3 ) * src( 2, 2 ) ) * t65;
      ( *this )( 0, 1 ) = -( src( 0, 1 ) * src( 2, 2 ) * src( 3, 3 ) - src( 0, 1 ) * src( 2, 3 ) * src( 3, 2 ) -
          t71 * src( 3, 3 ) + t73 * src( 3, 2 ) + t75 * src( 2, 3 ) - t77 * src( 2, 2 ) ) * t65;
      ( *this )( 0, 2 ) = ( t81 * src( 3, 3 ) - t83 * src( 3, 2 ) - t85 * src( 3, 3 ) + t87 * src( 3, 2 ) +
          t75 * src( 1, 3 ) - t77 * src( 1, 2 ) ) * t65;
      ( *this )( 0, 3 ) = -( t81 * src( 2, 3 ) - t83 * src( 2, 2 ) - t85 * src( 2, 3 ) + t87 * src( 2, 2 ) +
          t71 * src( 1, 3 ) - t73 * src( 1, 2 ) ) * t65;
      ( *this )( 1, 0 ) = -( t101 * src( 3, 3 ) - t103 * src( 3, 2 ) - t105 * src( 3, 3 ) + t107 * src( 3, 2 ) +
          t109 * src( 2, 3 ) - t111 * src( 2, 2 ) ) * t65;
      ( *this )( 1, 1 ) = ( t115 * src( 3, 3 ) - t117 * src( 3, 2 ) - t119 * src( 3, 3 ) + t121 * src( 3, 2 ) +
          t123 * src( 2, 3 ) - t125 * src( 2, 2 ) ) * t65;
      ( *this )( 1, 2 ) = -( t129 * src( 3, 3 ) - t131 * src( 3, 2 ) - t133 * src( 3, 3 ) + t135 * src( 3, 2 ) +
          t123 * src( 1, 3 ) - t125 * src( 1, 2 ) ) * t65;
      ( *this )( 1, 3 ) = ( t129 * src( 2, 3 ) - t131 * src( 2, 2 ) - t133 * src( 2, 3 ) + t135 * src( 2, 2 ) +
          t119 * src( 1, 3 ) - t121 * src( 1, 2 ) ) * t65;
      ( *this )( 2, 0 ) = ( t32 * src( 3, 3 ) - t103 * src( 3, 1 ) - t46 * src( 3, 3 ) + t107 * src( 3, 1 ) +
          t57 * src( 2, 3 ) - t111 * src( 2, 1 ) ) * t65;
      ( *this )( 2, 1 ) = -( t19 * src( 3, 3 ) - t117 * src( 3, 1 ) - t43 * src( 3, 3 ) + t121 * src( 3, 1 ) +
          t54 * src( 2, 3 ) - t125 * src( 2, 1 ) ) * t65;
      ( *this )( 2, 2 ) = ( t14 * src( 3, 3 ) - t131 * src( 3, 1 ) - t29 * src( 3, 3 ) + t135 * src( 3, 1 ) +
          t54 * src( 1, 3 ) - t125 * src( 1, 1 ) ) * t65;
      ( *this )( 2, 3 ) = -( t14 * src( 2, 3 ) - t131 * src( 2, 1 ) - t29 * src( 2, 3 ) + t135 * src( 2, 1 ) +
          t43 * src( 1, 3 ) - t121 * src( 1, 1 ) ) * t65;
      ( *this )( 3, 0 ) = -( t32 * src( 3, 2 ) - t101 * src( 3, 1 ) - t46 * src( 3, 2 ) + t105 * src( 3, 1 ) +
          t57 * src( 2, 2 ) - t109 * src( 2, 1 ) ) * t65;
      ( *this )( 3, 1 ) = ( t19 * src( 3, 2 ) - t115 * src( 3, 1 ) - t43 * src( 3, 2 ) + t119 * src( 3, 1 ) +
          t54 * src( 2, 2 ) - t123 * src( 2, 1 ) ) * t65;
      ( *this )( 3, 2 ) = -( t14 * src( 3, 2 ) - t129 * src( 3, 1 ) - t29 * src( 3, 2 ) + t133 * src( 3, 1 ) +
          t54 * src( 1, 2 ) - t123 * src( 1, 1 ) ) * t65;
      ( *this )( 3, 3 ) = ( t14 * src( 2, 2 ) - t129 * src( 2, 1 ) - t29 * src( 2, 2 ) + t133 * src( 2, 1 ) +
          t43 * src( 1, 2 ) - t119 * src( 1, 1 ) ) * t65;

      return;
    }

    default:
      {
      // Copy M in this
      this->m_values = src.m_values;

      // Declare workspace for permutations and scratch array
      lapack_int NN = integer_conversion<lapack_int>( this->get_nCols() );
      array1d<lapack_int> IPIV( NN );
      lapack_int INFO;
      array1d<double> INV_WORK( NN );

      // Call to LAPACK using LAPACKE
      // --- Compute LU factorization (LAPACK function DGETRF)
      INFO = LAPACKE_dgetrf( LAPACK_COL_MAJOR,
                             NN,
                             NN,
                             src.m_values.data(),
                             NN,
                             IPIV.data() );
      GEOS_ASSERT_MSG( INFO == 0, "LAPACKE_dgetrf: LU factorization failed" );

      // --- Invert (LAPACK function DGETRI)
      INFO = LAPACKE_dgetri( LAPACK_COL_MAJOR,
                             NN,
                             this->m_values.data(),
                             NN,
                             IPIV.data() );
      GEOS_ASSERT_MSG( INFO == 0, "LAPACKE_dgetri failed" );

//      std::vector<int> IPIV( NN );
//      int INFO;
//      std::vector<double> INV_WORK(NN);
//      // Call to BLAS using Epetra BLAS Wrapper Class (Epetra_BLAS)
//      Epetra_LAPACK Lapack;
//      // --- Compute LU factorization (LAPACK function DGETRF)
//      Lapack.GETRF(NN, NN, &(*this)(0,0), NN, &IPIV[0], &INFO);
//      ASSERT(INFO == 0, "LAPACK DGETRF error");
//      // --- Invert (LAPACK function DGETRI)
//      Lapack.GETRI(NN, &(*this)(0,0), NN, &IPIV[0], &INV_WORK[0], &NN, &INFO);
//      ASSERT(INFO == 0, "LAPACK DGETRI error");
      return;
    }
  }
}

// matrix-matrix sum (optional scaling)
void BlasMatrix::MatAdd(BlasMatrix const & A,
                        real64 const scalarA)
{
  unsigned int nRowA = A.get_nRows();
  unsigned int nColA = A.get_nCols();

  GEOS_ASSERT_MSG(A.get_nRows() == m_height,
                  "Matrix dimensions not compatible for sum");
  GEOS_ASSERT_MSG(A.get_nCols() == m_width,
                  "Matrix dimensions not compatible for sum");

  cblas_daxpy(static_cast<int>(m_values.size()),
              scalarA,
              A.m_values.data(),
              1,
              m_values.data(),
              1);
  return;

}

// matrix-matrix multiplication (optional scaling/accumulation)
void BlasMatrix::GEMM( BlasMatrix const & A,
                       BlasMatrix const & B,
                       real64 const scalarAB,
                       real64 const scalarThis )
{
  GEOS_ASSERT_MSG( A.get_nCols() == B.get_nRows(),
                   "Matrix dimensions not compatible for product" );
  GEOS_ASSERT_MSG( A.get_nRows() == m_height,
                   "Matrix dimensions not compatible for product" );
  GEOS_ASSERT_MSG( B.get_nCols() == m_width,
                   "Matrix dimensions not compatible for product" );

  int nRowAi = integer_conversion<int>( A.get_nRows() );
  int nColAi = integer_conversion<int>( A.get_nCols() );
  int nRowBi = integer_conversion<int>( B.get_nRows() );
  int nColBi = integer_conversion<int>( B.get_nCols() );

  cblas_dgemm( CblasColMajor,
               CblasNoTrans,
               CblasNoTrans,
               nRowAi,
               nColBi,
               nColAi,
               scalarAB,
               A.m_values.data(),
               nRowAi,
               B.m_values.data(),
               nRowBi,
               scalarThis,
               this->m_values.data(),
               nRowAi );
  return;

}

//
//// matrix-matrix multiplication (optional scaling/accumulation)
//void SerialDenseMatrix::MatMatMult(SerialDenseMatrix& A,
//                                   SerialDenseMatrix& B,
//                                   double scalarAB,
//                                   double scalarThis)
//{
//  unsigned int nRowA = A.get_nRows();
//  unsigned int nColA = A.get_nCols();
//  unsigned int nRowB = B.get_nRows();
//  unsigned int nColB = B.get_nCols();
//
//  ASSERT(nColA == nRowB, "Matrix dimensions not compatible for product");
//  ASSERT(nRowA == m_nRows, "Matrix dimensions not compatible for product");
//  ASSERT(nColB == m_nCols, "Matrix dimensions not compatible for product");
//
//#ifdef WITH_TRILINOS
//  // Call to BLAS using Epetra BLAS Wrapper Class (Epetra_BLAS)
//  int nRowAi = static_cast<int>(nRowA);
//  int nColAi = static_cast<int>(nColA);
//  int nRowBi = static_cast<int>(nRowB);
//  int nColBi = static_cast<int>(nColB);
//  Epetra_BLAS Blas;
//  Blas.GEMM('N', 'N', nRowAi, nColBi, nColAi,
//            scalarAB, &A(0,0), nRowAi, &B(0, 0), nRowBi,
//            scalarThis, &(*this)(0,0), m_nRows);
//  return;
//#else
//  // Loops are arranged without paying attention to accessing data in a
//  // contiguous way
//  for (unsigned int i = 0; i < m_nRows; ++i)
//    for (unsigned int j = 0; j < m_nCols; ++j)
//    {
//      double add_value = (abs(scalarThis) > 0.) ? (*this)(i,j) : 0.;
//      for (unsigned int k = 0; k < nColA; ++k)
//        add_value += scalarAB * A(i,k) * B(k,j);
//      (*this)(i,j) = add_value;
//    }
//#endif
//}
//
//// matrix-transpose(matrix) multiplication (optional scaling/accumulation)
//void SerialDenseMatrix::MatMatTMult(SerialDenseMatrix& A,
//                                    SerialDenseMatrix& B,
//                                    double scalarAB,
//                                    double scalarThis)
//{
//  unsigned int nRowA = A.get_nRows();
//  unsigned int nColA = A.get_nCols();
//  unsigned int nRowB = B.get_nRows();
//  unsigned int nColB = B.get_nCols();
//
//  ASSERT(nColA == nColB, "Matrix dimensions not compatible for product");
//  ASSERT(nRowA == m_nRows, "Matrix dimensions not compatible for product");
//  ASSERT(nRowB == m_nCols, "Matrix dimensions not compatible for product");
//
//#ifdef WITH_TRILINOS
//  // Call to BLAS using Epetra BLAS Wrapper Class (Epetra_BLAS)
//  int nRowAi = static_cast<int>(nRowA);
//  int nColAi = static_cast<int>(nColA);
//  int nRowBi = static_cast<int>(nRowB);
//  Epetra_BLAS Blas;
//  Blas.GEMM('N', 'T', nRowAi, nRowBi, nColAi,
//            scalarAB, &A(0,0), nRowAi, &B(0, 0), nRowBi,
//            scalarThis, &(*this)(0,0), m_nRows);
//  return;
//#else
//  // Loops are arranged without paying attention to accessing data in a
//  // contiguous way
//  for (unsigned int i = 0; i < m_nRows; ++i)
//    for (unsigned int j = 0; j < m_nCols; ++j)
//    {
//      double add_value = (abs(scalarThis) > 0.) ? (*this)(i,j) : 0.;
//      for (unsigned int k = 0; k < nColA; ++k)
//        add_value += scalarAB * A(i,k) * B(j,k);
//      (*this)(i,j) = add_value;
//    }
//#endif
//}
//
//// transpose(matrix)-matrix multiplication (optional scaling/accumulation)
//void SerialDenseMatrix::MatTMatMult(SerialDenseMatrix& A,
//                                    SerialDenseMatrix& B,
//                                    double scalarAB,
//                                    double scalarThis)
//{
//  unsigned int nRowA = A.get_nRows();
//  unsigned int nColA = A.get_nCols();
//  unsigned int nRowB = B.get_nRows();
//  unsigned int nColB = B.get_nCols();
//
//  ASSERT(nRowA == nRowB, "Matrix dimensions not compatible for product");
//  ASSERT(nColA == m_nRows, "Matrix dimensions not compatible for product");
//  ASSERT(nColB == m_nCols, "Matrix dimensions not compatible for product");
//
//#ifdef WITH_TRILINOS
//  // Call to BLAS using Epetra BLAS Wrapper Class (Epetra_BLAS)
//  int nRowAi = static_cast<int>(nRowA);
//  int nColAi = static_cast<int>(nColA);
//  int nRowBi = static_cast<int>(nRowB);
//  int nColBi = static_cast<int>(nColB);
//  Epetra_BLAS Blas;
//  Blas.GEMM('T', 'N', nColAi, nColBi, nRowAi,
//            scalarAB, &A(0,0), nRowAi, &B(0, 0), nRowBi,
//            scalarThis, &(*this)(0,0), m_nRows);
//  return;
//#else
//  // Loops are arranged without paying attention to accessing data in a
//  // contiguous way
//  for (unsigned int i = 0; i < m_nRows; ++i)
//    for (unsigned int j = 0; j < m_nCols; ++j)
//    {
//      double add_value = (abs(scalarThis) > 0.) ? (*this)(i,j) : 0.;
//      for (unsigned int k = 0; k < nRowA; ++k)
//        add_value += scalarAB * A(k,i) * B(k,j);
//      (*this)(i,j) = add_value;
//    }
//#endif
//}
//
//// transpose(Matrix)-transpose(matrix) multiplication
//// (optional scaling/accumulation)
//void SerialDenseMatrix::MatTMatTMult(SerialDenseMatrix& A,
//                                     SerialDenseMatrix& B,
//                                     double scalarAB,
//                                     double scalarThis)
//{
//  unsigned int nRowA = A.get_nRows();
//  unsigned int nColA = A.get_nCols();
//  unsigned int nRowB = B.get_nRows();
//  unsigned int nColB = B.get_nCols();
//
//  ASSERT(nRowA == nColB, "Matrix dimensions not compatible for product");
//  ASSERT(nColA == m_nRows, "Matrix dimensions not compatible for product");
//  ASSERT(nRowB == m_nCols, "Matrix dimensions not compatible for product");
//
//#ifdef WITH_TRILINOS
//  // Call to BLAS using Epetra BLAS Wrapper Class (Epetra_BLAS)
//  int nRowAi = static_cast<int>(nRowA);
//  int nColAi = static_cast<int>(nColA);
//  int nRowBi = static_cast<int>(nRowB);
//  Epetra_BLAS Blas;
//  Blas.GEMM('T', 'T', nColAi, nRowBi, nRowAi,
//            scalarAB, &A(0,0), nRowAi, &B(0, 0), nRowBi,
//            scalarThis, &(*this)(0,0), m_nRows);
//  return;
//#else
//  // Loops are arranged without paying attention to accessing data in a
//  // contiguous way
//  for (unsigned int i = 0; i  <m_nRows; ++i)
//    for (unsigned int j = 0; j < m_nCols; ++j)
//    {
//      double add_value = (abs(scalarThis) > 0.) ? (*this)(i,j) : 0.;
//      for (unsigned int k = 0; k < nRowA; ++k)
//        add_value += scalarAB * A(k,i) * B(j,k);
//      (*this)(i,j) = add_value;
//    }
//#endif
//}
//
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
  for( localIndex i = 0 ; i < m_height ; ++i )
  {
    for( localIndex j = 0 ; j < m_width ; ++j )
      printf( "%10.2e ", m_values[j * m_height + i] );
    printf( "\n" );
  }
}

} // end geosx namespace
