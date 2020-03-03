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
 * @file HypreMatrix.cpp
 */

#include "HypreMatrix.hpp"

#include "HYPRE.h"
#include "_hypre_IJ_mv.h"
#include "_hypre_parcsr_mv.h"
#include "HypreUtils.hpp"

#include <iomanip>

namespace geosx
{

// Helper function that performs the following sequence of IJMatrix
// call: Create, SetObjectType, Initialize.
static void initialize( MPI_Comm const & comm,
                        HYPRE_BigInt const & ilower,
                        HYPRE_BigInt const & iupper,
                        HYPRE_BigInt const & jlower,
                        HYPRE_BigInt const & jupper,
                        arraySlice1d< HYPRE_Int const > const & ncols,
                        HYPRE_IJMatrix & ij_matrix )
{
  GEOSX_LAI_CHECK_ERROR( HYPRE_IJMatrixCreate( comm,
                                               ilower,
                                               iupper,
                                               jlower,
                                               jupper,
                                               &ij_matrix ) );

  GEOSX_LAI_CHECK_ERROR( HYPRE_IJMatrixSetObjectType( ij_matrix, HYPRE_PARCSR ) );
  GEOSX_LAI_CHECK_ERROR( HYPRE_IJMatrixSetRowSizes( ij_matrix, ncols.data() ) );
  GEOSX_LAI_CHECK_ERROR( HYPRE_IJMatrixInitialize( ij_matrix ) );
}

HypreMatrix::HypreMatrix()
: LinearOperator(),
  MatrixBase(),
  m_ij_mat{},
  m_parcsr_mat{}
{}

HypreMatrix::HypreMatrix( HypreMatrix const & src )
: HypreMatrix()
{
#if 0
  GEOSX_LAI_ASSERT( src.ready() );

  HYPRE_BigInt ilower, iupper, jlower, jupper;

  GEOSX_LAI_CHECK_ERROR( HYPRE_IJMatrixGetLocalRange( src.m_ij_mat,
                                                      &ilower,
                                                      &iupper,
                                                      &jlower,
                                                      &jupper ) );

  // Get number of non-zeroes per row
  HYPRE_Int nrows = integer_conversion< HYPRE_Int >( iupper - ilower + 1 );
  array1d< HYPRE_BigInt > rows( nrows );
  array1d< HYPRE_Int > row_sizes( nrows );

  for( HYPRE_BigInt i = ilower ; i <= iupper ; ++i )
  {
    rows[i - ilower] = i;
  }

  GEOSX_LAI_CHECK_ERROR( HYPRE_IJMatrixGetRowCounts( src.unwrapped(),
                                                     nrows,
                                                     rows.data(),
                                                     row_sizes.data() ) );

  // Get number of non-zeroes (ncols) for row indeces in rows,
  // column indeces and values
  HYPRE_Int nnz = 0;
  for( HYPRE_Int i = 0 ; i < nrows ; ++i )
  {
    nnz += row_sizes[i];
  }

  array1d< HYPRE_BigInt > cols( nnz );
  array1d< HYPRE_Real > values( nnz );

  GEOSX_LAI_CHECK_ERROR( HYPRE_IJMatrixGetValues( src.m_ij_mat,
                                                  -nrows,
                                                  row_sizes.data(),
                                                  rows.data(),
                                                  cols.data(),
                                                  values.data() ) );

  initialize( src.getComm(),
              ilower,
              iupper,
              jlower,
              jupper,
              row_sizes,
              m_ij_mat );

  GEOSX_LAI_CHECK_ERROR( HYPRE_IJMatrixSetValues( m_ij_mat,
                                                  nrows,
                                                  row_sizes.data(),
                                                  rows.data(),
                                                  cols.data(),
                                                  values.data() ) );

  close();
#endif

  GEOSX_LAI_ASSERT( src.ready() );

  // Copy parcsr matrix
  HYPRE_ParCSRMatrix const dst_parcsr = hypre_ParCSRMatrixClone( src.m_parcsr_mat, 1 );

  // Create IJ layer (with matrix closed)
  parCSRtoIJ( dst_parcsr );
}

HypreMatrix::~HypreMatrix()
{
  reset();
}

void HypreMatrix::createWithGlobalSize( globalIndex const globalRows,
                                        globalIndex const globalCols,
                                        localIndex const maxEntriesPerRow,
                                        MPI_Comm const & comm )
{
  GEOSX_LAI_ASSERT( closed() );
  GEOSX_LAI_ASSERT_GE( globalRows, 0 );
  GEOSX_LAI_ASSERT_GE( globalCols, 0 );
  GEOSX_LAI_ASSERT_GE( maxEntriesPerRow, 0 );

  reset();

  HYPRE_Int const rank  = integer_conversion< HYPRE_Int >( MpiWrapper::Comm_rank( comm ) );
  HYPRE_Int const nproc = integer_conversion< HYPRE_Int >( MpiWrapper::Comm_size( comm ) );

  HYPRE_Int const localRowSize = integer_conversion< HYPRE_Int >( globalRows / nproc );
  HYPRE_Int const rowResidual = integer_conversion< HYPRE_Int >( globalRows % nproc );

  HYPRE_Int const localColSize = integer_conversion< HYPRE_Int >( globalCols / nproc );
  HYPRE_Int const colResidual = integer_conversion< HYPRE_Int >( globalCols % nproc );

  HYPRE_BigInt const ilower = rank * localRowSize + ( rank == 0 ? 0 : rowResidual );
  HYPRE_BigInt const iupper = ilower + localRowSize + ( rank == 0 ? rowResidual : 0 ) - 1;
  HYPRE_BigInt const jlower = rank * localColSize + ( rank == 0 ? 0 : colResidual );
  HYPRE_BigInt const jupper = jlower + localColSize + ( rank == 0 ? colResidual : 0 ) - 1;

  array1d< HYPRE_Int > row_sizes( integer_conversion<localIndex>( iupper - ilower + 1 ) );
  row_sizes = integer_conversion< HYPRE_Int >( maxEntriesPerRow );

  initialize( comm,
              ilower,
              iupper,
              jlower,
              jupper,
              row_sizes,
              m_ij_mat );
}

void HypreMatrix::createWithLocalSize( localIndex const localRows,
                                       localIndex const localCols,
                                       localIndex const maxEntriesPerRow,
                                       MPI_Comm const & comm )
{
  GEOSX_LAI_ASSERT_GE( localRows, 0 );
  GEOSX_LAI_ASSERT_GE( localCols, 0 );
  GEOSX_LAI_ASSERT_GE( maxEntriesPerRow, 0 );

  reset();

  HYPRE_BigInt const ilower = MpiWrapper::PrefixSum<HYPRE_BigInt>( localRows, comm ).first;
  HYPRE_BigInt const iupper = ilower + localRows - 1;
  HYPRE_BigInt const jlower = MpiWrapper::PrefixSum<HYPRE_BigInt>( localCols, comm ).first;
  HYPRE_BigInt const jupper = ilower + localCols - 1;

  array1d< HYPRE_Int > row_sizes( localRows );
  row_sizes = integer_conversion< HYPRE_Int >( maxEntriesPerRow );

  initialize( comm,
              ilower,
              iupper,
              jlower,
              jupper,
              row_sizes,
              m_ij_mat );
}

void HypreMatrix::set( real64 const value )
{
  GEOSX_LAI_ASSERT( ready() );
  open();
  GEOSX_LAI_CHECK_ERROR( HYPRE_IJMatrixSetConstantValues( m_ij_mat, static_cast< HYPRE_Real >( value ) ) );
  close();
}

void HypreMatrix::reset()
{
  MatrixBase::reset();
  if( m_ij_mat )
  {
    GEOSX_LAI_CHECK_ERROR( HYPRE_IJMatrixDestroy( m_ij_mat ) );
    m_ij_mat = nullptr;
    m_parcsr_mat = nullptr;
  }
}

void HypreMatrix::zero()
{
  GEOSX_LAI_ASSERT( ready() );
  open();
  GEOSX_LAI_CHECK_ERROR( hypre_IJMatrixSetConstantValuesParCSR( m_ij_mat, 0.0 ) );
  close();
}

void HypreMatrix::open()
{
  GEOSX_LAI_ASSERT( created() && closed() );
  GEOSX_LAI_CHECK_ERROR( HYPRE_IJMatrixInitialize( m_ij_mat ) );
  m_closed = false;
}

void HypreMatrix::close()
{
  GEOSX_LAI_ASSERT( !closed() );

  GEOSX_LAI_CHECK_ERROR( HYPRE_IJMatrixAssemble( m_ij_mat ) );

  // Get a reference to the constructed matrix object. Done only on the first
  // assembly call when the sparsity pattern of the matrix is defined.
  if( !m_assembled )
  {
    GEOSX_LAI_CHECK_ERROR( HYPRE_IJMatrixGetObject( m_ij_mat, (void * *) &m_parcsr_mat ) );
    // Compute column partitioning if needed
    if( hypre_IJMatrixRowPartitioning( m_ij_mat ) !=
        hypre_IJMatrixColPartitioning( m_ij_mat ) )
    {
      if( !hypre_ParCSRMatrixCommPkg( m_parcsr_mat ) )
      {
        GEOSX_LAI_CHECK_ERROR( hypre_MatvecCommPkgCreate( m_parcsr_mat ) );
      }
    }
  }

  m_closed = true;
  m_assembled = true;
}

bool HypreMatrix::created() const
{
  return m_ij_mat != nullptr;
}

void HypreMatrix::add( globalIndex const rowIndex,
                       globalIndex const colIndex,
                       real64 const value )
{
  GEOSX_LAI_ASSERT( modifiable() );

  HYPRE_Int ncols = 1;
  GEOSX_LAI_CHECK_ERROR( HYPRE_IJMatrixAddToValues( m_ij_mat,
                                                    1,
                                                    &ncols,
                                                    toHYPRE_BigInt( &rowIndex ),
                                                    toHYPRE_BigInt( &colIndex ),
                                                    &value ) );
}

void HypreMatrix::set( globalIndex const rowIndex,
                       globalIndex const colIndex,
                       real64 const value )
{
  GEOSX_LAI_ASSERT( modifiable() );
  GEOSX_LAI_ASSERT_GE( rowIndex, ilower() );
  GEOSX_LAI_ASSERT_GT( iupper(), rowIndex );

  HYPRE_Int ncols = 1;
  GEOSX_LAI_CHECK_ERROR( HYPRE_IJMatrixSetValues( m_ij_mat,
                                                  1,
                                                  &ncols,
                                                  toHYPRE_BigInt( &rowIndex ),
                                                  toHYPRE_BigInt( &colIndex ),
                                                  &value ) );

}

void HypreMatrix::insert( globalIndex const rowIndex,
                          globalIndex const colIndex,
                          real64 const value )
{
  GEOSX_LAI_ASSERT( insertable() );

  HYPRE_Int ncols = 1;
  GEOSX_LAI_CHECK_ERROR( HYPRE_IJMatrixAddToValues( m_ij_mat,
                                                    1,
                                                    &ncols,
                                                    toHYPRE_BigInt( &rowIndex ),
                                                    toHYPRE_BigInt( &colIndex ),
                                                    &value ) );
}

void HypreMatrix::add( globalIndex const rowIndex,
                       globalIndex const * colIndices,
                       real64 const * values,
                       localIndex size )
{
  GEOSX_LAI_ASSERT( modifiable() );

  HYPRE_Int ncols = integer_conversion< HYPRE_Int >( size );
  GEOSX_LAI_CHECK_ERROR( HYPRE_IJMatrixAddToValues( m_ij_mat,
                                                    1,
                                                    &ncols,
                                                    toHYPRE_BigInt( &rowIndex ),
                                                    toHYPRE_BigInt( colIndices ),
                                                    values ) );
}

void HypreMatrix::set( globalIndex const rowIndex,
                       globalIndex const * colIndices,
                       real64 const * values,
                       localIndex size )
{
  GEOSX_LAI_ASSERT( modifiable() );
  GEOSX_LAI_ASSERT_GE( rowIndex, ilower() );
  GEOSX_LAI_ASSERT_GT( iupper(), rowIndex );

  HYPRE_Int ncols = integer_conversion< HYPRE_Int >( size );
  GEOSX_LAI_CHECK_ERROR( HYPRE_IJMatrixSetValues( m_ij_mat,
                                                  1,
                                                  &ncols,
                                                  toHYPRE_BigInt( &rowIndex ),
                                                  toHYPRE_BigInt( colIndices ),
                                                  values ) );
}

void HypreMatrix::insert( globalIndex const rowIndex,
                          globalIndex const * colIndices,
                          real64 const * values,
                          localIndex size )
{
  GEOSX_LAI_ASSERT( insertable() );

  HYPRE_Int ncols = integer_conversion< HYPRE_Int >( size );
  GEOSX_LAI_CHECK_ERROR( HYPRE_IJMatrixAddToValues( m_ij_mat,
                                                    1,
                                                    &ncols,
                                                    toHYPRE_BigInt( &rowIndex ),
                                                    toHYPRE_BigInt( colIndices ),
                                                    values ) );
}

void HypreMatrix::add( globalIndex const rowIndex,
                       arraySlice1d< globalIndex const > const & colIndices,
                       arraySlice1d< real64 const > const & values )
{
  GEOSX_LAI_ASSERT( modifiable() );

  HYPRE_Int ncols = integer_conversion< HYPRE_Int >( colIndices.size() );
  GEOSX_LAI_CHECK_ERROR( HYPRE_IJMatrixAddToValues( m_ij_mat,
                                                    1,
                                                    &ncols,
                                                    toHYPRE_BigInt( &rowIndex ),
                                                    toHYPRE_BigInt( colIndices),
                                                    values ) );
}

void HypreMatrix::set( globalIndex const rowIndex,
                       arraySlice1d< globalIndex const > const & colIndices,
                       arraySlice1d< real64 const > const & values )
{
  GEOSX_LAI_ASSERT( modifiable() );
  GEOSX_LAI_ASSERT_GE( rowIndex, ilower() );
  GEOSX_LAI_ASSERT_GT( iupper(), rowIndex );

  HYPRE_Int ncols = integer_conversion< HYPRE_Int >( colIndices.size() );
  GEOSX_LAI_CHECK_ERROR( HYPRE_IJMatrixSetValues( m_ij_mat,
                                                  1,
                                                  &ncols,
                                                  toHYPRE_BigInt( &rowIndex ),
                                                  toHYPRE_BigInt( colIndices ),
                                                  values ) );
}

void HypreMatrix::insert( globalIndex const rowIndex,
                          arraySlice1d< globalIndex const > const & colIndices,
                          arraySlice1d< real64 const > const & values )
{
  GEOSX_LAI_ASSERT( insertable() );

  HYPRE_Int ncols = integer_conversion< HYPRE_Int >( colIndices.size() );
  GEOSX_LAI_CHECK_ERROR( HYPRE_IJMatrixAddToValues( m_ij_mat,
                                                    1,
                                                    &ncols,
                                                    toHYPRE_BigInt( &rowIndex ),
                                                    toHYPRE_BigInt( colIndices ),
                                                    values ) );
}

void HypreMatrix::add( arraySlice1d< globalIndex const > const & rowIndices,
                       arraySlice1d< globalIndex const > const & colIndices,
                       arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const & values )
{
  for( localIndex i = 0 ; i < rowIndices.size() ; ++i )
  {
    add( rowIndices[i], colIndices, values[i] );
  }
}

void HypreMatrix::set( arraySlice1d< globalIndex const > const & rowIndices,
                       arraySlice1d< globalIndex const > const & colIndices,
                       arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const & values )
{
  for( localIndex i = 0 ; i < integer_conversion< localIndex >( rowIndices.size() ) ; ++i )
  {
    set( rowIndices[i], colIndices, values[i] );
  }
}

void HypreMatrix::insert( arraySlice1d< globalIndex const > const & rowIndices,
                          arraySlice1d< globalIndex const > const & colIndices,
                          arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const & values )
{
  for( localIndex i = 0 ; i < rowIndices.size() ; ++i )
  {
    insert( rowIndices[i], colIndices, values[i] );
  }
}

template<typename T, int SRC_USD, int DST_USD>
static void convertArrayLayout( arraySlice2d<T const, SRC_USD> const & src,
                                arraySlice2d<T, DST_USD> const & dst )
{
  GEOSX_ASSERT( src.size( 0 ) == dst.size( 0 ) );
  GEOSX_ASSERT( src.size( 1 ) == dst.size( 1 ) );
  for( localIndex i = 0; i < src.size( 0 ); ++i )
  {
    for( localIndex j = 0; j < src.size( 1 ); ++j )
    {
      dst(i, j) = src(i, j);
    }
  }
}

void HypreMatrix::add( arraySlice1d< globalIndex const > const & rowIndices,
                       arraySlice1d< globalIndex const > const & colIndices,
                       arraySlice2d< real64 const, MatrixLayout::COL_MAJOR > const & values )
{
  array2d<real64, MatrixLayout::ROW_MAJOR_PERM> valuesRowMajor( rowIndices.size(), colIndices.size() );
  convertArrayLayout( values, valuesRowMajor.toSlice() );
  add( rowIndices, colIndices, valuesRowMajor );
}

void HypreMatrix::set( arraySlice1d< globalIndex const > const & rowIndices,
                       arraySlice1d< globalIndex const > const & colIndices,
                       arraySlice2d< real64 const, MatrixLayout::COL_MAJOR > const & values )
{
  array2d<real64, MatrixLayout::ROW_MAJOR_PERM> valuesRowMajor( rowIndices.size(), colIndices.size() );
  convertArrayLayout( values, valuesRowMajor.toSlice() );
  set( rowIndices, colIndices, valuesRowMajor );
}

void HypreMatrix::insert( arraySlice1d< globalIndex const > const & rowIndices,
                          arraySlice1d< globalIndex const > const & colIndices,
                          arraySlice2d< real64 const, MatrixLayout::COL_MAJOR > const & values )
{
  array2d<real64, MatrixLayout::ROW_MAJOR_PERM> valuesRowMajor( rowIndices.size(), colIndices.size() );
  convertArrayLayout( values, valuesRowMajor.toSlice() );
  insert( rowIndices, colIndices, valuesRowMajor );
}

void HypreMatrix::add( globalIndex const * rowIndices,
                       globalIndex const * colIndices,
                       real64 const * values,
                       localIndex const numRows,
                       localIndex const numCols )
{
  for( localIndex i = 0 ; i < numRows ; ++i )
  {
    add( rowIndices[i], colIndices, values + numCols * i, numCols );
  }
}

void HypreMatrix::set( globalIndex const * rowIndices,
                       globalIndex const * colIndices,
                       real64 const * values,
                       localIndex const numRows,
                       localIndex const numCols )
{
  for( localIndex i = 0 ; i < numRows ; ++i )
  {
    set( rowIndices[i], colIndices, values + numCols * i, numCols );
  }
}

void HypreMatrix::insert( globalIndex const * rowIndices,
                          globalIndex const * colIndices,
                          real64 const * values,
                          localIndex const numRows,
                          localIndex const numCols )
{
  for( localIndex i = 0 ; i < numRows ; ++i )
  {
    insert( rowIndices[i], colIndices, values + numCols * i, numCols );
  }
}

void HypreMatrix::apply( HypreVector const & src,
                         HypreVector & dst ) const
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT( src.ready() );
  GEOSX_LAI_ASSERT( dst.ready() );
  GEOSX_LAI_CHECK_ERROR( hypre_ParCSRMatrixMatvec( 1.0,
                                                   m_parcsr_mat,
                                                   src.unwrappedParVector(),
                                                   0.0,
                                                   dst.unwrappedParVector() ) );
}

void HypreMatrix::multiply( HypreMatrix const & src,
                            HypreMatrix & dst,
                            bool const closeResult ) const
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT( src.ready() );
  GEOSX_LAI_ASSERT_EQ( numGlobalCols(), src.numGlobalRows() );

  // Compute product
  HYPRE_ParCSRMatrix const dst_parcsr = hypre_ParMatmul( m_parcsr_mat, src.m_parcsr_mat );

  // Create IJ layer (with matrix closed)
  dst.parCSRtoIJ( dst_parcsr );

  // Reopen matrix if desired
  if( !closeResult )
  {
    dst.open();
  }
}

void HypreMatrix::leftMultiplyTranspose( HypreMatrix const & src,
                                         HypreMatrix & dst,
                                         bool const closeResult ) const
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT( src.ready() );
  GEOSX_LAI_ASSERT_EQ( numGlobalRows(), src.numGlobalRows() );

  // Compute product
  HYPRE_ParCSRMatrix const dst_parcsr = hypre_ParTMatmul( m_parcsr_mat, src.m_parcsr_mat );

  // Create IJ layer (with matrix closed)
  dst.parCSRtoIJ( dst_parcsr );

  // Reopen matrix if desired
  if( !closeResult )
  {
    dst.open();
  }

}

void HypreMatrix::rightMultiplyTranspose( HypreMatrix const & src,
                                          HypreMatrix & dst,
                                          bool const closeResult ) const
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT( src.ready() );
  GEOSX_LAI_ASSERT_EQ( numGlobalCols(), src.numGlobalCols() );

  // Transpose this
  HypreMatrix tmp;
  transpose( tmp );

  // Compute product
  src.multiply( tmp, dst );

  // Reopen matrix if desired
  if( !closeResult )
  {
    dst.open();
  }

}

void HypreMatrix::parCSRtoIJ( HYPRE_ParCSRMatrix const & parCSRMatrix )
{
  reset();
  m_closed = false;

  hypre_IJMatrix *ijmatrix;

  ijmatrix = hypre_CTAlloc( hypre_IJMatrix, 1, HYPRE_MEMORY_HOST );

  hypre_IJMatrixComm( ijmatrix ) = hypre_ParCSRMatrixComm( parCSRMatrix );

  hypre_IJMatrixObject( ijmatrix ) = parCSRMatrix;
  hypre_IJMatrixTranslator( ijmatrix ) = NULL;
  hypre_IJMatrixAssumedPart( ijmatrix ) = hypre_ParCSRMatrixAssumedPartition( parCSRMatrix );

  hypre_IJMatrixAssembleFlag( ijmatrix ) = 1;

  hypre_IJMatrixObjectType( ijmatrix ) = HYPRE_PARCSR;
#ifdef HYPRE_USING_OPENMP
  hypre_IJMatrixOMPFlag( ijmatrix ) = 1;
#else
  hypre_IJMatrixOMPFlag( ijmatrix ) = 0;
#endif
  hypre_IJMatrixPrintLevel( ijmatrix ) = 0;

  array1d< HYPRE_BigInt > info( 2 );
  if( MpiWrapper::Comm_rank( hypre_IJMatrixComm( ijmatrix ) ) == 0 )
  {
    info( 0 ) = hypre_ParCSRMatrixFirstRowIndex( parCSRMatrix );
    info( 1 ) = hypre_ParCSRMatrixFirstColDiag( parCSRMatrix );
  }
  MpiWrapper::bcast( info.data(), 2, 0, hypre_IJMatrixComm( ijmatrix ) );
  hypre_IJMatrixGlobalFirstRow( ijmatrix ) = info( 0 );
  hypre_IJMatrixGlobalFirstCol( ijmatrix ) = info( 1 );

  hypre_IJMatrixGlobalNumRows( ijmatrix ) = hypre_ParCSRMatrixGlobalNumRows( parCSRMatrix );
  hypre_IJMatrixGlobalNumCols( ijmatrix ) = hypre_ParCSRMatrixGlobalNumCols( parCSRMatrix );

  // Set row partitioning
  if( hypre_ParCSRMatrixOwnsRowStarts( parCSRMatrix ) )
  {
    hypre_IJMatrixRowPartitioning( ijmatrix ) = hypre_ParCSRMatrixRowStarts( parCSRMatrix );
  }
  else
  {
    HYPRE_BigInt *row_partitioning;
    HYPRE_BigInt *row_starts = hypre_ParCSRMatrixRowStarts( parCSRMatrix );
#ifdef HYPRE_NO_GLOBAL_PARTITION
    row_partitioning = hypre_CTAlloc( HYPRE_BigInt, 2, HYPRE_MEMORY_HOST );
    row_partitioning[0] = row_starts[0];
    row_partitioning[1] = row_starts[1];
#else
    GEOSX_ERROR( "HYPRE intended to be used only with HYPRE_NO_GLOBAL_PARTITION" )
#endif
    hypre_IJMatrixRowPartitioning( ijmatrix ) = row_partitioning;
    hypre_ParCSRMatrixRowStarts( parCSRMatrix ) = row_partitioning;
  }
  hypre_ParCSRMatrixOwnsRowStarts( parCSRMatrix ) = 0;

  if( hypre_IJMatrixGlobalNumRows( ijmatrix ) != hypre_IJMatrixGlobalNumCols( ijmatrix ) )
  {
    // Rectangular matrix
    // Set column partitioning
    if( hypre_ParCSRMatrixOwnsColStarts( parCSRMatrix ) )
    {
      hypre_IJMatrixColPartitioning( ijmatrix ) = hypre_ParCSRMatrixColStarts( parCSRMatrix );
    }
    else
    {
      HYPRE_BigInt *col_partitioning;
      HYPRE_BigInt *col_starts = hypre_ParCSRMatrixColStarts( parCSRMatrix );
#ifdef HYPRE_NO_GLOBAL_PARTITION
      col_partitioning = hypre_CTAlloc( HYPRE_BigInt, 2, HYPRE_MEMORY_HOST );
      col_partitioning[0] = col_starts[0];
      col_partitioning[1] = col_starts[1];
#else
      GEOSX_ERROR( "HYPRE intended to be used only with HYPRE_NO_GLOBAL_PARTITION" )
#endif
      hypre_IJMatrixColPartitioning( ijmatrix ) = col_partitioning;
      hypre_ParCSRMatrixColStarts( parCSRMatrix ) = col_partitioning;
    }
    hypre_ParCSRMatrixOwnsColStarts( parCSRMatrix ) = 0;
  }
  else
  {
    // Square matrix
    hypre_IJMatrixColPartitioning( ijmatrix ) = hypre_IJMatrixRowPartitioning( ijmatrix );
    hypre_ParCSRMatrixOwnsColStarts( parCSRMatrix ) = hypre_ParCSRMatrixOwnsRowStarts( parCSRMatrix );
  }

  m_ij_mat = (HYPRE_IJMatrix) ijmatrix;

  close();

}

void HypreMatrix::gemv( real64 const alpha,
                        HypreVector const & x,
                        real64 const beta,
                        HypreVector & y,
                        bool useTranspose ) const
{
  GEOSX_LAI_ASSERT( ready() );

  if( !useTranspose )
  {
    GEOSX_LAI_CHECK_ERROR( hypre_ParCSRMatrixMatvec( alpha,
                                                     m_parcsr_mat,
                                                     x.unwrappedParVector(),
                                                     beta,
                                                     y.unwrappedParVector() ) );
  }
  else
  {
    GEOSX_LAI_CHECK_ERROR( hypre_ParCSRMatrixMatvecT( alpha,
                                                      m_parcsr_mat,
                                                      x.unwrappedParVector(),
                                                      beta,
                                                      y.unwrappedParVector() ) );
  }
}

void HypreMatrix::scale( real64 const scalingFactor )
{
  GEOSX_LAI_ASSERT( ready() );

  hypre_CSRMatrix * const prt_diag_CSR = hypre_ParCSRMatrixDiag( m_parcsr_mat );
  HYPRE_Int const diag_nnz = hypre_CSRMatrixNumNonzeros( prt_diag_CSR );
  HYPRE_Real * const ptr_diag_data = hypre_CSRMatrixData( prt_diag_CSR );

  hypre_CSRMatrix * const prt_offdiag_CSR = hypre_ParCSRMatrixOffd( m_parcsr_mat );
  HYPRE_Int const offdiag_nnz = hypre_CSRMatrixNumNonzeros( prt_offdiag_CSR );
  HYPRE_Real * const ptr_offdiag_data = hypre_CSRMatrixData( prt_offdiag_CSR );

  HYPRE_Real alpha = static_cast< HYPRE_Real >( scalingFactor );
  for( HYPRE_Int i = 0 ; i < diag_nnz ; ++i )
  {
    ptr_diag_data[i] *= alpha;
  }
  for( HYPRE_Int i = 0 ; i < offdiag_nnz ; ++i )
  {
    ptr_offdiag_data[i] *= alpha;
  }
}

void HypreMatrix::leftScale( HypreVector const & vec )
{
  GEOSX_LAI_ASSERT( ready() );

  hypre_Vector * const ptr_vec = hypre_ParVectorLocalVector( vec.unwrappedParVector() );
  HYPRE_Real * ptr_vec_data = hypre_VectorData( ptr_vec );

  hypre_CSRMatrix * const prt_diag_CSR = hypre_ParCSRMatrixDiag( m_parcsr_mat );
  HYPRE_Real * const ptr_diag_data = hypre_CSRMatrixData( prt_diag_CSR );

  HYPRE_Int nrows = hypre_CSRMatrixNumRows( prt_diag_CSR );
  HYPRE_Int const * IA = hypre_CSRMatrixI( prt_diag_CSR );

  for( HYPRE_Int i = 0 ; i < nrows ; ++i )
  {
    for( HYPRE_Int j = IA[i] ; j < IA[i + 1] ; ++j )
    {
      ptr_diag_data[j] *= ptr_vec_data[i];
    }
  }

  hypre_CSRMatrix * const prt_offdiag_CSR = hypre_ParCSRMatrixOffd( m_parcsr_mat );
  HYPRE_Real * const ptr_offdiag_data = hypre_CSRMatrixData( prt_offdiag_CSR );

  nrows = hypre_CSRMatrixNumRows( prt_offdiag_CSR );
  IA = hypre_CSRMatrixI( prt_offdiag_CSR );

  for( HYPRE_Int i = 0 ; i < nrows ; ++i )
  {
    for( HYPRE_Int j = IA[i] ; j < IA[i + 1] ; ++j )
    {
      ptr_offdiag_data[j] *= ptr_vec_data[i];
    }
  }

}

localIndex HypreMatrix::maxRowLength() const
{
  GEOSX_LAI_ASSERT( assembled() );

  HYPRE_Int const nrows = integer_conversion< HYPRE_Int >( numLocalRows() );

  array1d< HYPRE_BigInt > rows( nrows );
  array1d< HYPRE_Int > ncols( nrows );

  for( HYPRE_Int i = 0 ; i < nrows ; ++i )
  {
    rows[i] = ilower() + integer_conversion< HYPRE_BigInt >( i );
  }

  GEOSX_LAI_CHECK_ERROR( HYPRE_IJMatrixGetRowCounts( m_ij_mat,
                                                     nrows,
                                                     rows.data(),
                                                     ncols.data() ) );

  HYPRE_Int const localMaxRowLength = *std::max_element( ncols.begin(), ncols.end() );
  return MpiWrapper::Max( localMaxRowLength, getComm() );
}

localIndex HypreMatrix::localRowLength( localIndex localRowIndex ) const
{
  return globalRowLength( getGlobalRowID( localRowIndex ) );
}

localIndex HypreMatrix::globalRowLength( globalIndex globalRowIndex ) const
{
  GEOSX_LAI_ASSERT( assembled() );

  HYPRE_BigInt row = integer_conversion< HYPRE_BigInt >( globalRowIndex );
  HYPRE_Int ncols;

  GEOSX_LAI_CHECK_ERROR( HYPRE_IJMatrixGetRowCounts( m_ij_mat,
                                                     1,
                                                     &row,
                                                     &ncols ) );

  return integer_conversion< localIndex >( ncols );
}

void HypreMatrix::getRowCopy( globalIndex globalRow,
                              arraySlice1d< globalIndex > const & colIndices,
                              arraySlice1d< real64 > const & values ) const
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT_GE( globalRow, ilower() );
  GEOSX_LAI_ASSERT_GT( iupper(), globalRow );

  HYPRE_Int numEntries;
  GEOSX_LAI_CHECK_ERROR( HYPRE_IJMatrixGetRowCounts( m_ij_mat,
                                                     1,
                                                     toHYPRE_BigInt( &globalRow ),
                                                     &numEntries ) );

  GEOSX_LAI_ASSERT_GE( colIndices.size(), numEntries );
  GEOSX_LAI_ASSERT_GE( values.size(), numEntries );

  GEOSX_LAI_CHECK_ERROR( hypre_IJMatrixGetValuesParCSR( m_ij_mat,
                                                        -1,
                                                        &numEntries,
                                                        toHYPRE_BigInt( &globalRow ),
                                                        toHYPRE_BigInt( colIndices ),
                                                        values ) );
}

real64 HypreMatrix::getDiagValue( globalIndex globalRow ) const
{
  GEOSX_LAI_ASSERT( assembled() );
  GEOSX_LAI_ASSERT_GE( globalRow, ilower());
  GEOSX_LAI_ASSERT_GT( iupper(), globalRow );

  // Get local row index
  HYPRE_Int const localRow = integer_conversion< HYPRE_Int >( getLocalRowID( globalRow ) );

  // Get diagonal block
  hypre_CSRMatrix const * const prt_CSR  = hypre_ParCSRMatrixDiag( m_parcsr_mat );
  HYPRE_Int const       * const IA       = hypre_CSRMatrixI( prt_CSR );
  HYPRE_Int const       * const JA       = hypre_CSRMatrixJ( prt_CSR );
  HYPRE_Real const      * const ptr_data = hypre_CSRMatrixData( prt_CSR );

  for( HYPRE_Int j = IA[localRow] ; j < IA[localRow + 1] ; ++j )
  {
    if( JA[j] == localRow )
    {
      return ptr_data[j];
    }
  }
  return 0.0;
}

void HypreMatrix::clearRow( globalIndex const globalRow,
                            real64 const diagValue )
{
  GEOSX_LAI_ASSERT( modifiable() );
  GEOSX_LAI_ASSERT_GE( globalRow, ilower() );
  GEOSX_LAI_ASSERT_GT( iupper(), globalRow );

  // Get local row index
  HYPRE_Int const localRow = integer_conversion< HYPRE_Int >( getLocalRowID( globalRow ) );

  // Clear row in diagonal block
  hypre_CSRMatrix * prt_CSR  = hypre_ParCSRMatrixDiag( m_parcsr_mat );
  HYPRE_Int *       IA       = hypre_CSRMatrixI( prt_CSR );
  HYPRE_Real *      ptr_data = hypre_CSRMatrixData( prt_CSR );

  for( HYPRE_Int j = IA[localRow] ; j < IA[localRow + 1] ; ++j )
  {
    ptr_data[j] = 0.0;
  }

  // Clear row in off-diagonal block
  prt_CSR  = hypre_ParCSRMatrixOffd( m_parcsr_mat );
  IA       = hypre_CSRMatrixI( prt_CSR );
  ptr_data = hypre_CSRMatrixData( prt_CSR );

  for( HYPRE_Int j = IA[localRow] ; j < IA[localRow + 1] ; ++j )
  {
    ptr_data[j] = 0.0;
  }

  // Set diagonal value
  if( std::fabs(diagValue) > 0.0 && numGlobalRows() == numGlobalCols() )
  {
    set( globalRow, globalRow, diagValue );
  }
}

HYPRE_IJMatrix const & HypreMatrix::unwrapped() const
{
  return m_ij_mat;
}

HYPRE_IJMatrix & HypreMatrix::unwrapped()
{
  return m_ij_mat;
}

HYPRE_ParCSRMatrix const & HypreMatrix::unwrappedParCSR() const
{
  return m_parcsr_mat;
}

HYPRE_ParCSRMatrix & HypreMatrix::unwrappedParCSR()
{
  return m_parcsr_mat;
}

localIndex HypreMatrix::getLocalRowID( globalIndex const index ) const
{
  GEOSX_LAI_ASSERT( created() );
  HYPRE_BigInt ilower, iupper, jlower, jupper;
  GEOSX_LAI_CHECK_ERROR( HYPRE_IJMatrixGetLocalRange( m_ij_mat,
                                                      &ilower, &iupper,
                                                      &jlower, &jupper ) );
  return (index >= ilower && index < iupper) ? integer_conversion< localIndex >( index - ilower ) : -1;
}

globalIndex HypreMatrix::getGlobalRowID( localIndex const index ) const
{
  GEOSX_LAI_ASSERT( created() );
  GEOSX_LAI_ASSERT_GE( index, 0 );
  GEOSX_LAI_ASSERT_GT( numLocalRows(), index );
  return ilower() + index;
}

globalIndex HypreMatrix::numGlobalRows() const
{
  GEOSX_LAI_ASSERT( created() );
  return hypre_IJMatrixGlobalNumRows( m_ij_mat );
}

globalIndex HypreMatrix::numGlobalCols() const
{
  GEOSX_LAI_ASSERT( created() );
  return hypre_IJMatrixGlobalNumCols( m_ij_mat );
}

localIndex HypreMatrix::numLocalRows() const
{
  GEOSX_LAI_ASSERT( created() );
  return integer_conversion< localIndex >( iupper() - ilower() );
}

localIndex HypreMatrix::numLocalCols() const
{
  GEOSX_LAI_ASSERT( created() );
  return integer_conversion< localIndex >( jupper() - jlower() );
}

globalIndex HypreMatrix::ilower() const
{
  GEOSX_LAI_ASSERT( created() );

  HYPRE_BigInt ilower, iupper, jlower, jupper;
  GEOSX_LAI_CHECK_ERROR( HYPRE_IJMatrixGetLocalRange( m_ij_mat,
                                                      &ilower, &iupper,
                                                      &jlower, &jupper ) );
  return integer_conversion< globalIndex >( ilower );
}

globalIndex HypreMatrix::iupper() const
{
  GEOSX_LAI_ASSERT( created() );

  HYPRE_BigInt ilower, iupper, jlower, jupper;
  GEOSX_LAI_CHECK_ERROR( HYPRE_IJMatrixGetLocalRange( m_ij_mat,
                                                      &ilower,
                                                      &iupper,
                                                      &jlower,
                                                      &jupper ) );
  return integer_conversion< globalIndex >( iupper + 1 );
}

globalIndex HypreMatrix::jlower() const
{
  GEOSX_LAI_ASSERT( created() );

  HYPRE_BigInt ilower, iupper, jlower, jupper;
  GEOSX_LAI_CHECK_ERROR( HYPRE_IJMatrixGetLocalRange( m_ij_mat,
                                                      &ilower,
                                                      &iupper,
                                                      &jlower,
                                                      &jupper ) );
  return integer_conversion< globalIndex >( jlower );
}

globalIndex HypreMatrix::jupper() const
{
  GEOSX_LAI_ASSERT( created() );

  HYPRE_BigInt ilower, iupper, jlower, jupper;
  GEOSX_LAI_CHECK_ERROR( HYPRE_IJMatrixGetLocalRange( m_ij_mat,
                                                      &ilower,
                                                      &iupper,
                                                      &jlower,
                                                      &jupper ) );
  return integer_conversion< globalIndex >( jupper + 1 );
}

localIndex HypreMatrix::numLocalNonzeros() const
{
  GEOSX_LAI_ASSERT( assembled() );

  hypre_CSRMatrix * prt_diag_CSR = hypre_ParCSRMatrixDiag( m_parcsr_mat );
  HYPRE_Int diag_nnz = hypre_CSRMatrixNumNonzeros( prt_diag_CSR );

  hypre_CSRMatrix * prt_offdiag_CSR = hypre_ParCSRMatrixOffd( m_parcsr_mat );
  HYPRE_Int offdiag_nnz = hypre_CSRMatrixNumNonzeros( prt_offdiag_CSR );

  return integer_conversion< localIndex >( diag_nnz + offdiag_nnz );
}

globalIndex HypreMatrix::numGlobalNonzeros() const
{
  return MpiWrapper::Sum( integer_conversion< globalIndex >( numLocalNonzeros() ), getComm() );
}

void HypreMatrix::print( std::ostream & os ) const
{
  GEOSX_LAI_ASSERT( ready() );

  int const this_mpi_process = MpiWrapper::Comm_rank( getComm() );
  int const n_mpi_process = MpiWrapper::Comm_size( getComm() );
  char str[77];

  if ( this_mpi_process == 0 )
  {
    os << "MPI_Process         GlobalRowID         GlobalColID                   Value" << std::endl;
  }

  for( int iRank = 0; iRank < n_mpi_process; iRank++ )
  {
	  MpiWrapper::Barrier( getComm() );
    if ( iRank == this_mpi_process )
    {
      globalIndex const firstRowID = ilower();
      globalIndex const firstDiagColID = jlower();

      hypre_CSRMatrix const * const prt_diag_CSR = hypre_ParCSRMatrixDiag( m_parcsr_mat );
      HYPRE_Int const * const diag_IA = hypre_CSRMatrixI( prt_diag_CSR );
      HYPRE_Int const * const diag_JA = hypre_CSRMatrixJ( prt_diag_CSR );
      HYPRE_Real const * const ptr_diag_data = hypre_CSRMatrixData( prt_diag_CSR );

      hypre_CSRMatrix const * const prt_offdiag_CSR = hypre_ParCSRMatrixOffd( m_parcsr_mat );
      HYPRE_Int const * const offdiag_IA = hypre_CSRMatrixI( prt_offdiag_CSR );
      HYPRE_Int const * const offdiag_JA = hypre_CSRMatrixJ( prt_offdiag_CSR );
      HYPRE_BigInt const * const col_map_offdiag = hypre_ParCSRMatrixColMapOffd( m_parcsr_mat );
      HYPRE_Real const * const ptr_offdiag_data = hypre_CSRMatrixData( prt_offdiag_CSR );

      for( HYPRE_Int i = 0 ; i < hypre_CSRMatrixNumRows( prt_diag_CSR ); ++i )
      {
        for( HYPRE_Int j = diag_IA[i] ; j < diag_IA[i + 1] ; ++j )
    	  {

          sprintf( str,
                   "%11i%20lli%20lli%24.10e\n",
                   iRank,
                   firstRowID + integer_conversion<globalIndex>(i),
                   firstDiagColID + integer_conversion<globalIndex>( diag_JA[j] ),
                   ptr_diag_data[j] );
          os << str;
        }
        for( HYPRE_Int j = offdiag_IA[i] ; j < offdiag_IA[i + 1] ; ++j )
    	  {
          sprintf( str,
                   "%11i%20lli%20lli%24.10e\n",
                   iRank,
                   firstRowID + integer_conversion<globalIndex>(i),
                   col_map_offdiag[ offdiag_JA[j] ],
                   ptr_offdiag_data[j] );
          os << str;
        }
      }
    }
  }
}

void HypreMatrix::write( string const & filename,
                         LAIOutputFormat const format ) const
{
  GEOSX_LAI_ASSERT( ready() );

  switch( format )
  {
    case LAIOutputFormat::NATIVE_ASCII:
    {
      GEOSX_LAI_CHECK_ERROR( hypre_ParCSRMatrixPrintIJ( m_parcsr_mat, 1, 1, filename.c_str() ) );
      break;
    }

    default:
      GEOSX_ERROR( "Unsupported matrix output format" );
  }
}

real64 HypreMatrix::norm1() const
{
  GEOSX_LAI_ASSERT( ready() );

  HypreMatrix matT;
  transpose( matT );
  return matT.normInf();
}

real64 HypreMatrix::normInf() const
{
  GEOSX_LAI_ASSERT( ready() );

  hypre_CSRMatrix const * const prt_diag_CSR = hypre_ParCSRMatrixDiag( m_parcsr_mat );
  HYPRE_Real const * const ptr_diag_data = hypre_CSRMatrixData( prt_diag_CSR );

  HYPRE_Int nrows = hypre_CSRMatrixNumRows( prt_diag_CSR );
  HYPRE_Int const * IA = hypre_CSRMatrixI( prt_diag_CSR );

  array1d< HYPRE_Real > row_abs_sum( nrows );
  row_abs_sum = 0.0;

  for( HYPRE_Int i = 0 ; i < nrows ; ++i )
  {
    for( HYPRE_Int j = IA[i] ; j < IA[i + 1] ; ++j )
    {
      row_abs_sum[i] += std::fabs( ptr_diag_data[j] );
    }
  }

  hypre_CSRMatrix const * const prt_offdiag_CSR = hypre_ParCSRMatrixOffd( m_parcsr_mat );
  HYPRE_Real const * const ptr_offdiag_data = hypre_CSRMatrixData( prt_offdiag_CSR );

  nrows = hypre_CSRMatrixNumRows( prt_offdiag_CSR );
  IA = hypre_CSRMatrixI( prt_offdiag_CSR );

  for( HYPRE_Int i = 0 ; i < nrows ; ++i )
  {
    for( HYPRE_Int j = IA[i] ; j < IA[i + 1] ; ++j )
    {
      row_abs_sum[i] += std::fabs( ptr_offdiag_data[j] );
    }
  }

  real64 const local_norm = *std::max_element( row_abs_sum.begin(), row_abs_sum.end() );
  return MpiWrapper::Max( local_norm, getComm() );

}

real64 HypreMatrix::normFrobenius() const
{
  GEOSX_LAI_ASSERT( ready() );
  return hypre_ParCSRMatrixFnorm( m_parcsr_mat );
}

void HypreMatrix::rightScale( HypreVector const & GEOSX_UNUSED_PARAM( vec ) )
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_ERROR( "Not implemented" );
}

void HypreMatrix::leftRightScale( HypreVector const & GEOSX_UNUSED_PARAM( vecLeft ),
                                  HypreVector const & GEOSX_UNUSED_PARAM( vecRight ) )
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_ERROR( "Not implemented" );
}

void HypreMatrix::transpose( HypreMatrix & dst ) const
{
  GEOSX_LAI_ASSERT( ready() );

  // Transpose this->m_parcsr_mat
  HYPRE_ParCSRMatrix dst_parcsr;
  GEOSX_LAI_CHECK_ERROR( hypre_ParCSRMatrixTranspose( m_parcsr_mat, &dst_parcsr, 1 ) );

  // Create IJ layer (with matrix closed)
  dst.parCSRtoIJ( dst_parcsr );
}

MPI_Comm HypreMatrix::getComm() const
{
  GEOSX_LAI_ASSERT( created() );
  return hypre_IJMatrixComm( m_ij_mat );
}

}// end namespace geosx

