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

// BEGIN_RST_NARRATIVE HypreMatrix.rst
// ==============================
// Epetra-based Matrix Object
// ==============================
// This class contains the ParallelMatrix wrappers based on Epetra_Crs Objects.
// The class contains a unique pointer to an Epetra_CrsMatrix as well as constructors,
// functions and accessors for Epetra objects.

// Include the corresponding header file.
#include "HypreMatrix.hpp"

// Include required Hypre headers
#include "HYPRE.h"
#include "_hypre_IJ_mv.h"
#include "_hypre_parcsr_mv.h"
#include "HypreUtils.hpp"

// Put everything under the geosx namespace.
namespace geosx
{

// Helper function that performs the following sequence of IJMatrix
// call: Create, SetObjectType, Initialize.
static void initialize( MPI_Comm const & comm,
                        HYPRE_BigInt const & ilower,
                        HYPRE_BigInt const & iupper,
                        HYPRE_BigInt const & jlower,
                        HYPRE_BigInt const & jupper,
                        array1d< HYPRE_Int > const & ncols,
                        HYPRE_IJMatrix & ij_matrix )
{
  HYPRE_IJMatrixCreate( comm,
                        ilower,
                        iupper,
                        jlower,
                        jupper,
                        &ij_matrix );
  HYPRE_IJMatrixSetObjectType( ij_matrix,
                               HYPRE_PARCSR );

  HYPRE_IJMatrixSetRowSizes( ij_matrix,
                             ncols.data() );

  HYPRE_IJMatrixInitialize( ij_matrix );
}

// ----------------------------
// Constructors
// ----------------------------

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Create an empty matrix (meant to be used for declaration)
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""

HypreMatrix::HypreMatrix()
{}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Copy constructor
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""

HypreMatrix::HypreMatrix( HypreMatrix const & src )
{
  GEOSX_ASSERT_MSG( src.unwrappedPointer() != nullptr,
                    "Input matrix looks empty" );
  GEOSX_ASSERT_MSG( src.isAssembled(),
                    "Input matrix hasn't been closed before copy" );

  this->reset();

  HYPRE_BigInt ilower, iupper, jlower, jupper;

  HYPRE_IJMatrixGetLocalRange( src.m_ij_mat,
                               &ilower,
                               &iupper,
                               &jlower,
                               &jupper );

  // Get number of non-zeroes per row
  HYPRE_Int nrows = integer_conversion< HYPRE_Int >( iupper - ilower + 1 );
  array1d< HYPRE_BigInt > rows( nrows );
  array1d< HYPRE_Int > row_sizes( nrows );

  for( HYPRE_BigInt i = ilower ; i <= iupper ; ++i )
  {
    rows[i - ilower] = i;
  }

  HYPRE_IJMatrixGetRowCounts( *src.unwrappedPointer(),
                              nrows,
                              rows.data(),
                              row_sizes.data() );

  // Get number of non-zeroes (ncols) for row indeces in rows,
  // column indeces and values
  HYPRE_Int nnz = 0;
  for( HYPRE_Int i = 0 ; i < nrows ; ++i )
  {
    nnz += row_sizes[i];
  }

  array1d< HYPRE_BigInt > cols( nnz );
  array1d< HYPRE_Real > values( nnz );

  HYPRE_IJMatrixGetValues( src.m_ij_mat,
                           -nrows,
                           row_sizes.data(),
                           rows.data(),
                           cols.data(),
                           values.data() );

  initialize( hypre_IJMatrixComm( src.m_ij_mat ),
              ilower,
              iupper,
              jlower,
              jupper,
              row_sizes,
              m_ij_mat );

  HYPRE_IJMatrixSetValues( m_ij_mat,
                           nrows,
                           row_sizes.data(),
                           rows.data(),
                           cols.data(),
                           values.data() );

  this->close();
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Destructor
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
HypreMatrix::~HypreMatrix()
{
  this->reset();
}

// -----------------------------
// Create
// -----------------------------
// Allocate matrix (prepare to be filled with data).

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Create a matrix from an HYPRE_IJMatrix.
// """""""""""""""""""""""""""""""""""""""""""""""

void HypreMatrix::create()
{
  //TODO:
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Create a matrix from number of elements
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""

// Note: When the global size is provided a near-even distribution
// of elements across processors is produced. All processors get the
// same number of rows/cols, except proc 0 which gets any remainder
// elements necessary when the number of processors does not divide
// evenly into the vector length.

void HypreMatrix::createWithGlobalSize( globalIndex const globalSize,
                                        localIndex const maxEntriesPerRow,
                                        MPI_Comm const & comm )
{
  this->reset();
  // just call general version
  this->createWithGlobalSize( globalSize,
                              globalSize,
                              maxEntriesPerRow,
                              comm );

}

void HypreMatrix::createWithGlobalSize( globalIndex const globalRows,
                                        globalIndex const globalCols,
                                        localIndex const maxEntriesPerRow,
                                        MPI_Comm const & comm )
{
  this->reset();

  int this_mpi_process = MpiWrapper::Comm_rank( comm );
  int n_mpi_process = MpiWrapper::Comm_size( comm );

  HYPRE_Int localRowSize = integer_conversion< HYPRE_Int >( globalRows / n_mpi_process );
  HYPRE_Int rowResidual = integer_conversion< HYPRE_Int >( globalRows % n_mpi_process );

  GEOSX_ASSERT_MSG( localRowSize >= 0,
                    "localRowSize is is lower than 0" );

  HYPRE_Int localColSize = integer_conversion< HYPRE_Int >( globalCols / n_mpi_process );
  HYPRE_Int colResidual = integer_conversion< HYPRE_Int >( globalCols % n_mpi_process );

  GEOSX_ASSERT_MSG( localColSize >= 0,
                    "localColSize is is lower than 0" );

  HYPRE_BigInt ilower, iupper, jlower, jupper;

  if( this_mpi_process == 0 )
  {
    ilower = 0;
    localRowSize = localRowSize + rowResidual;
    iupper = integer_conversion< HYPRE_BigInt >( localRowSize - 1 );
    jlower = 0;
    localColSize = localColSize + colResidual;
    jupper =  integer_conversion< HYPRE_BigInt >( localColSize - 1 );
  }
  else
  {
    ilower = integer_conversion< HYPRE_BigInt >( this_mpi_process * localRowSize + rowResidual );
    iupper = ilower + integer_conversion< HYPRE_BigInt >( localRowSize ) - 1;
    jlower = integer_conversion< HYPRE_BigInt >( this_mpi_process * localColSize + colResidual );
    jupper = jlower + integer_conversion< HYPRE_BigInt >( localColSize ) - 1;
  }

  array1d< HYPRE_Int > row_sizes( localRowSize );
  row_sizes = integer_conversion< HYPRE_Int >( integer_conversion< HYPRE_Int >( maxEntriesPerRow ) );

  initialize( comm,
              ilower,
              iupper,
              jlower,
              jupper,
              row_sizes,
              m_ij_mat );
}

void HypreMatrix::createWithLocalSize( localIndex const localSize,
                                       localIndex const maxEntriesPerRow,
                                       MPI_Comm const & comm )
{
  this->reset();
  // just call general version
  this->createWithLocalSize( localSize,
                             localSize,
                             maxEntriesPerRow,
                             comm );
}

void HypreMatrix::createWithLocalSize( localIndex const localRows,
                                       localIndex const localCols,
                                       localIndex const maxEntriesPerRow,
                                       MPI_Comm const & comm )
{
  GEOSX_ASSERT_MSG( localRows >= 0,
                    "local rows are lower than 0" );
  GEOSX_ASSERT_MSG( localCols >= 0,
                    "local columns are lower than 0" );
  this->reset();

  int this_mpi_process = MpiWrapper::Comm_rank( comm );
  int n_mpi_process = MpiWrapper::Comm_size( comm );

  array1d< HYPRE_Int > localSizeArray( n_mpi_process * 2 );
  array1d< HYPRE_Int > tmp( 2 );
  tmp[0] = integer_conversion< HYPRE_Int >( localRows );
  tmp[1] = integer_conversion< HYPRE_Int >( localCols );

  MpiWrapper::allGather( tmp.toViewConst(),
                         localSizeArray );

  HYPRE_BigInt ilower, iupper, jlower, jupper;

  ilower = 0;
  jlower = 0;
  for( int i = 0 ; i < this_mpi_process * 2 ; i += 2 )
  {
    ilower += integer_conversion< HYPRE_BigInt >( localSizeArray[i] );
    jlower += integer_conversion< HYPRE_BigInt >( localSizeArray[i + 1] );
  }
  iupper = ilower + integer_conversion< HYPRE_BigInt >( localRows ) - 1;
  jupper = jlower + integer_conversion< HYPRE_BigInt >( localCols ) - 1;

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

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Reinitialize.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Keeps the map and graph but sets all values to user-defined value.
void HypreMatrix::set( real64 const value )
{
  if( m_is_ready_to_use )
  {
    this->open();
    m_is_ready_to_use = false;
  }
  HYPRE_IJMatrixSetConstantValues( m_ij_mat, static_cast< HYPRE_Real >( value ) );
}


/**
 * @brief Reset the object
 *
 */
void HypreMatrix::reset()
{
  if( m_ij_mat )
  {
    HYPRE_IJMatrixDestroy( m_ij_mat );
    m_ij_mat = nullptr;
    m_parcsr_mat = nullptr;
  }
  m_is_pattern_fixed = false;
  m_is_ready_to_use = false;
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Reinitialize.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Keeps the map and graph but sets all values to 0. If the matrix is closed
// it will automatically re-open it.

void HypreMatrix::zero()
{
  if( m_is_ready_to_use )
  {
    this->open();
    m_is_ready_to_use = false;
  }
  hypre_IJMatrixSetConstantValuesParCSR( m_ij_mat, 0.0 );
  this->close();
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Open
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Empty open function (implemented for HYPRE compatibility).

void HypreMatrix::open()
{
  if( m_is_ready_to_use )
  {
    HYPRE_IJMatrixInitialize( m_ij_mat );
    m_is_ready_to_use = false;
  }
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Close
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Fix the sparsity pattern, make the data contiguous in memory and
// optimize storage.

void HypreMatrix::close()
{

  if( !m_is_ready_to_use )
  {
    HYPRE_IJMatrixAssemble( m_ij_mat );

    // Get a reference to the constructed matrix object. Done only on the first
    // assembly call when the sparsity pattern of the matrix is defined.
    if( !m_is_pattern_fixed )
    {
      HYPRE_IJMatrixGetObject( m_ij_mat, (void * *) &m_parcsr_mat );
      m_is_pattern_fixed = true;

      // Compute column partitioning if needed
      if( hypre_IJMatrixRowPartitioning( m_ij_mat ) !=
          hypre_IJMatrixColPartitioning( m_ij_mat ) )
      {
        if( !hypre_ParCSRMatrixCommPkg( m_parcsr_mat ))
        {
          hypre_MatvecCommPkgCreate( m_parcsr_mat );
        }
      }
    }

    m_is_ready_to_use = true;
  }
}

// -------------------------
// Add/Set
// -------------------------

// 1x1

void HypreMatrix::add( globalIndex const rowIndex,
                       globalIndex const colIndex,
                       real64 const value )
{
  GEOSX_ERROR_IF( m_ij_mat == nullptr,
                  "matrix appears to be empty (not created)" );
  HYPRE_Int ncols = 1;
  HYPRE_IJMatrixAddToValues( m_ij_mat,
                             1,
                             &ncols,
                             toHYPRE_BigInt( &rowIndex ),
                             toHYPRE_BigInt( &colIndex ),
                             toHYPRE_Real( &value ) );
}

void HypreMatrix::set( globalIndex const rowIndex,
                       globalIndex const colIndex,
                       real64 const value )
{
  GEOSX_ASSERT_MSG( m_ij_mat != nullptr,
                    "matrix appears to be empty (not created)" );
  GEOSX_ASSERT_MSG( this->getLocalRowID( rowIndex ) >= 0,
                    "HypreMatrix, it is not possible to set values on other processors" );

  HYPRE_Int ncols = 1;
  HYPRE_IJMatrixSetValues( m_ij_mat,
                           1,
                           &ncols,
                           toHYPRE_BigInt( &rowIndex ),
                           toHYPRE_BigInt( &colIndex ),
                           toHYPRE_Real( &value ) );

}

void HypreMatrix::insert( globalIndex const rowIndex,
                          globalIndex const colIndex,
                          real64 const value )
{
  this->add( rowIndex,
             colIndex,
             value );
}

// 1xN c-style

void HypreMatrix::add( globalIndex const rowIndex,
                       globalIndex const * colIndices,
                       real64 const * values,
                       localIndex size )
{
  GEOSX_ASSERT_MSG( m_ij_mat != nullptr,
                    "matrix appears to be empty (not created)" );
  HYPRE_Int ncols = integer_conversion< HYPRE_Int >( size );
  HYPRE_IJMatrixAddToValues( m_ij_mat,
                             1,
                             &ncols,
                             toHYPRE_BigInt( &rowIndex ),
                             toHYPRE_BigInt( colIndices ),
                             toHYPRE_Real( values ) );
}

void HypreMatrix::set( globalIndex const rowIndex,
                       globalIndex const * colIndices,
                       real64 const * values,
                       localIndex size )
{
  GEOSX_ASSERT_MSG( m_ij_mat != nullptr,
                    "matrix appears to be empty (not created)" );
  GEOSX_ASSERT_MSG( this->getLocalRowID( rowIndex ) >= 0,
                    "HypreMatrix, it is not possible to set values on other processors" );

  HYPRE_Int ncols = integer_conversion< HYPRE_Int >( size );
  HYPRE_IJMatrixSetValues( m_ij_mat,
                           1,
                           &ncols,
                           toHYPRE_BigInt( &rowIndex ),
                           toHYPRE_BigInt( colIndices ),
                           toHYPRE_Real( values ) );
}

void HypreMatrix::insert( globalIndex const rowIndex,
                          globalIndex const * colIndices,
                          real64 const * values,
                          localIndex size )
{
  this->add( rowIndex,
             colIndices,
             values,
             size );
}

// 1xN array1d style

void HypreMatrix::add( globalIndex const rowIndex,
                       array1d< globalIndex > const & colIndices,
                       array1d< real64 > const & values )
{
  GEOSX_ASSERT_MSG( m_ij_mat != nullptr,
                    "matrix appears to be empty (not created)" );
  HYPRE_Int ncols = integer_conversion< HYPRE_Int >( colIndices.size() );
  HYPRE_IJMatrixAddToValues( m_ij_mat,
                             1,
                             &ncols,
                             toHYPRE_BigInt( &rowIndex ),
                             toHYPRE_BigInt( colIndices.data() ),
                             toHYPRE_Real( values.data() ) );
}

void HypreMatrix::set( globalIndex const rowIndex,
                       array1d< globalIndex > const & colIndices,
                       array1d< real64 > const & values )
{
  GEOSX_ASSERT_MSG( m_ij_mat != nullptr,
                    "matrix appears to be empty (not created)" );
  GEOSX_ASSERT_MSG( this->getLocalRowID( rowIndex ) >= 0,
                    "HypreMatrix, it is not possible to set values on other processors" );

  HYPRE_Int ncols = integer_conversion< HYPRE_Int >( colIndices.size() );
  HYPRE_IJMatrixSetValues( m_ij_mat,
                           1,
                           &ncols,
                           toHYPRE_BigInt( &rowIndex ),
                           toHYPRE_BigInt( colIndices.data() ),
                           toHYPRE_Real( values.data() ) );
}

void HypreMatrix::insert( globalIndex const rowIndex,
                          array1d< globalIndex > const & colIndices,
                          array1d< real64 > const & values )
{
  this->add( rowIndex,
             colIndices,
             values );
}

//// MxN array2d style

void HypreMatrix::add( array1d< globalIndex > const & rowIndices,
                       array1d< globalIndex > const & colIndices,
                       array2d< real64 > const & values )
{
  GEOSX_ASSERT_MSG( m_ij_mat != nullptr,
                    "matrix appears to be empty (not created)" );
  localIndex nCols = integer_conversion< localIndex >( colIndices.size() );
  for( localIndex i = 0 ; i < integer_conversion< localIndex >( rowIndices.size() ) ; ++i )
  {
    this->add( rowIndices[i],
               colIndices.data(),
               values.begin() + nCols * i,
               nCols );
  }
}

void HypreMatrix::set( array1d< globalIndex > const & rowIndices,
                       array1d< globalIndex > const & colIndices,
                       array2d< real64 > const & values )
{
  GEOSX_ASSERT_MSG( m_ij_mat != nullptr,
                    "matrix appears to be empty (not created)" );
  GEOSX_ASSERT_MSG( this->getLocalRowID( *std::min_element( rowIndices.begin(), rowIndices.end() ) ) >= 0 &&
                    this->getLocalRowID( *std::max_element( rowIndices.begin(), rowIndices.end() ) ) >= 0,
                    "HypreMatrix, it is not possible to set values on other processors" );

  localIndex nCols = integer_conversion< localIndex >( colIndices.size() );
  for( localIndex i = 0 ; i < integer_conversion< localIndex >( rowIndices.size() ) ; ++i )
  {
    this->set( rowIndices[i],
               colIndices.data(),
               values.begin() + nCols * i,
               nCols );
  }
}

void HypreMatrix::insert( array1d< globalIndex > const & rowIndices,
                          array1d< globalIndex > const & colIndices,
                          array2d< real64 > const & values )
{
  this->add( rowIndices,
             colIndices,
             values );
}

void HypreMatrix::add( globalIndex const * rowIndices,
                       globalIndex const * colIndices,
                       real64 const * values,
                       localIndex const numRows,
                       localIndex const numCols )
{
  for( localIndex i = 0 ; i < numRows ; ++i )
  {
    this->add( rowIndices[i],
               colIndices,
               values + numCols * i,
               numCols );
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
    this->set( rowIndices[i],
               colIndices,
               values + numCols * i,
               numCols );
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
    this->add( rowIndices[i],
               colIndices,
               values + numCols * i,
               numCols );
  }
}

// -------------------------
// Linear Algebra
// -------------------------

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Matrix/vector multiplication
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Perform the matrix-vector product A*src = dst.

void HypreMatrix::multiply( HypreVector const & src,
                            HypreVector & dst ) const
{
  hypre_ParCSRMatrixMatvec( 1.0,
                            m_parcsr_mat,
                            *src.getHypreParVectorPointer(),
                            0.0,
                            *dst.getHypreParVectorPointer() );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Matrix/matrix multiplication
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Perform the matrix-matrix product this * src = dst.
void HypreMatrix::multiply( HypreMatrix const & src,
                            HypreMatrix & dst,
                            bool const closeResult ) const
{
  // Error check
  GEOSX_ASSERT_MSG( this->globalCols() == src.globalRows(),
                    "Incompatible matrix dimensions" );

  // Compute product
  HYPRE_ParCSRMatrix dst_parcsr;
  dst_parcsr = hypre_ParMatmul( m_parcsr_mat,
                                src.m_parcsr_mat );

  // Create IJ layer (with matrix closed)
  dst.parCSRtoIJ( dst_parcsr );

  // Reopen matrix if desired
  if( !closeResult )
  {
    dst.open();
  }
}

// Perform the matrix-matrix product transpose(this) * src = dst.
void HypreMatrix::leftMultiplyTranspose( HypreMatrix const & src,
                                         HypreMatrix & dst,
                                         bool const closeResult ) const
{
  // Error check
  GEOSX_ASSERT_MSG( this->globalRows() == src.globalRows(),
                    "Incompatible matrix dimensions" );

  // Compute product
  HYPRE_ParCSRMatrix dst_parcsr;

  dst_parcsr = hypre_ParTMatmul( m_parcsr_mat,
                                 src.m_parcsr_mat );

  // Create IJ layer (with matrix closed)
  dst.parCSRtoIJ( dst_parcsr );

  // Reopen matrix if desired
  if( !closeResult )
  {
    dst.open();
  }

}

// Perform the matrix-matrix product src * transpose(this) = dst.
void HypreMatrix::rightMultiplyTranspose( HypreMatrix const & src,
                                          HypreMatrix & dst,
                                          bool const closeResult ) const
{
  // Error check
  GEOSX_ASSERT_MSG( src.globalCols() == this->globalCols(),
                    "Incompatible matrix dimensions" );

  // Transpose this
  HYPRE_ParCSRMatrix tmp;
  hypre_ParCSRMatrixTranspose( src.m_parcsr_mat,
                               &tmp,
                               1 );

  // Compute product
  HYPRE_ParCSRMatrix dst_parcsr;
  dst_parcsr = hypre_ParMatmul( src.m_parcsr_mat,
                                tmp );

  // Create IJ layer (with matrix closed)
  dst.parCSRtoIJ( dst_parcsr );

  // Destroy temporary matrix
  hypre_ParCSRMatrixDestroy( tmp );

  // Reopen matrix if desired
  if( !closeResult )
  {
    dst.open();
  }

}


void HypreMatrix::parCSRtoIJ( HYPRE_ParCSRMatrix const & parCSRMatrix )
{
  this->reset();

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

  array1d< HYPRE_Int > info( 2 );
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

  this->close();

}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Compute residual.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Compute res = b - Ax (residual form).

void HypreMatrix::residual( HypreVector const & x,
                            HypreVector const & b,
                            HypreVector & r ) const
{

  HYPRE_ParVectorCopy( *b.getHypreParVectorPointer(),
                       *r.getHypreParVectorPointer() );
  hypre_ParCSRMatrixMatvec( -1.0,
                            m_parcsr_mat,
                            *x.getHypreParVectorPointer(),
                            1.0,
                            *r.getHypreParVectorPointer() );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Generalized matrix/vector product.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Compute gemv <tt>y = alpha*A*x + beta*y</tt>.

void HypreMatrix::gemv( real64 const alpha,
                        HypreVector const & x,
                        real64 const beta,
                        HypreVector & y,
                        bool useTranspose ) const
{

  if( !useTranspose )
  {
    hypre_ParCSRMatrixMatvec( static_cast< HYPRE_Real >( alpha ),
                              m_parcsr_mat,
                              *x.getHypreParVectorPointer(),
                              static_cast< HYPRE_Real >( beta ),
                              *y.getHypreParVectorPointer() );
  }
  else
  {
    hypre_ParCSRMatrixMatvecT( static_cast< HYPRE_Real >( alpha ),
                               m_parcsr_mat,
                               *x.getHypreParVectorPointer(),
                               static_cast< HYPRE_Real >( beta ),
                               *y.getHypreParVectorPointer() );
  }
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Scale.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Multiply all elements by scalingFactor.

void HypreMatrix::scale( real64 const scalingFactor )
{
  hypre_CSRMatrix * prt_diag_CSR = hypre_ParCSRMatrixDiag( m_parcsr_mat );
  HYPRE_Int diag_nnz = hypre_CSRMatrixNumNonzeros( prt_diag_CSR );
  HYPRE_Real * ptr_diag_data = hypre_CSRMatrixData( prt_diag_CSR );

  hypre_CSRMatrix * prt_offdiag_CSR = hypre_ParCSRMatrixOffd( m_parcsr_mat );
  HYPRE_Int offdiag_nnz = hypre_CSRMatrixNumNonzeros( prt_offdiag_CSR );
  HYPRE_Real * ptr_offdiag_data = hypre_CSRMatrixData( prt_offdiag_CSR );

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

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Left and right scaling
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""

void HypreMatrix::leftScale( HypreVector const & vec )
{
  hypre_Vector * ptr_vec = hypre_ParVectorLocalVector( *vec.getHypreParVectorPointer() );
  HYPRE_Real * ptr_vec_data = hypre_VectorData( ptr_vec );

  hypre_CSRMatrix * prt_diag_CSR = hypre_ParCSRMatrixDiag( m_parcsr_mat );
  HYPRE_Real * ptr_diag_data = hypre_CSRMatrixData( prt_diag_CSR );

  HYPRE_Int nrows = hypre_CSRMatrixNumRows( prt_diag_CSR );
  HYPRE_Int * IA = hypre_CSRMatrixI( prt_diag_CSR );

  for( HYPRE_Int i = 0 ; i < nrows ; ++i )
  {
    for( HYPRE_Int j = IA[i] ; j < IA[i + 1] ; ++j )
    {
      ptr_diag_data[j] *= ptr_vec_data[i];
    }
  }

  hypre_CSRMatrix * prt_offdiag_CSR = hypre_ParCSRMatrixOffd( m_parcsr_mat );
  HYPRE_Real * ptr_offdiag_data = hypre_CSRMatrixData( prt_offdiag_CSR );

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

  HYPRE_Int nrows = integer_conversion< HYPRE_Int >( this->localRows());

  array1d< HYPRE_BigInt > rows( nrows );
  array1d< HYPRE_Int > ncols( nrows );

  for( HYPRE_Int i = 0 ; i < nrows ; ++i )
  {
    rows[i] = this->ilower() + integer_conversion< HYPRE_BigInt >( i );
  }

  HYPRE_Int ierr;
  ierr = HYPRE_IJMatrixGetRowCounts( m_ij_mat,
                                     nrows,
                                     rows.data(),
                                     ncols.data() );
  GEOSX_UNUSED_VAR( ierr );
  GEOSX_ASSERT_MSG( ierr == 0,
                    "Error getting row counts - error code: " +
                    std::to_string( ierr ) );

  HYPRE_Int maxRowLength;
  HYPRE_Int localMaxRowLength = *std::max_element( ncols.begin(), ncols.end() );

  MpiWrapper::allReduce( &localMaxRowLength,
                         &maxRowLength,
                         1,
                         MPI_MAX,
                         hypre_IJMatrixComm( m_ij_mat ) );

  return integer_conversion< localIndex >( maxRowLength );
}

localIndex HypreMatrix::localRowLength( localIndex localRowIndex ) const
{
  return this->globalRowLength( this->getGlobalRowID( localRowIndex ) );
}

localIndex HypreMatrix::globalRowLength( globalIndex globalRowIndex ) const
{
  HYPRE_BigInt row = integer_conversion< HYPRE_BigInt >( globalRowIndex );
  HYPRE_Int ncols;

  HYPRE_Int ierr;
  ierr = HYPRE_IJMatrixGetRowCounts( m_ij_mat,
                                     1,
                                     &row,
                                     &ncols );
  GEOSX_UNUSED_VAR( ierr );
  GEOSX_ASSERT_MSG( ierr == 0,
                    "Error getting row counts - error code: " +
                    std::to_string( ierr ) );

  return integer_conversion< localIndex >( ncols );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get global row copy
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""

void HypreMatrix::getRowCopy( globalIndex globalRow,
                              array1d< globalIndex > & colIndices,
                              array1d< real64 > & values ) const
{
  GEOSX_ASSERT_MSG( m_is_pattern_fixed == true,
                    "matrix sparsity pattern not created" );
  GEOSX_ASSERT_MSG( this->getLocalRowID( globalRow ) >= 0,
                    std::string( "Row " ) +
                    std::to_string( globalRow ) +
                    std::string( " not on processor " ) +
                    std::to_string( MpiWrapper::Comm_rank( hypre_IJMatrixComm( m_ij_mat ) ) ) );

  HYPRE_Int n_entries;
  HYPRE_IJMatrixGetRowCounts( m_ij_mat,
                              1,
                              toHYPRE_BigInt( &globalRow ),
                              &n_entries );

  values.resize( n_entries );
  colIndices.resize( n_entries );
  hypre_IJMatrixGetValuesParCSR( m_ij_mat,
                                 -1,
                                 &n_entries,
                                 toHYPRE_BigInt( &globalRow ),
                                 toHYPRE_BigInt( colIndices.data() ),
                                 values.data() );
}

real64 HypreMatrix::getDiagValue( globalIndex globalRow ) const
{
  GEOSX_ASSERT_MSG( m_is_pattern_fixed == true,
                    "matrix sparsity pattern not created" );
  GEOSX_ASSERT_MSG( this->getLocalRowID( globalRow ) >= 0,
                    std::string( "Row " ) +
                    std::to_string( globalRow ) +
                    std::string( " not on processor " ) +
                    std::to_string( MpiWrapper::Comm_rank( hypre_IJMatrixComm( m_ij_mat ) ) ) );

  // Get local row index
  HYPRE_Int localRow = integer_conversion< HYPRE_Int >( globalRow - this->ilower());

  // Get diagonal block
  hypre_CSRMatrix * prt_CSR  = hypre_ParCSRMatrixDiag( m_parcsr_mat );
  HYPRE_Int *       IA       = hypre_CSRMatrixI( prt_CSR );
  HYPRE_Int *       JA       = hypre_CSRMatrixJ( prt_CSR );
  HYPRE_Real *          ptr_data = hypre_CSRMatrixData( prt_CSR );

  for( HYPRE_Int j = IA[localRow] ; j < IA[localRow + 1] ; ++j )
  {
    if( JA[j] == localRow )
    {
      return static_cast< real64 >( ptr_data[j] );
    }
  }
  return 0.0;
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Clear row
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Clear the row.  By default the diagonal value will be set
// to zero, but the user can pass a desired diagValue.

void HypreMatrix::clearRow( globalIndex const globalRow,
                            real64 const diagValue )
{
  GEOSX_ASSERT_MSG( m_is_pattern_fixed,
                    std::string( "HypreMatrix, sparsity pattern must be fixed" ) );
  GEOSX_ASSERT_MSG( !m_is_ready_to_use,
                    std::string( "HypreMatrix, matrix must be opened" ) );
  GEOSX_ASSERT_MSG( this->getLocalRowID( globalRow ) >= 0,
                    std::string( "HypreMatrix, it is not possible to clear" ) +
                    std::string( "rows on other processors" ) );

  // Get local row index
  HYPRE_Int localRow = integer_conversion< HYPRE_Int >( globalRow - this->ilower());

  // Clear row in diagonal block
  hypre_CSRMatrix * prt_CSR  = hypre_ParCSRMatrixDiag( m_parcsr_mat );
  HYPRE_Int *       IA       = hypre_CSRMatrixI( prt_CSR );
  HYPRE_Real *          ptr_data = hypre_CSRMatrixData( prt_CSR );

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
  this->set( globalRow, globalRow, diagValue );
}


// ---------------------------------------------------------
//  Accessors
// ---------------------------------------------------------

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get pointer.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Accessor for the pointer to the raw Epetra matrix

HYPRE_IJMatrix const * HypreMatrix::unwrappedPointer() const
{
  return &m_ij_mat;
}

HYPRE_IJMatrix * HypreMatrix::unwrappedPointer()
{
  return &m_ij_mat;
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// getLocalRowID
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Map a global row index to local row index
localIndex HypreMatrix::getLocalRowID( globalIndex const index ) const
{
  if( index < this->iupper() )
  {
    return std::max( integer_conversion< localIndex >( -1 ),
                     integer_conversion< localIndex >( index - this->ilower() ) );
  }
  else
  {
    return integer_conversion< localIndex >( -1 );
  }

}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// getGlobalRowID
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Map a local row index to global row index
globalIndex HypreMatrix::getGlobalRowID( localIndex const index ) const
{
  return this->ilower() + integer_conversion< globalIndex >( index );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get number of global rows.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Accessor for the number of global rows

globalIndex HypreMatrix::globalRows() const
{
  return hypre_IJMatrixGlobalNumRows( m_ij_mat );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get number of global columns.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Accessor for the number of global columns

globalIndex HypreMatrix::globalCols() const
{
  return hypre_IJMatrixGlobalNumCols( m_ij_mat );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get number of local rows.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""

localIndex HypreMatrix::localRows() const
{
  return integer_conversion< localIndex >( this->iupper() -
                                           this->ilower() );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get number of local columns.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""

localIndex HypreMatrix::localCols() const
{
  return integer_conversion< localIndex >( this->jupper() -
                                           this->jlower() );
}

//
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get the lower index owned by processor.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Accessor for the index of the first global row
globalIndex HypreMatrix::ilower() const
{
  HYPRE_BigInt ilower, iupper, jlower, jupper;

  HYPRE_IJMatrixGetLocalRange( m_ij_mat,
                               &ilower,
                               &iupper,
                               &jlower,
                               &jupper );
  return integer_conversion< globalIndex >( ilower );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get the upper index owned by processor.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Accessor for the index of the last global row
globalIndex HypreMatrix::iupper() const
{
  HYPRE_BigInt ilower, iupper, jlower, jupper;

  HYPRE_IJMatrixGetLocalRange( m_ij_mat,
                               &ilower,
                               &iupper,
                               &jlower,
                               &jupper );
  return integer_conversion< globalIndex >( iupper + 1 );
}

//
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get the lower index owned by processor.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Accessor for the index of the first global row
globalIndex HypreMatrix::jlower() const
{
  HYPRE_BigInt ilower, iupper, jlower, jupper;

  HYPRE_IJMatrixGetLocalRange( m_ij_mat,
                               &ilower,
                               &iupper,
                               &jlower,
                               &jupper );
  return integer_conversion< globalIndex >( jlower );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get the upper index owned by processor.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Accessor for the index of the last global row
globalIndex HypreMatrix::jupper() const
{
  HYPRE_BigInt ilower, iupper, jlower, jupper;

  HYPRE_IJMatrixGetLocalRange( m_ij_mat,
                               &ilower,
                               &iupper,
                               &jlower,
                               &jupper );
  return integer_conversion< globalIndex >( jupper + 1 );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get number of local nonzeros.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Accessor for the number of local nonzeros
localIndex HypreMatrix::localNonzeros() const
{
  hypre_CSRMatrix * prt_diag_CSR = hypre_ParCSRMatrixDiag( m_parcsr_mat );
  HYPRE_Int diag_nnz = hypre_CSRMatrixNumNonzeros( prt_diag_CSR );

  hypre_CSRMatrix * prt_offdiag_CSR = hypre_ParCSRMatrixOffd( m_parcsr_mat );
  HYPRE_Int offdiag_nnz = hypre_CSRMatrixNumNonzeros( prt_offdiag_CSR );

  return integer_conversion< localIndex >( diag_nnz + offdiag_nnz );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Print to terminal.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Wrapper to print the trilinos output of the matrix
void HypreMatrix::print( std::ostream & os ) const
{
  GEOSX_ERROR_IF( !m_is_ready_to_use,
                  "matrix appears to be empty (not created) or not finalized" );
  if( MpiWrapper::Comm_rank( hypre_IJMatrixComm( m_ij_mat ) ) == 0 )
  {
    os << "Hypre interface: no output on screen available\n";
    os << "                 use write method";
  }
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Write to matlab-compatible file
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Note: EpetraExt also supports a MatrixMarket format as well
void HypreMatrix::write( string const & filename,
                         bool const mtxFormat ) const
{
  GEOSX_ERROR_IF( !m_is_ready_to_use,
                  "matrix appears to be empty (not created) or not finalized" );
  if( mtxFormat )
  {
    GEOSX_WARNING( "MatrixMarket not available for HypreMatrix, default used" );
  }
  hypre_ParCSRMatrixPrintIJ( m_parcsr_mat,
                             1,
                             1,
                             filename.c_str() );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Inf-norm.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Returns the infinity norm of the matrix.
real64 HypreMatrix::normInf() const
{
  HypreMatrix matT;
  HYPRE_ParCSRMatrix parCSRT;
  hypre_ParCSRMatrixTranspose( m_parcsr_mat,
                               &parCSRT,
                               1 );
  matT.parCSRtoIJ( parCSRT );

  real64 normInf = matT.norm1();
  matT.reset();
  return normInf;
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// 1-norm.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Returns the one norm of the matrix.
real64 HypreMatrix::norm1() const
{
  hypre_CSRMatrix * prt_diag_CSR = hypre_ParCSRMatrixDiag( m_parcsr_mat );
  HYPRE_Real * ptr_diag_data = hypre_CSRMatrixData( prt_diag_CSR );

  HYPRE_Int nrows = hypre_CSRMatrixNumRows( prt_diag_CSR );
  HYPRE_Int * IA = hypre_CSRMatrixI( prt_diag_CSR );

  array1d< HYPRE_Real > row_abs_sum( nrows );
  row_abs_sum = 0.0;

  for( HYPRE_Int i = 0 ; i < nrows ; ++i )
  {
    for( HYPRE_Int j = IA[i] ; j < IA[i + 1] ; ++j )
    {
      row_abs_sum[i] += std::fabs( ptr_diag_data[j] );
    }
  }

  hypre_CSRMatrix * prt_offdiag_CSR = hypre_ParCSRMatrixOffd( m_parcsr_mat );
  HYPRE_Real * ptr_offdiag_data = hypre_CSRMatrixData( prt_offdiag_CSR );

  nrows = hypre_CSRMatrixNumRows( prt_offdiag_CSR );
  IA = hypre_CSRMatrixI( prt_offdiag_CSR );

  for( HYPRE_Int i = 0 ; i < nrows ; ++i )
  {
    for( HYPRE_Int j = IA[i] ; j < IA[i + 1] ; ++j )
    {
      row_abs_sum[i] += std::fabs( ptr_offdiag_data[j] );
    }
  }

  real64 norm1;
  real64 local_norm1 = static_cast< real64 >( *std::max_element( row_abs_sum.begin(), row_abs_sum.end() ) );
  MpiWrapper::allReduce( &local_norm1,
                         &norm1,
                         1,
                         MPI_MAX,
                         hypre_ParCSRMatrixComm( m_parcsr_mat ) );

  return norm1;

}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Frobenius-norm.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Returns the Frobenius norm of the matrix.
real64 HypreMatrix::normFrobenius() const
{
  hypre_CSRMatrix * prt_diag_CSR = hypre_ParCSRMatrixDiag( m_parcsr_mat );
  HYPRE_Int diag_nnz = hypre_CSRMatrixNumNonzeros( prt_diag_CSR );
  HYPRE_Real * ptr_diag_data = hypre_CSRMatrixData( prt_diag_CSR );

  hypre_CSRMatrix * prt_offdiag_CSR = hypre_ParCSRMatrixOffd( m_parcsr_mat );
  HYPRE_Int offdiag_nnz = hypre_CSRMatrixNumNonzeros( prt_offdiag_CSR );
  HYPRE_Real * ptr_offdiag_data = hypre_CSRMatrixData( prt_offdiag_CSR );

  HYPRE_Real normFrob = 0.0;
  HYPRE_Real local_normFrob = 0.0;

  for( HYPRE_Int i = 0 ; i < diag_nnz ; ++i )
    local_normFrob += ptr_diag_data[i] * ptr_diag_data[i];

  for( HYPRE_Int i = 0 ; i < offdiag_nnz ; ++i )
    local_normFrob += ptr_offdiag_data[i] * ptr_offdiag_data[i];

  MpiWrapper::allReduce( &local_normFrob,
                         &normFrob,
                         1,
                         MPI_SUM,
                         hypre_ParCSRMatrixComm( m_parcsr_mat ) );

  normFrob = sqrt( normFrob );
  return static_cast< real64 >( normFrob );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Is-assembled.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Boolean indicator. True = matrix assembled and ready to be used.
bool HypreMatrix::isAssembled() const
{
  return m_is_pattern_fixed;
}

// Boolean indicator. True = matrix assembled and ready to be used.
bool HypreMatrix::isClosed() const
{
  return m_is_ready_to_use;
}

std::ostream & operator<<( std::ostream & os,
                           HypreMatrix const & matrix )
{
  matrix.print( os );
  return os;
}

}// end geosx namespace

// END_RST_NARRATIVE
