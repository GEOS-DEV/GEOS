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

// Put everything under the geosx namespace.
namespace geosx
{

inline HYPRE_Int * toHYPRE_Int( globalIndex * const index )
{
  return reinterpret_cast<HYPRE_Int*>(index);
}

inline HYPRE_Int const * toHYPRE_Int( globalIndex const * const index )
{
  return reinterpret_cast<HYPRE_Int const*>(index);
}

// ----------------------------
// Constructors
// ----------------------------

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Create an empty matrix (meant to be used for declaration)
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""

HypreMatrix::HypreMatrix()
{
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Copy constructor
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""

HypreMatrix::HypreMatrix( HypreMatrix const &src )
{
  GEOS_ERROR_IF( src.unwrappedPointer() == nullptr,
                 "Input matrix looks empty" );
  GEOS_ERROR_IF( m_is_ready_to_use,
                 "Input matrix hasn't been closed before copy" );

  MPI_Comm comm = hypre_IJVectorComm( *src.unwrappedPointer() );
  HYPRE_Int ilower, iupper, jlower, jupper;
  HYPRE_Int objectType;

  HYPRE_IJMatrixGetLocalRange( *src.unwrappedPointer(),
                               &ilower,
                               &iupper,
                               &jlower,
                               &jupper );
  HYPRE_IJMatrixGetObjectType( *src.unwrappedPointer(),
                               &objectType );

  HYPRE_IJMatrixCreate( comm,
                        ilower,
                        iupper,
                        jlower,
                        jupper,
                        &m_ij_mat );
  HYPRE_IJMatrixSetObjectType( m_ij_mat, objectType );

  // Get number of non-zeroes per row
  HYPRE_Int nrows = iupper - ilower + 1;
  array1d<HYPRE_Int> rows( nrows );
  array1d<HYPRE_Int> ncols( nrows );

  for( HYPRE_Int i = ilower ; i <= iupper ; ++i )
    rows[i - ilower] = i;

  HYPRE_IJMatrixGetRowCounts( *src.unwrappedPointer(),
                              nrows,
                              rows.data(),
                              ncols.data() );

  HYPRE_Int size = 0;
  for( HYPRE_Int i = 0 ; i < nrows ; ++i )
    size += ncols[i];

  array1d<HYPRE_Int> cols( size );
  array1d<double> values( size );

  hypre_IJMatrixGetValuesParCSR( *src.unwrappedPointer(),
                                 -nrows,
                                 ncols.data(),
                                 rows.data(),
                                 cols.data(),
                                 values.data() );

  HYPRE_IJMatrixSetRowSizes( m_ij_mat,
                             ncols.data() );

  HYPRE_IJMatrixInitialize( m_ij_mat );

  HYPRE_IJMatrixSetValues( m_ij_mat,
                           nrows,
                           ncols.data(),
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
  if( m_ij_mat )
    HYPRE_IJMatrixDestroy( m_ij_mat );
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
  this->createWithGlobalSize( globalSize,
                              globalSize,
                              maxEntriesPerRow,
                              comm ); // just call general version

}

void HypreMatrix::createWithGlobalSize( globalIndex const globalRows,
                                        globalIndex const globalCols,
                                        localIndex const maxEntriesPerRow,
                                        MPI_Comm const & comm )
{
  int this_mpi_process;
  int n_mpi_process;
  MPI_Comm_rank( comm, &this_mpi_process );
  MPI_Comm_size( comm, &n_mpi_process );

  globalIndex localRowSize = globalRows / n_mpi_process;
  globalIndex rowResidual = globalRows % n_mpi_process;

  GEOS_ERROR_IF( localRowSize < 1,
                 "localRowSize is lower than 1: less than one processor per row" );

  globalIndex localColSize = globalCols / n_mpi_process;
  globalIndex colResidual = globalCols % n_mpi_process;

  GEOS_ERROR_IF( localColSize < 1,
                 "localColSize is lower than 1: less than one processor per column" );

  globalIndex ilower, iupper, jlower, jupper;

  if( this_mpi_process == 0 )
  {
    ilower = 0;
    localRowSize = localRowSize + rowResidual;
    iupper = localRowSize - 1;
    jlower = 0;
    localColSize = localColSize + colResidual;
    jupper = localColSize - 1;
  }
  else
  {
    ilower = this_mpi_process * localRowSize + rowResidual;
    iupper = ilower + localRowSize - 1;
    jlower = this_mpi_process * localColSize + colResidual;
    jupper = jlower + localColSize - 1;
  }

  HYPRE_IJMatrixCreate( comm,
                        integer_conversion<HYPRE_Int>( ilower ),
                        integer_conversion<HYPRE_Int>( iupper ),
                        integer_conversion<HYPRE_Int>( jlower ),
                        integer_conversion<HYPRE_Int>( jupper ),
                        &m_ij_mat );
  HYPRE_IJMatrixSetObjectType( m_ij_mat, HYPRE_PARCSR );

  array1d<HYPRE_Int> row_sizes( localRowSize );
  row_sizes = integer_conversion<HYPRE_Int>( maxEntriesPerRow );
  HYPRE_IJMatrixSetRowSizes( m_ij_mat,
                             row_sizes.data() );
  HYPRE_IJMatrixInitialize( m_ij_mat );
  m_is_ready_to_use = false;

}

void HypreMatrix::createWithLocalSize( localIndex const localSize,
                                       localIndex const maxEntriesPerRow,
                                       MPI_Comm const & comm )
{
  this->createWithLocalSize( localSize,
                             localSize,
                             maxEntriesPerRow,
                             comm ); // just call general version
}

void HypreMatrix::createWithLocalSize( localIndex const localRows,
                                       localIndex const localCols,
                                       localIndex const maxEntriesPerRow,
                                       MPI_Comm const & comm )
{
  GEOS_ERROR_IF( localRows < 1,
                 "local rows are lower than 1: less than one processor per row" );
  GEOS_ERROR_IF( localCols < 1,
                 "local columns are lower than 1: less than one processor per column" );

  int this_mpi_process;
  int n_mpi_process;
  MPI_Comm_rank( comm, &this_mpi_process );
  MPI_Comm_size( comm, &n_mpi_process );

  array1d<int> localSizeArray( n_mpi_process * 2 );
  array1d<int> tmp( 2 );
  tmp[0] = integer_conversion<int>( localRows );
  tmp[1] = integer_conversion<int>( localCols );

  MPI_Allgather( tmp.data(),
                 2,
                 MPI_INT,
                 localSizeArray.data(),
                 2,
                 MPI_INT,
                 comm );

  globalIndex ilower, iupper, jlower, jupper;

  ilower = 0;
  jlower = 0;
  for( int i = 0 ; i < this_mpi_process * 2 ; i += 2 )
  {
    ilower += integer_conversion<long long int>( localSizeArray[i] );
    jlower += integer_conversion<long long int>( localSizeArray[i + 1] );
  }
  iupper = ilower + integer_conversion<long long int>( localRows ) - 1;
  jupper = jlower + integer_conversion<long long int>( localCols ) - 1;

  HYPRE_IJMatrixCreate( comm,
                        integer_conversion<HYPRE_Int>( ilower ),
                        integer_conversion<HYPRE_Int>( iupper ),
                        integer_conversion<HYPRE_Int>( jlower ),
                        integer_conversion<HYPRE_Int>( jupper ),
                        &m_ij_mat );
  HYPRE_IJMatrixSetObjectType( m_ij_mat, HYPRE_PARCSR );

  array1d<HYPRE_Int> row_sizes( localRows );
  row_sizes = integer_conversion<HYPRE_Int>( maxEntriesPerRow );
  HYPRE_IJMatrixSetRowSizes( m_ij_mat,
                             row_sizes.data() );

  HYPRE_IJMatrixInitialize( m_ij_mat );
  m_is_ready_to_use = false;

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
  HYPRE_IJMatrixSetConstantValues( m_ij_mat, value );
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
  HYPRE_IJMatrixSetConstantValues( m_ij_mat, 0.0 );
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
      HYPRE_IJMatrixGetObject( m_ij_mat, (void**) &m_parcsr_mat );
      m_is_pattern_fixed = true;
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
  GEOS_ERROR_IF( m_ij_mat == nullptr,
	             "matrix appears to be empty (not created)" );
  HYPRE_Int ncols = 1;
  HYPRE_IJMatrixAddToValues( m_ij_mat,
                             1,
                             &ncols,
                             &rowIndex,
                             &colIndex,
                             &value );
}

void HypreMatrix::set( globalIndex const rowIndex,
                       globalIndex const colIndex,
                       real64 const value )
{
  GEOS_ERROR_IF( m_ij_mat == nullptr,
	             "matrix appears to be empty (not created)" );
  GEOS_ASSERT_MSG( this->ilower() <= rowIndex &&
		           rowIndex < this->iupper(),
	               "HypreMatrix, it is not possible to set values on other processors");


  HYPRE_Int ncols = 1;
  HYPRE_IJMatrixSetValues( m_ij_mat,
                           1,
                           &ncols,
                           &rowIndex,
                           &colIndex,
                           &value );

}

void HypreMatrix::insert( globalIndex const rowIndex,
                          globalIndex const colIndex,
                          real64 const value )
{
  this->set( rowIndex,
             colIndex,
             value );
}

// 1xN c-style

void HypreMatrix::add( globalIndex const rowIndex,
                       globalIndex const * colIndices,
                       real64 const * values,
                       localIndex size )
{
  GEOS_ERROR_IF( m_ij_mat == nullptr,
				 "matrix appears to be empty (not created)" );

  HYPRE_Int ncols = integer_conversion<HYPRE_Int>( size );
  HYPRE_IJMatrixAddToValues( m_ij_mat,
                             1,
                             &ncols,
                             &rowIndex,
                             colIndices,
                             values );
}

void HypreMatrix::set( globalIndex const rowIndex,
                       globalIndex const * colIndices,
                       real64 const * values,
                       localIndex size )
{
  GEOS_ERROR_IF( m_ij_mat == nullptr,
				 "matrix appears to be empty (not created)" );
  GEOS_ASSERT_MSG( this->ilower() <= rowIndex &&
		           rowIndex < this->iupper(),
	               "HypreMatrix, it is not possible to set values on other processors");

  HYPRE_Int ncols = integer_conversion<HYPRE_Int>( size );
  HYPRE_IJMatrixSetValues( m_ij_mat,
                           1,
                           &ncols,
                           &rowIndex,
                           colIndices,
                           values );
}

void HypreMatrix::insert( globalIndex const rowIndex,
                          globalIndex const * colIndices,
                          real64 const * values,
                          localIndex size )
{
  this->set( rowIndex,
             colIndices,
             values,
             size );
}

// 1xN array1d style

void HypreMatrix::add( globalIndex const rowIndex,
                       array1d<globalIndex> const &colIndices,
                       array1d<real64> const &values )
{
  GEOS_ERROR_IF( m_ij_mat == nullptr,
				 "matrix appears to be empty (not created)" );
  HYPRE_Int ncols = integer_conversion<HYPRE_Int>( colIndices.size() );
  HYPRE_IJMatrixAddToValues( m_ij_mat,
                             1,
                             &ncols,
                             &rowIndex,
                             colIndices.data(),
                             values.data() );
}

void HypreMatrix::set( globalIndex const rowIndex,
                       array1d<globalIndex> const &colIndices,
                       array1d<real64> const &values )
{
  GEOS_ERROR_IF( m_ij_mat == nullptr,
				 "matrix appears to be empty (not created)" );
  GEOS_ASSERT_MSG( this->ilower() <= rowIndex &&
		           rowIndex < this->iupper(),
	               "HypreMatrix, it is not possible to set values on other processors");


  HYPRE_Int ncols = integer_conversion<HYPRE_Int>( colIndices.size() );
  HYPRE_IJMatrixSetValues( m_ij_mat,
                           1,
                           &ncols,
                           &rowIndex,
                           colIndices.data(),
                           values.data() );
}

void HypreMatrix::insert( globalIndex const rowIndex,
                          array1d<globalIndex> const &colIndices,
                          array1d<real64> const &values )
{
  this->set( rowIndex,
             colIndices,
             values );
}

//// MxN array2d style

void HypreMatrix::add( array1d<globalIndex> const & rowIndices,
                       array1d<globalIndex> const & colIndices,
                       array2d<real64> const & values )
{
  GEOS_ERROR_IF( m_ij_mat == nullptr,
	             "vector appears to be empty (not created)" );
  globalIndex nCols = colIndices.size();
  for( globalIndex i = 0 ; i < rowIndices.size() ; ++i )
  {
    this->add( rowIndices[i],
               colIndices.data(),
               values.begin() + nCols * i,
               nCols );
  }
}

void HypreMatrix::set( array1d<globalIndex> const & rowIndices,
                       array1d<globalIndex> const & colIndices,
                       array2d<real64> const & values )
{
  GEOS_ERROR_IF( m_ij_mat == nullptr,
	             "vector appears to be empty (not created)" );
  GEOS_ASSERT_MSG( this->ilower() <= *std::min_element(rowIndices.data(),
	                                                   rowIndices.data() + rowIndices.size() ) &&
	                   *std::max_element(rowIndices.data(),
	                		             rowIndices.data() + rowIndices.size()) < this->iupper(),
	                   "HypreVector, it is not possible to set values on other processors");


  globalIndex nCols = colIndices.size();
  for( globalIndex i = 0 ; i < rowIndices.size() ; ++i )
  {
    this->set( rowIndices[i],
               colIndices.data(),
               values.begin() + nCols * i,
               nCols );
  }
}

void HypreMatrix::insert( array1d<globalIndex> const & rowIndices,
                          array1d<globalIndex> const & colIndices,
                          array2d<real64> const & values )
{
  this->set( rowIndices,
             colIndices,
             values );
}

// -------------------------
// Linear Algebra
// -------------------------

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Matrix/vector multiplication
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Perform the matrix-vector product A*src = dst.

void HypreMatrix::multiply( HypreVector const &src,
                            HypreVector &dst ) const
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
  GEOS_ASSERT_MSG( this->globalCols() == src.globalRows(),
		           "Incompatible matrix dimensions");
  hypre_ParCSRMatrix *dst_parcsr;
  dst_parcsr = hypre_ParMatmul( m_parcsr_mat,
		                        src.m_parcsr_mat );

  // PARCSR to IJ
  HYPRE_Int ierr = 0;
  HYPRE_Int first_local_row, last_local_row, local_num_rows;
  HYPRE_Int first_local_col, last_local_col;//, local_num_cols;
  HYPRE_Int M, N;
  HYPRE_Int i,j;
  HYPRE_Int                 local_row;
  HYPRE_Int                *diag_sizes;
  HYPRE_Int                *offdiag_sizes;
  HYPRE_Int  size;
  HYPRE_Int                *col_inds;
  HYPRE_Real *values;

  MPI_Comm comm = hypre_IJMatrixComm(m_ij_mat);

  ierr = HYPRE_ParCSRMatrixGetLocalRange( dst_parcsr,
            &first_local_row, &last_local_row ,
            &first_local_col, &last_local_col );

  local_num_rows = last_local_row - first_local_row + 1;

  ierr += HYPRE_ParCSRMatrixGetDims( dst_parcsr, &M, &N );

  ierr += HYPRE_IJMatrixCreate( comm, first_local_row, last_local_row,
                                first_local_col, last_local_col,
                                &dst.m_ij_mat );

  ierr += HYPRE_IJMatrixSetObjectType( 	dst.m_ij_mat, HYPRE_PARCSR );

  diag_sizes = hypre_CTAlloc(HYPRE_Int,  local_num_rows, HYPRE_MEMORY_HOST);
  offdiag_sizes = hypre_CTAlloc(HYPRE_Int,  local_num_rows, HYPRE_MEMORY_HOST);
  local_row = 0;
  for (i=first_local_row; i<= last_local_row; i++)
  {
    ierr += HYPRE_ParCSRMatrixGetRow( dst_parcsr, i, &size,
                                      &col_inds, &values );
    for (j=0; j < size; j++)
    {
      if (col_inds[j] < first_local_row || col_inds[j] > last_local_row)
      {
        offdiag_sizes[local_row]++;
      }
      else
      {
        diag_sizes[local_row]++;
      }
    }
    local_row++;
    ierr += HYPRE_ParCSRMatrixRestoreRow( dst_parcsr, i, &size,
                                          &col_inds, &values );
  }
  ierr += HYPRE_IJMatrixSetDiagOffdSizes( dst.m_ij_mat,
                                          (const HYPRE_Int *) diag_sizes,
                                          (const HYPRE_Int *) offdiag_sizes );
  hypre_TFree(diag_sizes, HYPRE_MEMORY_HOST);
  hypre_TFree(offdiag_sizes, HYPRE_MEMORY_HOST);

  ierr = HYPRE_IJMatrixInitialize( dst.m_ij_mat );

  for (i=first_local_row; i<= last_local_row; i++)
  {
    ierr += HYPRE_ParCSRMatrixGetRow( dst_parcsr, i, &size,
                                      &col_inds, &values );

           ierr += HYPRE_IJMatrixSetValues( dst.m_ij_mat, 1, &size, &i,
                                            (const HYPRE_Int *) col_inds,
                                            (const HYPRE_Real *) values );

           ierr += HYPRE_ParCSRMatrixRestoreRow( dst_parcsr, i, &size,
                                                 &col_inds, &values );
  }

  if (closeResult)
  {
	  dst.close();
	  ierr += HYPRE_ParCSRMatrixDestroy(dst_parcsr);


  }
  GEOS_ASSERT_MSG( ierr == 0,
				  "Error in driver building IJMatrix from parcsr matrix" );

}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Compute residual.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Compute res = b - Ax (residual form).

void HypreMatrix::residual( HypreVector const &x,
                            HypreVector const &b,
                            HypreVector &r ) const
                            {
//  m_matrix->Multiply( false, *x.unwrappedPointer(), *r.unwrappedPointer() );
//  r.axpby( 1.0, b, -1.0 );
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
                        HypreVector const &x,
                        real64 const beta,
                        HypreVector &y,
                        bool useTranspose )
{

  if( !useTranspose )
  {
    hypre_ParCSRMatrixMatvec( alpha,
                              m_parcsr_mat,
                              *x.getHypreParVectorPointer(),
                              beta,
                              *y.getHypreParVectorPointer() );
  }
  else
  {
    hypre_ParCSRMatrixMatvecT( alpha,
                               m_parcsr_mat,
                               *x.getHypreParVectorPointer(),
                               beta,
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
  double * ptr_diag_data = hypre_CSRMatrixData( prt_diag_CSR );

  hypre_CSRMatrix * prt_offdiag_CSR = hypre_ParCSRMatrixOffd( m_parcsr_mat );
  HYPRE_Int offdiag_nnz = hypre_CSRMatrixNumNonzeros( prt_offdiag_CSR );
  double * ptr_offdiag_data = hypre_CSRMatrixData( prt_offdiag_CSR );

  for( HYPRE_Int i = 0 ; i < diag_nnz ; ++i )
    ptr_diag_data[i] *= scalingFactor;

  for( HYPRE_Int i = 0 ; i < offdiag_nnz ; ++i )
    ptr_offdiag_data[i] *= scalingFactor;
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Left and right scaling
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""

void HypreMatrix::leftScale( HypreVector const &vec )
{
  hypre_Vector * ptr_vec = hypre_ParVectorLocalVector( *vec.getHypreParVectorPointer() );
  double * ptr_vec_data = hypre_VectorData( ptr_vec );

  hypre_CSRMatrix * prt_diag_CSR = hypre_ParCSRMatrixDiag( m_parcsr_mat );
  //HYPRE_Int diag_nnz = hypre_CSRMatrixNumNonzeros( prt_diag_CSR );
  double * ptr_diag_data = hypre_CSRMatrixData( prt_diag_CSR );

  HYPRE_Int nrows = hypre_CSRMatrixNumRows( prt_diag_CSR );
  HYPRE_Int * IA = hypre_CSRMatrixI( prt_diag_CSR );

  for( HYPRE_Int i = 0 ; i < nrows ; ++i )
    for( HYPRE_Int j = IA[i] ; j < IA[i + 1] ; ++j )
      ptr_diag_data[j] *= ptr_vec_data[i];

  hypre_CSRMatrix * prt_offdiag_CSR = hypre_ParCSRMatrixOffd( m_parcsr_mat );
  //HYPRE_Int offdiag_nnz = hypre_CSRMatrixNumNonzeros( prt_offdiag_CSR );
  double * ptr_offdiag_data = hypre_CSRMatrixData( prt_offdiag_CSR );

  nrows = hypre_CSRMatrixNumRows( prt_offdiag_CSR );
  IA = hypre_CSRMatrixI( prt_offdiag_CSR );

  for( HYPRE_Int i = 0 ; i < nrows ; ++i )
    for( HYPRE_Int j = IA[i] ; j < IA[i + 1] ; ++j )
      ptr_offdiag_data[j] *= ptr_vec_data[i];

}

//void EpetraMatrix::rightScale( EpetraVector const &vec )
//{
//  m_matrix->RightScale( *(*vec.unwrappedPointer())(0) );
//}
//
//void EpetraMatrix::leftRightScale( EpetraVector const &vecLeft,
//                                   EpetraVector const &vecRight )
//{
//  leftScale(vecLeft);
//  rightScale(vecRight);
//}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get global row copy
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""

void HypreMatrix::getRowCopy( globalIndex globalRow,
                              array1d<globalIndex> & colIndices,
                              array1d<real64> & values) const
{
	GEOS_ASSERT_MSG(  m_is_pattern_fixed == true,
			          "matrix sparsity pattern not created" );
    GEOS_ASSERT_MSG( this->ilower() <= globalRow &&
                     globalRow < this->iupper(),
					 std::string("Row ") +
					 std::to_string(globalRow) +
                     std::string(" not on processor ") +
					 std::to_string( MpiWrapper::Comm_rank( hypre_IJMatrixComm(m_ij_mat) ) ) );

    HYPRE_Int n_entries;
    HYPRE_IJMatrixGetRowCounts( m_ij_mat,
			                    1,
								toHYPRE_Int(&globalRow),
			                    &n_entries );

//	localIndex length = integer_conversion<localIndex>(n_entries);

	values.resize(n_entries);
	colIndices.resize(n_entries);
	hypre_IJMatrixGetValuesParCSR( m_ij_mat,
	                              -1,
	    						  &n_entries,
		                          toHYPRE_Int(&globalRow),
		         				  toHYPRE_Int(colIndices.data()),
								  values.data() );
}
//
//
//// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
//// Clear row
//// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
//// Clear the row.  By default the diagonal value will be set
//// to zero, but the user can pass a desired diagValue.
//
//void EpetraMatrix::clearRow( globalIndex const globalRow,
//                             real64 const diagValue )
//{
//  double *values = nullptr;
//  int length;
//
//  int err = m_matrix->ExtractGlobalRowView( globalRow, length, values );
//  GEOS_ERROR_IF(err != 0, "getRowView failed. This often happens if the requested global row is not local to this processor, or if close() hasn't been called.");
//
//  for(int j=0; j<length; ++j)
//    values[j] = 0.0;
//
//  set(globalRow,globalRow,diagValue);
//}
//
//
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
// Get number of global rows.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Accessor for the number of global rows

globalIndex HypreMatrix::globalRows() const
{
  return hypre_IJMatrixGlobalNumRows(m_ij_mat);
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get number of global columns.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Accessor for the number of global columns

globalIndex HypreMatrix::globalCols() const
{
  return hypre_IJMatrixGlobalNumCols(m_ij_mat);
}

//
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get the lower index owned by processor.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Accessor for the index of the first global row
globalIndex HypreMatrix::ilower() const
{
  HYPRE_Int ilower, iupper, jlower, jupper;

  HYPRE_IJMatrixGetLocalRange( m_ij_mat,
                               &ilower,
                               &iupper,
                               &jlower,
                               &jupper );
  return integer_conversion<globalIndex>( ilower );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get the upper index owned by processor.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Accessor for the index of the last global row
globalIndex HypreMatrix::iupper() const
{
  HYPRE_Int ilower, iupper, jlower, jupper;

  HYPRE_IJMatrixGetLocalRange( m_ij_mat,
                               &ilower,
                               &iupper,
                               &jlower,
                               &jupper );
  return integer_conversion<globalIndex>( iupper + 1);
}

//// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
//// Print to terminal.
//// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
//// Wrapper to print the trilinos output of the matrix
//void EpetraMatrix::print() const
//{
//  if( m_matrix.get() != nullptr )
//  {
//    std::cout << *m_matrix << std::endl;
//  }
//}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Write to matlab-compatible file
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Note: EpetraExt also supports a MatrixMarket format as well
void HypreMatrix::write( string const & filename ) const
                         {
  GEOS_ERROR_IF( !m_is_ready_to_use,
                 "matrix appears to be empty (not created) or not finalized" );
  HYPRE_IJMatrixPrint( m_ij_mat, filename.c_str() );
}

//// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
//// Inf-norm.
//// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
//// Returns the infinity norm of the matrix.
//real64 EpetraMatrix::normInf() const
//{
//  return m_matrix->NormInf();
//}
//
//// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
//// 1-norm.
//// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
//// Returns the one norm of the matrix.
//real64 EpetraMatrix::norm1() const
//{
//  return m_matrix->NormOne();
//}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Frobenius-norm.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Returns the Frobenius norm of the matrix.
real64 HypreMatrix::normFrobenius() const
{
  MPI_Comm comm = hypre_IJMatrixComm( m_ij_mat );

  hypre_CSRMatrix * prt_diag_CSR = hypre_ParCSRMatrixDiag( m_parcsr_mat );
  HYPRE_Int diag_nnz = hypre_CSRMatrixNumNonzeros( prt_diag_CSR );
  double * ptr_diag_data = hypre_CSRMatrixData( prt_diag_CSR );

  hypre_CSRMatrix * prt_offdiag_CSR = hypre_ParCSRMatrixOffd( m_parcsr_mat );
  HYPRE_Int offdiag_nnz = hypre_CSRMatrixNumNonzeros( prt_offdiag_CSR );
  double * ptr_offdiag_data = hypre_CSRMatrixData( prt_offdiag_CSR );

  double normFrob = 0.0;
  double local_normFrob = 0.0;

  for( HYPRE_Int i = 0 ; i < diag_nnz ; ++i )
    local_normFrob += ptr_diag_data[i] * ptr_diag_data[i];

  for( HYPRE_Int i = 0 ; i < offdiag_nnz ; ++i )
    local_normFrob += ptr_offdiag_data[i] * ptr_offdiag_data[i];

  MPI_Allreduce( &local_normFrob,
                 &normFrob,
                 1,
                 MPI_DOUBLE,
                 MPI_SUM,
                 comm );

  normFrob = sqrt( normFrob );
  return static_cast<real64>( normFrob );
}

//// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
//// Is-assembled.
//// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
//// Boolean indicator. True = matrix assembled and ready to be used.
//bool EpetraMatrix::isAssembled() const
//{
//  return assembled;
//}

}// end geosx namespace

// END_RST_NARRATIVE
