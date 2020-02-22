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
 * @file PetscSparseMatrix.cpp
 */

#include "PetscSparseMatrix.hpp"

#if !defined(PETSC_USE_64BIT_INDICES)
#define PETSC_USE_64BIT_INDICES
#endif

#include "linearAlgebra/interfaces/PetscUtils.hpp"
#include "mpiCommunications/MpiWrapper.hpp"

#include <petscvec.h>
#include <petscmat.h>

namespace geosx
{

PetscSparseMatrix::PetscSparseMatrix()
: LinearOperator(),
  MatrixBase(),
  m_mat{}
{}

PetscSparseMatrix::PetscSparseMatrix( PetscSparseMatrix const & src )
: PetscSparseMatrix()
{
  GEOSX_LAI_MATRIX_STATUS( src.ready() );
  MatDuplicate( src.m_mat, MAT_COPY_VALUES, &m_mat );
  m_assembled = true;
}

PetscSparseMatrix::~PetscSparseMatrix()
{
  reset();
}

void PetscSparseMatrix::createWithLocalSize( localIndex const localRows,
                                             localIndex const localCols,
                                             localIndex const maxEntriesPerRow,
                                             MPI_Comm const & comm )
{
  GEOSX_LAI_MATRIX_STATUS( closed() );
  GEOSX_LAI_ASSERT_GE( localRows, 0 );
  GEOSX_LAI_ASSERT_GE( localCols, 0 );
  GEOSX_LAI_ASSERT_GE( maxEntriesPerRow, 0 );

  reset();

  // set up matrix
  MatCreate( comm, &m_mat );
  MatSetType( m_mat, MATMPIAIJ );
  MatSetSizes( m_mat, localRows, localCols, PETSC_DETERMINE, PETSC_DETERMINE );
  MatMPIAIJSetPreallocation( m_mat, maxEntriesPerRow, nullptr, maxEntriesPerRow, nullptr );
  MatSetUp( m_mat );
}

void PetscSparseMatrix::createWithGlobalSize( globalIndex const globalRows,
                                              globalIndex const globalCols,
                                              localIndex const maxEntriesPerRow,
                                              MPI_Comm const & comm )
{
  GEOSX_LAI_MATRIX_STATUS( closed() );
  GEOSX_LAI_ASSERT_GE( globalRows, 0 );
  GEOSX_LAI_ASSERT_GE( globalCols, 0 );
  GEOSX_LAI_ASSERT_GE( maxEntriesPerRow, 0 );

  reset();

  // set up matrix
  MatCreate( comm, &m_mat );
  MatSetType( m_mat, MATMPIAIJ );
  MatSetSizes( m_mat, PETSC_DECIDE, PETSC_DECIDE, globalRows, globalCols );
  MatMPIAIJSetPreallocation( m_mat, maxEntriesPerRow, nullptr, maxEntriesPerRow, nullptr );
  MatSetUp( m_mat );
}

bool PetscSparseMatrix::created() const
{
  return m_mat != nullptr;
}

void PetscSparseMatrix::reset()
{
  MatrixBase::reset();
  MatDestroy( &m_mat );
}

void PetscSparseMatrix::set( real64 const value )
{
  GEOSX_LAI_MATRIX_STATUS( ready() );

  PetscInt firstrow;
  PetscInt lastrow;
  MatGetOwnershipRange( m_mat, &firstrow, &lastrow );

  PetscInt numEntries;
  const PetscInt * inds;
  const PetscScalar * vals;

  PetscInt numEntries_;
  array1d<PetscScalar> vals_;
  array1d<PetscInt> inds_;

  PetscInt const maxNumRows = MpiWrapper::Max( lastrow - firstrow, getComm() );

  // loop over rows
  for( PetscInt row = firstrow; row < lastrow; row++)
  {
    // get entries in row
    MatGetRow( m_mat, row, &numEntries, &inds, &vals );
    numEntries_ = numEntries;
    inds_.resize( numEntries_ );
    for ( int i = 0; i < numEntries_; i++ ) 
    {
      inds_[i] = inds[i];
    }   
    MatRestoreRow( m_mat, row, &numEntries, &inds, &vals );

    // set entries to value
    open();
    if( numEntries_ > 0 )
    {
      vals_.resize( numEntries_ );
      for ( int i = 0; i < numEntries_; i++ ) 
      {
        vals_[i] = value;
      }
      MatSetValues( m_mat, 1, &row, numEntries_, inds_, vals_, INSERT_VALUES );
    }
    close();
  }

  // TODO: is there a way to set local rows without global ops?
  // ensure all ranks call close() the same number of times
  PetscInt const numExtra = maxNumRows - (lastrow - firstrow);
  for( PetscInt i = 0; i < numExtra; ++i )
  {
    open();
    close();
  }
}

void PetscSparseMatrix::zero()
{
  GEOSX_LAI_MATRIX_STATUS( ready() );
  MatZeroEntries( m_mat );
}

void PetscSparseMatrix::open()
{
  GEOSX_LAI_MATRIX_STATUS( created() && closed() );
  m_closed = false;
}

void PetscSparseMatrix::close()
{
  GEOSX_LAI_MATRIX_STATUS( !closed() );
  MatAssemblyBegin( m_mat, MAT_FINAL_ASSEMBLY );
  MatAssemblyEnd( m_mat, MAT_FINAL_ASSEMBLY );
  MatSetOption( m_mat, MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE );
  MatSetOption( m_mat, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE );
  m_assembled = true;
  m_closed = true;
}

void PetscSparseMatrix::add( globalIndex const rowIndex,
                             globalIndex const colIndex,
                             real64 const value )
{
  GEOSX_LAI_MATRIX_STATUS( modifiable() );
  MatSetValue( m_mat, rowIndex, colIndex, value, ADD_VALUES );
}

void PetscSparseMatrix::set( globalIndex const rowIndex,
                             globalIndex const colIndex,
                             real64 const value )
{
  GEOSX_LAI_MATRIX_STATUS( modifiable() );
  MatSetValue( m_mat, rowIndex, colIndex, value, INSERT_VALUES );
}

void PetscSparseMatrix::insert( globalIndex const rowIndex,
                                globalIndex const colIndex,
                                real64 const value )
{
  GEOSX_LAI_MATRIX_STATUS( insertable() );
  MatSetValue( m_mat, rowIndex, colIndex, value, INSERT_VALUES );
}

void PetscSparseMatrix::add( globalIndex const rowIndex,
                             globalIndex const * colIndices,
                             real64 const * values,
                             localIndex size )
{
  GEOSX_LAI_MATRIX_STATUS( modifiable() );
  PetscInt rows[1] = {rowIndex};
  MatSetValues( m_mat, 1, rows, size, toPetscInt( colIndices ), values, ADD_VALUES );
}

void PetscSparseMatrix::set( globalIndex const rowIndex,
                             globalIndex const * colIndices,
                             real64 const * values,
                             localIndex size )
{
  GEOSX_LAI_MATRIX_STATUS( modifiable() );
  PetscInt rows[1] = {rowIndex};
  MatSetValues( m_mat, 1, rows, size, toPetscInt( colIndices ), values, INSERT_VALUES );
}

void PetscSparseMatrix::insert( globalIndex const rowIndex,
                                globalIndex const * colIndices,
                                real64 const * values,
                                localIndex size )
{
  GEOSX_LAI_MATRIX_STATUS( insertable() );
  PetscInt rows[1] = {rowIndex};
  MatSetValues( m_mat, 1, rows, size, toPetscInt( colIndices ), values, INSERT_VALUES );
}

void PetscSparseMatrix::add( globalIndex const rowIndex,
                             arraySlice1d<globalIndex const> const & colIndices,
                             arraySlice1d<real64 const> const &values )
{
  GEOSX_LAI_MATRIX_STATUS( modifiable() );
  PetscInt rows[1] = {rowIndex};
  MatSetValues( m_mat, 1, rows, values.size(), toPetscInt( colIndices ), values.data(), ADD_VALUES );
}

void PetscSparseMatrix::set( globalIndex const rowIndex,
                             arraySlice1d<globalIndex const> const &colIndices,
                             arraySlice1d<real64 const> const &values )
{
  GEOSX_LAI_MATRIX_STATUS( modifiable() );
  PetscInt rows[1] = {rowIndex};
  MatSetValues( m_mat, 1, rows, values.size(), toPetscInt( colIndices ), values.data(), INSERT_VALUES );
}

void PetscSparseMatrix::insert( globalIndex const rowIndex,
                                arraySlice1d<globalIndex const> const &colIndices,
                                arraySlice1d<real64 const> const &values )
{
  GEOSX_LAI_MATRIX_STATUS( insertable() );
  PetscInt rows[1] = {rowIndex};
  MatSetValues( m_mat, 1, rows, values.size(), toPetscInt( colIndices ), values.data(), INSERT_VALUES );
}

void PetscSparseMatrix::add( arraySlice1d<globalIndex const> const & rowIndices,
                             arraySlice1d<globalIndex const> const & colIndices,
                             arraySlice2d<real64 const, 1> const & values )
{
  GEOSX_LAI_MATRIX_STATUS( modifiable() );
  MatSetValues( m_mat,
                rowIndices.size(),
                toPetscInt(rowIndices),
                colIndices.size(),
                toPetscInt(colIndices),
                values.data(),
                ADD_VALUES );
}

void PetscSparseMatrix::set( arraySlice1d<globalIndex const> const & rowIndices,
                             arraySlice1d<globalIndex const> const & colIndices,
                             arraySlice2d<real64 const, 1> const & values )
{
  GEOSX_LAI_MATRIX_STATUS( modifiable() );
  MatSetValues( m_mat,
                rowIndices.size(),
                toPetscInt(rowIndices),
                colIndices.size(),
                toPetscInt(colIndices),
                values.data(),
                INSERT_VALUES );
}

void PetscSparseMatrix::insert( arraySlice1d<globalIndex const> const & rowIndices,
                                arraySlice1d<globalIndex const> const & colIndices,
                                arraySlice2d<real64 const, 1> const & values )
{
  GEOSX_LAI_MATRIX_STATUS( insertable() );
  MatSetValues( m_mat,
                rowIndices.size(),
                toPetscInt(rowIndices),
                colIndices.size(),
                toPetscInt(colIndices),
                values.data(),
                INSERT_VALUES );
}

void PetscSparseMatrix::add( arraySlice1d<globalIndex const> const & rowIndices,
                             arraySlice1d<globalIndex const> const & colIndices,
                             arraySlice2d<real64 const, 0> const & values )
{
  GEOSX_LAI_MATRIX_STATUS( modifiable() );
  MatSetOption( m_mat, MAT_ROW_ORIENTED, PETSC_FALSE );
  MatSetValues( m_mat,
                rowIndices.size(),
                toPetscInt(rowIndices),
                colIndices.size(),
                toPetscInt(colIndices),
                values.data(),
                ADD_VALUES );
  MatSetOption( m_mat, MAT_ROW_ORIENTED, PETSC_TRUE );
}

void PetscSparseMatrix::set( arraySlice1d<globalIndex const> const & rowIndices,
                             arraySlice1d<globalIndex const> const & colIndices,
                             arraySlice2d<real64 const, 0> const & values )
{
  GEOSX_LAI_MATRIX_STATUS( modifiable() );
  MatSetOption( m_mat, MAT_ROW_ORIENTED, PETSC_FALSE );
  MatSetValues( m_mat,
                rowIndices.size(),
                toPetscInt(rowIndices),
                colIndices.size(),
                toPetscInt(colIndices),
                values.data(),
                INSERT_VALUES );
  MatSetOption( m_mat, MAT_ROW_ORIENTED, PETSC_TRUE );
}

void PetscSparseMatrix::insert( arraySlice1d<globalIndex const> const & rowIndices,
                                arraySlice1d<globalIndex const> const & colIndices,
                                arraySlice2d<real64 const, 0> const & values )
{
  GEOSX_LAI_MATRIX_STATUS( insertable() );
  MatSetOption( m_mat, MAT_ROW_ORIENTED, PETSC_FALSE );
  MatSetValues( m_mat,
                rowIndices.size(),
                toPetscInt(rowIndices),
                colIndices.size(),
                toPetscInt(colIndices),
                values.data(),
                INSERT_VALUES );
  MatSetOption( m_mat, MAT_ROW_ORIENTED, PETSC_TRUE );
}

void PetscSparseMatrix::add( globalIndex const * rowIndices,
                             globalIndex const * colIndices,
                             real64 const * values,
                             localIndex const numRows,
                             localIndex const numCols )
{
  GEOSX_LAI_MATRIX_STATUS( modifiable() );
  MatSetValues( m_mat,
                numRows,
                toPetscInt(rowIndices),
                numCols,
                toPetscInt(colIndices),
                values,
                ADD_VALUES );
}

void PetscSparseMatrix::set( globalIndex const * rowIndices,
                             globalIndex const * colIndices,
                             real64 const * values,
                             localIndex const numRows,
                             localIndex const numCols )
{
  GEOSX_LAI_MATRIX_STATUS( modifiable() );
  MatSetValues( m_mat,
                numRows,
                toPetscInt(rowIndices),
                numCols,
                toPetscInt(colIndices),
                values,
                INSERT_VALUES );
}

void PetscSparseMatrix::insert( globalIndex const * rowIndices,
                                globalIndex const * colIndices,
                                real64 const * values,
                                localIndex const numRows,
                                localIndex const numCols )
{
  GEOSX_LAI_MATRIX_STATUS( insertable() );
  MatSetValues( m_mat,
                numRows,
                toPetscInt(rowIndices),
                numCols,
                toPetscInt(colIndices),
                values,
                INSERT_VALUES );
}

void PetscSparseMatrix::multiply( PetscVector const & src,
                                  PetscVector & dst ) const
{
  GEOSX_LAI_MATRIX_STATUS( ready() );
  GEOSX_LAI_ASSERT( src.ready() );
  GEOSX_LAI_ASSERT( dst.ready() );
  MatMult( m_mat, src.unwrapped(), dst.unwrapped() );
}

void PetscSparseMatrix::multiply( PetscSparseMatrix const & src, 
                                  PetscSparseMatrix & dst,
                                  bool const closeResult ) const
{
  GEOSX_LAI_MATRIX_STATUS( ready() );
  GEOSX_LAI_MATRIX_STATUS( src.ready() );
  GEOSX_LAI_ASSERT_EQ( numGlobalCols(), src.numGlobalRows() );

  MatReuse const reuse = dst.created() ? MAT_REUSE_MATRIX : MAT_INITIAL_MATRIX;
  if( !dst.created() )
  {
    dst.createWithLocalSize( numLocalRows(), src.numLocalCols(), 1, getComm() );
  }

  MatMatMult( m_mat, src.unwrapped(), reuse, PETSC_DEFAULT, &dst.unwrapped() );
  dst.m_assembled = closeResult;
  dst.m_closed = closeResult;
}

void PetscSparseMatrix::leftMultiplyTranspose( PetscSparseMatrix const & src,
                                               PetscSparseMatrix & dst,
                                               bool const closeResult ) const
{
  GEOSX_LAI_MATRIX_STATUS( ready() );
  GEOSX_LAI_MATRIX_STATUS( src.ready() );
  GEOSX_LAI_ASSERT_EQ( numGlobalRows(), src.numGlobalRows() );

  MatReuse const reuse = dst.created() ? MAT_REUSE_MATRIX : MAT_INITIAL_MATRIX;
  if( !dst.created() )
  {
    dst.createWithLocalSize( numLocalCols(), src.numLocalCols(), 1, getComm() );
  }

  MatTransposeMatMult( m_mat, src.unwrapped(), reuse, PETSC_DEFAULT, &dst.unwrapped() );
  dst.m_assembled = closeResult;
  dst.m_closed = closeResult;
}

void PetscSparseMatrix::rightMultiplyTranspose( PetscSparseMatrix const & src,
                                                PetscSparseMatrix & dst,
                                                bool const closeResult ) const
{
  GEOSX_LAI_MATRIX_STATUS( ready() );
  GEOSX_LAI_MATRIX_STATUS( src.ready() );
  GEOSX_LAI_ASSERT_EQ( numGlobalCols(), src.numGlobalCols() );

  MatReuse const reuse = dst.created() ? MAT_REUSE_MATRIX : MAT_INITIAL_MATRIX;
  if( !dst.created() )
  {
    dst.createWithLocalSize( numLocalRows(), src.numLocalRows(), 1, getComm() );
  }

  MatMatTransposeMult( m_mat, src.unwrapped(), reuse, PETSC_DEFAULT, &dst.unwrapped() );
  dst.m_assembled = closeResult;
  dst.m_closed = closeResult;
}

void PetscSparseMatrix::gemv( real64 const alpha,
                              PetscVector const & x,
                              real64 const beta,
                              PetscVector & y,
                              bool useTranspose ) const
{
  GEOSX_LAI_MATRIX_STATUS( ready() );

  PetscVector x_( x );
  PetscVector b_( x );

  x_.scale( alpha ); // alpha*x_
  y.scale( beta ); // beta*y

  if ( useTranspose ) 
  {
    MatMultTranspose( m_mat, x_.unwrapped(), b_.unwrapped() );
  }
  else
  {
    MatMult( m_mat, x_.unwrapped(), b_.unwrapped() ); // alpha*A*x_ = b_
  }
  VecAXPY( y.unwrapped(), 1, b_.unwrapped() ); // alpha*A*x_ + beta*y = y
}

void PetscSparseMatrix::scale( real64 const scalingFactor )
{
  GEOSX_LAI_MATRIX_STATUS( ready() );
  MatScale( m_mat, scalingFactor );
}

void PetscSparseMatrix::leftScale( PetscVector const &vec )
{
  GEOSX_LAI_MATRIX_STATUS( ready() );
  MatDiagonalScale( m_mat, vec.unwrapped(), nullptr );
}

void PetscSparseMatrix::rightScale( PetscVector const &vec )
{
  GEOSX_LAI_MATRIX_STATUS( ready() );
  MatDiagonalScale( m_mat, nullptr, vec.unwrapped() );
}

void PetscSparseMatrix::leftRightScale( PetscVector const &vecLeft,
                                        PetscVector const &vecRight )
{
  GEOSX_LAI_MATRIX_STATUS( ready() );
  MatDiagonalScale( m_mat, vecLeft.unwrapped(), vecRight.unwrapped() );
}

void PetscSparseMatrix::clearRow( globalIndex const globalRow,
                                  real64 const diagValue )
{
  GEOSX_LAI_MATRIX_STATUS( modifiable() );
  GEOSX_LAI_ASSERT_GE( globalRow, ilower() );
  GEOSX_LAI_ASSERT_GT( iupper(), globalRow );

  PetscInt rows[1] = {globalRow};
  MatZeroRows( m_mat, 1, rows, diagValue, nullptr, nullptr );
}

localIndex PetscSparseMatrix::maxRowLength() const
{
  GEOSX_LAI_MATRIX_STATUS( assembled() );
  localIndex maxLocalLength = 0;
  for( localIndex i = ilower(); i < iupper(); ++i )
  {
    maxLocalLength = std::max( maxLocalLength, globalRowLength( i ) );
  }
  return MpiWrapper::Max( maxLocalLength, getComm() );
}

localIndex PetscSparseMatrix::localRowLength( localIndex localRowIndex ) const
{
  return globalRowLength( getGlobalRowID( localRowIndex ) );
}

localIndex PetscSparseMatrix::globalRowLength( globalIndex globalRowIndex ) const
{
  GEOSX_LAI_MATRIX_STATUS( assembled() );
  PetscInt ncols;
  MatGetRow( m_mat, globalRowIndex, &ncols, nullptr, nullptr );
  localIndex const nnz = ncols;
  MatRestoreRow( m_mat, globalRowIndex, &ncols, nullptr, nullptr );
  return nnz;
}

void PetscSparseMatrix::getRowCopy( globalIndex globalRow,
                                    array1d<globalIndex> & colIndices,
                                    array1d<real64> & values ) const
{
  GEOSX_LAI_MATRIX_STATUS( ready() );
  GEOSX_LAI_ASSERT_GE( globalRow, ilower() );
  GEOSX_LAI_ASSERT_GT( iupper(), globalRow );

  PetscScalar const * vals;
  PetscInt const * inds;
  PetscInt numEntries;

  MatGetRow( m_mat, globalRow, &numEntries, &inds, &vals );

  values.resize( numEntries );
  colIndices.resize( numEntries );

  for ( int i = 0; i < numEntries; i++ ) 
  {
    colIndices[i] = inds[i];
  }

  for ( int i = 0; i < numEntries; i++ ) 
  {
    values[i] = vals[i];
  }

  MatRestoreRow( m_mat, globalRow, &numEntries, &inds, &vals );
}

real64 PetscSparseMatrix::getDiagValue( globalIndex globalRow ) const
{
  GEOSX_LAI_MATRIX_STATUS( assembled() );
  GEOSX_LAI_ASSERT_GE( globalRow, ilower());
  GEOSX_LAI_ASSERT_GT( iupper(), globalRow );

  PetscScalar const * vals = nullptr;
  PetscInt const * cols = nullptr;
  PetscInt ncols;

  MatGetRow( m_mat, globalRow, &ncols, &cols, &vals );
  for( PetscInt i = 0; i < ncols; i++ )
  {
    if( cols[i] == globalRow ) 
    {
      return vals[i];
    }
  }
  MatRestoreRow( m_mat, globalRow, &ncols, &cols, &vals );

  return 0.0;
}

Mat & PetscSparseMatrix::unwrapped()
{
  GEOSX_LAI_MATRIX_STATUS( created() );
  return m_mat;
}

const Mat & PetscSparseMatrix::unwrapped() const
{
  GEOSX_LAI_MATRIX_STATUS( created() );
  return m_mat;
}

globalIndex PetscSparseMatrix::numGlobalRows() const
{
  GEOSX_LAI_MATRIX_STATUS( created() );
  PetscInt num_rows;
  PetscInt num_cols;
  MatGetSize( m_mat, &num_rows, &num_cols );
  return integer_conversion<globalIndex>( num_rows );
}

globalIndex PetscSparseMatrix::numGlobalCols() const
{
  GEOSX_LAI_MATRIX_STATUS( created() );
  PetscInt num_rows;
  PetscInt num_cols;
  MatGetSize( m_mat, &num_rows, &num_cols );
  return integer_conversion<globalIndex>( num_cols );
}

globalIndex PetscSparseMatrix::ilower() const
{
  GEOSX_LAI_MATRIX_STATUS( created() );
  PetscInt firstrow;
  PetscInt lastrow;
  MatGetOwnershipRange( m_mat, &firstrow, &lastrow );
  return integer_conversion<globalIndex>( firstrow );
}

globalIndex PetscSparseMatrix::iupper() const
{
  GEOSX_LAI_MATRIX_STATUS( created() );
  PetscInt firstrow;
  PetscInt lastrow;
  MatGetOwnershipRange( m_mat, &firstrow, &lastrow );
  return integer_conversion<globalIndex>( lastrow );
}

localIndex PetscSparseMatrix::numLocalNonzeros() const
{
  GEOSX_LAI_MATRIX_STATUS( assembled() );
  PetscInt firstrow, lastrow;
  MatGetOwnershipRange( m_mat, &firstrow, &lastrow );

  PetscInt numEntries;
  localIndex result = 0;

  // loop over rows
  for( PetscInt row = firstrow; row < lastrow; ++row )
  {
    MatGetRow( m_mat, row, &numEntries, nullptr, nullptr );
    result += numEntries;
    MatRestoreRow( m_mat, row, &numEntries, nullptr, nullptr );
  }

  return result;
}

globalIndex PetscSparseMatrix::numGlobalNonzeros() const
{
  return MpiWrapper::Sum( integer_conversion<globalIndex>( numLocalNonzeros() ), getComm() );
}

real64 PetscSparseMatrix::normInf() const
{
  GEOSX_LAI_MATRIX_STATUS( ready() );
  real64 normInf;
  MatNorm( m_mat, NORM_INFINITY, &normInf );
  return normInf;
}

real64 PetscSparseMatrix::norm1() const
{
  GEOSX_LAI_MATRIX_STATUS( ready() );
  real64 norm1;
  MatNorm( m_mat, NORM_1, &norm1 );
  return norm1;
}

real64 PetscSparseMatrix::normFrobenius() const
{
  GEOSX_LAI_MATRIX_STATUS( ready() );
  real64 normFrob;
  MatNorm( m_mat, NORM_FROBENIUS, &normFrob );
  return normFrob;
}

localIndex PetscSparseMatrix::getLocalRowID( globalIndex const index ) const
{
  GEOSX_LAI_MATRIX_STATUS( created() );
  PetscInt low, high;
  MatGetOwnershipRange( m_mat, &low, &high);
  return (index >= low && index < high) ? integer_conversion< localIndex >( index - low ) : -1;
}

globalIndex PetscSparseMatrix::getGlobalRowID( localIndex const index ) const
{
  GEOSX_LAI_MATRIX_STATUS( created() );
  GEOSX_LAI_ASSERT_GE( index, 0 );
  GEOSX_LAI_ASSERT_GT( numLocalRows(), index );
  PetscInt low, high;
  MatGetOwnershipRange( m_mat, &low, &high );
  return integer_conversion<globalIndex>( index + low );
}

localIndex PetscSparseMatrix::numLocalCols() const
{
  GEOSX_LAI_MATRIX_STATUS( created() );
  PetscInt cols;
  MatGetSize( m_mat, nullptr, &cols );
  return integer_conversion<localIndex>( cols );
}

localIndex PetscSparseMatrix::numLocalRows() const
{
  GEOSX_LAI_MATRIX_STATUS( created() );
  PetscInt low, high;
  MatGetOwnershipRange( m_mat, &low, &high);
  return integer_conversion<localIndex >( high - low );
}

MPI_Comm PetscSparseMatrix::getComm() const
{
  GEOSX_LAI_MATRIX_STATUS( created() );
  MPI_Comm comm;
  PetscObjectGetComm( reinterpret_cast<PetscObject>( m_mat ), &comm );
  return comm;
}

void PetscSparseMatrix::print( std::ostream & os ) const
{
  GEOSX_LAI_MATRIX_STATUS( ready() );
  GEOSX_ERROR_IF( &os != &std::cout, "Only output to stdout currently supported" );
  MatView( m_mat, PETSC_VIEWER_STDOUT_( getComm() ) );
}

void PetscSparseMatrix::write( string const & filename,
                               LAIOutputFormat const format ) const
{
  GEOSX_LAI_MATRIX_STATUS( ready() );

  PetscViewer viewer;
  PetscViewerASCIIOpen( getComm(), filename.c_str(), &viewer );
  PetscViewerFormat petscFormat;

  switch( format )
  {
    case LAIOutputFormat::NATIVE_ASCII:
      petscFormat = PETSC_VIEWER_DEFAULT;
      break;
    case LAIOutputFormat::NATIVE_BINARY:
      petscFormat = PETSC_VIEWER_NATIVE;
      break;
    case LAIOutputFormat::MATLAB_ASCII:
      petscFormat = PETSC_VIEWER_ASCII_MATLAB;
      break;
    case LAIOutputFormat::MATLAB_BINARY:
      petscFormat = PETSC_VIEWER_BINARY_MATLAB;
      break;
    case LAIOutputFormat::MATRIX_MARKET:
      petscFormat = PETSC_VIEWER_ASCII_MATRIXMARKET;
      break;
    default:
      GEOSX_ERROR( "Unsupported matrix output format" );
  }

  PetscViewerPushFormat( viewer, petscFormat );
  MatView( m_mat, viewer );
  PetscViewerDestroy( &viewer );
}

} // end geosx namespace
