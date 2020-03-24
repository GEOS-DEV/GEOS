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
 * @file PetscMatrix.cpp
 */

#include "PetscMatrix.hpp"

#if !defined(PETSC_USE_64BIT_INDICES)
#define PETSC_USE_64BIT_INDICES
#endif

#include "PetscUtils.hpp"
#include "mpiCommunications/MpiWrapper.hpp"

#include <petscvec.h>
#include <petscmat.h>

namespace geosx
{

PetscMatrix::PetscMatrix()
  : LinearOperator(),
  MatrixBase(),
  m_mat{}
{}

PetscMatrix::PetscMatrix( PetscMatrix const & src )
  : PetscMatrix()
{
  GEOSX_LAI_ASSERT( src.ready() );
  GEOSX_LAI_CHECK_ERROR( MatDuplicate( src.m_mat, MAT_COPY_VALUES, &m_mat ) );
  m_assembled = true;
}

PetscMatrix::~PetscMatrix()
{
  reset();
}

void PetscMatrix::createWithLocalSize( localIndex const localRows,
                                       localIndex const localCols,
                                       localIndex const maxEntriesPerRow,
                                       MPI_Comm const & comm )
{
  GEOSX_LAI_ASSERT( closed() );
  GEOSX_LAI_ASSERT_GE( localRows, 0 );
  GEOSX_LAI_ASSERT_GE( localCols, 0 );
  GEOSX_LAI_ASSERT_GE( maxEntriesPerRow, 0 );

  reset();

  // set up matrix
  GEOSX_LAI_CHECK_ERROR( MatCreate( comm, &m_mat ) );
  GEOSX_LAI_CHECK_ERROR( MatSetType( m_mat, MATMPIAIJ ) );
  GEOSX_LAI_CHECK_ERROR( MatSetSizes( m_mat, localRows, localCols, PETSC_DETERMINE, PETSC_DETERMINE ) );
  GEOSX_LAI_CHECK_ERROR( MatMPIAIJSetPreallocation( m_mat, maxEntriesPerRow, nullptr, maxEntriesPerRow, nullptr ) );
  GEOSX_LAI_CHECK_ERROR( MatSetUp( m_mat ) );
  GEOSX_LAI_CHECK_ERROR( MatSetOption( m_mat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE ) );
}

void PetscMatrix::createWithGlobalSize( globalIndex const globalRows,
                                        globalIndex const globalCols,
                                        localIndex const maxEntriesPerRow,
                                        MPI_Comm const & comm )
{
  GEOSX_LAI_ASSERT( closed() );
  GEOSX_LAI_ASSERT_GE( globalRows, 0 );
  GEOSX_LAI_ASSERT_GE( globalCols, 0 );
  GEOSX_LAI_ASSERT_GE( maxEntriesPerRow, 0 );

  reset();

  // set up matrix
  GEOSX_LAI_CHECK_ERROR( MatCreate( comm, &m_mat ) );
  GEOSX_LAI_CHECK_ERROR( MatSetType( m_mat, MATMPIAIJ ) );
  GEOSX_LAI_CHECK_ERROR( MatSetSizes( m_mat, PETSC_DECIDE, PETSC_DECIDE, globalRows, globalCols ) );
  GEOSX_LAI_CHECK_ERROR( MatMPIAIJSetPreallocation( m_mat, maxEntriesPerRow, nullptr, maxEntriesPerRow, nullptr ) );
  GEOSX_LAI_CHECK_ERROR( MatSetUp( m_mat ) );
  GEOSX_LAI_CHECK_ERROR( MatSetOption( m_mat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE ) );
}

bool PetscMatrix::created() const
{
  return m_mat != nullptr;
}

void PetscMatrix::reset()
{
  MatrixBase::reset();
  GEOSX_LAI_CHECK_ERROR( MatDestroy( &m_mat ) );
}

void PetscMatrix::set( real64 const value )
{
  GEOSX_LAI_ASSERT( ready() );

  PetscInt firstrow;
  PetscInt lastrow;
  GEOSX_LAI_CHECK_ERROR( MatGetOwnershipRange( m_mat, &firstrow, &lastrow ) );

  PetscInt numEntries;
  const PetscInt * inds;
  const PetscScalar * vals;

  PetscInt numEntries_;
  array1d< PetscScalar > vals_;
  array1d< PetscInt > inds_;

  PetscInt const maxNumRows = MpiWrapper::Max( lastrow - firstrow, getComm() );

  // loop over rows
  for( PetscInt row = firstrow; row < lastrow; row++ )
  {
    // get entries in row
    GEOSX_LAI_CHECK_ERROR( MatGetRow( m_mat, row, &numEntries, &inds, &vals ) );
    numEntries_ = numEntries;
    inds_.resize( numEntries_ );
    for( int i = 0; i < numEntries_; i++ )
    {
      inds_[i] = inds[i];
    }
    GEOSX_LAI_CHECK_ERROR( MatRestoreRow( m_mat, row, &numEntries, &inds, &vals ) );

    // set entries to value
    if( numEntries_ > 0 )
    {
      vals_.resize( numEntries_ );
      for( int i = 0; i < numEntries_; i++ )
      {
        vals_[i] = value;
      }
      GEOSX_LAI_CHECK_ERROR( MatSetValues( m_mat, 1, &row, numEntries_, inds_, vals_, INSERT_VALUES ) );
    }
    GEOSX_LAI_CHECK_ERROR( MatAssemblyBegin( m_mat, MAT_FINAL_ASSEMBLY ) );
    GEOSX_LAI_CHECK_ERROR( MatAssemblyEnd( m_mat, MAT_FINAL_ASSEMBLY ) );
  }

  // TODO: is there a way to set local rows without global/collective calls?
  // ensure all ranks call MatAssemblyEnd the same number of times
  PetscInt const numExtra = maxNumRows - (lastrow - firstrow);
  for( PetscInt i = 0; i < numExtra; ++i )
  {
    GEOSX_LAI_CHECK_ERROR( MatAssemblyBegin( m_mat, MAT_FLUSH_ASSEMBLY ) );
    GEOSX_LAI_CHECK_ERROR( MatAssemblyEnd( m_mat, MAT_FLUSH_ASSEMBLY ) );
  }

  // call one final time
  GEOSX_LAI_CHECK_ERROR( MatAssemblyBegin( m_mat, MAT_FINAL_ASSEMBLY ) );
  GEOSX_LAI_CHECK_ERROR( MatAssemblyEnd( m_mat, MAT_FINAL_ASSEMBLY ) );
}

void PetscMatrix::zero()
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_CHECK_ERROR( MatZeroEntries( m_mat ) );
}

void PetscMatrix::open()
{
  GEOSX_LAI_ASSERT( created() && closed() );
  m_closed = false;
}

void PetscMatrix::close()
{
  GEOSX_LAI_ASSERT( !closed() );

  // Initiate global assembly
  GEOSX_LAI_CHECK_ERROR( MatAssemblyBegin( m_mat, MAT_FINAL_ASSEMBLY ) );
  GEOSX_LAI_CHECK_ERROR( MatAssemblyEnd( m_mat, MAT_FINAL_ASSEMBLY ) );

  // Clear any rows that have been marked for clearing
  if( assembled() )
  {
    GEOSX_LAI_CHECK_ERROR( MatZeroRows( m_mat, m_rowsToClear.size(), m_rowsToClear, 1.0, nullptr, nullptr ) );

    for( localIndex i = 0; i < m_rowsToClear.size(); ++i )
    {
      set( m_rowsToClear[i], m_rowsToClear[i], m_diagValues[i] );
    }
    GEOSX_LAI_CHECK_ERROR( MatAssemblyBegin( m_mat, MAT_FINAL_ASSEMBLY ) );
    GEOSX_LAI_CHECK_ERROR( MatAssemblyEnd( m_mat, MAT_FINAL_ASSEMBLY ) );
    m_rowsToClear.clear();
    m_diagValues.clear();
  }

  // If this is initial assembly, set some useful options on the matrix
  if( !assembled() )
  {
    GEOSX_LAI_CHECK_ERROR( MatSetOption( m_mat, MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE ) );
    GEOSX_LAI_CHECK_ERROR( MatSetOption( m_mat, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE ) );
    GEOSX_LAI_CHECK_ERROR( MatSetOption( m_mat, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE ) );
    GEOSX_LAI_CHECK_ERROR( MatSetOption( m_mat, MAT_NO_OFF_PROC_ZERO_ROWS, PETSC_TRUE ) );
  }
  m_assembled = true;
  m_closed = true;
}

void PetscMatrix::add( globalIndex const rowIndex,
                       globalIndex const colIndex,
                       real64 const value )
{
  GEOSX_LAI_ASSERT( modifiable() );
  GEOSX_LAI_CHECK_ERROR( MatSetValue( m_mat, rowIndex, colIndex, value, ADD_VALUES ) );
}

void PetscMatrix::set( globalIndex const rowIndex,
                       globalIndex const colIndex,
                       real64 const value )
{
  GEOSX_LAI_ASSERT( modifiable() );
  GEOSX_LAI_CHECK_ERROR( MatSetValue( m_mat, rowIndex, colIndex, value, INSERT_VALUES ) );
}

void PetscMatrix::insert( globalIndex const rowIndex,
                          globalIndex const colIndex,
                          real64 const value )
{
  GEOSX_LAI_ASSERT( insertable() );
  GEOSX_LAI_CHECK_ERROR( MatSetValue( m_mat, rowIndex, colIndex, value, INSERT_VALUES ) );
}

void PetscMatrix::add( globalIndex const rowIndex,
                       globalIndex const * colIndices,
                       real64 const * values,
                       localIndex size )
{
  GEOSX_LAI_ASSERT( modifiable() );
  PetscInt rows[1] = {rowIndex};
  GEOSX_LAI_CHECK_ERROR( MatSetValues( m_mat, 1, rows, size, toPetscInt( colIndices ), values, ADD_VALUES ) );
}

void PetscMatrix::set( globalIndex const rowIndex,
                       globalIndex const * colIndices,
                       real64 const * values,
                       localIndex size )
{
  GEOSX_LAI_ASSERT( modifiable() );
  PetscInt rows[1] = {rowIndex};
  GEOSX_LAI_CHECK_ERROR( MatSetValues( m_mat, 1, rows, size, toPetscInt( colIndices ), values, INSERT_VALUES ) );
}

void PetscMatrix::insert( globalIndex const rowIndex,
                          globalIndex const * colIndices,
                          real64 const * values,
                          localIndex size )
{
  GEOSX_LAI_ASSERT( insertable() );
  PetscInt rows[1] = {rowIndex};
  GEOSX_LAI_CHECK_ERROR( MatSetValues( m_mat, 1, rows, size, toPetscInt( colIndices ), values, INSERT_VALUES ) );
}

void PetscMatrix::add( globalIndex const rowIndex,
                       arraySlice1d< globalIndex const > const & colIndices,
                       arraySlice1d< real64 const > const & values )
{
  GEOSX_LAI_ASSERT( modifiable() );
  PetscInt rows[1] = {rowIndex};
  GEOSX_LAI_CHECK_ERROR( MatSetValues( m_mat,
                                       1,
                                       rows,
                                       colIndices.size(),
                                       toPetscInt( colIndices ),
                                       values,
                                       ADD_VALUES ) );
}

void PetscMatrix::set( globalIndex const rowIndex,
                       arraySlice1d< globalIndex const > const & colIndices,
                       arraySlice1d< real64 const > const & values )
{
  GEOSX_LAI_ASSERT( modifiable() );
  PetscInt rows[1] = {rowIndex};
  GEOSX_LAI_CHECK_ERROR( MatSetValues( m_mat,
                                       1,
                                       rows,
                                       colIndices.size(),
                                       toPetscInt( colIndices ),
                                       values,
                                       INSERT_VALUES ) );
}

void PetscMatrix::insert( globalIndex const rowIndex,
                          arraySlice1d< globalIndex const > const & colIndices,
                          arraySlice1d< real64 const > const & values )
{
  GEOSX_LAI_ASSERT( insertable() );
  PetscInt rows[1] = {rowIndex};
  GEOSX_LAI_CHECK_ERROR( MatSetValues( m_mat,
                                       1,
                                       rows,
                                       colIndices.size(),
                                       toPetscInt( colIndices ),
                                       values,
                                       INSERT_VALUES ) );
}

void PetscMatrix::add( arraySlice1d< globalIndex const > const & rowIndices,
                       arraySlice1d< globalIndex const > const & colIndices,
                       arraySlice2d< real64 const, 1 > const & values )
{
  GEOSX_LAI_ASSERT( modifiable() );
  GEOSX_LAI_CHECK_ERROR( MatSetValues( m_mat,
                                       rowIndices.size(),
                                       toPetscInt( rowIndices ),
                                       colIndices.size(),
                                       toPetscInt( colIndices ),
                                       values.dataIfContiguous(),
                                       ADD_VALUES ) );
}

void PetscMatrix::set( arraySlice1d< globalIndex const > const & rowIndices,
                       arraySlice1d< globalIndex const > const & colIndices,
                       arraySlice2d< real64 const, 1 > const & values )
{
  GEOSX_LAI_ASSERT( modifiable() );
  GEOSX_LAI_CHECK_ERROR( MatSetValues( m_mat,
                                       rowIndices.size(),
                                       toPetscInt( rowIndices ),
                                       colIndices.size(),
                                       toPetscInt( colIndices ),
                                       values.dataIfContiguous(),
                                       INSERT_VALUES ) );
}

void PetscMatrix::insert( arraySlice1d< globalIndex const > const & rowIndices,
                          arraySlice1d< globalIndex const > const & colIndices,
                          arraySlice2d< real64 const, 1 > const & values )
{
  GEOSX_LAI_ASSERT( insertable() );
  GEOSX_LAI_CHECK_ERROR( MatSetValues( m_mat,
                                       rowIndices.size(),
                                       toPetscInt( rowIndices ),
                                       colIndices.size(),
                                       toPetscInt( colIndices ),
                                       values.dataIfContiguous(),
                                       INSERT_VALUES ) );
}

void PetscMatrix::add( arraySlice1d< globalIndex const > const & rowIndices,
                       arraySlice1d< globalIndex const > const & colIndices,
                       arraySlice2d< real64 const, 0 > const & values )
{
  GEOSX_LAI_ASSERT( modifiable() );
  GEOSX_LAI_CHECK_ERROR( MatSetOption( m_mat, MAT_ROW_ORIENTED, PETSC_FALSE ) );
  GEOSX_LAI_CHECK_ERROR( MatSetValues( m_mat,
                                       rowIndices.size(),
                                       toPetscInt( rowIndices ),
                                       colIndices.size(),
                                       toPetscInt( colIndices ),
                                       values.dataIfContiguous(),
                                       ADD_VALUES ) );
  GEOSX_LAI_CHECK_ERROR( MatSetOption( m_mat, MAT_ROW_ORIENTED, PETSC_TRUE ) );
}

void PetscMatrix::set( arraySlice1d< globalIndex const > const & rowIndices,
                       arraySlice1d< globalIndex const > const & colIndices,
                       arraySlice2d< real64 const, 0 > const & values )
{
  GEOSX_LAI_ASSERT( modifiable() );
  GEOSX_LAI_CHECK_ERROR( MatSetOption( m_mat, MAT_ROW_ORIENTED, PETSC_FALSE ) );
  GEOSX_LAI_CHECK_ERROR( MatSetValues( m_mat,
                                       rowIndices.size(),
                                       toPetscInt( rowIndices ),
                                       colIndices.size(),
                                       toPetscInt( colIndices ),
                                       values.dataIfContiguous(),
                                       INSERT_VALUES ) );
  GEOSX_LAI_CHECK_ERROR( MatSetOption( m_mat, MAT_ROW_ORIENTED, PETSC_TRUE ) );
}

void PetscMatrix::insert( arraySlice1d< globalIndex const > const & rowIndices,
                          arraySlice1d< globalIndex const > const & colIndices,
                          arraySlice2d< real64 const, 0 > const & values )
{
  GEOSX_LAI_ASSERT( insertable() );
  GEOSX_LAI_CHECK_ERROR( MatSetOption( m_mat, MAT_ROW_ORIENTED, PETSC_FALSE ) );
  GEOSX_LAI_CHECK_ERROR( MatSetValues( m_mat,
                                       rowIndices.size(),
                                       toPetscInt( rowIndices ),
                                       colIndices.size(),
                                       toPetscInt( colIndices ),
                                       values.dataIfContiguous(),
                                       INSERT_VALUES ) );
  GEOSX_LAI_CHECK_ERROR( MatSetOption( m_mat, MAT_ROW_ORIENTED, PETSC_TRUE ) );
}

void PetscMatrix::add( globalIndex const * rowIndices,
                       globalIndex const * colIndices,
                       real64 const * values,
                       localIndex const numRows,
                       localIndex const numCols )
{
  GEOSX_LAI_ASSERT( modifiable() );
  GEOSX_LAI_CHECK_ERROR( MatSetValues( m_mat,
                                       numRows,
                                       toPetscInt( rowIndices ),
                                       numCols,
                                       toPetscInt( colIndices ),
                                       values,
                                       ADD_VALUES ) );
}

void PetscMatrix::set( globalIndex const * rowIndices,
                       globalIndex const * colIndices,
                       real64 const * values,
                       localIndex const numRows,
                       localIndex const numCols )
{
  GEOSX_LAI_ASSERT( modifiable() );
  GEOSX_LAI_CHECK_ERROR( MatSetValues( m_mat,
                                       numRows,
                                       toPetscInt( rowIndices ),
                                       numCols,
                                       toPetscInt( colIndices ),
                                       values,
                                       INSERT_VALUES ) );
}

void PetscMatrix::insert( globalIndex const * rowIndices,
                          globalIndex const * colIndices,
                          real64 const * values,
                          localIndex const numRows,
                          localIndex const numCols )
{
  GEOSX_LAI_ASSERT( insertable() );
  GEOSX_LAI_CHECK_ERROR( MatSetValues( m_mat,
                                       numRows,
                                       toPetscInt( rowIndices ),
                                       numCols,
                                       toPetscInt( colIndices ),
                                       values,
                                       INSERT_VALUES ) );
}

void PetscMatrix::apply( PetscVector const & src,
                         PetscVector & dst ) const
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT( src.ready() );
  GEOSX_LAI_ASSERT( dst.ready() );
  GEOSX_LAI_ASSERT_EQ( numGlobalRows(), dst.globalSize() );
  GEOSX_LAI_ASSERT_EQ( numGlobalCols(), src.globalSize() );

  GEOSX_LAI_CHECK_ERROR( MatMult( m_mat, src.unwrapped(), dst.unwrapped() ) );
}

void PetscMatrix::applyTranspose( Vector const & src,
                                  Vector & dst ) const
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT( src.ready() );
  GEOSX_LAI_ASSERT( dst.ready() );
  GEOSX_LAI_ASSERT_EQ( numGlobalCols(), dst.globalSize() );
  GEOSX_LAI_ASSERT_EQ( numGlobalRows(), src.globalSize() );

  GEOSX_LAI_CHECK_ERROR( MatMultTranspose( m_mat, src.unwrapped(), dst.unwrapped() ) );
}

void PetscMatrix::multiply( PetscMatrix const & src,
                            PetscMatrix & dst ) const
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT( src.ready() );
  GEOSX_LAI_ASSERT_EQ( numGlobalCols(), src.numGlobalRows() );

  dst.reset();
  GEOSX_LAI_CHECK_ERROR( MatMatMult( m_mat, src.m_mat, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &dst.m_mat ) );
  dst.m_assembled = true;
  dst.m_closed = true;
}

void PetscMatrix::leftMultiplyTranspose( PetscMatrix const & src,
                                         PetscMatrix & dst ) const
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT( src.ready() );
  GEOSX_LAI_ASSERT_EQ( numGlobalRows(), src.numGlobalRows() );

  dst.reset();
  GEOSX_LAI_CHECK_ERROR( MatTransposeMatMult( m_mat, src.m_mat, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &dst.m_mat ) );
  dst.m_assembled = true;
  dst.m_closed = true;
}

void PetscMatrix::rightMultiplyTranspose( PetscMatrix const & src,
                                          PetscMatrix & dst ) const
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT( src.ready() );
  GEOSX_LAI_ASSERT_EQ( numGlobalCols(), src.numGlobalCols() );

  dst.reset();
  GEOSX_LAI_CHECK_ERROR( MatMatTransposeMult( m_mat, src.m_mat, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &dst.m_mat ) );
  dst.m_assembled = true;
  dst.m_closed = true;
}

void PetscMatrix::gemv( real64 const alpha,
                        PetscVector const & x,
                        real64 const beta,
                        PetscVector & y,
                        bool useTranspose ) const
{
  GEOSX_LAI_ASSERT( ready() );

  PetscVector x_( x );
  PetscVector b_( x );

  x_.scale( alpha ); // alpha*x_
  y.scale( beta ); // beta*y

  if( useTranspose )
  {
    GEOSX_LAI_CHECK_ERROR( MatMultTranspose( m_mat, x_.unwrapped(), b_.unwrapped() ) );
  }
  else
  {
    GEOSX_LAI_CHECK_ERROR( MatMult( m_mat, x_.unwrapped(), b_.unwrapped() ) ); // alpha*A*x_ = b_
  }
  GEOSX_LAI_CHECK_ERROR( VecAXPY( y.unwrapped(), 1, b_.unwrapped() ) ); // alpha*A*x_ + beta*y = y
}

void PetscMatrix::scale( real64 const scalingFactor )
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_CHECK_ERROR( MatScale( m_mat, scalingFactor ) );
}

void PetscMatrix::leftScale( PetscVector const & vec )
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_CHECK_ERROR( MatDiagonalScale( m_mat, vec.unwrapped(), nullptr ) );
}

void PetscMatrix::rightScale( PetscVector const & vec )
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_CHECK_ERROR( MatDiagonalScale( m_mat, nullptr, vec.unwrapped() ) );
}

void PetscMatrix::leftRightScale( PetscVector const & vecLeft,
                                  PetscVector const & vecRight )
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_CHECK_ERROR( MatDiagonalScale( m_mat, vecLeft.unwrapped(), vecRight.unwrapped() ) );
}

void PetscMatrix::multiplyRAP( PetscMatrix const & R,
                               PetscMatrix const & P,
                               PetscMatrix & dst ) const
{
  // No builtin RAP function in PETSc, so use double product for now
  // TODO: research this better
  MatrixBase::multiplyRAP( R, P, dst );
}

void PetscMatrix::multiplyPtAP( PetscMatrix const & P,
                                PetscMatrix & dst ) const
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT( P.ready() );
  GEOSX_LAI_ASSERT_EQ( numGlobalRows(), P.numGlobalRows() );
  GEOSX_LAI_ASSERT_EQ( numGlobalCols(), P.numGlobalRows() );

  dst.reset();
  GEOSX_LAI_CHECK_ERROR( MatPtAP( m_mat, P.m_mat, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &dst.m_mat ) );
  dst.m_assembled = true;
}

void PetscMatrix::transpose( PetscMatrix & dst ) const
{
  GEOSX_LAI_ASSERT( ready() );

  dst.reset();
  GEOSX_LAI_CHECK_ERROR( MatTranspose( m_mat, MAT_INITIAL_MATRIX, &dst.m_mat ) );
  dst.m_assembled = true;
}

real64 PetscMatrix::clearRow( globalIndex const globalRow,
                              bool const keepDiag,
                              real64 const diagValue )
{
  GEOSX_LAI_ASSERT( modifiable() );
  GEOSX_LAI_ASSERT_GE( globalRow, ilower() );
  GEOSX_LAI_ASSERT_GT( iupper(), globalRow );

  PetscInt numEntries = 0;
  PetscInt const * inds;
  PetscReal const * vals;
  GEOSX_LAI_CHECK_ERROR( MatGetRow( m_mat, globalRow, &numEntries, &inds, &vals ) );

  bool const square = numGlobalRows() == numGlobalCols();

  real64 oldDiag = 0.0;
  if( square )
  {
    for( localIndex i = 0; i < numEntries; ++i )
    {
      if( inds[i] == globalRow )
      {
        oldDiag = vals[i];
      }
    }
  }
  GEOSX_LAI_CHECK_ERROR( MatRestoreRow( m_mat, globalRow, &numEntries, &inds, &vals ) );

  // There is no legitimate way in PETSc to zero out a row in a local fashion,
  // without any collective calls. Therefore, we store the rows to be cleared
  // and perform the clearing at the next close() call.
  m_rowsToClear.push_back( globalRow );
  m_diagValues.push_back( keepDiag ? oldDiag : diagValue );

  return oldDiag;
}

localIndex PetscMatrix::maxRowLength() const
{
  GEOSX_LAI_ASSERT( assembled() );
  localIndex maxLocalLength = 0;
  for( localIndex i = ilower(); i < iupper(); ++i )
  {
    maxLocalLength = std::max( maxLocalLength, globalRowLength( i ) );
  }
  return MpiWrapper::Max( maxLocalLength, getComm() );
}

localIndex PetscMatrix::localRowLength( localIndex localRowIndex ) const
{
  return globalRowLength( getGlobalRowID( localRowIndex ) );
}

localIndex PetscMatrix::globalRowLength( globalIndex globalRowIndex ) const
{
  GEOSX_LAI_ASSERT( assembled() );
  PetscInt ncols;
  GEOSX_LAI_CHECK_ERROR( MatGetRow( m_mat, globalRowIndex, &ncols, nullptr, nullptr ) );
  localIndex const nnz = ncols;
  GEOSX_LAI_CHECK_ERROR( MatRestoreRow( m_mat, globalRowIndex, &ncols, nullptr, nullptr ) );
  return nnz;
}

void PetscMatrix::getRowCopy( globalIndex globalRow,
                              arraySlice1d< globalIndex > const & colIndices,
                              arraySlice1d< real64 > const & values ) const
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT_GE( globalRow, ilower() );
  GEOSX_LAI_ASSERT_GT( iupper(), globalRow );

  PetscScalar const * vals;
  PetscInt const * inds;
  PetscInt numEntries;

  GEOSX_LAI_CHECK_ERROR( MatGetRow( m_mat, globalRow, &numEntries, &inds, &vals ) );

  GEOSX_LAI_ASSERT_GE( colIndices.size(), numEntries );
  GEOSX_LAI_ASSERT_GE( values.size(), numEntries );

  for( int i = 0; i < numEntries; i++ )
  {
    colIndices[i] = inds[i];
  }
  for( int i = 0; i < numEntries; i++ )
  {
    values[i] = vals[i];
  }

  GEOSX_LAI_CHECK_ERROR( MatRestoreRow( m_mat, globalRow, &numEntries, &inds, &vals ) );
}

real64 PetscMatrix::getDiagValue( globalIndex globalRow ) const
{
  GEOSX_LAI_ASSERT( assembled() );
  GEOSX_LAI_ASSERT_GE( globalRow, ilower());
  GEOSX_LAI_ASSERT_GT( iupper(), globalRow );

  PetscScalar const * vals = nullptr;
  PetscInt const * cols = nullptr;
  PetscInt ncols;
  real64 diagValue = 0.0;

  GEOSX_LAI_CHECK_ERROR( MatGetRow( m_mat, globalRow, &ncols, &cols, &vals ) );
  for( PetscInt i = 0; i < ncols; i++ )
  {
    if( cols[i] == globalRow )
    {
      diagValue = vals[i];
      break;
    }
  }
  GEOSX_LAI_CHECK_ERROR( MatRestoreRow( m_mat, globalRow, &ncols, &cols, &vals ) );

  return diagValue;
}

Mat & PetscMatrix::unwrapped()
{
  GEOSX_LAI_ASSERT( created() );
  return m_mat;
}

const Mat & PetscMatrix::unwrapped() const
{
  GEOSX_LAI_ASSERT( created() );
  return m_mat;
}

globalIndex PetscMatrix::numGlobalRows() const
{
  GEOSX_LAI_ASSERT( created() );
  PetscInt num_rows;
  PetscInt num_cols;
  GEOSX_LAI_CHECK_ERROR( MatGetSize( m_mat, &num_rows, &num_cols ) );
  return integer_conversion< globalIndex >( num_rows );
}

globalIndex PetscMatrix::numGlobalCols() const
{
  GEOSX_LAI_ASSERT( created() );
  PetscInt num_rows;
  PetscInt num_cols;
  GEOSX_LAI_CHECK_ERROR( MatGetSize( m_mat, &num_rows, &num_cols ) );
  return integer_conversion< globalIndex >( num_cols );
}

globalIndex PetscMatrix::ilower() const
{
  GEOSX_LAI_ASSERT( created() );
  PetscInt firstrow;
  PetscInt lastrow;
  GEOSX_LAI_CHECK_ERROR( MatGetOwnershipRange( m_mat, &firstrow, &lastrow ) );
  return integer_conversion< globalIndex >( firstrow );
}

globalIndex PetscMatrix::iupper() const
{
  GEOSX_LAI_ASSERT( created() );
  PetscInt firstrow;
  PetscInt lastrow;
  GEOSX_LAI_CHECK_ERROR( MatGetOwnershipRange( m_mat, &firstrow, &lastrow ) );
  return integer_conversion< globalIndex >( lastrow );
}

localIndex PetscMatrix::numLocalNonzeros() const
{
  GEOSX_LAI_ASSERT( assembled() );
  PetscInt firstrow, lastrow;
  GEOSX_LAI_CHECK_ERROR( MatGetOwnershipRange( m_mat, &firstrow, &lastrow ) );

  PetscInt numEntries;
  localIndex result = 0;

  // loop over rows
  for( PetscInt row = firstrow; row < lastrow; ++row )
  {
    GEOSX_LAI_CHECK_ERROR( MatGetRow( m_mat, row, &numEntries, nullptr, nullptr ) );
    result += numEntries;
    GEOSX_LAI_CHECK_ERROR( MatRestoreRow( m_mat, row, &numEntries, nullptr, nullptr ) );
  }

  return result;
}

globalIndex PetscMatrix::numGlobalNonzeros() const
{
  return MpiWrapper::Sum( integer_conversion< globalIndex >( numLocalNonzeros() ), getComm() );
}

real64 PetscMatrix::normInf() const
{
  GEOSX_LAI_ASSERT( ready() );
  real64 normInf;
  GEOSX_LAI_CHECK_ERROR( MatNorm( m_mat, NORM_INFINITY, &normInf ) );
  return normInf;
}

real64 PetscMatrix::norm1() const
{
  GEOSX_LAI_ASSERT( ready() );
  real64 norm1;
  GEOSX_LAI_CHECK_ERROR( MatNorm( m_mat, NORM_1, &norm1 ) );
  return norm1;
}

real64 PetscMatrix::normFrobenius() const
{
  GEOSX_LAI_ASSERT( ready() );
  real64 normFrob;
  GEOSX_LAI_CHECK_ERROR( MatNorm( m_mat, NORM_FROBENIUS, &normFrob ) );
  return normFrob;
}

localIndex PetscMatrix::getLocalRowID( globalIndex const index ) const
{
  GEOSX_LAI_ASSERT( created() );
  PetscInt low, high;
  GEOSX_LAI_CHECK_ERROR( MatGetOwnershipRange( m_mat, &low, &high ) );
  return (index >= low && index < high) ? integer_conversion< localIndex >( index - low ) : -1;
}

globalIndex PetscMatrix::getGlobalRowID( localIndex const index ) const
{
  GEOSX_LAI_ASSERT( created() );
  GEOSX_LAI_ASSERT_GE( index, 0 );
  GEOSX_LAI_ASSERT_GT( numLocalRows(), index );
  PetscInt low, high;
  GEOSX_LAI_CHECK_ERROR( MatGetOwnershipRange( m_mat, &low, &high ) );
  return integer_conversion< globalIndex >( index + low );
}

localIndex PetscMatrix::numLocalCols() const
{
  GEOSX_LAI_ASSERT( created() );
  PetscInt cols;
  GEOSX_LAI_CHECK_ERROR( MatGetSize( m_mat, nullptr, &cols ) );
  return integer_conversion< localIndex >( cols );
}

localIndex PetscMatrix::numLocalRows() const
{
  GEOSX_LAI_ASSERT( created() );
  PetscInt low, high;
  GEOSX_LAI_CHECK_ERROR( MatGetOwnershipRange( m_mat, &low, &high ) );
  return integer_conversion< localIndex >( high - low );
}

MPI_Comm PetscMatrix::getComm() const
{
  GEOSX_LAI_ASSERT( created() );
  MPI_Comm comm;
  GEOSX_LAI_CHECK_ERROR( PetscObjectGetComm( reinterpret_cast< PetscObject >( m_mat ), &comm ) );
  return comm;
}

void PetscMatrix::print( std::ostream & os ) const
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_ERROR_IF( &os != &std::cout, "Only output to stdout currently supported" );
  GEOSX_LAI_CHECK_ERROR( MatView( m_mat, PETSC_VIEWER_STDOUT_( getComm() ) ) );
}

void PetscMatrix::write( string const & filename,
                         LAIOutputFormat const format ) const
{
  GEOSX_LAI_ASSERT( ready() );

  PetscViewer viewer;
  GEOSX_LAI_CHECK_ERROR( PetscViewerASCIIOpen( getComm(), filename.c_str(), &viewer ) );
  PetscViewerFormat petscFormat = PETSC_VIEWER_DEFAULT;

  switch( format )
  {
    case LAIOutputFormat::NATIVE_ASCII:
    {
      petscFormat = PETSC_VIEWER_ASCII_COMMON;
    }
    break;
    case LAIOutputFormat::NATIVE_BINARY:
    {
      petscFormat = PETSC_VIEWER_NATIVE;
    }
    break;
    case LAIOutputFormat::MATLAB_ASCII:
    {
      petscFormat = PETSC_VIEWER_ASCII_MATLAB;
    }
    break;
    case LAIOutputFormat::MATLAB_BINARY:
    {
      petscFormat = PETSC_VIEWER_BINARY_MATLAB;
    }
    break;
    case LAIOutputFormat::MATRIX_MARKET:
    {
      petscFormat = PETSC_VIEWER_ASCII_MATRIXMARKET;
    }
    break;
    default:
      GEOSX_ERROR( "Unsupported matrix output format" );
  }

  GEOSX_LAI_CHECK_ERROR( PetscViewerPushFormat( viewer, petscFormat ) );
  GEOSX_LAI_CHECK_ERROR( MatView( m_mat, viewer ) );
  GEOSX_LAI_CHECK_ERROR( PetscViewerDestroy( &viewer ) );
}

} // end geosx namespace
