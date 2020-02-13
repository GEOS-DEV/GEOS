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

// BEGIN_RST_NARRATIVE PetscSparseMatrix.rst
// ==============================
// PETSc-based Matrix Object
// ==============================
// This class contains the ParallelMatrix wrappers based on PETSc Mat objects.
// The class contains a unique pointer to a Mat as well as constructors,
// functions and accessors for Mat objects.

// Include the corresponding header file.
#include "PetscSparseMatrix.hpp"

#if !defined(PETSC_USE_64BIT_INDICES)
#define PETSC_USE_64BIT_INDICES
#endif

#include <petscvec.h>
#include <petscmat.h>

#include "mpiCommunications/MpiWrapper.hpp"

// Put everything under the geosx namespace.
namespace geosx
{

static_assert( sizeof(PetscInt) == sizeof(globalIndex), "sizeof(PetscInt) != sizeof(localIndex)");
static_assert( std::is_same<PetscScalar, real64>::value, "PetscScalar != real64" );

inline PetscInt * toPetscInt( globalIndex * const index )
{
  return reinterpret_cast<PetscInt*>(index);
}

inline PetscInt const * toPetscInt( globalIndex const * const index )
{
  return reinterpret_cast<PetscInt const*>(index);
}

// ----------------------------
// Constructors
// ----------------------------

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Create an empty matrix (meant to be used for declaration)
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
PetscSparseMatrix::PetscSparseMatrix()
: Base(),
  m_mat{}
{}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Copy constructor
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
PetscSparseMatrix::PetscSparseMatrix( PetscSparseMatrix const & src )
: PetscSparseMatrix()
{
  GEOSX_ASSERT( !src.isOpen() );
  GEOSX_ASSERT( src.isAssembled() );

  MatDuplicate( src.m_mat, MAT_COPY_VALUES, &m_mat );
  m_assembled = true;
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Destructor
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
PetscSparseMatrix::~PetscSparseMatrix()
{
  MatDestroy( &m_mat );
}

// -----------------------------
// Create
// -----------------------------
// Allocate matrix (prepare to be filled with data).
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Create a matrix from number of elements
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
void PetscSparseMatrix::createWithLocalSize( localIndex const localSize,
                                             localIndex const maxEntriesPerRow,
                                             MPI_Comm const & comm )
{
  // call general version
  createWithLocalSize( localSize, localSize, maxEntriesPerRow, comm );
}

void PetscSparseMatrix::createWithGlobalSize( globalIndex const globalSize,
                                              localIndex const maxEntriesPerRow,
                                              MPI_Comm const & comm )
{
  // call general version
  createWithGlobalSize( globalSize, globalSize, maxEntriesPerRow, comm );
}

void PetscSparseMatrix::createWithLocalSize( localIndex const localRows,
                                             localIndex const localCols,
                                             localIndex const maxEntriesPerRow,
                                             MPI_Comm const & comm )
{
  GEOSX_ASSERT( !isOpen() );
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
  GEOSX_ASSERT( !isOpen() );
  reset();

  // set up matrix
  MatCreate( comm, &m_mat );
  MatSetType( m_mat, MATMPIAIJ );
  MatSetSizes( m_mat, PETSC_DECIDE, PETSC_DECIDE, globalRows, globalCols );
  MatMPIAIJSetPreallocation( m_mat, maxEntriesPerRow, nullptr, maxEntriesPerRow, nullptr );
  MatSetUp( m_mat );
}

bool PetscSparseMatrix::isCreated() const
{
  return m_mat != nullptr;
}

void PetscSparseMatrix::reset()
{
  MatrixBase::reset();
  MatDestroy( &m_mat );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Reinitialize.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Sets all values to user-defined value.
void PetscSparseMatrix::set( real64 const value )
{
  GEOSX_ASSERT( !isOpen() );
  GEOSX_ASSERT( isAssembled() );

  PetscInt firstrow;
  PetscInt lastrow;
  MatGetOwnershipRange( m_mat, &firstrow, &lastrow );

  PetscInt numEntries;
  const PetscInt * inds;
  const PetscScalar * vals;

  PetscInt numEntries_;
  array1d<PetscScalar> vals_;
  array1d<PetscInt> inds_;

  PetscInt const maxNumRows = MpiWrapper::Max( lastrow - firstrow );

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

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Reinitialize.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Sets all values to 0.
void PetscSparseMatrix::zero()
{
  GEOSX_ASSERT( !isOpen() );
  GEOSX_ASSERT( isAssembled() );
  MatZeroEntries( m_mat );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Open
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Empty open function (implemented for HYPRE compatibility).
void PetscSparseMatrix::open()
{
  GEOSX_ASSERT( !isOpen() );
  GEOSX_ASSERT_MSG( isCreated(), "The matrix has not been created" );
  m_open = true;
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Close
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// PETSc matrix assembly. Space for preallocated nonzeros that is not filled are compressed out by assembly.
void PetscSparseMatrix::close()
{
  GEOSX_ASSERT( isOpen() );
  MatAssemblyBegin( m_mat, MAT_FINAL_ASSEMBLY );
  MatAssemblyEnd( m_mat, MAT_FINAL_ASSEMBLY );
  MatSetOption( m_mat, MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE );
  MatSetOption( m_mat, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE );
  m_assembled = true;
  m_open = false;
}

// -------------------------
// Add/Set
// -------------------------

// 1x1
void PetscSparseMatrix::add( globalIndex const rowIndex,
                             globalIndex const colIndex,
                             real64 const value )
{
  GEOSX_ASSERT( isOpen() );
  GEOSX_ASSERT( isAssembled() );
  MatSetValue( m_mat, rowIndex, colIndex, value, ADD_VALUES );
}

void PetscSparseMatrix::set( globalIndex const rowIndex,
                             globalIndex const colIndex,
                             real64 const value )
{
  GEOSX_ASSERT( isOpen() );
  GEOSX_ASSERT( isAssembled() );
  MatSetValue( m_mat, rowIndex, colIndex, value, INSERT_VALUES );
}

void PetscSparseMatrix::insert( globalIndex const rowIndex,
                                globalIndex const colIndex,
                                real64 const value )
{
  GEOSX_ASSERT( isOpen() );
  GEOSX_ASSERT( !isAssembled() );
  MatSetValue( m_mat, rowIndex, colIndex, value, INSERT_VALUES );
}

 // 1xN c-style
void PetscSparseMatrix::add( globalIndex const rowIndex,
                             globalIndex const * colIndices,
                             real64 const * values,
                             localIndex size )
{
  GEOSX_ASSERT( isOpen() );
  GEOSX_ASSERT( isAssembled() );
  PetscInt rows[1] = {rowIndex};
  MatSetValues( m_mat, 1, rows, size, toPetscInt( colIndices ), values, ADD_VALUES );
}

void PetscSparseMatrix::set( globalIndex const rowIndex,
                             globalIndex const * colIndices,
                             real64 const * values,
                             localIndex size )
{
  GEOSX_ASSERT( isOpen() );
  GEOSX_ASSERT( isAssembled() );
  PetscInt rows[1] = {rowIndex};
  MatSetValues( m_mat, 1, rows, size, toPetscInt( colIndices), values, INSERT_VALUES );
}

void PetscSparseMatrix::insert( globalIndex const rowIndex,
                                globalIndex const * colIndices,
                                real64 const * values,
                                localIndex size )
{
  GEOSX_ASSERT( isOpen() );
  GEOSX_ASSERT( !isAssembled() );
  PetscInt rows[1] = {rowIndex};
  MatSetValues( m_mat, 1, rows, size, toPetscInt( colIndices), values, INSERT_VALUES );
}

// 1xN array1d style 
void PetscSparseMatrix::add( globalIndex const rowIndex,
                             arraySlice1d<globalIndex const> const &colIndices,
                             arraySlice1d<real64 const> const &values )
{
  GEOSX_ASSERT( isOpen() );
  GEOSX_ASSERT( isAssembled() );
  PetscInt rows[1] = {rowIndex};
  MatSetValues( m_mat, 1, rows, values.size(), toPetscInt( colIndices.data()), values.data(), ADD_VALUES );
}

void PetscSparseMatrix::set( globalIndex const rowIndex,
                             arraySlice1d<globalIndex const> const &colIndices,
                             arraySlice1d<real64 const> const &values )
{
  GEOSX_ASSERT( isOpen() );
  GEOSX_ASSERT( isAssembled() );
  PetscInt rows[1] = {rowIndex};
  MatSetValues( m_mat, 1, rows, values.size(), toPetscInt( colIndices.data()), values.data(), INSERT_VALUES );
}

void PetscSparseMatrix::insert( globalIndex const rowIndex,
                                arraySlice1d<globalIndex const> const &colIndices,
                                arraySlice1d<real64 const> const &values )
{
  GEOSX_ASSERT( isOpen() );
  GEOSX_ASSERT( !isAssembled() );
  PetscInt rows[1] = {rowIndex};
  MatSetValues( m_mat, 1, rows, values.size(), toPetscInt( colIndices.data()), values.data(), INSERT_VALUES );
}

// MxN array2d style
void PetscSparseMatrix::add( arraySlice1d<globalIndex const> const & rowIndices,
                             arraySlice1d<globalIndex const> const & colIndices,
                             arraySlice2d<real64 const, 1> const & values )
{
  GEOSX_ASSERT( isOpen() );
  GEOSX_ASSERT( isAssembled() );
  MatSetValues( m_mat,
                rowIndices.size(),
                toPetscInt(rowIndices.data()),
                colIndices.size(),
                toPetscInt(colIndices.data()),
                values.data(),
                ADD_VALUES );
}

void PetscSparseMatrix::set( arraySlice1d<globalIndex const> const & rowIndices,
                             arraySlice1d<globalIndex const> const & colIndices,
                             arraySlice2d<real64 const, 1> const & values )
{
  GEOSX_ASSERT( isOpen() );
  GEOSX_ASSERT( isAssembled() );
  MatSetValues( m_mat,
                rowIndices.size(),
                toPetscInt(rowIndices.data()),
                colIndices.size(),
                toPetscInt(colIndices.data()),
                values.data(),
                INSERT_VALUES );
}

void PetscSparseMatrix::insert( arraySlice1d<globalIndex const> const & rowIndices,
                                arraySlice1d<globalIndex const> const & colIndices,
                                arraySlice2d<real64 const, 1> const & values )
{
  GEOSX_ASSERT( isOpen() );
  GEOSX_ASSERT( !isAssembled() );
  MatSetValues( m_mat,
                rowIndices.size(),
                toPetscInt(rowIndices.data()),
                colIndices.size(),
                toPetscInt(colIndices.data()),
                values.data(),
                INSERT_VALUES );
}

void PetscSparseMatrix::add( arraySlice1d<globalIndex const> const & rowIndices,
                             arraySlice1d<globalIndex const> const & colIndices,
                             arraySlice2d<real64 const, 0> const & values )
{
  GEOSX_ASSERT( isOpen() );
  GEOSX_ASSERT( isAssembled() );
  MatSetOption( m_mat, MAT_ROW_ORIENTED, PETSC_FALSE );
  MatSetValues( m_mat,
                rowIndices.size(),
                toPetscInt(rowIndices.data()),
                colIndices.size(),
                toPetscInt(colIndices.data()),
                values.data(),
                ADD_VALUES );
  MatSetOption( m_mat, MAT_ROW_ORIENTED, PETSC_TRUE );
}

void PetscSparseMatrix::set( arraySlice1d<globalIndex const> const & rowIndices,
                             arraySlice1d<globalIndex const> const & colIndices,
                             arraySlice2d<real64 const, 0> const & values )
{
  GEOSX_ASSERT( isOpen() );
  GEOSX_ASSERT( isAssembled() );
  MatSetOption( m_mat, MAT_ROW_ORIENTED, PETSC_FALSE );
  MatSetValues( m_mat,
                rowIndices.size(),
                toPetscInt(rowIndices.data()),
                colIndices.size(),
                toPetscInt(colIndices.data()),
                values.data(),
                INSERT_VALUES );
  MatSetOption( m_mat, MAT_ROW_ORIENTED, PETSC_TRUE );
}

void PetscSparseMatrix::insert( arraySlice1d<globalIndex const> const & rowIndices,
                                arraySlice1d<globalIndex const> const & colIndices,
                                arraySlice2d<real64 const, 0> const & values )
{
  GEOSX_ASSERT( isOpen() );
  GEOSX_ASSERT( !isAssembled() );
  MatSetOption( m_mat, MAT_ROW_ORIENTED, PETSC_FALSE );
  MatSetValues( m_mat,
                rowIndices.size(),
                toPetscInt(rowIndices.data()),
                colIndices.size(),
                toPetscInt(colIndices.data()),
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
  GEOSX_ASSERT( isOpen() );
  GEOSX_ASSERT( isAssembled() );
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
  GEOSX_ASSERT( isOpen() );
  GEOSX_ASSERT( isAssembled() );
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
  GEOSX_ASSERT( isOpen() );
  GEOSX_ASSERT( !isAssembled() );
  MatSetValues( m_mat,
                numRows,
                toPetscInt(rowIndices),
                numCols,
                toPetscInt(colIndices),
                values,
                INSERT_VALUES );
}

// -------------------------
// Linear Algebra
// -------------------------

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Matrix/vector multiplication
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Perform the matrix-vector product A*src = dst.
//
// NOTE: src and dst must be different vectors
void PetscSparseMatrix::multiply( PetscVector const &src,
                                  PetscVector &dst ) const
{
  GEOSX_ASSERT( !isOpen() );
  GEOSX_ASSERT( isAssembled() );
  MatMult( m_mat, src.getConstVec(), dst.getVec() );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Matrix/matrix multiplication
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Perform the matrix-matrix product A*src = dst.
//
// NOTE: src and dst must be different vectors
void PetscSparseMatrix::multiply( PetscSparseMatrix const & src, 
                                  PetscSparseMatrix & dst,
                                  bool const closeResult ) const
{
  GEOSX_ASSERT( !isOpen() );
  GEOSX_ASSERT( isAssembled() );
  GEOSX_ASSERT( !src.isOpen() );
  GEOSX_ASSERT( src.isAssembled() );

  MatReuse const reuse = dst.isCreated() ? MAT_REUSE_MATRIX : MAT_INITIAL_MATRIX;
  if( !dst.isCreated() )
  {
    dst.createWithLocalSize( localRows(), src.localCols(), 1, getComm() );
  }

  MatMatMult( m_mat, src.unwrappedPointer(), reuse, PETSC_DEFAULT, &dst.unwrappedPointer() );
  dst.m_assembled = closeResult;
  dst.m_open = !closeResult;
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Matrix/matrix multiplication
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Perform the matrix-matrix product this^T * src = dst.
void PetscSparseMatrix::leftMultiplyTranspose( PetscSparseMatrix const & src,
                                               PetscSparseMatrix & dst,
                                               bool const GEOSX_UNUSED_ARG( closeResult ) ) const
{
  GEOSX_ASSERT( !isOpen() );
  GEOSX_ASSERT( isAssembled() );
  GEOSX_ASSERT( !src.isOpen() );
  GEOSX_ASSERT( src.isAssembled() );

  MatReuse const reuse = dst.isCreated() ? MAT_REUSE_MATRIX : MAT_INITIAL_MATRIX;
  if( !dst.isCreated() )
  {
    dst.createWithLocalSize( localCols(), src.localCols(), 1, getComm() );
  }

  MatTransposeMatMult( m_mat, src.unwrappedPointer(), reuse, PETSC_DEFAULT, &dst.unwrappedPointer() );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Matrix/matrix multiplication
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Perform the matrix-matrix product src * this^T  = dst.
void PetscSparseMatrix::rightMultiplyTranspose( PetscSparseMatrix const & src,
                                                PetscSparseMatrix & dst,
                                                bool const GEOSX_UNUSED_ARG( closeResult ) ) const
{
  GEOSX_ASSERT( !isOpen() );
  GEOSX_ASSERT( isAssembled() );
  GEOSX_ASSERT( !src.isOpen() );
  GEOSX_ASSERT( src.isAssembled() );

  MatReuse const reuse = dst.isCreated() ? MAT_REUSE_MATRIX : MAT_INITIAL_MATRIX;
  if( !dst.isCreated() )
  {
    dst.createWithLocalSize( localRows(), src.localRows(), 1, getComm() );
  }

  MatMatTransposeMult( m_mat, src.unwrappedPointer(), reuse, PETSC_DEFAULT, &dst.unwrappedPointer() );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Generalized matrix/vector product.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Compute gemv <tt>y = alpha*A*x + beta*y</tt>.
void PetscSparseMatrix::gemv( real64 const alpha,
                              PetscVector const & x,
                              real64 const beta,
                              PetscVector & y,
                              bool useTranspose ) const
{
  GEOSX_ASSERT( !isOpen() );
  GEOSX_ASSERT( isAssembled() );

  PetscVector x_( x );
  PetscVector b_( x );

  x_.scale( alpha ); // alpha*x_
  y.scale( beta ); // beta*y

  if ( useTranspose ) 
  {
    MatMultTranspose( m_mat, x_.getVec(), b_.getVec() );
  }
  else
  {
    MatMult( m_mat, x_.getVec(), b_.getVec() ); // alpha*A*x_ = b_
  }
  VecAXPY( y.getVec(), 1, b_.getVec() ); // alpha*A*x_ + beta*y = y
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Scale.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Multiply all elements by scalingFactor.
void PetscSparseMatrix::scale( real64 const scalingFactor )
{
  GEOSX_ASSERT( !isOpen() );
  GEOSX_ASSERT( isAssembled() );
  MatScale( m_mat, scalingFactor );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Left and right scaling
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
void PetscSparseMatrix::leftScale( PetscVector const &vec )
{
  GEOSX_ASSERT( !isOpen() );
  GEOSX_ASSERT( isAssembled() );
  MatDiagonalScale( m_mat, vec.getConstVec(), nullptr );
}

void PetscSparseMatrix::rightScale( PetscVector const &vec )
{
  GEOSX_ASSERT( !isOpen() );
  GEOSX_ASSERT( isAssembled() );
  MatDiagonalScale( m_mat, nullptr, vec.getConstVec() );
}

void PetscSparseMatrix::leftRightScale( PetscVector const &vecLeft,
                                        PetscVector const &vecRight )
{
  GEOSX_ASSERT( !isOpen() );
  GEOSX_ASSERT( isAssembled() );
  MatDiagonalScale( m_mat, vecLeft.getConstVec(), vecRight.getConstVec() );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get global row copy
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// The challenge here is that columns are stored with local, not global,
// indices, so we need to do conversions back and forth
void PetscSparseMatrix::getRowCopy( globalIndex globalRow,
                                    array1d<globalIndex> & colIndices,
                                    array1d<real64> & values ) const
{
  GEOSX_ASSERT( !isOpen() );
  GEOSX_ASSERT( isAssembled() );
  GEOSX_ASSERT_GE( globalRow, ilower() );
  GEOSX_ASSERT_GT( iupper(), globalRow );

  const PetscScalar* vals;
  const PetscInt* inds;
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
  GEOSX_ASSERT( isAssembled() );
  GEOSX_ASSERT_GE( globalRow, ilower());
  GEOSX_ASSERT_GT( iupper(), globalRow );

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

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Clear row.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Clear the row.  By default the diagonal value will be set
// to zero, but the user can pass a desired diagValue
void PetscSparseMatrix::clearRow( globalIndex const globalRow,
                                  real64 const diagValue )
{
  GEOSX_ASSERT( isOpen() );
  GEOSX_ASSERT( isAssembled() );
  GEOSX_ASSERT_GE( globalRow, ilower() );
  GEOSX_ASSERT_GT( iupper(), globalRow );

  PetscInt rows[1] = {globalRow};
  MatZeroRows( m_mat, 1, rows, diagValue, nullptr, nullptr );
}

// ----------------------------
//  Accessors
// ----------------------------

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get pointer.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""

// Accessor for the pointer to the matrix
Mat & PetscSparseMatrix::unwrappedPointer()
{
  return m_mat;
}

// Accessor for the pointer to the matrix
const Mat & PetscSparseMatrix::unwrappedPointer() const
{
  return m_mat;
}

// Accessor for the MPI communicator
MPI_Comm PetscSparseMatrix::getComm() const
{
  GEOSX_ASSERT_MSG( isCreated(), "Matrix has not been created" );
  MPI_Comm comm;
  PetscObjectGetComm( reinterpret_cast<PetscObject>( m_mat ), &comm );
  return comm;
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get number of global rows.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Accessor for the number of global rows
globalIndex PetscSparseMatrix::globalRows() const
{
  GEOSX_ASSERT_MSG( isCreated(), "Matrix has not been created" );
  PetscInt num_rows;
  PetscInt num_cols;
  MatGetSize( m_mat, &num_rows, &num_cols );
  return integer_conversion<globalIndex>( num_rows );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get number of global columns.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Accessor for the number of global columns
globalIndex PetscSparseMatrix::globalCols() const
{
  GEOSX_ASSERT_MSG( isCreated(), "Matrix has not been created" );
  PetscInt num_rows;
  PetscInt num_cols;
  MatGetSize( m_mat, &num_rows, &num_cols );
  return integer_conversion<globalIndex>( num_cols );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get the lower index owned by processor.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Accessor for the index of the first global row
globalIndex PetscSparseMatrix::ilower() const
{
  GEOSX_ASSERT_MSG( isCreated(), "Matrix has not been created" );
  PetscInt firstrow;
  PetscInt lastrow;
  MatGetOwnershipRange( m_mat, &firstrow, &lastrow );
  return integer_conversion<globalIndex>( firstrow );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get the upper index owned by processor.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Accessor for the index of the last global row
globalIndex PetscSparseMatrix::iupper() const
{
  GEOSX_ASSERT_MSG( isCreated(), "Matrix has not been created" );
  PetscInt firstrow;
  PetscInt lastrow;
  MatGetOwnershipRange( m_mat, &firstrow, &lastrow );
  return integer_conversion<globalIndex>( lastrow );
 } 

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Print to terminal.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Wrapper to print the petsc output of the matrix
void PetscSparseMatrix::print( std::ostream & os ) const
{
  GEOSX_ASSERT_MSG( isCreated(), "Matrix has not been created" );
  GEOSX_ERROR_IF( &os != &std::cout, "Only output to stdout currently supported" );
  MatView( m_mat, PETSC_VIEWER_STDOUT_( getComm() ) );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Write to matlab-compatible file
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
void PetscSparseMatrix::write( string const & filename,
                               MatrixOutputFormat const format ) const
{
  GEOSX_ASSERT( isCreated() );

  PetscViewer viewer;
  PetscViewerASCIIOpen( getComm(), filename.c_str(), &viewer );
  PetscViewerFormat petscFormat;

  switch( format )
  {
    case MatrixOutputFormat::NATIVE_ASCII:
      petscFormat = PETSC_VIEWER_DEFAULT;
      break;
    case MatrixOutputFormat::NATIVE_BINARY:
      petscFormat = PETSC_VIEWER_NATIVE;
      break;
    case MatrixOutputFormat::MATLAB_ASCII:
      petscFormat = PETSC_VIEWER_ASCII_MATLAB;
      break;
    case MatrixOutputFormat::MATLAB_BINARY:
      petscFormat = PETSC_VIEWER_BINARY_MATLAB;
      break;
    case MatrixOutputFormat::MATRIX_MARKET:
      petscFormat = PETSC_VIEWER_ASCII_MATRIXMARKET;
      break;
    default:
      GEOSX_ERROR( "Unsupported matrix output format" );
  }

  PetscViewerPushFormat( viewer, petscFormat );
  MatView( m_mat, viewer );
  PetscViewerDestroy( &viewer );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Inf-norm.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Returns the infinity norm of the matrix.
real64 PetscSparseMatrix::normInf() const
{
  GEOSX_ASSERT( !isOpen() );
  GEOSX_ASSERT( isAssembled() );
  real64 normInf;
  MatNorm( m_mat, NORM_INFINITY, &normInf );
  return normInf;
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// 1-norm.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Returns the one norm of the matrix.
real64 PetscSparseMatrix::norm1() const
{
  GEOSX_ASSERT( !isOpen() );
  GEOSX_ASSERT( isAssembled() );
  real64 norm1;
  MatNorm( m_mat, NORM_1, &norm1 );
  return norm1;
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Frobenius-norm.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Returns the Frobenius norm of the matrix.
real64 PetscSparseMatrix::normFrobenius() const
{
  GEOSX_ASSERT( !isOpen() );
  GEOSX_ASSERT( isAssembled() );
  real64 normFrob;
  MatNorm( m_mat, NORM_FROBENIUS, &normFrob );
  return normFrob;
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// getLocalRowID
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Map a global row index to local row index
localIndex PetscSparseMatrix::getLocalRowID( globalIndex const index ) const
{
  GEOSX_ASSERT( isCreated() );
  GEOSX_ASSERT_GE( index, ilower() );
  GEOSX_ASSERT_GT( iupper(), index );
  PetscInt low, high;
  MatGetOwnershipRange( m_mat, &low, &high);
  return integer_conversion<localIndex>( index - low );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// getGlobalRowID
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Map a local row index to global row index
globalIndex PetscSparseMatrix::getGlobalRowID( localIndex const index ) const
{
  GEOSX_ASSERT( isCreated() );
  GEOSX_ASSERT_GE( index, 0 );
  GEOSX_ASSERT_GT( localRows(), index );
  PetscInt low, high;
  MatGetOwnershipRange( m_mat, &low, &high);
  return integer_conversion<globalIndex>( index + low );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// localCols
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Return the local number of columns on each processor
//
// NOTE: PETSc MPI matrices are partitioned row-wise so that the local number
// of columns is the global number.
localIndex PetscSparseMatrix::localCols() const
{
  GEOSX_ASSERT( isCreated());
  PetscInt cols;
  MatGetSize( m_mat, nullptr, &cols );
  return integer_conversion<localIndex>( cols );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// localRows
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Return the local number of rows on each processor
localIndex PetscSparseMatrix::localRows() const
{
  GEOSX_ASSERT( isCreated() );
  PetscInt low, high;
  MatGetOwnershipRange( m_mat, &low, &high);
  return integer_conversion<localIndex >( high - low );
}

localIndex PetscSparseMatrix::localNonzeros() const
{
  GEOSX_ASSERT( isAssembled() );
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

globalIndex PetscSparseMatrix::globalNonzeros() const
{
  GEOSX_ASSERT( isAssembled() );
  return MpiWrapper::Sum( integer_conversion<globalIndex>( localNonzeros() ) );
}

} // end geosx namespace
