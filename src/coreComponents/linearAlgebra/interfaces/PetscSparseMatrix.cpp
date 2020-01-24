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
{}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Copy constructor
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
PetscSparseMatrix::PetscSparseMatrix( PetscSparseMatrix const &src )
{
  MatDuplicate( src.m_mat, MAT_COPY_VALUES, &m_mat );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Destructor
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
PetscSparseMatrix::~PetscSparseMatrix()
{
  //MatDestroy( &m_mat );
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
  // set up matrix
  MatCreate( comm, &m_mat );
  MatSetType( m_mat, MATMPIAIJ );
  MatSetSizes( m_mat, PETSC_DECIDE, PETSC_DECIDE, globalRows, globalCols );
  MatMPIAIJSetPreallocation( m_mat, maxEntriesPerRow, nullptr, maxEntriesPerRow, nullptr );
  MatSetUp( m_mat );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Reinitialize.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Sets all values to user-defined value.
void PetscSparseMatrix::set( real64 const value )
{
  PetscInt firstrow;
  PetscInt lastrow;
  MatGetOwnershipRange( m_mat, &firstrow, &lastrow );

  PetscInt numEntries;
  const PetscInt* inds;
  const PetscScalar* vals;

  PetscInt numEntries_;
  PetscScalar* vals_;
  PetscInt* inds_;

  // loop over rows
  for( PetscInt row = firstrow; row < lastrow; row++){

    // get entries in row
    close(); 
    MatGetRow( m_mat, row, &numEntries, &inds, &vals );
    numEntries_ = numEntries;
    inds_ = new PetscInt[numEntries_];
    for ( int i = 0; i < numEntries_; i++ ) 
    {
      inds_[i] = inds[i];
    }   
    MatRestoreRow( m_mat, row, &numEntries, &inds, &vals );
    close(); 

    // set entries to value
    if( numEntries_ > 0 ) {

      vals_ = new PetscScalar[numEntries_];
      for ( int i = 0; i < numEntries_; i++ ) 
      {
        vals_[i] = value;
      }

      PetscInt rows[1] = {row};
      MatSetValues( m_mat, 1, rows, numEntries_, inds_, vals_, INSERT_VALUES );
    }
  }
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Reinitialize.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Sets all values to 0.
void PetscSparseMatrix::zero()
{
  MatZeroEntries( m_mat );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Open
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Empty open function (implemented for HYPRE compatibility).
void PetscSparseMatrix::open()
{}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Close
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// PETSc matrix assembly. Space for preallocated nonzeros that is not filled are compressed out by assembly.
void PetscSparseMatrix::close()
{
  MatAssemblyBegin( m_mat, MAT_FINAL_ASSEMBLY );
  MatAssemblyEnd( m_mat, MAT_FINAL_ASSEMBLY );
  m_assembled = true;
}

// -------------------------
// Add/Set
// -------------------------

// 1x1
void PetscSparseMatrix::add( globalIndex const rowIndex,
                             globalIndex const colIndex,
                             real64 const value )
{
  MatSetValue( m_mat, rowIndex, colIndex, value, ADD_VALUES );
}

void PetscSparseMatrix::set( globalIndex const rowIndex,
                             globalIndex const colIndex,
                             real64 const value )
{
  MatSetValue( m_mat, rowIndex, colIndex, value, INSERT_VALUES );
}

void PetscSparseMatrix::insert( globalIndex const rowIndex,
                                globalIndex const colIndex,
                                real64 const value )
{
  MatSetValue( m_mat, rowIndex, colIndex, value, INSERT_VALUES );
}

 // 1xN c-style
void PetscSparseMatrix::add( globalIndex const rowIndex,
                             globalIndex const * colIndices,
                             real64 const * values,
                             localIndex size )
{
  PetscInt rows[1] = {rowIndex};
  MatSetValues( m_mat, 1, rows, size, toPetscInt( colIndices), values, ADD_VALUES );
}

void PetscSparseMatrix::set( globalIndex const rowIndex,
                             globalIndex const * colIndices,
                             real64 const * values,
                             localIndex size )
{ 
  PetscInt rows[1] = {rowIndex};
  MatSetValues( m_mat, 1, rows, size, toPetscInt( colIndices), values, INSERT_VALUES );
}

void PetscSparseMatrix::insert( globalIndex const rowIndex,
                                globalIndex const * colIndices,
                                real64 const * values,
                                localIndex size )
{
  PetscInt rows[1] = {rowIndex};
  MatSetValues( m_mat, 1, rows, size, toPetscInt( colIndices), values, INSERT_VALUES );
}

// 1xN array1d style 
void PetscSparseMatrix::add( globalIndex const rowIndex,
                             array1d<globalIndex> const &colIndices,
                             array1d<real64> const &values )
{
  PetscInt rows[1] = {rowIndex};
  MatSetValues( m_mat, 1, rows, values.size(), toPetscInt( colIndices.data()), values.data(), ADD_VALUES );
}

void PetscSparseMatrix::set( globalIndex const rowIndex,
                             array1d<globalIndex> const &colIndices,
                             array1d<real64> const &values )
{
  PetscInt rows[1] = {rowIndex};
  MatSetValues( m_mat, 1, rows, values.size(), toPetscInt( colIndices.data()), values.data(), INSERT_VALUES );
}

void PetscSparseMatrix::insert( globalIndex const rowIndex,
                                array1d<globalIndex> const &colIndices,
                                array1d<real64> const &values )
{
  PetscInt rows[1] = {rowIndex};
  MatSetValues( m_mat, 1, rows, values.size(), toPetscInt( colIndices.data()), values.data(), INSERT_VALUES );
}

// MxN array2d style
void PetscSparseMatrix::add( array1d<globalIndex> const & rowIndices,
                             array1d<globalIndex> const & colIndices,
                             array2d<real64> const & values )
{
  MatSetValues( m_mat,
                rowIndices.size(),
                toPetscInt(rowIndices.data()),
                colIndices.size(),
                toPetscInt(colIndices.data()),
                values.data(),
                ADD_VALUES );
}

void PetscSparseMatrix::set( array1d<globalIndex> const & rowIndices,
                             array1d<globalIndex> const & colIndices,
                             array2d<real64> const & values )
{
  MatSetValues( m_mat,
                rowIndices.size(),
                toPetscInt(rowIndices.data()),
                colIndices.size(),
                toPetscInt(colIndices.data()),
                values.data(),
                INSERT_VALUES );
}

void PetscSparseMatrix::insert( array1d<globalIndex> const & rowIndices,
                                array1d<globalIndex> const & colIndices,
                                array2d<real64> const & values )
{
  MatSetValues( m_mat,
                rowIndices.size(),
                toPetscInt(rowIndices.data()),
                colIndices.size(),
                toPetscInt(colIndices.data()),
                values.data(),
                INSERT_VALUES );
}

void PetscSparseMatrix::insert( globalIndex const * rowIndices,
                                globalIndex const * colIndices,
                                real64 const * values,
                                localIndex const numRows,
                                localIndex const numCols )
{
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
                                  bool const GEOSX_UNUSED_ARG( closeResult ) ) const
{
  MatMatMult( m_mat, src.getConstMat(), MAT_INITIAL_MATRIX, PETSC_DEFAULT, dst.unwrappedNonConstPointer() );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Matrix/matrix multiplication
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Perform the matrix-matrix product this^T * src = dst.
void PetscSparseMatrix::leftMultiplyTranspose( PetscSparseMatrix const & src,
                                               PetscSparseMatrix & dst,
                                               bool const GEOSX_UNUSED_ARG( closeResult ) ) const
{
  MatTransposeMatMult( m_mat, src.getConstMat(), MAT_INITIAL_MATRIX, PETSC_DEFAULT, dst.unwrappedNonConstPointer() );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Matrix/matrix multiplication
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Perform the matrix-matrix product src * this^T  = dst.
void PetscSparseMatrix::rightMultiplyTranspose( PetscSparseMatrix const & src,
                                                PetscSparseMatrix & dst,
                                                bool const GEOSX_UNUSED_ARG( closeResult ) ) const
{
  MatMatTransposeMult( m_mat, src.getConstMat(), MAT_INITIAL_MATRIX, PETSC_DEFAULT, dst.unwrappedNonConstPointer() );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Compute residual.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Compute res = b - Ax (residual form).
void PetscSparseMatrix::residual( PetscVector const & x,
                                  PetscVector const & b,
                                  PetscVector & r ) const
{
  MatMult( m_mat, x.getConstVec(), r.getConstVec() );
  VecAXPY( r.getVec(), -1, b.getConstVec() );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Generalized matrix/vector product.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Compute gemv <tt>y = alpha*A*x + beta*y</tt>.
void PetscSparseMatrix::gemv( real64 const alpha,
                              PetscVector const & x,
                              real64 const beta,
                              PetscVector & y,
                              bool useTranspose )
{
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
  MatScale( m_mat, scalingFactor );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Left and right scaling
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
void PetscSparseMatrix::leftScale( PetscVector const &vec )
{
  MatDiagonalScale( m_mat, vec.getConstVec(), nullptr );
}

void PetscSparseMatrix::rightScale( PetscVector const &vec )
{
  MatDiagonalScale( m_mat, nullptr, vec.getConstVec() );
}

void PetscSparseMatrix::leftRightScale( PetscVector const &vecLeft,
                                        PetscVector const &vecRight )
{
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
  GEOSX_ERROR_IF( !m_assembled, "Attempting to call " << __FUNCTION__ << " before close() is illegal" );

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
  PetscInt rows[1] = {globalRow};

  // zero row and set diagonal to diagValue
  MatZeroRows( m_mat, 1, rows, diagValue, nullptr, nullptr );
}

// ----------------------------
//  Accessors
// ----------------------------

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get pointer.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""

// Accessor for the pointer to the matrix
Mat * PetscSparseMatrix::unwrappedNonConstPointer()
{
  return &( m_mat );
}

// Accessor for the pointer to the matrix
const Mat * PetscSparseMatrix::unwrappedPointer() const
{
  return &( m_mat );
}

// Accessor for the MPI communicator
MPI_Comm PetscSparseMatrix::getComm() const
{
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
  PetscInt num_rows;
  PetscInt num_cols;
  MatGetSize( m_mat, &num_rows, &num_cols );
  return num_rows;
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get number of global columns.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Accessor for the number of global columns
globalIndex PetscSparseMatrix::globalCols() const
{
  PetscInt num_rows;
  PetscInt num_cols;
  MatGetSize( m_mat, &num_rows, &num_cols );
  return num_cols;
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get the lower index owned by processor.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Accessor for the index of the first global row
globalIndex PetscSparseMatrix::ilower() const
{
  PetscInt firstrow;
  PetscInt lastrow;
  MatGetOwnershipRange( m_mat, &firstrow, &lastrow );
  return firstrow; 
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get the upper index owned by processor.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Accessor for the index of the last global row
globalIndex PetscSparseMatrix::iupper() const
{
  PetscInt firstrow;
  PetscInt lastrow;
  MatGetOwnershipRange( m_mat, &firstrow, &lastrow );
  return lastrow; 
 } 

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Print to terminal.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Wrapper to print the petsc output of the matrix
void PetscSparseMatrix::print() const 
{
  MatView( m_mat, PETSC_VIEWER_STDOUT_WORLD );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Write to matlab-compatible file
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
void PetscSparseMatrix::write( string const & filename,
                               bool const mtxFormat ) const
{
  PetscViewer viewer;

  if( mtxFormat )
  {
    // ".mtx" extension
    string name( filename );
    if( filename.substr( filename.find_last_of( "." ) + 1 ) != "mtx" ) 
    {
      name = filename.substr( 0, filename.find_last_of( "." ) ) + ".mtx";
    }
    PetscViewerASCIIOpen( getComm(), name.c_str(), &viewer);
    PetscViewerPushFormat( viewer, PETSC_VIEWER_ASCII_MATRIXMARKET );
  }
  else
  {
    PetscViewerASCIIOpen( getComm(), filename.c_str(), &viewer);
    PetscViewerPushFormat( viewer, PETSC_VIEWER_ASCII_MATLAB );
  }

  MatView( m_mat, viewer );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Inf-norm.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Returns the infinity norm of the matrix.
real64 PetscSparseMatrix::normInf() const
{
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
  real64 normFrob;
  MatNorm( m_mat, NORM_FROBENIUS, &normFrob );
  return normFrob;
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Is-assembled.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Boolean indicator. True = matrix assembled and ready to be used.
bool PetscSparseMatrix::isAssembled() const
{
  return m_assembled;
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// getLocalRowID
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Map a global row index to local row index
localIndex PetscSparseMatrix::getLocalRowID( globalIndex const index ) const
{
  PetscInt low, high;
  MatGetOwnershipRange( m_mat, &low, &high);
  if ( index < low || high <= index ) 
  {
    GEOSX_ERROR( "getLocalRowID: processor does not own global row index" );
  } 
  return index - low; 
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// getGlobalRowID
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Map a local row index to global row index
localIndex PetscSparseMatrix::getGlobalRowID( localIndex const index ) const
{
  PetscInt low, high;
  MatGetOwnershipRange( m_mat, &low, &high);
  if ( high - low < index ) 
  {
    GEOSX_ERROR( "getGloballRowID: processor does not own this many rows" );
  } 
  return static_cast<localIndex>( index + low );  
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// localCols
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Return the local number of columns on each processor
//
// NOTE: PETSc MPI matrices are partitioned row-wise so that the local number
// of columns is the global number.
localIndex PetscSparseMatrix::localCols( ) const
{
  PetscInt cols;
  MatGetSize( m_mat, nullptr, &cols );
  return static_cast<localIndex>( cols ); 
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// localRows
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Return the local number of rows on each processor
localIndex PetscSparseMatrix::localRows() const
{
  PetscInt low, high;
  MatGetOwnershipRange( m_mat, &low, &high);
  return high - low;
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// printParallelMatrix
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Print the given parallel matrix in Matrix Market format (MTX file)
void PetscSparseMatrix::printParallelMatrix( string const & fileName ) const
{
  PetscViewer viewer;

  // ".mtx" extension
  string name( fileName );
  if( fileName.substr( fileName.find_last_of( "." ) + 1 ) != "mtx" ) 
  {
    name = fileName.substr( 0, fileName.find_last_of( "." ) ) + ".mtx";
  }

  PetscViewerASCIIOpen( getComm(), name.c_str(), &viewer);
  PetscViewerPushFormat( viewer, PETSC_VIEWER_ASCII_MATRIXMARKET );
  MatView( m_mat, viewer );
}

// Get PETSc object
Mat PetscSparseMatrix::getConstMat() const
{
  return m_mat;
}

/* get PETSc object */
Mat PetscSparseMatrix::getMat()
{
  return m_mat;
}

localIndex PetscSparseMatrix::localNonzeros() const
{
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

} // end geosx namespace
