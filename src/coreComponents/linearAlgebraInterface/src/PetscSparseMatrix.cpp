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

// Put everything under the geosx namespace.
namespace geosx
{

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
	MatDuplicate( src._mat, MAT_COPY_VALUES, &_mat );
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
	MatCreate( comm, &_mat );
	MatSetType( _mat, MATMPIAIJ );
	MatSetSizes( _mat, localRows, localCols, PETSC_DETERMINE, PETSC_DETERMINE );
	MatMPIAIJSetPreallocation( _mat, maxEntriesPerRow, nullptr, maxEntriesPerRow, nullptr );
	MatSetUp( _mat );
}

void PetscSparseMatrix::createWithGlobalSize( globalIndex const globalRows,
                                              globalIndex const globalCols,
                                              localIndex const maxEntriesPerRow,
                                              MPI_Comm const & comm )
{
	// set up matrix
	MatCreate( comm, &_mat );
	MatSetType( _mat, MATMPIAIJ );
	MatSetSizes( _mat, PETSC_DECIDE, PETSC_DECIDE, globalRows, globalCols );
	MatMPIAIJSetPreallocation( _mat, maxEntriesPerRow, nullptr, maxEntriesPerRow, nullptr );
	MatSetUp( _mat );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Reinitialize.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Sets all values to user-defined value.
void PetscSparseMatrix::set( real64 const value )
{
	PetscInt firstrow;
	PetscInt lastrow;
	MatGetOwnershipRange( _mat, &firstrow, &lastrow );

	PetscInt numEntries;
	const PetscInt* inds;
	const PetscScalar* vals;

	PetscInt numEntries_;
	PetscScalar* vals_;
	PetscInt* inds_;

	// loop over rows
	for( PetscInt row = firstrow; row < lastrow; row++){

		// get entries in row
		close(); // hannah: this is not efficient
		MatGetRow( _mat, row, &numEntries, &inds, &vals );
		numEntries_ = numEntries;
		inds_ = new PetscInt[numEntries_];
		for ( int i = 0; i < numEntries_; i++ ) {
			inds_[i] = inds[i];
		}	
  		MatRestoreRow( _mat, row, &numEntries, &inds, &vals );
		close(); // hannah: when do I need this?

		// set entries to value
		if( numEntries_ > 0 ){

			vals_ = new PetscScalar[numEntries_];
			for ( int i = 0; i < numEntries_; i++ ) {
				vals_[i] = value;
			}

			PetscInt rows[1] = {row};
			MatSetValues( _mat, 1, rows, numEntries_, inds_, vals_, INSERT_VALUES );
		}
		// hannah: delete allocated vals_
	}
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Reinitialize.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Sets all values to 0.
void PetscSparseMatrix::zero()
{
  MatZeroEntries( _mat );
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
	MatAssemblyBegin( _mat, MAT_FINAL_ASSEMBLY );
	MatAssemblyEnd( _mat, MAT_FINAL_ASSEMBLY );
	assembled = true;
}

// -------------------------
// Add/Set
// -------------------------

// 1x1
void PetscSparseMatrix::add( globalIndex const rowIndex,
														globalIndex const colIndex,
														real64 const value )
{
	MatSetValue( _mat, rowIndex, colIndex, value, ADD_VALUES );
}

void PetscSparseMatrix::set( globalIndex const rowIndex,
														globalIndex const colIndex,
														real64 const value )
{
	MatSetValue( _mat, rowIndex, colIndex, value, INSERT_VALUES );
}

void PetscSparseMatrix::insert( globalIndex const rowIndex,
																globalIndex const colIndex,
																real64 const value )
{
	MatSetValue( _mat, rowIndex, colIndex, value, INSERT_VALUES );
}

 // 1xN c-style
void PetscSparseMatrix::add( globalIndex const rowIndex,
                             globalIndex const * colIndices,
                             real64 const * values,
                             localIndex size )
{
	PetscInt rows[1] = {rowIndex};
 	MatSetValues( _mat, 1, rows, size, colIndices, values, ADD_VALUES );
}

void PetscSparseMatrix::set( globalIndex const rowIndex,
                             globalIndex const * colIndices,
                             real64 const * values,
                             localIndex size )
{ 
	PetscInt rows[1] = {rowIndex};
 	MatSetValues( _mat, 1, rows, size, colIndices, values, INSERT_VALUES );
}

void PetscSparseMatrix::insert( globalIndex const rowIndex,
                                globalIndex const * colIndices,
                                real64 const * values,
                                localIndex size )
{
	PetscInt rows[1] = {rowIndex};
 	MatSetValues( _mat, 1, rows, size, colIndices, values, INSERT_VALUES );
}

// 1xN array1d style 
void PetscSparseMatrix::add( globalIndex const rowIndex,
                             array1d<globalIndex> const &colIndices,
                             array1d<real64> const &values )
{
	PetscInt rows[1] = {rowIndex};
 	MatSetValues( _mat, 1, rows, values.size(), colIndices.data(), values.data(), ADD_VALUES );
}

void PetscSparseMatrix::set( globalIndex const rowIndex,
                             array1d<globalIndex> const &colIndices,
                             array1d<real64> const &values )
{
	PetscInt rows[1] = {rowIndex};
 	MatSetValues( _mat, 1, rows, values.size(), colIndices.data(), values.data(), INSERT_VALUES );
}

void PetscSparseMatrix::insert( globalIndex const rowIndex,
                                array1d<globalIndex> const &colIndices,
                                array1d<real64> const &values )
{
	PetscInt rows[1] = {rowIndex};
 	MatSetValues( _mat, 1, rows, values.size(), colIndices.data(), values.data(), INSERT_VALUES );
}

// MxN array2d style
void PetscSparseMatrix::add( array1d<globalIndex> const & rowIndices,
                             array1d<globalIndex> const & colIndices,
                             array2d<real64> const & values )
{
 	MatSetValues( _mat, rowIndices.size(), rowIndices.data(), colIndices.size(), colIndices.data(), values.data(), ADD_VALUES );
}

void PetscSparseMatrix::set( array1d<globalIndex> const & rowIndices,
                             array1d<globalIndex> const & colIndices,
                             array2d<real64> const & values )
{
 	MatSetValues( _mat, rowIndices.size(), rowIndices.data(), colIndices.size(), colIndices.data(), values.data(), INSERT_VALUES );
}

void PetscSparseMatrix::insert( array1d<globalIndex> const & rowIndices,
                                array1d<globalIndex> const & colIndices,
                                array2d<real64> const & values )
{
 	MatSetValues( _mat, rowIndices.size(), rowIndices.data(), colIndices.size(), colIndices.data(), values.data(), INSERT_VALUES );
}

// -------------------------
// Linear Algebra
// -------------------------

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Matrix/vector multiplication
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Perform the matrix-vector product A*src = dst.
void PetscSparseMatrix::multiply( PetscVector const &src,
            					            PetscVector &dst ) const
{
	MatMult( _mat, src.getConstVec(), dst.getVec() );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Matrix/matrix multiplication
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Perform the matrix-matrix product A*src = dst.
void PetscSparseMatrix::multiply( PetscSparseMatrix const & src, 
    							                PetscSparseMatrix & dst ) const
{
	MatMatMult( _mat, src.getConstMat(), MAT_INITIAL_MATRIX, PETSC_DEFAULT, dst.unwrappedNonConstPointer() );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Compute residual.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Compute res = b - Ax (residual form).
void PetscSparseMatrix::residual( PetscVector const &x,
            					            PetscVector const &b,
            					            PetscVector &r ) const
{
	MatMult( _mat, x.getConstVec(), r.getConstVec() );
	VecAXPY( r.getVec(), -1, b.getConstVec() );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Generalized matrix/vector product.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Compute gemv <tt>y = alpha*A*x + beta*y</tt>.
void PetscSparseMatrix::gemv( real64 const alpha,
         					            PetscVector  const &x,
         					            real64 const beta,
         					            PetscVector  &y,
         					            bool useTranspose)
{
	PetscVector x_( x );
	PetscVector b_( x );

	x_.scale( alpha ); // alpha*x_
	y.scale( beta ); // beta*y

	if ( useTranspose ){
		MatMultTranspose( _mat, x_.getVec(), b_.getVec() );
	} else {
		MatMult( _mat, x_.getVec(), b_.getVec() ); // alpha*A*x_ = b_
	}
	VecAXPY( y.getVec(), 1, b_.getVec() ); // alpha*A*x_ + beta*y = y
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Scale.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Multiply all elements by scalingFactor.
void PetscSparseMatrix::scale( real64 const scalingFactor )
{
	MatScale( _mat, scalingFactor );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Left and right scaling
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
void PetscSparseMatrix::leftScale( PetscVector const &vec )
{
	MatDiagonalScale( _mat, vec.getConstVec(), nullptr );
}

void PetscSparseMatrix::rightScale( PetscVector const &vec )
{
	MatDiagonalScale( _mat, nullptr, vec.getConstVec() );
}

void PetscSparseMatrix::leftRightScale( PetscVector const &vecLeft,
                  					            PetscVector const &vecRight )
{
	MatDiagonalScale( _mat, vecLeft.getConstVec(), vecRight.getConstVec() );
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

  // hannah: need to check ownership?
  MatGetRow( _mat, globalRow, &numEntries, &inds, &vals );

  values.resize( numEntries );
  colIndices.resize( numEntries );

  // hannah: create array1d from c array?
  for ( int i = 0; i < numEntries; i++ ) {
	colIndices[i] = inds[i];
  }

  for ( int i = 0; i < numEntries; i++ ) {
	values[i] = vals[i];
  }

  MatRestoreRow( _mat, globalRow, &numEntries, &inds, &vals );
}

real64 PetscSparseMatrix::getDiagValue( globalIndex globalRow ) const
{
  GEOS_ERROR_IF( !assembled, "Attempting to call " << __FUNCTION__ << " before close() is illegal" );

  const PetscScalar *vals = nullptr;
  const PetscInt *cols = nullptr;
  PetscInt ncols;

  MatGetRow( _mat, globalRow, &ncols, &cols, &vals );
  for( int i = 0; i < ncols; i++ ){
    if( cols[i] == globalRow )
    {
      return vals[i]; // hannah (*vals)[i]?
    }
  }
  MatRestoreRow( _mat, globalRow, &ncols, &cols, &vals );

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
	MatZeroRows( _mat, 1, rows, diagValue, nullptr, nullptr );
 }

// ----------------------------
//  Accessors
// ----------------------------

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get pointer.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""

// Accessor for the pointer to the matrix
Mat* PetscSparseMatrix::unwrappedNonConstPointer()
{
	return &( _mat );
}

// Accessor for the pointer to the matrix
const Mat* PetscSparseMatrix::unwrappedPointer() const
{
	return &( _mat );
}

// Accessor for the MPI communicator
MPI_Comm PetscSparseMatrix::getComm() const
{
	MPI_Comm comm;
	PetscObjectGetComm( reinterpret_cast<PetscObject>( _mat ), &comm );
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
	// hannah: type conversion
	MatGetSize( _mat, &num_rows, &num_cols );
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
	MatGetSize( _mat, &num_rows, &num_cols );
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
	MatGetOwnershipRange( _mat, &firstrow, &lastrow );
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
 	MatGetOwnershipRange( _mat, &firstrow, &lastrow );
 	return lastrow; 
 } 

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Print to terminal.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Wrapper to print the petsc output of the matrix
void PetscSparseMatrix::print() const 
{
	MatView( _mat, PETSC_VIEWER_STDOUT_WORLD );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Write to matlab-compatible file
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
void PetscSparseMatrix::write( string const & filename,
                               bool const mtxFormat ) const
{
	PetscViewer viewer;
	const char * filename_ = filename.c_str();

	
	if( mtxFormat ){

		// ".mtx" extension
    	string name( filename );
		if( filename.substr( filename.find_last_of( "." ) + 1 ) != "mtx" ){
			name = filename.substr( 0, filename.find_last_of( "." ) ) + ".mtx";
		}
		PetscViewerASCIIOpen( getComm(), name.c_str(), &viewer);
		PetscViewerPushFormat( viewer, PETSC_VIEWER_ASCII_MATRIXMARKET );

	} else {
		PetscViewerASCIIOpen( getComm(), filename.c_str(), &viewer);
		PetscViewerPushFormat( viewer, PETSC_VIEWER_ASCII_MATLAB );
	}
	MatView( _mat, viewer );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Inf-norm.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Returns the infinity norm of the matrix.
real64 PetscSparseMatrix::normInf() const
{
	real64 normInf;
	MatNorm( _mat, NORM_INFINITY, &normInf );
	return normInf;
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// 1-norm.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Returns the one norm of the matrix.
real64 PetscSparseMatrix::norm1() const
{
	real64 norm1;
	MatNorm( _mat, NORM_1, &norm1 );
	return norm1;
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Frobenius-norm.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Returns the Frobenius norm of the matrix.
real64 PetscSparseMatrix::normFrobenius() const
{
	real64 normFrob;
	MatNorm( _mat, NORM_FROBENIUS, &normFrob );
	return normFrob;
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Is-assembled.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Boolean indicator. True = matrix assembled and ready to be used.
bool PetscSparseMatrix::isAssembled() const
{
	return assembled;
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// MatrixMatrixMultiply
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Perform the matrix-matrix product A*src = dst.
void PetscSparseMatrix::MatrixMatrixMultiply( bool const transA,
                                         PetscSparseMatrix const &B,
                                         bool const transB,
                                         PetscSparseMatrix &C,
                                         bool const call_FillComplete ) const
{
	// hannah: fill complete?
	if( transA && transB ) {
		MatMatMult( _mat, B.getConstMat(), MAT_INITIAL_MATRIX, PETSC_DEFAULT, C.unwrappedNonConstPointer() );
		MatTranspose( C.getConstMat(), MAT_INPLACE_MATRIX, C.unwrappedNonConstPointer() );
	} else if (transA ) {
		MatTransposeMatMult( _mat, B.getConstMat(), MAT_INITIAL_MATRIX, PETSC_DEFAULT, C.unwrappedNonConstPointer() );	
	} else if (transB ) {
		Mat BT; // hannah: better way to do this?
		MatMatTransposeMult( _mat, B.getConstMat(), MAT_INITIAL_MATRIX, PETSC_DEFAULT, C.unwrappedNonConstPointer() ); // not for mpiaij matrices
	} else {
		MatMatMult( _mat, B.getConstMat(), MAT_INITIAL_MATRIX, PETSC_DEFAULT, C.unwrappedNonConstPointer() );	
	}
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// getLocalRowID
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Map a global row index to local row index
localIndex PetscSparseMatrix::getLocalRowID( globalIndex const index ) const
{
	PetscInt low, high;
	MatGetOwnershipRange( _mat, &low, &high);
	if ( index < low || high <= index ) {
		GEOS_ERROR( "getLocalRowID: processor does not own global row index" );	// hannah: error here?
	} 
	return index - low;	
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// getGlobalRowID
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Map a local row index to global row index
localIndex PetscSparseMatrix::getGlobalRowID( localIndex const index ) const // hannah: should this return globalIndex?
{
	PetscInt low, high;
	MatGetOwnershipRange( _mat, &low, &high);
	if ( high - low < index ) {
		GEOS_ERROR( "getGloballRowID: processor does not own this many rows" );	// hannah: error here?
	} 
	return static_cast<localIndex>( index + low );	
	// return index + low;	
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// numMyCols
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Return the local number of columns on each processor
localIndex PetscSparseMatrix::numMyCols( ) const
{
	// total number of columns
	// hannah: right now PETSc is partitioned row-wise, would be different if block partitioned
	PetscInt cols;
	MatGetSize( _mat, nullptr, &cols );
	return static_cast<localIndex>( cols );	
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
	if( fileName.substr( fileName.find_last_of( "." ) + 1 ) != "mtx" ){
		name = fileName.substr( 0, fileName.find_last_of( "." ) ) + ".mtx";
	}

	PetscViewerASCIIOpen( getComm(), name.c_str(), &viewer);
	PetscViewerPushFormat( viewer, PETSC_VIEWER_ASCII_MATRIXMARKET );
	MatView( _mat, viewer );
}

// Get PETSc object
Mat PetscSparseMatrix::getConstMat() const
{
	return _mat;
}

/* get PETSc object */
Mat PetscSparseMatrix::getMat()
{
	return _mat;
}

} // end geosx namespace
