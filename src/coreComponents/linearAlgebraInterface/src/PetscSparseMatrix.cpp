/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
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
{
	// do nothing
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Copy constructor
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
PetscSparseMatrix::PetscSparseMatrix( PetscSparseMatrix const &in_matrix )
{
	// Hannah: check that in_matrix is not empty
	MatDuplicate(in_matrix._mat, MAT_COPY_VALUES, &_mat);
}

// -----------------------------
// Create
// -----------------------------
// Allocate matrix (prepare to be filled with data).

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Create a matrix from an Epetra_FECrsGraph.
// """""""""""""""""""""""""""""""""""""""""""""""
void PetscSparseMatrix::create( Epetra_FECrsGraph const &graph )
{
	// Hannah: skip?
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Create a matrix from number of elements
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
void PetscSparseMatrix::createWithLocalSize( localIndex const localSize,
                                             localIndex const maxEntriesPerRow,
                                             MPI_Comm const & comm )
{
		// set up matrix
		MatCreate(comm, &_mat);
		MatSetType(_mat, MATMPIAIJ);
		MatSetSizes(_mat, localSize, localSize, PETSC_DETERMINE, PETSC_DETERMINE);
		// Hannah: being conservative here with maxEntriesPerRow
		MatMPIAIJSetPreallocation(_mat, maxEntriesPerRow, NULL, maxEntriesPerRow, NULL);
		MatSetOption(_mat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
		MatSetUp(_mat);

		// assemble
		// Hannah: assemble here?
		MatAssemblyBegin(_mat, MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(_mat, MAT_FINAL_ASSEMBLY);

		assembled = true;

		// Hannah: something wrong with number of columns

	}

  void PetscSparseMatrix::createWithGlobalSize(globalIndex const globalSize,
                                               localIndex const maxEntriesPerRow,
                                               MPI_Comm const & comm )
{
		// set up matrix
		MatCreate(comm, &_mat);
		MatSetType(_mat, MATMPIAIJ);
		MatSetSizes(_mat, PETSC_DECIDE, PETSC_DECIDE, globalSize, globalSize);
		MatMPIAIJSetPreallocation(_mat, maxEntriesPerRow, NULL, maxEntriesPerRow, NULL);
		MatSetOption(_mat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
		MatSetUp(_mat);

		// assemble
		MatAssemblyBegin(_mat, MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(_mat, MAT_FINAL_ASSEMBLY);

		assembled = true;
	}

  void PetscSparseMatrix::createWithLocalSize(localIndex const localRows,
                                              localIndex const localCols,
                                              localIndex const maxEntriesPerRow,
                                              MPI_Comm const & comm )
{
		// set up matrix
		MatCreate(comm, &_mat);
		MatSetType(_mat, MATMPIAIJ);
		MatSetSizes(_mat, localRows, localCols, PETSC_DETERMINE, PETSC_DETERMINE);
		MatMPIAIJSetPreallocation(_mat, maxEntriesPerRow, NULL, maxEntriesPerRow, NULL);
		MatSetOption(_mat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
		MatSetUp(_mat);

		// assemble
		MatAssemblyBegin(_mat, MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(_mat, MAT_FINAL_ASSEMBLY);

		assembled = true;
	}

  void PetscSparseMatrix::createWithGlobalSize(globalIndex const globalRows,
                                               globalIndex const globalCols,
                                               localIndex const maxEntriesPerRow,
                                               MPI_Comm const & comm )
{
		// set up matrix
		MatCreate(comm, &_mat);
		MatSetType(_mat, MATMPIAIJ);
		MatSetSizes(_mat, PETSC_DECIDE, PETSC_DECIDE, globalRows, globalCols);
		MatMPIAIJSetPreallocation(_mat, maxEntriesPerRow, NULL, maxEntriesPerRow, NULL);
		MatSetOption(_mat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
		MatSetUp(_mat);

		// assemble
		MatAssemblyBegin(_mat, MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(_mat, MAT_FINAL_ASSEMBLY);

		assembled = true;
	}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Reinitialize.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Keeps the map and graph but sets all values to 0.
 void PetscSparseMatrix::zero()
 {
 	MatZeroEntries(_mat);
 }

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Open
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Empty open function (implemented for HYPRE compatibility).
 void PetscSparseMatrix::open()
 {
 	// do nothing
 }

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Close
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Hannah: blurb here
 void PetscSparseMatrix::close()
 {
 	MatAssemblyBegin(_mat, MAT_FINAL_ASSEMBLY);
 	MatAssemblyEnd(_mat, MAT_FINAL_ASSEMBLY);
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
 	MatSetValue(_mat, rowIndex, colIndex, value, ADD_VALUES);
 }

 void PetscSparseMatrix::set( globalIndex const rowIndex,
                              globalIndex const colIndex,
                              real64 const value )
 {
 	MatSetValue(_mat, rowIndex, colIndex, value, INSERT_VALUES);
 }

 void PetscSparseMatrix::insert( globalIndex const rowIndex,
                              globalIndex const colIndex,
                              real64 const value )
 {
 	MatSetValue(_mat, rowIndex, colIndex, value, INSERT_VALUES);
 }

 // 1xN with arrays
void PetscSparseMatrix::add( globalIndex const rowIndex,
                       globalIndex const * colIndices,
                       real64 const * values,
                       localIndex size )
{
	int rows[1] = {rowIndex};
 	MatSetValues(_mat, 1, rows, size, colIndices, values, ADD_VALUES);
}

void PetscSparseMatrix::set( globalIndex const rowIndex,
                       globalIndex const * colIndices,
                       real64 const * values,
                       localIndex size )
{
	int rows[1] = {rowIndex};
 	MatSetValues(_mat, 1, rows, size, colIndices, values, INSERT_VALUES);
}

void PetscSparseMatrix::insert( globalIndex const rowIndex,
                          globalIndex const * colIndices,
                          real64 const * values,
                          localIndex size )
{
	int rows[1] = {rowIndex};
 	MatSetValues(_mat, 1, rows, size, colIndices, values, INSERT_VALUES);
	// Hannah: insert vs set?
}

// 1xN with array1d 
void PetscSparseMatrix::add( globalIndex const rowIndex,
                        array1d<globalIndex> const &colIndices,
                        array1d<real64> const &values )
{
	// Hannah: to do
}

void PetscSparseMatrix::set( globalIndex const rowIndex,
                        array1d<globalIndex> const &colIndices,
                        array1d<real64> const &values )
{
	// Hannah: to do
}

void PetscSparseMatrix::insert( globalIndex const rowIndex,
                           array1d<globalIndex> const &colIndices,
                           array1d<real64> const &values )
{
	// Hannah: to do
}

// MxN with array2d 
void PetscSparseMatrix::add( array1d<globalIndex> const & rowIndices,
                        array1d<globalIndex> const & colIndices,
                        array2d<real64> const & values )
{
	// Hannah: to do
}

void PetscSparseMatrix::set( array1d<globalIndex> const & rowIndices,
                        array1d<globalIndex> const & colIndices,
                        array2d<real64> const & values )
{
	// Hannah: to do
}

void PetscSparseMatrix::insert( array1d<globalIndex> const & rowIndices,
                           array1d<globalIndex> const & colIndices,
                           array2d<real64> const & values )
{
	// Hannah: to do
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
	// Hannah: check dst is not NULL
	MatMult(_mat, src.getConstVec(), dst.getVec());
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Matrix/matrix multiplication
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Perform the matrix-matrix product A*src = dst.
void PetscSparseMatrix::multiply(PetscSparseMatrix const & src, 
    														 PetscSparseMatrix & dst) const
{
	MatMatMult(_mat, src.getConstMat(), MAT_INITIAL_MATRIX, PETSC_DEFAULT, dst.getPointer_());
	// Hannah: getPointer_ name?
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Compute residual.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Compute res = b - Ax (residual form).
void PetscSparseMatrix::residual( PetscVector const &x,
            					  PetscVector const &b,
            					  PetscVector &res ) const
{
	MatMult(_mat, x.getConstVec(), res.getConstVec());
	VecAXPY(res.getVec(), -1, b.getConstVec());
	// Hannah: everything should be a pointer?
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
	PetscVector x_(x);
	PetscVector b_(x);

	x_.scale(alpha); // alpha*x_
	y.scale(beta); // beta*y

	if (useTranspose){
		MatMultTranspose(_mat, x_.getVec(), b_.getVec());
	} else {
		MatMult(_mat, x_.getVec(), b_.getVec()); // alpha*A*x_ = b_
	}
	VecAXPY(y.getVec(), 1, b_.getVec()); // alpha*A*x_ + beta*y = y
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Scale.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Multiply all elements by scalingFactor.
void PetscSparseMatrix::scale( real64 const scalingFactor )
{
	MatScale(_mat, scalingFactor);
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Left and right scaling
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
void PetscSparseMatrix::leftScale( PetscVector const &vec )
{
	MatDiagonalScale(_mat, vec.getConstVec(), NULL);
}

void PetscSparseMatrix::rightScale( PetscVector const &vec )
{
	MatDiagonalScale(_mat, NULL, vec.getConstVec());
}

void PetscSparseMatrix::leftRightScale( PetscVector const &vecLeft,
                  					 PetscVector const &vecRight )
{
	MatDiagonalScale(_mat, vecLeft.getConstVec(), vecRight.getConstVec());
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
  // Hannah: to do
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Clear row.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Clear the row.  By default the diagonal value will be set
// to zero, but the user can pass a desired diagValue
void PetscSparseMatrix::clearRow( int const row,
            				   	  PetscScalar const factor )
 {
	int rows[1] = {row};
	double diag = 1;

	// who does this row belong to
	int firstRow, lastRow;
	MatGetOwnershipRange(_mat, &firstRow, &lastRow);

	if(firstRow <= row && row < lastRow)
	{
		// get diagonal entry
		const double* vals;
		const int* inds;
		int numEntries;

		MatGetRow(_mat, row, &numEntries, &inds, &vals);

		for(int i = 0; i < numEntries; i++)
		{
			if(inds[i] == row) diag = vals[i];
		}
		MatRestoreRow(_mat, row, &numEntries, &inds, &vals);
	}

	// zero row and multiply diagonal by factor
	MatZeroRows(_mat, 1, rows, diag*factor, NULL, NULL);
 }

// ----------------------------
//  Accessors
// ----------------------------

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get pointer.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Accessor for the pointer to the matrix
const Mat* PetscSparseMatrix::unwrappedPointer() const
{
	return &(_mat);
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get number of global rows.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Accessor for the number of global rows
globalIndex PetscSparseMatrix::globalRows() const
{
	globalIndex num_rows;
	globalIndex num_cols;
	// Hannah: type conversion
	MatGetSize(_mat, &num_rows, &num_cols);
	return num_rows;
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get number of global columns.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Accessor for the number of global columns
globalIndex PetscSparseMatrix::globalCols() const
{
	globalIndex num_rows;
	globalIndex num_cols;
	// Hannah: type conversion
	MatGetSize(_mat, &num_rows, &num_cols);
	return num_cols;
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get the lower index owned by processor.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Accessor for the index of the first global row
globalIndex PetscSparseMatrix::ilower() const
{
	globalIndex firstrow;
	globalIndex lastrow;
	// Hannah: type conversion
	MatGetOwnershipRange(_mat, &firstrow, &lastrow);
	return firstrow; 
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get the upper index owned by processor.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Accessor for the index of the last global row
globalIndex PetscSparseMatrix::iupper() const
{
 	globalIndex firstrow;
 	globalIndex lastrow;
 	MatGetOwnershipRange(_mat, &firstrow, &lastrow);
	 // Hannah: off by one?
 	return lastrow; 
 } 

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Print to terminal.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Wrapper to print the petsc output of the matrix
void PetscSparseMatrix::print() const 
{
	MatView(_mat, PETSC_VIEWER_STDOUT_WORLD);
}


// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Write to matlab-compatible file
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Note: EpetraExt also supports a MatrixMarket format as well
void PetscSparseMatrix::write( string const & filename ) const
{
  //  strcpy(filename_char, filename.c_str()); 

	// 	// set up PETSc viewer
	// 	PetscViewer viewer;
	// 	PetscViewerCreate(PETSC_COMM_WORLD, &viewer);
	// 	PetscViewerSetType(viewer, PETSCVIEWERBINARY);
	// 	PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
	// 	PetscViewerFileSetMode(viewer, FILE_MODE_WRITE);
	// 	PetscViewerFileSetName(viewer, filename_char);
	// 	MatView(_mat, viewer);

	// Hannah: to do
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Inf-norm.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Returns the infinity norm of the matrix.
real64 PetscSparseMatrix::normInf() const
{
	real64 normInf;
	MatNorm(_mat, NORM_INFINITY, &normInf);
	return normInf;
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// 1-norm.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Returns the one norm of the matrix.
real64 PetscSparseMatrix::norm1() const
{
	real64 norm1;
	MatNorm(_mat, NORM_1, &norm1);
	return norm1;
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Frobenius-norm.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Returns the Frobenius norm of the matrix.
real64 PetscSparseMatrix::normFrobenius() const
{
	real64 normFrob;
	MatNorm(_mat, NORM_FROBENIUS, &normFrob);
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

/* get PETSc object */
Mat PetscSparseMatrix::getMat()
{
	return _mat;
}

} // end geosx namespace