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
	MatDuplicate(in_matrix._mat, MAT_COPY_VALUES, &_mat);
}

virtual ~PetscSparseMatrix()
{
	if(_mat) MatDestroy(_mat);
}

// -----------------------------
// Create/Finalize
// -----------------------------
// Allocate matrix (prepare to be filled with data).

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Create a matrix from number of elements
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
void PetscSparseMatrix::create( MPI_Comm const comm,
								int const m_nRowGlobal,
								int const nMaxEntriesPerRow )
{
	// set up matrix
	MatCreate(comm, &_mat);
	MatSetType(_mat, MATMPIAIJ);
	MatSetSizes(_mat, PETSC_DECIDE, PETSC_DECIDE, m_nRowGlobal, m_nRowGlobal);
	MatMPIAIJSetPreallocation(_mat, nMaxEntriesPerRow, NULL, nMaxEntriesPerRow, NULL);
	MatSetOption(_mat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
	MatSetUp(_mat);
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Create a rectangular matrix from number of elements
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
 void PetscSparseMatrix::create( MPI_Comm const comm,
              					 int const m_nRowGlobal,
              					 int const m_nColGlobal,
              					 int const nMaxEntriesPerRow )
 {
	// set up matrix
	MatCreate(comm, &_mat);
	MatSetType(_mat, MATMPIAIJ);
	MatSetSizes(_mat, PETSC_DECIDE, PETSC_DECIDE, m_nRowGlobal, m_nColGlobal);
	MatMPIAIJSetPreallocation(_mat, nMaxEntriesPerRow, NULL, nMaxEntriesPerRow, NULL);
	MatSetOption(_mat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
	MatSetUp(_mat);
}

  // """""""""""""""""""""""""""""""""""""""""""""""""""""""""
  // Create a square matrix from number of elements in each row
  // """""""""""""""""""""""""""""""""""""""""""""""""""""""""
  void PetscSparseMatrix::create( MPI_Comm const comm,
               					  int const m_nRowGlobal,
               					  std::vector<int> const nMaxEntriesPerRow )
  {
  	// set up matrix
	MatCreate(comm, &_mat);
	MatSetType(_mat, MATMPIAIJ);
	MatSetSizes(_mat, PETSC_DECIDE, PETSC_DECIDE, m_nRowGlobal, m_nRowGlobal);
	MatMPIAIJSetPreallocation(_mat, 0, nMaxEntriesPerRow.data(), 0, nMaxEntriesPerRow.data());
	MatSetOption(_mat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
	MatSetUp(_mat);
  }

  // """""""""""""""""""""""""""""""""""""""""""""""""""""""""
  // Create a rectangular matrix from number of elements in each row
  // """""""""""""""""""""""""""""""""""""""""""""""""""""""""
  void PetscSparseMatrix::create( MPI_Comm const comm,
               int const m_nRowGlobal,
               int const m_nColGlobal,
               std::vector<int> const nMaxEntriesPerRow )
  {
  	// set up matrix
	MatCreate(comm, &_mat);
	MatSetType(_mat, MATMPIAIJ);
	MatSetSizes(_mat, PETSC_DECIDE, PETSC_DECIDE, m_nRowGlobal, m_nColGlobal);
	MatMPIAIJSetPreallocation(_mat, 0, nMaxEntriesPerRow.data(), 0, nMaxEntriesPerRow.data());
	MatSetOption(_mat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
	MatSetUp(_mat);
  }

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Create a matrix from an PetscSparseMatrix.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
 void PetscSparseMatrix::create( PetscSparseMatrix &matrix )
 {
 	MatDuplicate(matrix._mat, MAT_COPY_VALUES, &_mat);
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
 	// MatAssemblyBegin(_mat, MAT_FINAL_ASSEMBLY);
 }

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Close
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Fix the sparsity pattern, make the data contiguous in memory and optimize storage.
 void PetscSparseMatrix::close()
 {
 	MatAssemblyBegin(_mat, MAT_FINAL_ASSEMBLY);
 	MatAssemblyEnd(_mat, MAT_FINAL_ASSEMBLY);
 	assembled = true;
 }

 // -------------------------
// Add/Set
// -------------------------

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Add row values.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Add row values at row iRow and columns cols (size nCols)
 void PetscSparseMatrix::add( int const iRow,
           					  int const nCols,
           					  real64 const *values,
           					  int const *cols )
 {
 	int rows[1] = {iRow};
 	MatSetValues(_mat, 1, rows, nCols, cols, values, ADD_VALUES);
 }

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Add single value.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Add a value at row iRow and column iCol
 void PetscSparseMatrix::add( int const iRow,
           					  int const iCol,
           					  real64 const value )
 {
 	MatSetValue(_mat, iRow, iCol, value, ADD_VALUES);
 }

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Set row values.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Add the value of elements (iRow,[nCols]) to values.
 void PetscSparseMatrix::set( int const iRow,
           					  int const nCols,
           					  real64 const *values,
           					  int const *cols )
  {
 	int rows[1] = {iRow};
 	MatSetValues(_mat, 1, rows, nCols, cols, values, INSERT_VALUES);
 }

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Set single value.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Set the value of the element (iRow,iCol) to value.
 void PetscSparseMatrix::set( int const iRow,
           					  int const iCol,
           					  real64 const value )
  {
 	MatSetValue(_mat, iRow, iCol, value, INSERT_VALUES);
 }

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Insert values.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
 void PetscSparseMatrix::insert( int const iRow,
              					 int const nCols,
              					 real64 const *values,
              					 int const *cols )
  {
 	int rows[1] = {iRow};
 	MatSetValues(_mat, 1, rows, nCols, cols, values, INSERT_VALUES);
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
	// need to make sure dst is not NULL
	MatMult(_mat, src.getConstVec(), dst.getVec());
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
// Left scale (diagonal scaling).
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Pre-multiplies (left) with diagonal matrix consisting of the values in vec.
void PetscSparseMatrix::leftScale( PetscVector const &vec )
{
	MatDiagonalScale(_mat, vec.getConstVec(), NULL);
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Right scale (diagonal scaling)
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Post-multiplies (right) with diagonal matrix consisting of the values in vec.
void PetscSparseMatrix::rightScale( PetscVector const &vec )
{
	MatDiagonalScale(_mat, NULL, vec.getConstVec());
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Left and Right scale (diagonal scalings)
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Pre-multiplies (left) with diagonal matrix consisting of the values in vecLeft and
// Post-multiplies (right) with diagonal matrix consisting of the values in vecRight.
void PetscSparseMatrix::leftRightScale( PetscVector const &vecLeft,
                  					 PetscVector const &vecRight )
{
	MatDiagonalScale(_mat, vecLeft.getConstVec(), vecRight.getConstVec());
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Clear row.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Clear row and multiply diagonal term by factor.
void PetscSparseMatrix::clearRow( int const row,
            				   	  real64 const factor )
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
// Get global row.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Extract the global row and output the number of elements, values and column indices.
void PetscSparseMatrix::getRow( int GlobalRow,
           						int &NumEntries,
           						real64* Values,
           						int* Indices ) const
{
 	const double* vals;
 	const int* inds;
 	int numEntries;
	MatGetRow(_mat, GlobalRow, &numEntries, &inds, &vals);

    for(int i = 0; i < numEntries; i++)
    {
    	Values[i] = vals[i];
    	Indices[i] = inds[i];
    }
    NumEntries = numEntries;
    
	MatRestoreRow(_mat, GlobalRow, &numEntries, &inds, &vals);
 }

/*
* Get global row myRow
* - numEntries: number of nonzeros 
* - vecValues: vector of values
* - vecIndices: vector of column indices */
void PetscSparseMatrix::getRow( int GlobalRow,
             					int &NumEntries,
             					std::vector<real64> &vecValues,
             					std::vector<int> &vecIndices ) const
{
	const double* vals;
	const int* inds;
	int numEntries;
	MatGetRow(_mat, GlobalRow, &numEntries, &inds, &vals);

	vecIndices.assign(inds, inds + numEntries);
	vecValues.assign(vals, vals + numEntries);
	NumEntries = numEntries;

	MatRestoreRow(_mat, GlobalRow, &numEntries, &inds, &vals);
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get local row.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Extract the local row and output the number of elements, values and column indices.
void PetscSparseMatrix::getLocalRow( int myRow,
               					 	 int & NumEntries,
               					 	 real64 * & Values,
               					 	 int * & Indices ) const
{
 	// myRow -> globalRow
 	int firstRow, lastRow;
 	int GlobalRow;
 	MatGetOwnershipRange(_mat, &firstRow, &lastRow);
 	GlobalRow = firstRow + myRow;

 	const double* vals;
 	const int* inds;
 	int numEntries;
	MatGetRow(_mat, GlobalRow, &numEntries, &inds, &vals);

    for(int i = 0; i < numEntries; i++)
    {
    	Values[i] = vals[i];
    	Indices[i] = inds[i];
    }
    NumEntries = numEntries;
    
	MatRestoreRow(_mat, GlobalRow, &numEntries, &inds, &vals);
 }

 /*
  * Get local row myRow
  * - numEntries: number of nonzeros 
  * - vecValues: vector of values
  * - vecIndices: vector of column indices */
 void PetscSparseMatrix::getLocalRow( int myRow,
                   					int &NumEntries,
                   					std::vector<real64> &vecValues,
                   					std::vector<int> &vecIndices ) const

 {
 	// myRow -> globalRow
 	int firstRow, lastRow;
 	int GlobalRow;
 	MatGetOwnershipRange(_mat, &firstRow, &lastRow);
 	GlobalRow = firstRow + myRow;

 	const double* vals;
 	const int* inds;
 	int numEntries;
	MatGetRow(_mat, GlobalRow, &numEntries, &inds, &vals);

	vecIndices.assign(inds, inds + numEntries);
    vecValues.assign(vals, vals + numEntries);
 	NumEntries = numEntries;

	MatRestoreRow(_mat, GlobalRow, &numEntries, &inds, &vals);
 }

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get number of global rows.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Accessor for the number of global rows
int PetscSparseMatrix::globalRows() const
{
	int num_rows;
	int num_cols;
	MatGetSize(_mat, &num_rows, &num_cols);
	return num_rows;
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get number of global columns.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Accessor for the number of global columns
int PetscSparseMatrix::globalCols() const
{
	int num_rows;
	int num_cols;
	MatGetSize(_mat, &num_rows, &num_cols);
	return num_cols;
}

/* Number of unique columns */
// int PetscSparseMatrix::uniqueCols() const;

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get the lower index owned by processor.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Accessor for the index of the first global row
int PetscSparseMatrix::ilower() const
{
	int firstrow;
	int lastrow;
	MatGetOwnershipRange(_mat, &firstrow, &lastrow);
	return firstrow; 
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get the upper index owned by processor.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Accessor for the index of the last global row
int PetscSparseMatrix::iupper() const
{
 	int firstrow;
 	int lastrow;
 	MatGetOwnershipRange(_mat, &firstrow, &lastrow);
 	return lastrow; 
 } 

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Test for ownership of global index.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Returns the local index if owned by this processor. Returns -1 if
// not owned.
int PetscSparseMatrix::rowMapLID( int const GID ) const
{
	int firstrow;
 	int lastrow;
 	MatGetOwnershipRange(_mat, &firstrow, &lastrow);

 	if (firstrow <= GID && GID < lastrow){
 		return GID - firstrow;
 	} else {
 		return -1;
 	}
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get the number of local rows.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Accessor for the number of local rows
int PetscSparseMatrix::myRows() const
{
	int firstrow;
	int lastrow;
	MatGetOwnershipRange(_mat, &firstrow, &lastrow);
	return lastrow - firstrow; 
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get the number of local columns.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Accessor for the number of local columns
int PetscSparseMatrix::myCols() const
{
	int firstcol;
	int lastcol;
	MatGetOwnershipRangeColumn(_mat, &firstcol, &lastcol);
	return lastcol - firstcol; 
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

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Print to terminal.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Wrapper to print the petsc output of the matrix
void PetscSparseMatrix::print() const 
{
	MatView(_mat, PETSC_VIEWER_STDOUT_WORLD);
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get pointer.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Accessor for the pointer to the matrix
const Mat* PetscSparseMatrix::getPointer() const
{
	return &(_mat);
}

/* get PETSc object */
Mat PetscSparseMatrix::getMat()
{
	return _mat;
}

