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
 * @file PetscVector.cpp
 *
 *  Created on: Jan 28, 2019
 *  Author: Hannah Morgan
 */

// BEGIN_RST_NARRATIVE PetscSparseVector.rst
// ==============================
// Petsc-based Vector Object
// ==============================
// This class contains the ParallelVector wrappers based on Petsc Vec Objects.
// The class contains a unique pointer to a Vec as well as constructors,
// functions and accessors for PETSc objects.

// Include the corresponding header file.
#include "PetscVector.hpp"

// Put everything under the geosx namespace.
namespace geosx
{

// ----------------------------
// Constructors
// ----------------------------

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Empty constructor
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Create an empty vector
PetscVector::PetscVector()
{
	// Do nothing
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Copy constructor
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Create a unique PetscVector from PetscVector.  
// The data from the input vector is copied to a new memory location. 
PetscVector::PetscVector(PetscVector const & vec)
{
	// Hannah: check if vec is empty
	VecDuplicate(vec._vec, &_vec); 
	VecCopy(vec._vec, _vec); 
}

// Create a unique PetscVector from a PETSc Vec.  
PetscVector::PetscVector(Vec vec)
{
	_vec = vec;
}

// ----------------------------
// Create
// ----------------------------

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Create with known size
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// There are two variants of this function.  In the first, the user knows
// the local size and wants the global size to be the sum of each processor's
// contributions.  In the second, the user knows the global size and wants a
// near-even distribution of elements across processors. All processors
// get the same number of elements, except proc 0 which gets any remainder
// elements necessary when the number of processors does not divide evenly
// into the vector length.
void PetscVector::createWithLocalSize( localIndex const localSize, MPI_Comm const & comm )
{
  VecCreateMPI(comm, localSize, PETSC_DETERMINE, &_vec);
}

void PetscVector::createWithGlobalSize( globalIndex const globalSize, MPI_Comm const & comm )
{
  VecCreateMPI(comm, PETSC_DECIDE, globalSize, &_vec);
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Create from array
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Create a vector from local array data.  The global vector contains
// local arrays stitched together.
void PetscVector::create( array1d<real64> const & localValues, MPI_Comm const & comm )
{
	// Hannah: to do
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Add/Set value(s)
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Add/set entries in the vector.  
// Hannah: might need type conversions

// single element 
void PetscVector::set( globalIndex const globalRow,
                        real64 const value 
{
	VecSetValue(_vec, globalRow, value, INSERT_VALUES);
	VecAssemblyBegin(_vec);
  	VecAssemblyEnd(_vec);
}

void PetscVector::add( globalIndex const globalRow,
                        real64 const value 
{
	VecSetValue(_vec, globalRow, value, ADD_VALUES);
	VecAssemblyBegin(_vec);
  	VecAssemblyEnd(_vec);
}

// multiple elements, arrays
void PetscVector::set( globalIndex const * globalIndices,
                        real64 const * values,
                        localIndex size )
{
  VecSetValues(_vec, size, globalIndices, values, INSERT_VALUES);
  VecAssemblyBegin(_vec);
  VecAssemblyEnd(_vec);
}

void PetscVector::add( globalIndex const * globalIndices,
                        real64 const * values,
                        localIndex size )
{
  VecSetValues(_vec, size, globalIndices, values, ADD_VALUES);
  VecAssemblyBegin(_vec);
  VecAssemblyEnd(_vec);
}

// multiple elements, array1d
void PetscVector::set( array1d<globalIndex> const & globalIndices,
                        array1d<real64> const & values )
{
	// Hannah: to do
}
void PetscVector::add( array1d<globalIndex> const & globalIndices,
                        array1d<real64> const & values )
{
	// Hannah: to do
}

// set all elements
void PetscVector::set(real64 value)
{
	VecSet(_vec, value);
	VecAssemblyBegin(_vec);
  	VecAssemblyEnd(_vec);
}

// zero all entries
void PetscVector::zero()
{
	VecZeroEntries(_vec);
}

// fill with random numbers
void PetscVector::rand()
{
	// Hannah: how random do we need this to be?
	PetscRandom ran;
	PetscRandomCreate(PETSC_COMM_WORLD, &ran);

	time_t seconds;
  	seconds = time (NULL);
	PetscRandomSetSeed(ran, seconds);
	PetscRandomSeed(ran);

	VecSetRandom(_vec, ran);
	PetscRandomDestroy(&ran);
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Open / close
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
void PetscVector::open()
{
	// do nothing
}

void PetscVector::close()
{
	VecAssemblyBegin(_vec);
  	VecAssemblyEnd(_vec);	
}

// ----------------------------
// Linear Algebra
// ----------------------------
// The following functions support basic linear algebra ops

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Scale
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Multiply all elements by scalingFactor.
void PetscVector::scale(real64 const scalingFactor)
{
	VecScale(_vec, scalingFactor);
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Dot
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Dot product with the vector vec.
real64 PetscVector::dot(PetscVector const &vec)
{
	real64 dot;
	VecDot(_vec, vec._vec, &dot);
	return dot;
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Copy
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Update vector as this = x.
void PetscVector::copy(PetscVector const &x)
{
	VecSet(_vec, 0);
	VecAXPY(_vec, 1.0, x._vec);
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Axpy
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Update vector as this = alpha*x + this.
void PetscVector::axpy(real64 const alpha, PetscVector const &x)
{
	VecAXPY(_vec, alpha, x._vec);
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Axpby
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Update vector as this = alpha*x + beta*this.
void PetscVector::axpby(real64 const alpha, PetscVector &x, real64 const beta)
{
	VecScale(_vec, beta);
	VecAXPY(_vec, alpha, x._vec);
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// 1-norm
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// 1-norm of the vector.
real64 PetscVector::norm1() const
{
	real64 result;
	VecNorm(_vec, NORM_1, &result);
	return result;
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// 2-norm
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// 2-norm of the vector.
real64 PetscVector::norm2() const
{
	real64 result;
	VecNorm(_vec, NORM_2, &result);
	return result;
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Inf-norm
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Inf-norm of the vector.
real64 PetscVector::normInf() const
{
	real64 result;
	VecNorm(_vec, NORM_INFINITY, &result);
	return result;
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Print
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Print vector to the terminal in PETSc format.
void PetscVector::print() const
{
	VecView(_vec, PETSC_VIEWER_STDOUT_WORLD);
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Write to matlab-compatible file
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
void PetscVector::write( string const & filename ) const
{
	// // need a char[] for PETSc
	// char filename_char[filename.length() + 1]; 
    // strcpy(filename_char, filename.c_str()); 

	// // set up PETSc viewer
	// PetscViewer viewer;
	// PetscViewerCreate(PETSC_COMM_WORLD, &viewer);
	// PetscViewerSetType(viewer, PETSCVIEWERBINARY);
	// PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
	// PetscViewerFileSetMode(viewer, FILE_MODE_WRITE);
	// PetscViewerFileSetName(viewer, filename_char);
	// VecView(_vec, viewer)

	// Hannah: to do
}

// ----------------------------
// Acessors
// ----------------------------

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get value
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get element globalRow
real64 PetscVector::get(globalIndex globalRow) const
{
	double value[1];
	int index[1] = {globalRow};
	VecGetValues(_vec, 1, index, value);
  	return value[0];
}

void PetscVector::get( array1d<globalIndex> const & globalIndices,
                        array1d<real64> & values ) const
{
  // Hannah: to do
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get unwrapped pointer
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get const pointer
const Vec* PetscVector::unwrappedPointer() const
{
	return &(_vec);
}

// Get non-const pointer
const Vec* PetscVector::unwrappedPointer()
{
	return &(_vec);
}

// Get PETSc object
Vec PetscVector::getConstVec() const
{
	return _vec;
}

// Get PETSc object
Vec PetscVector::getVec()
{
	return _vec;
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get the number of global elements
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Return the global size of the vector (total number of elements).
int PetscVector::globalSize() const
{
	int size;
	VecGetSize(_vec, &size);
	return size;
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get the number of local elements
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Return the local size of the vector (total number of local elements).
int PetscVector::localSize() const
{
	int size;
	VecGetLocalSize(_vec, &size);
	return size;
}

} // end geosx


// END_RST_NARRATIVE
