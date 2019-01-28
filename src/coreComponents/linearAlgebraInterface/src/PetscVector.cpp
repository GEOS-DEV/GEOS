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
// Create a unique PETSc Vec.  The data from the input vector is
// copied to a new memory location. 
PetscVector::PetscVector(PetscVector const & vec)
{
	VecDuplicate(vec._vec, &_vec); 
	VecCopy(vec._vec, _vec); 
}

PetscVector::PetscVector(Vec vec)
{
	_vec = vec;
}

PetscVector::~PetscVector()
{
  if (_vec) VecDestroy(&_vec);
}

// ----------------------------
// Create/finalize
// ----------------------------

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Create from array
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Create a vector from array
void PetscVector::create(const int size, real64 *vals)
{

  int indices[size];
  for (int i = 0; i < size; i++){
    indices[i] = i;
  }

  VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, size, &_vec);
  VecSetValues(_vec, size, indices, vals, INSERT_VALUES);
  VecAssemblyBegin(_vec);
  VecAssemblyEnd(_vec);
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Create from vector
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Create a vector from vector
void PetscVector::create(std::vector<real64> &vec)
{
  int size = vec.size();

  int indices[size];
  for (int i = 0; i < size; i++){
    indices[i] = i;
  }

  VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, size, &_vec);
  VecSetValues(_vec, size, indices, vec.data(), INSERT_VALUES);
  VecAssemblyBegin(_vec);
  VecAssemblyEnd(_vec);
}

/* set value of vector element */
void PetscVector::set(int element, real64 value)
{
	VecSetValue(_vec, element, value, INSERT_VALUES);
	VecAssemblyBegin(_vec);
  	VecAssemblyEnd(_vec);
}

// /* set values of vector elements */
// void PetscVector::set(array1d<int> elements, array1d<real64> values){}

/* add value to vector element */
void PetscVector::add(int element, real64 value)
{
	VecSetValue(_vec, element, value, ADD_VALUES);
	VecAssemblyBegin(_vec);
  	VecAssemblyEnd(_vec);
}

// /* add values to vector elements */
// void PetscVector::set(array1d<int> elements, array1d<real64> values){}

// ----------------------------
// Linear Algebra
// ----------------------------

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
void PetscVector::dot(PetscVector const vec, real64 *dst)
{
	VecDot(_vec, vec._vec, dst);
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
void PetscVector::norm1(real64 &result) const
{
	VecNorm(_vec, NORM_1, &result);
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// 2-norm
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// 2-norm of the vector.
void PetscVector::norm2(real64 &result) const
{
	VecNorm(_vec, NORM_2, &result);
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Inf-norm
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Inf-norm of the vector.
void PetscVector::normInf(real64 &result) const
{
	VecNorm(_vec, NORM_INFINITY, &result);
}

// ----------------------------
// Acessors
// ----------------------------

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

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get pointer (const)
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get const pointer
const Vec* PetscVector::getPointer() const
{
	return &(_vec);
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get pointer
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get non-const pointer
const Vec* PetscVector::getPointer()
{
	return &(_vec);
}

Vec PetscVector::getConstVec() const
{
	return _vec;
}

Vec PetscVector::getVec()
{
	return _vec;
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Print
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Print vector to the terminal in Trilinos format.
void PetscVector::print() const
{
	VecView(_vec, PETSC_VIEWER_STDOUT_WORLD);
}

// END_RST_NARRATIVE
