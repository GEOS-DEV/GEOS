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
 * @file PetscVector.cpp
 *
 *  Created on: Jan 28, 2019
 *  Author: Hannah Morgan
 */

// BEGIN_RST_NARRATIVE PetscVector.rst
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
{}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Copy constructor
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Create a unique PetscVector from PetscVector.  
// The data from the input vector is copied to a new memory location. 
PetscVector::PetscVector( PetscVector const & vec )
{
  VecDuplicate( vec._vec, &_vec ); 
  VecCopy( vec._vec, _vec ); 
}

// Create a unique PetscVector from a PETSc Vec.  
PetscVector::PetscVector( Vec vec )
{
  _vec = vec;
}

// ----------------------------
// Create
// ----------------------------


// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Create from PetscVector
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
void PetscVector::create( PetscVector const & src )
{
  VecDuplicate( src._vec, &_vec ); 
  VecCopy( src._vec, _vec );  
}

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
// 
// NOTE: creates a CPU MPI vector
void PetscVector::createWithLocalSize( localIndex const localSize, MPI_Comm const & comm )
{
  VecCreate( comm, &_vec );
  VecSetType( _vec, VECMPI );
  VecSetSizes( _vec, localSize, PETSC_DETERMINE);
}

void PetscVector::createWithGlobalSize( globalIndex const globalSize, MPI_Comm const & comm )
{
  VecCreate( comm, &_vec );
  VecSetType( _vec, VECMPI );
  VecSetSizes( _vec, PETSC_DECIDE, globalSize);
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Create from array
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Create a vector from local array data.  The global vector contains
// local arrays stitched together.
//
// NOTE: creates a CPU MPI vector, must assemble vector after use
void PetscVector::create( array1d<real64> const & localValues, MPI_Comm const & comm )
{
  PetscInt size = localValues.size(); 
  PetscScalar *values;

  VecCreate( comm, &_vec );
  VecSetType( _vec, VECMPI );
  VecSetSizes( _vec, size, PETSC_DETERMINE);
  VecGetArray( _vec, &values );

  // set vector values
  for (int i = 0; i < size; i++)
  {
    values[i] = localValues[i];
  }

  VecRestoreArray( _vec, &values );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Add/Set value(s)
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Add/set entries in the vector.  

// single element 
void PetscVector::set( globalIndex const globalRow,
                       real64 const value )
{
  VecSetValue( _vec, globalRow, value, INSERT_VALUES );
}

void PetscVector::add( globalIndex const globalRow,
                       real64 const value )
{
  VecSetValue( _vec, globalRow, value, ADD_VALUES );
}

// n-element, c-style options
void PetscVector::set( globalIndex const * globalIndices,
                       real64 const * values,
                       localIndex size )
{
  VecSetValues( _vec, size, toPetscInt(globalIndices), values, INSERT_VALUES );
}

void PetscVector::add( globalIndex const * globalIndices,
                       real64 const * values,
                       localIndex size )
{
  VecSetValues( _vec, size, toPetscInt(globalIndices), values, ADD_VALUES );
}

// n-element, array1d options
void PetscVector::set( array1d<globalIndex> const & globalIndices,
                       array1d<real64> const & values )
{
  VecSetValues( _vec, values.size(), toPetscInt(globalIndices.data()), values.data(), INSERT_VALUES );
}

void PetscVector::add( array1d<globalIndex> const & globalIndices,
                       array1d<real64> const & values )
{
  VecSetValues( _vec, values.size(), toPetscInt(globalIndices.data()), values.data(), ADD_VALUES );
}

// additional convenience options
void PetscVector::set( real64 value )
{
  VecSet( _vec, value );
}

void PetscVector::zero()
{
  VecZeroEntries( _vec );
}

void PetscVector::rand( unsigned long seed )
{
  // create random context
  PetscRandom ran;
  PetscRandomCreate( PETSC_COMM_WORLD, &ran );

  // set random seed
  PetscRandomSetSeed( ran, seed );    
  PetscRandomSeed( ran );

  // create random vector
  VecSetRandom( _vec, ran );
  PetscRandomDestroy( &ran );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Open / close
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
void PetscVector::open()
{}

void PetscVector::close()
{
  // assemble the vector after setting values
  VecAssemblyBegin( _vec );
  VecAssemblyEnd( _vec ); 
}

// ----------------------------
// Linear Algebra
// ----------------------------
// The following functions support basic linear algebra ops

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Scale
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Multiply all elements by scalingFactor.
void PetscVector::scale( real64 const scalingFactor )
{
  VecScale( _vec, scalingFactor );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Dot
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Dot product with the vector vec.
real64 PetscVector::dot( PetscVector const &vec )
{
  real64 dot;
  VecDot( _vec, vec._vec, &dot );
  return dot;
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Copy
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Update vector as this = x.
void PetscVector::copy( PetscVector const &x )
{
  VecSet( _vec, 0 );
  VecAXPY( _vec, 1.0, x._vec );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Axpy
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Update vector as this = alpha*x + this.
void PetscVector::axpy( real64 const alpha, PetscVector const &x )
{
  VecAXPY( _vec, alpha, x._vec );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Axpby
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Update vector as this = alpha*x + beta*this.
void PetscVector::axpby( real64 const alpha, 
                         PetscVector &x, 
                         real64 const beta )
{
  VecScale( _vec, beta );
  VecAXPY( _vec, alpha, x._vec );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// 1-norm
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// 1-norm of the vector.
real64 PetscVector::norm1() const
{
  real64 result;
  VecNorm( _vec, NORM_1, &result );
  return result;
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// 2-norm
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// 2-norm of the vector.
real64 PetscVector::norm2() const
{
  real64 result;
  VecNorm( _vec, NORM_2, &result );
  return result;
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Inf-norm
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Inf-norm of the vector.
real64 PetscVector::normInf() const
{
  real64 result;
  VecNorm( _vec, NORM_INFINITY, &result );
  return result;
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Print
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Print vector to the terminal in PETSc format.
void PetscVector::print() const
{
  VecView( _vec, PETSC_VIEWER_STDOUT_WORLD );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Write to matlab-compatible file
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
void PetscVector::write( string const & filename,
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
  } else 
  {
    PetscViewerASCIIOpen( getComm(), filename.c_str(), &viewer);
    PetscViewerPushFormat( viewer, PETSC_VIEWER_ASCII_MATLAB );
  }
  VecView( _vec, viewer );
}

// ----------------------------
// Acessors
// ----------------------------

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get value
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get element globalRow
real64 PetscVector::get( globalIndex globalRow ) const
{
  real64 value[1];
  PetscInt index[1] = {globalRow};
  VecGetValues( _vec, 1, index, value );
  return value[0];
}

void PetscVector::get( array1d<globalIndex> const & globalIndices,
                       array1d<real64> & values ) const
{
  GEOS_ERROR( "not yet implemented" );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get unwrapped pointer
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get const pointer
const Vec* PetscVector::unwrappedPointer() const
{
  return &( _vec );
}

// Get non-const pointer
Vec* PetscVector::unwrappedPointer()
{
  return &( _vec );
}

// Get const PETSc object
Vec PetscVector::getConstVec() const
{
  return _vec;
}

// Get PETSc object
Vec PetscVector::getVec()
{
  return _vec;
}

// Accessor for the MPI communicator
MPI_Comm PetscVector::getComm() const
{
  MPI_Comm comm;
  PetscObjectGetComm( reinterpret_cast<PetscObject>( _vec ), &comm );
  return comm;
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get the number of global elements
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Return the global size of the vector (total number of elements).
globalIndex PetscVector::globalSize() const
{
  PetscInt size;
  VecGetSize( _vec, &size );
  return size;
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get the number of local elements
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Return the local size of the vector (total number of local elements).
localIndex PetscVector::localSize() const
{
  PetscInt size;
  VecGetLocalSize( _vec, &size );
  return size;
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// getLocalRowID
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Map a global row index to local row index
//
// NOTE: error if requesting processor does not own row index
localIndex PetscVector::getLocalRowID( globalIndex const index ) const
{
  PetscInt low, high;
  VecGetOwnershipRange( _vec, &low, &high );
  if ( index < low || high <= index ) 
  {
    GEOS_ERROR( "getLocalRowID: processor does not own global row index" );
  } 
  return index - low; 
} 

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// extractLocalVector
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Extract a view of the local portion of the array
// User allocates memory for localVector
void PetscVector::extractLocalVector( real64 ** localVector ) const
{
  PetscScalar *avec;
  localIndex size = localSize();

  // get the vector
  VecGetArray( _vec, &avec );

  // copy
  for( int i = 0; i < size; i++ )
  {
    ( *localVector )[i] = avec[i];
  }

  // return local vector
  VecRestoreArray( _vec, &avec );
} 

} // end geosx

// END_RST_NARRATIVE
