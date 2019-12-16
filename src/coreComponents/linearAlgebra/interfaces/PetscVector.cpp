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
 * @file PetscVector.cpp
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

#include <petscvec.h>

// Put everything under the geosx namespace.
namespace geosx
{

static_assert( sizeof(PetscInt) == sizeof(globalIndex), "sizeof(PetscInt) != sizeof(globalIndex)");
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
// Empty constructor
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Create an empty vector
PetscVector::PetscVector()
: m_vec()
{}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Copy constructor
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Create a unique PetscVector from PetscVector.  
// The data from the input vector is copied to a new memory location. 
PetscVector::PetscVector( PetscVector const & vec )
: m_vec()
{
  VecDuplicate( vec.m_vec, &m_vec );
  VecCopy( vec.m_vec, m_vec );
}

// Create a unique PetscVector from a PETSc Vec.  
PetscVector::PetscVector( Vec vec )
{
  m_vec = vec;
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Destructor
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
PetscVector::~PetscVector()
{
  VecDestroy( &m_vec );
}

// ----------------------------
// Create
// ----------------------------


// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Create from PetscVector
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
void PetscVector::create( PetscVector const & src )
{
  VecDestroy( &m_vec );
  VecDuplicate( src.m_vec, &m_vec );
  VecCopy( src.m_vec, m_vec );
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
  VecDestroy( &m_vec );
  VecCreate( comm, &m_vec );
  VecSetType( m_vec, VECMPI );
  VecSetSizes( m_vec, localSize, PETSC_DETERMINE);
}

void PetscVector::createWithGlobalSize( globalIndex const globalSize, MPI_Comm const & comm )
{
  VecDestroy( &m_vec );
  VecCreate( comm, &m_vec );
  VecSetType( m_vec, VECMPI );
  VecSetSizes( m_vec, PETSC_DECIDE, globalSize);
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

  VecDestroy( &m_vec );
  VecCreate( comm, &m_vec );
  VecSetType( m_vec, VECMPI );
  VecSetSizes( m_vec, size, PETSC_DETERMINE);
  VecGetArray( m_vec, &values );

  // set vector values
  for (int i = 0; i < size; i++)
  {
    values[i] = localValues[i];
  }

  VecRestoreArray( m_vec, &values );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Add/Set value(s)
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Add/set entries in the vector.

// single element
void PetscVector::set( globalIndex const globalRow,
                       real64 const value )
{
  VecSetValue( m_vec, globalRow, value, INSERT_VALUES );
}

void PetscVector::add( globalIndex const globalRow,
                       real64 const value )
{
  VecSetValue( m_vec, globalRow, value, ADD_VALUES );
}

// n-element, c-style options
void PetscVector::set( globalIndex const * globalIndices,
                       real64 const * values,
                       localIndex size )
{
  VecSetValues( m_vec, size, toPetscInt( globalIndices), values, INSERT_VALUES );
}

void PetscVector::add( globalIndex const * globalIndices,
                       real64 const * values,
                       localIndex size )
{
  VecSetValues( m_vec, size, toPetscInt( globalIndices), values, ADD_VALUES );
}

// n-element, array1d options
void PetscVector::set( array1d<globalIndex> const & globalIndices,
                       array1d<real64> const & values )
{
  VecSetValues( m_vec, values.size(), toPetscInt( globalIndices.data()), values.data(), INSERT_VALUES );
}

void PetscVector::add( array1d<globalIndex> const & globalIndices,
                       array1d<real64> const & values )
{
  VecSetValues( m_vec, values.size(), toPetscInt( globalIndices.data()), values.data(), ADD_VALUES );
}

// additional convenience options
void PetscVector::set( real64 value )
{
  VecSet( m_vec, value );
}

void PetscVector::zero()
{
  set( 0.0 );
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
  VecSetRandom( m_vec, ran );
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
  VecAssemblyBegin( m_vec );
  VecAssemblyEnd( m_vec );
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
  VecScale( m_vec, scalingFactor );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Dot
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Dot product with the vector vec.
real64 PetscVector::dot( PetscVector const &vec )
{
  real64 dot;
  VecDot( m_vec, vec.m_vec, &dot );
  return dot;
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Copy
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Update vector as this = x.
void PetscVector::copy( PetscVector const &x )
{
  VecSet( m_vec, 0 );
  VecAXPY( m_vec, 1.0, x.m_vec );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Axpy
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Update vector as this = alpha*x + this.
void PetscVector::axpy( real64 const alpha, PetscVector const &x )
{
  VecAXPY( m_vec, alpha, x.m_vec );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Axpby
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Update vector as this = alpha*x + beta*this.
void PetscVector::axpby( real64 const alpha,
                         PetscVector &x,
                         real64 const beta )
{
  VecScale( m_vec, beta );
  VecAXPY( m_vec, alpha, x.m_vec );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// 1-norm
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// 1-norm of the vector.
real64 PetscVector::norm1() const
{
  real64 result;
  VecNorm( m_vec, NORM_1, &result );
  return result;
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// 2-norm
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// 2-norm of the vector.
real64 PetscVector::norm2() const
{
  real64 result;
  VecNorm( m_vec, NORM_2, &result );
  return result;
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Inf-norm
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Inf-norm of the vector.
real64 PetscVector::normInf() const
{
  real64 result;
  VecNorm( m_vec, NORM_INFINITY, &result );
  return result;
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Print
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Print vector to the terminal in PETSc format.
void PetscVector::print() const
{
  VecView( m_vec, PETSC_VIEWER_STDOUT_WORLD );
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
  VecView( m_vec, viewer );
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
  real64 value;
  VecGetValues( m_vec, 1, toPetscInt( &globalRow ), &value );
  return value;
}

void PetscVector::get( array1d<globalIndex> const & globalIndices,
                       array1d<real64> & values ) const
{
  values.resize( globalIndices.size() );
  VecGetValues( m_vec, globalIndices.size(), toPetscInt( globalIndices.data() ), values.data() );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get unwrapped pointer
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get const pointer
const Vec* PetscVector::unwrappedPointer() const
{
  return &( m_vec );
}

// Get non-const pointer
Vec* PetscVector::unwrappedPointer()
{
  return &( m_vec );
}

// Get const PETSc object
Vec PetscVector::getConstVec() const
{
  return m_vec;
}

// Get PETSc object
Vec PetscVector::getVec()
{
  return m_vec;
}

// Accessor for the MPI communicator
MPI_Comm PetscVector::getComm() const
{
  MPI_Comm comm;
  PetscObjectGetComm( reinterpret_cast<PetscObject>( m_vec ), &comm );
  return comm;
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get the number of global elements
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Return the global size of the vector (total number of elements).
globalIndex PetscVector::globalSize() const
{
  PetscInt size;
  VecGetSize( m_vec, &size );
  return size;
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get the number of local elements
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Return the local size of the vector (total number of local elements).
localIndex PetscVector::localSize() const
{
  PetscInt size;
  VecGetLocalSize( m_vec, &size );
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
  VecGetOwnershipRange( m_vec, &low, &high );
  GEOSX_ERROR_IF( index < low || high <= index, "getLocalRowID: processor does not own global row index" );
  return index - low;
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// extractLocalVector
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Extract a view of the local portion of the array
real64 const * PetscVector::extractLocalVector() const
{
  PetscScalar * avec;
  VecGetArray( m_vec, &avec );
  real64 const * localVector = avec;
  VecRestoreArray( m_vec, &avec );
  return localVector;
}

real64 * PetscVector::extractLocalVector()
{
  PetscScalar * avec;
  VecGetArray( m_vec, &avec );
  real64 * localVector = avec;
  VecRestoreArray( m_vec, &avec );
  return localVector;
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// ilower
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Returns the index of the first global row owned by that processor.
globalIndex PetscVector::ilower() const
{
  PetscInt low, high;
  VecGetOwnershipRange( m_vec, &low, &high );
  return low;
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// iupper
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// eturns the next index after last global row owned by that processor.
globalIndex PetscVector::iupper() const
{
  PetscInt low, high;
  VecGetOwnershipRange( m_vec, &low, &high );
  return high;
}

} // end geosx

// END_RST_NARRATIVE
