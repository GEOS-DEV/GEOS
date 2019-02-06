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
 * @file EpetraVector.cpp
 *
 *  Created on: Jul 24, 2018
 *  Author: Matthias Cremon
 */

// BEGIN_RST_NARRATIVE EpetraSparseVector.rst
// ==============================
// Epetra-based Vector Object
// ==============================
// This class contains the ParallelVector wrappers for Epetra_FEVector Objects.
// The class contains a unique pointer to an Epetra_FEVector as well as constructors,
// functions and accessors for Epetra objects.

// Include the corresponding header file.
#include "EpetraVector.hpp"

// Put everything under the geosx namespace.
namespace geosx
{

// ----------------------------
// Constructors
// ----------------------------

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Empty constructor
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Construct as an empty vector
EpetraVector::EpetraVector()
{}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Copy constructor
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Create a unique pointer to the raw vector.  The data from the input vector is
// copied to a new memory location. Checks if the vector to be copied is empty.
EpetraVector::EpetraVector( EpetraVector const &src )
{
  GEOS_ERROR_IF( src.unwrappedPointer() == nullptr, "source vector appears to be empty" );
  m_vector = std::unique_ptr<Epetra_FEVector>( new Epetra_FEVector( *src.unwrappedPointer()));
}

// ----------------------------
// Create
// ----------------------------

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Create from Epetra_Map
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Create a vector from an Epetra_Map.
void EpetraVector::create( Epetra_Map const &map )
{
  m_vector = std::unique_ptr<Epetra_FEVector>( new Epetra_FEVector( map ));
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
void EpetraVector::createWithLocalSize( localIndex const localSize, MPI_Comm const & comm )
{
  Epetra_Map map = Epetra_Map( -1, integer_conversion<int, localIndex>( localSize ), 0, Epetra_MpiComm( comm ) );
  create( map );
}

void EpetraVector::createWithGlobalSize( globalIndex const globalSize, MPI_Comm const & comm )
{
  Epetra_Map map = Epetra_Map( globalSize, 0, Epetra_MpiComm( comm ) );
  create( map );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Create from array
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Create a vector from local array data.  The global vector contains
// local arrays stitched together.
void EpetraVector::create( array1d<real64> const & localValues, MPI_Comm const & comm )
{
  int localSize = integer_conversion<int, localIndex>( localValues.size());
  Epetra_Map map = Epetra_Map( -1, localSize, 0, Epetra_MpiComm( comm ) );
  m_vector = std::unique_ptr<Epetra_FEVector>( new Epetra_FEVector( View, map, const_cast<double*>(localValues.data()), localSize, 1 ));
}


// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Add/Set value(s)
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Add/set entries in the vector.  Note that Epetra_Vector seems to have
// haphazard support for long ints.  This motivated using FEVector to be
// sure we have proper globalRow support.

// single element options
void EpetraVector::set( globalIndex const globalRow,
                        real64 const value )
{
  m_vector->ReplaceGlobalValues( 1, &globalRow, &value );
}

void EpetraVector::add( globalIndex const globalRow,
                        real64 const value )
{
  m_vector->SumIntoGlobalValues( 1, &globalRow, &value );
}

// n-element, c-style options
void EpetraVector::set( globalIndex const * globalIndices,
                        real64 const * values,
                        localIndex size )
{
  m_vector->ReplaceGlobalValues( integer_conversion<int, localIndex>( size ), globalIndices, values );
}

void EpetraVector::add( globalIndex const * globalIndices,
                        real64 const * values,
                        localIndex size )
{
  m_vector->SumIntoGlobalValues( integer_conversion<int, localIndex>( size ), globalIndices, values );
}

// n-element, array1d options
void EpetraVector::set( array1d<globalIndex> const & globalIndices,
                        array1d<real64> const & values )
{
  m_vector->ReplaceGlobalValues( integer_conversion<int, localIndex>( values.size()), globalIndices.data(), values.data() );
}
void EpetraVector::add( array1d<globalIndex> const & globalIndices,
                        array1d<real64> const & values )
{
  m_vector->SumIntoGlobalValues( integer_conversion<int, localIndex>( values.size()), globalIndices.data(), values.data() );
}

//additional convenience options:
void EpetraVector::set( real64 value )
{
  m_vector->PutScalar( value );
}

void EpetraVector::zero()
{
  set( 0.0 );
}

void EpetraVector::rand()
{
  m_vector->SetSeed( 1984 );
  m_vector->Random();
}


// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Open / close
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
void EpetraVector::open()
{
// ... nothing to do here ...
}

void EpetraVector::close()
{
  m_vector->GlobalAssemble();
}

// ---------------------------------------------------------
// Linear Algebra
// ---------------------------------------------------------
// The following functions support basic linear algebra ops

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Scale
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Multiply all elements by scalingFactor.
void EpetraVector::scale( real64 const scalingFactor )
{
  m_vector.get()->Scale( scalingFactor );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Dot
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Dot product with the vector vec.
real64 EpetraVector::dot( EpetraVector const &vec )
{
  real64 tmp;
  m_vector.get()->Dot( *vec.unwrappedPointer(), &tmp );
  return tmp;
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Copy
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Update vector as this = x.
void EpetraVector::copy( EpetraVector const &x )
{
  m_vector.get()->Update( 1., *x.unwrappedPointer(), 0. );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Axpy
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Update vector as this = alpha*x + this.
void EpetraVector::axpy( real64 const alpha,
                         EpetraVector const &x )
{
  m_vector.get()->Update( alpha, *x.unwrappedPointer(), 1. );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Axpby
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Update vector as this = alpha*x + beta*this.
void EpetraVector::axpby( real64 const alpha,
                          EpetraVector const &x,
                          real64 const beta )
{
  m_vector.get()->Update( alpha, *x.unwrappedPointer(), beta );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// 1-norm
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// 1-norm of the vector.
real64 EpetraVector::norm1() const
{
  real64 tmp;
  m_vector.get()->Norm1( &tmp );
  return tmp;
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// 2-norm
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// 2-norm of the vector.
real64 EpetraVector::norm2() const
{
  real64 tmp;
  m_vector.get()->Norm2( &tmp );
  return tmp;
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Inf-norm
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Inf-norm of the vector.
real64 EpetraVector::normInf() const
{
  real64 tmp;
  m_vector.get()->NormInf( &tmp );
  return tmp;
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Print
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Print vector to the std::cout in Trilinos format.
void EpetraVector::print() const
{
  GEOS_ERROR_IF( m_vector.get() == nullptr, "Vector appears to be empty" );
  std::cout << *m_vector.get() << std::endl;
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Write to matlab-compatible file
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Note: EpetraExt also supports a MatrixMarket format as well
//       if we prefer that.
void EpetraVector::write( string const & filename ) const
{
  EpetraExt::MultiVectorToMatlabFile( filename.c_str(), *m_vector );
}

// ----------------------------
// Acessors
// ----------------------------

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get value
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get element globalRow
// TODO: implementation not straightforward
real64 EpetraVector::get( globalIndex globalRow ) const
{
  GEOS_ERROR( "not yet implemented" );
  return std::numeric_limits<double>::quiet_NaN();
}

void EpetraVector::get( array1d<globalIndex> const & globalIndices,
                        array1d<real64> & values ) const
{
  GEOS_ERROR( "not yet implemented" );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get unwrapped pointer
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get pointer to raw Epetra object, with const and non-const versions.
Epetra_FEVector const * EpetraVector::unwrappedPointer() const
{
  return m_vector.get();
}

Epetra_FEVector* EpetraVector::unwrappedPointer()
{
  return m_vector.get();
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get the number of global elements
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Return the global size of the vector (total number of elements).
globalIndex EpetraVector::globalSize() const
{
  return m_vector.get()->GlobalLength64();
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get the number of local elements
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Return the local size of the vector (total number of local elements).
localIndex EpetraVector::localSize() const
{
  return m_vector.get()->MyLength();
}

} // end geosx

// END_RST_NARRATIVE
