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
 * @file EpetraVector.cpp
 *
 *  Created on: Jul 24, 2018
 *  Author: Matthias Cremon
 */

// BEGIN_RST_NARRATIVE EpetraSparseVector.rst
// ==============================
// Epetra-based Vector Object
// ==============================
// This class contains the ParallelVector wrappers based on Epetra_Vector Objects.
// The class contains a unique pointer to an Epetra_CrsVector as well as constructors,
// functions and accessors for Epetra objects.

// Include the corresponding header file.
#include "TrilinosVector.hpp"

#include <ostream>

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
EpetraVector::EpetraVector()
{}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Copy constructor
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Create a unique pointer to an Epetra_CrsMatrix.  The data from the input vector is
// copied to a new memory location. First checks if the vector to be copied is not empty.
EpetraVector::EpetraVector( EpetraVector const &in_vec )
{
  if( in_vec.getPointer() != nullptr )
  {
    m_vector = std::unique_ptr<Epetra_Vector>( new Epetra_Vector( *in_vec.getPointer()));
  }
}

// ----------------------------
// Create/finalize
// ----------------------------

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Create from array
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Create a vector from array
void EpetraVector::create( trilinosTypes::gid const size,
                           double *V )
{
  // Create an Epetra_Map of size size.
  Epetra_Map map = Epetra_Map( size, 0, Epetra_MpiComm( MPI_COMM_WORLD ));
  // Create a unique pointer to an Epetra_Vector defined from the Epetra_Map.
  m_vector = std::unique_ptr<Epetra_Vector>( new Epetra_Vector( View, map, V ));
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Create from Epetra_Map
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Create a vector from Epetra_Map
void EpetraVector::create( Epetra_Map const &Map,
                           double *V )
{
  // Create an Epetra_Map from the input map.
  Epetra_Map map = Epetra_Map( Map );
  // Create a unique pointer to an Epetra_Vector defined from the Epetra_Map.
  m_vector = std::unique_ptr<Epetra_Vector>( new Epetra_Vector( View, map, V ));
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Create from vector
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Create a vector from vector
void EpetraVector::create( std::vector<double> &vec )
{
  // Get the size of the vector
  trilinosTypes::gid m_size = vec.size();
  // Create an Epetra_Map from that size.
  Epetra_Map map = Epetra_Map( m_size, 0, Epetra_MpiComm( MPI_COMM_WORLD ));
  // Create a unique pointer to an Epetra_Vector defined from the Epetra_Map.
  m_vector = std::unique_ptr<Epetra_Vector>( new Epetra_Vector( View, map, vec.data()));
}

// Add into value (TODO This needs to use integers for some reason! No longlong).
void EpetraVector::add( integer const element,
                        real64 const value )
{
  m_vector->SumIntoGlobalValues( 1, &value, &element );
}

// Add into values (TODO This needs to use integers for some reason! No longlong).
void EpetraVector::add( array1d<integer> const elements,
                        array1d<real64> const values )
{
  m_vector->SumIntoGlobalValues( static_cast<integer>(elements.size()), values.data(), elements.data() );
}

// Add values into vector
void EpetraVector::add( lid const numIndices,
                        gid const * const elements,
                        real64 const * const values )
{
  m_vector->SumIntoGlobalValues( static_cast<int>(numIndices), values, elements );
}


// Set value
void EpetraVector::set( trilinosTypes::gid const element,
                        real64 const value )
{
  m_vector->ReplaceGlobalValues( 1, &value, &element );
}

// Set values
void EpetraVector::set( array1d<trilinosTypes::gid> const elements,
                        array1d<real64> const values )
{
  m_vector->ReplaceGlobalValues( static_cast<integer>(elements.size()), values.data(), elements.data() );
}

// ----------------------------
// Linear Algebra
// ----------------------------

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Scale
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Multiply all elements by scalingFactor.
void EpetraVector::scale( real64 const scalingFactor )
{
  // Scale every element in the vector.
  m_vector.get()->Scale( scalingFactor );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Dot
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Dot product with the vector vec.
void EpetraVector::dot( EpetraVector const &vec,
                        real64 &dst )
{
  // Compute a dot product with vector vec. Put the result in dst.
  m_vector.get()->Dot( *vec.getPointer(), &dst );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Copy
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Update vector as this = x.
void EpetraVector::copy( EpetraVector const &x )
{
  // Copy the input vector. Points to the same memory location but update the values.
  m_vector.get()->Update( 1., *x.getPointer(), 0. );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Axpy
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Update vector as this = alpha*x + this.
void EpetraVector::axpy( real64 const alpha,
                          EpetraVector const &x )
{
  // Update the vector. Points to the same memory location but update the values.
  m_vector.get()->Update( alpha, *x.getPointer(), 1. );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Axpby
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Update vector as this = alpha*x + beta*this.
void EpetraVector::axpby( real64 const alpha,
                          EpetraVector const &x,
                          real64 const beta )
{
  // Update the vector. Points to the same memory location but update the values.
  m_vector.get()->Update( alpha, *x.getPointer(), beta );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// 1-norm
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// 1-norm of the vector.
void EpetraVector::norm1( real64 &dst ) const
{
  // Compute the 1-norm of the vector and put the result in dst.
  m_vector.get()->Norm1( &dst );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// 2-norm
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// 2-norm of the vector.
void EpetraVector::norm2( real64 &dst ) const
{
  // Compute the 2-norm of the vector and put the result in dst.
  m_vector.get()->Norm2( &dst );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Inf-norm
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Inf-norm of the vector.
void EpetraVector::normInf( real64 &dst ) const
{
  // Compute the inf-norm of the vector and put the result in dst.
  m_vector.get()->NormInf( &dst );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Print
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Print vector to the terminal in Trilinos format.
void EpetraVector::print( std::ostream & outputStream ) const
{
  if( m_vector.get() != nullptr )
    outputStream << *m_vector.get() << std::endl;
}

// ----------------------------
// Acessors
// ----------------------------


real64 const * EpetraVector::getValues() const
{
  return m_vector->Values();
}

real64 * EpetraVector::getValues()
{
  return m_vector->Values();
}


// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get element
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get element i
real64 EpetraVector::getElement(trilinosTypes::gid i) const
{
  real64 * temp = m_vector->Values();
  return temp[i];
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get pointer (const)
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get const pointer
const Epetra_Vector* EpetraVector::getPointer() const
{
  return m_vector.get();
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get pointer
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get non-const pointer
Epetra_Vector* EpetraVector::getPointer()
{
  return m_vector.get();
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get the number of global elements
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Return the global size of the vector (total number of elements).
trilinosTypes::gid EpetraVector::globalSize() const
{
  return m_vector.get()->GlobalLength64();
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get the number of local elements
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Return the local size of the vector (total number of local elements).
trilinosTypes::lid EpetraVector::localSize() const
{
  return m_vector.get()->MyLength();
}

}

// END_RST_NARRATIVE
