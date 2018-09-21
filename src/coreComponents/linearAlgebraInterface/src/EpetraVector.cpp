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

#include "EpetraVector.hpp"

namespace geosx
{
// Create an empty vector
EpetraVector::EpetraVector()
{}

// Copy a vector
EpetraVector::EpetraVector( EpetraVector const &in_vec )
{
  if( in_vec.getPointer() != nullptr )
    m_vector = std::unique_ptr<Epetra_Vector>( new Epetra_Vector( *in_vec.getPointer()));
}

// Create a vector from array
void EpetraVector::create( const trilinosTypes::gid size, double *V )
{
  Epetra_Map map = Epetra_Map( size, 0, Epetra_MpiComm( MPI_COMM_WORLD ));
  m_vector = std::unique_ptr<Epetra_Vector>( new Epetra_Vector( View, map, V ));
}

// Create a vector from array
void EpetraVector::create( const Epetra_Map& Map, double *V )
{
  Epetra_Map map = Epetra_Map( Map );
  m_vector = std::unique_ptr<Epetra_Vector>( new Epetra_Vector( View, map, V ));
}

// Create a vector from vector
void EpetraVector::create( std::vector<double> &vec )
{
  trilinosTypes::gid m_size = vec.size();
  Epetra_Map map = Epetra_Map( m_size, 0, Epetra_MpiComm( MPI_COMM_WORLD ));
  m_vector = std::unique_ptr<Epetra_Vector>( new Epetra_Vector( View, map, vec.data()));
}

// Multiply all elements by scalingFactor.
void EpetraVector::scale( real64 const scalingFactor )
{
  m_vector.get()->Scale( scalingFactor );
}
// Dot product with the vector vec.
void EpetraVector::dot( EpetraVector const &vec,
                        real64 &dst )
{
  m_vector.get()->Dot( *vec.getPointer(), &dst );
}

// Update vector as this = x.
void EpetraVector::copy( EpetraVector const &x )
{
  m_vector.get()->Update( 1., *x.getPointer(), 0. );
}

// Update vector as this = alpha*x + beta*this.
void EpetraVector::axpby( real64 const alpha,
                          EpetraVector const &x,
                          real64 const beta )
{
  m_vector.get()->Update( alpha, *x.getPointer(), beta );
}

// 1-norm of the vector.
void EpetraVector::norm1( real64 &dst ) const
{
  m_vector.get()->Norm1( &dst );
}

// 2-norm of the vector.
void EpetraVector::norm2( real64 &dst ) const
{
  m_vector.get()->Norm2( &dst );
}

// Inf-norm of the vector.
void EpetraVector::normInf( real64 &dst ) const
{
  m_vector.get()->NormInf( &dst );
}


// Return the global size of the vector (total number of elements).
trilinosTypes::gid EpetraVector::globalSize() const
{
  return m_vector.get()->GlobalLength64();
}

// Return the local size of the vector (total number of local elements).
trilinosTypes::lid EpetraVector::localSize() const
{
  return m_vector.get()->MyLength();
}

// Print vector to the terminal in Trilinos format.
void EpetraVector::print() const
{
  if( m_vector.get() != nullptr )
    std::cout << *m_vector.get() << std::endl;
}

// Accessors
// Get const pointer
const Epetra_Vector* EpetraVector::getPointer() const
{
  return m_vector.get();
}

// Get non-const pointer
Epetra_Vector* EpetraVector::getPointer()
{
  return m_vector.get();
}

}
