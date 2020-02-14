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
 * @file EpetraVector.cpp
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

// Include required Epetra headers
#include <Epetra_FEVector.h>
#include <Epetra_Map.h>
#include <EpetraExt_MultiVectorOut.h>

#ifdef GEOSX_USE_MPI
#include <Epetra_MpiComm.h>
#else
#include<Epetra_SerialComm.h>
typedef Epetra_SerialComm Epetra_MpiComm;
#endif

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
EpetraVector::EpetraVector() = default;

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Copy constructor
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Create a unique pointer to the raw vector.  The data from the input vector is
// copied to a new memory location. Checks if the vector to be copied is empty.
EpetraVector::EpetraVector( EpetraVector const & src )
{
  GEOSX_ERROR_IF( src.unwrappedPointer() == nullptr, "source vector appears to be empty" );
  m_vector = std::make_unique< Epetra_FEVector >( *src.unwrappedPointer() );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Destructor
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Deletes the underlying Epetra vector
EpetraVector::~EpetraVector() = default;

// ----------------------------
// Create
// ----------------------------

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Create from EpetraVector
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
void EpetraVector::create( EpetraVector const & src )
{
  GEOSX_ERROR_IF( src.unwrappedPointer() == nullptr, "source vector appears to be empty" );
  m_vector = std::make_unique< Epetra_FEVector >( *src.unwrappedPointer() );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Create from Epetra_Map
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Create a vector from an Epetra_Map.
void EpetraVector::create( Epetra_Map const & map )
{
  m_vector = std::make_unique< Epetra_FEVector >( map );
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
void EpetraVector::createWithLocalSize( localIndex const localSize,
                                        MPI_Comm const & MPI_PARAM(comm) )
{
  Epetra_Map map = Epetra_Map( integer_conversion< globalIndex >( -1 ),
                               integer_conversion< int >( localSize ),
                               0,
                               Epetra_MpiComm( MPI_PARAM(comm) ) );
  create( map );
}

void EpetraVector::createWithGlobalSize( globalIndex const globalSize,
                                         MPI_Comm const & MPI_PARAM(comm) )
{
  Epetra_Map map = Epetra_Map( globalSize, 0, Epetra_MpiComm( MPI_PARAM(comm) ) );
  create( map );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Create from array
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Create a vector from local array data.  The global vector contains
// local arrays stitched together.
void EpetraVector::create( array1d< real64 > const & localValues,
                           MPI_Comm const & MPI_PARAM(comm) )
{
  int const localSize = integer_conversion< int >( localValues.size() );
  Epetra_Map map = Epetra_Map( -1, localSize, 0, Epetra_MpiComm( MPI_PARAM(comm) ) );
  m_vector = std::make_unique< Epetra_FEVector >( View, map, localValues.data(), localSize, 1 );
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
  m_vector->ReplaceGlobalValues( integer_conversion< int >( size ), globalIndices, values );
}

void EpetraVector::add( globalIndex const * globalIndices,
                        real64 const * values,
                        localIndex size )
{
  m_vector->SumIntoGlobalValues( integer_conversion< int >( size ), globalIndices, values );
}

// n-element, array1d options
void EpetraVector::set( array1d< globalIndex > const & globalIndices,
                        array1d< real64 > const & values )
{
  m_vector->ReplaceGlobalValues( integer_conversion< int >( values.size() ), globalIndices.data(), values.data() );
}

void EpetraVector::add( array1d< globalIndex > const & globalIndices,
                        array1d< real64 > const & values )
{
  m_vector->SumIntoGlobalValues( integer_conversion< int >( values.size() ), globalIndices.data(), values.data() );
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
real64 EpetraVector::dot( EpetraVector const & vec ) const
{
  real64 tmp;
  m_vector.get()->Dot( *vec.unwrappedPointer(), &tmp );
  return tmp;
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Copy
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Update vector as this = x.
void EpetraVector::copy( EpetraVector const & x )
{
  m_vector.get()->Update( 1., *x.unwrappedPointer(), 0. );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Axpy
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Update vector as this = alpha*x + this.
void EpetraVector::axpy( real64 const alpha,
                         EpetraVector const & x )
{
  m_vector.get()->Update( alpha, *x.unwrappedPointer(), 1. );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Axpby
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Update vector as this = alpha*x + beta*this.
void EpetraVector::axpby( real64 const alpha,
                          EpetraVector const & x,
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
void EpetraVector::print( std::ostream & os ) const
{
  if( m_vector )
  {
    m_vector->Print( os );
  }
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Write to matlab-compatible file
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Note: EpetraExt also supports a MatrixMarket format as well
//       if we prefer that.
void EpetraVector::write( string const & filename,
                          bool const mtxFormat ) const
{
  if( mtxFormat )
  {
    // Ensure the ".mtx" extension
    string name( filename );
    if( filename.substr( filename.find_last_of( "." ) + 1 ) != "mtx" )
    {
      name = filename.substr( 0, filename.find_last_of( "." ) ) + ".mtx";
    }
    EpetraExt::MultiVectorToMatrixMarketFile( name.c_str(), *m_vector );
  }
  else
  {
    EpetraExt::MultiVectorToMatlabFile( filename.c_str(), *m_vector );
  }
}

// ----------------------------
// Accessors
// ----------------------------

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get value
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get element globalRow
real64 EpetraVector::get( globalIndex globalRow ) const
{
  return extractLocalVector()[getLocalRowID( globalRow )];
}

void EpetraVector::get( array1d< globalIndex > const & globalIndices,
                        array1d< real64 > & values ) const
{
  real64 const * localVector = extractLocalVector();
  values.resize( globalIndices.size() );

  for( localIndex i = 0; i < globalIndices.size(); ++i )
  {
    values[i] = localVector[getLocalRowID( globalIndices[i] )];
  }
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get unwrapped pointer
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get pointer to raw Epetra object, with const and non-const versions.
Epetra_FEVector const * EpetraVector::unwrappedPointer() const
{
  return m_vector.get();
}

Epetra_FEVector * EpetraVector::unwrappedPointer()
{
  return m_vector.get();
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get the number of global elements
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Return the global size of the vector (total number of elements).
globalIndex EpetraVector::globalSize() const
{
  return m_vector->GlobalLength64();
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get the number of local elements
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Return the local size of the vector (total number of local elements).
localIndex EpetraVector::localSize() const
{
  return m_vector->MyLength();
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// getLocalRowID
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Map a global row index to local row index
localIndex EpetraVector::getLocalRowID( globalIndex const index ) const
{
  return m_vector->Map().LID( integer_conversion< long long >( index ) );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// getGlobalRowID
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Map a local row index to global row index
globalIndex EpetraVector::getGlobalRowID( localIndex const index ) const
{
  return m_vector->Map().GID64( integer_conversion< int >( index ) );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// extractLocalVector
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Extract a view of the local portion of the array
real64 const * EpetraVector::extractLocalVector() const
{
  int dummy;
  double * localVector;
  m_vector->ExtractView( &localVector, &dummy );
  return localVector;
}

real64 * EpetraVector::extractLocalVector()
{
  int dummy;
  double * localVector;
  m_vector->ExtractView( &localVector, &dummy );
  return localVector;
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// ilower
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Returns the index of the first global row owned by that processor.
globalIndex EpetraVector::ilower() const
{
  return m_vector->Map().MinMyGID64();
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// iupper
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// eturns the next index after last global row owned by that processor.
globalIndex EpetraVector::iupper() const
{
  return m_vector->Map().MaxMyGID64() + 1;
}

std::ostream & operator<<( std::ostream & os,
                           EpetraVector const & vec )
{
  vec.print( os );
  return os;
}

} // end geosx

// END_RST_NARRATIVE
