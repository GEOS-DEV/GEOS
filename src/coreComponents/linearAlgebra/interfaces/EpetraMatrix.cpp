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
 * @file EpetraMatrix.cpp
 */

// BEGIN_RST_NARRATIVE EpetraMatrix.rst
// ==============================
// Epetra-based Matrix Object
// ==============================
// This class contains the ParallelMatrix wrappers based on Epetra_Crs Objects.
// The class contains a unique pointer to an Epetra_CrsMatrix as well as constructors,
// functions and accessors for Epetra objects.

// Include the corresponding header file.
#include "EpetraMatrix.hpp"

// Include required Epetra headers
#include <Epetra_Map.h>
#include <Epetra_FECrsGraph.h>
#include <Epetra_FECrsMatrix.h>
#include <Epetra_FEVector.h>
#include <EpetraExt_MatrixMatrix.h>
#include <EpetraExt_RowMatrixOut.h>

#ifdef GEOSX_USE_MPI
#include <Epetra_MpiComm.h>
#else
#include<Epetra_SerialComm.h>
using Epetra_MpiComm =  Epetra_SerialComm;
#endif

// Put everything under the geosx namespace.
namespace geosx
{

// ----------------------------
// Constructors
// ----------------------------

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Create an empty matrix (meant to be used for declaration)
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
EpetraMatrix::EpetraMatrix()
: Base()
{ }

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Copy constructor
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
EpetraMatrix::EpetraMatrix( EpetraMatrix const & src )
: Base()
{
  GEOSX_ASSERT( !src.isOpen() );
  GEOSX_ASSERT( src.unwrappedPointer() != nullptr );

  m_matrix = std::make_unique< Epetra_FECrsMatrix >( *src.unwrappedPointer() );
  m_src_map = std::make_unique< Epetra_Map >( m_matrix->DomainMap() );
  m_dst_map = std::make_unique< Epetra_Map >( m_matrix->RangeMap() );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Destructor
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
EpetraMatrix::~EpetraMatrix() = default;

// -----------------------------
// Create
// -----------------------------
// Allocate matrix (prepare to be filled with data).

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Create a matrix from number of elements
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
void EpetraMatrix::createWithGlobalSize( globalIndex const globalSize,
                                         localIndex const maxEntriesPerRow,
                                         MPI_Comm const & comm )
{
  createWithGlobalSize( globalSize, globalSize, maxEntriesPerRow, comm ); // just call general version
}

void EpetraMatrix::createWithGlobalSize( globalIndex const globalRows,
                                         globalIndex const globalCols,
                                         localIndex const maxEntriesPerRow,
                                         MPI_Comm const & MPI_PARAM(comm) )
{
  m_dst_map = std::make_unique< Epetra_Map >( globalRows,
                                              0,
                                              Epetra_MpiComm( MPI_PARAM(comm) ) );
  m_src_map = std::make_unique< Epetra_Map >( globalCols,
                                              0,
                                              Epetra_MpiComm( MPI_PARAM(comm) ) );
  m_matrix = std::make_unique< Epetra_FECrsMatrix >( Copy,
                                                     *m_dst_map,
                                                     integer_conversion< int >( maxEntriesPerRow ),
                                                     false );
}

void EpetraMatrix::createWithLocalSize( localIndex const localSize,
                                        localIndex const maxEntriesPerRow,
                                        MPI_Comm const & comm )
{
  createWithLocalSize( localSize, localSize, maxEntriesPerRow, comm ); // just call general version
}

void EpetraMatrix::createWithLocalSize( localIndex const localRows,
                                        localIndex const localCols,
                                        localIndex const maxEntriesPerRow,
                                        MPI_Comm const & MPI_PARAM(comm) )
{
  m_dst_map = std::make_unique< Epetra_Map >( integer_conversion< globalIndex >( -1 ),
                                              integer_conversion< int >( localRows ),
                                              0,
                                              Epetra_MpiComm( MPI_PARAM(comm) ) );
  m_src_map = std::make_unique< Epetra_Map >( integer_conversion< globalIndex >( -1 ),
                                              integer_conversion< int >( localCols ),
                                              0,
                                              Epetra_MpiComm( MPI_PARAM(comm) ) );
  m_matrix = std::make_unique< Epetra_FECrsMatrix >( Copy,
                                                     *m_dst_map,
                                                     integer_conversion< int >( maxEntriesPerRow ),
                                                     false );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Reinitialize.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Keeps the map and graph but sets all values to user-defined value.
void EpetraMatrix::set( real64 const value )
{
  GEOSX_ASSERT( !isOpen() );
  m_matrix->PutScalar( value );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Reinitialize.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Keeps the map and graph but sets all values to 0.
void EpetraMatrix::zero()
{
  GEOSX_ASSERT( !isOpen() );
  m_matrix->PutScalar( 0 );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Open
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Empty open function (implemented for HYPRE compatibility).
void EpetraMatrix::open()
{
  GEOSX_ASSERT( !isOpen() );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Close
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Fix the sparsity pattern, make the data contiguous in memory and optimize storage.
void EpetraMatrix::close()
{
  GEOSX_ASSERT( isOpen() );
  m_matrix->GlobalAssemble( *m_src_map, *m_dst_map );
  m_assembled = true;
}

// -------------------------
// Add/Set
// -------------------------

// 1x1
void EpetraMatrix::add( globalIndex const rowIndex,
                        globalIndex const colIndex,
                        real64 const value )
{
  GEOSX_ASSERT( isOpen() );
  m_matrix->SumIntoGlobalValues( rowIndex, 1, &value, &colIndex );
}

void EpetraMatrix::set( globalIndex const rowIndex,
                        globalIndex const colIndex,
                        real64 const value )
{
  GEOSX_ASSERT( isOpen() );
  m_matrix->ReplaceGlobalValues( rowIndex, 1, &value, &colIndex );
}

void EpetraMatrix::insert( globalIndex const rowIndex,
                           globalIndex const colIndex,
                           real64 const value )
{
  GEOSX_ASSERT( isOpen() );
  m_matrix->InsertGlobalValues( rowIndex, 1, &value, &colIndex );
}

// 1xN c-style
void EpetraMatrix::add( globalIndex const rowIndex,
                        globalIndex const * colIndices,
                        real64 const * values,
                        localIndex size )
{
  GEOSX_ASSERT( isOpen() );
  m_matrix->SumIntoGlobalValues( rowIndex, integer_conversion< int >( size ), values, colIndices );
}

void EpetraMatrix::set( globalIndex const rowIndex,
                        globalIndex const * colIndices,
                        real64 const * values,
                        localIndex size )
{
  GEOSX_ASSERT( isOpen() );
  m_matrix->ReplaceGlobalValues( rowIndex, integer_conversion< int >( size ), values, colIndices );
}

void EpetraMatrix::insert( globalIndex const rowIndex,
                           globalIndex const * colIndices,
                           real64 const * values,
                           localIndex size )
{
  GEOSX_ASSERT( isOpen() );
  m_matrix->InsertGlobalValues( rowIndex, integer_conversion< int >( size ), values, colIndices );
}

// 1xN array1d style
void EpetraMatrix::add( globalIndex const rowIndex,
                        array1d< globalIndex > const & colIndices,
                        array1d< real64 > const & values )
{
  GEOSX_ASSERT( isOpen() );
  m_matrix->SumIntoGlobalValues( rowIndex, integer_conversion< int >( colIndices.size() ), values.data(),
                                 colIndices.data() );
}

void EpetraMatrix::set( globalIndex const rowIndex,
                        array1d< globalIndex > const & colIndices,
                        array1d< real64 > const & values )
{
  GEOSX_ASSERT( isOpen() );
  m_matrix->ReplaceGlobalValues( rowIndex, integer_conversion< int >( colIndices.size() ), values.data(),
                                 colIndices.data() );
}

void EpetraMatrix::insert( globalIndex const rowIndex,
                           array1d< globalIndex > const & colIndices,
                           array1d< real64 > const & values )
{
  GEOSX_ASSERT( isOpen() );
  m_matrix->InsertGlobalValues( rowIndex, integer_conversion< int >( colIndices.size() ), values.data(),
                                colIndices.data() );
}

// MxN array2d style
void EpetraMatrix::add( array1d< globalIndex > const & rowIndices,
                        array1d< globalIndex > const & colIndices,
                        array2d< real64 > const & values )
{
  GEOSX_ASSERT( isOpen() );
  m_matrix->SumIntoGlobalValues( integer_conversion< int >( rowIndices.size() ), rowIndices.data(),
                                 integer_conversion< int >( colIndices.size() ), colIndices.data(),
                                 values.data(), Epetra_FECrsMatrix::ROW_MAJOR );
}

void EpetraMatrix::set( array1d< globalIndex > const & rowIndices,
                        array1d< globalIndex > const & colIndices,
                        array2d< real64 > const & values )
{
  GEOSX_ASSERT( isOpen() );
  m_matrix->ReplaceGlobalValues( integer_conversion< int >( rowIndices.size() ), rowIndices.data(),
                                 integer_conversion< int >( colIndices.size() ), colIndices.data(),
                                 values.data(), Epetra_FECrsMatrix::ROW_MAJOR );
}

void EpetraMatrix::insert( array1d< globalIndex > const & rowIndices,
                           array1d< globalIndex > const & colIndices,
                           array2d< real64 > const & values )
{
  GEOSX_ASSERT( isOpen() );
  m_matrix->InsertGlobalValues( integer_conversion< int >( rowIndices.size() ), rowIndices.data(),
                                integer_conversion< int >( colIndices.size() ), colIndices.data(),
                                values.data(), Epetra_FECrsMatrix::ROW_MAJOR );
}

void EpetraMatrix::add( globalIndex const * rowIndices,
                        globalIndex const * colIndices,
                        real64 const * values,
                        localIndex const numRows,
                        localIndex const numCols )
{
  GEOSX_ASSERT( isOpen() );
  m_matrix->SumIntoGlobalValues( integer_conversion< int >( numRows ), rowIndices,
                                 integer_conversion< int >( numCols ), colIndices,
                                 values, Epetra_FECrsMatrix::ROW_MAJOR );
}

void EpetraMatrix::set( globalIndex const * rowIndices,
                        globalIndex const * colIndices,
                        real64 const * values,
                        localIndex const numRows,
                        localIndex const numCols )
{
  GEOSX_ASSERT( isOpen() );
  m_matrix->ReplaceGlobalValues( integer_conversion< int >( numRows ), rowIndices,
                                 integer_conversion< int >( numCols ), colIndices,
                                 values, Epetra_FECrsMatrix::ROW_MAJOR );
}

void EpetraMatrix::insert( globalIndex const * rowIndices,
                           globalIndex const * colIndices,
                           real64 const * values,
                           localIndex const numRows,
                           localIndex const numCols )
{
  GEOSX_ASSERT( isOpen() );
  m_matrix->InsertGlobalValues( integer_conversion< int >( numRows ), rowIndices,
                                integer_conversion< int >( numCols ), colIndices,
                                values, Epetra_FECrsMatrix::ROW_MAJOR );
}

// -------------------------
// Linear Algebra
// -------------------------

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Matrix/vector multiplication
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Perform the matrix-vector product this * src = dst.
void EpetraMatrix::multiply( EpetraVector const & src,
                             EpetraVector & dst ) const
{
  GEOSX_ASSERT( !isOpen() );
  m_matrix->Multiply( false, *src.unwrappedPointer(), *dst.unwrappedPointer() );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Matrix/matrix multiplication
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Perform the matrix-matrix product this * src = dst.
void EpetraMatrix::multiply( EpetraMatrix const & src,
                             EpetraMatrix & dst,
                             bool const closeResult ) const
{
  this->multiply( false, src, false, dst, closeResult );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Matrix/matrix multiplication
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Perform the matrix-matrix product this^T * src = dst.
void EpetraMatrix::leftMultiplyTranspose( EpetraMatrix const & src,
                                          EpetraMatrix & dst,
                                          bool const closeResult ) const
{
  this->multiply( true, src, false, dst, closeResult );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Matrix/matrix multiplication
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Perform the matrix-matrix product src * this^T = dst.
void EpetraMatrix::rightMultiplyTranspose( EpetraMatrix const & src,
                                           EpetraMatrix & dst,
                                           bool const closeResult ) const
{
  GEOSX_ASSERT_MSG( m_matrix, "Matrix has not been created" );
  GEOSX_ASSERT( !isOpen() );
  src.multiply( false, *this, true, dst, closeResult );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Generalized matrix/vector product.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Compute gemv <tt>y = alpha*A*x + beta*y</tt>.
void EpetraMatrix::gemv( real64 const alpha,
                         EpetraVector const & x,
                         real64 const beta,
                         EpetraVector & y,
                         bool useTranspose )
{
  GEOSX_ASSERT_MSG( m_matrix, "Matrix has not been created" );
  GEOSX_ASSERT( !isOpen() );

  EpetraVector Ax( y );
  m_matrix->Multiply( useTranspose, *x.unwrappedPointer(), *Ax.unwrappedPointer() );
  y.axpby( alpha, Ax, beta );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Scale.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Multiply all elements by scalingFactor.
void EpetraMatrix::scale( real64 const scalingFactor )
{
  GEOSX_ASSERT_MSG( m_matrix, "Matrix has not been created" );
  GEOSX_ASSERT( !isOpen() );
  m_matrix->Scale( scalingFactor );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Left and right scaling
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
void EpetraMatrix::leftScale( EpetraVector const & vec )
{
  GEOSX_ASSERT_MSG( m_matrix, "Matrix has not been created" );
  GEOSX_ASSERT( !isOpen() );
  m_matrix->LeftScale( *( *vec.unwrappedPointer() )( 0 ) );
}

void EpetraMatrix::rightScale( EpetraVector const & vec )
{
  GEOSX_ASSERT_MSG( m_matrix, "Matrix has not been created" );
  GEOSX_ASSERT( !isOpen() );
  m_matrix->RightScale( *( *vec.unwrappedPointer() )( 0 ) );
}

void EpetraMatrix::leftRightScale( EpetraVector const & vecLeft,
                                   EpetraVector const & vecRight )
{
  leftScale( vecLeft );
  rightScale( vecRight );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get global row copy
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// The challenge here is that columns are stored with local, not global,
// indices, so we need to do conversions back and forth
void EpetraMatrix::getRowCopy( globalIndex globalRow,
                               array1d< globalIndex > & colIndices,
                               array1d< real64 > & values ) const
{
  GEOSX_ASSERT_MSG( m_matrix, "Matrix has not been created" );
  GEOSX_ASSERT_GE( globalRow, ilower());
  GEOSX_ASSERT_GT( iupper(), globalRow );
  GEOSX_ASSERT_MSG( m_matrix->IndicesAreLocal(), "Matrix must be assembled by calling close() at least once" );

  int n_entries;
  int * indices_ptr;
  double * values_ptr;

  int const err = m_matrix->ExtractMyRowView( m_matrix->LRID( globalRow ), n_entries, values_ptr, indices_ptr );
  GEOSX_ERROR_IF( err != 0,
                 "getRowCopy failed. This often happens if the requested global row "
                 "is not local to this processor, or if close() hasn't been called." );

  localIndex const length = integer_conversion< localIndex >( n_entries );
  values.resize( length );
  colIndices.resize( length );

  for( localIndex i = 0; i < length; ++i )
  {
    colIndices[i] = m_matrix->GCID64( indices_ptr[i] );
    values[i] = values_ptr[i];
  }
}

real64 EpetraMatrix::getDiagValue( globalIndex globalRow ) const
{
  GEOSX_ASSERT_MSG( m_matrix, "Matrix has not been created" );
  GEOSX_ASSERT_GE( globalRow, ilower());
  GEOSX_ASSERT_GT( iupper(), globalRow );
  GEOSX_ASSERT_MSG( m_matrix->IndicesAreLocal(), "Matrix must be assembled by calling close() at least once" );

  double * values = nullptr;
  int * indices = nullptr;
  int length;

  int const err = m_matrix->ExtractMyRowView( m_matrix->LRID( globalRow ), length, values, indices );
  GEOSX_ERROR_IF_NE_MSG( err, 0, "EpetraCrsMatrix::ExtractGlobalRowView failed" );

  for( int j = 0; j < length; ++j )
  {
    if( m_matrix->GCID64( indices[j] ) == globalRow )
    {
      return values[j];
    }
  }

  return 0.0;
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Clear row
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Clear the row.  By default the diagonal value will be set
// to zero, but the user can pass a desired diagValue.
void EpetraMatrix::clearRow( globalIndex const globalRow,
                             real64 const diagValue )
{
  GEOSX_ASSERT_MSG( m_matrix, "Matrix has not been created" );
  GEOSX_ASSERT_GE( globalRow, ilower());
  GEOSX_ASSERT_GT( iupper(), globalRow );
  GEOSX_ASSERT_MSG( m_matrix->IndicesAreLocal(), "Matrix must be assembled by calling close() at least once" );

  double * values = nullptr;
  int length;

  int const err = m_matrix->ExtractMyRowView( m_matrix->LRID( globalRow ), length, values );
  GEOSX_ERROR_IF_NE_MSG( err, 0, "EpetraCrsMatrix::ExtractGlobalRowView failed" );

  for( int j = 0; j < length; ++j )
  {
    values[j] = 0.0;
  }

  set( globalRow, globalRow, diagValue );
}

// ---------------------------------------------------------
//  Accessors
// ---------------------------------------------------------

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get pointer.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Accessor for the pointer to the raw Epetra matrix
Epetra_FECrsMatrix * EpetraMatrix::unwrappedPointer() const
{
  GEOSX_ASSERT_MSG( m_matrix, "Matrix has not been created" );
  return m_matrix.get();
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get number of global rows.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Accessor for the number of global rows
globalIndex EpetraMatrix::globalRows() const
{
  GEOSX_ASSERT_MSG( m_matrix, "Matrix has not been created" );
  GEOSX_ASSERT( !isOpen() );
  return m_matrix->NumGlobalRows64();
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get number of global columns.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Accessor for the number of global columns
globalIndex EpetraMatrix::globalCols() const
{
  GEOSX_ASSERT_MSG( m_matrix, "Matrix has not been created" );
  GEOSX_ASSERT( !isOpen() );
  return m_matrix->NumGlobalCols64();
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get the lower index owned by processor.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Accessor for the index of the first global row
globalIndex EpetraMatrix::ilower() const
{
  GEOSX_ASSERT_MSG( m_matrix, "Matrix has not been created" );
  GEOSX_ASSERT( !isOpen() );
  return m_matrix->RowMap().MinMyGID64();
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get the next index after upper index owned by processor.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Accessor for the index of the last global row
globalIndex EpetraMatrix::iupper() const
{
  GEOSX_ASSERT_MSG( m_matrix, "Matrix has not been created" );
  GEOSX_ASSERT( !isOpen() );
  return m_matrix->RowMap().MaxMyGID64() + 1;
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get number of local nonzeros.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Accessor for the number of local nonzeros
localIndex EpetraMatrix::localNonzeros() const
{
  GEOSX_ASSERT_MSG( m_matrix, "Matrix has not been created" );
  GEOSX_ASSERT( !isOpen() );
  return m_matrix->NumMyNonzeros();
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get number of global nonzeros.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Accessor for the number of global nonzeros
globalIndex EpetraMatrix::globalNonzeros() const
{
  GEOSX_ASSERT_MSG( m_matrix, "Matrix has not been created" );
  GEOSX_ASSERT( !isOpen() );
  return m_matrix->NumGlobalNonzeros64();
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Print to terminal.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Wrapper to print the trilinos output of the matrix
void EpetraMatrix::print( std::ostream & os ) const
{
  GEOSX_ASSERT( !isOpen() );
  if( m_matrix )
  {
    m_matrix->Print( os );
  }
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Write to matlab-compatible file
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Note: EpetraExt also supports a MatrixMarket format as well
void EpetraMatrix::write( string const & filename,
                          bool const mtxFormat ) const
{
  GEOSX_ASSERT( !isOpen() );
  if( mtxFormat )
  {
    // Ensure the ".mtx" extension
    string name( filename );
    if( filename.substr( filename.find_last_of( "." ) + 1 ) != "mtx" )
    {
      name = filename.substr( 0, filename.find_last_of( "." ) ) + ".mtx";
    }
    EpetraExt::RowMatrixToMatrixMarketFile( name.c_str(), *m_matrix );
  }
  else
  {
    EpetraExt::RowMatrixToMatlabFile( filename.c_str(), *m_matrix );
  }
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Inf-norm.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Returns the infinity norm of the matrix.
real64 EpetraMatrix::normInf() const
{
  GEOSX_ASSERT_MSG( m_matrix, "Matrix has not been created" );
  GEOSX_ASSERT( !isOpen() );
  return m_matrix->NormInf();
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// 1-norm.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Returns the one norm of the matrix.
real64 EpetraMatrix::norm1() const
{
  GEOSX_ASSERT_MSG( m_matrix, "Matrix has not been created" );
  GEOSX_ASSERT( !isOpen() );
  return m_matrix->NormOne();
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Frobenius-norm.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Returns the Frobenius norm of the matrix.
real64 EpetraMatrix::normFrobenius() const
{
  GEOSX_ASSERT_MSG( m_matrix, "Matrix has not been created" );
  GEOSX_ASSERT( !isOpen() );
  return m_matrix->NormFrobenius();
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Is-assembled.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Boolean indicator. True = matrix assembled and ready to be used.
bool EpetraMatrix::isAssembled() const
{
  return m_assembled;
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// MatrixMatrixMultiply
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Perform the matrix-matrix product A*B = C.
void EpetraMatrix::multiply( bool const transA,
                             EpetraMatrix const & B,
                             bool const transB,
                             EpetraMatrix & C,
                             bool const closeResult ) const
{
  GEOSX_ASSERT_MSG( m_matrix, "Matrix has not been created" );
  GEOSX_ASSERT( !isOpen() );
  GEOSX_ASSERT_MSG( B.m_matrix, "Matrix has not been created" );
  GEOSX_ASSERT( !B.isOpen() );

  if( !C.m_matrix )
  {
    C.createWithLocalSize( transA ? localCols() : localRows(),
                           transB ? B.localRows() : B.localCols(),
                           1,
                           comm() )
  }

  int const err = EpetraExt::MatrixMatrix::Multiply( *m_matrix,
                                                     transA,
                                                     *B.unwrappedPointer(),
                                                     transB,
                                                     *C.unwrappedPointer(),
                                                     closeResult );

  GEOSX_ERROR_IF_NE_MSG( err, 0, "Error in EpetraExt::MatrixMatrix::Multiply()" );
  if( closeResult )
  {
    C.m_open = false;
  }
}



// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// getLocalRowID
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Map a global row index to local row index
localIndex EpetraMatrix::getLocalRowID( globalIndex const index ) const
{
  GEOSX_ASSERT_MSG( m_matrix, "Matrix has not been created" );
  return m_matrix->LRID( index );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// getGlobalRowID
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Map a local row index to global row index
globalIndex EpetraMatrix::getGlobalRowID( localIndex const index ) const
{
  GEOSX_ASSERT_MSG( m_matrix, "Matrix has not been created" );
  return m_matrix->GRID64( integer_conversion< int >( index ) );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// localCols
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Return the local number of columns on each processor
// NOTE: direct use of NumMyCols() counts also for overlays. To avoid those, DomainMap() is needed
localIndex EpetraMatrix::localCols() const
{
  GEOSX_ASSERT_MSG( m_matrix, "Matrix has not been created" );
  return m_matrix->DomainMap().NumMyElements();
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// localRows
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Return the local number of columns on each processor
localIndex EpetraMatrix::localRows() const
{
  GEOSX_ASSERT_MSG( m_matrix, "Matrix has not been created" );
  return m_matrix->RowMap().NumMyElements();
}

std::ostream & operator<<( std::ostream & os,
                           EpetraMatrix const & matrix )
{
  matrix.print( os );
  return os;
}

} // end geosx namespace

// END_RST_NARRATIVE

/* TODO: We should make a decision about thread safety in another
 * pull request.  Either we make Epetra threadsafe or we move to
 * Tpetra as an alternative.
 */

/* SCRATCH CODE - possible template for threadsafe assembly */
/* DELETE WHEN NO LONGER NEEDED */

/*
   // """""""""""""""""""""""""""""""""""""""""""""""""""""""""""
   // Ignore the following. Thread-safe prototype, commented out.
   // """""""""""""""""""""""""""""""""""""""""""""""""""""""""""
   //  template<typename int_type>
   //  int Epetra_CrsMatrix::TSumIntoGlobalValues(int_type Row,
   //              int NumEntries,
   //              const double * srcValues,
   //              const int_type *Indices)
   //  {
   int j;
   int ierr = 0;
   int Loc = 0;


   int locRow = Graph_.LRID( Row ); // Normalize row range

   if( locRow < 0 || locRow >= NumMyRows_ )
   {
    EPETRA_CHK_ERR( -1 ); // Not in Row range
   }

   if( StaticGraph() && !Graph_.HaveColMap())
   {
    EPETRA_CHK_ERR( -1 );
   }

   double * RowValues = Values( locRow );

   if( !StaticGraph())
   {
    for( j=0 ; j<NumEntries ; j++ )
    {
      int_type Index = Indices[j];
      if( Graph_.FindglobalIndexLoc( locRow, Index, j, Loc ))
   //  #ifdef EPETRA_HAVE_OMP
   //  #ifdef EPETRA_HAVE_OMP_NONASSOCIATIVE
   //  #pragma omp atomic
   //  #endif
   //  #endif
   //          RowValues[Loc] += srcValues[j];
        RAJA::atomicAdd<ATOMIC_POL2>( &RowValues [Loc], srcValues[j] );
      else
        ierr = 2;   // Value Excluded
    }
   }
   else
   {
    const Epetra_BlockMap& colmap = Graph_.ColMap();
    int NumColIndices = Graph_.NumMyIndices( locRow );
    const int* ColIndices = Graph_.Indices( locRow );

    if( Graph_.Sorted())
    {
      int insertPoint;
      for( j=0 ; j<NumEntries ; j++ )
      {
        int Index = colmap.LID( Indices[j] );

        // Check whether the next added element is the subsequent element in
        // the graph indices, then we can skip the binary search
        if( Loc < NumColIndices && Index == ColIndices[Loc] )
   //  #ifdef EPETRA_HAVE_OMP
   //  #ifdef EPETRA_HAVE_OMP_NONASSOCIATIVE
   //  #pragma omp atomic
   //  #endif
   //  #endif
   //            RowValues[Loc] += srcValues[j];
          RAJA::atomicAdd<ATOMIC_POL2>( &RowValues [Loc], srcValues[j] );
        else
        {
          Loc = Epetra_Util_binary_search( Index, ColIndices, NumColIndices, insertPoint );
          if( Loc > -1 )
   //  #ifdef EPETRA_HAVE_OMP
   //  #ifdef EPETRA_HAVE_OMP_NONASSOCIATIVE
   //  #pragma omp atomic
   //  #endif
   //  #endif
   //              RowValues[Loc] += srcValues[j];
            RAJA::atomicAdd<ATOMIC_POL2>( &RowValues [Loc], srcValues[j] );
          else
            ierr = 2;   // Value Excluded
        }
 ++Loc;
      }
    }
    else
      for( j=0 ; j<NumEntries ; j++ )
      {
        int Index = colmap.LID( Indices[j] );
        if( Graph_.FindMyIndexLoc( NumColIndices, ColIndices, Index, j, Loc ))
   //  #ifdef EPETRA_HAVE_OMP
   //  #ifdef EPETRA_HAVE_OMP_NONASSOCIATIVE
   //  #pragma omp atomic
   //  #endif
   //  #endif
   //            RowValues[Loc] += srcValues[j];
          RAJA::atomicAdd<ATOMIC_POL2>( &RowValues [Loc], srcValues[j] );
        else
          ierr = 2;   // Value Excluded
      }
   }

   NormOne_ = -1.0;   // Reset Norm so it will be recomputed.
   NormInf_ = -1.0;   // Reset Norm so it will be recomputed.
   NormFrob_ = -1.0;

   EPETRA_CHK_ERR( ierr );

   return(0);
   //  }

 */
