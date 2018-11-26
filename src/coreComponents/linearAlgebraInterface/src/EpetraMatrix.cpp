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
{}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Copy constructor
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""

EpetraMatrix::EpetraMatrix( EpetraMatrix const &src )
{
  GEOS_ERROR_IF(src.unwrappedPointer() == nullptr, "Input matrix looks empty" );
  //TODO GEOS_ERROR_IF( !src.isClosed(), "Input matrix hasn't been properly closed before copy");
  
  m_matrix  = std::unique_ptr<Epetra_FECrsMatrix>( new Epetra_FECrsMatrix( *src.unwrappedPointer() ) );
  m_src_map = std::unique_ptr<Epetra_Map>(new Epetra_Map( m_matrix->DomainMap()));
  m_dst_map = std::unique_ptr<Epetra_Map>(new Epetra_Map( m_matrix->RangeMap()));
}

// -----------------------------
// Create
// -----------------------------
// Allocate matrix (prepare to be filled with data).

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Create a matrix from an Epetra_FECrsGraph.
// """""""""""""""""""""""""""""""""""""""""""""""

void EpetraMatrix::create( Epetra_FECrsGraph const &graph )
{
  m_matrix  = std::unique_ptr<Epetra_FECrsMatrix>( new Epetra_FECrsMatrix( Copy, graph ) );
  m_src_map = std::unique_ptr<Epetra_Map>(new Epetra_Map( m_matrix->DomainMap()));
  m_dst_map = std::unique_ptr<Epetra_Map>(new Epetra_Map( m_matrix->RangeMap()));
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Create a matrix from number of elements
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""

void EpetraMatrix::createWithGlobalSize( trilinosTypes::gid const globalSize,
                                         trilinosTypes::lid const maxEntriesPerRow,
                                         MPI_Comm const & comm )
{
  createWithGlobalSize(globalSize,globalSize,maxEntriesPerRow,comm); // just call general version
}


void EpetraMatrix::createWithGlobalSize( trilinosTypes::gid const globalRows,
                                         trilinosTypes::gid const globalCols,
                                         trilinosTypes::lid const maxEntriesPerRow,
                                         MPI_Comm const & comm )
{
  m_dst_map = std::unique_ptr<Epetra_Map>(new Epetra_Map( globalRows, 0, Epetra_MpiComm( comm ) ));
  m_src_map = std::unique_ptr<Epetra_Map>(new Epetra_Map( globalCols, 0, Epetra_MpiComm( comm ) ));
  m_matrix = std::unique_ptr<Epetra_FECrsMatrix>( new Epetra_FECrsMatrix( Copy, *m_dst_map, maxEntriesPerRow, false ) );
}


void EpetraMatrix::createWithLocalSize( trilinosTypes::lid const localSize,
                                        trilinosTypes::lid const maxEntriesPerRow,
                                        MPI_Comm const & comm )
{
  createWithLocalSize(localSize,localSize,maxEntriesPerRow,comm); // just call general version
}


void EpetraMatrix::createWithLocalSize( trilinosTypes::lid const localRows,
                                        trilinosTypes::lid const localCols,
                                        trilinosTypes::lid const maxEntriesPerRow,
                                        MPI_Comm const & comm )
{
  m_dst_map = std::unique_ptr<Epetra_Map>(new Epetra_Map( -1, localRows, 0, Epetra_MpiComm( comm ) ));
  m_src_map = std::unique_ptr<Epetra_Map>(new Epetra_Map( -1, localCols, 0, Epetra_MpiComm( comm ) ));
  m_matrix = std::unique_ptr<Epetra_FECrsMatrix>( new Epetra_FECrsMatrix( Copy, *m_dst_map, maxEntriesPerRow, false ) );
}


// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Reinitialize.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Keeps the map and graph but sets all values to 0.

void EpetraMatrix::zero()
{
  m_matrix->PutScalar( 0 );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Open
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Empty open function (implemented for HYPRE compatibility).

void EpetraMatrix::open()
{}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Close
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Fix the sparsity pattern, make the data contiguous in memory and optimize storage.

void EpetraMatrix::close()
{

  m_matrix->GlobalAssemble(*m_src_map,*m_dst_map);
  assembled = true;
}

// -------------------------
// Add/Set
// -------------------------

// 1x1 

void EpetraMatrix::add( trilinosTypes::gid const rowIndex,
                        trilinosTypes::gid const colIndex,
                        real64 const value )
{
  m_matrix->SumIntoGlobalValues( rowIndex, 1, &value, &colIndex );
}


void EpetraMatrix::set( trilinosTypes::gid const rowIndex,
                        trilinosTypes::gid const colIndex,
                        real64 const value )
{
  m_matrix->ReplaceGlobalValues( rowIndex, 1, &value, &colIndex );
}


void EpetraMatrix::insert( trilinosTypes::gid const rowIndex,
                           trilinosTypes::gid const colIndex,
                           real64 const value )
{
  m_matrix->InsertGlobalValues( rowIndex, 1, &value, &colIndex );
}

// 1xN c-style

void EpetraMatrix::add( trilinosTypes::gid const rowIndex,
                        trilinosTypes::gid const * colIndices,
                        real64 const * values,
                        trilinosTypes::lid size )
{
  m_matrix->SumIntoGlobalValues( rowIndex, size, values, colIndices );
}

void EpetraMatrix::set( trilinosTypes::gid const rowIndex,
                        trilinosTypes::gid const * colIndices,
                        real64 const * values,
                        trilinosTypes::lid size )
{
  m_matrix->ReplaceGlobalValues( rowIndex, size, values, colIndices );
}


void EpetraMatrix::insert( trilinosTypes::gid const rowIndex,
                           trilinosTypes::gid const * colIndices,
                           real64 const * values,
                           trilinosTypes::lid size )
{
  m_matrix->InsertGlobalValues( rowIndex, size, values, colIndices );
}

// 1xN array1d style

void EpetraMatrix::add( trilinosTypes::gid const rowIndex,
                        array1d<trilinosTypes::gid> const &colIndices,
                        array1d<real64> const &values )
{
  // TODO: add integer_conversion
  m_matrix->SumIntoGlobalValues( rowIndex, colIndices.size(), values.data(), colIndices.data() );
}

void EpetraMatrix::set( trilinosTypes::gid const rowIndex,
                        array1d<trilinosTypes::gid> const &colIndices,
                        array1d<real64> const &values )
{
  // TODO: add integer_conversion
  m_matrix->ReplaceGlobalValues( rowIndex, colIndices.size(), values.data(), colIndices.data() );
}

void EpetraMatrix::insert( trilinosTypes::gid const rowIndex,
                           array1d<trilinosTypes::gid> const &colIndices,
                           array1d<real64> const &values )
{
  // TODO: add integer_conversion
  m_matrix->InsertGlobalValues( rowIndex, colIndices.size(), values.data(), colIndices.data() );
}

// MxN array2d style

void EpetraMatrix::add( array1d<trilinosTypes::gid> const & rowIndices,
                        array1d<trilinosTypes::gid> const & colIndices,
                        array2d<real64> const & values )
{
  // TODO: add integer_conversion
  m_matrix->SumIntoGlobalValues( rowIndices.size(), rowIndices.data(), colIndices.size(), colIndices.data(), values.data(), Epetra_FECrsMatrix::ROW_MAJOR );
}

void EpetraMatrix::set( array1d<trilinosTypes::gid> const & rowIndices,
                        array1d<trilinosTypes::gid> const & colIndices,
                        array2d<real64> const & values )
{
  // TODO: add integer_conversion
  m_matrix->ReplaceGlobalValues( rowIndices.size(), rowIndices.data(), colIndices.size(), colIndices.data(), values.data(), Epetra_FECrsMatrix::ROW_MAJOR );
}

void EpetraMatrix::insert( array1d<trilinosTypes::gid> const & rowIndices,
                           array1d<trilinosTypes::gid> const & colIndices,
                           array2d<real64> const & values )
{
  // TODO: add integer_conversion
  m_matrix->InsertGlobalValues( rowIndices.size(), rowIndices.data(), colIndices.size(), colIndices.data(), values.data(), Epetra_FECrsMatrix::ROW_MAJOR );
}


// -------------------------
// Linear Algebra
// -------------------------

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Matrix/vector multiplication
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Perform the matrix-vector product A*src = dst.

void EpetraMatrix::multiply( EpetraVector const &src,
                             EpetraVector &dst ) const
{
  m_matrix->Multiply( false, *src.unwrappedPointer(), *dst.unwrappedPointer() );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Compute residual.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Compute res = b - Ax (residual form).

void EpetraMatrix::residual( EpetraVector const &x,
                             EpetraVector const &b,
                             EpetraVector &r ) const
{
  m_matrix->Multiply( false, *x.unwrappedPointer(), *r.unwrappedPointer() );
  r.axpby( 1.0, b, -1.0 );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Generalized matrix/vector product.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Compute gemv <tt>y = alpha*A*x + beta*y</tt>.

void EpetraMatrix::gemv( real64 const alpha,
                         EpetraVector const &x,
                         real64 const beta,
                         EpetraVector &y,
                         bool useTranspose )
{
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
  m_matrix->Scale( scalingFactor );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Left and right scaling
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""

void EpetraMatrix::leftScale( EpetraVector const &vec )
{
  m_matrix->LeftScale( *(*vec.unwrappedPointer())(0) );
}

void EpetraMatrix::rightScale( EpetraVector const &vec )
{
  m_matrix->RightScale( *(*vec.unwrappedPointer())(0) );
}

void EpetraMatrix::leftRightScale( EpetraVector const &vecLeft,
                                   EpetraVector const &vecRight )
{
  leftScale(vecLeft);
  rightScale(vecRight);
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Clear row.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""

void EpetraMatrix::clearRow( trilinosTypes::gid const row,
                             real64 const factor )
{
  // Get the local index corresponding to the global index input. Returns -1
  // if the global index is not owned by this processor.
  int local_row = m_matrix->LRID( row );

  // If this processor owns the input row:
  if( local_row >= 0 )
  {
    double *values = nullptr;
    int *col_indices = nullptr;
    int num_entries;

    // Extract the row and get the number of entries, the values and the column indices.
    m_matrix->ExtractMyRowView( local_row, num_entries, values, col_indices );

    // TODO this if may be unnecessary?
    if( values != nullptr && col_indices != nullptr && num_entries > 0 )
    {
      int* diag_find = std::find( col_indices, col_indices+num_entries-1, local_row );
      long int diag_index = (diag_find - col_indices);

      // Set all non-diagonal entries to 0
      for( int j=0 ; j<num_entries ; ++j )
      {
        if( diag_index != j )
        {
          values[j] = 0.;
        }
      }
      // Scale the diagonal (existing) value with factor.
      values[diag_index] *= factor;
    }
  }
}


// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get global row.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Extract the global row and output the number of elements, values and column indices.

void EpetraMatrix::getRow( trilinosTypes::gid GlobalRow,
                           trilinosTypes::lid &NumEntries,
                           real64* Values,
                           trilinosTypes::gid* Indices ) const
{
  m_matrix->ExtractGlobalRowView( GlobalRow, NumEntries, Values, Indices );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get local row.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Extract the local row and output the number of elements, values and column indices.

void EpetraMatrix::getLocalRow( trilinosTypes::lid localRow,
                                trilinosTypes::lid & NumEntries,
                                real64 * & Values,
                                trilinosTypes::lid * & Indices ) const
{
  m_matrix->ExtractMyRowView( localRow, NumEntries, Values, Indices );
}

void EpetraMatrix::getRow( trilinosTypes::gid GlobalRow,
                                 trilinosTypes::lid &NumEntries,
                                 std::vector<real64> &vecValues,
                                 std::vector<trilinosTypes::gid> &vecIndices ) const
{
  real64* Values;
  trilinosTypes::gid* Indices;
  m_matrix->ExtractGlobalRowView( GlobalRow, NumEntries, Values, Indices );
  if( m_matrix->MyGRID( GlobalRow ))
  {
    vecIndices.assign( Indices, Indices+NumEntries );
    vecValues.assign( Values, Values+NumEntries );
  }
}

void EpetraMatrix::getLocalRow( trilinosTypes::lid localRow,
                                      trilinosTypes::lid &NumEntries,
                                      std::vector<real64> &vecValues,
                                      std::vector<trilinosTypes::lid> &vecIndices ) const
{
  real64* Values;
  trilinosTypes::lid* Indices;
  m_matrix->ExtractMyRowView( localRow, NumEntries, Values, Indices );
  vecIndices.assign( Indices, Indices+NumEntries );
  vecValues.assign( Values, Values+NumEntries );
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
  return m_matrix.get();
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get number of global rows.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Accessor for the number of global rows

trilinosTypes::gid EpetraMatrix::globalRows() const
{
  return m_matrix->NumGlobalRows64();
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get number of global columns.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Accessor for the number of global columns

trilinosTypes::gid EpetraMatrix::globalCols() const
{
  return m_matrix->NumGlobalCols64();
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get the number of local rows.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Accessor for the number of local rows
trilinosTypes::lid EpetraMatrix::localRows() const
{
  return m_matrix->NumMyRows();
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get the number of local columns.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Accessor for the number of local columns
trilinosTypes::lid EpetraMatrix::localCols() const
{
  return m_matrix->NumMyCols();
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get the lower index owned by processor.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Accessor for the index of the first global row
trilinosTypes::gid EpetraMatrix::ilower() const
{
  return m_matrix->RowMap().MyGlobalElements64()[0];
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get the upper index owned by processor.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Accessor for the index of the last global row
trilinosTypes::gid EpetraMatrix::iupper() const
{
  return m_matrix->RowMap().MyGlobalElements64()[0] + m_matrix->RowMap().NumMyElements();
}


// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Global/Local index mapping
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Given a global index return the local index if owned by this processor (or vice versa).
// Returns -1 if not owned.

trilinosTypes::lid EpetraMatrix::rowLID( trilinosTypes::gid const GID ) const
{
  return m_matrix->LRID( GID );
}

trilinosTypes::gid EpetraMatrix::rowGID( trilinosTypes::lid const LID ) const
{
  return m_matrix->GRID64( LID );
}


// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Print to terminal.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Wrapper to print the trilinos output of the matrix
void EpetraMatrix::print() const
{
  if( m_matrix.get() != nullptr )
  {
    std::cout << *m_matrix << std::endl;
  }
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Inf-norm.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Returns the infinity norm of the matrix.
real64 EpetraMatrix::normInf() const
{
  return m_matrix->NormInf();
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// 1-norm.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Returns the one norm of the matrix.
real64 EpetraMatrix::norm1() const
{
  return m_matrix->NormOne();
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Frobenius-norm.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Returns the Frobenius norm of the matrix.
real64 EpetraMatrix::normFrobenius() const
{
  return m_matrix->NormFrobenius();
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Is-assembled.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Boolean indicator. True = matrix assembled and ready to be used.
bool EpetraMatrix::isAssembled() const
{
  return assembled;
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
      if( Graph_.FindtrilinosTypes::gidLoc( locRow, Index, j, Loc ))
//  #ifdef EPETRA_HAVE_OMP
//  #ifdef EPETRA_HAVE_OMP_NONASSOCIATIVE
//  #pragma omp atomic
//  #endif
//  #endif
//          RowValues[Loc] += srcValues[j];
        RAJA::atomic::atomicAdd<ATOMIC_POL2>( &RowValues [Loc], srcValues[j] );
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
          RAJA::atomic::atomicAdd<ATOMIC_POL2>( &RowValues [Loc], srcValues[j] );
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
            RAJA::atomic::atomicAdd<ATOMIC_POL2>( &RowValues [Loc], srcValues[j] );
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
          RAJA::atomic::atomicAdd<ATOMIC_POL2>( &RowValues [Loc], srcValues[j] );
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
