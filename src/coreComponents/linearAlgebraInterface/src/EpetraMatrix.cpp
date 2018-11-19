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
EpetraMatrix::EpetraMatrix( EpetraMatrix const &in_mat )
{
  // Check if the vector to be copied is not empty
  if( in_mat.unwrappedPointer() != nullptr )
  {
    // Create a unique pointer to an Epetra_CrsMatrix. The data from the input matrix is
    // copied to a new memory location.
    m_matrix = std::unique_ptr<Epetra_CrsMatrix>( new Epetra_CrsMatrix( *in_mat.unwrappedPointer() ) );
  }
}

// -----------------------------
// Create/Finalize
// -----------------------------
// Allocate matrix (prepare to be filled with data).

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Create a matrix from number of elements
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
void EpetraMatrix::create( MPI_Comm const comm,
                                 trilinosTypes::gid const size,
                                 trilinosTypes::lid const nMaxEntriesPerRow )
{
  // Create an Epetra_Map of size size.
  Epetra_Map map = Epetra_Map( size, 0, Epetra_MpiComm( comm ) );
  // Create a unique pointer to an Epetra_CrsMatrix defined from the Epetra_Map.
  // Specifies the maximum number of non-zeros per row.
  m_matrix = std::unique_ptr<Epetra_CrsMatrix>( new Epetra_CrsMatrix( Copy, map, nMaxEntriesPerRow, false ) );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Create a rectangular matrix from number of elements
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// (TODO this has not been tested
// and from previous experience it may not be the correct way of getting a rectangular
// matrix).
void EpetraMatrix::create( MPI_Comm const comm,
                                 trilinosTypes::gid const in_m_nRowGlobal,
                                 trilinosTypes::gid const in_m_nColGlobal,
                                 trilinosTypes::lid const nMaxEntriesPerRow )
{
  // Create an Epetra_Map of size in_m_nRowGlobal.
  Epetra_Map rowMap = Epetra_Map( in_m_nRowGlobal, 0, Epetra_MpiComm( comm ));
  // Create an Epetra_Map of size in_m_nColGlobal.
  Epetra_Map colMap = Epetra_Map( in_m_nColGlobal, 0, Epetra_MpiComm( comm ));
  // Create a unique pointer to an Epetra_CrsMatrix defined from the two Epetra_Maps.
  // Specifies the maximum number of non-zeros per row.
  m_matrix = std::unique_ptr<Epetra_CrsMatrix>( new Epetra_CrsMatrix( Copy, rowMap, colMap, nMaxEntriesPerRow, false ) );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Create a matrix from Epetra_Map
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
void EpetraMatrix::create( Epetra_Map const &input_map,
                                 trilinosTypes::lid const nMaxEntriesPerRow )
{
  // Create a unique pointer to an Epetra_CrsMatrix defined from the input Epetra_Map.
  // Specifies the maximum number of non-zeros per row.
  m_matrix = std::unique_ptr<Epetra_CrsMatrix>( new Epetra_CrsMatrix( Copy, input_map, nMaxEntriesPerRow, false ) );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Create a matrix from two Epetra_Maps
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// (TODO this has not been tested
// and from previous experience it may not be the correct way of getting a rectangular
// matrix).
void EpetraMatrix::create( Epetra_Map const &input_row_map,
                                 Epetra_Map const &input_col_map,
                                 trilinosTypes::lid const nMaxEntriesPerRow )
{
  // Create a unique pointer to an Epetra_CrsMatrix defined from the two input Epetra_Map.
  // Specifies the maximum number of non-zeros per row.
  m_matrix = std::unique_ptr<Epetra_CrsMatrix>( new Epetra_CrsMatrix( Copy, input_row_map, input_col_map, nMaxEntriesPerRow, false ) );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Create a matrix from an Epetra_CrsGraph.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
void EpetraMatrix::create( Epetra_CrsGraph const &graph )
{
  // Create a unique pointer to an Epetra_CrsMatrix defined from the input Epetra_CrsGraph.
  // An Epetra_CrsGraph can be used as a sparsity pattern.
  m_matrix = std::unique_ptr<Epetra_CrsMatrix>( new Epetra_CrsMatrix( Copy, graph ) );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Create a matrix from an Epetra_CrsMatrix.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
void EpetraMatrix::create( Epetra_CrsMatrix &in_matrix )
{
  // Create a unique pointer to an Epetra_CrsMatrix defined from the input Epetra_CrsMatrix.
  // Similar to the copy constructor but does no copy the data, just points to the input matrix.
  m_matrix = std::unique_ptr<Epetra_CrsMatrix>( &in_matrix );
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

  m_matrix->FillComplete();
  assembled = true;
}

// -------------------------
// Add/Set
// -------------------------

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Add single value.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Add a value at row iRow and column iCol
void EpetraMatrix::add( trilinosTypes::gid const iRow,
                              trilinosTypes::gid const iCol,
                              real64 const value )
{
  m_matrix->SumIntoGlobalValues( iRow, 1, &value, &iCol );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Add row values.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Add row values at row iRow and columns cols (size nCols)
void EpetraMatrix::add( trilinosTypes::gid const iRow,
                              trilinosTypes::lid const nCols,
                              real64 const *values,
                              trilinosTypes::gid const *cols )
{

#if 1
  m_matrix->SumIntoGlobalValues( iRow, nCols, values, cols );
#else

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
#endif
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Add dense matrix.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// TODO THIS FUNCTION HAS NOT BEEN PROPERLY TESTED.
void EpetraMatrix::add( array1d<trilinosTypes::gid> const rowIndices,
                              array1d<trilinosTypes::gid> const colIndices,
                              array2d<real64> const values)
{
  // Dummy variable for a row of the matrix values.
//  array1d<real64> valuesInRow;

  // Loop over the rows
  for ( integer i = 0; i < rowIndices.size(); i++ )
  {

    // Loop over the columns to extract the values in the row (TODO there is probably
    // a better way to access a row of an array2d object).
//    for (integer j = 0; j < colIndices.size(); j++)
//    {
//      valuesInRow[j] = values[i][j];
//    }

//      Fill row i.
//    m_matrix->SumIntoGlobalValues( rowIndices[i], colIndices.size(), valuesInRow.data(), colIndices.data() );

    // This version uses the sliced array2d. TODO test it properly! But should be more efficient.
    m_matrix->SumIntoGlobalValues( rowIndices[i], static_cast<integer>(colIndices.size()), values[i], colIndices.data() );
  }
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Set single value.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Set the value of the element (iRow,iCol) to value.
void EpetraMatrix::set( trilinosTypes::gid const iRow,
                              trilinosTypes::gid const iCol,
                              real64 const value )
{
  m_matrix->ReplaceGlobalValues( iRow, 1, &value, &iCol );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Set row values.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Add the value of elements (iRow,[nCols]) to values.
void EpetraMatrix::set( trilinosTypes::gid const iRow,
                              trilinosTypes::lid const nCols,
                              real64 const *values,
                              trilinosTypes::gid const *cols )
{
  m_matrix->ReplaceGlobalValues( iRow, nCols, values, cols );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Insert values.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// TODO remove the possibility to dynamically construct the sparsity pattern.
// Insert the value to the element (iRow,iCol). This routine should not be used unless
// there is a guarantee that each element will only be accessed once. If that is not the
// case, it will allocate (way) too much memory and sum all the values allocated to the
// same element.
void EpetraMatrix::insert( trilinosTypes::gid const iRow,
                                 trilinosTypes::lid const nCols,
                                 real64 const *values,
                                 trilinosTypes::gid const *cols )
{
  m_matrix->InsertGlobalValues( iRow, nCols, values, cols );
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
// Left scale (diagonal scaling).
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Pre-multiplies (left) with diagonal matrix consisting of the values in vec.
void EpetraMatrix::leftScale( EpetraVector const &vec )
{
  m_matrix->LeftScale( *(*vec.unwrappedPointer())(0) );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Right scale (diagonal scaling)
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Post-multiplies (right) with diagonal matrix consisting of the values in vec.
void EpetraMatrix::rightScale( EpetraVector const &vec )
{
  m_matrix->RightScale( *(*vec.unwrappedPointer())(0) );
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Left and Right scale (diagonal scalings)
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Pre-multiplies (left) with diagonal matrix consisting of the values in vecLeft and
// Post-multiplies (right) with diagonal matrix consisting of the values in vecRight.
void EpetraMatrix::leftRightScale( EpetraVector const &vecLeft,
                                   EpetraVector const &vecRight )
{
  leftScale(vecLeft);
  rightScale(vecRight);
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Clear row.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Clear row and multiply diagonal term by factor.
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
  // Same as the other getRow function but outputs the values and column indices to
  // standard vectors.
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
  // Same as the other getLocalRow function but outputs the values and column indices to
  // standard vectors.
  real64* Values;
  trilinosTypes::lid* Indices;
  m_matrix->ExtractMyRowView( localRow, NumEntries, Values, Indices );
  vecIndices.assign( Indices, Indices+NumEntries );
  vecValues.assign( Values, Values+NumEntries );
}

// ----------------------------
//  Accessors
// ----------------------------

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get pointer.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Accessor for the pointer to the matrix
Epetra_CrsMatrix * EpetraMatrix::unwrappedPointer() const
{
  return m_matrix.get();
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get number of global rows.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Accessor for the number of global rows
trilinosTypes::gid EpetraMatrix::globalRows() const
{
  return m_matrix->RowMap().NumGlobalElements64();
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get number of global columns.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Accessor for the number of global columns
trilinosTypes::gid EpetraMatrix::globalCols() const
{
  return m_matrix->ColMap().NumGlobalElements64();
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get number of unique columns.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Accessor for the number of global columns
trilinosTypes::gid EpetraMatrix::uniqueCols() const
{
  return m_matrix->DomainMap().NumGlobalElements64();
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
// Get the Epetra row map.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Accessor for the row map
Epetra_Map const & EpetraMatrix::RowMap() const
{
  return m_matrix->RowMap();
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get the Epetra column map.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Accessor for the column map.
Epetra_Map const & EpetraMatrix::ColMap() const
{
  return m_matrix->ColMap();
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get the Epetra domain map.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Accessor for the domain map
Epetra_Map const & EpetraMatrix::DomainMap() const
{
  return m_matrix->DomainMap();
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get the number of local rows.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Accessor for the number of local rows
int EpetraMatrix::myRows() const
{
  return m_matrix->NumMyRows();
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Get the number of local columns.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Accessor for the number of local columns
int EpetraMatrix::myCols() const
{
  return m_matrix->NumMyCols();
}

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Test for ownership of global index.
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Returns the local index if owned by this processor. Returns -1 if
// not owned.
trilinosTypes::lid EpetraMatrix::rowMapLID( trilinosTypes::gid const GID ) const
{
  return m_matrix->RowMap().LID( GID );
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

}

// END_RST_NARRATIVE
